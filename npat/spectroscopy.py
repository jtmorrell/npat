from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sqlite3, os, json
import numpy as np
import pandas as pd
import datetime as dtm
import multiprocessing

from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.interpolate import interp1d

from .isotope import Isotope
from .plotter import colors
from .plotter import _init_plot
from .plotter import _close_plot

global DB_CONNECTIONS_DICT
DB_CONNECTIONS_DICT = {}

class Calibration(object):
	"""Calibration is a class for calibrating
	gamma ray data from High-Purity Germanium (HPGe)
	detectors.

	...

	Parameters
	----------
	filename : str (optional)
		Path to .json file with stored calibration data.

	Attributes
	----------
	engcal : np.ndarray
		2 or 3 parameter energy calibration.
	effcal : np.ndarray
		3 or 5 parameter efficiency calibration.
	unc_effcal : np.ndarray
		3x3 or 5x5 dimensional covariance matrix on the efficiency calibration.
	rescal : np.ndarray
		1 or 2 parameter resolution (shape) calibration.

	Methods
	-------
	plot(show=True, fit=True, saveas=None, zoom=None, logscale=True, 
			grayscale=False, labels=False, square_fig=False)
		Plots the spectrum. Various options for plotting.
	calibrate(spectra, saveas=None, db=None, auto_calibrate=False)
		Functionality depends on filetype.
	saveas(filename)
	open(filename)
	eng(idx, *cal)
	eff(energy, *c)
	unc_eff(energy, c=None, u=None)
	res(idx, *cal)
	map_idx(energy, *cal)

	Notes
	-----

	References
	----------

	Examples
	--------

	"""
	def __init__(self, filename=None):
		self.calib = {}
		self._calib_data = {}
		self._default = {'engcal': [0.0, 0.3],
						'effcal': [0.331, 0.158, 0.410, 0.001, 1.476],
						'unc_effcal': [[ 5.038e-02,  3.266e-02, -2.151e-02, -4.869e-05, -7.748e-03],
									 [ 3.266e-02,  2.144e-02, -1.416e-02, -3.416e-05, -4.137e-03],
									 [-2.151e-02, -1.416e-02,  9.367e-03,  2.294e-05,  2.569e-03],
									 [-4.869e-05, -3.416e-05,  2.294e-05,  5.411e-07, -1.165e-04],
									 [-7.748e-03, -4.137e-03,  2.569e-03, -1.165e-04,  3.332e-02]],
						'rescal': [2.0, 4e-4]}
		self.update()

		self._path, self._fnm = None, None
		if filename is not None:
			self._path, self._fnm = os.path.split(filename)
			if self._path in ['',' ']:
				self._path = os.getcwd()
			if self._fnm in ['',' ']:
				raise ValueError('Invalid Filename: {}'.format(filename))

			if os.path.exists(os.path.join(self._path, self._fnm)):
				self.open(os.path.join(self._path, self._fnm))
			else:
				raise ValueError('File does not exist: {}'.format(filename))

	def update(self, meta={}):
		for nm in meta:
			if nm.endswith('cal'):
				self.calib[nm] = meta[nm]
		for nm in ['engcal', 'effcal', 'unc_effcal', 'rescal']:
			if nm not in self.calib:
				self.calib[nm] = self._default[nm]

	@property
	def engcal(self):
		return self.calib['engcal']

	@engcal.setter
	def engcal(self, cal):
		self.calib['engcal'] = cal

	@property
	def effcal(self):
		return self.calib['effcal']

	@effcal.setter
	def effcal(self, cal):
		self.calib['effcal'] = cal

	@property
	def unc_effcal(self):
		return self.calib['unc_effcal']

	@unc_effcal.setter
	def unc_effcal(self, cal):
		self.calib['unc_effcal'] = cal

	@property
	def rescal(self):
		return self.calib['rescal']

	@rescal.setter
	def rescal(self, cal):
		self.calib['rescal'] = cal
	
	def eng(self, idx, *cal):
		"""Description

		...

		Parameters
		----------
		idx : array_like
			Description of parameter `x`.
		cal : array_like (len 2 or 3)
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		idx = np.asarray(idx)
		cal = cal if len(cal) else self.calib['engcal']
		if len(cal)==2:
			return cal[0]+cal[1]*idx
		elif len(cal)==3:
			return cal[0]+cal[1]*idx+cal[2]*idx**2

	def eff(self, energy, *c):
		"""Description

		...

		Parameters
		----------
		energy : array_like
			Description of parameter `x`.
		c : array_like (len 3 or 5)
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		energy = np.asarray(energy)
		c = c if len(c) else self.calib['effcal']
		if len(c)==3:
			return c[0]*np.exp(-c[1]*energy**c[2])
		elif len(c)==5:
			return c[0]*np.exp(-c[1]*energy**c[2])*(1.0-np.exp(-c[3]*energy**c[4]))
			
	def unc_eff(self, energy, c=None, u=None):
		"""Description

		...

		Parameters
		----------
		energy : array_like
			Description of parameter `x`.
		c : array_like (len 3 or 5)
			Description of parameter `x`.
		u : array_like (3x3 or 5x5)
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		energy = np.asarray(energy)
		if c is None or u is None:
			c, u = self.calib['effcal'], self.calib['unc_effcal']
		eps = 1E-8
		var= np.zeros(len(energy)) if energy.shape else 0.0
		for n in range(len(c)):
			for m in range(n, len(c)):
				if not np.isfinite(u[n][m]):
					return np.inf
				c_n, c_m = list(c), list(c)
				c_n[n], c_m[m] = c_n[n]+eps, c_m[m]+eps
				par_n = (self.eff(energy, *c_n)-self.eff(energy, *c))/eps
				par_m = (self.eff(energy, *c_m)-self.eff(energy, *c))/eps
				var += u[n][m]*par_n*par_m*(2.0 if n!=m else 1.0)
		return np.sqrt(var)

	def res(self, idx, *cal):
		"""Description

		...

		Parameters
		----------
		idx : array_like
			Description of parameter `x`.
		cal : array_like (len 1 or 2)
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		idx = np.asarray(idx)
		cal = cal if len(cal) else self.calib['rescal']
		if len(cal)==1:
			return cal[0]*np.sqrt(idx)
		elif len(cal)==2:
			return cal[0]+cal[1]*idx

	def map_idx(self, energy, *cal):
		"""Description

		...

		Parameters
		----------
		energy : array_like
			Description of parameter `x`.
		cal : array_like (len 2 or 3)
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		energy = np.asarray(energy)
		cal = cal if len(cal) else self.calib['engcal']
		if len(cal)==3:
			if cal[2]!=0.0:
				return np.array(np.rint(0.5*(np.sqrt(cal[1]**2-4.0*cal[2]*(cal[0]-energy))-cal[1])/cal[2]), dtype=np.int32)
		return np.array(np.rint((energy-cal[0])/cal[1]), dtype=np.int32)

	def calibrate(self, spectra, saveas=None, db=None, auto_calibrate=False):
		"""Description

		...

		Parameters
		----------
		spectra: list
			Description of parameter `x`.
		saveas: 


		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		shelves = []
		sps = []
		
		for sp in spectra:
			if type(sp)==str:
				sp = Spectrum(sp, db)
			sps.append(sp)
			if ('A0' in sp.meta and 'ref_date' in sp.meta):
				if 'shelf' not in sp.meta:
					sp.meta = {'shelf':None}
				shelves.append(sp.meta['shelf'])
		cal_spec = [sp for sp in spectra if ('A0' in sp.meta and 'ref_date' in sp.meta)]
		if len(cal_spec)==0:
			raise ValueError('Cannot calibrate: A0 and ref_date required for calibration.')

		if auto_calibrate:
			for sp in cal_spec:
				sp.auto_calibrate(**{'_cb_cal':False})
		p0_engcal = np.average([s.cb.engcal for s in cal_spec if len(s.cb.engcal)==len(cal_spec[0].cb.engcal)], axis=0)
		engcal = self._calibrate_energy(cal_spec, p0_engcal)
		rescal = self._calibrate_resolution(cal_spec, np.average([s.cb.rescal for s in cal_spec], axis=0))
		for sp in sps:
			sp.meta = {'engcal':engcal, 'rescal':rescal}
		self.update({'engcal':engcal, 'rescal':rescal})

		shelves = {sh:[s for s in cal_spec if s.meta['shelf']==sh] for sh in set(shelves)}
		for shelf in shelves:
			if 'effcal' not in self._calib_data:
				self._calib_data['effcal'] = []

			shelf_specs = shelves[shelf]
			p0_effcal = np.average([s.cb.effcal for s in shelf_specs if len(s.cb.effcal)==len(shelf_specs[0].cb.effcal)], axis=0)
			effcal, unc_effcal = self._calibrate_efficiency(shelf_specs, p0_effcal)
			for sp in sps:
				if shelf==sp.meta['shelf']:
					sp.meta = {'effcal':effcal, 'unc_effcal':unc_effcal}
					self.update({'effcal':effcal, 'unc_effcal':unc_effcal})
					if saveas is not None:
						if type(saveas)==str:
							saveas = [saveas]
						sp.save(*saveas)

	def _calibrate_energy(self, spectra, p0):
		E, chn, unc_chn = [], [], []
		for sp in spectra:
			sp.meta = {'engcal':p0}
			for f in sp.fits:
				E += f['gm'][:,0].tolist()
				chn += [f['fit'][i] for i in f['ix']['mu']]
				unc_chn += [np.sqrt(f['unc'][i][i]) for i in f['ix']['mu']]

		E, chn, unc_chn = np.array(E), np.array(chn), np.array(unc_chn)
		idx = np.where((chn>unc_chn)&(np.isfinite(unc_chn))&(unc_chn>1E-3))
		E, chn, unc_chn = E[idx], chn[idx], unc_chn[idx]

		p0 = p0 if len(p0)==3 else [p0[0], p0[1], 0.0]
		fit, unc = curve_fit(self.eng, chn, E, sigma=self.eng(unc_chn), p0=p0)

		self._calib_data['engcal'] = {'fit':fit, 'unc':unc, 'x':E, 'y':chn, 'yerr':unc_chn}
		return fit

	def _calibrate_efficiency(self, spectra, p0):
		E, eff, unc_eff = [], [], []
		lms = {}
		for sp in spectra:
			sp.meta = {'effcal':p0}
			p = sp.peaks
			E += p['energy'].tolist()
			if type(sp.meta['A0'])==list:
				A0 = np.array([sp.meta['A0'][sp.meta['istp'].index(i)] for i in p['isotope']])
			else:
				A0 = np.array([sp.meta['A0'] for i in p['isotope']])
			if type(sp.meta['ref_date'])==list:
				rd = [sp.meta['ref_date'][sp.meta['istp'].index(i)] for i in p['isotope']]
				d0 = [dtm.datetime.strptime(r, '%m/%d/%Y %H:%M:%S') for r in rd]
				t_d = np.array([(sp.meta['start_time']-d).total_seconds() for d in d0])
			else:
				rd = sp.meta['ref_date']
				d0 = dtm.datetime.strptime(rd, '%m/%d/%Y %H:%M:%S')
				t_d = np.array([(sp.meta['start_time']-d0).total_seconds() for i in p['isotope']])

			lm = []
			for i in p['isotope']:
				if i not in lms:
					lms[i] = Isotope(i).decay_const(unc=True)
				lm.append(lms[i])
			lm = np.array(lm)

			ef = p['counts']*lm[:,0]/((1.0-np.exp(-lm[:,0]*sp.meta['real_time']))*np.exp(-lm[:,0]*t_d)*p['intensity']*A0*(sp.meta['live_time']/sp.meta['real_time']))
			unc_eff += (ef*np.sqrt((p['unc_counts']/p['counts'])**2+(p['unc_intensity']/p['intensity'])**2+(lm[:,1]/lm[:,0])**2+0.01**2)).tolist()
			eff += ef.tolist()

		E, eff, unc_eff = np.array(E), np.array(eff), np.array(unc_eff)
		idx = np.where((eff>unc_eff)&(unc_eff>0.0)&(np.isfinite(unc_eff)))
		E, eff, unc_eff = E[idx], eff[idx], unc_eff[idx]
		p0 = p0.tolist() if len(p0)==5 else p0.tolist()+[0.001, 1.476]
		p0[0] = max([min([p0[0]*np.average(np.array(eff)/self.eff(np.array(E), *p0), weights=(self.eff(np.array(E), *p0)/np.array(unc_eff))**2),99.9]),0.001])

		bounds = ([0.0, 0.0, -1.0, 0.0, -2.0], [100.0, 3.0, 3.0, 0.5, 3.0])
		if any([sp.fit_config['xrays'] for sp in spectra]):
			try:
				fit3, unc3 = curve_fit(self.eff, E, eff, sigma=unc_eff, p0=p0[:3], bounds=(bounds[0][:3], bounds[1][:3]))
				fit5, unc5 = curve_fit(self.eff, E, eff, sigma=unc_eff, p0=fit3.tolist()+p0[3:], bounds=bounds)
				
				chi5 = np.sum((eff-self.eff(E, *fit5))**2/unc_eff**2)
				chi3 = np.sum((eff-self.eff(E, *fit3))**2/unc_eff**2)
				## Invert to find which is closer to one
				chi5 = chi5 if chi5>1.0 else 1.0/chi5
				chi3 = chi3 if chi3>1.0 else 1.0/chi3
				fit, unc = (fit3, unc3) if chi3<=chi5 else (fit5, unc5)
			except:
				fit, unc = curve_fit(self.eff, E, eff, sigma=unc_eff, p0=p0[:3], bounds=(bounds[0][:3], bounds[1][:3]))

		else:
			fit, unc = curve_fit(self.eff, E, eff, sigma=unc_eff, p0=p0[:3], bounds=(bounds[0][:3], bounds[1][:3]))

		self._calib_data['effcal'].append({'fit':fit, 'unc':unc, 'x':E, 'y':eff, 'yerr':unc_eff, 'shelf':spectra[0].meta['shelf']})
		return fit, unc

	def _calibrate_resolution(self, spectra, p0):
		mu, sig, unc_sig = [], [], []
		for sp in spectra:
			sp.meta = {'rescal':p0}
			for f in sp.fits:
				mu += [f['fit'][i] for i in f['ix']['mu']]
				sig += [f['fit'][i] for i in f['ix']['sig']]
				u_sig = [f['unc'][i][i] for i in f['ix']['sig']]
				u_mu = [f['unc'][i][i] for i in f['ix']['mu']]
				unc_sig += [np.sqrt(s+0.25*p0[0]*u_mu[n]/mu[-1*(len(u_mu)-n)]) for n,s in enumerate(u_sig)]
		mu, sig, unc_sig = np.array(mu), np.array(sig), np.array(unc_sig)
		idx = np.where((sig>unc_sig)&(np.isfinite(unc_sig))&(unc_sig>1E-3))
		mu, sig, unc_sig = mu[idx], sig[idx], unc_sig[idx]
		fit, unc = curve_fit(self.res, mu, sig, sigma=unc_sig, p0=p0, bounds=([-np.inf, 0.0],[np.inf, np.inf]))

		self._calib_data['rescal'] = {'fit':fit, 'unc':unc, 'x':mu, 'y':sig, 'yerr':unc_sig}
		return fit

	def open(self, filename):
		"""Description

		...

		Parameters
		----------
		x : type
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		if filename.endswith('.json'):
			with open(filename) as f:
				js = json.load(f)

			for c in ['engcal','rescal','effcal','unc_effcal']:
				if c in js:
					self.update({c:np.array(js[c])})
			if '_calib_data' in js:
				for c in ['engcal','rescal']:
					if c in js['_calib_data']:
						self._calib_data[c] = {str(i):np.array(js['_calib_data'][c][i]) for i in js['_calib_data'][c]}
				if 'effcal' in js['_calib_data']:
					self._calib_data['effcal'] = []
					for e in js['_calib_data']['effcal']:
						self._calib_data['effcal'].append({str(i):(np.array(e[i]) if i!='shelf' else (str(e[i]) if e[i] is not None else None)) for i in e})


	def saveas(self, filename):
		"""Description

		...

		Parameters
		----------
		x : type
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		if filename.endswith('.json'):
			js = {}
			for c in ['engcal','rescal','effcal','unc_effcal']:
				if type(self.calib[c])==list:
					js[c] = self.calib[c]
				else:
					js[c] = self.calib[c].tolist()
			js['_calib_data'] = {}
			for cl in ['engcal','rescal']:
				if cl in self._calib_data:
					js['_calib_data'][cl] = {i:self._calib_data[cl][i].tolist() for i in self._calib_data[cl]}
			if 'effcal' in self._calib_data:
				js['_calib_data']['effcal'] = []
				for e in self._calib_data['effcal']:
					js['_calib_data']['effcal'].append({i:(e[i].tolist() if i!='shelf' else e[i]) for i in e})

			with open(filename, 'w') as f:
				json.dump(js, f, indent=4)

	def plot(self, **kwargs):
		"""Description

		...

		Parameters
		----------
		x : type
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		cm = colors(aslist=True)
		cm_light = colors(shade='light', aslist=True)
		f, ax = _init_plot(N_plots=3, figsize=(12.8, 4.8))

		if 'engcal' in self._calib_data:
			d = self._calib_data['engcal']
			ax[0].errorbar(d['x'], d['y'], yerr=d['yerr'], ls='None', marker='o')
			x = np.arange(min(d['x']), max(d['x']), 0.1)
			ax[0].plot(x, self.map_idx(x, *d['fit']))
			ax[0].set_xlabel('Energy (keV)')
			ax[0].set_ylabel('ADC Channel')
			ax[0].set_title('Energy Calibration')

		if 'rescal' in self._calib_data:
			d = self._calib_data['rescal']
			ax[1].errorbar(d['x'], d['y'], yerr=d['yerr'], ls='None', marker='o')
			x = np.arange(min(d['x']), max(d['x']), 0.1)
			ax[1].plot(x, self.res(x, *d['fit']))
			ax[1].set_xlabel('ADC Channel')
			ax[1].set_ylabel('Peak Width')
			ax[1].set_title('Resolution Calibration')
		
		if 'effcal' in self._calib_data:
			for N,d in enumerate(sorted(self._calib_data['effcal'], key=lambda h: h['shelf'])):
				ax[2].errorbar(d['x'], d['y'], yerr=d['yerr'], ls='None', marker='o', color=cm[N%len(cm)])
				x = np.arange(min(d['x']), max(d['x']), 0.1)
				ax[2].plot(x, self.eff(x, *d['fit']), color=cm[N%len(cm)], label=d['shelf'])
				low = self.eff(x, *d['fit'])-self.unc_eff(x, d['fit'], d['unc'])
				high = self.eff(x, *d['fit'])+self.unc_eff(x, d['fit'], d['unc'])
				ax[2].fill_between(x, low, high, facecolor=cm_light[N%len(cm)], alpha=0.5)
			ax[2].set_xlabel('Energy (keV)')
			ax[2].set_ylabel('Efficiency')
			ax[2].set_title('Efficiency Calibration')
			if any([d['shelf'] for d in self._calib_data['effcal']]):
				ax[2].legend(loc=0)

		return _close_plot(f, None, default_log=False, **kwargs)







def _parallel_fit(dat):
	x, y, fit_config, snip, p0, bounds = dat['x'], dat['y'], dat['fit_config'], dat['snip'], dat['p0'], dat['bounds']

	def _multiplet(x, *args):

		def _peak(x, A, mu, sig, R, alpha, step):
			r2 = 1.41421356237
			return A*np.exp(-0.5*((x-mu)/sig)**2)+R*A*np.exp((x-mu)/(alpha*sig))*erfc((x-mu)/(r2*sig)+1.0/(r2*alpha))+step*A*erfc((x-mu)/(r2*sig))

		if fit_config['bg_fit']:
			if fit_config['quad_bg']:
				b, peak = 3, args[0]+args[1]*x+args[2]*x**2
			else:
				b, peak = 2, args[0]+args[1]*x
		else:
			b, peak = 0, snip.copy()
		R, alpha, step = fit_config['R'], fit_config['alpha'], fit_config['step']
		if fit_config['skew_fit']:
			if fit_config['step_fit']:
				for n in range(int((len(args)-b)/6)):
					peak += _peak(x,*args[6*n+b:6*n+6+b])
			else:
				for n in range(int((len(args)-b)/5)):
					peak += _peak(x,*(args[5*n+b:5*n+5+b]+(step,)))
		else:
			if fit_config['step_fit']:
				for n in range(int((len(args)-b)/4)):
					peak += _peak(x,*(args[4*n+b:4*n+b+3]+(R, alpha)+(args[4*n+3+b])))
			else:
				for n in range(int((len(args)-b)/3)):
					peak += _peak(x,*(args[3*n+b:3*n+b+3]+(R,alpha,step)))

		return peak

	try:
		fit, unc = curve_fit(_multiplet, x, y, p0=p0, bounds=bounds, sigma=np.sqrt(y+0.1))
		converged = True
	except:
		fit, unc = np.array(p0), np.inf*np.ones((len(p0), len(p0)))
		converged = False
	dat['fit'], dat['unc'], dat['converged'] = fit, unc, converged
	return dat




class Spectrum(object):
	"""Spectrum is a class for fitting and plotting
	gamma ray data from High-Purity Germanium (HPGe)
	detectors.

	...


	Parameters
	----------
	filename : str
		Path to .Spe file.
	db : str, optional
		Path to sqlite database

	Attributes
	----------
	fits : list
		List of peaks (PeakFit class) found in spectrum fit.
	meta : dict
		Metadata about spectrum. 
	hist : np.ndarray
		1D histogram of pulse amplitudes (energy).

	Methods
	-------
	plot(show=True, fit=True, saveas=None, zoom=None, logscale=True, 
			grayscale=False, labels=False, square_fig=False)
		Plots the spectrum. Various options for plotting.
	save(*saveas)
		Functionality depends on filetype.
	summarize(printout=True, saveas=None)
		Prints and/or saves summary of peaks/isotopes.

	Notes
	-----

	References
	----------

	Examples
	--------

	"""

	def __init__(self, filename=None, db=None):
		self._default_params()

		if filename is not None:
			self._path, self._fnm = os.path.split(filename)
			if self._path in ['',' ']:
				self._path = os.getcwd()
			if self._fnm in ['',' ']:
				raise ValueError('Invalid Spectrum Filename: {}'.format(filename))

			if os.path.exists(os.path.join(self._path, self._fnm)):
				print('Reading Spectrum {}'.format(self._fnm))
				if self._fnm.endswith('.Spe'):
					self._from_Spe()
				elif self._fnm.endswith('.Chn'):
					self._from_Chn()
			else:
				raise ValueError('File does not exist: {}'.format(filename))

		self._check_db(db)
		self.cb.update(self._meta)


	def _default_params(self):
		self._path, self._fnm = None, None
		self.db, self.db_connection, self.db_fnm = None, None, None
		self._hist = np.zeros(2, dtype=np.int32)
		self._meta = {'fit_config':{'snip_adj':1.0, 'R':0.1, 'alpha':0.9, 
									'step':0.0, 'bg_fit':False, 'skew_fit':False, 
									'step_fit':False, 'quad_bg':False, 'SNR_cut':4.0, 
									'A_bound':10.0, 'mu_bound':1.5, 'sig_bound':1.5,
									'xrays':False, 'pk_width':7.5, 'E_min':75.0, 
									'I_min':0.05,'threads':multiprocessing.cpu_count()},
						'istp':[], 'shelf':None, 'spec_id':1}
		self._cb = Calibration()
		for nm in self.cb.calib:
			self._meta[nm] = self.cb.calib[nm]
		self._gamma_list = None
		self._fits, self._peaks = None, None

	def _get_db_params(self):
		if self.db is not None:
			from ast import literal_eval
			tables = [str(i[0]) for i in self.db.execute('SELECT name FROM sqlite_master WHERE type="table"')]

			if 'spectra' not in tables:
				q = 'CREATE TABLE `spectra` (`spec_id` INTEGER, `filename` TEXT, `file_path` TEXT, '
				q += '`start_time` TEXT, `live_time` REAL, `real_time` REAL, `meta` TEXT, `fit_config` TEXT, '
				q += '`engcal` TEXT, `effcal` TEXT, `unc_effcal` TEXT, `rescal` TEXT);'
				self.db.execute(q)
				self.db_connection.commit()

			q = list(self.db.execute('SELECT * FROM spectra WHERE filename LIKE ? AND file_path LIKE ?', ('%'+self._fnm+'%', '%'+self._path+'%')))
			if len(q):
				self.meta = {'engcal':literal_eval(str(q[0][8])), 'effcal':literal_eval(str(q[0][9])), 
							'unc_effcal':np.array(literal_eval(str(q[0][10]).replace('inf',"'inf'")), dtype=np.float64), 
							'rescal':literal_eval(str(q[0][11]))}
				self._meta['spec_id'] = q[0][0]
			else:
				ids = [i[0] for i in self.db.execute('SELECT spec_id FROM spectra')]
				self._meta['spec_id'] = max(ids)+1 if len(ids) else 1
				
				vals = [self.meta['spec_id'], self._fnm, self._path, dtm.datetime.strftime(self.meta['start_time'], '%m/%d/%Y %H:%M:%S'),
						self.meta['live_time'], self.meta['real_time'],
						str({i:self.meta[i] for i in self.meta if i not in ['start_time', 'fit_config','engcal','effcal','unc_effcal','rescal']}),
						str({i:self.fit_config[i] for i in self.fit_config if i not in []}),
						str(list(self.meta['engcal'])), str(list(self.meta['effcal'])), 
						(str(self.meta['unc_effcal']) if type(self.meta['unc_effcal'])==list else str(self.meta['unc_effcal'].tolist())),
						str(list(self.meta['rescal']))]
				self.db.execute('INSERT INTO spectra VALUES('+','.join(12*['?'])+')', tuple(vals))
				self.db_connection.commit()

	def _check_db(self, db=None):
		if db is not None:
			path, fnm = os.path.split(db)
			if path in ['',' ']:
				path = os.getcwd()
			if fnm in ['',' ']:
				raise ValueError('Invalid db Filename: {}'.format(db))
			self.db_fnm = os.path.join(path, fnm)

			global DB_CONNECTIONS_DICT
			if os.path.exists(self.db_fnm):
				if self.db_fnm not in DB_CONNECTIONS_DICT:
					DB_CONNECTIONS_DICT[self.db_fnm] = sqlite3.connect(self.db_fnm)
				self.db_connection = DB_CONNECTIONS_DICT[self.db_fnm]
				self.db = self.db_connection.cursor()
			else:
				print('WARNING: DB {} does not exist, creating new file.'.format(fnm))
				from sqlite3 import Error
				try:
					self.db_connection = sqlite3.connect(self.db_fnm)
					DB_CONNECTIONS_DICT[self.db_fnm] = self.db_connection
					self.db = self.db_connection.cursor()
				except Error as e:
					print(e)
			self._get_db_params()

	def _update_db(self, meta):
		if self.db is not None:
			if meta:
				meta_str = str({i:self.meta[i] for i in self.meta if i not in ['start_time','fit_config','engcal','effcal','unc_effcal','rescal']})
				self.db.execute('UPDATE spectra SET meta=? WHERE spec_id=?', (meta_str, self.meta['spec_id']))
				fit_config_str = str({i:self.fit_config[i] for i in self.fit_config if i not in []})
				self.db.execute('UPDATE spectra SET fit_config=? WHERE spec_id=?', (fit_config_str, self.meta['spec_id']))
				vals = [str(list(self.meta['engcal'])), str(list(self.meta['effcal'])), 
						(str(self.meta['unc_effcal']) if type(self.meta['unc_effcal'])==list else str(self.meta['unc_effcal'].tolist())),
						str(list(self.meta['rescal'])), self.meta['spec_id']]
				self.db.execute('UPDATE spectra SET engcal=?, effcal=?, unc_effcal=?, rescal=? WHERE spec_id=?', tuple(vals))
			elif self._fits is not None:
				if 'peaks' in [str(i[0]) for i in self.db.execute('SELECT name FROM sqlite_master WHERE type="table"')]:
					self.db.execute('DELETE FROM peaks WHERE spec_id=?',(self.meta['spec_id'],))
				self.peaks.to_sql('peaks', self.db_connection, if_exists='append', index=False)
			self.db_connection.commit()

	@property
	def cb(self):
		return self._cb

	@cb.setter
	def cb(self, _cb):
		import copy
		self._cb = copy.deepcopy(_cb)
		if 'effcal' in _cb._calib_data:
			shelf = self.meta['shelf'] if 'shelf' in self.meta else None
			for sh in _cb._calib_data['effcal']:
				if sh['shelf']==shelf:
					self.meta = {'effcal':sh['fit'], 'unc_effcal':sh['unc']}	

	@property
	def hist(self):
		return self._hist

	@hist.setter
	def hist(self, hist_array):
		self._hist = np.asarray(hist_array, dtype=np.int32)
		self.channels = np.arange(len(self._hist))
		self._snip = self.snip_bg()
		self.snip = interp1d(self.channels, self._snip)
		self._fits, self._peaks = None, None


	@property
	def meta(self):
		return self._meta

	@meta.setter
	def meta(self, meta_dict):
		for nm in meta_dict:
			if nm.endswith('cal'):
				self.cb.calib[nm] = meta_dict[nm]
			self._meta[nm] = meta_dict[nm]
			if nm=='istp':
				self._fits, self._peaks, self._gamma_list = None, None, None
		self._update_db(True)
	
	@property
	def fit_config(self):
		return self._meta['fit_config']
	
	@fit_config.setter
	def fit_config(self, config_dict):
		self._fits, self._peaks, self._gamma_list = None, None, None
		for nm in config_dict:
			self._meta['fit_config'][nm] = config_dict[nm]

	def __add__(self, other):
		"""Description

		...

		Parameters
		----------
		x : type
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""
		
		if not type(other)==Spectrum:
			raise ValueError('Cannot add spectrum to type {}'.format(type(other)))

		if self.meta['start_time']==other.meta['start_time']:
			alpha = np.sum(self.hist)/float(np.sum(other.hist))
			dead_time = alpha*(self.meta['real_time']-self.meta['live_time'])+(1.0-alpha)*(other.meta['real_time']-other.meta['live_time'])
			self.meta['real_time'] = alpha*self.meta['real_time']+(1.0-alpha)*other.meta['real_time']
			self.meta['live_time'] = self.meta['real_time']-dead_time
		else:
			self.meta['live_time'] += other.meta['live_time']
			self.meta['real_time'] += other.meta['real_time']

		if len(self.hist)==len(other.hist):
			if len(self.cb.engcal)==len(other.cb.engcal):
				if len(np.where(self.cb.engcal==other.cb.engcal)[0])==len(self.cb.engcal):
					self.hist += other.hist
					return self

		other_bins = other.cb.eng(np.arange(-0.5, len(other.hist)+0.5, 1.0))
		dNdE = np.array(other.hist, dtype=np.float64)/(other_bins[1:]-other_bins[:-1])
		f = interp1d(other.cb.eng(np.arange(len(other.hist))), dNdE, bounds_error=False, fill_value=0.0)

		bins = self.cb.eng(np.arange(-0.5, len(self.hist)+0.5, 1.0))
		edges = np.append(bins, bins).reshape((2,len(bins))).T.flatten()[1:-1].reshape((len(bins)-1, 2)).T
		if np.__version__>='1.16.0':
			e_grid = np.linspace(edges[0], edges[1], num=10).T
		else:
			e_grid = np.array([np.linspace(edges[0][n], edges[1][n], num=10) for n in range(len(edges[0]))])
		N = np.asarray(np.trapz(f(e_grid), e_grid, axis=1), dtype=np.int64)
		self.hist += np.random.poisson(N)

		return self

	def rebin(self, N_bins):
		"""Description

		...

		Parameters
		----------
		x : type
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		L = len(self.hist)
		if N_bins>L:
			raise ValueError('N_bins: {0} must not be greater than current value: {1}'.format(N_bins, L))
		r = int(round(L/float(N_bins)))
		self.hist = np.sum(self.hist.reshape((int(L/r), r)), axis=1)
		ec = self.cb.engcal
		self.meta = {'engcal':[ec[0], ec[1]*r]+([ec[2]*r**2] if len(ec)==3 else [])}
		rc = self.cb.rescal
		self.meta = {'rescal':([rc[0]*np.sqrt(r)] if len(rc)==1 else [rc[0], rc[1]*r])}

	def exp_smooth(self, x, alpha=0.3):
		N = int(2.0/alpha)-1
		wts = np.array([alpha*(1.0-alpha)**abs(N-i) for i in range(2*N+1)])
		y = np.concatenate((np.ones(N)*x[0], x, np.ones(N)*x[-1]))
		return np.dot(wts/np.sum(wts), np.take(y, [np.arange(i,len(x)+i) for i in range(2*N+1)]))

	def snip_bg(self):
		adj = self.fit_config['snip_adj']

		x, dead = self.channels, int(7.5*adj*self.cb.res(len(self._hist)))
		V_i, L = np.log(np.log(np.sqrt(self._hist+1.0)+1.0)+1.0), len(x)-dead
		while self._hist[dead]==0:
			dead += 1
		
		for M in np.linspace(0, 7.5*adj, 10):
			l, h = np.array(x-M*self.cb.res(x), dtype=np.int32), np.array(x+M*self.cb.res(x), dtype=np.int32)
			V_i[dead:L] = np.minimum(V_i[dead:L],0.5*(V_i[l[dead:L]]+V_i[h[dead:L]]))

		snip = (np.exp(np.exp(V_i)-1.0)-1.0)**2-1.0
		snip += adj*1.5*np.sqrt(snip+0.01)

		R1, R2 = 0.3/(adj*self.cb.res(0.1*L)), 0.3/(adj*self.cb.res(0.9*L))
		wt1, wt2 = np.flip(x, 0)/x[-1], x/x[-1]
		return wt1*self.exp_smooth(snip, R1) + wt2*self.exp_smooth(snip, R2)

	@property
	def gamma_list(self):
		if self._gamma_list is None:
			if 'istp' in self._meta:
				xrays, E_min, I_min = self.fit_config['xrays'], self.fit_config['E_min'], self.fit_config['I_min']
				gammas = [Isotope(i).gammas(E_lim=[E_min, self.cb.eng(0.99*len(self._hist))], I_lim=[I_min, None], xrays=xrays) for i in self._meta['istp']]
				self._gamma_list = [(np.array([g['E'], g['I'], g['dI']]).T)*np.array([1.0, 0.01, 0.01]) for g in gammas]
			else:
				self._gamma_list = []
		return self._gamma_list

	def _guess_forward_activity(self, *engcal):
		x, Y, p0_A = self.channels, [], []
		R, alpha, step = self.fit_config['R'], self.fit_config['alpha'], self.fit_config['step']
		clip = np.where(self._hist>self._snip, self._hist-self._snip, 0.0)
		W = self.fit_config['pk_width']
		for gm in self.gamma_list:
			L = len(self._hist)
			y = np.zeros(L)
			if len(gm):
				idx = self.cb.map_idx(gm[:,0], *engcal)
				sig = self.cb.res(idx)
				norm = self.cb.eff(gm[:,0])*gm[:,1]/sig
				idx, sig, norm = idx[idx<L], sig[idx<L], norm[idx<L]
				p0_A.append(np.exp(np.average(np.log((clip[idx]/norm)+1.0), weights=np.sqrt(norm)))-1.0)
				if p0_A[-1]<0:
					p0_A[-1] = 0.0
				l, h = np.array(idx-W*sig, dtype=np.int32), np.array(idx+W*sig, dtype=np.int32)
				for m,i in enumerate(idx):
					if l[m]>0 and h[m]<len(x):
						y[l[m]:h[m]] += self._peak(x[l[m]:h[m]], norm[m], i, sig[m], R, alpha, step)
				Y.append(y)
		return p0_A, np.array(Y)

	def _fit_forward_activity(self):
		p0_A, Y = self._guess_forward_activity()
		try:
			return curve_fit(lambda x, *A: self._snip+np.dot(A, Y), self.channels, self._hist, p0=p0_A)
		except:
			return p0_A, np.ones((len(p0_A), len(p0_A)))*np.inf

	def auto_calibrate(self, guess=[], data=[], **param):
		"""Description

		...

		Parameters
		----------
		x : type
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		from scipy.optimize import differential_evolution
		guess = guess if len(guess) else list(self.cb.engcal)
		if len(data):
			data = np.asarray(data, dtype=np.float64)
			if len(data)==1:
				guess = [0.0, data[0][1]/data[0][0]]
			else:
				N = 2 if len(data)<5 else 3
				M = np.column_stack([data[:,0]**m for m in range(N)])
				b = np.array([np.sum(data[:,0]**m*data[:,1]) for m in range(N)])
				M_inv = np.linalg.inv(np.dot(M.T, M))
				guess = np.dot(M_inv, b).tolist()
		guess = guess if len(guess)==3 else list(guess)+[0.0]
		obj = lambda m: np.average(np.sqrt(np.absolute(self._hist-self._snip-np.dot(*self._guess_forward_activity(*[guess[0], m[0], guess[2]])))))
		self.meta = {'engcal': [guess[0], differential_evolution(obj, [(0.995*guess[1], 1.005*guess[1])]).x[0], guess[2]]}

		if ('A0' in self.meta and 'ref_date' in self.meta) and (param['_cb_cal'] if '_cb_cal' in param else True):
			eff, u_eff, res = self.cb.effcal, self.cb.unc_effcal, self.cb.rescal
			try:
				self.cb.calibrate([self])
			except Exception as e:
				print('Error:', e)
				self.meta = {'engcal':guess, 'effcal':eff, 'unc_effcal':u_eff, 'rescal':res}
				try:
					self.cb.calibrate([self])
				except Exception as e:
					print('Error:', e)
					self.meta = {'engcal':guess, 'effcal':eff, 'unc_effcal':u_eff, 'rescal':res}
					print('WARNING: auto calibrate failed. Setting engcal to guess.')

		elif len(self.meta['istp']) and (param['_cb_cal'] if '_cb_cal' in param else True):
			eng, res = self.cb.engcal, self.cb.rescal
			fit_config = self.fit_config.copy()
			self.fit_config = {'bg_fit':False, 'skew_fit':False}
			calib_data = self.cb._calib_data.copy()
			try:
				self.cb._calib_data['engcal'] = {}
				self.cb._calib_data['rescal'] = {}
				self.meta = {'engcal':self.cb._calibrate_energy([self],eng), 'rescal':self.cb._calibrate_resolution([self],res)}
			except Exception as e:
				print('Error:', e)
				self.meta = {'engcal':guess, 'rescal':res}
				self.cb._calib_data = calib_data
				print('WARNING: auto calibrate failed. Setting engcal to guess.')
			self.fit_config = fit_config
			self._fits, self._peaks = None, None


	def _chi2(self, fit, lh):
		non_zero = np.where(self.hist[lh[0]:lh[1]]>0)
		dof = float(len(non_zero[0])-len(fit)-1)
		if dof==0:
			return np.inf
		resid = self.hist[lh[0]:lh[1]][non_zero]-self.multiplet(self.channels[lh[0]:lh[1]][non_zero], *fit)
		return np.sum(resid**2/self.hist[lh[0]:lh[1]][non_zero])/dof
		
	def _unc_calc(self, fn, p0, cov):
		var, eps = 0.0, 1E-8
		for n in range(len(p0)):
			for m in range(n, len(p0)):
				c_n, c_m = list(p0), list(p0)
				c_n[n], c_m[m] = c_n[n]+eps, c_m[m]+eps
				par_n = (fn(*c_n)-fn(*p0))/eps
				par_m = (fn(*c_m)-fn(*p0))/eps
				var += cov[n][m]*par_n*par_m*(2.0 if n!=m else 1.0)
		return var
	
	def _counts(self, fit, cov, lh):
		cfg = self.fit_config
		L = 3+2*int(cfg['skew_fit'])+int(cfg['step_fit'])
		M = int(cfg['bg_fit'])*(2+int(cfg['quad_bg']))
		N_cts, unc_N = [], []

		skew_fn = lambda A, R, alpha, sig: 2*A*R*alpha*sig*np.exp(-0.5/alpha**2)
		min_skew_fn = lambda A, sig: 2*A*cfg['R']*cfg['alpha']*sig*np.exp(-0.5/cfg['alpha']**2)
		pk_fn = lambda A, sig: 2.506628*A*sig

		for m in range(int((len(fit)-M)/L)):
			i = L*m+M
			p_i = [i, i+2]
			s_i = [i, i+3, i+4, i+2]
			if cfg['skew_fit']:
				N_cts.append(pk_fn(*fit[p_i])+skew_fn(*fit[s_i]))
			else:
				N_cts.append(pk_fn(*fit[p_i])+min_skew_fn(*fit[p_i]))
			if cov is not None:
				if not np.isinf(cov[i][i]):
					if cfg['skew_fit']:
						skew_unc = self._unc_calc(skew_fn, fit[s_i], cov[np.ix_(s_i, s_i)])
					else:
						skew_unc = self._unc_calc(min_skew_fn, fit[p_i], cov[np.ix_(p_i, p_i)])
					unc_N.append(np.sqrt(self._unc_calc(pk_fn, fit[p_i], cov[np.ix_(p_i, p_i)])+skew_unc))
				else:
					unc_N.append(np.inf)
			else:
				unc_N.append(np.inf)
		return np.array(N_cts), np.array(unc_N)
	
	def _decays(self, N, unc_N, gm):
		if 'live_time' not in self.meta or 'real_time' not in self.meta:
			raise ValueError('Need real-time/live-time to compute decay-rate.')

		D = N/(gm[:,1]*self.cb.eff(gm[:,0])*(self.meta['live_time']/self.meta['real_time']))
		unc_D = np.sqrt((unc_N/N)**2+(self.cb.unc_eff(gm[:,0])/self.cb.eff(gm[:,0]))**2+(gm[:,2]/gm[:,1])**2)*D

		A, unc_A = D/self.meta['real_time'], unc_D/self.meta['real_time']
		return D, unc_D, A, unc_A

	def _peak(self, x, A, mu, sig, R, alpha, step):
		r2 = 1.41421356237
		return A*np.exp(-0.5*((x-mu)/sig)**2)+R*A*np.exp((x-mu)/(alpha*sig))*erfc((x-mu)/(r2*sig)+1.0/(r2*alpha))+step*A*erfc((x-mu)/(r2*sig))

	def multiplet(self, x, *args):
		"""Description

		...

		Parameters
		----------
		x : type
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		if self.fit_config['bg_fit']:
			if self.fit_config['quad_bg']:
				b, peak = 3, args[0]+args[1]*x+args[2]*x**2
			else:
				b, peak = 2, args[0]+args[1]*x
		elif len(x)==len(self._snip):
			b, peak = 0, self._snip.copy()
		else:
			b, peak = 0, self.snip(x)
		R, alpha, step = self.fit_config['R'], self.fit_config['alpha'], self.fit_config['step']
		if self.fit_config['skew_fit']:
			if self.fit_config['step_fit']:
				for n in range(int((len(args)-b)/6)):
					peak += self._peak(x,*args[6*n+b:6*n+6+b])
			else:
				for n in range(int((len(args)-b)/5)):
					peak += self._peak(x,*(args[5*n+b:5*n+5+b]+(step,)))
		else:
			if self.fit_config['step_fit']:
				for n in range(int((len(args)-b)/4)):
					peak += self._peak(x,*(args[4*n+b:4*n+b+3]+(R, alpha)+(args[4*n+3+b])))
			else:
				for n in range(int((len(args)-b)/3)):
					peak += self._peak(x,*(args[3*n+b:3*n+b+3]+(R,alpha,step)))
		return peak

	def _get_p0(self, lh, px, gm, istp):
		ix = {'A':[],'mu':[],'sig':[],'R':[],'alpha':[],'step':[]}
		N = int(self.fit_config['bg_fit'])*(2+int(self.fit_config['quad_bg']))
		if self.fit_config['bg_fit']:
			M = np.column_stack([self.channels[lh[0]:lh[1]]**m for m in range(N)])
			b = np.array([np.sum(self.channels[lh[0]:lh[1]]**m*self._snip[lh[0]:lh[1]]) for m in range(N)])
			M_inv = np.linalg.inv(np.dot(M.T, M))
			p0 = np.dot(M_inv, b).tolist()
			resid = self._snip[lh[0]:lh[1]]-np.dot(M, p0)
			unc = np.sqrt(np.abs(M_inv*np.dot(resid.T, resid)/max([float((lh[1]-lh[0])-N),1.0])))
			bounds = [[i-10*unc[n][n] for n,i in enumerate(p0)],[i+10*unc[n][n] for n,i in enumerate(p0)]]
		else:
			p0, bounds = [],[[],[]]

		R, alpha, step = self.fit_config['R'], self.fit_config['alpha'], self.fit_config['step']
		bA, bm, bs = self.fit_config['A_bound'], self.fit_config['mu_bound'], self.fit_config['sig_bound']
		for n,g in enumerate(gm):
			p0 += px[n]
			bounds[0] += [0.0, p0[-2]-p0[-1]*bm, p0[-1]/bs]
			bounds[1] += [p0[-3]*bA, p0[-2]+p0[-1]*bm, p0[-1]*bs]
			ix['A'].append(N)
			ix['mu'].append(N+1)
			ix['sig'].append(N+2)
			N += 3
			if self.fit_config['skew_fit']:
				p0 += [R, alpha]
				bounds[0] += [0.0, 0.5]
				bounds[1] += [1.0, max((2.5, alpha))]
				ix['R'].append(N)
				ix['alpha'].append(N+1)
				N += 2
			if self.fit_config['step_fit']:
				p0.append(step)
				bounds[0].append(0.0)
				bounds[1].append(0.1)
				ix['step'].append(N)
				N += 1

		if self.fit_config['threads']>1:
			return {'p0':p0, 'lh':lh, 'bounds':bounds, 'gm':gm, 'istp':istp, 'ix':ix, 'x':self.channels[lh[0]:lh[1]],
					 'y':self.hist[lh[0]:lh[1]], 'snip':self.snip(self.channels[lh[0]:lh[1]]), 'fit_config':self.fit_config.copy()}
		else:
			return {'p0':p0, 'lh':lh, 'bounds':bounds, 'gm':gm, 'istp':istp, 'ix':ix}

	def _get_multiplets(self):
		p0_A, cov_p0_A = self._fit_forward_activity()
		pks, p0 = [], {'idx':[],'sig':[],'A':[]}
		W, o = self.fit_config['pk_width'], 0
		for n,gm in enumerate(self.gamma_list):
			if len(gm):
				idx = self.cb.map_idx(gm[:,0])
				sig = self.cb.res(idx)
				A = p0_A[n-o]*self.cb.eff(gm[:,0])*gm[:,1]/sig
				SNR = A/np.sqrt(self._snip[idx])
				cut = np.where((SNR>self.fit_config['SNR_cut'])&(idx>W*sig)&(idx<(len(self._hist)-W*sig)))
				idx, sig, gm_idx, A = idx[cut], sig[cut], np.arange(len(gm))[cut], A[cut]
				p0['idx'].append(idx)
				p0['sig'].append(sig)
				p0['A'].append(A)
				l, h = np.array(idx-W*sig, dtype=np.int32), np.array(idx+W*sig, dtype=np.int32)
				pks.append(np.array([l, h, np.repeat(n, len(l)), gm_idx, np.repeat(o, len(l)), np.arange(len(l)), np.repeat(n, len(l))], dtype=np.int32).T)
			else:
				o += 1
		

		if 'p0' in self.fit_config:
			mp0 = self.fit_config['p0']
			if type(mp0)==dict:
				mp0 = [{i:mp0[i][n] for i in mp0} for n in range(len([mp0[i] for i in mp0][0]))]

			for n,p in enumerate(mp0):
				if 'mu' in p:
					p['mu'] = int(p['mu'])
					if 'E' not in p:
						p['E'] = self.cb.eng(p['mu'])
				elif 'E' in p:
					p['mu'] = self.cb.map_idx(p['E'])
				else:
					raise ValueError('E or mu must be specified in p0.')
				if 'I' not in p:
					p['I'] = -1.0
				if 'unc_I' not in p:
					p['unc_I'] = 1.0
				if 'A' not in p:
					p['A'] = self._hist[p['mu']]-self._snip[p['mu']]
				if 'sig' not in p:
					p['sig'] = self.cb.res(p['mu'])
				if 'istp' not in p:
					p['istp'] = None

				if p['istp'] not in self.meta['istp']:
					self.meta['istp'].append(p['istp'])
					self.gamma_list.append(np.array([[p['E'],p['I'],p['unc_I']]]))
					ig = len(self.gamma_list)-1
				else:
					ig = self.meta['istp'].index(p['istp'])
					self.gamma_list[ig] = np.append(self.gamma_list[ig], [[p['E'],p['I'],p['unc_I']]], axis=0)
				p0['idx'].append([p['mu']])
				p0['sig'].append([p['sig']])
				p0['A'].append([p['A']])
				l, h = p['mu']-W*p['sig'], p['mu']+W*p['sig']
				pks.append(np.array([[l], [h], [len(p0['A'])-1], [len(self.gamma_list[ig])-1], [o], [0], [ig]], dtype=np.int32).T)

		if len(pks)==0:
			return pks

		pks = np.concatenate(pks)
		pks = pks[np.argsort(pks[:,0])]
		pairs = np.where(pks[:-1,1]>pks[1:,0])[0]
		pairs = np.unique(np.concatenate((pairs, pairs+1)))
		if len(pairs):
			groups = [[pairs[0], pairs[1]]]
			for p in pairs[2:]:
				if (pks[p][0]<pks[groups[-1][-1]][1]) and len(groups[-1])<8:
					groups[-1].append(p)
				else:
					groups.append([p])
		pk_pairs = [[pks[i] for i in p] for p in groups] if len(pairs) else []
		pks = [[p] for p in np.delete(pks, pairs, axis=0)]+pk_pairs

		multiplets = []
		for p in pks:
			lh = [p[0][0], p[-1][1]]
			N = [(i[2]-i[4], i[3], i[5], i[6]) for i in p]
			px = [[p0['A'][n][j], p0['idx'][n][j], p0['sig'][n][j]] for n,m,j,k in N]
			gm = np.asarray([self.gamma_list[k][m] for n,m,j,k in N])
			istp = [self.meta['istp'][k] for n,m,j,k in N]
			multiplets.append(self._get_p0(lh, px, gm, istp))

		return multiplets

	def _multi_fit(self, multi):
		p0, lh, bounds, gm = multi['p0'], multi['lh'], multi['bounds'], multi['gm']
		if 'fit' not in multi:
			try:
				fit, unc = curve_fit(self.multiplet, self.channels[lh[0]:lh[1]], 
					self.hist[lh[0]:lh[1]], p0=p0, bounds=bounds, sigma=np.sqrt(self.hist[lh[0]:lh[1]]+0.1))
				converged = True
			except:
				fit, unc = np.array(p0), np.inf*np.ones((len(p0), len(p0)))
				converged = False
		else:
			fit, unc, converged = multi['fit'], multi['unc'], multi['converged']
		
		N, unc_N = self._counts(fit, unc, lh)
		D, unc_D, A, unc_A = self._decays(N, unc_N, gm)
		chi2 = self._chi2(fit, lh)

		cols = ['spec_id','filename','isotope','energy','counts','unc_counts',
				'intensity','unc_intensity','efficiency','unc_efficiency',
				'decays','unc_decays','decay_rate','unc_decay_rate','chi2','converged']

		return [{'fit':fit,'unc':unc,'l':lh[0],'h':lh[1],'gm':gm,'converged':converged,'ix':multi['ix'],'istp':multi['istp']}, 
				pd.DataFrame({'spec_id':self.meta['spec_id'], 'filename':self._fnm, 'isotope':multi['istp'],
								'energy':gm[:,0], 'counts':N, 'unc_counts':unc_N,
								'intensity':gm[:,1], 'unc_intensity':gm[:,2], 'efficiency':self.cb.eff(gm[:,0]),
								'unc_efficiency':self.cb.unc_eff(gm[:,0]), 'decays':D, 'unc_decays':unc_D,
								'decay_rate':A, 'unc_decay_rate':unc_A, 'chi2':chi2, 'converged':converged}, columns=cols)]

	@property
	def fits(self):
		if self._fits is None:
			p0 = self._get_multiplets()
			if len(p0):
				if self.fit_config['threads']>1:
					pool = multiprocessing.Pool(processes=self.fit_config['threads'])
					p0 = pool.map(_parallel_fit, p0)
					pool.close()
					pool.join()
				multiplets = sorted(list(map(self._multi_fit, p0)), key=lambda h:h[0]['l'])
				self._fits = [i[0] for i in multiplets]
				self._peaks = pd.concat([i[1] for i in multiplets], ignore_index=True)
			else:
				self._fits = []
				self._peaks = pd.DataFrame([])
			self._update_db(False)
		return self._fits


	@property
	def peaks(self):
		if self._peaks is None:
			fits = self.fits
		return self._peaks

	def _split_fits(self):
		fits, istp, lh, N_sub = [], [], [], []
		cfg = self.fit_config
		for ft in self.fits:			
			if not ft['converged']:
				N_sub.append(0)
				continue
			f = ft['fit']
			B = int(cfg['bg_fit'])*(2+int(cfg['quad_bg']))
			p0 = f[:B].tolist()
			L = 3+2*int(cfg['skew_fit'])+int(cfg['step_fit'])
			N_sub.append(int((len(f)-B)/L))
			for n in range(N_sub[-1]):
				istp.append(ft['istp'][n])
				fits.append(p0+f[B+n*L:B+(n+1)*L].tolist())
				mu, sig = f[B+n*L+1], f[B+n*L+2]
				lh.append([mu-cfg['pk_width']*sig, mu+cfg['pk_width']*sig])
		itp_set = sorted(list(set(istp)))
		return {i:[[fits[n],lh[n]] for n,ip in enumerate(istp) if ip==i] for i in itp_set}, itp_set, N_sub

	

	def saveas(self, *fnms):
		"""Description

		...

		Parameters
		----------
		x : type
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		for fl in fnms:
			if fl.startswith('*'):
				fl = os.path.join(self._path, '.'.join(self._fnm.split('.')[:-1])+'.'+fl.split('.')[-1])

			if any([fl.endswith(e) for e in ['.png','.pdf','.eps','.pgf','.ps','.raw','.rgba','.svg','.svgz']]):
				self.plot(saveas=fl, show=False)

			if fl.endswith('.Spe'):
				### Maestro ASCII .Spe ###
				self._meta['DATE_MEA'] = [dtm.datetime.strftime(self.meta['start_time'], '%m/%d/%Y %H:%M:%S')]
				self._meta['MEAS_TIM'] = ['{0} {1}'.format(int(self.meta['live_time']), int(self.meta['real_time']))]
				self._meta['ENER_FIT'] = ['{0} {1}'.format(self.cb.engcal[0], self.cb.engcal[1])]
				self._meta['MCA_CAL'] = ['3','{0} {1} {2} keV'.format(self.cb.engcal[0], self.cb.engcal[1], (self.cb.engcal[2] if len(self.cb.engcal)>2 else 0.0))]
				defaults = {'ROI':['0'],'SPEC_REM':['DET# 0','DETDESC# None','AP# Maestro Version 7.01'],'PRESETS':['0'],
							'SHAPE_CAL':['3','0E+00 0E+00 0E+00'],'SPEC_ID':['No sample description was entered.']}
				for d in defaults:
					if d not in self.meta:
						self._meta[d] = defaults[d]
				ss = '\n'.join(['${}:\n'.format(sc)+'\n'.join(self.meta[sc]) for sc in ['SPEC_ID','SPEC_REM','DATE_MEA','MEAS_TIM']])
				ss += '\n$DATA:\n0 {}\n'.format(len(self.hist)-1)
				ss += '\n'.join([' '*(8-len(i))+i for i in map(str, self.hist)])+'\n'
				ss += '\n'.join(['${}:\n'.format(sc)+'\n'.join(self.meta[sc]) for sc in ['ROI','PRESETS','ENER_FIT','MCA_CAL','SHAPE_CAL']])+'\n'
				with open(fl,'w') as f:
					f.write(ss)

			if fl.endswith('.Chn'):
				### Maestro integer .Chn ###
				if 'SPEC_REM' in self.meta:
					det_no = int(self.meta['SPEC_REM'][0].split(' ')[1].strip())
				else:
					det_no = 0
				ss = np.array([-1, det_no, 1], dtype='i2').tobytes()
				st = self.meta['start_time']
				ss += np.array(('0' if st.second<10 else '')+str(st.second), dtype='S2').tobytes()
				ss += np.array([self.meta['real_time']/0.02, self.meta['live_time']/0.02], dtype='i4').tobytes()
				months = {1:'Jan',2:'Feb',3:'Mar',4:'Apr',5:'May',6:'Jun',7:'Jul',8:'Aug',9:'Sep',10:'Oct',11:'Nov',12:'Dec'}
				ss += np.array(('0' if st.day<10 else '')+str(st.day)+months[st.month]+str(st.year)[-2:]+('1' if st.year>1999 else '0'), dtype='S8').tobytes()
				ss += np.array(('0' if st.hour<10 else '')+str(st.hour)+('0' if st.minute<10 else '')+str(st.minute), dtype='S4').tobytes()
				ss += np.array([0, len(self.hist)], dtype='i2').tobytes()
				ss += np.array(self.hist, dtype='i4').tobytes()
				ss += np.array([-102, 0], dtype='i2').tobytes()
				ss += np.array(self.cb.engcal, dtype='f4').tobytes()
				if 'SHAPE_CAL' in self.meta:
					ss += np.array(self.meta['SHAPE_CAL'][1].split(' '), dtype='f4').tobytes()
				else:
					ss += np.array([0,0,0], dtype='f4').tobytes()
				ss += np.zeros(228, dtype='i1').tobytes()
				if 'SPEC_REM' in self.meta:
					L = len(self.meta['SPEC_REM'][1].split('# ')[1])
					L = min((63, L))
					ss += np.array(L, dtype='i1').tobytes()
					ss += np.array(self.meta['SPEC_REM'][1].split('# ')[1][:L], dtype='S{}'.format(L)).tobytes()
				else:
					ss += np.array(0, dtype='i1').tobytes()
				if 'SPEC_ID' in self.meta:
					if self.meta['SPEC_ID'][0]!='No sample description was entered.':
						L = len(''.join(self.meta['SPEC_ID']))
						L = min((63, L))
						ss += np.array(L, dtype='i1').tobytes()
						ss += np.array(''.join(self.meta['SPEC_ID'])[:L])
					else:
						ss += np.array(0, dtype='i1').tobytes()
						ss += np.array('\x00'*63, dtype='S63').tobytes()
				else:
					ss += np.array(0, dtype='i1').tobytes()
					ss += np.array('\x00'*63, dtype='S63').tobytes()
				ss += np.array('\x00'*128, dtype='S128').tobytes()
				with open(fl, 'wb') as f:
					f.write(ss)

			if fl.endswith('.spe'):
				### Radware gf3 .spe ###
				name = self._fnm[:8]+' '*(8-len(self._fnm))
				with open(fl, 'wb') as f:
					ss = np.array(24, dtype=np.uint32).tobytes()
					ss += np.array(name, dtype='c').tobytes()
					ss += np.array([len(self.hist), 1, 1, 1, 24, 4*len(self.hist)], dtype=np.uint32).tobytes()
					ss += np.array(self.hist, dtype=np.float32).tobytes()
					ss += np.array(4*len(self.hist), dtype=np.uint32).tobytes()
					f.write(ss)

			if fl.endswith('.csv'):
				self.peaks.to_csv(fl, index=False)

			if fl.endswith('.db'):
				meta = self.meta
				self._check_db(fl)
				self.meta = {i:meta[i] for i in meta if i.endswith('cal')}
				self._update_db(True)
				self._update_db(False)


	def _from_Spe(self, filename=None):
		if filename is None:
			filename = os.path.join(self._path, self._fnm)
		with open(filename) as f:
			ln = f.readline()
			while ln:
				if ln.startswith('$'):
					section = ln.strip()[1:-1]
					if section=='DATA':
						L = int(f.readline().strip().split(' ')[1])+1
						self.hist = np.fromfile(f, dtype=np.int64, count=L, sep='\n')
					else:
						self._meta[section] = []
				else:
					self._meta[section].append(ln.strip())
				ln = f.readline()

		self._meta['start_time'] = dtm.datetime.strptime(self._meta['DATE_MEA'][0], '%m/%d/%Y %H:%M:%S')
		self._meta['live_time'], self._meta['real_time'] = tuple(map(float, self._meta['MEAS_TIM'][0].split(' ')))
		self._meta['engcal'] = list(map(float, self._meta['MCA_CAL'][-1].split(' ')[:-1]))


	def _from_Chn(self, filename=None):
		if filename is None:
			filename = os.path.join(self._path, self._fnm)
		with open(filename, 'rb') as f:
			det_no = np.frombuffer(f.read(6), dtype='i2')[1]
			sts = np.frombuffer(f.read(2), dtype='S2')[0].decode('utf-8')
			self._meta['real_time'], self._meta['live_time'] = tuple(map(float, 0.02*np.frombuffer(f.read(8), dtype='i4')))
			st = np.frombuffer(f.read(12), dtype='S12')[0].decode('utf-8')
			months = {'jan':'01','feb':'02','mar':'03','apr':'04',
						'may':'05','jun':'06','jul':'07','aug':'08',
						'sep':'09','oct':'10','nov':'11','dec':'12'}
			start_time = '{0}/{1}/{2} {3}:{4}:{5}'.format(months[st[2:5].lower()],st[:2],('20' if st[7]=='1' else '19')+st[5:7],st[8:10],st[10:],sts)
			self._meta['start_time'] = dtm.datetime.strptime(start_time, '%m/%d/%Y %H:%M:%S')
			L = np.frombuffer(f.read(4), dtype='i2')[1]
			self.hist = np.asarray(np.frombuffer(f.read(4*L), dtype='i4'), dtype=np.int64)
			f.read(4)
			self._meta['engcal'] = np.frombuffer(f.read(12),dtype='f4').tolist()

			shape = np.frombuffer(f.read(12), dtype='f4').tolist()
			self._meta['SHAPE_CAL'] = ['3', ' '.join(map(str, shape))]
			f.read(228)
			L = np.frombuffer(f.read(1), dtype='i1')[0]
			if L:
				det_desc = np.frombuffer(f.read(L), dtype='S{}'.format(L))[0].decode('utf-8')
				self._meta['SPEC_REM'] = ['DET# '+str(det_no), 'DETDESC# '+det_desc, 'AP# Maestro Version 7.01']
			if L<63:
				f.read(63-L)
			L = np.frombuffer(f.read(1), dtype='i1')[0]
			if L:
				sample_desc = np.frombuffer(f.read(L),dtype='S{}'.format(L))[0].decode('utf-8')
				self._meta['SPEC_ID'] = [sample_desc]

	def summarize(self):
		"""Description

		...

		Parameters
		----------
		x : type
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		print(self.__str__())

	def __str__(self):
		ss = ''
		for n,p in self.peaks.iterrows():
			ln1 = [p['energy'], p['isotope'], p['intensity']*100.0]
			ln0 = '{1} - {0} keV (I = {2}%)'.format(*ln1)
			ss += ln0 + '\n'
			ss += ''.join(['-']*len(ln0)) + '\n'
			ln2 = [int(p['counts']), (int(p['unc_counts']) if np.isfinite(p['unc_counts']) else np.inf),
					format(p['decays'], '.3e'), format(p['unc_decays'], '.3e'), 
					format(p['decay_rate'], '.3e'), format(p['unc_decay_rate'], '.3e'), 
					format(p['decay_rate']/3.7E4, '.3e'), format(p['unc_decay_rate']/3.7E4, '.3e'), 
					round(p['chi2'],3)]
			ss += 'counts: {0} +/- {1}'.format(ln2[0], ln2[1]) + '\n'
			ss += 'decays: {0} +/- {1}'.format(ln2[2], ln2[3]) + '\n'
			ss += 'activity (Bq): {0} +/- {1}'.format(ln2[4], ln2[5]) + '\n'
			ss += 'activity (uCi): {0} +/- {1}'.format(ln2[6], ln2[7]) + '\n'
			ss += 'chi2/dof: {}'.format(ln2[8]) + '\n'
			ss += '\n'
		return ss


	def plot(self, fit=True, labels=True, snip=False, xcalib=True, **kwargs):
		"""Description

		...

		Parameters
		----------
		x : type
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		
		xgrid = np.array([self.channels-0.5, self.channels+0.5]).T.flatten()
		spec = np.array([self.hist, self.hist]).T.flatten()
		erange = self.cb.eng(xgrid) if xcalib else xgrid

		f, ax = _init_plot(figsize=(12.8, 4.8), **kwargs)
		spec_label = self._fnm if self._fnm is not None else 'Spectrum'
		ax.plot(erange, spec, lw=1.2, zorder=1, label=(spec_label if labels else None))

		if snip:
			erange = self.cb.eng(self.channels) if xcalib else self.channels
			ax.plot(erange, self._snip)

		lbs = []
		if fit:
			cm, cl = colors(), colors(aslist=True)
			cl = [c for c in cl if c not in [cm['k'], cm['gy']]]
			ls = ['-','-.','--']
			sub, itp, N_sub = self._split_fits()
			
			for n,p in enumerate(self.fits):
				if not p['converged']:
					xgrid = np.arange(p['l'], p['h'], 0.1)
					pk_fit = self.multiplet(xgrid, *p['fit'])
					erange = self.cb.eng(xgrid) if xcalib else xgrid
					ax.plot(erange, np.where(pk_fit>0.1, pk_fit, 0.1), ls=':', lw=1.2, color=cm['gy'])
				elif N_sub[n]>1:
					xgrid = np.arange(p['l'], p['h'], 0.1)
					pk_fit = self.multiplet(xgrid, *p['fit'])
					erange = self.cb.eng(xgrid) if xcalib else xgrid
					ax.plot(erange, np.where(pk_fit>0.1, pk_fit, 0.1), ls='--', lw=1.4, color=cm['gy'])

			
			for n,i in enumerate(itp):
				c = cl[n%len(cl)]
				ilbl = Isotope(i).TeX
				for p,lh in sub[i]:
					xgrid = np.arange(lh[0], lh[1], 0.1)
					pk_fit = self.multiplet(xgrid, *p)
					lb = ilbl if (labels and ilbl not in lbs and len(lbs)<20) else None
					erange = self.cb.eng(xgrid) if xcalib else xgrid
					ax.plot(erange, np.where(pk_fit>0.1, pk_fit, 0.1), lw=1.4, color=c, label=lb, ls=ls[int(n/len(cl))%len(ls)])
					if lb is not None:
						lbs.append(lb)
		
		if xcalib:
			ax.set_xlabel('Energy (keV)')
		else:
			ax.set_xlabel('ADC Channel')
		ax.set_ylabel('Counts')
		if labels:
			ncol = min([max([int(len(lbs)/5),1]), 3])
			ax.legend(loc=0, ncol=ncol)

		return _close_plot(f, ax, **kwargs)
