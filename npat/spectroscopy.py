from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sqlite3, os
import numpy as np
import matplotlib.pyplot as plt
import datetime as dtm

from matplotlib import gridspec
from scipy.optimize import curve_fit
from scipy.special import erfc
from scipy.interpolate import interp1d

from .isotope import Isotope
from .plotter import colors

DB_CONNECTIONS_DICT = {}

class Calibration(object):
	def __init__(self, meta={}):
		self.calib = {}
		self._default = {'engcal': [0.0, 0.3],
						'effcal': [0.331, 0.158, 0.410, 0.001, 1.476],
						'unc_effcal': [[ 5.038e-02,  3.266e-02, -2.151e-02, -4.869e-05, -7.748e-03],
									 [ 3.266e-02,  2.144e-02, -1.416e-02, -3.416e-05, -4.137e-03],
									 [-2.151e-02, -1.416e-02,  9.367e-03,  2.294e-05,  2.569e-03],
									 [-4.869e-05, -3.416e-05,  2.294e-05,  5.411e-07, -1.165e-04],
									 [-7.748e-03, -4.137e-03,  2.569e-03, -1.165e-04,  3.332e-02]],
						'rescal': [0.05]}
		self.update(meta)

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

	@property
	def unc_effcal(self):
		return self.calib['unc_effcal']

	@property
	def rescal(self):
		return self.calib['rescal']
	
	def eng(self, idx, *cal):
		cal = cal if len(cal) else self.calib['engcal']
		if len(cal)==2:
			return cal[0]+cal[1]*idx
		elif len(cal)==3:
			return cal[0]+cal[1]*idx+cal[2]*idx**2

	def eff(self, energy, *c):
		c = c if len(c) else self.calib['effcal']
		if len(c)==3:
			return c[0]*np.exp(-c[1]*energy**c[2])
		elif len(c)==5:
			return c[0]*np.exp(-c[1]*energy**c[2])*(1.0-np.exp(-c[3]*energy**c[4]))
			
	def unc_eff(self, energy, c=None, u=None):
		if c is None or u is None:
			c, u = self.calib['effcal'], self.calib['unc_effcal']
		var, eps = np.zeros(len(energy)), 1E-8
		for n in range(len(c)):
			for m in range(n, len(c)):
				c_n, c_m = list(c), list(c)
				c_n[n], c_m[m] = c_n[n]+eps, c_m[m]+eps
				par_n = (self.eff(energy, *c_n)-self.eff(energy, *c))/eps
				par_m = (self.eff(energy, *c_m)-self.eff(energy, *c))/eps
				var += u[n][m]*par_n*par_m*(2.0 if n!=m else 1.0)
		return np.sqrt(var)

	def res(self, idx, *cal):
		cal = cal if len(cal) else self.calib['rescal']
		if len(cal)==1:
			return cal[0]*np.sqrt(idx)

	def map_idx(self, energy, *cal):
		cal = cal if len(cal) else self.calib['engcal']
		if len(cal)==2:
			return np.array(np.rint((energy-cal[0])/cal[1]), dtype=np.int32)
		elif len(cal)==3:
			return np.array(np.rint(0.5*(np.sqrt(cal[1]**2-4.0*cal[2]*(cal[0]-energy))-cal[1])/cal[2]), dtype=np.int32)

	def calibrate(self, spectra, saveas=None, db=None, auto_calibrate=False):
		shelves = []
		sps = []
		self._calib_data = {'effcal':[]}
		for sp in spectra:
			if type(sp)==str:
				sp = Spectrum(sp, db)
			sps.append(sp)
			if 'shelf' in sp.meta:
				shelves.append(sp.meta['shelf'])
			else:
				sp.meta = {'shelf':None}
				shelves.append(None)
		cal_spec = [sp for sp in spectra if ('A0' in sp.meta and 'ref_date' in sp.meta)]
		if len(cal_spec)==0:
			raise ValueError('Cannot calibrate: A0 and ref_date required for calibration.')

		if auto_calibrate:
			for sp in cal_spec:
				sp.auto_calibrate(**{'_cb_cal':False})
		p0_engcal = np.average([s.cb.engcal for s in cal_spec if len(s.cb.engcal)==len(cal_spec[0].cb.engcal)], axis=0)
		engcal = self._calibrate_energy(cal_spec, p0_engcal)
		rescal = self._calibrate_resolution(cal_spec, np.average([s.cb.rescal for s in cal_spec], axis=0))

		shelves = {sh:[s for s in cal_spec if s.meta['shelf']==sh] for sh in set(shelves)}
		for shelf in shelves:
			spectra = shelves[shelf]
			p0_effcal = np.average([s.cb.effcal for s in spectra if len(s.cb.effcal)==len(spectra[0].cb.effcal)], axis=0)
			effcal, unc_effcal = self._calibrate_efficiency(spectra, p0_effcal)
			for sp in sps:
				if shelf==sp.meta['shelf']:
					sp.meta = {'engcal':engcal, 'effcal':effcal, 'unc_effcal':unc_effcal, 'rescal':rescal}
					if saveas is not None:
						if type(saveas)==str:
							saveas = [saveas]
						sp.save(saveas)

	def _calibrate_energy(self, spectra, p0):
		E, chn, unc_chn = [], [], []
		for sp in spectra:
			sp.meta = {'engcal':p0}
			for f in sp.fits:
				E += f.gm[:,0].tolist()
				chn += f.fit_param('mu').tolist()
				unc_chn += np.sqrt(f.cov_param('mu')).tolist()

		E, chn, unc_chn = np.array(E), np.array(chn), np.array(unc_chn)
		idx = np.where(chn>unc_chn)
		E, chn, unc_chn = E[idx], chn[idx], unc_chn[idx]

		p0 = p0 if len(p0)==3 else [p0[0], p0[1], 0.0]
		fit, unc = curve_fit(self.map_idx, E, chn, sigma=unc_chn, p0=p0)

		self._calib_data['engcal'] = {'fit':fit, 'unc':unc, 'x':E, 'y':chn, 'yerr':unc_chn}
		return fit

	def _calibrate_efficiency(self, spectra, p0):
		### NEED TO UPDATE WITH MORE PRECISE EFFICIENCY ###
		E, eff, unc_eff = [], [], []
		for sp in spectra:
			sp.meta = {'effcal':p0}
			for f in sp.fits:
				E += f.gm[:,0].tolist()
				if type(sp.meta['A0'])==list:
					A = np.array([sp.meta['A0'][sp.meta['istp'].index(i)] for i in f.istp])
				else:
					A = sp.meta['A0']
				eff += (f.activity[0]*f.cb.eff(f.gm[:,0])/A).tolist()
				unc_eff += (f.activity[1]*f.cb.eff(f.gm[:,0])/A).tolist()

		E, eff, unc_eff = np.array(E), np.array(eff), np.array(unc_eff)
		idx = np.where(eff>unc_eff)
		E, eff, unc_eff = E[idx], eff[idx], unc_eff[idx]
		p0 = p0 if len(p0)==5 else p0.tolist()+[0.001, 1.476]
		p0[0] = p0[0]*np.average(np.array(eff)/self.eff(np.array(E), *p0), weights=(self.eff(np.array(E), *p0)/np.array(unc_eff))**2)

		try:		
			bounds = ([0.0, 0.0, -1.0, 0.0, -2.0], [10.0, 2.0, 3.0, 0.5, 3.0])
			fit, unc = curve_fit(self.eff, E, eff, sigma=unc_eff, p0=p0, bounds=bounds)
		except:
			fit, unc = curve_fit(self.eff, E, eff, sigma=unc_eff, p0=p0[:3], bounds=(bounds[0][:3], bounds[1][:3]))

		self._calib_data['effcal'].append({'fit':fit, 'unc':unc, 'x':E, 'y':eff, 'yerr':unc_eff, 'shelf':spectra[0].meta['shelf']})
		return fit, unc

	def _calibrate_resolution(self, spectra, p0):
		mu, sig, unc_sig = [], [], []
		for sp in spectra:
			sp.meta = {'rescal':p0}
			for f in sp.fits:
				mu += f.fit_param('mu').tolist()
				sig += f.fit_param('sig').tolist()
				unc_sig += np.sqrt(f.cov_param('sig')+0.25*p0[0]*f.cov_param('mu')/f.fit_param('mu')).tolist()
		mu, sig, unc_sig = np.array(mu), np.array(sig), np.array(unc_sig)
		idx = np.where(sig>unc_sig)
		mu, sig, unc_sig = mu[idx], sig[idx], unc_sig[idx]
		fit, unc = curve_fit(self.res, mu, sig, sigma=unc_sig, p0=p0)

		self._calib_data['rescal'] = {'fit':fit, 'unc':unc, 'x':mu, 'y':sig, 'yerr':unc_sig}
		return fit

	def plot(self, show=True, saveas=None):
		cm = colors()
		f = plt.figure(figsize=(10, 6))
		gs = gridspec.GridSpec(2, 10)
		ax0 = plt.subplot(gs[0, 0:4])
		ax1 = plt.subplot(gs[1, 0:4])
		ax2 = plt.subplot(gs[0:2, 4:10])

		d = self._calib_data['engcal']
		ax0.errorbar(d['x'], d['y'], yerr=d['yerr'], ls='None', marker='o')
		x = np.arange(min(d['x']), max(d['x']), 0.1)
		ax0.plot(x, self.map_idx(x, *d['fit']))
		ax0.set_xlabel('Energy (keV)')
		ax0.set_ylabel('ADC Channel')

		d = self._calib_data['rescal']
		ax1.errorbar(d['x'], d['y'], yerr=d['yerr'], ls='None', marker='o')
		x = np.arange(min(d['x']), max(d['x']), 0.1)
		ax1.plot(x, self.res(x, *d['fit']))
		ax1.set_xlabel('ADC Channel')
		ax1.set_ylabel('Peak Width')
		
		for d in self._calib_data['effcal']:
			ax2.errorbar(d['x'], d['y'], yerr=d['yerr'], ls='None', marker='o')
			x = np.arange(min(d['x']), max(d['x']), 0.1)
			ax2.plot(x, self.eff(x, *d['fit']), color=cm['k'], label=d['shelf'])
			low = self.eff(x, *d['fit'])-self.unc_eff(x, d['fit'], d['unc'])
			high = self.eff(x, *d['fit'])+self.unc_eff(x, d['fit'], d['unc'])
			ax2.fill_between(x, low, high, facecolor=cm['gy'], alpha=0.5)
		ax2.set_xlabel('Energy (keV)')
		ax2.set_ylabel('Efficiency')
		if len([d['shelf'] for d in self._calib_data['effcal']]):
			ax2.legend(loc=0)

		f.tight_layout()
		if saveas is not None:
			f.savefig(saveas)
		if show:
			plt.show()
		plt.close()


	


class PeakFit(object):
	def __init__(self,*args):
		if len(args):
			self.fit_params(*args)
	def fit_params(self, x, hist, gm, istp, A0, cov_A0, meta, cb, snip):
		self.x = x
		self.hist = hist
		self.gm = np.asarray(gm)
		self.istp = istp
		self.A0 = np.array(A0)
		self.cov_A0 = np.array(cov_A0)
		D_inv = 1.0/np.sqrt(np.diag(cov_A0))
		self.cor_A0 = cov_A0*D_inv.reshape((len(D_inv),1))*D_inv
		self.fit_config = meta['fit_config']
		self.meta = meta
		self.cb = cb
		self.snip = interp1d(x, snip)
		self._snip = snip
		self._fit, self._cov = None, None
		self._chi2, self._counts, self._activity = None, None, None
		self._fit_param, self._cov_param = None, None

	def p0_bg(self, quadratic=False):
		N = 3 if quadratic else 2
		M = np.column_stack([self.x**m for m in range(N)])
		b = np.array([np.sum(self.x**m*self._snip) for m in range(N)])
		M_inv = np.linalg.inv(np.dot(M.T, M))
		p = np.dot(M_inv, b)
		resid = self._snip-np.dot(M, p)
		chi2 = np.dot(resid.T, resid)/float(len(self.x)-N)
		return p.tolist(), M_inv*chi2

	def peak(self, x, A, mu, sig, R, alpha, step):
		r2 = 1.41421356237
		return A*np.exp(-0.5*((x-mu)/sig)**2)+R*A*np.exp((x-mu)/(alpha*sig))*erfc((x-mu)/(r2*sig)+1.0/(r2*alpha))+step*A*erfc((x-mu)/(r2*sig))

	def Npeak(self, x, *args):
		if self.fit_config['bg_fit']:
			if self.fit_config['quad_bg']:
				L, peak = 3, args[0]+args[1]*x+args[2]*x**2
			else:
				L, peak = 2, args[0]+args[1]*x
		elif len(x)==len(self._snip):
			L, peak = 0, self._snip.copy()
		else:
			L, peak = 0, self.snip(x)
		R, alpha, step = self.fit_config['R'], self.fit_config['alpha'], self.fit_config['step']
		if self.fit_config['skew_fit']:
			if self.fit_config['step_fit']:
				for n in range(int((len(args)-L)/6)):
					peak += self.peak(x,*args[6*n+L:6*n+6+L])
			else:
				for n in range(int((len(args)-L)/5)):
					peak += self.peak(x,*(args[5*n+L:5*n+5+L]+(step,)))
		else:
			if self.fit_config['step_fit']:
				for n in range(int((len(args)-L)/4)):
					peak += self.peak(x,*(args[4*n+L:4*n+L+3]+(R, alpha)+(args[4*n+3+L])))
			else:
				for n in range(int((len(args)-L)/3)):
					peak += self.peak(x,*(args[3*n+L:3*n+L+3]+(R,alpha,step)))
		return peak

	@property	
	def fit(self):
		if self._fit is not None:
			return self._fit
		B_f = self.fit_config['fit_bound'] if self.fit_config['fit_bound'] is not None else 1.0
		bounds = [[],[]]
		if self.fit_config['bg_fit']:
			p0, cov_p0 = self.p0_bg(self.fit_config['quad_bg'])
			bounds[0] += (p0-B_f*np.sqrt(np.diag(cov_p0))).tolist()
			bounds[0][0] -= 10.0*B_f*np.sqrt(cov_p0[0][0])
			bounds[1] += (p0+B_f*np.sqrt(np.diag(cov_p0))).tolist()
			bounds[1][1] += 10.0*B_f*np.sqrt(cov_p0[0][0])
		else:
			p0 = []
		R, alpha, step = self.fit_config['R'], self.fit_config['alpha'], self.fit_config['step']
		for n,g in enumerate(self.gm):
			idx = float(self.cb.map_idx(g[0]))
			sig = float(self.cb.res(idx))
			p0 += [self.A0[n]*self.cb.eff(g[0])*g[1]/sig, idx, sig]
			bounds[0] += [0.0, idx-B_f*idx/sig, sig-sig/B_f] ### bounds[0][0] = 0.1*B_f*self.A0[n]*self.cb.eff(g[0])*g[1]/sig
			bounds[1] += [p0[-3]+10.0*B_f*self.A0[n]*self.cb.eff(g[0])*g[1]/sig, idx+B_f*idx/sig, sig+0.5*sig*B_f]
			if self.fit_config['skew_fit']:
				p0 += [R, alpha]
				bounds[0] += [0.0, 0.5]
				bounds[1] += [1.0, max((2.5, alpha))]
			if self.fit_config['step_fit']:
				p0.append(step)
				bounds[0].append(0.0)
				bounds[1].append(0.1)
		if self.fit_config['fit_bound'] is None:
			bounds = ([-np.inf for i in p0],[np.inf for i in p0])
		try:
			self._fit, self._cov = curve_fit(self.Npeak, self.x, self.hist, sigma=np.sqrt(self.hist+0.1), p0=p0, bounds=tuple(bounds))
		except Exception, err:
			print('Error on peak at {} (keV)'.format(round(self.gm[0][0],1)))
			print(Exception, err)
			self._fit, self._cov = np.array(p0), None
		return self._fit

	@property
	def cov(self):
		if self._fit is None:
			fit = self.fit
		return self._cov

	def fit_param(self, param='mu'):
		if self._fit_param is None:
			self._fit_param = {'A':[], 'mu':[], 'sig':[]}
			cfg = self.fit_config
			if cfg['skew_fit']:
				self._fit_param['R'] = []
				self._fit_param['alpha'] = []
			if cfg['step_fit']:
				self._fit_param['step'] = []
			f = self.fit
			if cfg['bg_fit']:
				self._fit_param['a'] = f[0]
				self._fit_param['b'] = f[1]
				if cfg['quad_bg']:
					self._fit_param['c'] = f[2]
			L = 3+2*int(cfg['skew_fit'])+int(cfg['step_fit'])
			for m in range(int(len(f)/L)):
				i = L*m+2*int(cfg['bg_fit'])+int(cfg['bg_fit'])*int(cfg['quad_bg'])
				self._fit_param['A'].append(f[i])
				self._fit_param['mu'].append(f[i+1])
				self._fit_param['sig'].append(f[i+2])
				if cfg['skew_fit']:
					self._fit_param['R'].append(f[i+3])
					self._fit_param['alpha'].append(f[i+4])
				if cfg['step_fit']:
					self._fit_param['step'].append(f[i+5])
			self._fit_param = {p:np.array(self._fit_param[p]) for p in self._fit_param}
		return self._fit_param[param]

	def cov_param(self, param='mu'):
		if self._cov_param is None:
			if self.cov is None:
				return np.inf*np.ones(len(self.gm))
			self._cov_param = {'A':[], 'mu':[], 'sig':[]}
			cfg = self.fit_config
			if cfg['skew_fit']:
				self._cov_param['R'] = []
				self._cov_param['alpha'] = []
			if cfg['step_fit']:
				self._cov_param['step'] = []
			f = self.cov
			if cfg['bg_fit']:
				self._cov_param['a'] = f[0][0]
				self._cov_param['b'] = f[1][1]
				if cfg['quad_bg']:
					self._cov_param['c'] = f[2][2]
			L = 3+2*int(cfg['skew_fit'])+int(cfg['step_fit'])
			for m in range(int(len(f)/L)):
				i = L*m+2*int(cfg['bg_fit'])+int(cfg['bg_fit'])*int(cfg['quad_bg'])
				self._cov_param['A'].append(f[i][i])
				self._cov_param['mu'].append(f[i+1][i+1])
				self._cov_param['sig'].append(f[i+2][i+2])
				if cfg['skew_fit']:
					self._cov_param['R'].append(f[i+3][i+3])
					self._cov_param['alpha'].append(f[i+4][i+4])
				if cfg['step_fit']:
					self._cov_param['step'].append(f[i+5][i+5])
			self._cov_param = {p:np.array(self._cov_param[p]) for p in self._cov_param}
		return self._cov_param[param]


	@property
	def chi2(self):
		if self._chi2 is not None:
			return self._chi2
		non_zero = np.where(self.hist>0)
		dof = float(len(non_zero[0])-len(self.fit)-1)
		if dof==0:
			self._chi2 = np.inf
		self._chi2 = np.sum((self.hist[non_zero]-self.Npeak(self.x[non_zero], *self.fit))**2/self.hist[non_zero])/dof
		return self._chi2


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
	
	@property
	def counts(self):
		if self._counts is not None:
			return self._counts
		cfg = self.fit_config
		L = 3+2*int(cfg['skew_fit'])+int(cfg['step_fit'])
		f, u = self.fit, self._cov
		N_cts, unc_N = [], []

		skew_fn = lambda A, R, alpha, sig: 2*A*R*alpha*sig*np.exp(-0.5/alpha**2)
		min_skew_fn = lambda A, sig: 2*A*cfg['R']*cfg['alpha']*sig*np.exp(-0.5/cfg['alpha']**2)
		pk_fn = lambda A, sig: 2.506628*A*sig

		for m in range(int(len(self.fit)/L)):
			i = L*m+2*int(cfg['bg_fit'])+int(cfg['bg_fit'])*int(cfg['quad_bg'])
			p_i = [i, i+2]
			s_i = [i, i+3, i+4, i+2]
			if cfg['skew_fit']:
				N_cts.append(pk_fn(*f[p_i])+skew_fn(*f[s_i]))
			else:
				N_cts.append(pk_fn(*f[p_i])+min_skew_fn(*f[p_i]))
			if u is not None:
				if not np.isinf(u[i][i]):
					if cfg['skew_fit']:
						skew_unc = self._unc_calc(skew_fn, f[s_i], u[np.ix_(s_i, s_i)])
					else:
						skew_unc = self._unc_calc(min_skew_fn, f[p_i], u[np.ix_(p_i, p_i)])
					unc_N.append(np.sqrt(self._unc_calc(pk_fn, f[p_i], u[np.ix_(p_i, p_i)])+skew_unc))
				else:
					unc_N.append(np.inf)
			else:
				unc_N.append(np.inf)
		self._counts = (np.array(N_cts), np.array(unc_N))
		return self._counts
	
	@property
	def activity(self):
		if 'live_time' not in self.meta or 'real_time' not in self.meta:
			print('WARNING: Need real-time/live-time to compute activity.')
			return np.zeros(len(gm)), np.zeros(len(gm))
		if self._activity is not None:
			return self._activity
		N, unc_N = self.counts
		A = N*(self.meta['real_time']/self.meta['live_time']**2)/(self.gm[:,1]*self.cb.eff(self.gm[:,0]))
		unc_A = np.sqrt((unc_N/N)**2+(self.cb.unc_eff(self.gm[:,0])/self.cb.eff(self.gm[:,0]))**2+(self.gm[:,2]/self.gm[:,1])**2)*A
		self._activity = (A, unc_A)
		return self._activity
	

	def plot(self, show=True, saveas=None, logscale=False, grayscale=False, labels=False):
		f, ax = plt.subplots()

		erange = self.cb.eng(np.array([self.x-0.5, self.x+0.5])).T.flatten()
		spec = np.array([self.hist, self.hist]).T.flatten()

		ax.plot(erange, spec, lw=1.2)
		xgrid = np.arange(self.x[0], self.x[-1], 0.1)
		if self.fit is not None:
			ax.plot(self.cb.eng(xgrid), self.Npeak(xgrid, *self.fit), lw=1.2)

		if logscale:
			ax.set_yscale('log')

		ax.set_xlabel('Energy (keV)')
		ax.set_ylabel('Counts')
		# ax.legend(loc=0)
		f.tight_layout()

		if saveas is not None:
			f.savefig(saveas)
		if show:
			plt.show()
		plt.close()




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
				if self._fnm.endswith('.Spe'):
					self._from_Spe()
			else:
				raise ValueError('File does not exist: {}'.format(filename))

			self._check_db(db)
		self.cb.update(self._meta)

	def _default_params(self):
		self._path, self._fnm = None, None
		self._hist = np.zeros(2, dtype=np.int32)
		self._meta = {'fit_config':{'snip_sig':8.5, 'snip_adj':1.5, 'snip_alpha':0.15,
									'R':0.1, 'alpha':0.9, 'step':0.0, 'bg_fit':False,
									'skew_fit':False, 'step_fit':False, 'quad_bg':False,
									'iter_fit':True, 'SNR_cut':3.0, 'fit_bound':2.5, 
									'xrays':False, 'pk_width':7.5, 'E_min':75.0,
									'I_min':0.05}, 'istp':[], 'shelf':None}
		self.cb = Calibration()
		for nm in self.cb.calib:
			self._meta[nm] = self.cb.calib[nm]
		self._gamma_list = None
		self._fits = None

	def _get_db_params(self):
		pass

	def _check_db(self, db=None):
		if db is not None:
			path, fnm = os.path.split(db)
			if path in ['',' ']:
				path = os.getcwd()
			if fnm in ['',' ']:
				raise ValueError('Invalid db Filename: {}'.format(db))
			self.db_fnm = os.path.join(path, fnm)
		else:
			self.db_fnm = os.path.join(self._path, 'npat_analysis.db')
		if os.path.exists(self.db_fnm):
			global DB_CONNECTIONS_DICT
			if self.db_fnm in DB_CONNECTIONS_DICT:
				self.db_connection = DB_CONNECTIONS_DICT[self.db_fnm]
				self.db = self.db_connection.cursor()
			else:
				self.db_connection = sqlite3.connect(self.db_fnm)
				DB_CONNECTIONS_DICT[self.db_fnm] = self.db_connection
				self.db = self.db_connection.cursor()
			self._get_db_params()
		else:
			self.db_connection = None
			self.db = None

	@property
	def hist(self):
		return self._hist

	@hist.setter
	def hist(self, hist_array):
		self._hist = np.asarray(hist_array, dtype=np.int32)
		self._snip = self.snip_bg()

	def _from_Spe(self, filename=None):
		if filename is None:
			filename = os.path.join(self._path, self._fnm)
		with open(filename) as f:
			fl = map(str.strip, f.readlines())

		idx = fl.index('$DATA:')
		head = fl[:idx]
		self.hist = fl[idx+2:idx+3+int(fl[idx+1].split(' ')[1])]
		foot = fl[idx+2+len(self._hist):]

		start_time = dtm.datetime.strptime(head[head.index('$DATE_MEA:')+1], '%m/%d/%Y %H:%M:%S')
		LT, RT = tuple(map(float, head[head.index('$MEAS_TIM:')+1].split(' ')))
		engcal = map(float, foot[foot.index('$MCA_CAL:')+2].split(' ')[:-1])

		self.meta = {'start_time':start_time, 'live_time':LT, 'real_time':RT, 'engcal':engcal,
						'Spe_head':'\n'.join(head), 'Spe_foot':'\n'.join(foot)}

	@property
	def meta(self):
		return self._meta

	@meta.setter
	def meta(self, meta_dict):
		for nm in meta_dict:
			if nm.endswith('cal'):
				self.cb.calib[nm] = meta_dict[nm]
			self._meta[nm] = meta_dict[nm]
	
	@property
	def fit_config(self):
		return self._meta['fit_config']
	
	@fit_config.setter
	def fit_config(self, config_dict):
		for nm in config_dict:
			self._meta['fit_config'][nm] = config_dict[nm]

	def exp_smooth(self, x, alpha=0.3):
		N = int(2.0/alpha)-1
		wts = np.array([alpha*(1.0-alpha)**abs(N-i) for i in range(2*N+1)])
		y = np.concatenate((np.ones(N)*x[0], x, np.ones(N)*x[-1]))
		return np.dot(wts/np.sum(wts), np.take(y, [np.arange(i,len(x)+i) for i in range(2*N+1)]))

	def snip_bg(self):
		sig, adj, alpha = self.fit_config['snip_sig'], self.fit_config['snip_adj'], self.fit_config['snip_alpha']
		r, x, dead = self.cb.rescal[0], np.arange(len(self._hist)), int(self.cb.rescal[0]*sig*len(self._hist)**0.5)
		V_i, L = np.log(np.log(np.sqrt(self._hist+1.0)+1.0)+1.0), len(x)-dead
		while self._hist[dead]==0:
			dead += 1
		
		for M in np.linspace(0, sig, 10):
			l, h = np.array(x-M*r*np.sqrt(x), dtype=np.int32), np.array(x+M*r*np.sqrt(x), dtype=np.int32)
			V_i[dead:L] = np.minimum(V_i[dead:L],0.5*(V_i[l[dead:L]]+V_i[h[dead:L]]))

		snip = self.exp_smooth((np.exp(np.exp(V_i)-1.0)-1.0)**2-1.0, alpha)
		return snip+adj*np.sqrt(snip+0.1)

	@property
	def gamma_list(self):
		if self._gamma_list is None:
			xrays, E_min, I_min = self.fit_config['xrays'], self.fit_config['E_min'], self.fit_config['I_min']
			gammas = [Isotope(i).gammas(E_lim=[E_min, self.cb.eng(0.99*len(self._hist))], I_lim=[I_min, None], xrays=xrays) for i in self._meta['istp']]
			no511 = [np.where(np.absolute(np.array(i['E'])-511.0)>5.0) for i in gammas]
			self._gamma_list = [(np.array([i['E'], i['I'], i['dI']]).T)[no511[n]]*np.array([1.0, 0.01, 0.01]) for n,i in enumerate(gammas)]
		return self._gamma_list

	def _guess_p0_A(self, *engcal):
		x, Y, p0_A, pk = np.arange(len(self._hist)), [], [], PeakFit().peak
		R, alpha, step = self.fit_config['R'], self.fit_config['alpha'], self.fit_config['step']
		clip = np.where(self._hist>self._snip, self._hist-self._snip, 0.0)
		W = self.fit_config['pk_width']
		for gm in self.gamma_list:
			y = np.zeros(len(self._hist))
			if len(gm):
				idx = self.cb.map_idx(gm[:,0], *engcal)
				sig = self.cb.res(idx)
				norm = self.cb.eff(gm[:,0])*gm[:,1]/sig
				p0_A.append((np.average((clip[idx]/norm), weights=np.sqrt(norm))))
				if p0_A[-1]<0:
					p0_A[-1] = 0.0
				l, h = np.array(idx-W*sig, dtype=np.int32), np.array(idx+W*sig, dtype=np.int32)
				for m,i in enumerate(idx):
					if l[m]>0 and h[m]<len(x):
						y[l[m]:h[m]] += pk(x[l[m]:h[m]], norm[m], i, sig[m], R, alpha, step)
				Y.append(y)
		return p0_A, np.array(Y)

	def _fit_p0_A(self):
		p0_A, Y = self._guess_p0_A()
		try:
			return curve_fit(lambda x, *A: self._snip+np.dot(A, Y), np.arange(len(self._hist)), self._hist, p0=p0_A)
		except:
			return p0_A, np.ones((len(p0_A), len(p0_A)))*np.inf

	def auto_calibrate(self, *guess, **param):
		from scipy.optimize import differential_evolution
		guess = guess if len(guess) else self.cb.engcal
		obj = lambda m: np.average(np.sqrt(np.absolute(self._hist-self._snip-np.dot(*self._guess_p0_A(*[guess[0], m[0]])))))
		self.meta = {'engcal': [guess[0], differential_evolution(obj, [(0.975*guess[1], 1.025*guess[1])]).x[0]]}
		if ('A0' in self.meta and 'ref_date' in self.meta) and (param['_cb_cal'] if '_cb_cal' in param else True):
			try:
				self.cb.calibrate([self])
			except:
				self.meta = {'engcal':guess}
				try:
					self.cb.calibrate([self])
				except:
					self.meta = {'engcal':guess}
					print('WARNING: auto calibrate failed. Setting engcal to guess.')



	def _get_p0(self):
		p0_A, cov_p0_A = self._fit_p0_A()
		pks = []
		W, o = self.fit_config['pk_width'], 0
		for n,gm in enumerate(self.gamma_list):
			if len(gm):
				idx = self.cb.map_idx(gm[:,0])
				sig = self.cb.res(idx)
				SNR = (p0_A[n-o]*self.cb.eff(gm[:,0])*gm[:,1]/sig)/np.sqrt(self._snip[idx])
				cut = np.where((SNR>self.fit_config['SNR_cut'])&(idx>W*sig)&(idx<(len(self._hist)-W*sig)))
				idx, sig, gm_idx = idx[cut], sig[cut], np.arange(len(gm))[cut]
				l, h = np.array(idx-W*sig, dtype=np.int32), np.array(idx+W*sig, dtype=np.int32)
				pks.append(np.array([l, h, np.repeat(n, len(l)), gm_idx, np.repeat(o, len(l))], dtype=np.int32).T)
			else:
				o += 1
		if len(pks)==0:
			return pks
		pks = np.concatenate(pks)
		pks = pks[np.argsort(pks[:,0])]
		pairs = np.where(pks[:-1,1]>pks[1:,0])[0]
		pairs = np.unique(np.concatenate((pairs, pairs+1)))
		if len(pairs):
			groups = [[pairs[0],pairs[1]]]
			for p in pairs[2:]:
				if (p-groups[-1][-1])==1 and len(groups[-1])<8:
					groups[-1].append(p)
				else:
					groups.append([p])
		pk_pairs = [[pks[i] for i in p] for p in groups] if len(pairs) else []
		x = np.arange(len(self._hist))
		p0 = [PeakFit(x[g[0][0]:g[-1][1]], self._hist[g[0][0]:g[-1][1]], [self.gamma_list[p[2]][p[3]] for p in g], [self._meta['istp'][p[2]] for p in g],
				[p0_A[p[2]-p[4]] for p in g], [[cov_p0_A[i[2]-i[4],p[2]-p[4]] for i in g] for p in g], self._meta, self.cb, self._snip[g[0][0]:g[-1][1]]) for g in pk_pairs]
		pks = np.delete(pks, pairs, axis=0)
		p0 += [PeakFit(x[p[0]:p[1]], self._hist[p[0]:p[1]], [self.gamma_list[p[2]][p[3]]], [self._meta['istp'][p[2]]],
				[p0_A[p[2]-p[4]]], [[cov_p0_A[p[2]-p[4],p[2]-p[4]]]], self._meta, self.cb, self._snip[p[0]:p[1]]) for p in pks]
		return p0

	@property
	def fits(self):
		"""List of peak fit objects
		"""
		if self._fits is None:
			self._fits = self._get_p0()
		return self._fits

	def save(self, *saveas):
		for fl in saveas:
			print(fl)

	def summarize(self, printout=True, saveas=None):
		if printout:
			for f in self.fits:
				for n in range(len(f.gm)):
					ln1 = [f.gm[n][0], f.istp[n], f.gm[n][1]*100.0]
					print('{0} keV peak in {1} ({2}% intensity)'.format(*ln1))
					ln2 = [int(f.counts[0][n]), (int(f.counts[1][n]) if np.isfinite(f.counts[1][n]) else np.inf),
							round(f.activity[0][n],1),round(f.activity[1][n],1), 
							round(f.activity[0][n]/3.7E4,3), round(f.activity[1][n]/3.7E4,3), 
							round(f.chi2,3)]
					print('N:{0}({1}), A(Bq):{2}({3}), A(uCi):{4}({5}), chi2/dof:{6}'.format(*ln2))
					print('')

	def plot(self, show=True, fit=True, saveas=None, zoom=None, logscale=True, 
				grayscale=False, labels=False, square_fig=False):
		
		erange = self.cb.eng(np.array([np.arange(len(self._hist))-0.5, np.arange(len(self._hist))+0.5])).T.flatten()
		spec = np.array([self.hist, self.hist]).T.flatten()

		f, ax = plt.subplots(figsize=(None if square_fig else (12.8, 4.8)))
		spec_label = self._fnm if self._fnm is not None else 'Spectrum'
		ax.plot(erange, spec, lw=1.2, zorder=1, label=spec_label)

		if fit:
			for p in self.fits:
				xgrid = np.arange(p.x[0], p.x[-1], 0.1)
				ax.plot(self.cb.eng(xgrid), p.Npeak(xgrid, *p.fit), lw=1.2)

		if logscale:
			ax.set_yscale('log')
			ax.set_ylim((0.5, ax.get_ylim()[1]))

		ax.set_xlabel('Energy (keV)')
		ax.set_ylabel('Counts')
		ax.legend(loc=0)
		f.tight_layout()

		if saveas is not None:
			f.savefig(saveas)
		if show:
			plt.show()
		plt.close()

