from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import re
import os
import numpy as np
import sqlite3
import pandas as pd

from .plotter import colors
from .plotter import _init_plot
from .plotter import _close_plot
from .dbmgr import get_cursor
from .isotope import Isotope

global ZG_CONNECTIONS_DICT
ZG_CONNECTIONS_DICT = {}

#### TODO: allow new compound definitions as {'El':wt_pct} and accept .json stack

class Ziegler(object):
	"""Method for solving energy loss within stacked-target foil experiment

	...

	Parameters
	----------
	x : type
		Description of parameter `x`.

	Attributes
	----------

	Methods
	-------

	Notes
	-----

	References
	----------

	Examples
	--------

	"""

	def __init__(self, stack, beam=None, **kwargs):
		### stack is list of dicts, which must have 'compound' specified.
		### Areal density specified by either 'ad' (mg/cm^2), both mass (g) and area (cm^2),
		### 'thickness' (mm) or both 'density' (g/cm^3) and 'thickness' (mm)

		zdb = get_cursor('ziegler')
		
		self.protons = {int(i[0]):list(map(float,i[1:])) for i in zdb.execute('SELECT * FROM protons')}
		self.helium = {int(i[0]):list(map(float,i[1:])) for i in zdb.execute('SELECT * FROM helium')}
		self.ionization = {int(i[0]):list(map(float,i[1:])) for i in zdb.execute('SELECT * FROM ionization')}
		self.weights = {int(i[0]):list(map(float,i[1:3])) for i in zdb.execute('SELECT * FROM weights')}
		self.compounds = {str(i[0]):[[int(h.split(':')[0]), float(h.split(':')[1])] for h in i[2].split(',')] for i in zdb.execute('SELECT * FROM compounds')}
		self.compounds = {cm:[[i[0],i[1]/sum([m[1] for m in self.compounds[cm]])] for i in self.compounds[cm]] for cm in self.compounds}
		self.densities = {str(i[0]):float(i[1]) for i in zdb.execute('SELECT * FROM compounds')}
		self.elements = sorted([[str(i[0]), int(i[2].split(':')[0])] for i in zdb.execute('SELECT * FROM compounds') if len(i[2].split(':'))==2], key=lambda h:len(h[0]), reverse=True)

		self._meta = {}
		self.meta = {'beam_istp':'1H', 'E0':33.0, 'dE0':0.3, 'N':10000, 'dp':1.0,
						'chunk_size':1E7, 'threads':1,
						'solved':False, 'accuracy':0.01, 'min_steps':2, 'max_steps':50}
		self._stack = []
		if beam is not None:
			print('Keyword `beam` deprecated: see documentation for proper input.')
			self.meta = beam
		self.meta = kwargs
		if type(stack)==str:
			stack = pd.read_csv(stack)
			stack = [{i:r[i] for i in stack.columns if not (np.isnan(r[i]) if type(r[i])==float else False)} for n,r in stack.iterrows()]
		self.stack = stack

	def check_db(self, db=None):
		if db is not None:
			path, fnm = os.path.split(db)
			if path in ['',' ']:
				path = os.getcwd()
			if fnm in ['',' ']:
				raise ValueError('Invalid db Filename: {}'.format(db))
			db_fnm = os.path.join(path, fnm)

			global ZG_CONNECTIONS_DICT
			if os.path.exists(db_fnm):
				if db_fnm not in ZG_CONNECTIONS_DICT:
					ZG_CONNECTIONS_DICT[db_fnm] = sqlite3.connect(db_fnm)
				self.db_connection = ZG_CONNECTIONS_DICT[db_fnm]
				self.db = self.db_connection.cursor()
			else:
				print('WARNING: DB {} does not exist, creating new file.'.format(fnm))
				from sqlite3 import Error
				try:
					self.db_connection = sqlite3.connect(db_fnm)
					ZG_CONNECTIONS_DICT[db_fnm] = self.db_connection
					self.db = self.db_connection.cursor()
				except Error as e:
					print(e)

	@property
	def meta(self):
		return self._meta

	@meta.setter
	def meta(self, meta_dict):
		for nm in meta_dict:
			self._meta[nm] = meta_dict[nm]
			if nm in ['beam_istp', 'istp']:
				_ITP = Isotope(meta_dict[nm])
				self._meta['Z'] = _ITP.Z
				self._meta['amu'] = _ITP.mass
		if self._meta['min_steps']>self._meta['max_steps']:
			self._meta['max_steps'] = self._meta['min_steps']+1
			print('WARNING: min_steps > max_steps, setting max_steps to {}'.format(self._meta['max_steps']))
		self._meta['solved'] = False

	def __getitem__(self, key):
		if type(key)==str:
			for s in self.stack:
				if s['name']==key:
					return s
			return None
		else:
			return self.stack[int(key)]

	@property
	def stack(self):
		if not self.meta['solved']:
			self._solve()
		return self._stack

	@stack.setter
	def stack(self, _stack):
		self._stack = list(_stack)
		self._meta['solved'] = False
		for s in self._stack:
			if 'name' not in s:
				s['name'] = None
			if 'compound' not in s:
				raise ValueError('compound must be specified')

			if type(s['compound'])==dict:
				cs = ''
				for c in s['compound']:
					cs = c
					self.compounds[c] = s['compound'][c]
				s['compound'] = cs
				if s['compound']=='':
					raise ValueError('compound must be specified')

			if 'ad' not in s:
				if 'area' in s and 'mass' in s:
					s['ad'] = 1e3*s['mass']/s['area']
				elif 'density' in s and 'thickness' in s:
					s['ad'] = 100.0*s['density']*s['thickness']
				elif s['compound'] in self.densities and 'thickness' in s:
					s['ad'], s['density'] = 100.0*self.densities[s['compound']]*s['thickness'], self.densities[s['compound']]

			if 'density' not in s:
				if s['compound'] in self.densities:
					s['density'] = self.densities[s['compound']]
					if 'thickness' not in s and 'ad' in s:
						s['thickness'] = s['ad']/(100.0*s['density'])

			if 'ad' not in s:
				raise ValueError('Areal density either not specified or not computable: {}'.format(s))

			if s['compound'] not in self.compounds:
				cm = s['compound']
				nums = [str(i) for i in range(10)]
				cm_list = []
				for el,Z in self.elements:
					if len(cm)==0:
						break
					f = cm.split(el)
					if len(f)>1:
						c = f[0]
						for e in f[1:]:
							if len(e):
								if e[0] in nums:
									cm_list.append([Z, float(e[0])])
									if len(e)>1:
										c += e[1:]
								else:
									cm_list.append([Z, 1.0])
									c += e
							else:
								cm_list.append([Z, 1.0])
						cm = c
				if len(cm)==0 and len(cm_list):
					self.compounds[s['compound']] = [[i[0], i[1]/sum([i[1] for i in cm_list])] for i in cm_list] 

			if s['compound'] not in self.compounds:
				raise ValueError('compound {} not known.'.format(s['compound']))


	def get_S(self, E, cm):
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

		# energy E in MeV , stopping power in MeV/(mg/cm2)
		E, r0 = np.asarray(E), False
		if not E.shape:
			E, r0 = np.array([E]), True
		S = np.zeros(len(E))
		A_ave = sum([self.weights[z2][0]*w for z2,w in self.compounds[cm]])
		for z2, w in self.compounds[cm]:
			S_nucl = self.get_S_nucl(E, self.meta['Z'], self.meta['amu'], z2, self.weights[z2][0])
			if self.meta['Z']==1:
				S += w*(S_nucl+self.get_S_p(E, z2, self.meta['amu']))
			elif self.meta['Z']==2:
				S += w*(S_nucl+self.get_S_He(E, z2, self.meta['amu']))
			else:
				S += w*(S_nucl+self.get_S_elec(E, z2, self.meta['amu'], self.meta['Z']))
		return S[0]*0.6022140857/A_ave if r0 else S*0.6022140857/A_ave

	def get_S_nucl(self, E, z1, m1, z2, m2):
		RM = (m1+m2)*np.sqrt((z1**(2/3.0)+z2**(2/3.0)))
		ER = 32.53*m2*1E3*E/(z1*z2*RM)
		return (0.5*np.log(1.0+ER)/(ER+0.10718+ER**0.37544))*8.462*z1*z2*m1/RM

	def get_S_p(self, eng, z2, M1=1.00727647):
		S = np.zeros(len(eng))
		E = 1E3*eng/M1
		A = self.protons[z2]
		beta_sq = np.where(E>=1E3,1.0-1.0/(1.0+E/931478.0)**2,0.9)
		B0 = np.where(E>=1E3,np.log(A[6]*beta_sq/(1.0-beta_sq))-beta_sq,0.0)
		Y = np.log(E[(1E3<=E)&(E<=5E4)])
		B0[np.nonzero(np.where((1E3<=E)&(E<=5E4),B0,0))] -= A[7]+A[8]*Y+A[9]*Y**2+A[10]*Y**3+A[11]*Y**4
		S[E>=1E3] = (A[5]/beta_sq[E>=1E3])*B0[E>=1E3]
		S_low = A[1]*E[(10<=E)&(E<1E3)]**0.45
		S_high = (A[2]/E[(10<=E)&(E<1E3)])*np.log(1.0+(A[3]/E[(10<=E)&(E<1E3)])+A[4]*E[(10<=E)&(E<1E3)])
		S[(10<=E)&(E<1E3)] = S_low*S_high/(S_low+S_high)
		S[(0<E)&(E<10)] = A[0]*E[(0<E)&(E<10)]**0.5
		return S

	def get_S_He(self,eng,z2,M1=4.003):
		S = np.zeros(len(eng))
		E = eng*4.0015/M1
		E = np.where(E>=0.001,E,0.001)
		A = self.helium[z2]
		S_low = A[0]*(1E3*E[E<=10])**A[1]
		S_high = (A[2]/E[E<=10])*np.log(1.0+(A[3]/E[E<=10])+A[4]*E[E<=10])
		S[E<=10] = S_low*S_high/(S_low+S_high)
		Y = np.log(1.0/E[E>10])
		S[E>10] = np.exp(A[5]+A[6]*Y+A[7]*Y**2+A[8]*Y**3)
		return S

	def get_S_elec(self, eng, z2, M1, z1):
		S = np.zeros(len(eng))
		E_keV = 1E3*eng
		S[E_keV/M1<1000] = self.get_eff_Z_ratio(E_keV[E_keV/M1<1000],z1,M1)**2*self.get_S_p(eng[E_keV/M1<1000],z2,M1)
		Y = E_keV[E_keV/M1>=1000]/M1
		beta_sq = 1.0-1.0/(1.0+Y/931478.0)**2
		FX = np.log(2E6*0.511003*beta_sq/(1.0-beta_sq))-beta_sq
		ZHY = 1.0-np.exp(-0.2*np.sqrt(Y)-0.0012*Y-0.00001443*Y**2)
		Z1EFF = self.get_eff_Z_ratio(E_keV[E_keV/M1>=1000],z1,M1)*ZHY
		S[E_keV/M1>=1000] = 4E-1*np.pi*(1.9732857/137.03604)**2*Z1EFF**2*z2*(FX-np.log(self.ionization[z2][0]))/(0.511003*beta_sq)
		return S

	def get_eff_Z_ratio(self, E_keV, z1, M1):
		if z1==1:
			return np.ones(len(eng))
		elif z1==2:
			Y = np.log(E_keV/M1)
			return z1*(1.0-np.exp(-0.7446-0.1429*Y-0.01562*Y**2+0.00267*Y**3-0.000001325*Y**8))
		elif z1==3:
			Y = E_keV/M1
			return z1*(1.0-np.exp(-0.7138-0.002797*Y-0.000001348*Y**2))
		BB = -0.886*np.sqrt(0.04*E_keV/M1)/z1**(2/3.0)
		return z1*(1.0-np.exp(BB-0.0378*np.sin(0.5*np.pi*BB))*(1.034-0.1777*np.exp(-0.08114*z1)))

	def _calc_bins(self):
		return np.arange(0.0, self.meta['E0']+10.0*self.meta['dE0'], min([0.1, self.meta['E0']/500.0]))

	def _solve_chunk(self, N):
		E0 = self.meta['E0']+self.meta['dE0']*np.random.normal(size=int(N))
		bins = self._calc_bins()
		hists = []
		dp = self.meta['dp']
		for n, sm in enumerate(self._stack):
			E_bar = [E0]
			if np.average(E0)<=0.0:
				hists.append(np.concatenate([[N],np.zeros(len(bins)-2)]))
			else:
				steps = int((1.0/self.meta['accuracy'])*sm['ad']*dp*self.get_S(np.average(E0), sm['compound'])/np.average(E0))
				steps = min([max([self.meta['min_steps'], steps]), self.meta['max_steps']])
				dr = (1.0/float(steps))
				for i in range(steps):
					S1 = self.get_S(E0, sm['compound'])
					E1 = E0 - dr*dp*sm['ad']*S1
					E1 = np.where(E1>0, E1, 0.0)
					E0 = E0 - dr*0.5*dp*sm['ad']*(S1+self.get_S(E1, sm['compound']))
					E0 = np.where(E0>0, E0, 0.0)
					E_bar.append(E0)
				hists.append(np.histogram(np.concatenate(E_bar), bins=bins)[0])
		return hists

	def _solve(self):
		if self.meta['solved']:
			return
		self.meta['solved'] = True
		dN = np.linspace(0, self.meta['N'], int(np.ceil(self.meta['N']/float(self.meta['chunk_size'])))+1, dtype=int)
		histos = list(map(self._solve_chunk, dN[1:]-dN[:-1]))

		bins = self._calc_bins()
		energy = 0.5*(bins[1:]+bins[:-1])
		for n,sm in enumerate(self._stack):
			sm['flux'] = np.sum([h[n] for h in histos], axis=0)
			sm['flux'] = sm['flux']/np.sum(sm['flux'])
			sm['mu_E'] = np.sum(sm['flux']*energy)
			sm['sig_E'] = np.sqrt(np.sum(sm['flux']*(energy-sm['mu_E'])**2))
			lh = np.where(sm['flux']>0)[0]
			if lh.size:
				if lh[0]==0:
					nm = sm['name'] if sm['name'] is not None else sm['compound']+str(n+1)
					print('WARNING: Beam stopped in foil {}'.format(nm))
			sm['flux'] = sm['flux'][lh[0]:lh[-1]]
			sm['energy'] = energy[lh[0]:lh[-1]]
			sm['bins'] = bins[lh[0]:lh[-1]+1]

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

		cols = ['name','compound','thickness','density','ad','mu_E','sig_E']
		stack = pd.DataFrame([{c:(sm[c] if c in sm else None) for c in cols} for sm in self.stack], columns=cols)
		cols = ['name','energy','flux']
		fluxes = pd.concat([pd.DataFrame({c:sm[c] for c in cols}, columns=cols) for sm in self.stack if sm['name'] is not None], ignore_index=True)
		for fl in fnms:
			if any([fl.endswith(e) for e in ['.png','.pdf','.eps','.pgf','.ps','.raw','.rgba','.svg','.svgz']]):
				self.plot(saveas=fl, show=False)
			if fl.endswith('.csv'):
				stack.to_csv(fl.replace('.csv','_stack.csv'), index=False)
				fluxes.to_csv(fl.replace('.csv','_fluxes.csv'), index=False)
			if fl.endswith('.db'):
				self.check_db(fl)
				stack.to_sql('stack', self.db_connection, if_exists='replace', index=False)
				fluxes.to_sql('fluxes', self.db_connection, if_exists='replace', index=False)

	def summarize(self, samples=None):
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

		for n,sm in enumerate(self.stack):
			if sm['name'] is not None:
				if samples is not None:
					if not any([re.match(s, sm['name']) for s in samples]):
						continue
				nm = sm['name'] if sm['name'] is not None else sm['compound']+str(n+1)
				print(nm+': '+str(round(sm['mu_E'], 2))+' +/- '+str(round(sm['sig_E'], 2))+' (MeV)')

	def plot_S(self, compound, energy=None, **kwargs):
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

		if energy is None:
			energy = 10.0**np.arange(0.1,2.8,0.1)
		f,ax = _init_plot(**kwargs)
		ax.plot(energy, self.get_S(energy, compound), label=compound.title())
		ax.set_xlabel('Energy (MeV)')
		ax.set_ylabel(r'Stopping Power (MeV$\cdot$mg$^{-1}\cdot$cm$^{-2}$)')
		ax.set_xscale('log')
		ax.legend(loc=0)
		return _close_plot(f, ax, **kwargs)

	def plot(self, samples=None,  **kwargs):
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

		if type(samples)==str:
			samples = [samples]
		f,ax = _init_plot(**kwargs)
		for sm in self.stack:
			if sm['name'] is not None:
				if samples is not None:
					if not any([re.match(s, sm['name']) for s in samples]):
						continue
				x, y = np.array([sm['bins'][:-1],sm['bins'][1:]]).T.flatten(), np.array([sm['flux'],sm['flux']]).T.flatten()
		
				ax.plot(x,y,label=sm['name'])

		ax.set_xlabel('Energy (MeV)')
		ax.set_ylabel('Flux (a.u.)')
		ax.legend(loc=0)
		return _close_plot(f, ax, **kwargs)
