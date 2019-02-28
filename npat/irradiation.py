from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import re
import numpy as np
import matplotlib.pyplot as plt
import datetime as dtm
from scipy.interpolate import interp1d

from .plotter import colors
from .plotter import _init_plot
from .plotter import _close_plot
from .dbmgr import get_cursor
from .isotope import Isotope


class Irradiation(object):
	"""Super class for activation experiments.

	...

	Parameters
	----------

	Attributes
	----------

	Methods
	-------

	"""
	def __init__(self):
		self._samples = []

	@property
	def samples(self):
		return self._samples
	
	def read_csv(self, filename):
		with open(filename) as f:
			lines = f.read().split('\n')

		dat = [[i if i!='' else None for i in l.split(',')] for l in lines if not l.startswith('#')]
		dat = [l for l in dat if len(l) and any(l)]
		
		if lines[0].startswith('#'):
			head = lines[0].replace('#','').split(',')
		else:
			head = [None for i in dat[0]]
		
		csv = []
		for n in range(len(head)):
			col = [d[n] for d in dat]
			try:
				csv.append([float(c) if c is not None else None for n,c in enumerate(col)])
			except:
				csv.append(col)

		return csv, head

	def save_csv(self, filename, csv, head=None):
		with open(filename, 'w+') as f:
			hs = '' if head is None else '#'+','.join(head)+'\n'
			f.write(hs+'\n'.join(map(','.join, zip(*[map(str,i) for i in csv]))))





class Sample(object):
	"""Material properties of samples in experiment

	...

	Parameters
	----------

	Attributes
	----------

	Methods
	-------

	"""
	def __init__(self):
		self.name = None
		self.energy = None
		self.flux = None
		self.bins = None

	def _set_hist(self, hist):
		self.flux = hist[0]
		self.bins = hist[1]
		self.energy = 0.5*(self.bins[1:]+self.bins[:-1])
		self.mu_E = np.sum(self.flux*self.energy)/np.sum(self.flux)
		self.sig_E = np.sqrt(np.sum(self.flux*(self.energy-self.mu_E)**2)/np.sum(self.flux))
		if len(self.energy)<2:
			self._flux_interp = interp1d([0, 1], [0, 0], bounds_error=False, fill_value=0.0)
		else:
			self._flux_interp = interp1d(self.energy, self.flux, bounds_error=False, fill_value=0.0)

	def reaction_integral(self, sig):
		# Trapezoidal Riemann sum
		phisig = self.flux*np.asarray(sig)
		return np.sum(0.5*(self.energy[1:]-self.energy[:-1])*(phisig[:-1]+phisig[1:]))

	def interp_flux(self, E=None):
		if E is None:
			return self._flux_interp
		return self._flux_interp(E)
	
	def plot(self, ax=None, **kwargs):
		x, y = np.array([self.bins[:-1],self.bins[1:]]).T.flatten(), np.array([self.flux,self.flux]).T.flatten()
		if ax is None:
			f, ax = _init_plot(**kwargs)
			ax.plot(x, y, label=self.name)
			ax.set_xlabel('Energy (MeV)')
			ax.set_ylabel('Flux (a.u.)')
			return _close_plot(f, ax, **kwargs)
		else:
			ax.plot(x, y, label=self.name)
		


class Ziegler(Irradiation):
	"""Method for solving energy loss within stacked-target foil experiment

	...

	Parameters
	----------

	Attributes
	----------

	Methods
	-------

	"""
	def __init__(self, stack=[], beam={}):
		### stack is list of dicts, which must have 'compound' specified.
		### Areal density specified by either 'ad' (mg/cm^2), both mass (g) and area (cm^2),
		### 'thickness' (mm) or both 'density' (g/cm^3) and 'thickness' (mm)

		db = get_cursor('ziegler')
		
		self.protons = {int(i[0]):list(map(float,i[1:])) for i in db.execute('SELECT * FROM protons')}
		self.helium = {int(i[0]):list(map(float,i[1:])) for i in db.execute('SELECT * FROM helium')}
		self.ionization = {int(i[0]):list(map(float,i[1:])) for i in db.execute('SELECT * FROM ionization')}
		self.weights = {int(i[0]):list(map(float,i[1:3])) for i in db.execute('SELECT * FROM weights')}
		self.compounds = {str(i[0]):[[int(h.split(':')[0]), float(h.split(':')[1])] for h in i[2].split(',')] for i in db.execute('SELECT * FROM compounds')}
		self.compounds = {cm:[[i[0],i[1]/sum([m[1] for m in self.compounds[cm]])] for i in self.compounds[cm]] for cm in self.compounds}
		self.densities = {str(i[0]):float(i[1]) for i in db.execute('SELECT * FROM compounds')}
		self.elements = sorted([[str(i[0]), int(i[2].split(':')[0])] for i in db.execute('SELECT * FROM compounds') if len(i[2].split(':'))==2], key=lambda h:len(h[0]), reverse=True)

		Irradiation.__init__(self)
		self.beam = beam
		self.stack = stack

	@property
	def beam(self):
		return self._beam

	@beam.setter
	def beam(self, _beam):
		self._beam = _beam
		default = {'istp':'1H', 'E0':33.0, 'dE0':0.3, 'N':10000}
		for k in default:
			if k not in _beam:
				self._beam[k] = default[k]
		_ITP = Isotope(self._beam['istp'])
		self._beam['Z'] = _ITP.Z
		self._beam['amu'] = _ITP.mass
		if len(self._samples):
			self._solve()

	@property
	def stack(self):
		return self._samples

	@stack.setter
	def stack(self, _stack):
		self._samples = []
		for s in _stack:
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
				s['density'] = self.densities[s['compound']]
				if 'thickness' not in s and 'ad' in s:
					s['thickness'] = s['ad']/(100.0*s['density'])

			sm = Sample()
			if 'name' in s:
				sm.name = s['name']
			if 'ad' not in s:
				raise ValueError('Areal density either not specified or not computable: {}'.format(s))
			sm.ad = s['ad']

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

			sm.compound = s['compound']
			sm.meta = s
			self._samples.append(sm)
		self._solve()

	def get_S(self, E, cm):
		# energy E in MeV , stopping power in MeV/(mg/cm2)
		E = np.asarray(E)
		S = 0.0
		A_ave = sum([self.weights[z2][0]*w for z2,w in self.compounds[cm]])
		for z2, w in self.compounds[cm]:
			S_nucl = self.get_S_nucl(E, self.beam['Z'], self.beam['amu'], z2, self.weights[z2][0])
			if self.beam['Z']==1:
				S += w*(S_nucl+self.get_S_p(E, z2, self.beam['amu']))
			elif self.beam['Z']==2:
				S += w*(S_nucl+self.get_S_He(E, z2, self.beam['amu']))
			else:
				S += w*(S_nucl+self.get_S_elec(E, z2, self.beam['amu'], self.beam['Z']))
		return S*0.6022140857/A_ave

	def get_S_nucl(self,E,z1,m1,z2,m2):
		RM = (m1+m2)*np.sqrt((z1**(2/3.0)+z2**(2/3.0)))
		ER = 32.53*m2*1E3*E/(z1*z2*RM)
		return (0.5*np.log(1.0+ER)/(ER+0.10718+ER**0.37544))*8.462*z1*z2*m1/RM

	def get_S_p(self,eng,z2,M1=1.00727647):
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

	def get_S_elec(self,eng,z2,M1,z1):
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

	def get_eff_Z_ratio(self,E_keV,z1,M1):
		if z1==1:
			return np.ones(len(E_keV))
		elif z1==2:
			Y = np.log(E_keV/M1)
			return z1*(1.0-np.exp(-0.7446-0.1429*Y-0.01562*Y**2+0.00267*Y**3-0.000001325*Y**8))
		elif z1==3:
			Y = E_keV/M1
			return z1*(1.0-np.exp(-0.7138-0.002797*Y-0.000001348*Y**2))
		BB = -0.886*np.sqrt(0.04*E_keV/M1)/z1**(2/3.0)
		return z1*(1.0-np.exp(BB-0.0378*np.sin(0.5*np.pi*BB))*(1.034-0.1777*np.exp(-0.08114*z1)))

	def _solve(self, dp=1.0):
		E0 = self.beam['E0']+self.beam['dE0']*np.random.normal(size=int(self.beam['N']))
		steps = 2
		warn = True
		for n, sm in enumerate(self._samples):
			E_bar = [E0]
			ds = (1.0/float(steps))
			for i in range(steps):
				S1 = self.get_S(E0, sm.compound)
				E1 = E0 - ds*dp*sm.ad*S1
				E1 = np.where(E1>0, E1, 0.0)
				E0 = E0 - ds*0.5*dp*sm.ad*(S1+self.get_S(E1, sm.compound))
				if np.average(np.abs(E0-E1)/np.where(E0>0, E0, np.ones(len(E0))))>1E-3:
					steps += 1
				E0 = np.where(E0>0, E0, 0.0)
				if not np.all(E0>0) and warn:
					print('WARNING: Beam stopped in foil {}'.format((sm.name if sm.name is not None else sm.compound+str(n+1))))
					warn = False
				E_bar.append(E0)
			self.samples[n]._set_hist(np.histogram(np.average(E_bar, axis=0), bins='auto'))

	def summarize(self, samples=None, printout=True, saveas=None, incl_no_names=False):
		names, mu_E, sig_E = [], [], []
		for sm in self.samples:
			if sm.name is not None:
				if samples is not None or incl_no_names:
					if not any([re.match(s, sm.name) for s in samples]):
						continue
				print(sm.name+': '+str(round(sm.mu_E, 2))+' +/- '+str(round(sm.sig_E, 2))+' (MeV)')
				if saveas is not None:
					names.append(sm.name)
					mu_E.append(sm.mu_E)
					sig_E.append(sm.sig_E)
		if saveas is not None:
			self.save_csv(saveas, [names, mu_E, sig_E], ['Name', 'mu_E (MeV)', 'sig_E (MeV)'])

	def plot_S(self, compound, energy=None, **kwargs):
		if energy is None:
			energy = 10.0**np.arange(0.1,2.8,0.1)
		f,ax = _init_plot(**kwargs)
		ax.plot(energy, self.get_S(energy, compound), label=compound.title())
		ax.set_xlabel('Energy (MeV)')
		ax.set_ylabel(r'Stopping Power (MeV$\cdot$mg$^{-1}\cdot$cm$^{-2}$)')
		ax.set_xscale('log')
		ax.legend(loc=0)
		return _close_plot(f, ax, **kwargs)

	def plot(self, samples=None, incl_no_names=False, **kwargs):
		if type(samples)==str:
			samples = [samples]
		f,ax = _init_plot(**kwargs)
		for sm in self.samples:
			if sm.name is not None or incl_no_names:
				if samples is not None:
					if not any([re.match(s, sm.name) for s in samples]):
						continue
				sm.plot(ax)

		ax.set_xlabel('Energy (MeV)')
		ax.set_ylabel('Flux (a.u.)')
		ax.legend(loc=0)
		return _close_plot(f, ax, **kwargs)
