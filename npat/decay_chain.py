from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import datetime as dtm
from scipy.optimize import curve_fit

from .isotope import Isotope
from .plotter import colors
from .plotter import _init_plot
from .plotter import _close_plot

class DecayChain(object):
	"""Calculating activities within a decay chain using Bateman equations

	...

	Parameters
	----------

	Attributes
	----------

	Methods
	-------

	"""
	def __init__(self, parent, units='s', R=None, A0=None, time=None):
		self.units = units

		istps = [Isotope(parent)]
		self.isotopes = [istps[0].name]
		self.chain = [[istps[0].decay_const(units), [], []]]
		stable_chain = [False]

		while not all(stable_chain):
			for n in [n for n,ch in enumerate(stable_chain) if not ch]:
				stable_chain[n] = True
				for prod,br in istps[n].decay_products():
					I = Isotope(prod)
					if I.name in self.isotopes:
						self.chain[self.isotopes.index(I.name)][1].append(br)
						self.chain[self.isotopes.index(I.name)][2].append(n)
					else:
						istps.append(I)
						self.isotopes.append(I.name)
						self.chain.append([I.decay_const(units), [br], [n]])
						stable_chain.append(self.chain[-1][0]<1E-12)
		self.chain = np.array(self.chain, dtype=object)

		self.set_R(R)
		self.A0 = np.zeros(len(self.chain))
		self.update_A0(A0)

		self.time = time

		self._counts = [[] for i in self.chain]
		self._others = []
		self._prev = []
		self._R_fit = None
		self._A0_fit = None
		self._EoB = None

	def filter_name(self, istp):
		if istp[-1] in 'g0123456789':
			return istp
		return istp+('1' if istp.endswith('m') else 'g')

	def index(self, istp):
		return self.isotopes.index(self.filter_name(istp))

	def _get_branches(self, istp):
		if self.filter_name(istp) not in self.isotopes:
			return [], []
		m = self.index(istp)
		BR = [[0.0]]+[[r] for r in self.chain[m,1]]
		CH = [[m]]+[[n] for n in self.chain[m,2]]
		while not all([c[-1]==0 for c in CH]):
			BR = [BR[n]+[i] for n,c in enumerate(CH) for i in self.chain[c[-1],1]]
			CH = [c+[i] for c in CH for i in self.chain[c[-1],2]]
		BR = [np.array(r)[::-1] for n,r in enumerate(BR)]
		CH = [np.array(c)[::-1] for n,c in enumerate(CH)]
		return BR, CH

	def activity(self, istp, t=None, units=None, _R=None, _A0=None):
		t = t if t is not None else self.time
		if t is None:
			raise ValueError('Time must be specified.')
		t = np.asarray(t)
		A = np.zeros(len(t)) if t.shape else 0.0
		R_lm = 1.0
		if units is not None:
			conv = {'ns':1e-9,'us':1e-6,'ms':1e-3,'s':1.0,'m':60.0,'h':3600.0,'d':86400.0,'y':31557.6E3,'ky':31557.6E6}
			R_lm = conv[units]/conv[self.units]

		for m,(BR, chain) in enumerate(zip(*self._get_branches(istp))):
			lm = R_lm*self.chain[chain,0]
			A0 = self.A0[chain] if _A0 is None else _A0[chain]
			R = self.R[chain] if _R is None else _R[chain]
			n = len(chain)
			for i in range(n):
				if i==n-1 and m>0:
					continue

				A_i = lm[-1]*(A0[i]/lm[i])
				B_i = np.prod(lm[i:-1]*BR[i:-1])

				for j in range(i, n):
					K = np.arange(i, n)
					C_j = np.prod(lm[K[K!=j]]-lm[j])
					if np.any((lm[K[K!=j]]-lm[j])==0):
						C_j = np.where((lm[K[K!=j]]-lm[j])!=0,lm[K[K!=j]]-lm[j],1E-12)
						C_j = np.prod(C_j)
					A += A_i*B_i*np.exp(-lm[j]*t)/C_j
					if lm[j]>1E-12:
						A += R[i]*lm[-1]*B_i*(1.0-np.exp(-lm[j]*t))/(lm[j]*C_j)
					else:
						A += R[i]*lm[-1]*B_i*t/C_j

		for ot in self._others:
			if self.filter_name(istp) in ot.isotopes:
				A += ot.activity(istp, t, units=(self.units if units is None else units))
		return A


	def decays(self, istp, t_start, t_stop, units=None, _A0=None):
		if np.any(self.R > 0):
			print('WARNING: decays during production not implemented.')
		t_start, t_stop = np.asarray(t_start), np.asarray(t_stop)
		D = np.zeros(len(t_start)) if t_start.shape else (np.zeros(len(t_stop)) if t_stop.shape else 0.0)
		R_lm = 1.0
		conv = {'ns':1e-9,'us':1e-6,'ms':1e-3,'s':1.0,'m':60.0,'h':3600.0,'d':86400.0,'y':31557600.0}
		if units is not None:
			R_lm = conv[units]/conv[self.units]

		for m,(BR, chain) in enumerate(zip(*self._get_branches(istp))):
			lm = R_lm*self.chain[chain,0]
			A0 = self.A0[chain] if _A0 is None else _A0[chain]
			n = len(chain)
			for i in range(n):
				if i==n-1 and m>0:
					continue

				A_i = lm[-1]*(A0[i]/lm[i])
				B_i = np.prod(lm[i:-1]*BR[i:-1])

				for j in range(i, len(chain)):
					K = np.arange(i, len(chain))
					C_j = np.prod(lm[K[K!=j]]-lm[j])
					if np.any((lm[K[K!=j]]-lm[j])==0):
						C_j = np.where((lm[K[K!=j]]-lm[j])!=0,lm[K[K!=j]]-lm[j],1E-12)
						C_j = np.prod(C_j)
					if lm[j]>1E-12:
						D += A_i*B_i*(np.exp(-lm[j]*t_start)-np.exp(-lm[j]*t_stop))/(lm[j]*C_j)
					else:
						D += A_i*B_i*(t_stop-t_start)/C_j
		D = D*conv[(self.units if units is None else units)]
		for ot in self._others:
			if self.filter_name(istp) in ot.isotopes:
				D += ot.decays(istp, t_start, t_stop, units=(self.units if units is None else units))
		return D

	def set_R(self, R):
		self.R = np.zeros(len(self.isotopes))
		if R is not None:
			if type(R)==dict:
				for istp in R:
					if self.filter_name(istp) in self.isotopes:
						self.R[self.index(istp)] += R[istp]
			else:
				self.R[0] += float(R)

	def update_A0(self, A0):
		if A0 is not None:
			if type(A0)==dict:
				for istp in A0:
					if self.filter_name(istp) in self.isotopes:
						self.A0[self.index(istp)] += A0[istp]
			else:
				self.A0[0] += float(A0)

	@property
	def EoB(self):
		return self._EoB

	@EoB.setter
	def EoB(self, eob):
		if type(eob)==str:
			eob = dtm.datetime.strptime(eob, '%m/%d/%Y %H:%M:%S')
		self._EoB = eob
	

	@property
	def counts(self):
		return self._counts

	@counts.setter
	def counts(self, N_c):
		if N_c is not None:
			if type(N_c)!=dict:
				N_c = {self.isotopes[0]:np.array(N_c)}
			for istp in N_c:
				if self.filter_name(istp) in self.isotopes:
					x = np.array(N_c[istp])
					i = self.index(istp)
					if x.shape:
						if len(x.shape)==1:
							x = np.array([x])
						if len(self._counts[i]):
							self._counts[i] = np.append(self._counts[i], x, axis=0)
						else:
							self._counts[i] = x

	@property
	def time(self):
		return self._time

	@time.setter
	def time(self, time):
		self._time = None
		if time is not None:
			t = np.asarray(time)
			if t.shape:
				self._time = t
			else:
				self._time = np.linspace(0.0, float(t), 1000)

	@property
	def A_meas(self):
		A_meas = []
		for n,ct in enumerate(self.counts):
			itp = self.isotopes[n]
			if len(ct):
				meas = np.array([self.calc_M(itp, itp, c[0], c[1])*c[2] for c in ct])
				if len(ct[0])==4:
					A_meas.append(np.column_stack((ct[:,0], meas, ct[:,3]*meas/ct[:,2])))
				else:
					A_meas.append(np.column_stack((ct[:,0], meas)))
			else:
				A_meas.append([])
		return A_meas
	
	def calc_L(self, i, t_m):
		a_0 = self.activity(i, 0.0)
		a_m = self.activity(i, t_m)
		if a_m>0:
			return a_0/a_m
		return 0.0

	def calc_M(self, i, j, t_start, t_stop):
		a_n = self.activity(i, t_start)
		d_m = self.decays(j, t_start, t_stop)
		if d_m>0:
			return a_n/d_m
		return 0.0

	def calc_P(self, i):
		a_0 = self.activity(i, 0.0)
		R_norm = self.R_norm[self.index(i)]
		if a_0>0:
			return R_norm/a_0
		return 0.0

	def calc_Q(self, i, t_m):
		a_m = self.activity(i, t_m)
		R_norm = self.R_norm[self.index(i)]
		if a_m>0:
			return R_norm/a_m
		return 0.0

	@property
	def R_norm(self):
		R_p,t = [p.R for p in self._prev if p.R.any()],[p.time[-1] for p in self._prev if p.R.any()]
		return np.dot(np.array(R_p).T,np.array(t))/np.sum(t)

	def fit_R(self, istp=None, _update=True):
		if self._R_fit is not None:
			if istp is None:
				return self._R_fit
			return self._R_fit[self.index(istp)]

		R_norm = self.R_norm
		nz = np.where(R_norm>0)[0]
		Y = []
		dY = []
		X = []
		time = []
		itp = []
		for n,ct in enumerate(self.counts):
			if len(ct):
				Y += ct[:,2].tolist()
				if len(ct[0])==4:
					dY += ct[:,3].tolist()
				X += [np.zeros(len(nz)) for i in ct]
				time += ct[:,:2].tolist()
				itp += (n*np.ones(len(ct),dtype=int)).tolist()

		Y = np.array(Y)
		dY = np.array(dY) if len(dY) else np.ones(len(Y))
		X = np.asarray(X)

		I = np.arange(len(self.isotopes))
		for m,z in enumerate(nz):
			A0 = np.copy(self._prev[0].A0) if len(self._prev) else np.zeros(len(self.isotopes))
			for p in self._prev:
				R_p = p.R/np.where(R_norm>0, R_norm, 1.0)
				R_p[I!=z] = 0.0
				A_old = np.copy(A0)
				for n in I:
					A0[n] = self.activity(self.isotopes[n], p.time[-1], _R=R_p, _A0=A_old)
			X[:,m] = np.array([self.decays(self.isotopes[itp[n]], time[n][0], time[n][1], _A0=A0) for n in range(len(X))])
		
		func = lambda X_f, *R_f: np.dot(X_f, np.asarray(R_f).T)
		fit, unc = curve_fit(func, X, Y, sigma=dY, p0=R_norm[nz], bounds=([0.0 for i in nz],[np.inf for i in nz]))
		self._R_fit = np.zeros(len(R_norm))
		self._R_fit[nz] = fit

		if _update:
			A0 = np.copy(self._prev[0].A0) if len(self._prev) else np.zeros(len(self.isotopes))
			for p in self._prev:
				p.A0 = np.copy(A0)
				p.R = p.R*self._R_fit/np.where(R_norm>0, R_norm, 1.0)
				for n in I:
					A0[n] = p.activity(self.isotopes[n], p.time[-1])
			self.A0 = np.copy(A0)

		if istp is None:
			return self._R_fit
		return self._R_fit[self.index(istp)]

	def fit_A0(self, istp=None, _update=True):
		if self._A0_fit is not None:
			if istp is None:
				return self._A0_fit
			return self._A0_fit[self.index(istp)]

		if len(self._prev):
			print('WARNING: Previous activities and production rates will not be fit.')

		nz = np.where(self.A0>0)[0]

		Y = []
		dY = []
		X = []
		time = []
		itp = []
		for n,ct in enumerate(self.counts):
			if len(ct):
				Y += ct[:,2].tolist()
				if len(ct[0])==4:
					dY += ct[:,3].tolist()
				X += [np.zeros(len(nz)) for i in ct]
				time += ct[:,:2].tolist()
				itp += (n*np.ones(len(ct),dtype=int)).tolist()

		Y = np.array(Y)
		dY = np.array(dY) if len(dY) else np.ones(len(Y))
		X = np.asarray(X)

		I = np.arange(len(self.isotopes))
		for z in nz:
			A0 = np.where(I!=z, np.zeros(len(self.isotopes)), 1.0)
			X[:,z] = np.array([self.decays(self.isotopes[itp[n]], time[n][0], time[n][1], _A0=A0) for n in range(len(X))])

		func = lambda X_f, *R_f: np.dot(X_f, np.asarray(R_f).T)
		fit, unc = curve_fit(func, X, Y, sigma=dY, p0=self.A0[nz], bounds=([0.0 for i in nz],[np.inf for i in nz]))
		self._A0_fit = np.copy(self.A0)
		self._A0_fit[self.A0>0] = fit

		if _update:
			self.A0 = np.copy(self._A0_fit)

		if istp is None:
			return self._A0_fit
		return self._A0_fit[self.index(istp)]


	def __add__(self, other):
		self._others.append(other)
		return self

	def append(self, other, time=None):
		if time is not None:
			self.time = time
		if self.time is None:
			raise ValueError('Decay Chain time is not set.')
		if self.isotopes[0]!=other.isotopes[0]:
			raise ValueError('Decay chains must have same parent.')
		if self.units!=other.units:
			raise ValueError('Units must be the same.')

		tmp = [self.R, self.A0, self.time, self._others, self.counts]
		new_A0 = {i:self.activity(i, self.time[-1]) for n,i in enumerate(self.isotopes) if self.chain[n,0]>1E-12}
		self.R = other.R
		self.A0 = other.A0
		self.update_A0(new_A0)
		self.time = other.time
		self._others = other._others
		self._counts = other.counts
		other.R = tmp[0]
		other.A0 = tmp[1]
		other.time = tmp[2]
		other._others = tmp[3]
		other._counts = tmp[4]
		self._prev += other._prev
		self._prev.append(other)
		other._prev = []

	def plot(self, time=None, N_plot=None, N_label=10, **kwargs):
		f, ax = _init_plot(**kwargs)
		if time is not None:
			self.time = time
		if self.time is not None:
			t = self.time
		else:
			t = np.linspace(0, 5.0*np.log(2)/self.chain[0,0], 1000)
		if N_plot is None:
			N_plot = len(self.isotopes)
		A_meas = self.A_meas
		ordr = int(np.floor(np.log10(np.average(self.A0))/3.0))
		lb_or = {-4:'p',-3:'n',-2:r'$\mu$',-1:'m',0:'',1:'k',2:'M',3:'G',4:'T'}[ordr]
		mult = 10.0**(-3*ordr)
		for n,istp in enumerate(self.isotopes):
			if self.chain[n,0]>1E-12 and n<N_plot:
				dt = 0.0
				tm, A = np.array([]),np.array([])
				for p in self._prev:
					dt += p.time[-1]
					t_p, A_p = p.time, p.activity(istp)
					tm, A = np.append(tm, t_p+(tm[-1] if len(tm) else 0.0)), np.append(A, A_p)
				tm, A = np.append(tm-dt, t), np.append(A, self.activity(istp, t))
				label = Isotope(istp).TeX if n<N_label else None
				line, = ax.plot(tm, mult*A, label=label)
				if len(A_meas[n]):
					am = A_meas[n]
					if len(am[0])==3:
						ax.errorbar(am[:,0], mult*am[:,1], yerr=mult*am[:,2], ls='None', marker='o', color=line.get_color())
					else:
						ax.plot(am[:,0], mult*am[:,1], ls='None', marker='o', color=line.get_color())

		ax.set_xlabel('Time ({})'.format(self.units))
		ax.set_ylabel('Activity ({}Bq)'.format(lb_or))
		ax.legend(loc=0)
		return _close_plot(f, ax, default_log=False, **kwargs)

