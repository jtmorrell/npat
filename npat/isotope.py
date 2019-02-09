from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import re
import numpy as np

from .dbmgr import get_cursor

class Isotope(object):
	"""Isotopic data

	...

	Parameters
	----------

	Attributes
	----------

	Methods
	-------

	"""
	def __init__(self, istp):
		if istp=='1n' or istp=='1ng':
			self.element, self.A, self.isomer = 'n', 1, 'g'
		else:
			self.element = ''.join(re.findall('[A-Z]+',istp))
			self.A = int(istp.split(self.element)[0])
			self.isomer = istp.split(self.element)[1]
		self.isotope = str(self.A)+self.element
		if self.isomer=='':
			self.isomer = 'g'
		if self.isomer=='m':
			self.isomer = 'm1'
		self.name = self.isotope+self.isomer
		self._meta = None
			
	@property
	def meta(self):
		if self._meta is None:
			self._meta = {}
			self.db = get_cursor('decay')
			for i in ['SFY','gm','el','bm','bp','al']:
				self._meta[i] = None
			Q = list(self.db.execute('SELECT * FROM chart WHERE isotope=?',(self.isotope,)))
			q = [i for i in Q if i[3]==self.isomer][0]
			self._meta['E_level'] = q[2]
			self._meta['J_pi'] = str(q[4]) if q[4] else '?'
			self._meta['Z'] = q[6]
			self._meta['N'] = q[7]
			self._meta['stable'] = bool(q[9])
			self._meta['t_half'] = q[10]
			self._meta['unc_t_half'] = q[11]
			self._meta['abundance'] = q[12]
			self._meta['unc_abundance'] = q[13]
			self._meta['mass'] = q[14]
			self._meta['Delta'] = q[15]
			self._meta['decay_mode'] = list(map(lambda i:[float(p) if n==2 else p for n,p in enumerate(i.split(':'))], str(q[16]).split(',')))
			state = '' if len(Q)==1 else (self.isomer[0] if len(Q)==2 else self.isomer)
			self._meta['TeX'] = r'$^{'+str(self.A)+state+r'}$'+self.element.title()
		return self._meta
	
	@property
	def E_level(self):
		return self.meta['E_level']

	@property
	def J_pi(self):
		return self.meta['J_pi']

	@property
	def Z(self):
		return self.meta['Z']

	@property
	def N(self):
		return self.meta['N']

	@property
	def stable(self):
		return self.meta['stable']

	@property
	def mass(self):
		return self.meta['mass']

	@property
	def Delta(self):
		return self.meta['Delta']

	@property
	def TeX(self):
		return self.meta['TeX']

	def __str__(self):
		return self.name

	def half_life(self, units='s', unc=False):
		if self.stable:
			return (np.inf,0.0) if unc else np.inf
		half_conv = {'ns':1e-9,'us':1e-6,'ms':1e-3,'s':1.0,'m':60.0,'h':3600.0,'d':86400.0,'y':31557.6E3,'ky':31557.6E6}[units]
		if unc:
			return self.meta['t_half']/half_conv,self.meta['unc_t_half']/half_conv
		return self.meta['t_half']/half_conv

	def decay_const(self, units='s', unc=False):
		if self.stable:
			return (0.0,0.0) if unc else 0.0
		if unc:
			T2,uT2 = self.half_life(units, True)
			return np.log(2.0)/T2,np.log(2.0)*uT2/T2**2
		return np.log(2.0)/self.half_life(units)

	def optimum_units(self):
		opt = ['ns']
		for units in ['us','ms','s','m','h','d','y']:
			if self.half_life(units)>1.0:
				opt.append(units)
		return opt[-1]

	def abundance(self, unc=False):
		if unc:
			return self.meta['abundance'], self.meta['unc_abundance']
		return self.meta['abundance']

	def get_SFY(self, unc=False, closest_SFY=False):
		if self.meta['SFY'] is None:
			self.meta['SFY'] = [(str(i[1]),i[2],i[3]) for i in self.db.execute('SELECT * FROM SFY WHERE parent=?',(self.name,))]
		SFY = list(self.meta['SFY'])
		if len(SFY)==0 and closest_SFY:
			itps = list(set([str(i[0]) for i in self.db.execute('SELECT parent FROM SFY')]))
			dA = [abs(self.A-int(i[:3])) for i in itps]
			SFY = [(str(i[1]),i[2],i[3]) for i in self.db.execute('SELECT * FROM SFY WHERE parent=?',(itps[dA.index(min(dA))],))]
		if unc:
			return [[i[0], i[1], i[2]] for i in SFY]
		return [[i[0], i[1]] for i in SFY]

	def decay_products(self, closest_SFY=False):
		prods = []
		for (mode, product, br) in self.meta['decay_mode']:
			if product=='SFY':
				prods += [[i, br*y] for i,y in self.get_SFY(False, closest_SFY)]
			else:
				prods.append([product, br])
		return [[i, br] for i,br in prods if br>1E-8]

	def gammas(self,I_lim=[None,None],E_lim=[None,None],xrays=False):
		if self.meta['gm'] is None:
			self.meta['gm'] = [[float(i[3]),float(i[4]),float(i[5]),str(i[6])] for i in self.db.execute('SELECT * FROM gammas WHERE isotope=? AND isomer=?',(self.isotope,self.isomer))]
		gammas = list(self.meta['gm'])
		for n,L in enumerate([E_lim,I_lim]):
			if L[0] is not None:
				gammas = [g for g in gammas if g[n]>=L[0]]
			if L[1] is not None:
				gammas = [g for g in gammas if g[n]<=L[1]]
		if not xrays:
			gammas = [g for g in gammas if g[3]=='' and abs(g[0]-511.0)>4.0]
		return {l:[g[n] for g in gammas] for n,l in enumerate(['E','I','dI','notes'])}

	def electrons(self,I_lim=(None,None),E_lim=(None,None),CE_only=False,Auger_only=False):
		if self.meta['el'] is None:
			self.meta['el'] = [[float(i[3]),float(i[4]),float(i[5]),str(i[6])] for i in self.db.execute('SELECT * FROM electrons WHERE isotope=? AND isomer=?',(self.isotope,self.isomer))]
		electrons = list(self.meta['el'])
		for n,L in enumerate([E_lim,I_lim]):
			if L[0] is not None:
				electrons = [e for e in electrons if e[n]>=L[0]]
			if L[1] is not None:
				electrons = [e for e in electrons if e[n]<=L[1]]
		if CE_only:
			electrons = [e for e in electrons if e[3].startswith('CE')]
		if Auger_only:
			electrons = [e for e in electrons if e[3].startswith('Aug')]
		return {l:[e[n]for e in electrons] for n,l in enumerate(['E','I','dI','notes'])}

	def beta_minus(self,I_lim=(None,None),Endpoint_lim=(None,None)):
		if self.meta['bm'] is None:
			self.meta['bm'] = [[float(i[3]),float(i[4]),float(i[5]),float(i[6])] for i in self.db.execute('SELECT * FROM beta_minus WHERE isotope=? AND isomer=?',(self.isotope,self.isomer))]
		betas = list(self.meta['bm'])
		for n,L in zip([3,1],[Endpoint_lim,I_lim]):
			if L[0] is not None:
				betas = [b for b in betas if b[n]>=L[0]]
			if L[1] is not None:
				betas = [b for b in betas if b[n]<=L[1]]
		return {l:[b[n] for b in betas] for n,l in enumerate(['muE','I','dI','endE'])}

	def beta_plus(self,I_lim=(None,None),Endpoint_lim=(None,None)):
		if self.meta['bp'] is None:
			self.meta['bp'] = [[float(i[3]),float(i[4]),float(i[5]),float(i[6])] for i in self.db.execute('SELECT * FROM beta_plus WHERE isotope=? AND isomer=?',(self.isotope,self.isomer))]
		betas = list(self.meta['bp'])
		for n,L in zip([3,1],[Endpoint_lim,I_lim]):
			if L[0] is not None:
				betas = [b for b in betas if b[n]>=L[0]]
			if L[1] is not None:
				betas = [b for b in betas if b[n]<=L[1]]
		return {l:[b[n] for b in betas] for n,l in enumerate(['muE','I','dI','endE'])}

	def alphas(self,I_lim=(None,None),E_lim=(None,None)):
		if self.meta['al'] is None:
			self.meta['al'] = [[float(i[3]),float(i[4]),float(i[5])] for i in self.db.execute('SELECT * FROM alphas WHERE isotope=? AND isomer=?',(self.isotope,self.isomer))]
		alphas = list(self.meta['al'])
		for n,L in enumerate([E_lim,I_lim]):
			if L[0] is not None:
				alphas = [a for a in alphas if a[n]>=L[0]]
			if L[1] is not None:
				alphas = [a for a in alphas if a[n]<=L[1]]
		return {l:[a[n] for a in alphas] for n,l in enumerate(['E','I','dI'])}

