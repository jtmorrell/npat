from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import re
import numpy as np

from .dbmgr import get_cursor

class Element(object):
	"""Elemental data

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

	def __init__(self, element):
		self.element = element.title()
		self._meta = None

	@property
	def meta(self):
		if self._meta is None:
			from scipy.interpolate import interp1d

			self._meta = {}
			zg = get_cursor('ziegler')
			self._meta['Z'] = [str(i[2]).split(':')[0] for i in zg.execute('SELECT * FROM compounds WHERE compound=?',(self.element,))][0]
			self._meta['mass'], self._meta['density'] = [(i[1],i[2]) for i in zg.execute('SELECT * FROM weights WHERE Z=?',(self._meta['Z'],))][0]
			coeff = np.array([i[1:] for i in zg.execute('SELECT * FROM mass_coeff WHERE Z=? ORDER BY energy',(self._meta['Z'],))])
			self._meta['mass_coeff'] = interp1d(coeff[:,0], coeff[:,1], bounds_error=False, fill_value=0.0)
			self._meta['mass_coeff_en'] = interp1d(coeff[:,0], coeff[:,2], bounds_error=False, fill_value=0.0)
			db = get_cursor('decay')
			self._meta['abundances'] = {str(i[1]):i[12] for i in db.execute('SELECT * FROM chart WHERE element=? AND abundance>0.0',(self.element,))}
			self._meta['isotopes'] = sorted([i for i in self._meta['abundances']])

		return self._meta

	@property
	def Z(self):
		return self.meta['Z']

	@property
	def mass(self):
		return self.meta['mass']

	@property
	def isotopes(self):
		return self.meta['isotopes']

	@property
	def abundances(self):
		return self.meta['abundances']

	@property
	def mass_coeff(self):
		return self.meta['mass_coeff']

	@property
	def mass_coeff_en(self):
		return self.meta['mass_coeff_en']

	def attenuation(self, E, x=1.0):
		return np.exp(-self.mass_coeff(E)*self.density*x)

	def transmission(self, E, x=1.0):
		return 1.0-np.exp(-self.mass_coeff(E)*self.density*x)
	
	@property
	def density(self):
		return self.meta['density']
	
	
	


class Isotope(object):
	"""Isotopic data

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

	def __init__(self, istp):
		if istp=='1n' or istp=='1ng':
			self.element, self.A, self.isomer = 'n', 1, 'g'
		else:
			self.element = ''.join(re.findall('[A-Z]+',istp))
			if istp.startswith('nat'):
				self.A = 'nat'
			else:
				self.A = int(istp.split(self.element)[0])
			self.isomer = istp.split(self.element)[1]
		self.isotope = str(self.A)+self.element
		if self.isomer=='' and type(self.A)==int:
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
			q = [i for i in Q if i[3]==self.isomer]
			if len(q):
				q = q[0]
			elif len(Q):
				q = Q[0]
			else:
				q = [0,0,0,0,'',0,None,None,None,None,None,None,None,None,None,None,'']
			self._meta['E_level'] = q[2]
			self._meta['J_pi'] = str(q[4]) if q[4] else '?'
			self._meta['Z'] = q[6]
			self._meta['N'] = q[7]
			self._meta['stable'] = bool(q[9])
			self._meta['t_half'] = q[10]
			self._meta['unc_t_half'] = q[11]
			self._meta['abundance'] = q[12]
			self._meta['unc_abundance'] = q[13]
			if self._meta['stable'] and self._meta['abundance'] is None:
				self._meta['abundance'] = 100.0
				self._meta['unc_abundance'] = 0.0
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

		if self.stable:
			return (np.inf,0.0) if unc else np.inf
		half_conv = {'ns':1e-9,'us':1e-6,'ms':1e-3,'s':1.0,'m':60.0,'h':3600.0,'d':86400.0,'y':31557.6E3,'ky':31557.6E6}[units]
		if unc:
			return self.meta['t_half']/half_conv,self.meta['unc_t_half']/half_conv
		return self.meta['t_half']/half_conv

	def decay_const(self, units='s', unc=False):
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

		if self.stable:
			return (0.0,0.0) if unc else 0.0
		if unc:
			T2,uT2 = self.half_life(units, True)
			return np.log(2.0)/T2,np.log(2.0)*uT2/T2**2
		return np.log(2.0)/self.half_life(units)

	def optimum_units(self):
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

		opt = ['ns']
		for units in ['us','ms','s','m','h','d','y']:
			if self.half_life(units)>1.0:
				opt.append(units)
		return opt[-1]

	def abundance(self, unc=False):
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

		if unc:
			return self.meta['abundance'], self.meta['unc_abundance']
		return self.meta['abundance']

	def get_SFY(self, unc=False, closest_SFY=False):
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

		prods = []
		for (mode, product, br) in self.meta['decay_mode']:
			if product=='SFY':
				prods += [[i, br*y] for i,y in self.get_SFY(False, closest_SFY)]
			else:
				prods.append([product, br])
		return [[i, br] for i,br in prods if br>1E-8]

	def gammas(self, I_lim=[None,None], E_lim=[None,None], xrays=False, dE_511=3.5):
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

		if self.meta['gm'] is None:
			self.meta['gm'] = [[float(i[3]),float(i[4]),float(i[5]),str(i[6])] for i in self.db.execute('SELECT * FROM gammas WHERE isotope=? AND isomer=?',(self.isotope,self.isomer))]
		gammas = list(self.meta['gm'])
		for n,L in enumerate([E_lim,I_lim]):
			if L[0] is not None:
				gammas = [g for g in gammas if g[n]>=L[0]]
			if L[1] is not None:
				gammas = [g for g in gammas if g[n]<=L[1]]
		if not xrays:
			gammas = [g for g in gammas if g[3]=='']
		return {l:[g[n] for g in gammas if abs(g[0]-511.0)>=dE_511] for n,l in enumerate(['E','I','dI','notes'])}

	def electrons(self,I_lim=(None,None),E_lim=(None,None),CE_only=False,Auger_only=False):
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

		if self.meta['bm'] is None:
			self.meta['bm'] = [[float(i[3]),float(i[4]),float(i[5]),float(i[6])] for i in self.db.execute('SELECT * FROM beta_minus WHERE isotope=? AND isomer=?',(self.isotope,self.isomer))]
		betas = list(self.meta['bm'])
		for n,L in zip([3,1],[Endpoint_lim,I_lim]):
			if L[0] is not None:
				betas = [b for b in betas if b[n]>=L[0]]
			if L[1] is not None:
				betas = [b for b in betas if b[n]<=L[1]]
		return {l:[b[n] for b in betas] for n,l in enumerate(['muE','I','dI','endE'])}

	def beta_plus(self, I_lim=(None,None), Endpoint_lim=(None,None)):
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

		if self.meta['bp'] is None:
			self.meta['bp'] = [[float(i[3]),float(i[4]),float(i[5]),float(i[6])] for i in self.db.execute('SELECT * FROM beta_plus WHERE isotope=? AND isomer=?',(self.isotope,self.isomer))]
		betas = list(self.meta['bp'])
		for n,L in zip([3,1],[Endpoint_lim,I_lim]):
			if L[0] is not None:
				betas = [b for b in betas if b[n]>=L[0]]
			if L[1] is not None:
				betas = [b for b in betas if b[n]<=L[1]]
		return {l:[b[n] for b in betas] for n,l in enumerate(['muE','I','dI','endE'])}

	def alphas(self, I_lim=(None,None), E_lim=(None,None)):
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

		if self.meta['al'] is None:
			self.meta['al'] = [[float(i[3]),float(i[4]),float(i[5])] for i in self.db.execute('SELECT * FROM alphas WHERE isotope=? AND isomer=?',(self.isotope,self.isomer))]
		alphas = list(self.meta['al'])
		for n,L in enumerate([E_lim,I_lim]):
			if L[0] is not None:
				alphas = [a for a in alphas if a[n]>=L[0]]
			if L[1] is not None:
				alphas = [a for a in alphas if a[n]<=L[1]]
		return {l:[a[n] for a in alphas] for n,l in enumerate(['E','I','dI'])}

	def dose_rate(self, activity=1.0, distance=30.0, units='R/hr'):
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

		def beta2(E_MeV, m_amu):
			return 1.0-(1.0/(1.0+(E_MeV/m_amu))**2)

		def e_range(E_keV):
			dEdx = lambda b2, t: 0.17*((np.log(3.61E5*t*np.sqrt(t+2)))+(0.5*(1-b2)*(1+t**2/8-(2*t+1)*np.log(2)))-4.312)/b2
			E = np.linspace(1E-12, E_keV*1E-3, 100)
			return np.trapz(1.0/((1.0+(E*7.22/800.0))*dEdx(beta2(E, 0.5109), E/0.5109)), E)

		def pos_range(E_keV):
			dEdx = lambda b2, t: 0.17*((np.log(3.61E5*t*np.sqrt(t+2)))+(np.log(2)-(b2/24)*(23+14/(t+2)+10/(t+2)**2+4/(t+2)**3))-4.312)/b2
			E = np.linspace(1E-12, E_keV*1E-3, 100)
			return np.trapz(1.0/((1.0+(E*7.22/800.0))*dEdx(beta2(E, 0.5109), E/0.5109)), E)

		def alpha_range(E_keV):
			dEdx = lambda b2: 0.17*4.0*((np.log(1.02E6*b2/(1.0-b2))-b2)-4.312)/b2
			return np.trapz(1.0/dEdx(beta2(np.linspace(1E-9, E_keV*1E-3, 100), 3.7284E3)), np.linspace(1E-9, E_keV*1E-3, 100))

		dose = {}
		gm = self.gammas(xrays=True, dE_511=0.0)
		al = self.alphas()
		bm = self.beta_minus()
		bp = self.beta_plus()
		el = self.electrons()

		dose['gammas'] = 1.4042E-12*np.sum(np.array(gm['E'])*np.array(gm['I']))*activity/distance**2
		dose['alphas'] = 5.2087E-14*np.sum([al['I'][n]*e/alpha_range(e) for n,e in enumerate(al['E'])])*activity/distance**2
		dose['beta_minus'] = 5.2087E-14*np.sum([bm['I'][n]*e/e_range(e) for n,e in enumerate(bm['muE'])])*activity/distance**2
		dose['beta_plus'] = 5.2087E-14*np.sum([bp['I'][n]*e/pos_range(e) for n,e in enumerate(bp['muE'])])*activity/distance**2
		dose['electrons'] = 5.2087E-14*np.sum([el['I'][n]*e/e_range(e) for n,e in enumerate(el['E'])])*activity/distance**2
		return dose


