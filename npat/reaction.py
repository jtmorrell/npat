from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import re
import numpy as np
from scipy.interpolate import interp1d

from .dbmgr import get_cursor
from .isotope import Isotope
from .plotter import colors
from .plotter import _init_plot
from .plotter import _close_plot

class Library(object):
	"""Interface with cross section library

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

	def __init__(self, name='tendl'):
		name = name.lower()
		if name in ['endf']:
			self.db_name = 'endf'
		elif name in ['tendl']:
			self.db_name = 'tendl'
		elif name in ['tendl_n_rp','tendl_nrp','tendl_n','nrp','rpn']:
			self.db_name = 'tendl_n_rp'
		elif name in ['tendl_p_rp','tendl_prp','tendl_p','prp','rpp']:
			self.db_name = 'tendl_p_rp'
		elif name in ['tendl_d_rp','tendl_drp','tendl_d','drp','rpd']:
			self.db_name = 'tendl_d_rp'
		elif name in ['irdff']:
			self.db_name = 'irdff'
		elif name in ['iaea','iaea-cpr','iaea-monitor','cpr','iaea_cpr','iaea_monitor','medical','iaea-medical','iaea_medical']:
			self.db_name = 'iaea_medical'
		else:
			raise ValueError('Library {} not recognized.'.format(name))
		self.db = get_cursor(self.db_name)
		self.name = {'endf':'ENDF/B-VII.1','tendl':'TENDL-2015','irdff':'IRDFF-II','iaea':'IAEA CP-Reference (2017)'}[self.db_name.split('_')[0]]
	
	def __str__(self):
		return self.name

	def check(self, target=None, incident=None, outgoing=None, product=None):
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

		return len(self.search(target, incident, outgoing, product))==1

	def search(self, target=None, incident=None, outgoing=None, product=None, _label=False):
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

		ss = 'SELECT * FROM all_reactions'
		if incident is not None:
			incident = incident.lower()
		if outgoing is not None:
			outgoing = outgoing.lower()

		if self.db_name in ['endf','tendl','irdff']:
			if incident is not None:
				if incident!='n':
					return []
			q = [(i+'%' if not n%2 else i) for n,i in enumerate([target, outgoing, product]) if i]
			ss += ' WHERE ' if len(q) else ''
			ss += ' AND '.join([i for i in [('target LIKE ?' if target else ''),('outgoing=?' if outgoing else ''),('product LIKE ?' if product else '')] if i])
			reacs = [list(map(str, i)) for i in self.db.execute(ss, tuple(q))]
			fmt = '{0}(n,{1}){2}'

		elif self.db_name in ['tendl_n_rp','tendl_p_rp','tendl_d_rp']:
			if incident is not None:
				if incident!=self.db_name.split('_')[1]:
					return []
			if product:
				if 'm' not in product and 'g' not in product:
					product += 'g'
					print('WARNING: Product isomeric state not specified, ground state assumed.')
			q = [i+'%' for i in [target, product] if i]
			ss += ' WHERE ' if len(q) else ''
			ss += ' AND '.join([i for i in [('target LIKE ?' if target else ''),('product LIKE ?' if product else '')] if i])
			reacs = [list(map(str, i)) for i in self.db.execute(ss, tuple(q))]
			fmt = '{0}('+self.db_name.split('_')[1]+',x){1}'

		elif self.db_name=='iaea_medical':
			q = [(i+'%' if n in [0,3] else i) for n,i in enumerate([target, incident, outgoing, product]) if i]
			ss += ' WHERE ' if len(q) else ''
			ss += ' AND '.join([i for i in [('target LIKE ?' if target else ''),('incident=?' if incident else ''),('outgoing=?' if outgoing else ''),('product LIKE ?' if product else '')] if i])
			reacs = [list(map(str, i)) for i in self.db.execute(ss, tuple(q))]
			fmt = '{0}({1},{2}){3}'

		if target:
			reacs = [i for i in reacs if i[0].lower()==target.lower()]
		if _label:
			return [i[-1] for i in reacs]
		return [fmt.format(*i) for i in reacs]

	def query(self, target=None, incident=None, outgoing=None, product=None):
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

		labels = self.search(target, incident, outgoing, product, _label=True)
		if not len(labels)==1:
			raise ValueError('{0}({1},{2}){3}'.format(target, incident, outgoing, product)+' is not a unique reaction.')
		if not target:
			raise ValueError('Target Must be specified.')
		if self.db_name in ['endf','tendl','tendl_n_rp','tendl_p_rp','tendl_d_rp']:
			table = ''.join(re.findall('[A-Z]+', target))+'_'+''.join(re.findall('[0-9]+', target))+('m' if 'm' in target else '')
			return np.array(list(self.db.execute('SELECT energy,{0} FROM {1}'.format(labels[0], table))))*(np.array([1E-6, 1E3]) if self.db_name=='endf' else np.ones(2))

		elif self.db_name=='irdff':
			return np.array(list(self.db.execute('SELECT * FROM {}'.format(labels[0]))))*np.array([1E-6, 1E3, 1E3])

		elif self.db_name=='iaea_medical':
			if incident is None:
				raise ValueError('Incident particle must be specified.')
			table = {'n':'neutron','p':'proton','d':'deuteron','h':'helion','a':'alpha','g':'gamma'}[incident]
			return np.array(list(self.db.execute('SELECT energy,cross_section,unc_cross_section FROM {0} WHERE target LIKE ? AND product=?'.format(table),('%'+target+'%', labels[0]))))







class Reaction(object):
	"""Retrieves cross sections for a given nuclear reaction

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

	def __init__(self, reaction_name, library='best'):
		self.target, p = tuple(reaction_name.split('('))
		p, self.product = tuple(p.split(')'))
		self.incident, self.outgoing = tuple(p.split(','))
		self.incident, self.outgoing = self.incident.lower(), self.outgoing.lower()
		self._rx = [self.target, self.incident, self.outgoing, self.product]
		self.name = reaction_name

		if library.lower()=='best':
			if self.incident=='n':
				for lb in ['irdff','endf','iaea','tendl','tendl_n_rp']:
					self.library = Library(lb)
					if lb=='tendl_n_rp':
						self._check(True)
					elif self._check():
						break
			elif self.incident in ['p','d']:
				for lb in ['iaea','tendl_'+self.incident+'_rp']:
					self.library = Library(lb)
					if lb=='tendl_d_rp':
						self._check(True)
					elif self._check():
						break
			else:
				self.library = Library('iaea')
				self._check(True)
		else:
			self.library = Library(library)
			self._check(True)

		self.name = self.library.search(*self._rx)[0]
		q = self.library.query(*self._rx)
		self.eng = q[:,0]
		self.xs = q[:,1]
		if q.shape[1]==3:
			self.unc_xs = q[:,2]
		else:
			self.unc_xs = np.zeros(len(self.xs))
		self._interp = None
		self._interp_unc = None
		self._tex = None

	def _check(self, err=False):
		c = self.library.check(*self._rx)
		if err and not c:
			raise ValueError('Reaction '+self.name+' not found or not unique.')
		return c

	def __str__(self):
		return self.name

	@property
	def TeX(self):
		if self._tex is None:
			if 'nat' not in self.target:
				target = Isotope(self.target).TeX
			else:
				target = r'$^{nat}$'+self.target[3:].title()
			product = Isotope(self.product).TeX if self.product else ''
			self._tex = '{0}({1},{2}){3}'.format(target, self.incident, self.outgoing, product)
		return self._tex
	
	@property
	def interp(self):
		if self._interp is None:
			self._interp = interp1d(self.eng, self.xs, bounds_error=False, fill_value=0.0)
		return self._interp

	@property
	def interp_unc(self):
		if self._interp_unc is None:
			self._interp_unc = interp1d(self.eng, self.unc_xs, bounds_error=False, fill_value=0.0)
		return self._interp_unc
	
	def integrate(self, energy, flux, unc=False):
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

		# Trapezoidal Riemann sum
		E = np.asarray(energy)
		phisig = np.asarray(flux)*self.interp(E)
		if unc:
			unc_phisig = np.asarray(flux)*self.interp_unc(E)
			return np.sum(0.5*(E[1:]-E[:-1])*(phisig[:-1]+phisig[1:])), np.sum(0.5*(E[1:]-E[:-1])*(unc_phisig[:-1]+unc_phisig[1:]))
		return np.sum(0.5*(E[1:]-E[:-1])*(phisig[:-1]+phisig[1:]))

	def average(self, energy, flux, unc=False):
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

		E, phi = np.asarray(energy), np.asarray(flux)
		phisig = phi*self.interp(E)
		dE = E[1:]-E[:-1]
		if unc:
			unc_phisig = np.asarray(flux)*self.interp_unc(E)
			return np.sum(0.5*dE*(phisig[:-1]+phisig[1:]))/np.sum(0.5*dE*(phi[:-1]+phi[1:])), np.sum(0.5*dE*(unc_phisig[:-1]+unc_phisig[1:]))/np.sum(0.5*dE*(phi[:-1]+phi[1:]))
		return np.sum(0.5*dE*(phisig[:-1]+phisig[1:]))/np.sum(0.5*dE*(phi[:-1]+phi[1:]))

	def plot(self, label=None, title=True, E_lim=None, **kwargs):
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

		f, ax = _init_plot(**kwargs)

		if title:
			ax.set_title(self.TeX)

		if label is not None:
			if label.lower() in ['both','library','reaction']:
				label = {'both':'{0}\n({1})'.format(self.TeX, self.library.name),'library':self.library.name,'reaction':self.TeX}[label.lower()]

		lh = [0, len(self.eng)] if E_lim is None else np.where((self.eng>=E_lim[0])&(self.eng<=E_lim[1]))[0][[0,-1]]
		l,h = lh[0], lh[1]
		
		line, = ax.plot(self.eng[l:h], self.xs[l:h], label=label)
		if np.any(self.unc_xs[l:h]>0):
			ax.fill_between(self.eng[l:h], self.xs[l:h]+self.unc_xs[l:h], self.xs[l:h]-self.unc_xs[l:h], facecolor=line.get_color(), alpha=0.5)

		ax.set_xlabel('Incident Energy (MeV)')
		ax.set_ylabel('Cross Section (mb)')

		if label:
			ax.legend(loc=0)

		return _close_plot(f, ax, default_log=False, **kwargs)


