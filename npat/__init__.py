# import matplotlib
# matplotlib.use('Agg')
from .dbmgr import *
from .spectroscopy import *
from .isotope import *
from .decay_chain import *
from .plotter import *
from .reaction import *
from .irradiation import *

import sys
if sys.version_info[0]<3:
	import copy_reg
	import types

	def _pickle_method(m):
		if m.im_self is None:
			return getattr, (m.im_class, m.im_func.func_name)
		else:
			return getattr, (m.im_self, m.im_func.func_name)

	copy_reg.pickle(types.MethodType, _pickle_method)


__version__ = '0.2.8'
__all__ = ['get_cursor', 'get_connection', 'colors',
			'set_style', 'Spectrum', 'Calibration', 
			'Isotope', 'DecayChain', 'Reaction',
			'Ziegler', 'Library']