from .dbmgr import *
from .spectroscopy import *
from .isotope import *
from .decay_chain import *
from .plotter import *
from .reaction import *
from .irradiation import *

import copy_reg
import types

def _pickle_method(m):
    if m.im_self is None:
        return getattr, (m.im_class, m.im_func.func_name)
    else:
        return getattr, (m.im_self, m.im_func.func_name)

copy_reg.pickle(types.MethodType, _pickle_method)

# def _pickle_method(method):
# 	func_name = method.im_func.__name__
# 	obj = method.im_self
# 	cls = method.im_class
# 	return _unpickle_method, (func_name, obj, cls)

# def _unpickle_method(func_name, obj, cls):
# 	for cls in cls.mro():
# 		try:
# 			func = cls.__dict__[func_name]
# 		except KeyError:
# 			pass
# 		else:
# 			break
# 		return func.__get__(obj, cls)

# copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

__version__ = '0.2.5'
__all__ = ['get_cursor', 'get_connection', 'colors',
			'set_style', 'Spectrum', 'Calibration', 
			'Isotope', 'DecayChain', 'Reaction',
			'Ziegler', 'Irradiation', 'Sample',
			'Library']