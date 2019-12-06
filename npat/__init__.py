"""
npat
=====

Provides
  1. Classes for plotting and fitting gamma ray spectra
  2. Access to data from evaluated nuclear reaction libraries
  3. Parsing capabilities for MVME listfiles
  4. Charged particle stopping power calculations
  5. Generalized Bateman equation solver
  6. Atomic and nuclear structure/decay data

How to use the documentation
----------------------------
Documentation is available in two forms: docstrings provided
with the code, and a users guide located on
`the npat homepage <https://jtmorrell.github.io/npat/build/html/index.html>`_.

Code snippets are indicated by three greater-than signs::

  >>> sp = npat.Spectrum('spectrum.Spe')
  >>> sp.plot()

Use the built-in ``help`` function to view a function or class's docstring::

  >>> help(npat.Spectrum)




"""

from .dbmgr import *
from .spectroscopy import *
from .isotope import *
from .decay_chain import *
from .plotter import *
from .reaction import *
from .irradiation import *
from .listfiles import *

__version__ = '0.3.10'
__all__ = ['get_cursor', 'get_connection', 'colors',
			'set_style', 'Spectrum', 'Calibration', 
			'Isotope', 'DecayChain', 'Reaction',
			'Ziegler', 'Library', 'Element', 
			'download', 'MVME']