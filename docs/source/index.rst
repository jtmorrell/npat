
================================
Welcome to NPAT's documentation!
================================

NPAT (nuclear physics analysis tools) is a python toolkit to aid in the analysis of experimental nuclear data.  

The primary application for NPAT is activation analysis, with specific utilities developed for stacked-target, charged-particle activation analysis.
However, the convenient python interface to nuclear decay data and a range of cross section libraries makes NPAT more generally useful in the nuclear sciences.

--------
Features
--------

NPAT features the following classes to aid in data analysis:

* Spectrum - Peak fitting for HPGe detector data.
* Calibration - Energy & efficiency calibration tool (for HPGe detectors).
* Ziegler - Stacked-target energy loss characterization.
* DecayChain - General purpose Bateman equation solver.
* Isotope - Isotopic and decay data.
* Reaction - Cross sections from multiple libraries.
* Library - Tool for searching and retrieving cross sections from multiple libraries.

--------
Contents
--------

.. toctree::
   :maxdepth: 1
   
   quickinstall
   methods
   usersguide/index
   api/index
   license

