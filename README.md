# NPAT

NPAT (nuclear physics analysis tools) is a python toolkit developed to aid in the analysis of experimental nuclear data.

The primary application for NPAT is activation analysis, with specific utilities developed for stacked-target, charged-particle activation analysis.
However, the convenient python interface to nuclear decay data and a range of cross section libraries makes NPAT more generally useful in the nuclear sciences.

## Features

NPAT features the following classes to aid in data analysis:

* Spectrum - Peak fitting for HPGe detector data.
* Calibration - Energy & efficiency calibration tool (for HPGe detectors).
* Ziegler - Stacked-target energy loss characterization.
* DecayChain - General purpose Bateman equation solver.
* Isotope - Isotopic and decay data.
* Reaction - Cross sections from multiple libraries.
* Library - Tool for searching and retrieving cross sections from multiple libraries.

More detail can be found in the [User's Guide](https://jtmorrell.github.io/npat/build/html/usersguide/index.html) and the [API](https://jtmorrell.github.io/npat/build/html/api/index.html).

## Quick Install

NPAT is available through the [Python Package index](https://pypi.org/), which allows installation using Python's standard command line utility [pip](https://pip.pypa.io/en/stable).  Assuming Python and Pip are installed already, you can install NPAT with the command:

```
pip install --user npat
```

or:

```
python -m pip install --user npat
```

Detailed installation instructions and troubleshooting can be found on the NPAT [documentation site](https://jtmorrell.github.io/npat/build/html/quickinstall.html). 



