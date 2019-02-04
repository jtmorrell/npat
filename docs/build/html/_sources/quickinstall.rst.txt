.. _quickinstall:

==================
Installation Guide
==================

NPAT is available through the `Python Package index`_, which allows installation using Python's standard command line utility `pip`_.  Assuming Python and Pip are installed already, you can install NPAT with the command::

	pip install --user npat


or::

	python -m pip install --user npat


If you'd like to install the most recent (unreleased) version of NPAT, you can clone the NPAT `github`_, and install by either running ``python -m setup.py install --user`` or ``make install``.

.. _Python Package index: https://pypi.org/
.. _pip: https://pip.pypa.io/en/stable
.. _github: https://github.com/jtmorrell/npat


Troubleshooting
---------------

NPAT should download about 1GB of nuclear data during setup, however if this download fails for some reason the databases can be downloaded by downloading the `setup`_ file and running it with the command::

	python setup.py install --user

.. _setup: https://github.com/jtmorrell/npat/blob/master/setup.py