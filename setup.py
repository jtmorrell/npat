from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from setuptools import setup, find_packages
import os
from distutils.core import setup
from distutils.command.install import install as _install


def _post_install(dir):
	path = lambda db: os.path.realpath(os.path.join(dir,'npat','data',db))

	if not os.path.isdir(path('')):
		os.mkdir(path(''))

	for fl in ['wwd6b1gk2ge5tgt/decay.db','tkndjqs036piojm/endf.db','zkoi6t2jicc9yqs/tendl.db','x2vfjr7uv7ffex5/tendl_d_rp.db','n0jjc0dv61j9of9/tendl_n_rp.db',
				'ib2a5lrhiwkcro5/tendl_p_rp.db','kq07684wtp890v5/ziegler.db','lzn8zs6y8zu3v0s/iaea_monitors.db','34sgcvt8n57b0aw/IRDFF.db']:
		fnm = fl.split('/')[1]
		if not os.path.isfile(path(fnm)):
			try:
				import urllib2
			except:
				import urllib.request as urllib2
			try:
				print('Downloading {}'.format(fnm))
				with open(path(fnm),'wb') as f:
					f.write(urllib2.urlopen('https://www.dropbox.com/s/{}?dl=1'.format(fl)).read())
			except:
				print('ERROR: Unable to download {}, please reinstall NPAT'.format(fnm))


class install(_install):
	def run(self):
		_install.run(self)
		self.execute(_post_install, (self.install_lib,),
					 msg="Downloading nuclear data files...")

setup(name='npat',
	  version='0.2.5',
	  description='Nuclear Physics Analysis Tools (NPAT) is a library written in python to assist in analysis of the physics of nuclear reactions and spectroscopy.',
	  url='https://github.com/jtmorrell/npat',
	  author='Jonathan Morrell',
	  author_email='jmorrell@berkeley.edu',
	  license='MIT',
	  packages=find_packages(),
	  include_package_data=True,
	  cmdclass={'install': install})
# , install_requires=['numpy>=1.11', 'matplotlib>=1.3', 'scipy>=0.19.1'])