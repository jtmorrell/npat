from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import sqlite3
import npat

path = lambda db: os.path.realpath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'data',db))

DECAY_connection = None
ZIEGLER_connection = None
ENDF_connection = None
TENDL_connection = None
TENDL_rpn_connection = None
TENDL_rpp_connection = None
TENDL_rpd_connection = None
IRDFF_connection = None
CPR_connection = None

def download(db='decay',f=False):
	"""Some nice color maps.

	...

	Parameters
	----------

	Returns
	-------

	Notes
	-----

	References
	----------

	Examples
	--------

	"""

	db = db.lower()
	if db in ['all','*']:
		d = ['decay','ziegler','endf','tendl','tendl_n_rp','tendl_p_rp','tendl_d_rp','IRDFF','iaea_monitors']
	elif db in ['decay']:
		d = ['decay']
	elif db in ['ziegler']:
		d = ['ziegler']
	elif db in ['endf']:
		d = ['endf']
	elif db in ['tendl_n_rp','tendl_nrp','tendl_n','nrp','rpn']:
		d = ['tendl_n_rp']
	elif db in ['tendl_p_rp','tendl_prp','tendl_p','prp','rpp']:
		d = ['tendl_p_rp']
	elif db in ['tendl_d_rp','tendl_drp','tendl_d','drp','rpd']:
		d = ['tendl_d_rp']
	elif db in ['irdff']:
		d = ['IRDFF']
	elif db in ['iaea','iaea-cpr','iaea-monitor','cpr','iaea_cpr','iaea_monitor','medical','iaea-medical','iaea_medical']:
		d = ['iaea_monitors']
	else:
		print('db={} not recognized.'.format(db))
		return

	addr = {'decay':'wwd6b1gk2ge5tgt', 'endf':'tkndjqs036piojm', 'tendl':'zkoi6t2jicc9yqs', 'tendl_d_rp':'x2vfjr7uv7ffex5', 'tendl_n_rp':'n0jjc0dv61j9of9',
				'tendl_p_rp':'ib2a5lrhiwkcro5', 'ziegler':'kq07684wtp890v5','iaea_monitors':'lzn8zs6y8zu3v0s','IRDFF':'34sgcvt8n57b0aw'}
	
	if not os.path.isdir(path('')):
		os.mkdir(path(''))

	try:
		import urllib2
	except:
		import urllib.request as urllib2

	for i in d:
		fnm = i+'.db'
		if (not os.path.isfile(path(fnm))) or f:
			
			try:
				print('Downloading {}'.format(fnm))
				with open(path(fnm),'wb') as f:
					f.write(urllib2.urlopen('https://www.dropbox.com/s/{0}/{1}?dl=1'.format(addr[i],fnm)).read())
			except Exception as e:
					print(e)
		else:
			print("{0}.db already installed. Run npat.download(db='{0}',f=True) to override.".format(i))



def get_connection(db='decay'):
	"""Some nice color maps.

	...

	Parameters
	----------

	Returns
	-------

	Notes
	-----

	References
	----------

	Examples
	--------

	"""

	db = db.lower().replace('.db','')

	if db in ['decay']:
		global DECAY_connection
		if DECAY_connection is None:
			DECAY_connection = sqlite3.connect(path('decay.db'))

		return DECAY_connection

	elif db in ['ziegler']:
		global ZIEGLER_connection
		if ZIEGLER_connection is None:
			ZIEGLER_connection = sqlite3.connect(path('ziegler.db'))

		return ZIEGLER_connection

	elif db in ['endf']:
		global ENDF_connection
		if ENDF_connection is None:
			ENDF_connection = sqlite3.connect(path('endf.db'))

		return ENDF_connection

	elif db in ['tendl']:
		global TENDL_connection
		if TENDL_connection is None:
			TENDL_connection = sqlite3.connect(path('tendl.db'))

		return TENDL_connection

	elif db in ['tendl_n_rp','tendl_nrp','tendl_n','nrp','rpn']:
		global TENDL_rpn_connection
		if TENDL_rpn_connection is None:
			TENDL_rpn_connection = sqlite3.connect(path('tendl_n_rp.db'))

		return TENDL_rpn_connection

	elif db in ['tendl_p_rp','tendl_prp','tendl_p','prp','rpp']:
		global TENDL_rpp_connection
		if TENDL_rpp_connection is None:
			TENDL_rpp_connection = sqlite3.connect(path('tendl_p_rp.db'))

		return TENDL_rpp_connection

	elif db in ['tendl_d_rp','tendl_drp','tendl_d','drp','rpd']:
		global TENDL_rpd_connection
		if TENDL_rpd_connection is None:
			TENDL_rpd_connection = sqlite3.connect(path('tendl_d_rp.db'))

		return TENDL_rpd_connection	

	elif db in ['irdff']:
		global IRDFF_connection
		if IRDFF_connection is None:
			IRDFF_connection = sqlite3.connect(path('IRDFF.db'))

		return IRDFF_connection

	elif db in ['iaea','iaea-cpr','iaea-monitor','cpr','iaea_cpr','iaea_monitor','medical','iaea-medical','iaea_medical']:
		global CPR_connection
		if CPR_connection is None:
			CPR_connection = sqlite3.connect(path('iaea_monitors.db'))

		return CPR_connection	


def get_cursor(db='decay'):
	"""Some nice color maps.

	...

	Parameters
	----------

	Returns
	-------

	Notes
	-----

	References
	----------

	Examples
	--------

	"""
	
	conn = get_connection(db)
	return conn.cursor()
