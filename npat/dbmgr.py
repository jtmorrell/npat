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

def get_connection(db='decay'):
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
	conn = get_connection(db)
	return conn.cursor()
