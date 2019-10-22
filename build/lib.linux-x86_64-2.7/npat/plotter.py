from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import matplotlib.pyplot as plt
from cycler import cycler


def colors(style='default', shade='dark', aslist=False):
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
	default = {'light':{'k':'#34495e','gy':'#95a5a6','r':'#e74c3c','b':'#3498db','g':'#2ecc71',
	'aq':'#1abc9c','o':'#e67e22','y':'#f1c40f','p':'#9b59b6','w':'#ecf0f1'},
	'dark':{'k':'#2c3e50','gy':'#7f8c8d','r':'#c0392b','b':'#2980b9','g':'#27ae60',
	'aq':'#16a085','o':'#d35400','y':'#f39c12','p':'#8e44ad','w':'#bdc3c7'}}

	american = {'light':{'k':'#636e72','gy':'#dfe6e9','r':'#ff7675','b':'#74b9ff','g':'#55efc4',
	'aq':'#81ecec','o':'#fab1a0','y':'#ffeaa7','p':'#a29bfe','w':'#ecf0f1'},
	'dark':{'k':'#2d3436','gy':'#b2bec3','r':'#d63031','b':'#0984e3','g':'#00b894',
	'aq':'#00cec9','o':'#e17055','y':'#fdcb6e','p':'#6c5ce7','w':'#bdc3c7'}}

	aussie = {'light':{'k':'#95afc0','gy':'#dff9fb','r':'#ff7979','b':'#686de0','g':'#badc58',
	'aq':'#7ed6df','o':'#ffbe76','y':'#f6e58d','p':'#e056fd','w':'#ecf0f1'},
	'dark':{'k':'#535c68','gy':'#c7ecee','r':'#eb4d4b','b':'#4834d4','g':'#6ab04c',
	'aq':'#22a6b3','o':'#f0932b','y':'#f9ca24','p':'#be2edd','w':'#bdc3c7'}}

	british = {'light':{'k':'#353b48','gy':'#7f8fa6','r':'#e84118','b':'#273c75','g':'#4cd137',
	'aq':'#1abc9c','o':'#e67e22','y':'#fbc531','p':'#9c88ff','w':'#f5f6fa'},
	'dark':{'k':'#2f3640','gy':'#718093','r':'#c23616','b':'#192a56','g':'#44bd32',
	'aq':'#16a085','o':'#d35400','y':'#e1b12c','p':'#8c7ae6','w':'#dcdde1'}}

	canadian = {'light':{'k':'#576574','gy':'#c8d6e5','r':'#ff6b6b','b':'#54a0ff','g':'#1dd1a1',
	'aq':'#00d2d3','o':'#feca57','y':'#f1c40f','p':'#5f27cd','w':'#ecf0f1'},
	'dark':{'k':'#222f3e','gy':'#8395a7','r':'#ee5253','b':'#2e86de','g':'#10ac84',
	'aq':'#01a3a4','o':'#ff9f43','y':'#f39c12','p':'#341f97','w':'#bdc3c7'}}

	chinese = {'light':{'k':'#57606f','gy':'#a4b0be','r':'#ff6b81','b':'#70a1ff','g':'#7bed9f',
	'aq':'#1abc9c','o':'#ff7f50','y':'#eccc68','p':'#5352ed','w':'#ffffff'},
	'dark':{'k':'#2f3542','gy':'#747d8c','r':'#ff4757','b':'#1e90ff','g':'#2ed573',
	'aq':'#16a085','o':'#ff6348','y':'#ffa502','p':'#3742fa','w':'#f1f2f6'}}

	german = {'light':{'k':'#778ca3','gy':'#d1d8e0','r':'#fc5c65','b':'#45aaf2','g':'#26de81',
	'aq':'#2bcbba','o':'#fd9644','y':'#fed330','p':'#a55eea','w':'#ecf0f1'},
	'dark':{'k':'#4b6584','gy':'#a5b1c2','r':'#eb3b5a','b':'#2d98da','g':'#20bf6b',
	'aq':'#0fb9b1','o':'#fa8231','y':'#f7b731','p':'#8854d0','w':'#bdc3c7'}}

	spanish = {'light':{'k':'#34495e','gy':'#d1ccc0','r':'#ff5252','b':'#34ace0','g':'#33d9b2',
	'aq':'#1abc9c','o':'#ff793f','y':'#ffb142','p':'#706fd3','w':'#f7f1e3'},
	'dark':{'k':'#2c3e50','gy':'#84817a','r':'#b33939','b':'#227093','g':'#218c74',
	'aq':'#16a085','o':'#cd6133','y':'#cc8e35','p':'#474787','w':'#aaa69d'}}

	swedish = {'light':{'k':'#485460','gy':'#d2dae2','r':'#ff5e57','b':'#4bcffa','g':'#0be881',
	'aq':'#34e7e4','o':'#ffc048','y':'#ffdd59','p':'#575fcf','w':'#ecf0f1'},
	'dark':{'k':'#1e272e','gy':'#808e9b','r':'#ff3f34','b':'#0fbcf9','g':'#05c46b',
	'aq':'#00d8d6','o':'#ffa801','y':'#ffd32a','p':'#3c40c6','w':'#bdc3c7'}}

	cm = {'default':default,'american':american,'aussie':aussie,'british':british,'canadian':canadian,'chinese':chinese,'german':german,'spanish':spanish,'swedish':swedish}
	if aslist:
		return [cm[style][shade][c] for c in ['k','gy','r','b','g','o','aq','y','p']]
	return cm[style][shade]


def set_style(sty='show'):
	"""Preset styles for various purposes.

	...

	Parameters
	----------

	Returns
	-------

	"""
	cm = colors()
	plt.rcParams['font.family']='sans-serif'
	plt.rcParams['axes.prop_cycle'] = cycler(color=[cm['k'],cm['r'],cm['b'],cm['gy'],cm['g'],cm['p'],cm['o'],cm['aq'],cm['y']])
	plt.rcParams['figure.autolayout'] = 'True'
	plt.rcParams['xtick.minor.visible']='True'
	plt.rcParams['ytick.minor.visible']='True'
	plt.rcParams['axes.titlepad']='8'
	if sty=='show':
		plt.rcParams['font.size']='14'
		plt.rcParams['xtick.major.pad']='8'
		plt.rcParams['ytick.major.pad']='8'
		plt.rcParams['legend.fontsize']='12'
		plt.rcParams['lines.markersize']='4.5'
		plt.rcParams['errorbar.capsize']='5.0'
		plt.rcParams['lines.linewidth'] = '1.8'
	elif sty=='paper':
		plt.rcParams['font.size']='14'
		plt.rcParams['xtick.major.pad']='8'
		plt.rcParams['ytick.major.pad']='8'
		plt.rcParams['legend.fontsize']='12'
		plt.rcParams['lines.markersize']='3.5'
		plt.rcParams['errorbar.capsize']='4.0'
		plt.rcParams['lines.linewidth'] = '2.4'
	elif sty=='poster':
		plt.rcParams['font.size']='18'
		plt.rcParams['font.weight']='bold'
		plt.rcParams['axes.titleweight']='bold'
		plt.rcParams['axes.titlepad']='14'
		plt.rcParams['xtick.major.pad']='6'
		plt.rcParams['ytick.major.pad']='6'
		plt.rcParams['xtick.major.size']='5'
		plt.rcParams['xtick.minor.size']='3.5'
		plt.rcParams['xtick.major.width']='1.5'
		plt.rcParams['xtick.minor.width']='1.2'
		plt.rcParams['ytick.major.size']='5'
		plt.rcParams['ytick.minor.size']='3.5'
		plt.rcParams['ytick.major.width']='1.5'
		plt.rcParams['ytick.minor.width']='1.2'
		plt.rcParams['axes.linewidth']='1.2'
		plt.rcParams['legend.fontsize']='14'
		plt.rcParams['lines.markersize']='4.5'
		plt.rcParams['errorbar.capsize']='5.0'
		plt.rcParams['lines.linewidth'] = '2.8'
	elif sty=='presentation':
		plt.rcParams['font.size']='14'
		plt.rcParams['xtick.major.pad']='8'
		plt.rcParams['ytick.major.pad']='8'
		plt.rcParams['legend.fontsize']='12'
		plt.rcParams['lines.markersize']='3.5'
		plt.rcParams['errorbar.capsize']='4.0'
		plt.rcParams['lines.linewidth'] = '1.8'


def _init_plot(**kwargs):
	f, ax = None, None
	if 'f' in kwargs and 'ax' in kwargs:
		f, ax = kwargs['f'], kwargs['ax']
	
	N = kwargs['N_plots'] if 'N_plots' in kwargs else 1

	if f is None or ax is None:
		if 'figsize' in kwargs:
			f, ax = plt.subplots(1, N, figsize=kwargs['figsize'])
		else:
			f, ax = plt.subplots(1, N)

	if 'style' in kwargs:
		set_style(kwargs['style'])

	return f, ax


def _close_plot(fig, axis, **kwargs):
	f, ax = fig, axis

	if 'default_log' in kwargs:
		if kwargs['default_log']:
			ax.set_yscale('log')
	else:
		ax.set_yscale('log')

	if 'scale' in kwargs:
		s = kwargs['scale'].lower()
		if s in ['log','logy']:
			ax.set_yscale('log')
		elif s in ['lin','liny']:
			ax.set_yscale('linear')
		elif s=='logx':
			ax.set_xscale('log')
		elif s=='linx':
			ax.set_xscale('linear')
		elif s in ['linlog','loglog','loglin','linlin']:
			ax.set_yscale(s[3:].replace('lin','linear'))
			ax.set_xscale(s[:3].replace('lin','linear'))

	if 'logscale' in kwargs:
		if kwargs['logscale']:
			ax.set_yscale('log')
		else:
			ax.set_yscale('linear')

	if 'logx' in kwargs:
		if kwargs['logx']:
			ax.set_xscale('log')
		else:
			ax.set_xscale('linear')

	if 'logy' in kwargs:
		if kwargs['logy']:
			ax.set_yscale('log')
		else:
			ax.set_yscale('linear')

	f.tight_layout()

	if 'saveas' in kwargs:
		if kwargs['saveas'] is not None:
			f.savefig(kwargs['saveas'])

	if 'show' in kwargs:
		if kwargs['show']:
			plt.show()
	else:
		plt.show()

	if 'f' in kwargs and 'ax' in kwargs:
		return f, ax

	plt.close()

set_style()
