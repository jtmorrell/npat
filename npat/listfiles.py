from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, re, zipfile, json
import numpy as np
import datetime as dtm

from npat import Spectrum


class MVME(object):
	"""Spectrum is a class for fitting and plotting
	gamma ray data from High-Purity Germanium (HPGe)
	detectors.

	...


	Parameters
	----------
	filename : str
		Path to .Spe file.
	db : str, optional
		Path to sqlite database

	Attributes
	----------
	fits : list
		List of peaks (PeakFit class) found in spectrum fit.
	meta : dict
		Metadata about spectrum. 
	hist : np.ndarray
		1D histogram of pulse amplitudes (energy).

	Methods
	-------
	plot(show=True, fit=True, saveas=None, zoom=None, logscale=True, 
			grayscale=False, labels=False, square_fig=False)
		Plots the spectrum. Various options for plotting.
	save(*saveas)
		Functionality depends on filetype.
	summarize(printout=True, saveas=None)
		Prints and/or saves summary of peaks/isotopes.

	Notes
	-----

	References
	----------

	Examples
	--------

	"""

	def __init__(self, filename):
		self.zipfilename = filename
		self.zp = zipfile.ZipFile(filename, 'r')
		_path, _fnm = os.path.split(filename)
		print('Reading {}'.format(_fnm))
		self._meta = {}
		self.meta = {'tdc_resolution':24E-12, 'time_bins':1, 'time_bin_length':None}

		self._parsed = False

	@property
	def meta(self):
		return self._meta

	@meta.setter
	def meta(self, meta_dict):
		for nm in meta_dict:
			self._meta[nm] = meta_dict[nm]

	def _parse_log(self, fl):
		self._real_time = 0.0
		for ln in fl.splitlines():
			if re.search('readout starting', ln):
				self._start_time = dtm.datetime.strptime(ln.split(' ')[-1], '%Y-%m-%dT%H:%M:%S')
			elif re.search('readout stopped', ln):
				self._real_time = float((dtm.datetime.strptime(ln.split(' ')[-1], '%Y-%m-%dT%H:%M:%S')-self._start_time).total_seconds())

	def _parse_analysis(self, fl):
		fl = json.loads(fl)
		modules = [str(c['moduleName']) for c in fl['AnalysisNG']['properties']['ModuleProperties'] if c['moduleTypeName'].startswith('mdpp16')]
		mod_id = [str(c['moduleId']) for c in fl['AnalysisNG']['properties']['ModuleProperties'] if c['moduleTypeName'].startswith('mdpp16')]

		self._cb = {i:[] for i in modules}
		for c in fl['AnalysisNG']['operators']:
			if c['class'].endswith('CalibrationMinMax'):
				if c['name'].endswith('.amplitude'):
					cal_id = None
					for cn in fl['AnalysisNG']['connections']:
						if cn['dstId']==c['id']:
							cal_id = cn['srcId']
					m_id = None
					if cal_id is not None:
						for cn in fl['AnalysisNG']['sources']:
							if cn['id']==cal_id:
								m_id = cn['moduleId']
					if m_id is not None:
						if m_id in mod_id:
							mod = modules[mod_id.index(m_id)]
							cl = [[0.0, 1.0] for i in range(16)]
							for j,cb in enumerate(c['data']['calibrations']):
								mn, mx = float(cb['unitMin']), float(cb['unitMax'])
								cl[j] = [mn, (mx-mn)/float(2**16)]
							self._cb[mod] = cl

	def _parse_header(self, head):
		self.modules = [str(m['name']) for m in head['DAQConfig']['events'][0]['modules']]

	def _default_fmap(self, amplitude, millis, board, channel, overflow, pileup):
		if self.spectra is None:
			self.spectra = [[[None for m in range(len(self._time_bins)-1)] for i in range(16)] for n in range(self._N_boards)]
		for bd in range(self._N_boards):
			for ch in range(16):
				for tm in range(len(self._time_bins)-1):
					start, stop = self._time_bins[tm], self._time_bins[tm+1]


					ix = np.where((channel==ch)&(board==bd)&(millis>=start*1E3)&(millis<stop*1E3))[0]
					if not len(ix):
						continue

					dead_time = 300E-9*len(amplitude[ix])
					pu_ix = np.where(pileup[ix])[0]
					if len(pu_ix):
						if pu_ix[-1]==len(ix)-1:
							pu_ix = pu_ix[:-1]
						diff = (millis[ix][pu_ix+1]-millis[ix][pu_ix])*1E3
						dead_time += np.sum(np.where((diff>0)&(diff<1), diff, 0.0))

					hist = np.histogram(amplitude[ix][np.where((pileup[ix]==0)&(overflow[ix]==0)&(amplitude[ix]>0))], bins=np.arange(-0.5,2**16+0.5,1))[0]
					if np.any(hist[500:]>0):
						if self.spectra[bd][ch][tm] is None:
							self.spectra[bd][ch][tm] = Spectrum()
							st = self._start_time + dtm.timedelta(seconds=int(start))
							rt = self._real_time-start if stop==self._time_bins[-1] else stop-start
							if stop==self._time_bins[-1]:
								rt = max((rt, 1E-3*(millis[ix][-1]-millis[ix][0])))
							self.spectra[bd][ch][tm].meta = {'start_time':st,'real_time':rt,'live_time':rt-dead_time}
							if len(self.modules)>bd:
								if self.modules[bd] in self._cb:
									if len(self._cb[self.modules[bd]]):
										self.spectra[bd][ch][tm].meta = {'engcal':self._cb[self.modules[bd]][ch]}

							self.spectra[bd][ch][tm].hist = hist
						else:
							self.spectra[bd][ch][tm].hist += hist
							self.spectra[bd][ch][tm].meta['live_time'] -= dead_time


	def parse(self, fmap=None):
		"""Description

		...

		Parameters
		----------
		x : type
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""


		#### 1.8 GiB/min on 10/2/19

		self.spectra = None
		self._cb = None
		self._N_boards = 1
		self._start_time = None
		self._real_time = None
		self._time_bins = None

		for fdat in self.zp.infolist():
			if fdat.filename.endswith('.mvmelst'):
				fnm = fdat.filename
				size = fdat.file_size
			elif fdat.filename.endswith('.log'):
				self._parse_log(self.zp.read(fdat.filename).decode('utf-8'))
			elif fdat.filename.endswith('.analysis'):
				self._parse_analysis(self.zp.read(fdat.filename).decode('utf-8'))

		if self._cb is None:
			raise ValueError('Analysis file not found in zipfile')

		if self._real_time is None:
			raise ValueError('Log file not found in zipfile')

		if self.meta['time_bin_length'] is not None:
			self._time_bins = np.arange(0.0, self._real_time+self.meta['time_bin_length'], self.meta['time_bin_length'])
			self._time_bins[-1] = self._real_time
		else:
			if type(self.meta['time_bins'])==int:
				self._time_bins = np.linspace(0.0, self._real_time, self.meta['time_bins']+1)
			else:
				self._time_bins = np.array(self.meta['time_bins'], dtype=np.float64)

		fl = self.zp.open(fnm)
		if fl.read(4).decode('utf-8')!='MVME' or np.frombuffer(fl.read(4), dtype='u2')[0]!=1:
			raise ValueError('Only listfile version 1 supported.')

		head_len = np.frombuffer(fl.read(4), dtype='u2')[0]
		self._parse_header(json.loads(fl.read(4*head_len).decode('utf-8')))
		size -= (12+4*head_len)
		full_size = float(size)

		tails = np.array([], dtype='u2').reshape(0,2)
		end_times = None
		while size>0:
			read_size = min((int(size/4), int(1E8)))
			size -= read_size*4
			print('Processing {0} out of {1} MiB ({2}%)'.format(int((full_size-size)/1048576.0), int(full_size/1048576.0), round(100*(full_size-size)/full_size),1))
			fl16 = np.concatenate((tails, np.frombuffer(fl.read(read_size*4), dtype='u2').reshape(read_size, 2)))

			section_idx = np.where(fl16[:,1]==8192)[0]
			if size>0:
				tails = fl16[section_idx[-1]:]
				fl16 = fl16[:section_idx[-1]]

			dat_idx = np.where((fl16[:,1]>=4096)&(fl16[:,1]<8192))[0]
			
			flags = fl16[:,1][dat_idx]
			channel = np.bitwise_and(flags,63)
			overflow = np.bitwise_and(flags,64)
			pileup = np.bitwise_and(flags,128)
			
			good_events = np.where((channel[:-1]-channel[1:]==16))[0]
			data = fl16[:,0][dat_idx]
			adc = data[good_events+1]
			tdc = data[good_events]
			ch = channel[good_events+1]
			overflow = overflow[good_events+1]
			pileup = pileup[good_events+1]

			mod_num = np.zeros(len(fl16), dtype='u1')
			n, non_ix = 0, dat_idx[good_events]

			while(len(non_ix)):
				n += 1
				x = fl16[:,1][non_ix-n]
				mod_num[non_ix[np.where((x==1024)|(x==2048))]] += 1
				non_ix = non_ix[np.where(x!=8192)]

			mod_num = mod_num[dat_idx][good_events]-1
			if len(mod_num):
				self._N_boards = max((np.max(mod_num)+1, self._N_boards))


			long_time = np.zeros(len(fl16), dtype='u4')
			n, non_ix = 0, dat_idx[good_events]

			while(len(non_ix)):
				n += 1
				x = fl16[:,1][non_ix+n]
				eoe_ix = np.where(x>=49152, True, False)
				long_time[non_ix[eoe_ix]] = np.left_shift(np.array(x[eoe_ix]-49152, dtype='u4'), 16)+np.array(fl16[:,0][non_ix+n][eoe_ix], dtype='u4')
				non_ix = non_ix[np.logical_not(eoe_ix)]

			long_time = np.array(long_time[dat_idx][good_events], dtype='i8')
			millis = np.zeros(len(long_time))
			if end_times is None:
				end_times = np.zeros((self._N_boards, 16), dtype='i8')
			for bd in range(self._N_boards):
				tdc_res = (self.meta['tdc_resolution'][bd] if type(self.meta['tdc_resolution'])==list else self.meta['tdc_resolution'])*1E3
				for c in range(16):
					ix = np.where((ch==c)&(mod_num==bd))[0]
					t = long_time[ix]
					t[1:] += np.cumsum(np.where((t[:-1]-t[1:])>1E9, 1073741823, 0), dtype='i8')
				
					dt = np.where((t[1:]-t[:-1])>1E7, -16776704, 0)
					dt = np.where((t[:-1]-t[1:])>1E7, 16776704, dt)
					t[1:] += np.cumsum(dt, dtype='i8')
					if len(t):
						if (end_times[bd][c]-t[0])>1E9:
							t += 1073741823*int(round((end_times[bd][c]-t[0])/1073741823))
						end_times[bd][c] = t[-1]

					millis[ix] = t/16E3 + tdc[ix]*tdc_res

			if fmap is not None:
				fmap(adc, millis, mod_num, ch, overflow, pileup)
			else:
				self._default_fmap(adc, millis, mod_num, ch, overflow, pileup)

		self._parsed = True

	def save(self, resolution=2**13):
		"""Description

		...

		Parameters
		----------
		x : type
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		self.save_to_dir(self.zipfilename.replace('.zip',''))

	def save_to_dir(self, directory='', resolution=2**13):
		"""Description

		...

		Parameters
		----------
		x : type
			Description of parameter `x`.

		Returns
		-------

		Notes
		-----

		References
		----------

		Examples
		--------

		"""

		if not self._parsed:
			self.parse()

		directory = os.path.abspath(directory)
		if not os.path.exists(directory):
			os.mkdir(directory)
		path, fnm = os.path.split(self.zipfilename)
		for n,bd in enumerate(self.spectra):
			for ch,chan in enumerate(bd):
				for tm,sp in enumerate(chan):
					if sp is None:
						continue
					if resolution!=2**16:
						sp.rebin(resolution)
					if len(self._time_bins)==2:
						sp.saveas(os.path.join(directory, '{0}_b{1}_ch{2}.Spe'.format(fnm.replace('.zip',''), n, ch)))
					else:
						sp.saveas(os.path.join(directory, '{0}_b{1}_ch{2}_t{3}.Spe'.format(fnm.replace('.zip',''), n, ch, tm)))
						
		print('Listfile saved to {}'.format(directory))

