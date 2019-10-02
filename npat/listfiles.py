from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, re, zipfile, json
import numpy as np

from npat import Spectrum


class MVME(object):
	def __init__(self, filename):
		self.zipfilename = filename
		self.zp = zipfile.ZipFile(filename, 'r')
		_path, _fnm = os.path.split(filename)
		print('Reading {}'.format(_fnm))
		self.N_boards = 1
		self.spectra = None
		self.cal = None
		self.start_time = None
		self.real_time = None
		self.tdc_resolution = 24E-12

	def _parse_log(self, fl):
		import datetime as dtm

		self.real_time = 0.0
		for ln in fl.splitlines():
			if re.search('readout starting', ln):
				self.start_time = dtm.datetime.strptime(ln.split(' ')[-1], '%Y-%m-%dT%H:%M:%S')
			elif re.search('readout stopped', ln):
				self.real_time = float((dtm.datetime.strptime(ln.split(' ')[-1], '%Y-%m-%dT%H:%M:%S')-self.start_time).total_seconds())

	def _parse_analysis(self, fl):
		fl = json.loads(fl)
		modules = [str(c['moduleName']) for c in fl['AnalysisNG']['properties']['ModuleProperties']]
		self.cal = {i:[] for i in modules}
		for n,m in enumerate(modules):
			for c in fl['AnalysisNG']['operators']:
				if c['class'].endswith('CalibrationMinMax'):
					if c['name']==m+'.amplitude' or c['name']==m+'cal.amplitude':
						cl = [[0.0, 1.0] for i in range(16)]
						for j,cb in enumerate(c['data']['calibrations']):
							mn, mx = float(cb['unitMin']), float(cb['unitMax'])
							cl[j] = [mn, (mx-mn)/float(2**16)]
						self.cal[modules[n]] = cl

	def _parse_header(self, head):
		self.modules = [str(m['name']) for m in head['DAQConfig']['events'][0]['modules']]

	def _default_fmap(self, amplitude, millis, board, channel, overflow, pileup):
		if self.spectra is None:
			self.spectra = [[None for i in range(16)] for n in range(self.N_boards)]
		for bd in range(self.N_boards):
			for ch in range(16):
				ix = np.where((channel==ch)&(board==bd))[0]

				dead_time = 300E-9*len(amplitude[ix])
				pu_ix = np.where(pileup[ix])[0]
				if len(pu_ix):
					if pu_ix[-1]==len(ix)-1:
						pu_ix = pu_ix[:-1]
					diff = (millis[ix][pu_ix+1]-millis[ix][pu_ix])*1E3
					dead_time += np.sum(np.where((diff>0)&(diff<1), diff, 0.0))
			
				hist = np.histogram(amplitude[ix][np.where((pileup[ix]==0)&(overflow[ix]==0)&(amplitude[ix]>0))], bins=np.arange(-0.5,2**16+0.5,1))[0]
				if np.any(hist[500:]>0):
					if self.spectra[bd][ch] is None:
						self.spectra[bd][ch] = Spectrum()
						self.spectra[bd][ch].meta = {'start_time':self.start_time,'real_time':self.real_time,'live_time':self.real_time-dead_time}
						if len(self.modules)>bd:
							if self.modules[bd] in self.cal:
								if len(self.cal[self.modules[bd]]):
									self.spectra[bd][ch].meta = {'engcal':self.cal[self.modules[bd]][ch]}
						self.spectra[bd][ch].hist = hist
					else:
						self.spectra[bd][ch].hist += hist
						self.spectra[bd][ch].meta['live_time'] -= dead_time


	def _parse(self, fmap=None):
		#### 1.8 GiB/min on 10/2/19
		for fdat in self.zp.infolist():
			if fdat.filename.endswith('.mvmelst'):
				fnm = fdat.filename
				size = fdat.file_size
			elif fdat.filename.endswith('.log'):
				self._parse_log(self.zp.read(fdat.filename).decode('utf-8'))
			elif fdat.filename.endswith('.analysis'):
				self._parse_analysis(self.zp.read(fdat.filename).decode('utf-8'))

		if self.cal is None:
			raise ValueError('Analysis file not found in zipfile')

		if self.real_time is None:
			raise ValueError('Log file not found in zipfile')

		fl = self.zp.open(fnm)
		if fl.read(4).decode('utf-8')!='MVME' or np.frombuffer(fl.read(4), dtype='u2')[0]!=1:
			raise ValueError('Only listfile version 1 supported.')

		head_len = np.frombuffer(fl.read(4), dtype='u2')[0]
		self._parse_header(json.loads(fl.read(4*head_len).decode('utf-8')))
		size -= (12+4*head_len)
		full_size = float(size)

		tails = np.array([], dtype='u2').reshape(0,2)
		while size>0:
			### Reads approx. 10MB/s
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
				self.N_boards = max((np.max(mod_num)+1, self.N_boards))


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
			for bd in range(self.N_boards):
				for c in range(16):
					ix = np.where((ch==c)&(mod_num==bd))[0]
					t = long_time[ix]
					t[1:] += np.cumsum(np.where((t[:-1]-t[1:])>1E9, 1073741823, 0), dtype='i8')
				
					dt = np.where((t[1:]-t[:-1])>1E7, -16776704, 0)
					dt = np.where((t[:-1]-t[1:])>1E7, 16776704, dt)
					t[1:] += np.cumsum(dt, dtype='i8')

					millis[ix] = t/16E9 + tdc[ix]*(self.tdc_resolution*1E-3)


			if fmap is not None:
				fmap(adc, millis, mod_num, ch, overflow, pileup)
			else:
				self._default_fmap(adc, millis, mod_num, ch, overflow, pileup)

	def save(self, resolution=2**13):
		self._parse()
		if not os.path.exists(self.zipfilename.replace('.zip','')):
			os.mkdir(self.zipfilename.replace('.zip',''))
		path, fnm = os.path.split(self.zipfilename)
		for n,bd in enumerate(self.spectra):
			for ch,sp in enumerate(bd):
				if sp is None:
					continue
				if resolution!=2**16:
					sp.rebin(resolution)
				sp.saveas('{0}/{1}/{1}_b{2}_ch{3}.Spe'.format(path, fnm.replace('.zip',''), n, ch))




if __name__=="__main__":
	fl = MVME('/home/jmorrell/Documents/mvme_testing/listfiles/Fluffy_001_190420_131634.zip')
	fl.save()
	# for fnm in list(os.walk('/home/jmorrell/Documents/mvme_testing/listfiles'))[0][2]:
	# 	fl = MVME('/home/jmorrell/Documents/mvme_testing/listfiles/'+fnm)
	# 	fl.save()