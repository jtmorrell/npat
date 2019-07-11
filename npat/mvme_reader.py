from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os, re, zipfile
import numpy as np
import datetime as dtm
import matplotlib.pyplot as plt
from npat import Spectrum

class Clover_MVME(object):
	def __init__(self, filename, spec_length=None):

		self.filename = filename
		self.spec_length = spec_length
		self.sums = None
		zp = zipfile.ZipFile(filename, 'r')
		for fdat in zp.infolist():
			if fdat.filename.endswith('.mvmelst'):
				fnm = fdat.filename
				size = fdat.file_size
			elif fdat.filename.endswith('.log'):
				self.parse_logfile(zp.read(fdat.filename).decode('utf-8'))
			elif fdat.filename.endswith('.analysis'):
				self.parse_analysisfile(zp.read(fdat.filename).decode('utf-8'))
		self.parse_listfile(zp.open(fnm), size)
	
	def parse_logfile(self, logfile):
		lines = logfile.splitlines()
		self.start_time = dtm.datetime.strptime([ln for ln in lines if re.search('readout starting', ln)][0].split(' ')[-1], '%Y-%m-%dT%H:%M:%S')
		self.start_time_str = self.start_time.strftime('%m/%d/%Y %H:%M:%S')

	def parse_analysisfile(self, analysisfile):
		for classes in analysisfile.split('"class":'):
			f1, f2 = False, False
			for line in classes.split('\n'):
				if re.search(r'.amplitude"', line):
					f1 = True
				if re.search('CalibrationMinMax', line):
					f2 = True
			if f1 and f2:
				mx, mn = [],[]
				for ln in classes.split('\n'):
					if re.search('unitMax', ln):
						mx.append(ln.split(':')[1].strip().replace(',',''))
					if re.search('unitMin', ln):
						mn.append(ln.split(':')[1].strip().replace(',',''))
				mx, mn = np.array(list(map(float, mx))), np.array(list(map(float, mn)))
				self.calibration = np.array([mn, (mx-mn)/float(2**16)]).T

	def parse_listfile(self, listfile, size):
		self.channels = [[[np.zeros(2**16),0.0,0.0]] for n in range(16)]
		self.prev = [[[0.0,0.0]] for n in range(16)]

		self.filtered_spectrum = np.zeros(2**16)
		leftovers = np.array([], dtype=np.uint32)
		last_overflow = 0
		end_real_time = 0
		s0, si = float(size), 0

		while size>0:

			read_size = min((int(size/4), int(1E8)))
			si += read_size
			print(round(100.0*4.0*si/s0,1), '% Complete')
			size -= read_size*4

			### End of Event ###
			fl32 = np.concatenate((leftovers,np.frombuffer(listfile.read(read_size*4), dtype=np.uint32)))
			e_idx = np.where(fl32>=3221225472)
			


			if len(e_idx[0])==0:
				break

			real_time = np.array(fl32[e_idx], dtype='u8')-3221225472
			leftovers = fl32[e_idx[0][-1]+1:]
			fl32 = fl32[:e_idx[0][-1]+1]
			
			x8 = fl32.view(dtype=np.dtype((np.uint32, {'0':(np.uint8,0), '1':(np.uint8,1), '2':(np.uint8,2), '3':(np.uint8,3)})))
			fl8 = np.array([x8['0'],x8['1'],x8['2'],x8['3']], dtype=np.uint8).T
			x16 = fl32.view(dtype=np.dtype((np.uint32, {'0':(np.uint16,0), '1':(np.uint16,2)})))
			fl16 = np.array([x16['0'],x16['1']], dtype=np.uint16).T

			# x8 = fl32.astype(np.dtype((np.uint32, {'0':(np.uint8,0), '1':(np.uint8,1), '2':(np.uint8,2), '3':(np.uint8,3)})))
			# fl8 = np.array([x8['0'],x8['1'],x8['2'],x8['3']], dtype=np.uint8).T
			# x16 = fl32.astype(np.dtype((np.uint32, {'0':(np.uint16,0), '1':(np.uint16,2)})))
			# fl16 = np.array([x16['0'],x16['1']], dtype=np.uint16).T



			### Header ###
			h_idx = np.where((fl8[:,3]>=64)&(fl8[:,3]<128))

			# print(h_idx[0][:20])

			# print(np.right_shift(np.left_shift(fl16[h_idx[0],0][:20],6),6))

			# for n,ln in enumerate(fl8[h_idx][:1000]):
			# 	print([''.join(map(str, np.unpackbits(l))) for l in ln[::-1]])

			module_id = fl8[h_idx][:,2]
			# print(len(h_idx[0]), len(e_idx[0]))
			# print(np.unpackbits(np.array([module_id[:40]]).T, axis=1))
			# print(set(module_id))
			# plt.hist(module_id)
			# plt.yscale('log')
			# plt.show()


			### Data ###
			d_idx = np.where((fl8[:,3]<32)&(fl8[:,3]>=16))
			bits, fl8 = fl8[d_idx][:,2], None
			pileup = np.array(np.right_shift(bits,7),dtype='b')
			overflow = np.array(np.right_shift(np.left_shift(bits,1),7),dtype='b')
			trigger = np.array(np.right_shift(np.left_shift(bits,2),7),dtype='b')
			channel = np.right_shift(np.left_shift(bits,3),3)
			print(set(channel))
			tg_idx = np.where((trigger==0)&(overflow==0))
			d_idx, pileup, channel = d_idx[0][tg_idx], pileup[tg_idx], channel[tg_idx]
			while channel[-1]>15:
				d_idx, pileup, channel = d_idx[:-1], pileup[:-1], channel[:-1]

			fl16 = fl16[d_idx][:,0]
			chn_pairs = np.where((channel[channel>15]-channel[np.where(channel>15)[0]+1])==16)
			tdc_idx = np.where(channel>15)[0][chn_pairs]

			TDC = fl16[tdc_idx]
			ADC = fl16[tdc_idx+1]
			pileup = pileup[tdc_idx+1]
			ADC_ch = channel[tdc_idx]-16

			ov_idx = np.append((np.where(real_time[:-1]>real_time[1:]+int(1E6))[0]+1),[len(real_time)-1])
			real_time += last_overflow*1073741823
			if real_time[0]+int(1E6)<end_real_time:
				last_overflow += 1
				real_time += 1073741823
			for n,idx in enumerate(ov_idx[1:]):
				real_time[ov_idx[n]:idx+1] += (n+1)*1073741823
			last_overflow += len(ov_idx)-1
			end_real_time = real_time[-1]

			e_idx = np.insert(e_idx[0],0,0)
			real_time = np.repeat(real_time,e_idx[1:]-e_idx[:-1])[d_idx][tdc_idx]
			TDC = (TDC*96E-12)+(real_time/16E6)

			pileup = [pileup[ADC_ch==i] for i in range(16)]

			for n in range(16):
				dat = [ADC[ADC_ch==n],TDC[ADC_ch==n]]
				dead = np.append(dat[1],(dat[1][-1] if len(dat[1])>0 else 0.0))
				A = dat[0][(pileup[n]==0)&(dat[0]>100)]
				LT = dat[1][(pileup[n]==0)&(dat[0]>100)]-np.cumsum(np.where(pileup[n],dead[1:]-dead[:-1],300E-9))[(pileup[n]==0)&(dat[0]>100)]
				RT = dat[1][(pileup[n]==0)&(dat[0]>100)]

				self.update_channels(A, LT, RT, n)

	def update_channels(self, ADC, LT, RT, channel):
		if self.spec_length is None:
			self.channels[channel][0][0] += np.histogram(ADC, bins=np.arange(-0.5,2**16+0.5,1))[0]
			if len(LT)>0:
				self.channels[channel][0][1] = LT[-1]
				self.channels[channel][0][2] = RT[-1]
		else:
			while len(RT)>0:
				if self.channels[channel][-1][2]<self.spec_length:
					w = np.where(RT<=(self.prev[channel][-1][1]+self.spec_length),True,False)
					if len(LT[w])>0:
						self.channels[channel][-1][0] += np.histogram(ADC[w], bins=np.arange(-0.5,2**16+0.5,1))[0]
						self.channels[channel][-1][1] = LT[w][-1]-self.prev[channel][-1][0]
						self.channels[channel][-1][2] = RT[w][-1]-self.prev[channel][-1][1]
						w = np.where(w,False,True)
						ADC = ADC[w]
						LT = LT[w]
						RT = RT[w]
					else:
						self.channels[channel][-1][2] = self.spec_length-(self.channels[channel][-1][1]-self.channels[channel][-1][2])
						self.channels[channel][-1][1] = self.spec_length
						self.prev[channel].append([LT[0], RT[0]])
						self.channels[channel].append([np.zeros(2**16),0.0,0.0])
				else:
					self.prev[channel].append([LT[0], RT[0]])
					self.channels[channel].append([np.zeros(2**16),0.0,0.0])
					w = np.where(RT<=(self.prev[channel][-1][1]+self.spec_length),True,False)
					self.channels[channel][-1][0] += np.histogram(ADC[w], bins=np.arange(-0.5,2**16+0.5,1))[0]
					self.channels[channel][-1][1] = LT[w][-1]-self.prev[channel][-1][0]
					self.channels[channel][-1][2] = RT[w][-1]-self.prev[channel][-1][1]
					w = np.where(w,False,True)
					ADC = ADC[w]
					LT = LT[w]
					RT = RT[w]

	def rebin(self, N_bits=14):
		self.channels = [[[np.sum(c[0].reshape((2**N_bits,int(2**16/2**N_bits))),axis=1),c[1],c[2]] for c in ch] for ch in self.channels]
		self.calibration = np.array([[cb[0], (2**16/2**N_bits)*cb[1]] for cb in self.calibration])


	def summation(self, groups):
		self.sums = []
		for group in groups:
			self.sums.append([])
			L = len(self.channels[group[0]][0][0])
			
			x = np.arange(L, dtype=np.int32)
			cb = [self.calibration[g] for g in group]
			N_cb = sorted([[n,cb[n][0]+L*cb[n][1]] for n in range(len(cb))], key=lambda h:h[1])[-1][0]
			M_lu = [[cb[u][1]/cb[l][1] for u in range(len(cb))] for l in range(len(cb))]
			B_lu = [[(cb[u][0]-cb[l][0])/cb[l][1] for u in range(len(cb))] for l in range(len(cb))]

			for t in range(len(self.channels[group[N_cb]])):
				sm = np.zeros(L)
				for m,g in enumerate(group):
					if t>=len(self.channels[g]):
						continue
					spec = self.channels[g][t][0]
					for o in np.arange(0.0,1.0,0.005):
						a = np.array(np.round(B_lu[N_cb][m]+(x+o)*M_lu[N_cb][m], 0), dtype=np.int32)
						a = np.where(a<L, a, L-1)
						sm[a] += 0.005*spec
				self.sums[-1].append([np.array(sm, dtype=np.int32), self.channels[group[N_cb]][t][1], self.channels[group[N_cb]][t][2], cb[N_cb], group])


	def save(self):
		if not os.path.exists(self.filename.replace('.zip','')):
			os.mkdir(self.filename.replace('.zip',''))
		path, fnm = os.path.split(self.filename)
		for n,ch in enumerate(self.channels):
			for m,c in enumerate(ch):
				if np.all(c[0]==0):
					continue
				sp = Spectrum()
				sp.hist = c[0]
				sp.meta = {'start_time':self.start_time+dtm.timedelta(seconds=self.prev[n][m][1]),'live_time':c[1], 'real_time':c[2], 'engcal':self.calibration[n]}
				if self.spec_length is None:
					sp.saveas('{0}/{1}/{1}_ch{2}.Spe'.format(path, fnm.replace('.zip',''), n))
				else:
					sp.saveas('{0}/{1}/{1}_ch{2}_t{3}.Spe'.format(path, fnm.replace('.zip',''), n, m))
		if self.sums is not None:
			for grp in self.sums:
				for m,sm in enumerate(grp):
					if np.all(sm[0]==0):
						continue
					sp = Spectrum()
					sp.hist = sm[0]
					sp.meta = {'start_time':self.start_time,'live_time':sm[1], 'real_time':sm[2], 'engcal':sm[3]}
					if self.spec_length is None:
						sp.saveas('{0}/{1}/{1}_sum_ch{2}.Spe'.format(path, fnm.replace('.zip',''), '-'.join(list(map(str, sm[4])))))
					else:
						sp.saveas('{0}/{1}/{1}_sum_ch{2}_t{3}.Spe'.format(path, fnm.replace('.zip',''), '-'.join(list(map(str, sm[4]))), m))


if __name__=="__main__":
	### arguments: filename, (optional) time bins in seconds - can also be None
	lst = Clover_MVME('/home/jmorrell/Documents/mvme_testing/listfiles/Fluffy_001_190420_131634.zip')

	### arguments: (optional) number of bits for rebinning - default 14
	lst.rebin()

	### arguments: groups of channels to sum together
	lst.summation([[0,1,2,3],[4,5,6,7]])

	### saves as ASCII .Spe in new directory with same name as .zip file
	lst.save()

