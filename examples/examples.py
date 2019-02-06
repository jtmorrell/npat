from npat import *
import matplotlib.pyplot as plt
import datetime as dtm
import numpy as np

if __name__=='__main__':

	######### Spectroscopy Tests #########

	sp = Spectrum('eu_calib_7cm.Spe','test.db')
	sp.meta = {'istp':['152EU'], 'A0':3.7E4, 'ref_date':'01/01/2009 12:00:00'}

	# sp.fit_config = {'xrays':True, 'E_min':20.0}
	sp.auto_calibrate()
	sp.cb.plot()
	sp.summarize()
	sp.plot()

	sp = Spectrum('AH20190131_Na22_10cm_IDM4.Chn','test.db')
	sp = Spectrum('AH20190131_Na22_10cm_IDM4.Spe')

	sp.saveas('AH20190131_Na22_10cm_IDM4 (copy).Chn')
	sp.saveas('AH20190131_Na22_10cm_IDM4 (copy).spe')

	sp = Spectrum('AH20190131_Na22_10cm_IDM4 (copy).Chn')
	sp.plot()



	########## Decay Chain Tests ##########

	dc = DecayChain('225RA','d',R={'225RA':9.0,'225AC':1.0},time=10.0/24.0)
	dc.append(DecayChain('225RA','d',R={'225RA':2.0,'225AC':1.0},time=33.0/24.0))
	dc.append(DecayChain('225RA','d',R={'225RA':5.0,'225AC':1.0},time=113.0/24.0))
	dc.append(DecayChain('225RA','d',time=21))
	dc.counts = {'225AC':[[5.0, 6.0, 6.0, 0.2],[6.0, 7.0, 6.2, 0.3]],'221FR':[5.5, 6.5, 6.0, 0.2]}
	dc.fit_R()
	dc.plot(N_plot=10, logscale=False)

	######## Ziegler Tests #############

	zg = Ziegler(stack=[{'compound':'Ni','name':'Ni01','thickness':0.025},{'compound':'Ti','name':'Ti01','thickness':0.025},{'compound':'U','name':'U01','thickness':0.7,'density':20.0}], beam={'istp':'2H', 'N':1E5})
	zg.plot(['U','Ni','Ti'])
	zg.summarize()


	####### Isotope Tests ############

	ISTP = Isotope('133BA')
	print(ISTP.TeX)
	print(ISTP.mass)
	print(ISTP.half_life(ISTP.optimum_units(),unc=True),ISTP.optimum_units())
	print(ISTP.gammas())


	####### Reaction Tests ##########

	f, ax = None, None
	for lb in ['irdff','endf','tendl']:
		rx = Reaction('90ZR(n,2n)89ZR', lb)
		f, ax = rx.plot(f=f, ax=ax, show=False, label='library', title=True)

	plt.show()
	

	f, ax = None, None
	for lb in ['endf','tendl']:
		rx = Reaction('226RA(n,2n)225RA', lb)
		f, ax = rx.plot(f=f, ax=ax, show=False, label='both')
		rx = Reaction('226RA(n,3n)224RA', lb)
		f, ax = rx.plot(f=f, ax=ax, show=False, label='both', title=True)

	plt.show()

	lb = Library('tendl_n')
	print(lb.search(target='226RA',product='225RAg'))