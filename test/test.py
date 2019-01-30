import numpy as np
from npat import *
import matplotlib.pyplot as plt

if __name__=='__main__':

	######### Spectroscopy Tests #########

	# sp = Spectrum('eu_calib_7cm.Spe')
	# sp.meta = {'istp':['152EU'], 'A0':3.7E4, 'ref_date':None}
	# sp.fit_config = {'xrays':True, 'E_min':20.0}
	# # cb = Calibration()
	# # cb.calibrate([sp], auto_calibrate=True)
	# # cb.plot()
	# sp.auto_calibrate()
	# sp.cb.plot()
	# sp.summarize()
	# sp.plot()

	# 
	# sp = Spectrum('La01.Spe')
	# sp.meta = {'istp': ['135CE','134CE','134LA','137CE','137CEm','139CE','133BA','133BAm','132CS','135LA','22NA','24NA']}
	sp = Spectrum('/home/jmorrell/Documents/Radium_Bernstein_Oct2018/data/count_room_data/experiment/225Ac_separated_20cm_000.Spe')
	sp.meta = {'istp':['221FR','221RN','213BI','209TL','226RA','214PB','214BI','228AC','210TL','210BI','212PB','212BI','212PO','208TL','40K']}
	# sp.fit_config = {'skew_fit':True}
	# sp.auto_calibrate()
	# sp = Spectrum('iridium_7cm.Spe')
	# sp.meta = {'istp':['192IR','194IR','193IRm']}
	sp.summarize()
	sp.plot()



	########## Decay Chain Tests ##########

	# dc = DecayChain('225RA','d',R={'225RA':9.0,'225AC':1.0},time=10.0/24.0)
	# dc.append(DecayChain('225RA','d',R={'225RA':2.0,'225AC':1.0},time=33.0/24.0))
	# dc.append(DecayChain('225RA','d',R={'225RA':5.0,'225AC':1.0},time=113.0/24.0))
	# dc.append(DecayChain('225RA','d',time=21))
	# dc.counts = {'225AC':[[5.0, 6.0, 6.0, 0.2],[6.0, 7.0, 6.2, 0.3]],'221FR':[5.5, 6.5, 6.0, 0.2]}
	# dc.fit_R()
	# dc.plot(N_plot=10, logscale=False)

	######## Ziegler Tests #############

	# zg = Ziegler(stack=[{'compound':'Ni','name':'Ni01','thickness':0.025},{'compound':'Ti','name':'Ti01','thickness':0.025},{'compound':'U','name':'U01','thickness':1.0,'density':20.0}], beam={'istp':'2H', 'N':1E5})
	# zg.plot(['U','Ni','Ti'])
	# zg.summarize(saveas='test.csv')


	####### Isotope Tests ############

	# ISTP = Isotope('133CS')
	# print ISTP.name
	# print ISTP.element
	# print ISTP.A
	# print ISTP.isomer
	# print ISTP.isotope
	# print ISTP.E_level
	# print ISTP.decay_mode
	# print ISTP.TeX
	# print ISTP.amu
	# print ISTP.half_life(ISTP.optimum_units(),unc=True),ISTP.optimum_units()
	# print ISTP.gammas()
	# print ISTP.electrons()
	# print ISTP.beta_minus()
	# print ISTP.beta_plus()
	# print ISTP.alphas()


	####### Reaction Tests ##########3

	# lb = Library('irdff')
	# print lb.search()

	# f, ax = None, None
	# for lb in ['irdff','endf','tendl']:
	# 	rx = Reaction('89Y(n,2n)88Y', lb)
	# 	f, ax = rx.plot(f=f, ax=ax, show=False, label='library')

	# plt.show()
	
	# f, ax = None, None
	# for lb in ['irdff','endf','tendl']:
	# 	rx = Reaction('90ZR(n,2n)89ZR', lb)
	# 	f, ax = rx.plot(f=f, ax=ax, show=False, label='library')

	# plt.show()
	

	# f, ax = None, None
	# for lb in ['endf','tendl']:
	# 	rx = Reaction('226RA(n,2n)225RA', lb)
	# 	f, ax = rx.plot(f=f, ax=ax, show=False, label='both')
	# 	rx = Reaction('226RA(n,3n)224RA', lb)
	# 	f, ax = rx.plot(f=f, ax=ax, show=False, label='both')

	# plt.show()

	# print rx.target, rx.product
	# for lbr in ['endf','tendl','tendl_n_rp','irdff']:
	# 	lb = Library(lbr)
	# 	print len(lb.query(target='58FE',outgoing='g',product='59FE'))
	# for lbr in ['iaea','tendl_p_rp']:
	# 	lb = Library(lbr)
	# 	print len(lb.query(target='63CU',incident='p',product='62ZN'))
	# for lbr in ['iaea','tendl_d_rp']:
	# 	lb = Library(lbr)
	# 	print len(lb.query(target='63CU',incident='d',product='62ZN'))
	# lb = Library('tendl')
	# print lb.query(target='226RA', outgoing='2n')
	# print lb.search(product='63ZN',incident='p')
	# print lb.check(product='115INm',target='115IN')