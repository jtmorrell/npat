import numpy as np
from npat import *
import matplotlib.pyplot as plt

if __name__=='__main__':

	########## Decay Chain Tests ##########

	ch = DecayChain('135CE', units='h')

	# ch = DecayChain(['225RA'],units='d')
	# t_i = np.arange(0,6,0.1)
	# t_d = np.arange(0,50,0.1)
	# A_r = ch.production_activity('225RA',t_i)
	# A_a = ch.production_activity('225AC',t_i)
	# A_r0 = A_r[-1]
	# A_a0 = A_a[-1]
	# A_rd = ch.decay_activity('225RA',t_d,A_r0)
	# A_ad = ch.decay_activity('225AC',t_d,A_r0)+A_a0*np.exp(-np.log(2.0)*t_d/10.0)
	# plt.plot((t_i-t_i[-1]).tolist()+t_d.tolist(),A_r.tolist()+A_rd.tolist(),label=r'$^{225}$Ra Activity')
	# plt.plot((t_i-t_i[-1]).tolist()+t_d.tolist(),A_a.tolist()+A_ad.tolist(),label=r'$^{225}$Ac Activity')
	# # Peak at 15.2 days for 5 day irradiation, 0.1003
	# # Peak at 14.7 days for 6 day irradiation, 0.1216
	# # Peak at 14.3 days for 7 day irradiation, 0.1413
	# # Peak at 13.8 days for 8 day irradiation, 0.1615
	# plt.xlabel('Time Since EOB (d)')
	# plt.ylabel('Activity (rel.)')
	# plt.legend(loc=0)
	# plt.show()

	# ch.plot('production_activity', logscale=False, **{'t':np.arange(0,700,0.1)})
	# ch.plot('decay_activity',logscale=False,**{'t':np.arange(0,50,0.01)})

	######## Ziegler Tests #############

	# zg = Ziegler(stack=[{'compound':'Ni','name':'Ni01','thickness':0.025},{'compound':'Ti','name':'Ti01','thickness':0.025},{'compound':'U','name':'U01','thickness':1.0}], beam={'istp':'2H', 'N':1E5})
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