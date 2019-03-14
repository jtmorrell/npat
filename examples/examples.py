### npat examples

def spectroscopy_examples():
	from npat import Spectrum, Calibration
	### Load and plot the spectrum
	sp = Spectrum('eu_calib_7cm.Spe')
	sp.plot()

	### Fit Europium Spectrum
	sp.meta = {'istp':['152EU']}
	sp.plot()

	### Perform efficiency calibration
	sp.meta = {'A0':3.7E4, 'ref_date':'01/01/2009 12:00:00'}
	sp.auto_calibrate()
	sp.cb.plot()

	### Save peak information
	sp.saveas('test.db')
	sp.saveas('test.csv')

	### Print out peaks
	sp.summarize()

	### Save as .Chn format
	sp.saveas('eu_calib_7cm.Chn')

	### Load with database
	sp = Spectrum('eu_calib_7cm.Chn', 'test.db')
	sp.meta = {'istp':['152EU'], 'A0':3.7E4, 'ref_date':'01/01/2009 12:00:00'}

	### Plot ADC channels instead of energy
	sp.plot(xcalib=False)

	### Pick out a few peaks for manual calibration
	cb_data = [[664.5, 121.8],
				[1338.5, 244.7],
				[1882.5, 344.3],
				[2428, 444],
				[7698, 1408]]

	sp.auto_calibrate(data=cb_data)

	### Efficiency calibration using the Calibration class
	cb = Calibration()
	cb.calibrate([sp])
	cb.plot()

	### Custom peaks
	sp.fit_config = {'p0':[{'E':1460.82, 'I':0.1066, 'dI':0.0017, 'istp':'40K'}]}
	sp.summarize()
	sp.plot()

	### More detailed fits
	sp.fit_config = {'xrays':True, 'E_min':20.0, 'bg_fit':True, 'quad_bg':True}
	### Save and show the plot
	sp.plot(saveas='europium.png')


def ziegler_examples():
	from npat import Ziegler
	
	zg = Ziegler(stack=[{'compound':'Ni','name':'Ni01','thickness':0.025},{'compound':'Ti','name':'Ti01','thickness':1.025},{'compound':'SrCO3','name':'SrCO3','thickness':0.7,'density':3.5}], beam_istp='2H', N=1E4, max_steps=100)
	print(zg.stack)
	zg.summarize()
	zg.plot(['Sr','Ni','Ti'])

	zg = Ziegler(stack=[{'compound':'Be','name':'Be Breakup','thickness':6.0}])
	zg.meta = {'istp':'2H','E0':33.0}
	zg.summarize()
	zg.plot()
	zg.saveas('breakup.csv', 'breakup.db', 'breakup.png')

	zg = Ziegler(stack='test_stack.csv')
	zg.meta = {'istp':'2H'}
	zg.plot()
	

def decay_chain_examples():
	from npat import DecayChain

	dc = DecayChain('225RA','d',R={'225RA':9.0,'225AC':1.0},time=10.0/24.0)
	dc.append(DecayChain('225RA','d',R={'225RA':2.0,'225AC':1.0},time=33.0/24.0))
	dc.append(DecayChain('225RA','d',R={'225RA':5.0,'225AC':1.0},time=113.0/24.0))
	dc.append(DecayChain('225RA','d',time=21))
	dc.counts = {'225AC':[[5.0, 6.0, 6.0, 0.2],[6.0, 7.0, 6.2, 0.3]],'221FR':[5.5, 6.5, 6.0, 0.2]}
	dc.fit_R()
	dc.plot(N_plot=10, logscale=False)

def isotope_examples():
	from npat import Isotope

	i = Isotope('133BA')
	print(i.TeX)
	print(i.mass)
	print(i.half_life(i.optimum_units(),unc=True), i.optimum_units())
	print(i.gammas())

def reaction_examples():
	from npat import Reaction, Library
	import matplotlib.pyplot as plt

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


if __name__=='__main__':

	spectroscopy_examples()
	decay_chain_examples()
	ziegler_examples()
	isotope_examples()
	reaction_examples()