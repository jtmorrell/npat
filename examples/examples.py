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
	sp.cb.saveas('eu_calib.json')
	sp.cb.open('eu_calib.json')
	sp.cb.plot()
	cb = Calibration('eu_calib.json')
	cb.plot()

	### Save peak information
	sp.saveas('test.db', 'test.csv')

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

def listfile_examples():
	from npat import MVME
	
	fl = MVME('mvmelst_007.zip')
	fl.save()


def ziegler_examples():
	from npat import Ziegler
	
	### Set up a stack with different input options.
	zg = Ziegler(stack=[{'compound':'Ni', 'name':'Ni01', 'thickness':0.025},  # Thickness only (mm)
						{'compound':'Kapton', 'thickness':0.05},				# No name - will not be tallied
						{'compound':'Ti', 'name':'Ti01', 'thickness':1.025},  # Very thick: should see straggle
						{'compound':{'inconel':[[26,33.0],[28,55.0]]},'ad':1.0,'name':'test'},
						{'compound':'SrCO3', 'name':'SrCO3', 'area':0.785, 'mass':4.8E-3}],  # Mass (g) and area (cm^2)
						beam_istp='2H', N=1E5, max_steps=100, E0=33.0)  ## 33 MeV deuteron beam

	### zg.stack contains all information, both input and calculated
	print(zg.stack)

	### Print mean energies on samples
	zg.summarize()

	### Plot only strontium and titanium fluxes
	zg.plot(['Sr', 'Ti'])

	### Find out if 6mm of Be will stop a deuteron beem
	zg = Ziegler(stack=[{'compound':'Be', 'name':'Be Breakup','thickness':6.0}])
	### Set beam options with zg.meta
	zg.meta = {'istp':'2H', 'E0':33.0}

	### Summarize, plot and save
	zg.summarize()
	zg.plot()
	zg.saveas('breakup.csv', 'breakup.db', 'breakup.png')

	### Import stack design from .csv file
	zg = Ziegler(stack='test_stack.csv')
	zg.meta = {'istp':'4HE','E0':70.0}
	zg.plot()
	

def decay_chain_examples():
	from npat import DecayChain

	### 225RA decay chain, units of days, 9.0/day production rate, for 0.5 days
	dc = DecayChain('225RA', 'd', R=9.0, time=0.5)
	dc.plot()

	### Additional production of 225AC, with production rate of 225RA fluctuating
	dc.append(DecayChain('225RA', 'd', R={'225RA':2.0, '225AC':1.0}, time=1.5))
	dc.append(DecayChain('225RA', 'd', R={'225RA':5.0, '225AC':1.0}, time=4.5))

	### 21 day decay time
	dc.append(DecayChain('225RA', 'd', time=21))

	### Measured counts: [start_time (d), stop_time (d), decays, unc_decays]
	### Times relative to last appended DecayChain, i.e. EoB time
	dc.counts = {'225AC':[[5.0, 5.1, 6E5, 2E4],
						  [6.0, 6.1, 7E5, 3E4]],
				'221FR':[5.5, 5.6, 6E5, 2E4]}

	### Find the scaled production rate that gives us these counts
	dc.fit_R()
	### Only plot the 5 most active isotopes in the decay chain
	dc.plot(N_plot=5)

def isotope_examples():
	from npat import Isotope

	i = Isotope('60CO')
	### Get LaTeX formatted name
	print(i.TeX)
	### Get isotope mass in amu
	print(i.mass)
	### Get half life in optimum units
	print(i.half_life(i.optimum_units(),unc=True), i.optimum_units())
	### Print list of the decay gammas
	print(i.gammas()['E'])

def reaction_examples():
	from npat import Reaction, Library
	import matplotlib.pyplot as plt

	### We will plot the same reaction from three different libraries
	### Passing f,ax to rx.plot allows multiple plots on the same figure
	f, ax = None, None
	for lb in ['irdff','endf','tendl']:
		rx = Reaction('90ZR(n,2n)89ZR', lb)
		f, ax = rx.plot(f=f, ax=ax, show=False, label='library', title=True)

	plt.show()
	
	### Compare (n,2n) and (n,3n) for endf vs tendl
	f, ax = None, None
	for lb in ['endf','tendl']:
		rx = Reaction('226RA(n,2n)225RA', lb)
		f, ax = rx.plot(f=f, ax=ax, show=False, label='both', E_lim=[0,30], logscale=True)
		rx = Reaction('226RA(n,3n)224RA', lb)
		f, ax = rx.plot(f=f, ax=ax, show=False, label='both', title=True, E_lim=[0,40])

	plt.show()

	### Search the TENDL-2015 neutron library for reactions producing 225RA from 226RA
	lb = Library('tendl_n')
	print(lb.search(target='226RA',product='225RAg'))


if __name__=='__main__':

	# spectroscopy_examples()
	listfile_examples()
	# decay_chain_examples()
	# ziegler_examples()
	# isotope_examples()
	# reaction_examples()