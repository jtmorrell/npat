from npat import Ziegler, Isotope

if __name__=='__main__':

	######## Ziegler Tests #############

	## Compound must be specified, and enough info to determine areal density
	## Units:
	## thickness: mm
	## mass: g
	## area: cm^2
	## ad (areal density): mg/cm^2
	## density: g/cm^3

	stack = [{'compound':'Ni', 'name':'Ni01', 'thickness':0.025},
			{'compound':'Ti', 'name':'Ti01', 'area':10.0, 'mass':0.100},
			{'compound':'U', 'name':'U01', 'ad':1200.0}]

	## must specify stack and beam...
	## stack is list of dicts, which specify composition of each
	## element of the stack
	##
	## beam is dict, specifying the isotope and incident energy
	## of the incoming beam.  dE0 (1-sigma width) and N (number of particles)
	## can optionally be included
	zg = Ziegler(stack=stack, beam={'istp':'2H', 'E0':33.0, 'N':1E5})

	## name from stack is used for plotting - either incl_no_names=True must 
	## be set, or name must be specified
	## if zg.plot() is called, all foils in stack are plotted that have a name.
	## optionally the first argument to plot is a list of samples that will be 
	## regex matched to the names in the stack
	zg.plot(['U','Ni','Ti'])

	## summarize either prints out the mean and 1-sigma energies of each
	## foil in the stack (with a name), or saves this info to a .csv using
	## the saveas argument.
	zg.summarize(saveas='test.csv')
	
	## The results of the stack calculation can be found in the zg.stack
	## variable, which is a list of Sample objects which store the energy, 
	## flux, and other info about each sample in the stack calculation.
	print zg.stack[2].meta
	print zg.stack[1].energy
	print zg.stack[1].flux


	####### Isotope Tests ############

	# ISTP = Isotope('133BA')
	# print ISTP.name
	# print ISTP.element
	# print ISTP.A
	# print ISTP.isomer
	# print ISTP.isotope
	# print ISTP.E_level
	# print ISTP.TeX
	# print ISTP.mass
	# print ISTP.half_life(ISTP.optimum_units(),unc=True),ISTP.optimum_units()
	# print ISTP.gammas()
	# print ISTP.electrons()
	# print ISTP.beta_minus()
	# print ISTP.beta_plus()
	# print ISTP.alphas()
