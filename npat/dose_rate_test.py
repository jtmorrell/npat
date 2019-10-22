import numpy as np

def dose_rate(activity=1.0, distance=30.0, units='R'):

	def beta2(E_MeV, m_amu):
		return 1.0-(1.0/(1.0+(E_MeV/m_amu))**2)

	def e_range(E_MeV):
		dEdx = lambda b2, t: 0.17*((np.log(3.61E5*t*np.sqrt(t+2)))+(0.5*(1-b2)*(1+t**2/8-(2*t+1)*np.log(2)))-4.312)/b2
		E = np.linspace(1E-9, E_MeV, 100)
		return np.trapz(1.0/((1.0+(E*7.22/800.0))*dEdx(beta2(E, 0.5109), E/0.5109)), E)

	def pos_range(E_MeV):
		dEdx = lambda b2, t: 0.17*((np.log(3.61E5*t*np.sqrt(t+2)))+(np.log(2)-(b2/24)*(23+14/(t+2)+10/(t+2)**2+4/(t+2)**3))-4.312)/b2
		E = np.linspace(1E-9, E_MeV, 100)
		return np.trapz(1.0/((1.0+(E*7.22/800.0))*dEdx(beta2(E, 0.5109), E/0.5109)), E)

	def alpha_range(E_MeV):
		dEdx = lambda b2: 0.17*4.0*((np.log(1.02E6*b2/(1.0-b2))-b2)-4.312)/b2
		return np.trapz(1.0/dEdx(beta2(np.linspace(1E-9, E_MeV, 100), 3.7284E3)), np.linspace(1E-9, E_MeV, 100))

