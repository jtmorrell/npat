from npat import Spectrum, Calibration


sp = Spectrum('eu_calib_7cm.Spe')
sp.meta = {'istp':['152EU'], 'A0':3.7E4, 'ref_date':None}
sp.fit_config = {'xrays':True, 'E_min':20.0}
# cb = Calibration()
# cb.calibrate([sp], auto_calibrate=True)
# cb.plot()
sp.auto_calibrate()
sp.cb.plot()
sp.summarize()
sp.plot()

# 
# sp = Spectrum('La01.Spe')
# sp.meta = {'istp': ['135CE','134CE','134LA','137CE','137CEm','139CE','133BA','133BAm','132CS','135LA','22NA','24NA']}
# sp = Spectrum('/home/jmorrell/Documents/Radium_Bernstein_Oct2018/data/count_room_data/experiment/225Ac_separated_20cm_000.Spe')
# sp.meta = {'istp':['221FR','217AT','213BI','213PO','209PB','224RA','220RN','216PO','212PB','212BI','214BI']}
# sp.fit_config = {'skew_fit':True}
# sp.auto_calibrate()
# sp = Spectrum('iridium_7cm.Spe')
# sp.meta = {'istp':['192IR','194IR','193IRm']}




