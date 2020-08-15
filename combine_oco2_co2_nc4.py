import numpy as np
import h5py
import glob
import read_oco2_co2_nc4

years = [2017,2018,2019,2020]
meses = np.arange(1,13)

for year in years:
	for mes in meses:
		datapath ="LITE/oco2_LtCO2_%i%02i*" % (year-2000,mes)
		lista = glob.glob(datapath)
		lista.sort()
		h5filename='MEXICO/mexico_co2_oco2_%i%02i.h5'% (year,mes)
		
		if len(lista) != 0:
			read_oco2_co2_nc4.processlist(lista,h5filename)
			#print len(lista)
		else:
			print "No hay archivos para %i-%02i" % (year,mes)

		
