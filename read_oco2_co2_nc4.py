import numpy as np
import h5py
import glob
import copy
import datetime as dt
import pytz

latmin=13.0
lonmin=-119
latmax=34.0
lonmax=-73

to=dt.datetime.utcfromtimestamp(0.0)

### STRING PARA FORMAR UN DATA TYPE QUE PERMITA LEER LOS CAMPOS RELEVANTES DE LOS ARCHIVOS DE OCO-2 AUTOMATICAMENTE, SE AGREGO AL PRINCIPIO tepoch, lat y lon POR QUE NO ESTA EN LOS ARCHIVOS 
str_type = \
"tepoch, float\lat, float\lon, float\/Retrieval/psurf, float\/Retrieval/psurf_apriori_b8, float\/Retrieval/xco2_raw, float\/Sounding/airmass, float\/Sounding/altitude, float\/Sounding/land_water_indicator, int\/co2_profile_apriori, float, 20\/date, int, 7\/latitude, float\/longitude, float\/pressure_levels, float, 20\/pressure_weight, float, 20\/solar_zenith_angle, float\/time, float\/xco2, float\/xco2_apriori, float\/xco2_averaging_kernel, float, 20\/xco2_uncertainty, float"

fields = []
for ele in str_type.split("\\"):
	trozos = tuple(ele.split(","))
	tupla = []
	for ii,ele in enumerate(trozos):
		if ii == 2:
			tupla.append(int(ele.strip()))
		else:
			tupla.append(ele.strip())
	fields.append(tuple(tupla))

oco2_co2_type=np.dtype(fields)
#print oco2_co2_type.names

def onebyone(mat_oco2,fh5):
	"""
	Recibe la matriz de mediciones filtradas por processlist y el archivo h5 donde se va a escribir las mediciones
	Escribe las mediciones por latitud y longitud en el archivo h5
	"""

	for ii,line in enumerate(mat_oco2):

		lat = line["/latitude"]
		lon = line["/longitude"]

		try:
			dset=fh5['%iN%iW' % (int(lat),int(lon))]
			n=dset.shape[0]
			dset.resize((n+1,))
			dset[n]=line
			fh5.flush()
		except:
			datalatlon=[line]
			try:
				dset=fh5.create_dataset('%iN%iW' % (int(lat),int(lon)),data=datalatlon,maxshape=(None,))

				fh5.flush()
			except:
				print 'nada works'

def processlist(lista,h5filename):
	"""
	Recibe una lista de archivos nc4 de oco2 para filtrarlos de acuerdo a latmin, latmax, lonmin y lonmax y el nombre del archivo h5 de salida
	Llama a onebyone
	"""

	try:
		fh5=h5py.File(h5filename,'r+')
		print 'file exist'
	except:
		fh5=h5py.File(h5filename,'w')
		print 'new file'

	cont = 0
	for filename in lista:
		#print filename
		datos=h5py.File(filename,'r')
		lat = datos["/latitude"][()]
		lon = datos["/longitude"][()]
		xco2_qf = datos["/xco2_quality_flag"][()]
		cond_latlon = ((lat > latmin) & (lat < latmax) & (lon > lonmin) & (lon < lonmax) & (xco2_qf == 0))	### FILTRO PARA LATITUD, LONGITUD Y xco2_quality_flag = 0 (GOOD)
		#print filename,len(lat),len(lat[cond_latlon])
		if len(lat[cond_latlon]) > 0:	### SOLO CUANDO AL APLICAR EL FILTRO SE OBTIENE AL MENOS UNA MEDICION
			mat_inter = np.empty(len(lat[cond_latlon]),dtype=oco2_co2_type)
			for name in oco2_co2_type.names[3:]:		### NO CONSIDERA tepoch, lat y lon POR QUE NO ESTAN EN EL HDF DE oco2
				mat_inter[name] = datos[name][()][cond_latlon]
			mat_inter["tepoch"] = np.array([(dt.datetime(year=huh[0],month=huh[1],day=huh[2],hour=huh[3],minute=huh[4],second=huh[5],microsecond=huh[6]*1000) - to).total_seconds() for huh in mat_inter["/date"]])
			mat_inter["lat"] = mat_inter["/latitude"]
			mat_inter["lon"] = mat_inter["/longitude"]
			#print "mat_inter",mat_inter.shape
			if cont == 0:
				mat_oco2 = copy.copy(mat_inter)	### SE CREA LA MATRIZ DONDE SE CONCATENAN LOS DIFERENTES DIAS AL INICIO, CUANDO EL CONTADOR ES 0
			else:
				mat_oco2 = np.concatenate((mat_oco2,mat_inter),axis=0)	### SE CONCATENA LOS VALORES FILTRADOS DE CADA DIA DEL MES A LA MATRIZ mat_oco2
			cont = cont+len(lat[cond_latlon])
		else:
			continue
			
		datos.close()
	#print "Contador",cont
	print "mat_oco2",h5filename,mat_oco2.shape
	onebyone(mat_oco2,fh5)	### SE MANDA LLAMAR onebyone
	#zlevels=np.arange(20)
	#fh5.attrs.create('zlevels', zlevels, dtype=zlevels.dtype )
	fh5.close()	### CIERRA EL ARCHIVO h5 DONDE SE ESCRIBE




