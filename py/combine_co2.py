##### combine_co2.py:
##### Este codigo fue tomado de EPR y modificado para hacerlo un poco mas general
##### Genera una lista de direcciones y nombres de archivo fuente
##### Genera una lista de nombres de archivos a escribir/generar
##### Llama a un script que toma estos nombres y genera los archivos HDF con formato EPR

# importamos los paquetes
import numpy as np
import h5py
import glob
import format_co2

#################################################################################################################
### Esta seccion requiere configuracion por el usuario

# genera la lista de anyos y el arreglo de meses
years = [2019, 2020]          
meses = np.arange(1,13)

# informacion sobre el path de los archivos
readpath = 'f:\\CCA\\NASAGESDISC\\OCO3L2Lite\\LITE\\'
readformat = '.nc4'

# use read_seedname si los nombres de la fuente tienen un formato establecido, de lo contrario use *
# NOTA: asegurese que el formato generado solo se repita en la seccion de fecha del archivo y no en otro identificador
# (e.g. orbita)
read_seedname = 'oco3_LtCO2_'

writepath = 'f:\\CCA\\NASAGESDISC\\OCO3L2Lite\\mexico\\'
write_seedname = 'mexico_oco3_lite_'

# yearformat: 0 si es dos digitos, 1 si es 4 digitos

yearformat = 0

### Fin de seccion que requiere configuracion por el usuario
#################################################################################################################

# testing
#for year in years:
#   for mes in meses:
#        if yearformat == 1:
#            filedate = read_seedname + '%i%02i*' % (year,mes)
#        else:
#            filedate = read_seedname + '%i%02i*' % (year-2000,mes)
#        datapath = readpath + filedate + readformat
#        #print(datapath)
#        lista = glob.glob(datapath)
#        lista.sort()
#        #print(lista)
#        h5filename = writepath + write_seedname + '%i%02i.h5' % (year,mes)
#        #print(h5filename)
# end of testing

# genera la lista de archivos fuente y el nombre del archivo hdf que se creara, llama a la funcion de procesamiento
for year in years:
    for mes in meses:
        if yearformat == 1:
            filedate = read_seedname + '%i%02i*' % (year,mes)
        else:
            filedate = read_seedname + '%i%02i*' % (year-2000,mes)
        datapath = readpath + filedate + readformat
        lista = glob.glob(datapath)
        lista.sort()
        h5filename = writepath + write_seedname + '%i%02i.h5' % (year,mes)

        if len(lista) != 0:
            format_co2.processlist(lista,h5filename)
            #print len(lista)
        else:
            print("No hay archivos para %i-%02i" % (year,mes))
