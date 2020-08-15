#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# importamos los paquetes
import numpy as np
import h5py
import glob
import format_oco3_target


# In[ ]:


# genera la lista de años y el arreglo de meses
years = [2019,2020]
meses = np.arange(1,13)


# In[ ]:


# genera la lista de archivos fuente y el nombre del archivo hdf que se creara, llama a la funcion de procesamiento
for year in years:
    for mes in meses:
        datapath ="d:\jupyter\oco3_mexico_city\h5\*%i%02i*.h5" % (year,mes)
        lista = glob.glob(datapath)
        lista.sort()
        h5filename='d:\jupyter\oco3_mexico_city\mexico\mexico_co3_target_%i%02i.h5' % (year,mes)

        if len(lista) != 0:
            format_oco3_target.processlist(lista,h5filename)
            #print len(lista)
        else:
            print("No hay archivos para %i-%02i" % (year,mes))

