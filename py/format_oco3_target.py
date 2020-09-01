#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# importamos todas las librerias
import numpy as np
import h5py
import copy
import glob
import datetime as dt
import pytz


# In[ ]:


# definimos el grid para la region
latmin = 13
latmax = 34
lonmin = -119
lonmax = -73


# In[ ]:


# inicializacion del tiempo para convertir TAI93 a epoch
t0=dt.datetime.utcfromtimestamp(0.0)                      # 1970-1-1 00:00:00
t93=dt.datetime(1993,1,1,0,0,0,0)                         # 1993-1-1 00:00:00
t93secs=(t93-t0).total_seconds()                          # segundos de 1970 a 1993


# In[ ]:


# obtiene las variables del archivo de texto (abre, copia a lista, y cierra)
vartxt = open('targetvariables.txt','r')
variables = vartxt.read().replace('\n','')
vartxt.close()
# print(variables)


# In[ ]:


# tomado de read_oco2_co2_nc4.py de D3_SATELITE en EPR03 y modificado un poco

fields = []                                      # se inicializa la lista fields
for ele in variables.split("\\"):                # divide variables en elementos separados por linea (\)
    trozos = tuple(ele.split(","))               # genera tuples que contienen como elemento trozos de cada linea 
    #print(trozos)
    tupla = []                                   # inicia la lista que se agregara como tuple a fields
    for ii,ele in enumerate(trozos):            
        if ii == 2:
            tupla.append(int(ele.strip()))       # si trozos tiene tres elementos, agrega el tercero como int
        else:                                    # strip() elimina los espacios vacios
            tupla.append(ele.strip())            # los dos primeros elementos del tuple trozos se agregan a la lista tupla 
    fields.append(tuple(tupla))                  # agrega el la lista tupla como tuple a fields
# print(fields)
target_hdf_tipos=np.dtype(fields)                # genera la estructura para crear la matriz
# print(tipos)


# In[ ]:


def onebyone(mat_target,fh5):
    """
    Recibe la matriz de mediciones filtradas por processlist y el archivo h5 donde se va a escribir las mediciones
    Escribe las mediciones por latitud y longitud en el archivo h5
    Original de read_oco2_co2_nc4.py 
    """

    for ii,line in enumerate(mat_target):

        lat = line["/latitude"]                                  # toma el valor de latitud de la linea
        lon = line["/longitude"]                                 # toma el valor de longitud de la linea

        # try/except intenta abrir la seccion xN-yW del hdf creado y agrega la linea, si no hay seccion entonces crea una
        try:
            dset=fh5['%iN%iW' % (int(lat),int(lon))]             # toma el grid del hdf casero donde debería ir el dato
            n=dset.shape[0]                                      # obtiene el tamaño
            dset.resize((n+1,))                                  # le suma un elemento al tamaño 
            dset[n]=line                                         # le pone la linea al ultimo espacio 
            fh5.flush()                                          # escribe en disco 
        except:
            datalatlon=[line]                                    # crea una lista con los datos de la linea
            # try/except intenta crear un dataset, si fall entonces muestra un mensaje
            try:
                # crea un dataset xN-yW sin limite (None) tomando la parte entera de lat y lon con los datos de la linea
                dset=fh5.create_dataset('%iN%iW' % (int(lat),int(lon)),data=datalatlon,maxshape=(None,))
                fh5.flush()                                      # escribe en disco
            except:
                print('nada works')


# In[ ]:


def processlist(lista,h5filename):
    """
    Recibe una lista de archivos hdf del oco3 target para filtrarlos de acuerdo a 
    latmin, latmax, lonmin y lonmax y el nombre del archivo h5 de salida
    Llama a onebyone
    Modificada por Mixtli Campos, Ago-2020
    Original de read_oco2_co2_nc4.py 
    """

    # esta seccion try/except intenta abrir un archivo hdf existente y si no puede crea uno nuevo
    try:
        fh5=h5py.File(h5filename,'r+')
        print('file exist')
    except:
        fh5=h5py.File(h5filename,'w')
        print('new file')

    # este contador ayuda a ver cuantas lineas se agregan
    cont = 0
    
    # esta seccion se encarga de generar una matriz para que la funcion onebyone escriba el hdf 
    for filename in lista:
        #print(filename)
        datos=h5py.File(filename,'r')                              # abre el archivo hdf del target
        lat = datos["/latitude"][()]                                # vector de latitud
        lon = datos["/longitude"][()]                               # vector de longitud
        #xco2_qf = datos["/xco2_quality_flag"][()]                  # vector de quality flag (sin usar para target)
        
        # Lo siguiente es el filtro para latitud, longitud
        cond_latlon = ((lat > latmin) & (lat < latmax) & (lon > lonmin) & (lon < lonmax)) 
        #print(filename,len(lat),len(lat[cond_latlon]))
        
        # el siguiente if/else crea la matriz mat_inter solo si se cumple la condicion
        if len(lat[cond_latlon]) > 0:                 
            mat_inter = np.empty(len(lat[cond_latlon]),dtype=target_hdf_tipos) # mat_inter es estructura target_hdf_tipos
            for name in target_hdf_tipos.names[1:]:                            # empieza desde 1 porque agregamos tepoch
                mat_inter[name] = datos[name][()][cond_latlon]
            ### Se suman los segundos de 1970 a 1993 al valor TAI93 del hdf para obtener epoch
            mat_inter["tepoch"] = np.array([huh+t93secs for huh in mat_inter["/time_tai93"]])
            #mat_inter["lat"] = mat_inter["/latitude"]                         # lat y lon adicionales no se necesitan
            #mat_inter["lon"] = mat_inter["/longitude"]                        # para los target hdf
            #print("mat_inter",mat_inter.shape)
            if cont == 0:
                ### SE CREA LA MATRIZ DONDE SE CONCATENAN LOS DIFERENTES DIAS AL INICIO, CUANDO EL CONTADOR ES 0
                mat_target = copy.copy(mat_inter)
            else:
                ### SE CONCATENA LOS VALORES FILTRADOS DE CADA DIA DEL MES A LA MATRIZ mat_target
                mat_target = np.concatenate((mat_target,mat_inter),axis=0)
            cont = cont+len(lat[cond_latlon])
        else:
            continue
            
        datos.close()
    #print("Contador",cont)
    print("mat_target",h5filename,mat_target.shape)
    onebyone(mat_target,fh5)  ### SE MANDA LLAMAR onebyone
    #zlevels=np.arange(20)
    #fh5.attrs.create('zlevels', zlevels, dtype=zlevels.dtype )
    fh5.close() ### CIERRA EL ARCHIVO h5 DONDE SE ESCRIBE

