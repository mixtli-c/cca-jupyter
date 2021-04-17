##### format_co2.py:
##### Script para combine_co2.py que se encarga de generar una estructura con la informacion de las
##### variables que se extraeran el archivo fuente, extrae los datos del archivo fuente y crea un
##### archivo HDF separado por grid.

# importamos todas las librerias
import numpy as np
import h5py
import copy
import glob
import datetime as dt
import pytz

#####################################################################################################
### Esta seccion requiere configuracion por parte del usuario

# definimos el grid para la region
latmin = 13
latmax = 34
lonmin = -119
lonmax = -73

# nombre del archivo .txt con las variables
archivo_variables = 'f:\\gitCCA\\cca-jupyter\\data\\oco3lite_variables.txt'

# nombre de variables de latitud, longitud (e.g. /RetrievalGeometry/latitude) tiempo (e.g /time_tai93)
latname = "/latitude"
lonname = "/longitude"
timename = "/time"

# tipo de formato de tiempo en el archivo fuente: 0 - posix, 1 - tai , 2 - date
whattime = 0

# existe variable de quality flag? 0 - No, 1 - Si
qfexists = 1
qfname = "/xco2_quality_flag"

### NOTA IMPORTANTE: Las variable target_id y target_name necesitan convertirse a numpy unicode ('U')
### Si usas estas variables necesitas revisar la funcion processlist
### Fin de la seccion que requiere configuracion por parte del usuario
####################################################################################################

# inicializacion del tiempo para convertir TAI93 a epoch
t0=dt.datetime.utcfromtimestamp(0.0)                      # 1970-1-1 00:00:00
t93=dt.datetime(1993,1,1,0,0,0,0)                         # 1993-1-1 00:00:00
t93secs=(t93-t0).total_seconds()                          # segundos de 1970 a 1993

# obtiene las variables del archivo de texto (abre, copia a lista, y cierra)
vartxt = open(archivo_variables,'r')
variables = vartxt.read().replace('\n','')
vartxt.close()
# print(variables)

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
#print(fields)
read_tipos=np.dtype(fields)                # genera la estructura para crear la matriz
print(read_tipos)


def onebyone(mat_datos,fh5):
    """
    Recibe la matriz de mediciones filtradas por processlist y el archivo h5 donde se va a escribir las mediciones
    Escribe las mediciones por latitud y longitud en el archivo h5
    Original de read_oco2_co2_nc4.py 
    """

    for ii,line in enumerate(mat_datos):

        lat = line["lat"]                                  # toma el valor de latitud de la linea
        lon = line["lon"]                                 # toma el valor de longitud de la linea

        # try/except intenta abrir la seccion xN-yW del hdf creado y agrega la linea, si no hay seccion entonces crea una
        try:
            dset=fh5['%iN%iW' % (int(lat),int(lon))]             # toma el grid del hdf casero donde debería ir el dato
            n=dset.shape[0]                                      # obtiene el tamaño
            dset.resize((n+1,))                                  # le suma un elemento al tamaño 
            dset[n]=line                                         # le pone la linea al ultimo espacio 
            fh5.flush()                                          # escribe en disco 
        except:
            datalatlon=[line]                                    # crea una lista con los datos de la linea
            #dset=fh5.create_dataset('%iN%iW' % (int(lat),int(lon)),data=datalatlon,maxshape=(None,))
            #fh5.flush() 
            #print("except works")
            # try/except intenta crear un dataset, si fall entonces muestra un mensaje
            try:
                # crea un dataset xN-yW sin limite (None) tomando la parte entera de lat y lon con los datos de la linea
                dset=fh5.create_dataset('%iN%iW' % (int(lat),int(lon)),data=datalatlon,maxshape=(None,))
                fh5.flush()                                      # escribe en disco
            except:
                print('nada works')

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
        print(filename)
        datos=h5py.File(filename,'r')                              # abre el archivo hdf del target
        lat = datos[latname][()]                                # vector de latitud
        lon = datos[lonname][()]                                # vector de longitud
        
        if qfexists == 1:
            xco2_qf = datos[qfname][()]                  # vector de quality flag (sin usar para target)
        
        if qfexists == 1:
            # Lo siguiente es el filtro para latitud, longitud, y xco2 quality flag
            cond_latlon = ((lat > latmin) & (lat < latmax) & (lon > lonmin) & (lon < lonmax) & (xco2_qf == 0))
            
        if qfexists == 0:
            # Lo siguiente es el filtro para latitud, longitud
            cond_latlon = ((lat > latmin) & (lat < latmax) & (lon > lonmin) & (lon < lonmax)) 
            # print(filename,len(lat),len(lat[cond_latlon]))
        
        # el siguiente if/else crea la matriz mat_inter solo si se cumple la condicion
        if len(lat[cond_latlon]) > 0:                 
            
            mat_inter = np.empty(len(lat[cond_latlon]),dtype=read_tipos) # mat_inter es estructura target_hdf_tipos
            
            for name in read_tipos.names[3:]:                            # empieza desde 3 porque agregamos tepoch, lat, lon
                if (name=='/Sounding/target_id') or (name=='/Sounding/target_name'):
                    #print(name)
                    datosobj = datos[name][()][cond_latlon]
                    datosstr = datosobj.astype('<S28')
                    mat_inter[name] = datosstr
                    #print(datosstr)
                    #print(mat_inter[name])
                    #print('string')
                else: 
                    mat_inter[name] = datos[name][()][cond_latlon]
                    #print(mat_inter[name])
                    #print('notstring')
            
            # Los siguientes 3 IF seleccionan como se creara "tepoch" dependiendo del formato de tiempo del archivo fuente
            if whattime == 0:
                ### POSIX: El tiempo solamente se copia
                mat_inter["tepoch"] = mat_inter[timename]
            if whattime == 1:
                ### TAI93: Se suman los segundos de 1970 a 1993 al valor TAI93 del hdf para obtener epoch
                mat_inter["tepoch"] = np.array([huh+t93secs for huh in mat_inter[timename]])
            if whattime == 2:
                ### DATE: Se convierte el formato de fecha a segundos totales a partir de epoch
                mat_inter["tepoch"] = np.array([(dt.datetime(year=huh[0],month=huh[1],day=huh[2],hour=huh[3],minute=huh[4],second=huh[5],microsecond=huh[6]*1000) - to).total_seconds() for huh in mat_inter[timename]])
                
            mat_inter["lat"] = mat_inter[latname]                         # lat y lon adicionales para preservar un formato
            mat_inter["lon"] = mat_inter[lonname]                         # para todos los tipos de archivos
            #print(mat_inter["lat"])
            #print(mat_inter["lon"])
            #wait = input("PRESS ENTER TO CONTINUE.")
            if cont == 0:
                ### SE CREA LA MATRIZ DONDE SE CONCATENAN LOS DIFERENTES DIAS AL INICIO, CUANDO EL CONTADOR ES 0
                mat_datos = copy.copy(mat_inter)
            
            else:
                ### SE CONCATENA LOS VALORES FILTRADOS DE CADA DIA DEL MES A LA MATRIZ mat_target
                mat_datos = np.concatenate((mat_datos,mat_inter),axis=0)
            
            cont = cont+len(lat[cond_latlon])
            #wait = input("PRESS ENTER TO CONTINUE.")
        
        else:
            continue
            
        datos.close()
    #print("Contador",cont)
    #print(mat_datos)
    print("mat_datos",h5filename,mat_datos.shape)
    #wait = input("PRESS ENTER TO CONTINUE.")
    onebyone(mat_datos,fh5)  ### SE MANDA LLAMAR onebyone
    #zlevels=np.arange(20)
    #fh5.attrs.create('zlevels', zlevels, dtype=zlevels.dtype )
    fh5.close() ### CIERRA EL ARCHIVO h5 DONDE SE ESCRIBE
