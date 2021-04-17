import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import h5py 
import datetime
import copy
from calendar import monthrange
from cartopy.io.shapereader import Reader
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter,
                                LatitudeLocator)


def yfitdata(dataset):
    x=dataset['tepoch']
    x=np.array(x)
    xv=x-1469992008.04
    print(np.average(x),xv)

    a=7.8131227887e-08
    b=401.937969974
    c=-1.26751651147
    d=1.70977382721
    e=0.961498971762
    f=-0.151533144835
    g=-0.0514855180354
    h=0.064624982476
    i=-0.00939599602783
    j=-0.112117468199
    
    hb=a*xv+b
    omega_anual=2*np.pi/(365*24*60*60)

    hhb=(c*np.cos(omega_anual*x))+(d*np.sin(omega_anual*x))+(e*np.cos(2*omega_anual*x))+(f*np.sin(2*omega_anual*x))+(g*np.cos(3*omega_anual*x))+(h*np.sin(3*omega_anual*x))+(i*np.cos(4*omega_anual*x))+(j*np.sin(4*omega_anual*x))

    yfit=hb+hhb
    
    return yfit

###########################################################################################################################################################################

def compilealldata(lats,lons,filelist,footprintradio,verbosity):
    """ Takes a list of HDF files and scans for rows that meet the latitude, longitude, and targetname conditions """
    minlat = (min(lats)-footprintradio)
    maxlat = (max(lats)+footprintradio)
    minlon = (min(lons)-footprintradio)
    maxlon = (max(lons)+footprintradio)
    
    ilats = np.arange(int(maxlat)-int(minlat)+1,dtype=int)+int(minlat)
    ilons = np.arange(int(maxlon)-int(minlon)+1,dtype=int)+int(minlon)


    for ifile,filename in enumerate(filelist):
        datos = h5py.File(filename,'r')
        #print(ifile)
        for lat in ilats:
            for lon in ilons:
                #print(lat,lon)
                try:
                    gridname = '%iN%iW' % (lat,lon)
                    #print(gridname)
                    lt = datos[gridname]['lat'][()]
                    ln = datos[gridname]['lon'][()]
                    #print(ln)
                    mask = ((lt <= maxlat) & (lt >= minlat) & \
                           (ln <= maxlon) & (ln >= minlon))
                    #print(mask)
                    data_inter = datos[gridname][mask]
                    #print(data_inter)
                    
                    try:
                        data = np.append(data,data_inter)
                        #print('matrix exists')
                    except:
                        data = copy.copy(data_inter)
                        #print(data)
                        #print('new matrix')
                    if verbosity == 1:
                        print('Data shape for current loop is: ', data.shape)
                except:
                    print('exception: no dataset')
                #print 'data:',len(data)
                #print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])

        datos.close()
    #plt.plot(data['lon'],data['lat'],'.')
    #plt.show()
    return data

###########################################################################################################################################################################

def compilealldatat(lats,lons,filelist,footprintradio,verbosity):
    """ Takes a list of HDF files and scans for rows that meet the latitude, longitude, and targetname conditions """
    minlat = (min(lats)-footprintradio/110)
    maxlat = (max(lats)+footprintradio/110)
    minlon = (min(lons)-footprintradio/110)
    maxlon = (max(lons)+footprintradio/110)
    
    ilats = np.arange(int(maxlat)-int(minlat)+1,dtype=int)+int(minlat)
    ilons = np.arange(int(maxlon)-int(minlon)+1,dtype=int)+int(minlon)


    for ifile,filename in enumerate(filelist):
        datos = h5py.File(filename,'r')
        #print(ifile)
        for lat in ilats:
            for lon in ilons:
                #print(lat,lon)
                try:
                    gridname = '%iN%iW' % (lat,lon)
                    #print(gridname)
                    lt = datos[gridname]['lat'][()]
                    ln = datos[gridname]['lon'][()]
                    #print(ln)
                    mask = ((lt <= maxlat) & (lt >= minlat) & \
                           (ln <= maxlon) & (ln >= minlon))
                    #print(mask)
                    data_inter = datos[gridname][mask]
                    #print(data_inter)
                    
                    try:
                        data = np.append(data,data_inter)
                        #print('matrix exists')
                    except:
                        data = copy.copy(data_inter)
                        #print(data)
                        #print('new matrix')
                    if verbosity == 1:
                        print('Data shape for current loop is: ', data.shape)
                except:
                    print('exception: no dataset')
                #print 'data:',len(data)
                #print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])

        datos.close()
    #plt.plot(data['lon'],data['lat'],'.')
    #plt.show()
    return data


###########################################################################################################################################################################

def compilealldatatarget(lats,lons,filelist,footprintradio,targetname1,targetname2,verbosity):
    """ Takes a list of HDF files and scans for rows that meet the latitude, longitude, and targetname conditions """
    minlat = (min(lats)-footprintradio)
    maxlat = (max(lats)+footprintradio)
    minlon = (min(lons)-footprintradio)
    maxlon = (max(lons)+footprintradio)
    
    ilats = np.arange(int(maxlat)-int(minlat)+1,dtype=int)+int(minlat)
    ilons = np.arange(int(maxlon)-int(minlon)+1,dtype=int)+int(minlon)


    for ifile,filename in enumerate(filelist):
        datos = h5py.File(filename,'r')
        #print(ifile)
        for lat in ilats:
            for lon in ilons:
                #print(lat,lon)
                try:
                    gridname = '%iN%iW' % (lat,lon)
                    #print(gridname)
                    lt = datos[gridname]['lat'][()]
                    ln = datos[gridname]['lon'][()]
                    #print(ln)
                    tid = datos[gridname]['/Sounding/target_id'][()] 
                    tname = datos[gridname]['/Sounding/target_name'][()]
                    tidstr = tid.astype('U')
                    tnamestr = tname.astype('U')
                    mask = ((lt <= maxlat) & (lt >= minlat) & \
                           (ln <= maxlon) & (ln >= minlon) & \
                           ((tidstr == targetname1) | (tidstr == targetname2)))
                    #print(mask)
                    data_inter = datos[gridname][mask]
                    #print(data_inter)
                    
                    try:
                        data = np.append(data,data_inter)
                        #print('matrix exists')
                    except:
                        data = copy.copy(data_inter)
                        #print(data)
                        #print('new matrix')
                    if verbosity == 1:
                        print('Data shape for current loop is: ', data.shape)
                except:
                    print('exception: no dataset')
                #print 'data:',len(data)
                #print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])

        datos.close()
    #plt.plot(data['lon'],data['lat'],'.')
    #plt.show()
    return data


###########################################################################################################################################################################

def compilealldatatargett(lats,lons,filelist,footprintradio,targetname1,targetname2,verbosity):
    """ Takes a list of HDF files and scans for rows that meet the latitude, longitude, and targetname conditions """
    minlat = (min(lats)-footprintradio/110)
    maxlat = (max(lats)+footprintradio/110)
    minlon = (min(lons)-footprintradio/110)
    maxlon = (max(lons)+footprintradio/110)
    
    ilats = np.arange(int(maxlat)-int(minlat)+1,dtype=int)+int(minlat)
    ilons = np.arange(int(maxlon)-int(minlon)+1,dtype=int)+int(minlon)


    for ifile,filename in enumerate(filelist):
        datos = h5py.File(filename,'r')
        #print(ifile)
        for lat in ilats:
            for lon in ilons:
                #print(lat,lon)
                try:
                    gridname = '%iN%iW' % (lat,lon)
                    #print(gridname)
                    lt = datos[gridname]['lat'][()]
                    ln = datos[gridname]['lon'][()]
                    #print(ln)
                    tid = datos[gridname]['/Sounding/target_id'][()] 
                    tname = datos[gridname]['/Sounding/target_name'][()]
                    tidstr = tid.astype('U')
                    tnamestr = tname.astype('U')
                    mask = ((lt <= maxlat) & (lt >= minlat) & \
                           (ln <= maxlon) & (ln >= minlon) & \
                           ((tidstr == targetname1) | (tidstr == targetname2)))
                    #print(mask)
                    data_inter = datos[gridname][mask]
                    #print(data_inter)
                    
                    try:
                        data = np.append(data,data_inter)
                        #print('matrix exists')
                    except:
                        data = copy.copy(data_inter)
                        #print(data)
                        #print('new matrix')
                    if verbosity == 1:
                        print('Data shape for current loop is: ', data.shape)
                except:
                    print('exception: no dataset')
                #print 'data:',len(data)
                #print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])

        datos.close()
    #plt.plot(data['lon'],data['lat'],'.')
    #plt.show()
    return data

###########################################################################################################################################################################

def datacorrection_ml(dset,data,name,psurf,verbosity):
    """ Takes a dataset with reference pressure and xCO2 (e.g. Altzomoni) and the compiled data from compilealldatatarget,
    makes a copy of the compiled data, selects a day-1 to day+1 interval and corrects for mixing layer """
    
    cont = 0
    t0 = datetime.datetime.utcfromtimestamp(0.0)
    datacorr = copy.copy(data)
    if verbosity == 1:
        print("Datos seleccionados:")
        print(data[name])
    for ii,ele in enumerate(dset['year']):
        
        try:
            tmin = (datetime.datetime(year=dset['year'][ii],month=dset['month'][ii],day=dset['day'][ii]-1,hour=0,minute=0,second=0)-t0).total_seconds()
        except:
            maxdays = monthrange(dset['year'][ii],dset['month'][ii]-1)
            tmin = (datetime.datetime(year=dset['year'][ii],month=dset['month'][ii]-1,day=maxdays,hour=0,minute=0,second=0)-t0).total_seconds()
        try:
            tmax = (datetime.datetime(year=dset['year'][ii],month=dset['month'][ii],day=dset['day'][ii]+1,hour=23,minute=59,second=59)-t0).total_seconds() 
        except:
            tmax = (datetime.datetime(year=dset['year'][ii],month=dset['month'][ii]+1,day=1,hour=23,minute=59,second=59)-t0).total_seconds()
        
        mask = (tmin <= data['tepoch']) & (tmax >= data['tepoch'])
      
        
        datacorr[name][mask] = (data[name][mask]*data[psurf][mask] - dset['xCO2altz'][ii]*dset['psurfaltz'][ii])/(data[psurf][mask]-dset['psurfaltz'][ii])
                
        cont = cont + len(data[name][mask])
        #print(len(data[name][mask]))

    if verbosity == 1:
        print("Number of elements modified") 
        print(cont)
        print("original data")
        print(data[name].shape)
        print(data[name])
        print("corrected")
        print(datacorr[name].shape)
        print(datacorr[name])

    return datacorr
###########################################################################################################################################################################

def datacorrection2_ml(daylist,dset,data,name,psurf,verbosity):
    """ Takes a dataset with reference pressure and xCO2 (e.g. Altzomoni) and the compiled data from compilealldatatarget,
    makes a copy of the compiled data, selects a day-1 to day+1 interval and corrects for mixing layer
    This one is  redone of datacorrection to scan an auxiliary daylist file and mask the data according to each day, the reference will then
    be looked in the dataset under three different conditions (same day referece, same month reference, no reference)"""
    
    cont = 0
    t0 = datetime.datetime.utcfromtimestamp(0.0)
    datacorr = copy.copy(data)
    if verbosity == 1:
        print("Datos seleccionados:")
        print(data[name])
    for ii,ele in enumerate(daylist['tepoch']):
        tmin = ele
        tmax = ele + 86400        
        mask = (data['tepoch'] >= tmin) & (data['tepoch'] < tmax)
        datemin = datetime.datetime.utcfromtimestamp(tmin)

        if tmin in dset['tepoch']:
            if verbosity == 1:
                print(datemin.strftime("%b/%d"),"Hay referecia para el dia")
            indexdset = np.where(dset['tepoch']==tmin)
            xco2altz = dset['xco2avg'][indexdset]
        elif datemin.month in dset['month']:
            if verbosity == 1:
                print(datemin.strftime("%b/%d"),"Hay refencia para el mes")
            indexmonth = np.where(dset['month']==datemin.month)
            xco2altz = np.average(dset['xco2avg'][indexmonth])
        else:
            if verbosity == 1:
                print(datemin.strftime("%b/%d"),"No hay referencia, utilizando 408 para xco2altz (hardcoded)")
            xco2altz = 408.0


        datacorr[name][mask] = (data[name][mask]*data[psurf][mask] - xco2altz*639.25)/(data[psurf][mask]-639.25)

                
        cont = cont + len(data[name][mask])
        #print(len(data[name][mask]))

    if verbosity == 1:
        print("Number of elements modified") 
        print(cont)
        print("original data")
        print(data[name].shape)
        print(data[name])
        print("corrected")
        print(datacorr[name].shape)
        print(datacorr[name])

    return datacorr
###########################################################################################################################################################################

def datacorrection3_ml(daylist,dset,data,name,psurf,verbosity):
    """ Takes a dataset with reference pressure and xCO2 (e.g. Altzomoni) and the compiled data from compilealldatatarget,
    makes a copy of the compiled data, selects a day-1 to day+1 interval and corrects for mixing layer
    This one is  redone of datacorrection to scan an auxiliary daylist file and mask the data according to each day, the reference will then
    be looked in the dataset (with a xco2 fit for each day)"""
    
    datacorr = copy.copy(data)
    if verbosity == 1:
        print("Datos seleccionados:")
        print(data[name])
    yfit = yfitdata(data)
    datacorr[name] = (data[name]*data[psurf] - yfit*639.25)/(data[psurf]-639.25)

    if verbosity == 1:
        print("original data")
        print(data[name].shape)
        print(data[name])
        print("corrected")
        print(datacorr[name].shape)
        print(datacorr[name])

    return datacorr

###########################################################################################################################################################################

def matrixcorrection_ml(dset,data,psurfmat,name,psurf,lats,lons,verbosity):
    """ Takes a dataset with reference xCO2 (e.g. Altzomoni), the compiled data from compilealldatatarget, and a matrix with surface pressure
    makes a copy of the compiled data, selecst a day-1 to day+1 interval and corrects for mixing layer using the pressure from psurfmat. 
    
    Correction is done element-wise rather than array-wise because it has to index the cell from psurfmat that must be used, this is done by
    selecting the last index from arrays that comply a condition, creates a whole vector for each corrected element of the subset of data that 
    complies with the time interval and replaces the whole subset instead of replacing element wise due to python's limitation of fancy indexing
    """
    
    cont = 0
    t0 = datetime.datetime.utcfromtimestamp(0.0)
    datacorr = copy.copy(data)
    if verbosity == 1:
        print("Datos seleccionados:")
        print(data[name])

    for ii,ele in enumerate(dset['year']):
        
        try:
            tmin = (datetime.datetime(year=dset['year'][ii],month=dset['month'][ii],day=dset['day'][ii]-1,hour=0,minute=0,second=0)-t0).total_seconds()
        except:
            maxdays = monthrange(dset['year'][ii],dset['month'][ii]-1)
            tmin = (datetime.datetime(year=dset['year'][ii],month=dset['month'][ii]-1,day=maxdays,hour=0,minute=0,second=0)-t0).total_seconds()
        try:
            tmax = (datetime.datetime(year=dset['year'][ii],month=dset['month'][ii],day=dset['day'][ii]+1,hour=23,minute=59,second=59)-t0).total_seconds() 
        except:
            tmax = (datetime.datetime(year=dset['year'][ii],month=dset['month'][ii]+1,day=1,hour=23,minute=59,second=59)-t0).total_seconds()
        
        mask = (tmin <= data['tepoch']) & (tmax >= data['tepoch'])
        
        if verbosity == 1:
            print("Length of masked data:")
            print(len(datacorr[psurf][mask]))
        
        x=np.zeros((len(datacorr[psurf][mask])),dtype=np.float64)
        for jj,item in enumerate(datacorr[name][mask]):
            latcells = np.where(datacorr['lat'][mask][jj]>lats)
            loncells = np.where(datacorr['lon'][mask][jj]>lons)
            try:
                indexlon = loncells[0][-1]
            except:
                try:
                    indexlat = latcells[0][-1]
                    psurfval = psurfmat[0][indexlat]
                except:
                    psurfval = psurfmat[0][0]
            try:
                indexlat = latcells[0][-1]
                psurfval = psurfmat[indexlon][indexlat]
            except:
                psurfval = psurfmat[indexlon][0]
            
            a=data[name][mask][jj]*psurfval
            b=dset['xco2avg'][ii]*dset['psurfaltz'][ii]
            c=psurfval-dset['psurfaltz'][ii]
            datacorr[name][mask][jj] = (data[name][mask][jj]*psurfval - dset['xco2avg'][ii]*dset['psurfaltz'][ii])/(psurfval-dset['psurfaltz'][ii])
            x[jj] = (a-b)/c
            
            if verbosity == 1:
                print("x:",x[jj],"datacorr copy: ",datacorr[name][mask][jj], "original data:", data[name][mask][jj])
        
        datacorr[name][mask] = x
        
        cont = cont + len(data[name][mask])
        #print(len(data[name][mask]))
    if verbosity == 1:
        print("Number of elements modified") 
        print(cont)
        print("original data")
        print(data[name].shape)
        print(data[name])
        print("corrected")
        print(datacorr[name].shape)
        print(datacorr[name])
    return datacorr
###########################################################################################################################################################################

def matrixcorrection2_ml(daylist,dset,data,psurfmat,name,psurf,lats,lons,verbosity):
    """ Takes a dataset with reference xCO2 (e.g. Altzomoni), the compiled data from compilealldatatarget, and a matrix with surface pressure
    makes a copy of the compiled data, selecst a day-1 to day+1 interval and corrects for mixing layer using the pressure from psurfmat. 
    
    Correction is done element-wise rather than array-wise because it has to index the cell from psurfmat that must be used, this is done by
    selecting the last index from arrays that comply a condition, creates a whole vector for each corrected element of the subset of data that 
    complies with the time interval and replaces the whole subset instead of replacing element wise due to python's limitation of fancy indexing
    
    This one is  redone of matrixcorrection to scan an auxiliary daylist file and mask the data according to each day, the reference will then
    be looked in the dataset under three different conditions (same day referece, same month reference, no reference)"""

    
    cont = 0
    t0 = datetime.datetime.utcfromtimestamp(0.0)
    datacorr = copy.copy(data)
    if verbosity == 1:
        print("Datos seleccionados:")
        print(data[name])

    for ii,ele in enumerate(daylist['tepoch']):
        tmin = ele
        tmax = ele + 86400        
        mask = (data['tepoch'] >= tmin) & (data['tepoch'] < tmax)
        datemin = datetime.datetime.utcfromtimestamp(tmin)

        if tmin in dset['tepoch']:
            if verbosity == 1:
                print(datemin.strftime("%b/%d"),"Hay referecia para el dia")
            indexdset = np.where(dset['tepoch']==tmin)
            xco2altz = dset['xco2avg'][indexdset]
        elif datemin.month in dset['month']:
            if verbosity == 1:
                print(datemin.strftime("%b/%d"),"Hay refencia para el mes")
            indexmonth = np.where(dset['month']==datemin.month)
            xco2altz = np.average(dset['xco2avg'][indexmonth])
        else:
            if verbosity == 1:
                print(datemin.strftime("%b/%d"),"No hay referencia, utilizando 408 para xco2altz (hardcoded)")
            xco2altz = 408

        if verbosity == 1:
            print("Length of masked data:")
            print(len(datacorr[psurf][mask]))
        
        x=np.zeros((len(datacorr[psurf][mask])),dtype=np.float64)
        for jj,item in enumerate(datacorr[name][mask]):
            latcells = np.where(datacorr['lat'][mask][jj]>lats)
            loncells = np.where(datacorr['lon'][mask][jj]>lons)
            try:
                indexlon = loncells[0][-1]
            except:
                try:
                    indexlat = latcells[0][-1]
                    psurfval = psurfmat[0][indexlat]
                except:
                    psurfval = psurfmat[0][0]
            try:
                indexlat = latcells[0][-1]
                psurfval = psurfmat[indexlon][indexlat]
            except:
                psurfval = psurfmat[indexlon][0]
            
            a=data[name][mask][jj]*psurfval
            b=xco2altz*639.25
            c=psurfval-639.25
            datacorr[name][mask][jj] = (data[name][mask][jj]*psurfval - xco2altz*639.25)/(psurfval-639.25)
            x[jj] = (a-b)/c
            
            if verbosity == 1:
                print("x:",x[jj],"datacorr copy: ",datacorr[name][mask][jj], "original data:", data[name][mask][jj])
        
        datacorr[name][mask] = x
        
        cont = cont + len(data[name][mask])
        #print(len(data[name][mask]))
    if verbosity == 1:
        print("Number of elements modified") 
        print(cont)
        print("original data")
        print(data[name].shape)
        print(data[name])
        print("corrected")
        print(datacorr[name].shape)
        print(datacorr[name])
    return datacorr
###########################################################################################################################################################################

def matrixcorrection3_ml(daylist,dset,data,psurfmat,name,psurf,lats,lons,verbosity):
    """ Takes a dataset with reference xCO2 (e.g. Altzomoni), the compiled data from compilealldatatarget, and a matrix with surface pressure
    makes a copy of the compiled data, selecst a day-1 to day+1 interval and corrects for mixing layer using the pressure from psurfmat. 
    
    Correction is done element-wise rather than array-wise because it has to index the cell from psurfmat that must be used, this is done by
    selecting the last index from arrays that comply a condition, creates a whole vector for each corrected element of the subset of data that 
    complies with the time interval and replaces the whole subset instead of replacing element wise due to python's limitation of fancy indexing
    
    This one is  redone of datacorrection to scan an auxiliary daylist file and mask the data according to each day, the reference will then
    be looked in the dataset (with a xco2 fit for each day)"""

    
    datacorr = copy.copy(data)
    if verbosity == 1:
        print("Datos seleccionados:")
        print(data[name])

    yfit = yfitdata(data)
    x=np.zeros((len(datacorr[psurf])),dtype=np.float64)
    for jj,item in enumerate(datacorr[name]):
        latcells = np.where(datacorr['lat'][jj]>lats)
        loncells = np.where(datacorr['lon'][jj]>lons)
        print(datacorr['lon'][jj],datacorr['lat'][jj])
        try:
            indexlon = loncells[0][-1]
        except:
            try:
                indexlat = latcells[0][-1]
                psurfval = psurfmat[0][indexlat]
                print(min(lons),lats[indexlat])
            except:
                psurfval = psurfmat[0][0]
                print(min(lons),min(lats))
        try:
            indexlat = latcells[0][-1]
            psurfval = psurfmat[indexlon][indexlat]
            print(lons[indexlon],lats[indexlat])
        except:
            psurfval = psurfmat[indexlon][0]
            print(lons[indexlon],min(lats))
        a=data[name][jj]*psurfval
        b=yfit[jj]*639.25
        c=psurfval-639.25
        datacorr[name][jj] = (data[name][jj]*psurfval - yfit[jj]*639.25)/(psurfval-639.25)
        x[jj] = (a-b)/c
            
        if verbosity == 1:
            print(psurfval, yfit[jj])
            print("x:",x[jj],"datacorr copy: ",datacorr[name][jj], "original data:", data[name][jj])
        
    datacorr[name] = x
    
    if verbosity == 1:
        print("original data")
        print(data[name].shape)
        print(data[name])
        print("corrected")
        print(datacorr[name].shape)
        print(datacorr[name])
    return datacorr

###########################################################################################################################################################################

def makematrixfromcompileddata(lats,lons,data,name,footprintradio):
    
    """ Takes the latitude and longitude arrays for the n x m grid, the data to average, the variable name, the radius for oversampling and a return flag.
    Averages the data into an n x m grid by using an oversampling defined by a radius, generates the average matrix, matrix of elements in the grid, stddev and
    std error.
    
    Transposes the resulting matrices due to the way the matrices are generated by the 'for' loops """

    #print 'data:'
    #print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])
    matrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
    countermatrix=np.zeros((len(lons),len(lats)),dtype=int)
    stdmatrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
    for ilat,lat in enumerate(lats):
        for ilon,lon in enumerate(lons):
            distances2=np.array((data['lat']-lat)**2+(data['lon']-lon)**2)
            conditions=(footprintradio**2 > distances2)
            #print(conditions)
            index=np.where(conditions)[0]
            if len(index) > 0:
                vec=data[name][index]
                #print(vec)
                matrix[ilon,ilat]=np.average(vec)
                #print(matrix[ilon,ilat],ilon,ilat)
                stdmatrix[ilon,ilat]=np.std(vec)
                countermatrix[ilon,ilat]=len(vec)
            else:
                #print(matrix[ilon,ilat],ilon,ilat)
                pass
                #plt.plot(data['lon'],data['lat'],'b.')
                #plt.plot([lon],[lat],'ro')
                #plt.show()
    errmatrix=stdmatrix/np.sqrt(countermatrix)
    #print(matrix)
    
    return matrix.T,stdmatrix.T, errmatrix.T
    

###########################################################################################################################################################################

def matrixforinterval(tmin,tmax,dataorg,lats,lons,name,footprintradio):

    """ Takes a time interval (tmin, tmax), compiled data to crop, and return flag.
    Generates a set of matrices (average, stdev, sterror) from data that complies with the time interval filter """

    index=np.where(np.logical_and(tmin <= dataorg['tepoch'],tmax >= dataorg['tepoch']))[0]
    datacrop=dataorg[index]
    matrix,stdmat,errmat=makematrixfromcompileddata(lats,lons,datacrop,name,footprintradio)
    
    return matrix,stdmat, errmat
    
###########################################################################################################################################################################

def makehrmatrixfromcompileddata(lats,lons,data,name,footprintradio):
    
    """ Takes the latitude and longitude arrays for the n x m grid, the data to average, the variable name, the radius for oversampling and a return flag.
    Averages the data into an n x m grid by using an oversampling defined by a radius, generates the average matrix, matrix of elements in the grid, stddev and
    std error.
    
    Transposes the resulting matrices due to the way the matrices are generated by the 'for' loops """

    #print 'data:'
    #print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])
    matrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
    countermatrix=np.zeros((len(lons),len(lats)),dtype=int)
    stdmatrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
    hours = []
    for t in data['tepoch']:
        date = datetime.datetime.utcfromtimestamp(t)
        hours.append(date.hour)
    hrs = np.array(hours)
    for ilat,lat in enumerate(lats):
        for ilon,lon in enumerate(lons):
            distances2=np.array((data['lat']-lat)**2+(data['lon']-lon)**2)
            conditions=(footprintradio**2 > distances2)
            #print(conditions)
            index=np.where(conditions)[0]
            if len(index) > 0:
                vec=hrs[index]
                #print(vec)
                matrix[ilon,ilat]=np.average(vec)
                #print(matrix[ilon,ilat],ilon,ilat)
                stdmatrix[ilon,ilat]=np.std(vec)
                countermatrix[ilon,ilat]=len(vec)
            else:
                #print(matrix[ilon,ilat],ilon,ilat)
                pass
                #plt.plot(data['lon'],data['lat'],'b.')
                #plt.plot([lon],[lat],'ro')
                #plt.show()
    errmatrix=stdmatrix/np.sqrt(countermatrix)
    #print(matrix)
    
    return matrix.T,stdmatrix.T, errmatrix.T
    

###########################################################################################################################################################################

def hourlymatrixforinterval(tmin,tmax,dataorg,lats,lons,name,footprintradio):

    """ Takes a time interval (tmin, tmax), compiled data to crop, and return flag.
    Generates a set of matrices (average, stdev, sterror) from data that complies with the time interval filter """

    index=np.where(np.logical_and(tmin <= dataorg['tepoch'],tmax >= dataorg['tepoch']))[0]
    datacrop=dataorg[index]
    matrix,stdmat,errmat=makehrmatrixfromcompileddata(lats,lons,datacrop,name,footprintradio)
    
    return matrix,stdmat, errmat
    
###########################################################################################################################################################################

def countermatrixfromcompileddata(lats,lons,data,name,footprintradio):
    
    """ Takes the latitude and longitude arrays for the n x m grid, the data to average, the variable name, the radius for oversampling and a return flag.
    Averages the data into an n x m grid by using an oversampling defined by a radius, generates the average matrix, matrix of elements in the grid, stddev and
    std error.
    
    Transposes the resulting matrices due to the way the matrices are generated by the 'for' loops """

    #print 'data:'
    #print max(data['lat']), min(data['lat']), max(data['lon']), min(data['lon'])
    matrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
    countermatrix=np.zeros((len(lons),len(lats)),dtype=int)
    stdmatrix=np.zeros((len(lons),len(lats)),dtype=np.float64)
    hours = []
    for t in data['tepoch']:
        date = datetime.datetime.utcfromtimestamp(t)
        hours.append(date.hour)
    hrs = np.array(hours)
    for ilat,lat in enumerate(lats):
        for ilon,lon in enumerate(lons):
            distances2=np.array((data['lat']-lat)**2+(data['lon']-lon)**2)
            conditions=(footprintradio**2 > distances2)
            #print(conditions)
            index=np.where(conditions)[0]
            if len(index) > 0:
                vec=hrs[index]
                #print(vec)
                matrix[ilon,ilat]=np.average(vec)
                #print(matrix[ilon,ilat],ilon,ilat)
                stdmatrix[ilon,ilat]=np.std(vec)
                countermatrix[ilon,ilat]=len(vec)
            else:
                #print(matrix[ilon,ilat],ilon,ilat)
                pass
                #plt.plot(data['lon'],data['lat'],'b.')
                #plt.plot([lon],[lat],'ro')
                #plt.show()
    errmatrix=stdmatrix/np.sqrt(countermatrix)
    #print(matrix)
    
    return countermatrix.T
    

###########################################################################################################################################################################

def countermatrixforinterval(tmin,tmax,dataorg,lats,lons,name,footprintradio):

    """ Takes a time interval (tmin, tmax), compiled data to crop, and return flag.
    Generates a set of matrices (average, stdev, sterror) from data that complies with the time interval filter """

    index=np.where(np.logical_and(tmin <= dataorg['tepoch'],tmax >= dataorg['tepoch']))[0]
    datacrop=dataorg[index]
    countermat=countermatrixfromcompileddata(lats,lons,datacrop,name,footprintradio)
    
    return countermat    


###########################################################################################################################################################################

def mappernxm(latmin,latmax,lonmin,lonmax,lons,lats,rows,cols,titles,fontsize,ticksize,lvls,matrices,figname,savefig,cbarname):
    
    """ Mapper N x M (N=rows, M=columns): """

    width = cols * 10
    height = rows * 10
    matricesplot = matrices[:(rows*cols)]
    plt.figure(figsize=(width,height))   # figsize goes width, height
    for ii,matrix in enumerate(matricesplot):
        mxmap = plt.subplot(rows, cols, ii+1, projection=ccrs.PlateCarree()) 
        mxmap.add_feature(cfeature.LAND)
        mxmap.add_feature(cfeature.OCEAN)
        mxmap.add_feature(cfeature.COASTLINE)
        mxmap.add_feature(cfeature.BORDERS, linestyle=':')
        mxmap.add_feature(cfeature.LAKES, alpha=0.5)
        mxmap.add_feature(cfeature.RIVERS)
        mxmap.add_feature(cfeature.STATES.with_scale('10m'))
        mxmap.set_extent((lonmin,lonmax,latmin,latmax))
        #mxmap.xaxis.set_visible(True)
        #mxmap.yaxis.set_visible(True)
    
        glines = mxmap.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
        glines.xformatter = LongitudeFormatter()
        glines.yformatter = LatitudeFormatter()
        glines.top_labels = None
        glines.right_labels = None
        
        if (ii+cols) % cols != 0:
            glines.left_labels = None
        
        glines.xlabel_style = {'size': 16}
        glines.ylabel_style = {'size': 16}
    
        plt.contourf(lons, lats, matrices[ii], levels = lvls[ii],cmap=plt.get_cmap("jet"),transform=ccrs.PlateCarree())
        plt.title(titles[ii],size=fontsize)
        plt.xlabel("Longitude",size=fontsize)
        plt.ylabel("Latitude",size=fontsize)
        cbar= plt.colorbar(fraction=0.041, pad=.02)
        
        if (ii+1) % cols == 0:
            cbar.set_label(cbarname, size=fontsize, labelpad=0.8)
        
        cbar.ax.tick_params(labelsize=ticksize)

    #plt.tight_layout()
    if savefig==1:
        plt.savefig(figname)
    plt.show()

###########################################################################################################################################################################
