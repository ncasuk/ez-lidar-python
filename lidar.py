from netCDF4 import Dataset,Variable
import numpy as np
import time
import os.path
import struct
import glob
import json
import matplotlib.pyplot as plt
import scipy.misc
from collections import OrderedDict
from lidar_aux import aux_file
from lidar_raw import lidar_raw,rebuild_raw
from lidar_aux import aux_file

class lidar(object):
    """
    Class encapsulating lidar data
    """
    rawfolder=r"D:\Leosphere\EZAeroData"
    ncfolder=r"D:\NetCDF"
    
    def __init__(self,data=None,aux='HTTP',trigger=2054,**kwargs):
        """
        Initialise from 
            data: Path to Netcdf file, or path to raw data, or extant lidar data
            aux:  Auxilliary ( location ) data, as Core netcdf, or "Horace" text file, or HTTP to live aicraft data
            trigger: Point in raw data where laser fired...
        """
        for k in kwargs:
            if(k in dir(self)):
                self.__dict__[k]=kwargs[k]

        self.aux=aux
        if(not(os.path.isdir(self.ncfolder))):
            self.ncfolder=""
        if(type(data)==str):
            if(data.endswith(".nc")):
                self.datapath=data
                self.ncfolder=os.path.dirname(data)
                self.data=Dataset(data,**kwargs)
            elif(os.path.isdir(data)):
                self.datapath,self.data=self.create(data,filename=self.ncfolder,**kwargs)
                self.rawfolder=data
                self.add_raw()
            else:
                raise ValueError('Data not recognised:"'+data+'"')            
        elif(data):
            self.data=data
        else:
            date=self.aux.date.replace("_","-")
            self.rawfolder=os.path.join(self.rawfolder,date)
            self.datapath,self.data=self.create(self.rawfolder,filename=self.ncfolder,**kwargs)
            self.add_raw()

        self.variables=self.data.variables
        try:
            self.whereblind,=np.where(~self.data.variables["Blind_offset0"][:].mask)
        except AttributeError:
            self.whereblind=np.arange(len(self.data.variables['Blind_offset0'][:]))
        self.bind=np.full(self["Time"].shape,-1,dtype=int)
        add=0
        wb=self.whereblind[:]
        while(-1 in self.bind):
            wb=wb[(wb+add)<len(self.bind)]
            wb=wb[self.bind[wb+add]==-1]
            self.bind[wb+add]=wb
            add+=1
        self.profile=[lidar.getprofile(self.get_prof,self,chan=0),
                      lidar.getprofile(self.get_prof,self,chan=1)]
        self.range_corrected=[lidar.getprofile(self.get_rc,self,chan=0),
                              lidar.getprofile(self.get_rc,self,chan=1)]
        self.image=[lidar.getprofile(self.make_img,self,chan=0),
                    lidar.getprofile(self.make_img,self,chan=1)]
        self.curtain=[lidar.getprofile(self.make_curtain,self,chan=0),
                    lidar.getprofile(self.make_curtain,self,chan=1)]
        self.trigger=trigger

    @property
    def aux(self):
        return self._aux
    @aux.setter
    def aux(self,aux):
        try:
            for v in aux.columns[1:]:
                self.__setattr__(v,lidar.getprofile(self.get_aux,self,para=v))
            self._aux=aux
        except AttributeError:
            self.aux=aux_file(aux)
          
    @property
    def trigger(self):
        return self._trigger
    @trigger.setter
    def trigger(self,trigger):
        self._trigger=trigger
        self.distance=np.arange(-trigger,self['rawSignal_0'].shape[0]-trigger)*self.getncattr('RawResolution (m)')
        self.distance[:trigger]=np.nan
        
    
    def get_raw_indexes(self):
        b=self.whereblind
        for i in range(len(b)):            
            try:
                yield (b[i],b[i+1])
            except IndexError:
                yield (b[i],len(self.data.variables["Time"]))
            
    def rebuild_raw(self,folder=''):
        formats={'Altitude (m)':'{:.6f}','Longitude (deg)':'{:.6f}','Latitude (deg)':'{:.6f}',
                 'Pressure (hPa)':'{:.1f}','Temperature (degC)':'{:.1f}','Humidity (%)':'{:.1f}',
                 'AngleAzimuth':'{:.1f}','AngleZenith':'{:.1f}','AnglesNB AA':'{:.0f}',
                 'AnglesNB ZA':'{:.0f}','NumberOfShot':'{:.0f}','Wave length (nm)':'{:.0f}',
                 'PRF (Hz)':'{:.0f}','Port':'{:.0f}','LineSeparator':'{:.0f}',
                 'DecimalSeparator':'{:.0f}','NbOfProfilesPerFile':'{:.0f}','DataCodage':'{:.0f}',
                 'WritingPosition (byte)':'{:.0f}','NumberOfSignal':'{:.0f}',
                 'HeaderSize':'{:.0f}','ID ALS':'{:.0f}'}
        t=self['Time'][:]
        basetime=86400*(t[0]//86400)
        nprof=self.getncattr('NbOfProfilesPerFile')
        variables=[u'Altitude (m)', u'Longitude (deg)', u'Latitude (deg)', u'Pressure (hPa)', u'Temperature (degC)', u'AngleAzimuth', u'AngleZenith']
        for start,stop in self.get_raw_indexes():
            filename=time.strftime("_%Y-%m-%d_%H-%M-%S",time.gmtime(t[start]))+time.strftime("_%H-%M-%S.raw",time.gmtime(t[stop-1]))
            filename=os.path.join(folder,filename)
            with open(filename,"wb") as f:
                try:
                    nwrite=self.getncattr("WritingPosition (byte)")
                except AttributeError:
                    nwrite=-1
                f.write("[ConfigSoftware]\r\n")
                for att in self.ncattrs():
                    line=att+"="
                    try:
                        form=formats[att]
                    except KeyError:
                        form="{:.9f}"
                    data=self.getncattr(att)
                    if(att=='NbOfProfilesPerFile'):
                        data=(stop-start)
                    if(att=='DateRun'):
                        data=time.strftime('%Y-%m-%d',time.gmtime(t[start]))
                    if(str(data)==data):
                        line+=str(data)
                    else:
                        try:
                            line+="\t".join([form.format(d) for d in data])
                        except TypeError:
                            line+=form.format(data)
                    line+="\r\n"
                    line=(line.replace('deg',u'\xb0')).encode("latin-1")
                    f.write(line)
                    if(att=="VARIABLES"):
                        for v in variables:
                            if(v in formats):
                                form=formats[v].replace("{:","%").replace("}","")
                            else:
                                form="%.9f"
                            f.write((v.replace('deg',u'\xb0')).encode("latin-1")+"=")
                            try:
                                self.variables[v][start:stop].tofile(f,sep="\t",format=form)
                            except NotImplementedError:
                                self.variables[v][start].tofile(f,sep="\t",format=form)
                            f.write("\r\n")
                        
                for sect,prefix in [("InfoBlindRef","Blind_"),("infoRaw","Raw_")]:
                    f.write("["+sect+"]\r\n")
                    for v in self.variables:
                        if(v.startswith(prefix)):
                            line=v.replace(prefix,"")
                            if(line in formats):
                                form=formats[line]
                            else:
                                form="{:.9f}"
                            line+="="+form.format(self.variables[v][start])+"\r\n"
                            f.write(line)
                if(nwrite>0):
                    f.seek(nwrite)
                dim1=self.variables['rawSignal_0'].shape[0]
                f.write(self.write_dims((2,dim1)))
                self.variables['rawBlind_0'][:,start].astype(">i4").tofile(f,"")
                self.variables['rawBlind_1'][:,start].astype(">i4").tofile(f,"")
                for tx in range(start,stop):
                    f.write(self.write_time(t[tx]))
                    f.write(self.write_dims((4,dim1)))
                    for var in ['rawSignal_0','rawSignal_1','rawPhoton_0','rawPhoton_1']:
                        self.variables[var][:,tx].astype(">i4").tofile(f,"")
                
                    


    def merge_aux(self,aux=None):
        if(aux):
            self.aux=aux
        if(not(self.aux)):
            raise TypeError("No auxilliary data set")
        keys=[('Altitude (m)','ALT_GIN'),
              ('Longitude (deg)','LON_GIN'),
              ('Latitude (deg)','LAT_GIN'),
              ('Pressure (hPa)','PALT_RVS') ]
        for k,j in keys:
            ans=self.__getattribute__(j)[:]
            if(j=='PALT_RVS'):
                ans=heightpress(ans)
            self.data.variables[k][:]=ans
           

    def write_dims(self,dims):
        return struct.pack('>II',*dims) 
                   
    def write_time(self,t):
        return time.strftime('%H-%M-%S',time.gmtime(t))           

    def __getattr__(self,att):
        try:
            return super.__getattr__(self,att)
        except AttributeError:
            try:
                return self.data.__getattribute__(att)
            except AttributeError:
                return self.variables[att]
        
    def make_curtain(self,n,chan=0,heights='ALT_GIN',maxheight=0):
        try:
            h=self.__getattribute__(heights)[n]
            w=self['Raw_NumberOfSignal'][self.bind[n]]/self.getncattr('PRF (Hz)')
        except AttributeError:
            raise AttributeError('No height data')
        h/=1.5
        if(maxheight==0):
            maxheight=np.nanmax(h)
        if(maxheight!=maxheight):
            print(Warning("Invalid height - NaN"))
            maxheight=0
        mxh=maxheight
        if(mxh<1):
            mxh=1
        #im=np.zeros((mxh,self.nprof))
        rc=self.range_corrected[chan][n][:] # [self.trigger:,:]

        if(len(rc.shape)<2):
           rc=rc.reshape(rc.shape+(1,))
        
        im=np.zeros((mxh,rc.shape[1]))
        
        for prof in range(rc.shape[1]):
            h1=h[prof]
            if(h1==h1):  # Check if NaN
                h1=int(h1)
                if(h1>0):
                   im[mxh-h1:,prof]=rc[:h1,prof]
        return im[::-1,:]

    def make_img(self,n,chan=0,heights='ALT_GIN',vs='Time',maxheight=0,reduction=10):
        try:
            h=self.__getattribute__(heights)[n]
            w=self['Raw_NumberOfSignal'][self.bind[n]]/self.getncattr('PRF (Hz)')
        except AttributeError:
            raise AttributeError('No height data')
        h/=1.5
        if(maxheight==0):
            maxheight=np.nanmax(h)
        if(maxheight!=maxheight):
            print(Warning("Invalid height - NaN"))
            maxheight=0
        mxh=maxheight/reduction
        if(mxh<1):
            mxh=1
        #im=np.zeros((mxh,self.nprof))
        rc=self.range_corrected[chan][n][self.trigger:,:]
        x=self[vs][n]
        if(len(rc.shape)<2):
           rc=rc.reshape(rc.shape+(1,))
        
        x1=np.min(x)
        im=np.full((mxh,(np.max(x)-x1)+1),-1000)
        
        for prof in range(rc.shape[1]):
            h1=(h[prof]/reduction)
            if(h1==h1):  # Check if NaN
                h1=int(h1)
                if(h1>0):
                   im[mxh-h1:,(x[prof]-x1):(x[prof]-x1+w[prof])]=np.mean(rc[:h1*reduction,prof].reshape(h1,reduction,1),axis=1)
                else:
                   im[-1,(x[prof]-x1):(x[prof]-x1+w[prof])]=np.mean(rc[:h1,rof])
        return im
 
    def get_prof(self,n,chan=0):
        s=self['Raw_gain%i' % chan][self.bind[n]]*self['rawSignal_%i' % chan][:,n]/self['Raw_NumberOfSignal'][self.bind[n]]
        blind=self['Blind_gain%i' % chan][self.bind[n]]*self['rawBlind_%i' % chan][:,self.bind[n]]/self['Blind_NumberOfSignal'][self.bind[n]]
        s-=blind
        sky=np.mean(s[:self.trigger-5],axis=0)
        s=(s-sky)
        return s
        

    def get_rc(self,n,chan=0,div=166.0):
        s=self.get_prof(n,chan=chan)
        if(len(s.shape)>1):
            d=self.distance.reshape(self.distance.shape+(1,))
        else:
            d=self.distance
        return (s*(1.0+d*2.0/div)**2)[self.trigger:]

    def get_rc_test(self,n,chan=0,div=166.0):
        s=self.get_prof(n,chan=chan)
        if(len(s.shape)>1):
            d=self.distance.reshape(self.distance.shape+(1,))
        else:
            d=self.distance
        d=(1.0+d[self.trigger:]*2.0/div)
        l=d.shape[0]
        z=np.identity(l)
        m=((l-1)/2.0-np.abs((l-1)/2.0-np.arange(l))).astype(int)
        st=np.arange(l)-m
        ml=1+2*m
        import scipy.signal
        for i in range(l):
            sig=d[i]/((np.pi*2)**0.5)
            j=scipy.signal.gaussian(ml[i],sig)
            z[i,st[i]:st[i]+ml[i]]=j*d[i]*d[i]/np.sum(j)
        
        return (np.mat(s[self.trigger:].T)*z).T.A


    def get_rc_old(self,n,chan=0):
        s=self.get_prof(n,chan=chan)
        if(len(s.shape)>1):
            d=self.distance.reshape(self.distance.shape+(1,))
        else:
            d=self.distance
        return (s*d**2)[self.trigger:]


    def get_aux(self,n,para='PALT_RVS'):
        return self.aux.get_values(self['Time'][n],para=para)


    def get_img(self,n,chan=0):
        rc=self.get_rc(n,chan=chan)
        
        
    def make_jpg(self,chan,filename='',heights=[],maxheight=0,**kwargs):
        im=self.make_img(chan,heights=heights,maxheight=maxheight)
        if not(filename) or os.path.isdir(filename):
            fn=('lidar_%10.10i.jpg') % self.times[0]
            filename=os.path.join(filename,fn)
        plt.imsave(filename,im,**kwargs)
        return filename
        
    def __getitem__(self,item):
        return self.variables[item]

    class getprofile:
        def __init__(self,funct,data,**kwargs):
            self.funct=funct
            self.data=data
            self.kwargs=kwargs
        def __getitem__(self,n):
            return self.funct(n,**self.kwargs)
        def __len__(self):
            return len(self.data['Time'])


    def create(self,folder,**kwargs):
        fs=glob.glob(os.path.join(folder,'*.raw'))
        try:
            l=lidar_raw(sorted(fs)[0])
            print(l)
            ncpath,nc=l.createrawNetCDF(**kwargs)
            l.addData(nc)
            return ncpath,nc
        except IndexError:
            raise IOError("No Raw data in "+folder)
        

    def add_raw(self,folder="",files=[]):
        if(folder):
            self.rawfolder=folder
        if(not(files)):
            files=glob.glob(os.path.join(self.rawfolder,'*.raw'))
        last=self['Time'][-1]
        added=False
        for f in sorted(files):
            t=time.mktime(time.strptime(f[-32:-13]+"-UTC","%Y-%m-%d_%H-%M-%S-%Z"))
            if(t>last):     
                lidar_raw(f).addData(self)
                added=True
        return added
                 
    def rebuild_raw(self,folder=''):
        rebuild_raw(self,folder)        


    def createCurtainNC(self,filename='',fltno='XXXX',revision=0):
        """ Opens a raw netcdf file and creates 
        variables and attibutes.
        The global attributes are based on the "ConfigSoftware" header info
        The variables are from InfoBlindRef and infoRaw as well
        as the raw signal, photon count and blind reference values
        """
        date=time.strftime('%Y%m%d',time.gmtime(self['Time'][0]))
        if(not(filename) or os.path.isdir(filename)):
            fn=('metoffice-lidar_faam_'+date+'_r%1.1i_'+fltno+'_level1.nc') % revision
            filename=os.path.join(filename,fn)
        nc=Dataset(filename,"w",clobber=True)
        for att in self.ncattrs():
            nc.setncattr(att,self.getncattr(att))
        print('Extracting curtains...')
        curtain=[self.curtain[i][:] for i in range(2)]
        print('Create dataset...')
        nc.createDimension('Time',None)
        nc.createDimension('Altitude',curtain[0].shape[0])
        t=nc.createVariable('Time',float,('Time'))
        t.setncattr("units","seconds since 1970-01-01 00:00:00 +0000")
        t.setncattr("long_name","time of measurement")
        t.setncattr("standard_name","time")
        h=nc.createVariable('Altitude',float,('Altitude'))
        h.setncattr("units","metres")
        h.setncattr("long_name","Altitude of measurement")
        h.setncattr("standard_name","altitude")
        lat=nc.createVariable('Latitude',float,('Time'))
        lat.setncattr("units","deg")
        lat.setncattr("long_name","Latitude of measurement")
        lat.setncattr("standard_name","latitude")
        lon=nc.createVariable('Longitude',float,('Time'))
        lon.setncattr("units","deg")
        lon.setncattr("long_name","Longitude of measurement")
        lon.setncattr("standard_name","longitude")
        t[:]=self['Time']
        lat[:]=self['Latitude (deg)'][:]
        lon[:]=self['Longitude (deg)'][:]
        h[:]=np.arange(curtain[0].shape[0],dtype=float)*1.5
        v=[nc.createVariable('rangeCorrected_%1.1i' % i,float,('Altitude','Time'),zlib=True) for i in range(2)]
        v[0][:]=curtain[0]
        v[1][:]=curtain[1]
        

        return nc    
        
        
            

            
def create(folder,**kwargs):
    fs=glob.glob(os.path.join(folder,'*.raw'))
    first=True
    for f in sorted(fs):
        l=lidar_raw(f)
        if(first):
            ncpath,nc=l.createrawNetCDF(**kwargs)
            first=False
        l.addData(nc)
    return lidar(nc)
    
def pressheight(press,qnh=1013.25):
    if(qnh==1013.25):
        return (1-(press/1013.25)**0.190284)*44307.69396
    else:
        return pressheight(press)-pressheight(qnh)

def heightpress(height,qnh=1013.25):
    if(qnh==1013.25):
        return 1013.25*((1-(height/44307.69396))**(1/0.190284))
    else:
        return heightpress(height+pressheight(qnh))


