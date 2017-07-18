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


class lidar_raw:
    """
    Class for reading in raw lidar files,
    and writing raw netcdf lidar data
    """
    def __init__(self,filename):
        self.filename=filename
        self.file=open(filename,'rb')
        self.header=OrderedDict()
        self.blind_smoothing=None
        section=''
        description=None
        lines=0
        headerlines=1000
        while(lines<headerlines):
            z=self.file.readline().decode('latin1').replace(u'\xb0','deg').strip().split("=")
            lines+=1
            if(len(z)>1):
                d=z[1].split("\t")
                for i in range(len(d)):
                    try:
                        d[i]=int(d[i])
                    except ValueError:
                        try:
                            d[i]=float(d[i])
                        except ValueError:
                            pass
                if len(d)==1:
                    d=d[0]
                if z[0]=="HeaderSize": headerlines=d
                if z[0]==z[0].upper():
                    if z[0]=="VARIABLES":
                        self.header[section][z[0]]=d
                        section="VARIABLES"
                        self.header[section]=OrderedDict()
                    else:
                        section="ConfigSoftware"
                        self.header[section][z[0]]=d
                else:                  
                    self.header[section][z[0]]=d
            else:
                if z[0].startswith('['):
                    section=z[0][1:-1]
                    self.header[section]=OrderedDict()
                    if(description):
                        self.header[section]["Description"]=description
                        description=None
                elif section=='':
                    description=z[0] # Old flight !!
        version=self.header['ConfigSoftware']['Version'].split('.')
        version=10000*float(version[0])+100*float(version[1])+float(version[2])
        try:
            self.pos=self.header['ConfigSoftware']['WritingPosition (byte)']
        except KeyError:
            self.pos=self.file.tell()
            #self.header['ConfigSoftware']['WritingPosition (byte)']=self.pos
            if(version>=11200):
                print("Warning: Versions >= 1.12.0 should include Writing Position")
            else:
                print("Info: Version < 1.12.0 ")
            
            
        self.it = self.header['ConfigSoftware']['NumberOfShot'] / self.header['ConfigSoftware']['PRF (Hz)']
        self.nprof=self.header['ConfigSoftware']['NbOfProfilesPerFile'] 
        self.get_basetime()
        self.file.seek(self.pos)
        (self.bdims,self.blindraw)=self.get_raw()
        #pos=self.pos+8+self.bdims[0]*self.bdims[1]*4
        #self.file.seek(pos)
        self.times=[self.get_time()]  
        (self.dims,raw)=self.get_raw()
        self.raw=np.empty((self.nprof,)+self.dims,dtype=raw.dtype)
        self.raw[0]=raw
        for n in range(self.nprof):
            if(n!=0):
                #pos=pos+16+self.dims[0]*self.dims[1]*4
                #self.file.seek(pos)
                self.times.append(self.get_time())
                (self.dims,raw)=self.get_raw()
                self.raw[n]=raw
        self.file.close()

    def get_time(self):
        t=self.file.read(8)
        return self.basetime+int(t[:2])*3600+int(t[3:5])*60+int(t[6:])
        
    def get_dims(self):
        return struct.unpack('>II',self.file.read(8))

    def get_raw(self):
        dims=self.get_dims()
        data=np.fromfile(self.file,count=dims[0]*dims[1],dtype='>i4').reshape(dims)
        return dims,data

    def getdate(self):
        dte=self.header['ConfigSoftware']['DateRun']
        return dte[:4]+dte[5:7]+dte[8:]

    def get_basetime(self):
        self.basetime=time.mktime(time.strptime(self.header['ConfigSoftware']['DateRun']+"-UTC","%Y-%m-%d-%Z"))

    def createrawNetCDF(self,filename='',fltno='XXXX',revision=0,**kwargs):
        if(not(filename) or os.path.isdir(filename)):
            fn=('metoffice-lidar_faam_'+self.getdate()+'_r%1.1i_'+fltno+'_raw.nc') % revision
            filename=os.path.join(filename,fn)
        return filename,self.openrawNetCDF(Dataset(filename,"w",clobber=True))

    def openrawNetCDF(self,nc):
        """ Opens a raw netcdf file and creates 
        variables and attibutes.
        The global attributes are based on the "ConfigSoftware" header info
        The variables are from InfoBlindRef and infoRaw as well
        as the raw signal, photon count and blind reference values
        """
        for att in self.header["ConfigSoftware"]:
            nc.setncattr(att,self.header["ConfigSoftware"][att])
        nc.createDimension('Time',None)
        nc.createDimension('Range',self.dims[1])
        t=nc.createVariable('Time',float,('Time'))
        t.setncattr("units","seconds since 1970-01-01 00:00:00 +0000")
        t.setncattr("long_name","time of measurement")
        t.setncattr("standard_name","time")
        
        for sect,prefix in [("InfoBlindRef","Blind_"),("infoRaw","Raw_"),("VARIABLES","")]:
            for att in self.header[sect]:
                nc.createVariable(prefix+att,float,('Time'))
        
        for i in range(2):
            nc.createVariable('rawSignal_%1.1i' % i,'i4',('Range','Time'),zlib=True)

        for i in range(2,self.dims[0]):
            nc.createVariable('rawPhoton_%1.1i' % (i-2),'i4',('Range','Time'),zlib=True)


        for i in range(self.bdims[0]):
            nc.createVariable('rawBlind_%1.1i' % i,'i4',('Range','Time'),zlib=True)

        return nc    
        
        
    def addData(self,nc):
        """
        Add data to a netcdf file
        """
        n=len(nc.variables['Time'])
        n2=n+len(self.times)
        nc.variables['Time'][n:]=self.times
        for i in range(2):
            nc.variables['rawSignal_%1.1i' % i][:,n:]=self.raw.T[:,i,:]

        for i in range(2,self.dims[0]):
            nc.variables['rawPhoton_%1.1i' % (i-2)][:,n:]=self.raw.T[:,i,:]

        """
        #Fill all blind refs
        for i in range(self.bdims[0]):
            for nn in range(n,n2):
                nc.variables['rawBlind_%1.1i' % i][:,nn]=self.blindraw[i]
        """
        #fill only relavent blindrefs up
        for i in range(self.bdims[0]):
            nc.variables['rawBlind_%1.1i' % i][:,n]=self.blindraw[i]
        
        
        for sect,prefix in [("InfoBlindRef","Blind_"),("infoRaw","Raw_"),("VARIABLES","")]:
            for att in self.header[sect]:
                try:
                    if(len(self.header[sect][att])>1):
                        nc.variables[prefix+att][n:]=self.header[sect][att]
                    else:
                        nc.variables[prefix+att][n]=self.header[sect][att]
                except TypeError:       
                    nc.variables[prefix+att][n]=self.header[sect][att]
            

def rebuild_raw(ncdata,folder=''):
    """
    Function to write lidar data as
    leosphere style raw files
    """
    formats={'Altitude (m)':'{:.6f}','Longitude (deg)':'{:.6f}','Latitude (deg)':'{:.6f}',
             'Pressure (hPa)':'{:.1f}','Temperature (degC)':'{:.1f}','Humidity (%)':'{:.1f}',
             'AngleAzimuth':'{:.1f}','AngleZenith':'{:.1f}','AnglesNB AA':'{:.0f}',
             'AnglesNB ZA':'{:.0f}','NumberOfShot':'{:.0f}','Wave length (nm)':'{:.0f}',
             'PRF (Hz)':'{:.0f}','Port':'{:.0f}','LineSeparator':'{:.0f}',
             'DecimalSeparator':'{:.0f}','NbOfProfilesPerFile':'{:.0f}','DataCodage':'{:.0f}',
             'WritingPosition (byte)':'{:.0f}','NumberOfSignal':'{:.0f}',
             'HeaderSize':'{:.0f}','ID ALS':'{:.0f}'}
    t=ncdata['Time'][:]
    basetime=86400*(t[0]//86400)
    nprof=ncdata.getncattr('NbOfProfilesPerFile')
    variables=[u'Altitude (m)', u'Longitude (deg)', u'Latitude (deg)', u'Pressure (hPa)', u'Temperature (degC)', u'AngleAzimuth', u'AngleZenith']
    for start,stop in ncdata.get_raw_indexes():
        filename=time.strftime("_%Y-%m-%d_%H-%M-%S",time.gmtime(t[start]))+time.strftime("_%H-%M-%S.raw",time.gmtime(t[stop-1]))
        filename=os.path.join(folder,filename)
        with open(filename,"wb") as f:
            try:
                nwrite=ncdata.getncattr("WritingPosition (byte)")
            except AttributeError:
                nwrite=-1
            f.write("[ConfigSoftware]\r\n")
            for att in ncdata.ncattrs():
                line=att+"="
                try:
                    form=formats[att]
                except KeyError:
                    form="{:.9f}"
                data=ncdata.getncattr(att)
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
                            ncdata.variables[v][start:stop].tofile(f,sep="\t",format=form)
                        except NotImplementedError:
                            ncdata.variables[v][start].tofile(f,sep="\t",format=form)
                        f.write("\r\n")
                    
            for sect,prefix in [("InfoBlindRef","Blind_"),("infoRaw","Raw_")]:
                f.write("["+sect+"]\r\n")
                for v in ncdata.variables:
                    if(v.startswith(prefix)):
                        line=v.replace(prefix,"")
                        if(line in formats):
                            form=formats[line]
                        else:
                            form="{:.9f}"
                        line+="="+form.format(ncdata.variables[v][start])+"\r\n"
                        f.write(line)
            if(nwrite>0):
                f.seek(nwrite)
            dim1=ncdata.variables['rawSignal_0'].shape[0]
            f.write(write_dims((2,dim1)))
            ncdata.variables['rawBlind_0'][:,start].astype(">i4").tofile(f,"")
            ncdata.variables['rawBlind_1'][:,start].astype(">i4").tofile(f,"")
            for tx in range(start,stop):
                f.write(write_time(t[tx]))
                f.write(write_dims((4,dim1)))
                for var in ['rawSignal_0','rawSignal_1','rawPhoton_0','rawPhoton_1']:
                    ncdata.variables[var][:,tx].astype(">i4").tofile(f,"")


def write_dims(dims):
    return struct.pack('>II',*dims) 
               
def write_time(t):
    return time.strftime('%H-%M-%S',time.gmtime(t))           

