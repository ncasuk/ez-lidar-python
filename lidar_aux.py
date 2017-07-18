import threading
import multiprocessing
import urllib2
import cookielib
import time
import numpy as np
from netCDF4 import Dataset
import time
import os.path
import socket
import struct
import json

class aux_file(object):
    """
    Class to encapsulate the auxilliary data required by the lidar ( basically position data )
    """
    columns=['Time','PALT_RVS','ALT_GIN','LON_GIN','LAT_GIN','PTCH_GIN','ROLL_GIN','HGT_RADR']
    paras=['time_since_midnight','pressure_height_m','gin_altitude','gin_longitude','gin_latitude','gin_pitch','gin_roll','radar_height']
    format="%7.1f %10.3f %10.3f %11.5f %11.5f %7.2f %7.2f %10.3f"
    default_tank="192.168.101.108"
    file_prefix="horace_"
    folder=r"C:\Horace"
    dataurl="/live/jsondata?"
    initial_data="para=flight_number&para=time_since_midnight&para=utc_time"
    timeout=2
    flt_no='XXXX'
    def __init__(self,source='HTTP',path=r"C:\Horace",**kwargs):
        """
        Reads or creates a "Horace" text file for compatibility with IDL lidar code
        this is optionally kept up to date with live data from decades, or 
        recreated post flight from core netCDF
        """ 
        for k in kwargs:
            if(k in dir(self)):
                self.__dict__[k]=kwargs[k]
        self.dtype=zip(self.columns,['f8']*len(self.columns))

        if(source.upper().startswith("HTTP")):
            if(source.upper()=="HTTP"):
                self.tank=self.default_tank
            else:
                s=source.split(":")
                self.tank=s[1].replace("/","")
            self.url="http://"+self.tank+self.dataurl
            self.paraurl=self.url+"para="+"&para=".join(self.paras)
            self.timeout=2
            self.initialise()                
            self.filename=os.path.join(path,self.file_prefix+self.date+".dat")
            try:
                self.read()
            except IOError:
                pass
            
        elif(source.endswith(".nc")):
            if(path):
                self.folder=path
            self.read_nc(source,**kwargs)
        elif(source.endswith(".dat")):
            self.filename=source
            self.folder=os.path.dirname(source)
            self.date=source
            self.read()


    def initialise(self):
        try:
            response = urllib2.urlopen(self.url+self.initial_data,timeout=self.timeout)
            js=response.read()
            if(js):
                dat=json.loads(js)
                self.basetime=dat['utc_time'][0]-dat['time_since_midnight'][0]
                self.flt_no=dat['flight_number'][0]
            else:
                self.basetime=None
        except(urllib2.HTTPError):
            self.basetime=None

    @property
    def date(self):
        return self._date
    @date.setter
    def date(self,d):
        if(d.endswith(".dat")):
            self._date=d[-14:-4]
        elif(d.startswith("seconds since")):
            self._date=d[14:24].replace("-","_")
        self._basetime=time.mktime(time.strptime(self._date+"UTC","%Y_%m_%d%Z"))

    @property
    def basetime(self):
        return self._basetime
    @basetime.setter
    def basetime(self,b):
        try:
            t=time.gmtime(b)
            self._basetime=time.mktime(t)
            self._date=time.strftime("%Y_%m_%d",t)
        except TypeError:
            pass

    @property
    def data(self):
        return self._data
    @data.setter
    def data(self,d):
        self._data=d
        self.times=d['Time']+self.basetime

    def start(self):
        self.thread=self.HTTP_thread(self)
        self.thread.start()
    def stop(self):
        self.thread.running=False
        
    class HTTP_thread(threading.Thread):
        def __init__(self,caller):
            self.caller=caller
            threading.Thread.__init__(self)
        def run(self):   
            self.running=True
            while self.running:
                time.sleep(1.0)
                self.caller.add_latest()
                      
        
    def add_latest(self):
        try:
            last=int(self.times[-1]+1)
        except AttributeError:
            last=0
        #print(self.paraurl+"&frm={:d}".format(last))
        response = urllib2.urlopen(self.paraurl+"&frm={:d}".format(last),timeout=self.timeout)
        js=response.read()
        if(js):
            dat=json.loads(js)
            if(len(dat['time_since_midnight'])):
                data=np.empty(len(dat['time_since_midnight']),dtype=self.dtype)
                for c,p in zip(self.columns,self.paras):
                    data[c]=dat[p]
                if(self.file_prefix):
                    with open(self.filename,"a") as f:
                        np.savetxt(f,data,fmt=self.format)
                try:
                    self.data=np.append(self.data,data)
                except AttributeError:
                    self.data=data
        

    def read_nc(self,path,fill=np.nan,**kwargs):
        data=Dataset(path)
        self.date=data.variables['Time'].units
        try:
            self.flt_no=data.getncattr('FLIGHT').upper()
        except:
            self.flt_no=os.path.basename(path)[27:31].upper()
                           
        self.filename=os.path.join(self.folder,self.file_prefix+self.date+".dat")
        d=np.empty(data.variables['Time'].shape,dtype=self.dtype)
        for c in self.columns:
            if(len(data.variables[c].shape)>1):
                try:
                    d[c]=data.variables[c][:,0].filled(fill)
                except AttributeError:
                    d[c]=data.variables[c][:,0]
                d[c][data.variables[c+"_FLAG"][:,0]>1]=fill
            else:
                try:
                    d[c]=data.variables[c][:].filled(fill)
                except AttributeError:
                    d[c]=data.variables[c][:]
                try:
                    d[c][data.variables[c+"_FLAG"][:]>1]=fill
                except KeyError:
                    pass
        data.close()
        self.data=d
        del data
        
    def get_values(self,times,para='ALT_GIN'):
        ind=self.get_indexes(times)
        d=self.data[para][ind]
        d[ind<0]=np.nan
        return d
        
    def get_indexes(self,times):
        return np.digitize(times-1,self.times)
        """
        i=[]
        try:
            for t in times:
                if(t in self.times):
                    i.append(np.where(self.times==t)[0][0])
                else:
                    i.append(-1)
        except TypeError:
            i.append(np.where(self.times==times)[0][0])
        return np.array(i,dtype=int)
        """

    def __getitem__(self,item):
        return self.data[item]

    def write(self):
        np.savetxt(self.filename,self.data,fmt=self.format)
                
    def read(self):
        self.data=np.genfromtxt(self.filename,names=self.columns)

