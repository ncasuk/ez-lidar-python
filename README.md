# EZ-LIDAR-PYTHON

Project to create a community resource using data from the various EZ-LIDARs at the Met Office, and NCAS.

## Getting Started

git clone https://github.com/ncasuk/ez-lidar-python.git

### Prerequisites

Must have Python installed and Numpy at the very least.

### Example use
```
import lidar
```
To create the file - do something like this
```
l_b920=lidar.lidar('2015-08-07_B920.zip',ncfolder='',aux='core_faam_20150807_v004_r3_b920.nc',clobber=True)
```

Merge in the auxilliary data
```
l_b920.merge_aux()
```

Make sure they are closed
```
l_b920.close()
```


To open a lidar file

```
l_b920=lidar.lidar('metoffice-lidar_faam_20150807_r0_B920_raw.nc')
```
### Parameters available
```
Time
Blind_gain0
Blind_offset0
Blind_gain1
Blind_offset1
Blind_NumberOfSignal
Blind_NoiseMean0
Blind_NoiseStd0
Blind_NoiseMean1
Blind_NoiseStd1
Raw_offset0
Raw_gain0
Raw_offset1
Raw_gain1
Raw_NumberOfSignal
Raw_gainPct
Raw_NoiseMean0
Raw_NoiseStd0
Raw_NoiseMean1
Raw_NoiseStd1
AngleAzimuth
AngleZenith
rawSignal_0
rawSignal_1
rawPhoton_0
rawPhoton_1
rawBlind_0
rawBlind_1

Altitude (m)
Longitude (deg)
Latitude (deg)
Pressure (hPa)
Temperature (degC)
Humidity (%)
```
To get some raw data out as an array 
```
Alt=l_b920['Altitude (m)'][:]
```
There are some 'special' parameters which are not directly from the NetCDF
```
l_b920.bind # This is indexing for the Blind measurements which is less frequent than the measurements..


l_b920.profile[chan][n] # profile n channel chan
                       # where chan can be 0, 1 or 2 ( which is the ratio of 1 to 0)
                       # n is a profile number or a slice eg. l_920.profile[0][0] or l_920.profile[0][10:20]
```
Similarly
```
l_b920.range_corrected[chan][n] # for range corrected profiles

l_b920.curtain[chan][n] # for a curtain ( range corrected and altitude taken into account )

l_b920.image[chan][n] # for an image like array ? 

l_b920.trigger # trigger point in profile ( when laser fired - set to 2054 )
```





There are also lots of attributes mostly taken directly from the raw file...

Accessed via
```
l_b920.getncattr('attribute name')

Description
GENERAL INFORMATIONS
Version
HeaderSize
ID ALS
DateRun
Location
User
COMMENTS
Setup
General Comments
Meteo
InternalComment
VARIABLES
ACQUISITION PARAMETERS
Acquisition Mode
Acquisition duration (s)
Profile per sample
Sampling period
Range (m)
RawResolution (m)
OutResolution (m)
DepolarResolution (m)
SCANNING PARAMETERS
AngleZenithInitial (deg)
AngleZenithFinal (deg)
AngleAzimuthInitial (deg)
AngleAzimuthFinal (deg)
AnglesNB AA
AnglesNB ZA
ACCU PARAMETERS
NumberOfShot
LASER PARAMETERS
Wave length (nm)
PRF (Hz)
ALGO PARAMETERS
Ca
DAQ PARAMETERS
voltage0 (V)
voltage1 (V)
SOFT PARAMETERS
IP
Port
LineSeparator
DecimalSeparator
NbOfProfilesPerFile
DataCodage
WritingPosition (byte)
PHOTON COUNTING
OTHER PARAMETERS
Dft Data Directory
Dft ASCII Directory
Azimuth offset with North (deg)
```
## Authors

* **Dave Tiddeman** - *Initial work* - Met Office.

## License

ez-lidar-python may be freely distributed, modified and used commercially under the terms
of its [GNU LGPLv3 license](COPYING.LESSER).

