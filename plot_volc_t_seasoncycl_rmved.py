#################################################
###### INTRO ####################################
#################################################
import os, string, gc, sys,getopt
#sys.path.append('/export/bonfils2/PYFORT/bonfils/Functions/')
#sys.path.append('/home/bonfils2/PYFORT/bonfils/Functions/')
import vcs, cdutil, cdms2 as cdms, MV2 as MV, genutil, cdtime
from genutil import statistics
#import Lynch1
import pyclimate.svdeofs
from Scientific.IO.NetCDF import *
from pyclimate.svdeofs import *
from pyclimate.ncstruct import *
#sys.path.insert(0,"/home/bonfils2/NEWPYFORT/TRANSFO/build/lib.linux-i686-2.5")
#import Lynch
import numpy.oldnumeric as Numeric
import time
from Scientific.IO import FortranFormat
import numpy.oldnumeric.ma as MA
import sys,getopt  # for external loop

f15=open('trend_trenderr_values_AUTRESPDOP_revised.txt','a')

###################################################
args=sys.argv[1:]
letters='e:i:f:r'
keywords=['exper=','ice=','filt=','rang=']
oexpt='default'
oice='default'
lowpass='default'
Pdateclimo='default'
opts,pargs=getopt.getopt(args,letters,keywords)
for o,p in opts:
    if o in ['-e','--exper']:
        oexpt=p
    if o in ['-i','--ice']:
        oice=p
    if o in ['-f','--filt']:
        lowpass=p
    if o in ['-r','--rang']:
        Pdateclimo=p
## ###################################


def timereg(s1):
    y1time=str(cdtime.reltime(s1.getTime()[0], s1.getTime().units).tocomp(s1.getTime().getCalendar()))
    y2time=str(cdtime.reltime(s1.getTime()[-1], s1.getTime().units).tocomp(s1.getTime().getCalendar()))
    return y1time, y2time


x=vcs.init()
tpl1 = x.gettemplate( 'ASD1_of_4')
tpl2 = x.gettemplate( 'ASD2_of_4')
tpl3 = x.gettemplate( 'ASD3_of_4')
tpl4 = x.gettemplate( 'ASD4_of_4')
tpl1.mean.priority=0; tpl2.mean.priority=0; tpl3.mean.priority=0; tpl4.mean.priority=0
tpl1.max.priority=0; tpl2.max.priority=0; tpl3.max.priority=0; tpl4.max.priority=0
tpl1.min.priority=0; tpl2.min.priority=0; tpl3.min.priority=0; tpl4.min.priority=0
tpl1.comment1.x = 0.159090995789
tpl1.comment1.y = 0.979090995789
tpl1.comment2.x = 0.159090995789
tpl1.comment2.y = 0.949090995789
tpl2.comment1.x = 0.659090995789
tpl2.comment1.y = 0.979090995789
tpl2.comment2.x = 0.659090995789
tpl2.comment2.y = 0.949090995789
tpl3.comment1.x = 0.109090995789
tpl3.comment1.y = 0.209090995789
tpl4.comment1.x = 0.659090995789
tpl4.comment1.y = 0.209090995789
tpl1.line1.priority=1
tpl1.line1.x1=tpl1.box1.x1
tpl1.line1.x2=tpl1.box1.x2
tpl1.line1.y1=(tpl1.box1.y1+tpl1.box1.y2)/2.
tpl1.line1.y2=(tpl1.box1.y1+tpl1.box1.y2)/2.
tpl2.line1.priority=1
tpl2.line1.x1=tpl2.box1.x1
tpl2.line1.x2=tpl2.box1.x2
tpl2.line1.y1=(tpl2.box1.y1+tpl2.box1.y2)/2.
tpl2.line1.y2=(tpl2.box1.y1+tpl2.box1.y2)/2.
tpl3.line1.priority=1
tpl3.line1.x1=tpl3.box1.x1
tpl3.line1.x2=tpl3.box1.x2
tpl3.line1.y1=(tpl3.box1.y1+tpl3.box1.y2)/2.
tpl3.line1.y2=(tpl3.box1.y1+tpl3.box1.y2)/2.
tpl4.line1.priority=1
tpl4.line1.x1=tpl4.box1.x1
tpl4.line1.x2=tpl4.box1.x2
tpl4.line1.y1=(tpl4.box1.y1+tpl4.box1.y2)/2.
tpl4.line1.y2=(tpl4.box1.y1+tpl4.box1.y2)/2.
tpl2.source.priority = 0
tyw1=x.createyxvsx('new01agg')
tyw2=x.createyxvsx('new01aggg')
tyw3=x.createyxvsx('new01agggg')
tyw1.datawc_y1=-2.5 ; tyw1.datawc_y2=2.5
tyw2.datawc_y1=-2.5 ; tyw2.datawc_y2=2.5
tyw3.datawc_y1=-2.5 ; tyw3.datawc_y2=2.5
tyw1.linecolor=241 ; tyw1.line='solid' ; tyw1.linewidth=2 # black
tyw2.linecolor=242 ; tyw2.line='solid' ; tyw2.linewidth=2 # red
tyw3.linecolor=53  ; tyw3.line='solid' ; tyw3.linewidth=2 # red


#OBSPROJ='HADISST1'
#OBSPROJ='NOAA2'
#oice='_lrmIce'
#oice='_rmdIce'
#oice='_Ice'
#ONOAA3='NOAA3'
#OHADISST1='HADISST1'
ONOAA3='NOAA3b'
OHADISST1='HADISST1.1'

OBSPROJ3='ncdc'

if lowpass=='_lynch119': lowpass2='_lynch119'
if lowpass=='_nofiltr': lowpass2=''

#oregion1='PDO' #'PCR'
oregion2='IPO'
oregion1='60S60N'
oregion3='12N60N'
OBSPROJ=oexpt
gg='pc1Arec'

ICEPROJ='Ice'
oice2='_lrmIce'
Palors='8models_weight'

if Pdateclimo=='19001993': peri=''
if Pdateclimo=='18571906': peri='_18571906'
if Pdateclimo=='19071956': peri='_19071956'
if Pdateclimo=='19572006': peri='_19572006'

oregion='60S60N'
if oregion=='60S60N': #-60.,  -20.,   6., 18.
    regionIPDO=cdutil.region.domain(latitude=(-62.,62.,'ccb'))

#####################################################
## f5 = cdms.open('/export/bonfils2/DATA/CA/SAOD/IVI2LoadingLatHeight501-2000_Oct2012.nc' )
## Loading=f5('Loading')
## f5.close()
## cdutil.times.setTimeBoundsMonthly(Loading)
## #Loading2 = cdutil.averager(Loading,axis='yz',weight='equal')
## Loading = cdutil.averager(Loading,axis='y')
## cumul=0
## count=0
## for aa in range(43):
##     cumul=Loading[:,aa]+cumul
##     count=count+1
## print cumul.shape
## print count
## Loading=cumul/count
## Loading=Loading[4188:16188]
## #x.plot(Loading,tpl1)
## #raise



models=['MPI-ESM-P','CNRM-CM5']

for m in models:
    if m=='MPI-ESM-P': varss=['tmt','tls','tlt','tts']; oexpt='past1000'; mm='mpi_esm_p'
    if m=='CNRM-CM5': varss=['tls']; oexpt='piControl'
    for var in varss:
        #f5 = cdms.open('/big_disk_1/pcmdi/ipcc/data/hist_rcp85/atm/mo/'+var+'/'+mm+'/'+oexpt+'/'+var+'_both___'+m+'.r1i1p1.'+oexpt+'.mon.nc' )
        #f5 = cdms.open('/big_disk_1/pcmdi/ipcc/data/hist_rcp85/atm/mo/'+var+'/'+mm+'/'+oexpt+'/'+var+'_both___'+m+'.r1i1p1.'+oexpt+'.mon.poub.nc' )
        f5 = cdms.open('./'+var+'_both___'+m+'.r1i1p1.'+oexpt+'.mon.poub.nc' )
        t1=f5('eqmsu_'+var)
        #m0=fs('eqmsu_'+var,regionIPDO)
        cdutil.times.setTimeBoundsMonthly(t1)
        #wgt = cdutil.area_weights(t1[:,:,:])
        #print t1.shape
        #t1 = cdutil.averager(t1,axis='xy',weight=wgt)
        #acts=cdutil.ANNUALCYCLE.climatology(t1(time=(stclimo[0],enclimo[0])))
        acts=cdutil.ANNUALCYCLE.climatology(t1)
        print acts.shape
        t11=cdutil.ANNUALCYCLE.departures(t1,ref=acts)
        raise
        for i in range (int(t1.shape[0]/12)):
            t1[i*12:i*12+12]=cdutil.ANNUALCYCLE.departures(t1[i*12:i*12+12],ref=acts)
        fout=cdms.open('/big_disk_1/pcmdi/ipcc/data/hist_rcp85/atm/mo/'+var+'/'+mm+'/'+oexpt+'/'+var+'_both___'+m+'.r1i1p1.'+oexpt+'.mon.rmSC.nc','w')
        fout.write(t1,id='eqmsu_'+var,typecode='f')
        fout.close()
raise


#tlspi.setAxis(0,tmt[500:500+10200].getTime())

