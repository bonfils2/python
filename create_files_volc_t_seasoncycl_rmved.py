#################################################
###### INTRO ####################################
#################################################
#/usr/local/uvcdat/2012-08-24/bin/python
import os, string, gc, sys,getopt
sys.path.append('/export/bonfils2/PYFORT/bonfils/Functions/')
#sys.path.append('/home/bonfils2/PYFORT/bonfils/Functions/')
import vcs, cdutil, cdms2 as cdms, MV2 as MV, genutil, cdtime
from genutil import statistics
import Lynch1
import pyclimate.svdeofs
from Scientific.IO.NetCDF import *
from pyclimate.svdeofs import *
from pyclimate.ncstruct import *
sys.path.insert(0,"/export/bonfils2/NEWPYFORT/TRANSFO/build/lib.linux-i686-2.5")
#import Lynch
import numpy.oldnumeric as Numeric
import time
from Scientific.IO import FortranFormat
import numpy.oldnumeric.ma as MA
import sys,getopt  # for external loop


value=0
cdms.setNetcdfShuffleFlag(value) ## where value is either 0 or 1
cdms.setNetcdfDeflateFlag(value) ## where value is either 0 or 1
cdms.setNetcdfDeflateLevelFlag(value) ## where value is a integer between 0 and 9 included

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
oregionBen='BEN'
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

oregion='825S825N'
if oregion=='825S825N': #-60.,  -20.,   6., 18.
    regionBEN=cdutil.region.domain(latitude=(-82.5,82.5,'ccb'))

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


#models=['MPI-ESM-P']
#,'CCSM4','CNRM-CM5']
#models=['CCSM4']
#models=['bcc-csm1-1']
#models=['MIROC-ESM']
models=['GISS-E2-R','bcc-csm1-1','MPI-ESM-P','CCSM4','MIROC-ESM','FGOALS-gl','IPSL-CM5A-LR']
#models=['IPSL-CM5A-LR']
#models=['FGOALS-gl']
#models=['CCSM4','MIROC-ESM']
#models=['MIROC-ESM']
models=['MPI-ESM-P']

for m in models:
    if m=='MPI-ESM-P': varss=['tls','tls','tlt','tts']; oexpt='past1000'; mm='mpi-esm-p'; pp='1'
    if m=='CCSM4':     varss=['tmt','tls','tlt','tts']; oexpt='past1000'; mm='ccsm4'; pp='1'
    if m=='CNRM-CM5':  varss=['tmt','tls','tlt','tts']; oexpt='past1000'; mm='cnrm'; pp='1'
    if m=='bcc-csm1-1':varss=['tmt','tls','tlt','tts']; oexpt='past1000'; mm='bcc'; pp='1'
    if m=='MIROC-ESM': varss=['tmt','tls','tlt','tts']; oexpt='past1000'; mm='miroc-esm'; pp='1'
    if m=='GISS-E2-R':    varss=['tmt','tls','tlt','tts']; oexpt='past1000'; mm='giss-e2-r'; pp='121'
    if m=='FGOALS-gl':    varss=['tmt','tls','tlt','tts']; oexpt='past1000'; mm='fgoals-gl'; pp='1'
    if m=='IPSL-CM5A-LR': varss=['tmt','tls','tlt','tts']; oexpt='past1000'; mm='ipsl'; pp='1'

    if m=='FGOALS-gl': yr1='1000-1-1'
    if m!='FGOALS-gl': yr1='850-1-1' 
    if m=='IPSL-CM5A-LR': yr2='1699-12-1'
    if m!='IPSL-CM5A-LR': yr2='1699-12-31'
    #if m=='IPSL-CM5A-LR': yr2='1805-12-1'
    #if m!='IPSL-CM5A-LR': yr2='1805-12-31'

    oexpt='past1000'   ### because!
    for var in varss:
        print m,var,mm, oexpt
        #f5 = cdms.open('/big_disk_1/pcmdi/ipcc/data/hist_rcp85/atm/mo/'+var+'/'+mm+'/'+oexpt+'/'+var+'_both___'+m+'.r1i1p1.'+oexpt+'.mon.nc' )
        f5 = cdms.open('/network_storage/painter1/MSU/'+mm+'_r1i1p'+pp+'/'+var+'_both___'+m+'.r1i1p'+pp+'.'+oexpt+'.mon.nc' )
        #f5 = cdms.open('/big_disk_1/pcmdi/ipcc/data/hist_rcp85/atm/mo/'+var+'/'+mm+'/'+oexpt+'/'+var+'_both___'+m+'.r1i1p1.'+oexpt+'.mon.poub.nc' )
        #t11=f5('eqmsu_'+var)
        t11=f5('eqmsu_'+var,regionBEN,time=(yr1,yr2))
        #t11=f5('eqmsu_'+var,time=(yr1,yr2))
        if string.split(timereg(t11)[0],'-')[1]!='1': raise
        if string.split(timereg(t11)[1],'-')[1]!='12': raise
        if m=='FGOALS-gl':
            if string.split(timereg(t11)[0],'-')[0]!='1000':
                raise 
        else:
            if string.split(timereg(t11)[0],'-')[0]!='850':
                raise
        if string.split(timereg(t11)[1],'-')[0]!='1699': raise
        #if string.split(timereg(t11)[1],'-')[0]!='1805': raise
        cdutil.times.setTimeBoundsMonthly(t11)
        print t11.shape
        t1=MV.zeros((t11.shape[0]),typecode=MV.float32)
        t2=MV.zeros((t11.shape[0]),typecode=MV.float32)
        t1z=MV.zeros((t11.shape[0],t11.shape[1]),typecode=MV.float32)
        for i in range(t1.shape[0]):
            print i
            wgt = cdutil.area_weights(t11[i,:,:])
            t1[i] = cdutil.averager(t11[i,:,:],axis='xy',weight=wgt)
            t1z[i] = cdutil.averager(t11[i,:,:],axis='x',weight=wgt)
        t1.setAxis(0,t11.getAxis(0))
        t2.setAxis(0,t11.getAxis(0))
        t1z.setAxis(0,t11.getAxis(0))
        t1z.setAxis(1,t11.getAxis(1))
            
        ##############
        ##############
        ##############
        ##############
        # comment one or the other    
        ##############
        #MOY GLOBAL
        ##############
        print timereg(t1)
        acts=cdutil.ANNUALCYCLE.climatology(t1)
        print timereg(t1)
        print acts.shape
        #############acts=cdutil.ANNUALCYCLE.climatology(t1(time=(stclimo[0],enclimo[0])))  # does not work... maybe fixed by Charles
        t11=cdutil.ANNUALCYCLE.departures(t1,ref=acts)# does not work... maybe fixed by Charles
        for i in range (int(t1.shape[0]/12)):
            t2[i*12:i*12+12]=cdutil.ANNUALCYCLE.departures(t1[i*12:i*12+12],ref=acts)
        x.plot(t1)
        x.plot(t2)
        x.plot(t11)
        raise
        #####fout=cdms.open('/big_disk_1/pcmdi/ipcc/data/hist_rcp85/atm/mo/'+var+'/'+mm+'/'+oexpt+'/'+var+'_both___'+m+'.r1i1p1.'+oexpt+'.mon.rmSC.nc','w')
        #fout=cdms.open('/bonfils2_storage_tmp/MSU/'+var+'_both___'+m+'.r1i1p1.'+oexpt+'.mon.rmSC.nc','w')
        #fout=cdms.open('/bonfils2_storage_tmp/MSU/'+var+'_both___'+m+'.r1i1p1.'+oexpt+'.mon.rmSC_BEN2.nc','w')
        fout=cdms.open(var+'_both___'+m+'.r1i1p1.'+oexpt+'.mon.rmSC_BEN2_wrong.nc','w')
        fout.write(t2,id='eqmsu_'+var,typecode='f')
        fout.close()

        
        ## ##############
        ## #MOY ZONAL
        ## #############acts=cdutil.ANNUALCYCLE.climatology(t1(time=(stclimo[0],enclimo[0])))
        ## print timereg(t1z), 'timereg 1'
        ## actsz=cdutil.ANNUALCYCLE.climatology(t1z)
        ## print actsz.shape
        ## #
        ## #
        ## for i in range (int(t1z.shape[0]/12)):
        ##     t1z[i*12:i*12+12]=cdutil.ANNUALCYCLE.departures(t1z[i*12:i*12+12],ref=actsz)
        ## #fout=cdms.open('/big_disk_1/pcmdi/ipcc/data/hist_rcp85/atm/mo/'+var+'/'+mm+'/'+oexpt+'/'+var+'_both___'+m+'.r1i1p1.'+oexpt+'.mon.rmSC_moyzonal.nc','w')
        ## fout=cdms.open('/bonfils2_storage_tmp/MSU/'+var+'_both___'+m+'.r1i1p1.'+oexpt+'.mon.rmSC_moyzonal_BEN.nc','w')
        ## fout.write(t1z,id='eqmsu_'+var,typecode='f')
        ## fout.close()
        ## print timereg(t1z), 'timereg 2'
        ## ##############
        ## ##############
        ## ##############
        ## ##############

        ## ##############
        ## # test to see if time step is unique or not
        ## tim=t1.getTime()[:]
        ## import numpy
        ## if len(tim)!=len(numpy.unique(tim)):
        ##     raise
        ## ##############





#tlspi.setAxis(0,tmt[500:500+10200].getTime())

