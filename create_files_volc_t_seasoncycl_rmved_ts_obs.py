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


dicmodels={'mpi_esm_lr': 'MPI-ESM-LR', 'cmcc_cesm': 'CMCC-CESM', 'cesm1_fastchem': 'CESM1-FASTCHEM', 'mri_cgcm3': 'MRI-CGCM3', 'ccsm4': 'CCSM4', 'cmcc_cms': 'CMCC-CMS', 'cesm1_waccm': 'CESM1-WACCM', 'cancm4': 'CanCM4', 'cesm1_cam5': 'CESM1-CAM5', 'cnrm_cm5': 'CNRM-CM5', 'csiro_mk3_6_0': 'CSIRO-Mk3-6-0', 'ipsl_cm5a_lr': 'IPSL-CM5A-LR', 'giss_e2_h_cc': 'GISS-E2-H-CC', 'noresm1-me': 'NorESM1-ME','gfdl_esm2g': 'GFDL-ESM2G', 'cesm1_bgc': 'CESM1-BGC', 'mpi_esm_mr': 'MPI-ESM-MR', 'gfdl_cm3': 'GFDL-CM3', 'noresm1-m': 'NorESM1-M', 'cmcc_cm': 'CMCC-CM', 'hadcm3': 'HadCM3', 'ipsl_cm5b_lr': 'IPSL-CM5B-LR', 'gfdl_esm2m': 'GFDL-ESM2M', 'bcc_csm1_1_m': 'bcc-csm1-1-m', 'inmcm4': 'inmcm4', 'access1_3': 'ACCESS1-3', 'access1_0': 'ACCESS1-0', 'miroc5': 'MIROC5', 'giss_e2_r_cc': 'GISS-E2-R-CC', 'bcc_csm1_1': 'bcc-csm1-1', 'mpi_esm_p': 'MPI-ESM-P', 'hadgem2_cc': 'HadGEM2-CC', 'ipsl_cm5a_mr': 'IPSL-CM5A-MR', 'giss_e2_h': 'GISS-E2-H', 'miroc_esm': 'MIROC-ESM', 'miroc-esm-chem': 'MIROC-ESM-CHEM', 'miroc_esm_lr': 'MIROC-ESM-LR', 'miroc4h': 'MIROC4h', 'fio_esm': 'FIO-ESM', 'canesm2': 'CanESM2', 'giss_e2_r': 'GISS-E2-R', 'bnu_esm': 'BNU-ESM', 'hadgem2_ao': 'HadGEM2-AO', 'hadgem2_es': 'HadGEM2-ES', 'gfdl_cm2_p1': 'GFDL-CM2-P1', 'ec_earth': 'EC-EARTH', 'gfdl_cm2p1': 'GFDL-CM2P1', 'fgoals_s2': 'FGOALS-s2', 'fgoals_g2': 'FGOALS-g2'}

dicmodels2={'bcc_csm1_1':'bcc_r1i1p1',
'canesm2':'canesm2_r1i1p1',
'ccsm4':'ccsm4_r1i1p1',
'cnrm_cm5':'cnrm_r1i1p1',
'csiro_mk3_6_0':'csiro_r1i1p1',
'fgoals-g2_r1i1p1':'fgoals-g2_r1i1p1',
'gfdl_cm3':'gfdl_cm3_r1i1p1',
'gfdl_esm2m':'gfdl_esm2m_r1i1p1',
'giss_e2_h':'giss-e2-h_r1i1p1',
'giss_e2_r':'giss-e2-r_r1i1p1',
'giss_e2_r3':'giss-e2-r_r1i1p3',
'hadgem2_es':'hadgem2-es_r1i1p1',
'ipsl_cm5a_lr':'ipsl-cm5a-lr_r1i1p1',
'miroc-esm-chem':'miroc-esm-chem_r1i1p1',
'miroc_esm':'miroc-esm_r1i1p1',
'mri_cgcm3':'mri-cgcm3_r1i1p1',
'noresm1-m':'noresm1-m_r1i1p1'}


#var3='SPI_12'
#var3='tos'
#var3='ts'
var3='sst'
#var3='tas'
#var3='tls'
#var3='pr'
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
    regionBEN=cdutil.region.domain(latitude=(-82.5,82.5,'ccb')); endfi=''
    regionBEN=cdutil.region.domain(latitude=(-20.,20.,'ccb')); endfi='20S20N'
    #regionBEN=cdutil.region.domain(latitude=(0.,82.5,'ccb')); endfi='NH'
    #regionBEN=cdutil.region.domain(latitude=(-82.5,0.,'ccb')); endfi='SH'
    #regionBEN=cdutil.region.domain(latitude=(30.,45.,'ccb'),longitude=(-122.,-90.,'ccb')); endfi='_NA'

#####################################################
#models=['GISS-E2-R','bcc-csm1-1','MPI-ESM-P','CCSM4','MIROC-ESM','FGOALS-gl','IPSL-CM5A-LR']
#models=['MPI-ESM-P','CCSM4','MIROC-ESM','FGOALS-gl','IPSL-CM5A-LR']
#models=['FGOALS-gl']
#models=['MPI-ESM-P']
#models=['GISS-E2-R','bcc-csm1-1']

models=os.listdir( '/big_disk_1/pcmdi/ipcc/data/hist_rcp85/atm/mo/ts')
rips=['r1i1p','r2i1p','r3i1p','r4i1p','r5i1p','r6i1p','r7i1p','r8i1p','r9i1p','r10i1p','r11i1p','r12i1p']
#rips=['r2i1p','r3i1p']
#rips=['r3i1p','r4i1p','r5i1p','r6i1p','r7i1p','r8i1p']
models.remove('PET0.ESMF_LogFile')
models.remove('ts_hist_rcp85_all_runs_rmclimo0009_removed_lrmIcepoub.nc')

#models=['giss_e2_r','bnu_esm','miroc_esm']

models=MV.arange(1,101,1)
models=list(models)

models=[0]

#oexpt='past1000'   ### because!
oexpt='hist_rcp85'   ### because!
#oexpt='historicalNat'   ### because!

rips=['1']

for m in models:
 for rip in rips:
    varss=[var3]; mmm=m; pp='1'
    if rip=='r2i1p' and m=='giss_e2_h' and oexpt=='hist_rcp85': rip='r1i1p'; pp='2'
    if rip=='r3i1p' and m=='giss_e2_h' and oexpt=='hist_rcp85': rip='r1i1p'; pp='3'
    if rip=='r2i1p' and m=='giss_e2_r' and oexpt=='hist_rcp85': rip='r1i1p'; pp='2'
    if rip=='r3i1p' and m=='giss_e2_r' and oexpt=='hist_rcp85': rip='r1i1p'; pp='3'
    print m,mm

    #if m=='FGOALS-gl': yr1='1000-1-1'
    #if m!='FGOALS-gl': yr1='850-1-1' 
    #if m=='IPSL-CM5A-LR': yr2='1699-12-1'
    #if m!='IPSL-CM5A-LR': yr2='1699-12-31'
    #if m=='IPSL-CM5A-LR': yr2='1805-12-1'
    #if m!='IPSL-CM5A-LR': yr2='1805-12-31'
    special='NO'
    for var in varss:
        print m,var,mm, oexpt
        if var=='ts'  or var=='sst' or var=='tos' or var=='tas' or var=='tls':
            #try: f5 = cdms.open('/work/durack1/Shared/obs_data/temperature/HadCRUT3/santerTimeFormat/HadSST.3.1.0.0.anomalies.'+str(m)+'_santerTimeFormat.nc' )
            #try: f5 = cdms.open('/big_disk_1/pcmdi/ipcc/data/20c3m/atm/mo/ts/BEN/HADISST1.1/20c3m_run1/Xy/ts_picntrl_HADISST1.1_onPCMgrid_rmclimo0009_removed_Ice.nc' )
            #try: f5 = cdms.open('/big_disk_1/pcmdi/ipcc/data/20c3m/atm/mo/ts/BEN/NOAA3b/20c3m_run1/Xy/ts_picntrl_NOAA3b_onPCMgrid_rmclimo0009_removed_Ice.nc' )
            try: f5 = cdms.open('/work/durack1/Shared/obs_data/temperature/ERSST_V3b/130205_ERSST_V3b_185401-201301.nc' )
            except: continue
            #f5 = cdms.open('/work/durack1/Shared/obs_data/temperature/HadCRUT3/santerTimeFormat/HadSST.3.1.0.0.anomalies.'+str(m)+'_santerTimeFormat.nc' )
            #f5 = cdms.open('/big_disk_1/pcmdi/ipcc/data/20c3m/atm/mo/ts/BEN/HADISST1.1/20c3m_run1/Xy/ts_picntrl_HADISST1.1_onPCMgrid_rmclimo0009_removed_Ice.nc' )
            #f5 = cdms.open('/big_disk_1/pcmdi/ipcc/data/20c3m/atm/mo/ts/BEN/NOAA3b/20c3m_run1/Xy/ts_picntrl_NOAA3b_onPCMgrid_rmclimo0009_removed_Ice.nc' )
            f5 = cdms.open('/work/durack1/Shared/obs_data/temperature/ERSST_V3b/130205_ERSST_V3b_185401-201301.nc' ); special='YES'
        if var=='SPI_12':
            if oexpt=='historicalNat': 
                try: f5 = cdms.open('/bonfils2_storage_tmp/cmip5/'+oexpt+'/cmip5'+'.'+mm+'.'+oexpt+'.'+rip+pp+'.mo.atm.Amon.pr.ver-1.SPI_ALL_12.Base_same_1961-1990.nc' )
                except: continue
            if oexpt=='hist_rcp85': 
                try: f5 = cdms.open('/bonfils2_storage_tmp/cmip5/'+oexpt+'/pr'+'_'+mmm+'.'+oexpt+'.'+rip+pp+'.mo.SPI_ALL_12.Base_same_1961-1990.nc' )
                except: continue
            if oexpt=='historicalNat':f5 = cdms.open('/bonfils2_storage_tmp/cmip5/'+oexpt+'/cmip5'+'.'+mm+'.'+oexpt+'.'+rip+pp+'.mo.atm.Amon.pr.ver-1.SPI_ALL_12.Base_same_1961-1990.nc' )
            if oexpt=='hist_rcp85':f5 = cdms.open('/bonfils2_storage_tmp/cmip5/'+oexpt+'/pr'+'_'+mmm+'.'+oexpt+'.'+rip+pp+'.mo.SPI_ALL_12.Base_same_1961-1990.nc' )
        #t11=f5(var3,regionBEN,time=(yr1,yr2))
        if special=='YES': t11=f5(var3,regionBEN,time=('1854','2013','con'))
        else: t11=f5(var3,regionBEN)
        print timereg(t11)
        x.plot(t11)
        #fm=cdms.open(str.split(mod_filesTSH[0][:-1],'ts_')[0]+'sftland_'+str.split(mod_filesTSH[0][:-1],'ts_')[1][:-17]+'nc')
        #tm=fm('ts')
        #t11=f5(var3) ; endfi='_global'
        if string.split(timereg(t11)[0],'-')[1]!='1': raise
        if string.split(timereg(t11)[1],'-')[1]!='12': raise
        #if m=='FGOALS-gl':
        #    if string.split(timereg(t11)[0],'-')[0]!='1000':
        #        raise 
        #else:
        #    if string.split(timereg(t11)[0],'-')[0]!='850':
        #        raise
        #if string.split(timereg(t11)[1],'-')[0]!='1699': raise
        #if string.split(timereg(t11)[1],'-')[0]!='1805': raise
        cdutil.times.setTimeBoundsMonthly(t11)
        print t11.shape
        t1=MV.zeros((t11.shape[0]),typecode=MV.float32)
        t2=MV.zeros((t11.shape[0]),typecode=MV.float32)
        #t1z=MV.zeros((t11.shape[0],t11.shape[1]),typecode=MV.float32)
        #wgt = cdutil.area_weights(t11[0,:,:])
        wgt2 = cdutil.area_weights(t11)
        #for i in range(t1.shape[0]):
        #    #print i
        #    t1[i] = cdutil.averager(t11[i,:,:],axis='xy',weight=wgt)
        #    #t1z[i] = cdutil.averager(t11[i,:,:],axis='x',weight=wgt)
        t1=cdutil.averager(t11,axis='xy',weight=wgt2)
        #raise
        #t1.setAxis(0,t11.getAxis(0))
        t2.setAxis(0,t11.getAxis(0))
        #t1z.setAxis(0,t11.getAxis(0))
        #t1z.setAxis(1,t11.getAxis(1))
            
        ##############
        ##############
        ##############
        ##############
        # comment one or the other    
        ##############
        #MOY GLOBAL
        ##############
        print timereg(t1)
        acts=cdutil.ANNUALCYCLE.climatology(t1(time=('1960','2000','con')))
        print timereg(t1(time=('1960','2000','con')))
        print timereg(t1)
        print acts.shape
        #############acts=cdutil.ANNUALCYCLE.climatology(t1(time=(stclimo[0],enclimo[0])))  # does not work... maybe fixed by Charles
        t2=cdutil.ANNUALCYCLE.departures(t1,ref=acts)# does not work... maybe fixed by Charles
        #for i in range (int(t1.shape[0]/12)):
        #    t11[i*12:i*12+12]=cdutil.ANNUALCYCLE.departures(t1[i*12:i*12+12],ref=acts)
        ###x.plot(t1)
        ##x.plot(t11)
        x.plot(t2)
        if special=='YES': 
            ztax=t2.getTime()
            ztax.toRelativeTime('months since 1800')
        #####
        if var3=='SPI_12': fout=cdms.open('/bonfils2_storage_tmp/cmip5/'+oexpt+'/'+var+'_'+mmm+'_'+oexpt+'.'+rip+pp+'_rmclimo0009_removed_lrmIce'+'.mon.rmSC'+endfi+'.nc','w')
        #else: fout=cdms.open('/bonfils2_storage_tmp/cmip5/obs_data/HadSST.3.1.0.0.anomalies.'+str(m)+'_santerTimeFormat_rmclimo0009_removed_lrmIce'+'.mon.rmSC'+endfi+'.nc','w')
        #else: fout=cdms.open('/bonfils2_storage_tmp/cmip5/obs_data/ts_picntrl_HADISST1.1_onPCMgrid_rmclimo0009_removed_Ice.mon.rmSC'+endfi+'.nc','w')
        #else: fout=cdms.open('/bonfils2_storage_tmp/cmip5/obs_data/ts_picntrl_NOAA3b_onPCMgrid_rmclimo0009_removed_Ice.mon.rmSC'+endfi+'.nc','w')
        else: fout=cdms.open('/bonfils2_storage_tmp/cmip5/obs_data/ts_130205_ERSST_V3b_185401-201301_Ice.mon.rmSC'+endfi+'.nc','w')
        fout.write(t2,id=var3,typecode='f')
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
        ## #fout=cdms.open('/big_disk_1/pcmdi/ipcc/data/hist_rcp85/atm/mo/'+var+'/'+mm+'/'+oexpt+'/'+var+'_both___'+m+'.r2i1p1.'+oexpt+'.mon.rmSC_moyzonal.nc','w')
        ## fout=cdms.open('/bonfils2_storage_tmp/MSU/'+var+'_both___'+m+'.r2i1p1.'+oexpt+'.mon.rmSC_moyzonal_BEN.nc','w')
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

