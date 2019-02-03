# naming - set for output file naming - see at end of program
fn               = 'ctepyfs-test-'
# input files with settings and initial values
inputcsv         = 'CTEPYinputdf.csv'
b1fn             = 'b1bunches.csv'
b2fn             = 'b2bunches.csv'

import numpy as np
import pandas as pd
import datetime
import time
from collections import OrderedDict
import functools
from scipy import constants as const
from random import randint
import ast

# importing FORTRAN modules
import CTEsetparameters as setparam
import getemittance as fgetemit
import CTEwrite as cwrite
import rfupdate as frfupdate
import sptranschrom as fsptranschrom
import keeptransvprof as fkeeptransprof
import CTEraddamping as setrad
import raddamp as frad
import ibslong as fibs
import collimation as fcollim
import blowup as fblowup
import collision as fcollision
import addparticles as faddpart
import ctepyfunctions as ctef


# dictionary keys that are used to read from the input files and construct the output files
switcheskeys   = ['RFswitch','betatronSwitch','raddampSwitch','IBSswitch','collimationSwitch','blowupSwitch','collisionSwitch','writeAllCoordSwitch', 'writeMountSwitch']
ringparamskeys = ['aatomb2', 'qatomb2', 'aatomb2', 'qatomb2', 'proton', 'ion', 'thib', 'nturns', 'nMacro', 'timeRatio', 'nwrite', 'gammat', 'circ', 'gamma1', 
                  'gamma2', 'vrfb2', 'nharmb2', 'vrf2b2', 'nharm2b2', 'vrfb2', 'nharmb2', 'vrf2b2', 'nharm2b2', 'tunex1', 'tuney1', 'chromx1', 'chromy1', 
                  'tunex2', 'tuney2', 'chromx2', 'chromy2', 'dqmin1', 'dqmin2', 'k2L1', 'k2Lskew2', 'k2L2', 'k2Lskew1', 'tfsb2', 'tfsb2']
startingconditionskeys = ['ex1', 'ey1', 'npart1', 'blen1', 'rmsDelta1', 'tauhat1', 'ex2', 'ey2', 'npart2', 'blen2', 'rmsDelta2', 'tauhat2', 
                          'longcoor', 'transcoor', 'bunchLenPrecis', 'power', 'alint', 'coordfn1', 'coordfn2', 'blen3', 'ex3', 'ey3', 'blen4', 'ex4', 'ey4', 
                          'npart3', 'npart4']
radiationdampingkeys   = ['radMethod', 'tradlong', 'tradperp', 'siglong', 'sigperp', 'rho0']
ibsdictkeys            = ['ibsMethod', 'coupleIBS', 'coulombLog', 'fracibstot', 'nbins', 'piwinski', 'bane', 'intfile']
colldictkeys           = ['refEmxy', 'cutoffAmpl', 'collimAvgSwitch', 'emitMethod', 'nSigCutBeta', 'nSigCutMom', 'betaxMom', 'dispxMom']
blowupdictkeys         = ['pxKickFac', 'pyKickFac', 'blowupMethod']
collisiondictkeys      = ['collRoutine', 'sigI', 'longIntBins', 'ips_leveling', 'ips_leveling_values']

dfintcols    = ['sim.turn','t(hours)','N1_macro','N1_real','NlostLum1','SumL1','NlostDebunch1','SumD1','NLostBet1','Sumb2','NlostMom1','SumM1']
dfemitcols   = ['sim.turn','t(hours)', 'ex1(m)','ey1(m)','el1(eV/s/charge)','sig1_T1','sig1_dP/P_2']
dfibscols    = ['sim.turn','t(hours)','Tp(hours)','Tx(hours)','Ty(hours)']
dfcollimcols = ['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)']

# -----------------------------------------------------------------------------------
# starting from the main program
# -----------------------------------------------------------------------------------

st = time.clock()
print 'start at:',time.ctime()

# constants
protonmass      = const.physical_constants['proton mass energy equivalent in MeV'][0]/1000 # GeV
electroncharge  = const.physical_constants['elementary charge'][0]

# reading in the data from the input files and constructing dictionaries for the input data
switches                                  = ctef.getswitches(inputcsv)
ringparams,dftfsb2head,dftfsb2head,energy = ctef.getringparams(inputcsv)
radiationdamping                          = ctef.getradiationdamping(inputcsv)
ibsdict                                   = ctef.getibsdict(inputcsv)
colldict                                  = ctef.getcolldict(inputcsv)
blowupdict                                = ctef.getblowupdict(inputcsv)
collisiondict                             = ctef.getcollisiondict(inputcsv,ringparams)
startingconditions,df1,df2                = ctef.getstartingconditions(inputcsv,b1fn,b2fn)

# converting string representation of dict to dict
collisiondict['ips_leveling']= ast.literal_eval(collisiondict['ips_leveling'])
collisiondict['ips_leveling_values']= ast.literal_eval(collisiondict['ips_leveling_values'])

# transforming normalized emittances to geometric ones
df1['ex'] = df1['exn'] / ringparams['gamma1']
df1['ey'] = df1['eyn'] / ringparams['gamma1']
df2['ex'] = df2['exn'] / ringparams['gamma2']
df2['ey'] = df2['eyn'] / ringparams['gamma2']

# converting the harmonic number to floats - to use in fortran functions
fnharmb1  = float(ringparams['nharmb1'])
fnharm2b1 = float(ringparams['nharm2b1'])
fnharmb2  = float(ringparams['nharmb2'])
fnharm2b2 = float(ringparams['nharm2b2'])

# initialisation of variables
hglassfacOld,betaSMin,betar1,betar2,vrev1,vrev2,trev1,trev2,trf1,trf2,frf1,frf2,omegarf1,omegarf2,eta1,eta2,eqtime1,eqtime2=setparam.set_parameters1(ringparams,collisiondict)


# slip factors
eta1 = 1/ringparams['gammat'] **2-1/ringparams['gamma1']**2
eta2 = 1/ringparams['gammat'] **2-1/ringparams['gamma2']**2

# calculating momentum spread and adding to the bunch data
df1 = df1.assign(rmsDelta = lambda row: ctef.getsigpfromvrfandsigs(row['blen'],trev1,ringparams['nharm2b1'],-1*ringparams['vrf2b1'],eta1,betar1,energy))
df2 = df2.assign(rmsDelta = lambda row: ctef.getsigpfromvrfandsigs(row['blen'],trev2,ringparams['nharm2b2'],-1*ringparams['vrf2b2'],eta2,betar2,energy))

# create lists of bucket number for both beams, these values are in the beam input files (csv files) 
# This construction is chosen to be compatible with data retrieved from timber. In a later stage
# initial setup can be extracted from timber data and directly compared with simulations (automation)
b1 = list(df1['bucket'])
b2 = list(df2['bucket'])

# initialization of a set of parameters
betax1,betay1,betax2,betay2,pcoeff1,tcoeff1,rcoeff1,pcoeff2,tcoeff2,rcoeff2,dgamma_hat1,dgamma_hat2,tunes1,tunes2,tunesexact1,tunesexact2,dftfs1,dftfs2=setparam.set_parameters2(ringparams,startingconditions,omegarf1,trev1,eta1,betar1,omegarf2,trev2,eta2,betar2)

# read bane and interpolate datafiles or create files containing zeros
dfbane,dfinterpolate = ctef.readibsinputfiles(switches,ibsdict)

# raddamping initialization
if (switches['raddampSwitch']==True)&(radiationdamping['radMethod']!='manual'):
    #       calcualte radiation damping times and equilibrium emittances if not given manually
    #       Reference: Chao, Tigner: Handbook of Accelerator physics and engineering, (1998) page 186
    if (radiationdamping['radMethod']=="approx"):
        I21,I31,I4x1,I4y1,I5x1,I5y1 = setrad.raddampingapprox(ringparams,radiationdamping,betax1,betay1)
        I22,I32,I4x2,I4y2,I5x2,I5y2 = setrad.raddampingapprox(ringparams,radiationdamping,betax2,betay2)
    elif (radiationdamping['radMethod']=='lattic'):
    #   calculate radiation integrals over lattice using twiss file
        I21,I31,I4x1,I4y1,I5x1,I5y1 = setrad.raddampinglattic(ringparams,radiationdamping,dftfs1)
        I22,I32,I4x2,I4y2,I5x2,I5y2 = setrad.raddampinglattic(ringparams,radiationdamping,dftfs2)
        
	print 'I21={:.6e}, I31={:.6e}, I4x1={:.6e}, I4y1={:.6e}, I5x1={:.6e}, I5y1={:.6e}'.format(I21,I31,I4x1,I4y1,I5x1,I5y1)
	print 'I22={:.6e}, I32={:.6e}, I4x2={:.6e}, I4y2={:.6e}, I5x2={:.6e}, I5y2={:.6e}'.format(I22,I32,I4x2,I4y2,I5x2,I5y2)
    
    # initialisation of some more parameters
    rIon1      = (ringparams['qatomb1']**2/ringparams['aatomb1'])*1.54e-18
    rIon2      = (ringparams['qatomb2']**2/ringparams['aatomb2'])*1.54e-18
    CalphaE3C1 = rIon1*ringparams['gamma1']**3/(3*trev1)  # C_alpha*E^3/C in handbook formulas p 221 eq 11
    CalphaE3C2 = rIon2*ringparams['gamma2']**3/(3*trev2)

    #print CalphaE3C1,CalphaE3C2
    tradperp1,tradlong1,Cq1,sigE0E01,siglong1,Jx1,Jy1,sigperp1 = setrad.raddampingconstfromI(I21,I31,I4x1,I4y1,I5x1,I5y1,1,CalphaE3C1,ringparams)
    tradperp2,tradlong2,Cq2,sigE0E02,siglong2,Jx2,Jy2,sigperp2 = setrad.raddampingconstfromI(I22,I32,I4x2,I4y2,I5x2,I5y2,2,CalphaE3C2,ringparams)
 
    print 'calculating radiation damping times and equilibrium beam sizes:'
    print 'Tx=',tradperp1/3600,'h,  Ts=',tradlong1/3600,'h,  sigPt=',siglong1,'sigperp=',sigperp1,'m'
    print 'Tx=',tradperp2/3600,'h,  Ts=',tradlong2/3600,'h,  sigPt=',siglong2,'sigperp=',sigperp2,'m'

# setting up collision stuff - which bunchec collide where and so on
tdc,lumip1,lumip2,lumip5,lumip8 = ctef.genlumikeys(b1,b2)
enc1,enc2,colls                 = ctef.getIPenc(np.array(b1),np.array(b2))

hglassfacOld = OrderedDict()
hglassfacOld['ip1'] = np.ones(colls[0])
hglassfacOld['ip2'] = np.ones(colls[1])
hglassfacOld['ip5'] = np.ones(colls[2])
hglassfacOld['ip8'] = np.ones(colls[3])

# initialisation of the variables to track evolution
collisiontimingdict,lumip1dict,lumip2dict,lumip5dict,lumip8dict = ctef.genlumikeys(b1,b2)

b1intdict  = OrderedDict()
b2intdict  = OrderedDict()
b1emitdict = OrderedDict()
b2emitdict = OrderedDict()
b1ibsdict  = OrderedDict()
b2ibsdict  = OrderedDict()
b1dfcolimpxbeta = OrderedDict()
b1dfcolimpybeta = OrderedDict()
b1dfcolimpxmom  = OrderedDict()
b1dfcolimpymom  = OrderedDict()
b2dfcolimpxbeta = OrderedDict()
b2dfcolimpybeta = OrderedDict()
b2dfcolimpxmom  = OrderedDict()
b2dfcolimpymom  = OrderedDict()

for b in b1:
    b1intdict[b]  = pd.DataFrame(columns=dfintcols)
    b1emitdict[b] = pd.DataFrame(columns=dfemitcols)
    b1ibsdict[b]  = pd.DataFrame(columns=dfibscols)
    b1dfcolimpxbeta[b] = pd.DataFrame(columns=dfcollimcols)
    b1dfcolimpybeta[b] = pd.DataFrame(columns=dfcollimcols)
    b1dfcolimpxmom[b]  = pd.DataFrame(columns=dfcollimcols)
    b1dfcolimpymom[b]  = pd.DataFrame(columns=dfcollimcols)
for b in b2:
    b2intdict[b]  = pd.DataFrame(columns=dfintcols)
    b2emitdict[b] = pd.DataFrame(columns=dfemitcols)
    b2ibsdict[b]  = pd.DataFrame(columns=dfibscols)
    b2dfcolimpxbeta[b] = pd.DataFrame(columns=dfcollimcols)
    b2dfcolimpybeta[b] = pd.DataFrame(columns=dfcollimcols)
    b2dfcolimpxmom[b]  = pd.DataFrame(columns=dfcollimcols)
    b2dfcolimpymom[b]  = pd.DataFrame(columns=dfcollimcols)

v001,v002 = ctef.setv00(ringparams)

# generating the distributions
b1distrdict = OrderedDict()
b2distrdict = OrderedDict()

b1avgline   = OrderedDict()
b1avglinex  = OrderedDict()
b1avgliney  = OrderedDict()
b2avgline   = OrderedDict()
b2avglinex  = OrderedDict()
b2avgliney  = OrderedDict()

b1nLostLum        = OrderedDict()
b1nLostLumSum     = OrderedDict()
b1nLostDebunch    = OrderedDict()
b1nLostDebunchSum = OrderedDict()
b1nLostMom        = OrderedDict()
b1nLostMomSum     = OrderedDict()
b1nLostBeta       = OrderedDict()
b1nLostBetaSum    = OrderedDict()

b2nLostLum        = OrderedDict()
b2nLostLumSum     = OrderedDict()
b2nLostDebunch    = OrderedDict()
b2nLostDebunchSum = OrderedDict()
b2nLostMom        = OrderedDict()
b2nLostMomSum     = OrderedDict()
b2nLostBeta       = OrderedDict()
b2nLostBetaSum    = OrderedDict()

b1nLostLum_ip1    = OrderedDict()
b1nLostLum_ip2    = OrderedDict()
b1nLostLum_ip5    = OrderedDict()
b1nLostLum_ip8    = OrderedDict()
b2nLostLum_ip1    = OrderedDict()
b2nLostLum_ip2    = OrderedDict()
b2nLostLum_ip5    = OrderedDict()
b2nLostLum_ip8    = OrderedDict()

b1dfparticles     = OrderedDict()
b2dfparticles     = OrderedDict()

b1dfwritemount    = OrderedDict()
b2dfwritemount    = OrderedDict()
b1dfwritemounttrans= OrderedDict()
b2dfwritemounttrans= OrderedDict()

# storing the initial number of real particles
b1pnumber = OrderedDict()
b2pnumber = OrderedDict()

# for output ibs
b1ibsrow = OrderedDict()
b2ibsrow = OrderedDict()

for b in b1:
    t,pt,x,px,y,py = ctef.addparticles(df1[df1['bucket']==b],startingconditions['tauhat1'],tcoeff1,trev1,v001,
                 fnharmb1,fnharm2b1,ringparams['vrfb1'],ringparams['vrf2b1'],ringparams['gamma1'],
                 startingconditions['bunchLenPrecis'],ringparams['nMacro'],3,colldict['cutoffAmpl'],betax1,betay1)
    b1distrdict[b]       = OrderedDict() 
    b1distrdict[b]['t']  = t
    b1distrdict[b]['pt'] = pt
    b1distrdict[b]['x']  = x
    b1distrdict[b]['px'] = px
    b1distrdict[b]['y']  = y
    b1distrdict[b]['py'] = py
    b1avgline[b]         = np.zeros(int(ibsdict['nbins']))
    b1avglinex[b]        = np.zeros(int(ibsdict['nbins']))
    b1avgliney[b]        = np.zeros(int(ibsdict['nbins']))
    b1nLostLum[b]        = 0
    b1nLostLumSum[b]     = 0
    b1nLostDebunch[b]    = 0
    b1nLostDebunchSum[b] = 0
    b1nLostMom[b]        = 0
    b1nLostMomSum[b]     = 0
    b1nLostBeta[b]       = 0
    b1nLostBetaSum[b]    = 0
    b1nLostLum_ip1[b]    = 0
    b1nLostLum_ip2[b]    = 0
    b1nLostLum_ip5[b]    = 0
    b1nLostLum_ip8[b]    = 0
    b1dfparticles[b]     = pd.DataFrame(columns=['x','px','y','py','t','pt'])
    b1dfparticles[b]['x']  = x
    b1dfparticles[b]['px'] = px
    b1dfparticles[b]['y']  = y
    b1dfparticles[b]['py'] = py
    b1dfparticles[b]['t']  = t
    b1dfparticles[b]['pt'] = pt
    b1pnumber[b]           = df1[df1['bucket']==b]['npart']
    
for b in b2:
    t,pt,x,px,y,py = ctef.addparticles(df2[df2['bucket']==b],startingconditions['tauhat2'],tcoeff2,trev2,v002,
                 fnharmb2,fnharm2b2,ringparams['vrfb2'],ringparams['vrf2b2'],ringparams['gamma2'],
                 startingconditions['bunchLenPrecis'],ringparams['nMacro'],3,colldict['cutoffAmpl'],betax2,betay2)
    b2distrdict[b]       = OrderedDict() 
    b2distrdict[b]['t']  = t
    b2distrdict[b]['pt'] = pt
    b2distrdict[b]['x']  = x
    b2distrdict[b]['px'] = px
    b2distrdict[b]['y']  = y
    b2distrdict[b]['py'] = py
    b2avgline[b]         = np.zeros(int(ibsdict['nbins']))
    b2avglinex[b]        = np.zeros(int(ibsdict['nbins']))
    b2avgliney[b]        = np.zeros(int(ibsdict['nbins']))
    b2nLostLum[b]        = 0
    b2nLostLumSum[b]     = 0
    b2nLostDebunch[b]    = 0
    b2nLostDebunchSum[b] = 0
    b2nLostMom[b]        = 0
    b2nLostMomSum[b]     = 0
    b2nLostBeta[b]       = 0
    b2nLostBetaSum[b]    = 0
    b2nLostLum_ip1[b]    = 0
    b2nLostLum_ip2[b]    = 0
    b2nLostLum_ip5[b]    = 0
    b2nLostLum_ip8[b]    = 0
    b2dfparticles[b]     = pd.DataFrame(columns=['x','px','y','py','t','pt'])
    b2dfparticles[b]['x']  = x
    b2dfparticles[b]['px'] = px
    b2dfparticles[b]['y']  = y
    b2dfparticles[b]['py'] = py
    b2dfparticles[b]['t']  = t
    b2dfparticles[b]['pt'] = pt
    b2pnumber[b]           = df2[df2['bucket']==b]['npart']

# setting lumi levelling parameters
if collisiondict['ips_leveling']['IP1']==True:
    lumimax1 = collisiondict['ips_leveling_values']['IP1']
else: 
    lumimax1 = 0.0
if collisiondict['ips_leveling']['IP2']==True:
    lumimax2 = collisiondict['ips_leveling_values']['IP2']
else: 
    lumimax2 = 0.0
if collisiondict['ips_leveling']['IP5']==True:
    lumimax5 = collisiondict['ips_leveling_values']['IP5']
else: 
    lumimax5 = 0.0
if collisiondict['ips_leveling']['IP8']==True:
    lumimax8 = collisiondict['ips_leveling_values']['IP8']
else: 
    lumimax8 = 0.0
    
# crossing angles
theta1=abs(dftfs1[dftfs1['NAME']=='IP1']['PY'].values[0])
theta2=abs(dftfs1[dftfs1['NAME']=='IP2']['PY'].values[0])
theta5=abs(dftfs1[dftfs1['NAME']=='IP5']['PX'].values[0])
theta8=abs(dftfs1[dftfs1['NAME']=='IP8']['PX'].values[0])

# defining collimation cuts 
xcut1 = colldict['cutoffAmpl'] * np.sqrt(betax1*colldict['refEmxy'])
ycut1 = colldict['cutoffAmpl'] * np.sqrt(betay1*colldict['refEmxy'])
xcut2 = colldict['cutoffAmpl'] * np.sqrt(betax2*colldict['refEmxy'])
ycut2 = colldict['cutoffAmpl'] * np.sqrt(betay2*colldict['refEmxy'])

# beta stars
bs1=dftfs1[dftfs1['NAME']=='IP1']['BETX'].values[0]
bs2=dftfs1[dftfs1['NAME']=='IP2']['BETX'].values[0]
bs5=dftfs1[dftfs1['NAME']=='IP5']['BETX'].values[0]
bs8=dftfs1[dftfs1['NAME']=='IP8']['BETX'].values[0]

# bane ibs and interpolat
if ibsdict['ibsMethod'] != 'baneApprox':
    gtab = np.array([])
    xmin,xmax,ymin,ymax,zmin,zmax = (0,0,0,0,0,0)
    pnuminterp = 0
else:
    print 'under construction'
    quit()

if ibsdict['ibsMethod'] != 'interpolat':
    ap,ax,ay = (np.array([]),np.array([]),np.array([]))
else:
    print 'under construction'
    quit()
    
print 'current time:',time.ctime()
print'Starting initialization of emittances ...'

b1exdict = OrderedDict()
b1eydict = OrderedDict()
b2exdict = OrderedDict()
b2eydict = OrderedDict()

# initialization of emittances for all the bunches
for b in b1:
    b1exdict[b], b1eydict[b]= fgetemit.getemittance(b1dfparticles[b]['x'].values,b1dfparticles[b]['y'].values,b1dfparticles[b]['px'].values,b1dfparticles[b]['py'].values,
                                                   betax1,betay1,xcut1,ycut1,colldict['cutoffAmpl'],colldict['emitMethod'],colldict['refEmxy'])
    print 'bunch= {:4d}, ex={:.6e}, ey={:.6e}'.format(b,b1exdict[b],b1eydict[b])
for b in b2:
    b2exdict[b], b2eydict[b]= fgetemit.getemittance(b2dfparticles[b]['x'].values,b2dfparticles[b]['y'].values,b2dfparticles[b]['px'].values,b2dfparticles[b]['py'].values,
                                                   betax2,betay2,xcut2,ycut2,colldict['cutoffAmpl'],colldict['emitMethod'],colldict['refEmxy'])
    print 'bunch= {:4d}, ex={:.6e}, ey={:.6e}'.format(b,b2exdict[b],b2eydict[b])

#  write first turn
kturns  = 0
iwrite  = 1
eqtime  = ringparams['nturns'] * trev1 * ringparams['timeRatio']/3600.

# writing the first turn data to output files
print 'current time:',time.ctime()
print'Starting writing of first turn data ...'

# write avgline distributions
if (switches['writeMountSwitch']):
    ctef.writedistr(b1avgline,b1avglinex,b1avgliney,b1exdict,b1eydict,ringparams['nMacro'],startingconditions['npart1'],ringparams['qatomb1'],
                   ringparams['nwrite'],kturns,ringparams['thib'],ibsdict['nbins'],betax1,betay1)
    ctef.writedistr(b2avgline,b2avglinex,b2avgliney,b2exdict,b2eydict,ringparams['nMacro'],startingconditions['npart2'],ringparams['qatomb2'],
                   ringparams['nwrite'],kturns,ringparams['thib'],ibsdict['nbins'],betax2,betay2)
    
# write emittances
for b in b1:
    b1emitdict[b] = cwrite.writeEmi(b1exdict[b],b1eydict[b],b1dfparticles[b]['t'].values,b1dfparticles[b]['pt'].values,betar1,ringparams['gamma1'],
                              ringparams['aatomb1'],ringparams['qatomb1'],eqtime,kturns,ringparams['nturns'],b1emitdict[b])  
for b in b2:
    b2emitdict[b] = cwrite.writeEmi(b2exdict[b],b2eydict[b],b2dfparticles[b]['t'].values,b2dfparticles[b]['pt'].values,betar2,ringparams['gamma2'],
                              ringparams['aatomb2'],ringparams['qatomb2'],eqtime,kturns,ringparams['nturns'],b2emitdict[b])
    
#     write all starting conditions to file
if (switches['writeAllCoordSwitch']==1): 
    for b in b1:
        b1dfparticles[b].to_csv('all-coord-b1-'+str(b)+ '_' + str(kturn)+'.csv',index=False)
    for b in b2:
        b2dfparticles[b].to_csv('all-coord-b2-'+str(b)+ '_' + str(kturn)+'.csv',index=False)

# calculation of parameters
psix1,psiy1,coeffchromx1,coeffchromy1,tpdqmin1 = setparam.settunechroma(1,ringparams)
psix2,psiy2,coeffchromx2,coeffchromy2,tpdqmin2 = setparam.settunechroma(2,ringparams)

# extra variables for convenience - reducing of the typing
nwrite = ringparams['nwrite']
nturns = ringparams['nturns']
fmix1  = 1

# start of the loop 
print 'start at:',time.ctime(),' took: ',time.clock()-st
print 'Starting main loop...'

b1np = OrderedDict()
b2np = OrderedDict()

df = b2dfparticles.copy()
for kturns in range(1,int(ringparams['nturns'])+1):
    kcheck = kturns % int(nwrite)
    #kcheck = kcheck*nwrite
    
    if(kcheck==0):
        iwrite = 1
        print 'took :', time.clock()-st
    else:
        iwrite = 0
    
    if (kturns==nturns+1):
        print ' writing last turn'
        nwrite = -1
        iwrite =  1

# calculating new emittances
    for b in b1:
        b1exdict[b], b1eydict[b]= fgetemit.getemittance(b1dfparticles[b]['x'].values,b1dfparticles[b]['y'].values,
                                                        b1dfparticles[b]['px'].values,b1dfparticles[b]['py'].values,
                                                        betax1,betay1,xcut1,ycut1,colldict['cutoffAmpl'],colldict['emitMethod'],colldict['refEmxy'])
    for b in b2:
        b2exdict[b], b2eydict[b]= fgetemit.getemittance(b2dfparticles[b]['x'].values,b2dfparticles[b]['y'].values,
                                                        b2dfparticles[b]['px'].values,b2dfparticles[b]['py'].values,
                                                        betax2,betay2,xcut2,ycut2,colldict['cutoffAmpl'],colldict['emitMethod'],colldict['refEmxy'])
        
#   synchrotron motion
    if(switches['RFswitch']):
        for b in b1:
            b1np[b],x,px,y,py,t,pt,b1nLostDebunch[b],kept = frfupdate.rfupdate(b1dfparticles[b]['x'].values,b1dfparticles[b]['px'].values,
                                                                               b1dfparticles[b]['y'].values,b1dfparticles[b]['py'].values,
                                                                               b1dfparticles[b]['t'].values,b1dfparticles[b]['pt'].values,
                                                                               fmix1,b1nLostDebunch[b],2*np.pi/trev1,v001,ringparams['thib'],
                                                                               tcoeff1,fnharmb1,fnharm2b1,ringparams['vrfb1'],ringparams['vrf2b1'])

            b1dfparticles[b] =  b1dfparticles[b].head(kept)
            
        for b in b2:
            b2np[b],x,px,y,py,t,pt,b2nLostDebunch[b],kept = frfupdate.rfupdate(b2dfparticles[b]['x'].values,b2dfparticles[b]['px'].values,
                                                                               b2dfparticles[b]['y'].values,b2dfparticles[b]['py'].values,
                                                                               b2dfparticles[b]['t'].values,b2dfparticles[b]['pt'].values,
                                                                               fmix1,b2nLostDebunch[b],2*np.pi/trev2,v002,ringparams['thib'],
                                                                               tcoeff2,fnharmb2,fnharm2b2,ringparams['vrfb2'],ringparams['vrf2b2'])

            b2dfparticles[b] = b2dfparticles[b].head(kept)
                
#  betatron motion
    if(switches['betatronSwitch']):
        for b in b1:
            fsptranschrom.sptranschrom(b1dfparticles[b]['x'].values,b1dfparticles[b]['px'].values,b1dfparticles[b]['y'].values,
                                                      b1dfparticles[b]['py'].values,b1dfparticles[b]['pt'].values,
                                                      psix1,psiy1,coeffchromx1,coeffchromy1,tpdqmin1,ringparams['k2L1'],ringparams['k2Lskew1'])

            b1avglinex[b],b1avgliney[b] = fkeeptransprof.keeptransvprof(b1avglinex[b],b1avgliney[b],x,y,b1exdict[b],b1eydict[b],betax1,betay1)
        
        for b in b2:
            
            fsptranschrom.sptranschrom(b2dfparticles[b]['x'].values,b2dfparticles[b]['px'].values,b2dfparticles[b]['y'].values,
                                                      b2dfparticles[b]['py'].values,b2dfparticles[b]['pt'].values,
                                                      psix2,psiy2,coeffchromx2,coeffchromy2,tpdqmin2,ringparams['k2L2'],ringparams['k2Lskew2'])

            
            b2avglinex[b],b2avgliney[b] = fkeeptransprof.keeptransvprof(b2avglinex[b],b2avgliney[b],x,y,b2exdict[b],b2eydict[b],betax2,betay2)

#  radiation damping and quantum excitation
    if(switches['raddampSwitch']):
        for b in b1:
            x,px,y,py,pt,iseed,b1np[b] = frad.raddamp(b1dfparticles[b]['x'].values,b1dfparticles[b]['px'].values,b1dfparticles[b]['y'].values,
                                                    b1dfparticles[b]['py'].values,b1dfparticles[b]['pt'].values,
                                                    tradlong1,tradperp1,trev1,siglong1,sigperp1,ringparams['timeRatio'],12763)
        for b in b2:
            x,px,y,py,pt,iseed,b2np[b] = frad.raddamp(b2dfparticles[b]['x'].values,b2dfparticles[b]['px'].values,b2dfparticles[b]['y'].values,
                                        b2dfparticles[b]['py'].values,b2dfparticles[b]['pt'].values,
                                        tradlong2,tradperp2,trev2,siglong2,sigperp2,ringparams['timeRatio'],12763)
            
# ibs  
    if (switches['IBSswitch']):
        for b in b1:
            px,py,t,pt,ibsrow,iseed = fibs.ibslong(b1pnumber[b],b1dfparticles[b]['px'].values,b1dfparticles[b]['py'].values,b1dfparticles[b]['t'].values,
                                                        b1dfparticles[b]['pt'].values,b1exdict[b],b1eydict[b],ringparams['nMacro'],ringparams['circ'],
                                                        iwrite,ringparams['thib'],b1avgline[b],betar1,ringparams['gamma1'],
                                                        ringparams['gammat'],betax1,betay1,ringparams['aatomb1'],ringparams['qatomb1'],
                                                        trev1,np.array(dftfs1['BETX'].values,order='F'),np.array(dftfs1['BETY'].values,order='F'),
                                                        np.array(dftfs1['DX'].values,order='F'),np.array(dftfs1['DY'].values,order='F'),
                                                        np.array(dftfs1['DPX'].values,order='F'),np.array(dftfs1['DPY'].values,order='F'),
                                                        np.array(dftfs1['ALFX'].values,order='F'),np.array(dftfs1['ALFY'].values,order='F'),
                                                        np.array(dftfs1['L'].values,order='F'),gtab,ibsdict['coulombLog'],xmin,xmax,ymin,ymax,zmin,zmax,ap,ax,ay,
                                                        pnuminterp,ibsdict['ibsMethod'],ibsdict['coupleIBS'],nturns,kturns,eqtime,
                                                        ibsdict['fracibstot'],const.c,np.pi,ringparams['timeRatio'] ,12673)
            if (iwrite==1):
                b1ibsdict[b] = b1ibsdict[b].append(pd.DataFrame(ibsrow.reshape((1,5)),columns=b1ibsdict[b].columns)) 
        for b in b2:
            px,py,t,pt,ibsrow,iseed = fibs.ibslong(b2pnumber[b],b2dfparticles[b]['px'].values,b2dfparticles[b]['py'].values,b2dfparticles[b]['t'].values,
                                                        b2dfparticles[b]['pt'].values,b2exdict[b],b2eydict[b],ringparams['nMacro'],ringparams['circ'],
                                                        iwrite,ringparams['thib'],b2avgline[b],betar2,ringparams['gamma2'],
                                                        ringparams['gammat'],betax2,betay2,ringparams['aatomb2'],ringparams['qatomb2'],
                                                        trev2,np.array(dftfs2['BETX'].values,order='F'),np.array(dftfs2['BETY'].values,order='F'),
                                                        np.array(dftfs2['DX'].values,order='F'),np.array(dftfs2['DY'].values,order='F'),
                                                        np.array(dftfs2['DPX'].values,order='F'),np.array(dftfs2['DPY'].values,order='F'),
                                                        np.array(dftfs2['ALFX'].values,order='F'),np.array(dftfs2['ALFY'].values,order='F'),
                                                        np.array(dftfs2['L'].values,order='F'),gtab,ibsdict['coulombLog'],xmin,xmax,ymin,ymax,zmin,zmax,ap,ax,ay,
                                                        pnuminterp,ibsdict['ibsMethod'],ibsdict['coupleIBS'],nturns,kturns,eqtime,
                                                        ibsdict['fracibstot'],const.c,np.pi,ringparams['timeRatio'] ,12673)
            if (iwrite==1):
                b2ibsdict[b] = b2ibsdict[b].append(pd.DataFrame(ibsrow.reshape((1,5)),columns=b2ibsdict[b].columns))
# collimation
    if (switches['collimationSwitch']):
        for b in b1:
            collx,colly,momx,kept= fcollim.collimation(b1dfparticles[b]['x'].values,b1dfparticles[b]['px'].values,b1dfparticles[b]['y'].values,
                                                       b1dfparticles[b]['py'].values,b1dfparticles[b]['t'].values,b1dfparticles[b]['pt'].values,
                                                       b1nLostMom[b],b1nLostBeta[b],eqtime,betax1,betay1,colldict['refEmxy'],colldict['collimAvgSwitch'],
                                                       colldict['nSigCutBeta'],colldict['nSigCutMom'],colldict['betaxMom'],betar1,ringparams['gamma1'],
                                                       colldict['dispxMom'],kturns,nturns)
            if len(b1dfparticles[b])!=kept-1:
                print kept,
            b1dfparticles[b]       = b1dfparticles[b].head(kept-1)
            collx = collx[~np.any(collx==0, axis=1)]
            colly = colly[~np.any(colly==0, axis=1)]
            momx  = momx[~np.any(momx==0, axis=1)]
            if (len(collx) >0):
                b1dfcolimpxbeta[b] = b1dfcolimpxbeta[b].append(pd.DataFrame(collx,columns=b1dfcolimpxbeta[b].columns))
            if (len(colly) >0):
                b1dfcolimpybeta[b] = b1dfcolimpybeta[b].append(pd.DataFrame(colly,columns=b1dfcolimpxbeta[b].columns))
            if (len(momx) >0):
                b1dfcolimpxmom[b] = b1dfcolimpxmom[b].append(pd.DataFrame(momx,columns=b1dfcolimpxmom[b].columns))
 
        for b in b2:
            collx,colly,momx,kept= fcollim.collimation(b2dfparticles[b]['x'].values,b2dfparticles[b]['px'].values,b2dfparticles[b]['y'].values,
                                                       b2dfparticles[b]['py'].values,b2dfparticles[b]['t'].values,b2dfparticles[b]['pt'].values,
                                                       b2nLostMom[b],b2nLostBeta[b],eqtime,betax2,betay2,colldict['refEmxy'],colldict['collimAvgSwitch'],
                                                       colldict['nSigCutBeta'],colldict['nSigCutMom'],colldict['betaxMom'],betar2,ringparams['gamma2'],
                                                       colldict['dispxMom'],kturns,nturns)
            b2dfparticles[b]       = b2dfparticles[b].head(kept-1)
            collx = collx[~np.any(collx==0, axis=1)]
            colly = colly[~np.any(colly==0, axis=1)]
            momx  = momx[~np.any(momx==0, axis=1)]
            if (len(collx) >0):
                b2dfcolimpxbeta[b] = b2dfcolimpxbeta[b].append(pd.DataFrame(collx,columns=b2dfcolimpxbeta[b].columns))
            if (len(colly) >0):
                b2dfcolimpybeta[b] = b2dfcolimpybeta[b].append(pd.DataFrame(colly,columns=b2dfcolimpxbeta[b].columns))
            if (len(momx) >0):
                b2dfcolimpxmom[b] = b2dfcolimpxmom[b].append(pd.DataFrame(momx,columns=b2dfcolimpxmom[b].columns))
# blowup
    if (switches['blowupSwitch']):
        for b in b1:
            px,py,iseed = fblowup.blowup(b1dfparticles[b]['px'].values,b1dfparticles[b]['py'].values,12673,
                                           blowupdict['blowupMethod'],betax1,betay1,colldict['refEmxy'],blowupdict['pxKickFac'],blowupdict['pyKickFac'])
        for b in b2:
            px,py,iseed = fblowup.blowup(b2dfparticles[b]['px'].values,b2dfparticles[b]['py'].values,12673,
                                           blowupdict['blowupMethod'],betax2,betay2,colldict['refEmxy'],blowupdict['pxKickFac'],blowupdict['pyKickFac'])
# collision 
# code is quite long and semi-duplicate, maybe a function can be constructed - under investigation
    if (switches['collisionSwitch']):
        if (collisiondict['collRoutine']=='6a'):
            # counters for selecting correct hglassfacOld
            ip1counter = 0
            ip2counter = 0
            ip5counter = 0
            ip8counter = 0
            betass1 = OrderedDict()
            betass2 = OrderedDict()
            betass5 = OrderedDict()
            betass8 = OrderedDict()
            lum1 = OrderedDict()
            lum2 = OrderedDict()
            lum5 = OrderedDict()
            lum8 = OrderedDict()
            for k in tdc.keys():
                # for each ip
                for i in range(4):
                    # if there is a collision
                    if tdc[k][i] is not None:
                        # if IP1
                        if i == 0:
                            # get the corresponding hglassfac
                            hgold = hglassfacOld['ip1'][ip1counter]
                            # get the bunches that are colliding
                            b1bunch = b1[tdc[k][i][0]]
                            b2bunch = b2[tdc[k][i][1]]
                            # perform the collision
                            len1 = len(b1dfparticles[b1bunch]['x'].values)
                            len2 = len(b2dfparticles[b2bunch]['x'].values)
                            hglassfac,betass1['b1-'+str(b1bunch)+'-b2-'+str(b2bunch)],kept1,kept2,lum1out =\
                                fcollision.collision6a(collisiondict['sigI'],ringparams['thib'],collisiondict['longIntBins'],
                                                                               ringparams['timeRatio'],ringparams['circ'],np.pi,const.c,12673,hgold,
                                                                               b1nLostLum_ip1[b1bunch],b2nLostLum_ip1[b2bunch],
                                                                               ringparams['nMacro'],ringparams['nMacro'],
                                                                               b1pnumber[b1bunch].values[0],
                                                                               b1dfparticles[b1bunch]['x'].values,
                                                                               b1dfparticles[b1bunch]['px'].values,
                                                                               b1dfparticles[b1bunch]['y'].values,
                                                                               b1dfparticles[b1bunch]['py'].values,
                                                                               b1dfparticles[b1bunch]['t'].values,
                                                                               b1dfparticles[b1bunch]['pt'].values,
                                                                               b1exdict[b1bunch],b1eydict[b1bunch],betax1,betay1,bs1,lumimax1,theta1,
                                                                               b2pnumber[b2bunch].values[0],
                                                                               b2dfparticles[b2bunch]['x'].values,
                                                                               b2dfparticles[b2bunch]['px'].values,
                                                                               b2dfparticles[b2bunch]['y'].values,
                                                                               b2dfparticles[b2bunch]['py'].values,
                                                                               b2dfparticles[b2bunch]['t'].values,
                                                                               b2dfparticles[b2bunch]['pt'].values,
                                                                               b2exdict[b2bunch],b2eydict[b2bunch],betax2,betay2)

                            b1nLostLum_ip1[b1bunch] = b1nLostLum_ip1[b1bunch] + len1 - kept1
                            b2nLostLum_ip1[b2bunch] = b2nLostLum_ip1[b2bunch] + len2 - kept2
                            
                            b1nLostLum[b1bunch] = b1nLostLum[b1bunch] + len1 - kept1
                            b2nLostLum[b2bunch] = b2nLostLum[b2bunch] + len2 - kept2
                            b1dfparticles[b1bunch] = b1dfparticles[b1bunch].head(kept1)
                            b2dfparticles[b2bunch] = b2dfparticles[b2bunch].head(kept2)
                            lum1['b1-'+str(b1bunch)+'-b2-'+str(b2bunch)] = lum1out
                            hglassfacOld['ip1'][ip1counter] = hglassfac
                            ip1counter += ip1counter

                        # IP2
                        if i == 1:
                            # get the corresponding hglassfac
                            hgold = hglassfacOld['ip2'][ip2counter]
                            # get the bunches that are colliding
                            b1bunch = b1[tdc[k][i][0]]
                            b2bunch = b2[tdc[k][i][1]]
                            # perform the collision
                            len1 = len(b1dfparticles[b1bunch]['x'].values)
                            len2 = len(b2dfparticles[b2bunch]['x'].values)
                            hglassfac,betass2['b1-'+str(b1bunch)+'-b2-'+str(b2bunch)],kept1,kept2,lum2out =\
                                fcollision.collision6a(collisiondict['sigI'],ringparams['thib'],collisiondict['longIntBins'],
                                                                               ringparams['timeRatio'],ringparams['circ'],np.pi,const.c,12673,hgold,
                                                                               b1nLostLum_ip2[b1bunch],b2nLostLum_ip2[b2bunch],
                                                                               ringparams['nMacro'],ringparams['nMacro'],
                                                                               b1pnumber[b1bunch].values[0],
                                                                               b1dfparticles[b1bunch]['x'].values,
                                                                               b1dfparticles[b1bunch]['px'].values,
                                                                               b1dfparticles[b1bunch]['y'].values,
                                                                               b1dfparticles[b1bunch]['py'].values,
                                                                               b1dfparticles[b1bunch]['t'].values,
                                                                               b1dfparticles[b1bunch]['pt'].values,
                                                                               b1exdict[b1bunch],b1eydict[b1bunch],betax1,betay1,bs2,lumimax2,theta2,
                                                                               b2pnumber[b2bunch].values[0],
                                                                               b2dfparticles[b2bunch]['x'].values,
                                                                               b2dfparticles[b2bunch]['px'].values,
                                                                               b2dfparticles[b2bunch]['y'].values,
                                                                               b2dfparticles[b2bunch]['py'].values,
                                                                               b2dfparticles[b2bunch]['t'].values,
                                                                               b2dfparticles[b2bunch]['pt'].values,
                                                                               b2exdict[b2bunch],b2eydict[b2bunch],betax2,betay2)

                            b1nLostLum_ip2[b1bunch] = b1nLostLum_ip2[b1bunch] + len1 - kept1
                            b2nLostLum_ip2[b2bunch] = b2nLostLum_ip2[b2bunch] + len2 - kept2
                            b1nLostLum[b1bunch] = b1nLostLum[b1bunch] + len1 - kept1
                            b2nLostLum[b2bunch] = b2nLostLum[b2bunch] + len2 - kept2
                            b1dfparticles[b1bunch] = b1dfparticles[b1bunch].head(kept1)
                            b2dfparticles[b2bunch] = b2dfparticles[b2bunch].head(kept2)
                            lum2['b1-'+str(b1bunch)+'-b2-'+str(b2bunch)] = lum2out
                            hglassfacOld['ip2'][ip2counter] = hglassfac
                            ip2counter += ip2counter
                        # IP5
                        if i == 2:
                            # get the corresponding hglassfac
                            hgold = hglassfacOld['ip5'][ip5counter]
                            # get the bunches that are colliding
                            b1bunch = b1[tdc[k][i][0]]
                            b2bunch = b2[tdc[k][i][1]]
                            # perform the collision
                            len1 = len(b1dfparticles[b1bunch]['x'].values)
                            len2 = len(b2dfparticles[b2bunch]['x'].values)
                            hglassfac,betass5['b1-'+str(b1bunch)+'-b2-'+str(b2bunch)],kept1,kept2,lum5out =\
                                fcollision.collision6a(collisiondict['sigI'],ringparams['thib'],collisiondict['longIntBins'],
                                                                               ringparams['timeRatio'],ringparams['circ'],np.pi,const.c,12673,hgold,
                                                                               b1nLostLum_ip5[b1bunch],b2nLostLum_ip5[b2bunch],
                                                                               ringparams['nMacro'],ringparams['nMacro'],
                                                                               b1pnumber[b1bunch].values[0],
                                                                               b1dfparticles[b1bunch]['x'].values,
                                                                               b1dfparticles[b1bunch]['px'].values,
                                                                               b1dfparticles[b1bunch]['y'].values,
                                                                               b1dfparticles[b1bunch]['py'].values,
                                                                               b1dfparticles[b1bunch]['t'].values,
                                                                               b1dfparticles[b1bunch]['pt'].values,
                                                                               b1exdict[b1bunch],b1eydict[b1bunch],betax1,betay1,bs5,lumimax5,theta5,
                                                                               b2pnumber[b2bunch].values[0],
                                                                               b2dfparticles[b2bunch]['x'].values,
                                                                               b2dfparticles[b2bunch]['px'].values,
                                                                               b2dfparticles[b2bunch]['y'].values,
                                                                               b2dfparticles[b2bunch]['py'].values,
                                                                               b2dfparticles[b2bunch]['t'].values,
                                                                               b2dfparticles[b2bunch]['pt'].values,
                                                                               b2exdict[b2bunch],b2eydict[b2bunch],betax2,betay2)
                                                       
                            b1nLostLum_ip5[b1bunch] = b1nLostLum_ip5[b1bunch] + len1 - kept1
                            b2nLostLum_ip5[b2bunch] = b2nLostLum_ip5[b2bunch] + len2 - kept2
                            b1nLostLum[b1bunch] = b1nLostLum[b1bunch] + len1 - kept1
                            b2nLostLum[b2bunch] = b2nLostLum[b2bunch] + len2 - kept2
                            b1dfparticles[b1bunch] = b1dfparticles[b1bunch].head(kept1)
                            b2dfparticles[b2bunch] = b2dfparticles[b2bunch].head(kept2)
                            lum5['b1-'+str(b1bunch)+'-b2-'+str(b2bunch)] = lum5out
                            hglassfacOld['ip5'][ip5counter] = hglassfac
                            ip5counter += ip5counter
                        # IP8
                        if i == 3:
                            # get the corresponding hglassfac
                            hgold = hglassfacOld['ip8'][ip8counter]
                            # get the bunches that are colliding
                            b1bunch = b1[tdc[k][i][0]]
                            b2bunch = b2[tdc[k][i][1]]
                            # perform the collision
                            len1 = len(b1dfparticles[b1bunch]['x'].values)
                            len2 = len(b2dfparticles[b2bunch]['x'].values)
                            hglassfac,betass8['b1-'+str(b1bunch)+'-b2-'+str(b2bunch)],kept1,kept2,lum8out =\
                                fcollision.collision6a(collisiondict['sigI'],ringparams['thib'],collisiondict['longIntBins'],
                                                                               ringparams['timeRatio'],ringparams['circ'],np.pi,const.c,12673,hgold,
                                                                               b1nLostLum_ip8[b1bunch],b2nLostLum_ip8[b2bunch],
                                                                               ringparams['nMacro'],ringparams['nMacro'],
                                                                               b1pnumber[b1bunch].values[0],
                                                                               b1dfparticles[b1bunch]['x'].values,
                                                                               b1dfparticles[b1bunch]['px'].values,
                                                                               b1dfparticles[b1bunch]['y'].values,
                                                                               b1dfparticles[b1bunch]['py'].values,
                                                                               b1dfparticles[b1bunch]['t'].values,
                                                                               b1dfparticles[b1bunch]['pt'].values,
                                                                               b1exdict[b1bunch],b1eydict[b1bunch],betax1,betay1,bs8,lumimax8,theta8,
                                                                               b2pnumber[b2bunch].values[0],
                                                                               b2dfparticles[b2bunch]['x'].values,
                                                                               b2dfparticles[b2bunch]['px'].values,
                                                                               b2dfparticles[b2bunch]['y'].values,
                                                                               b2dfparticles[b2bunch]['py'].values,
                                                                               b2dfparticles[b2bunch]['t'].values,
                                                                               b2dfparticles[b2bunch]['pt'].values,
                                                                               b2exdict[b2bunch],b2eydict[b2bunch],betax2,betay2)

                            b1nLostLum_ip8[b1bunch] = b1nLostLum_ip8[b1bunch] + len1 - kept1
                            b2nLostLum_ip8[b2bunch] = b2nLostLum_ip8[b2bunch] + len2 - kept2
                            b1nLostLum[b1bunch] = b1nLostLum[b1bunch] + len1 - kept1
                            b2nLostLum[b2bunch] = b2nLostLum[b2bunch] + len2 - kept2
                            b1dfparticles[b1bunch] = b1dfparticles[b1bunch].head(kept1)
                            b2dfparticles[b2bunch] = b2dfparticles[b2bunch].head(kept2)
                            lum8['b1-'+str(b1bunch)+'-b2-'+str(b2bunch)] = lum8out
                            hglassfacOld['ip8'][ip8counter] = hglassfac
                            ip8counter += ip8counter
        if (collisiondict['collRoutine']=='1d'):
            print 'under construction'
            quit()
    #first turn already written
    if ((kturns==1)&(switches['collisionSwitch']==1.0)):
        iwrite=0
        # ip1
        ip1counter = 0
        for k in lumip1.keys():
            lumip1[k] = cwrite.writeLumi(lum1[k],hglassfacOld['ip1'][ip1counter],lumip1[k],kturns,nturns,eqtime,betass1[k])
            ip1counter += ip1counter
            
        # ip2
        ip2counter = 0
        for k in lumip2.keys():
            lumip2[k] = cwrite.writeLumi(lum2[k],hglassfacOld['ip2'][ip2counter],lumip2[k],kturns,nturns,eqtime,betass2[k])
            ip2counter += ip2counter
        # ip5
        ip5counter = 0
        for k in lumip5.keys():
            lumip5[k] = cwrite.writeLumi(lum5[k],hglassfacOld['ip5'][ip5counter],lumip5[k],kturns,nturns,eqtime,betass5[k])
            ip5counter += ip5counter
        # ip8
        ip8counter = 0
        for k in lumip8.keys():
            lumip8[k] = cwrite.writeLumi(lum8[k],hglassfacOld['ip8'][ip8counter],lumip8[k],kturns,nturns,eqtime,betass8[k])
            ip8counter += ip8counter
#     iwrite=1
    if (iwrite==1):
        #call writemomentsshort(np1,y1,py1,t1,pt1,x1,px1,1)
        # call writemomentsshort(np2,y2,py2,t2,pt2,x2,px2,2)
        if switches['writeAllCoordSwitch']:   #     write all coordinates to file
            for b in b1:
                b1dfparticles[b].to_csv('all-coord-b1-bunch-'+str(b)+'-turn-'+str(kturns)+'.csv',index=False)
            for b in b2:
                b2dfparticles[b].to_csv('all-coord-b2-bunch-'+str(b)+'-turn-'+str(kturns)+'.csv',index=False)
        if switches['writeMountSwitch']:
            for b in b1:
                b1dfwritemount[b] =cwrite.mountainr(ringparams['nMacro'],b1pnumber[b],b1avgline[b],ringparams['qatomb1'],
                                         ringparams['nwrite'],kturns,ringparams['thib'],ibsdict['nbins'],b1dfwritemount[b])
                 
                b1dfwritemounttrans[b] = cwrite.mountTransv(b1avglinex[b],b1avgliney[b],b1exdict[b],b1eydict[b],betax1,betay1,ibsdict['nbins'],
                                                kturns,b1dfwritemounttrans[b])
            for b in b2:
                b2dfwritemount[b] =cwrite.mountainr(ringparams['nMacro'],b2pnumber[b],b2avgline[b],ringparams['qatomb2'],
                                         ringparams['nwrite'],kturns,ringparams['thib'],ibsdict['nbins'],b2dfwritemount[b])
                b2dfwritemounttrans[b] = cwrite.mountTransv(b2avglinex[b],b2avgliney[b],b2exdict[b],b2eydict[b],betax2,betay2,ibsdict['nbins'],
                                                kturns,b2dfwritemounttrans[b])
        for b in b1:
            b1emitdict[b] = cwrite.writeEmi(b1exdict[b],b1eydict[b],b1dfparticles[b]['t'],b1dfparticles[b]['pt'],betar1,ringparams['gamma1'],
                              ringparams['aatomb1'],ringparams['qatomb1'],eqtime,kturns,ringparams['nturns'],b1emitdict[b])
            b1intdict[b],b1nLostLumSum[b],b1nLostDebunchSum[b],b1nLostBetaSum[b],b1nLostMomSum[b]=\
                            cwrite.writeNb(b1dfparticles[b],ringparams['nMacro'],b1pnumber[b].values[0], b1nLostLum[b],b1nLostLumSum[b],
                                           b1nLostDebunch[b],b1nLostDebunchSum[b],b1nLostBeta[b],b1nLostBetaSum[b],
                                           b1nLostMom[b],b1nLostMomSum[b],kturns,ringparams['nturns'],eqtime,b1intdict[b])
        for b in b2:
            b2emitdict[b] = cwrite.writeEmi(b2exdict[b],b2eydict[b],b2dfparticles[b]['t'],b2dfparticles[b]['pt'],betar2,ringparams['gamma2'],
                              ringparams['aatomb2'],ringparams['qatomb2'],eqtime,kturns,ringparams['nturns'],b2emitdict[b])
            b2intdict[b],b2nLostLumSum[b],b2nLostDebunchSum[b],b2nLostBetaSum[b],b2nLostMomSum[b]=\
                            cwrite.writeNb(b2dfparticles[b],ringparams['nMacro'],b2pnumber[b].values[0], b2nLostLum[b],b2nLostLumSum[b],
                                           b2nLostDebunch[b],b2nLostDebunchSum[b],b2nLostBeta[b],b2nLostBetaSum[b],
                                           b2nLostMom[b],b2nLostMomSum[b],kturns,ringparams['nturns'],eqtime,b2intdict[b])
        if (switches['collisionSwitch']):
            # ip1
            ip1counter = 0
            for k in lumip1.keys():
                lumip1[k] = cwrite.writeLumi(lum1[k],hglassfacOld['ip1'][ip1counter],lumip1[k],kturns,nturns,eqtime,betass1[k])
                ip1counter += ip1counter

            # ip2
            ip2counter = 0
            for k in lumip2.keys():
                lumip2[k] = cwrite.writeLumi(lum2[k],hglassfacOld['ip2'][ip2counter],lumip2[k],kturns,nturns,eqtime,betass2[k])
                ip2counter += ip2counter
            # ip5
            ip5counter = 0
            for k in lumip5.keys():
                lumip5[k] = cwrite.writeLumi(lum5[k],hglassfacOld['ip5'][ip5counter],lumip5[k],kturns,nturns,eqtime,betass5[k])
                ip5counter += ip5counter
            # ip8
            ip8counter = 0
            for k in lumip8.keys():
                lumip8[k] = cwrite.writeLumi(lum8[k],hglassfacOld['ip8'][ip8counter],lumip8[k],kturns,nturns,eqtime,betass8[k])
                ip8counter += ip8counter
        print b1nLostLum,b2nLostLum
        # reset counters
        for k in b1nLostLum.keys():
            b1nLostLum[k]=0
        for k in b1nLostDebunch.keys():
            b1nLostDebunch[k]=0
        for k in b1nLostMom.keys():
            b1nLostMom[k]=0
        for k in b1nLostBeta.keys():
            b1nLostBeta[k]=0
        for k in b1nLostLum_ip1.keys():
            b1nLostLum_ip1[k]=0
        for k in b1nLostLum_ip2.keys():
            b1nLostLum_ip2[k]=0
        for k in b1nLostLum_ip5.keys():
            b1nLostLum_ip5[k]=0
        for k in b1nLostLum_ip8.keys():
            b1nLostLum_ip8[k]=0
        for k in b2nLostLum.keys():
            b2nLostLum[k]=0
        for k in b2nLostDebunch.keys():
            b2nLostDebunch[k]=0
        for k in b2nLostMom.keys():
            b2nLostMom[k]=0
        for k in b2nLostBeta.keys():
            b2nLostBeta[k]=0
        for k in b2nLostLum_ip1.keys():
            b2nLostLum_ip1[k]=0
        for k in b2nLostLum_ip2.keys():
            b2nLostLum_ip2[k]=0
        for k in b2nLostLum_ip5.keys():
            b2nLostLum_ip5[k]=0
        for k in b2nLostLum_ip8.keys():
            b2nLostLum_ip8[k]=0
            
print 'Writing output files'
for b in b1:
    b1intdict[b].to_csv(fn+'-b1-'+str(b)+'-int.csv')
    b1emitdict[b].to_csv(fn+'-b1-'+str(b)+'-emit.csv')
    b1ibsdict[b].to_csv(fn+'-b1-'+str(b)+'-ibs.csv')
print 'b1 written'
for b in b2:
    b2intdict[b].to_csv(fn+'-b2-'+str(b)+'-int.csv')
    b2emitdict[b].to_csv(fn+'-b2-'+str(b)+'-emit.csv')
    b2ibsdict[b].to_csv(fn+'-b2-'+str(b)+'-ibs.csv')
print 'b2 written'

for k in lumip1.keys():
    lumip1[k].to_csv(fn+'-ip1-'+str(k)+'-lumi.csv')
for k in lumip2.keys():
    lumip2[k].to_csv(fn+'-ip2-'+str(k)+'-lumi.csv')
for k in lumip5.keys():
    lumip5[k].to_csv(fn+'-ip5-'+str(k)+'-lumi.csv')
for k in lumip8.keys():
    lumip8[k].to_csv(fn+'-ip8-'+str(k)+'-lumi.csv')
    
print 'took :', time.clock()-st
print 'end at:',time.ctime()
