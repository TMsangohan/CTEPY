# coding: utf-8
import pandas as pd
import numpy as np
from collections import OrderedDict
from scipy import constants as const

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

switcheskeys = ['RFswitch','betatronSwitch','raddampSwitch','IBSswitch','collimationSwitch','blowupSwitch','collisionSwitch','writeAllCoordSwitch', 'writeMountSwitch']
ringparamskeys = ['aatomb1','aatomb2','thib', 'nturns', 'nMacro', 'timeRatio', 'nwrite',
                  'vrfb1', 'nharmb1', 'vrf2b1', 'nharm2b1', 'vrfb2', 'nharmb2', 'vrf2b2', 'nharm2b2', 
                  'dqmin1', 'dqmin2', 'k2L1', 'k2Lskew2', 'k2L2', 'k2Lskew1', 'tfsb1', 'tfsb2']
startingconditionskeys = ['tauhat1','tauhat2', 'longcoor', 'transcoor', 'bunchLenPrecis', 'power', 'alint']
radiationdampingkeys = ['radMethod', 'tradlong', 'tradperp', 'siglong', 'sigperp', 'rho0']
ibsdictkeys = ['ibsMethod', 'coupleIBS', 'coulombLog', 'fracibstot', 'nbins', 'piwinski', 'bane', 'intfile']
colldictkeys = ['refEmxy', 'cutoffAmpl', 'collimAvgSwitch', 'emitMethod', 'nSigCutBeta', 'nSigCutMom', 'betaxMom', 'dispxMom']
blowupdictkeys = ['pxKickFac', 'pyKickFac', 'blowupMethod']
collisiondictkeys = ['collRoutine', 'longIntBins', 'ips_leveling', 'ips_leveling_values']

def getswitches(infn):
    inputdf  = pd.read_csv(infn,index_col=False)
    switches = inputdf.ix[0][switcheskeys].to_dict()
    return switches

def getringparams(infn):
    inputdf    = pd.read_csv(infn,index_col=False)
    ringparams = inputdf.ix[0][ringparamskeys].to_dict()
    print ringparams['tfsb1']
    protonmass                      = const.physical_constants['proton mass energy equivalent in MeV'][0]/1000 # GeV
    ringparams['proton']            = protonmass
    ringparams['ion']               = 193.7291748489224
    
    dftfsb1head                     = pd.read_csv(ringparams['tfsb1'],nrows=44,delim_whitespace=True)
    ringparams['tunex1']            = float(dftfsb1head[dftfsb1head['NAME']=='Q1']['TWISS'].values[0])
    ringparams['tuney1']            = float(dftfsb1head[dftfsb1head['NAME']=='Q2']['TWISS'].values[0])
    ringparams['chromx1']           = float(dftfsb1head[dftfsb1head['NAME']=='DQ1']['TWISS'].values[0])
    ringparams['chromy1']           = float(dftfsb1head[dftfsb1head['NAME']=='DQ2']['TWISS'].values[0])

    dftfsb2head                     = pd.read_csv(ringparams['tfsb2'],nrows=44,delim_whitespace=True)
    ringparams['tunex2']            = float(dftfsb2head[dftfsb2head['NAME']=='Q1']['TWISS'].values[0])
    ringparams['tuney2']            = float(dftfsb2head[dftfsb2head['NAME']=='Q2']['TWISS'].values[0])
    ringparams['chromx2']           = float(dftfsb2head[dftfsb2head['NAME']=='DQ1']['TWISS'].values[0])
    ringparams['chromy2']           = float(dftfsb2head[dftfsb2head['NAME']=='DQ2']['TWISS'].values[0])

    energy                          = 1e9 * float(dftfsb1head[dftfsb1head['NAME']=='ENERGY']['TWISS'].values[0])/float(dftfsb1head[dftfsb1head['NAME']=='CHARGE']['TWISS'].values[0])
    ringparams['gammat']            = float(dftfsb1head[dftfsb1head['NAME']=='GAMMATR']['TWISS'].values[0])
    ringparams['circ']              = float(dftfsb1head[dftfsb1head['NAME']=='LENGTH']['TWISS'].values[0])
    ringparams['qatomb1']           = float(dftfsb1head[dftfsb1head['NAME']=='CHARGE']['TWISS'].values[0])
    ringparams['qatomb2']           = float(dftfsb2head[dftfsb2head['NAME']=='CHARGE']['TWISS'].values[0])
    ringparams['gamma1']            = float(dftfsb1head[dftfsb1head['NAME']=='GAMMA']['TWISS'].values[0])
    ringparams['gamma2']            = float(dftfsb2head[dftfsb2head['NAME']=='GAMMA']['TWISS'].values[0])

    return ringparams,dftfsb1head,dftfsb2head,energy

def getradiationdamping(infn):
    inputdf    = pd.read_csv(infn,index_col=False)
    radiationdamping = inputdf.ix[0][radiationdampingkeys].to_dict()
    return radiationdamping

def getibsdict(infn):
    inputdf    = pd.read_csv(infn,index_col=False)
    ibsdict = inputdf.ix[0][ibsdictkeys].to_dict()
    return ibsdict

def getcolldict(infn):
    inputdf    = pd.read_csv(infn,index_col=False)
    colldict = inputdf.ix[0][colldictkeys].to_dict()
    return colldict

def getblowupdict(infn):
    inputdf    = pd.read_csv(infn,index_col=False)
    blowupdict = inputdf.ix[0][blowupdictkeys].to_dict()
    return blowupdict

def getcollisiondict(infn,ringparams):
    inputdf    = pd.read_csv(infn,index_col=False)
    collisiondict = inputdf.ix[0][collisiondictkeys].to_dict()
    # settings for collision routine
    # note : 1d is not implemented yet
    # VERY IMPORTANT : DON'T FORGET TO SET THE CORRECT CROSS SECTION
    collisiondict['collRoutine']       = '6a'           #(1d is slow but without assumptions on distributions, 6a is fast with assumed Gaussian transverse)
    if ringparams['aatomb1'] == 208:
        if ringparams['aatomb2'] == 208:
            collisiondict['sigI'] = 515.          # (cross section for particle removal in collisions) 
        if ringparams['aatomb2'] == 1:
            collisiondict['sigI'] = 2.0
    else:
        if ringparams['aatomb2'] == 208:
            collisiondict['sigI'] = 2.0          # (cross section for particle removal in collisions) 
        if ringparams['aatomb2'] == 1:
            collisiondict['sigI'] = 1.0

    return collisiondict

def setv00(ringparams):
    if (ringparams['aatomb'+str(1)] == 208):
        v001 = ringparams['ion']/  ringparams['qatomb'+str(1)] * 1e9
    else:
            v001 = ringparams['proton']/  ringparams['qatomb'+str(1)] *1e9

    if (ringparams['aatomb'+str(2)] == 208):
            v002 = ringparams['ion']/  ringparams['qatomb'+str(2)] * 1e9
    else:
            v002 = ringparams['proton']/  ringparams['qatomb'+str(2)] *1e9
    return v001,v002


def getstartingconditions(infn,b1fn,b2fn):
    cols       = ['bucket','exn','eyn','npart','blen']
    inputdf    = pd.read_csv(infn,index_col=False)
    
    startingconditions = inputdf.ix[0][startingconditionskeys].to_dict()
    
    b1df = pd.read_csv(b1fn,delim_whitespace=True,names=cols)
    b2df = pd.read_csv(b2fn,delim_whitespace=True,names=cols)
    return startingconditions,b1df,b2df

# function to get energyspread from bunchlength and RF voltage
def getsigpfromvrfandsigs(sigs,trev,h,vrf,eta,beta,eb):
    return 2*np.pi*sigs/trev*np.sqrt(h*vrf*eta/(2*np.pi*beta*eb))/eta/const.c

# reading bane or interpolat files for ibs if necessary
def readibsinputfiles(switches,ibsdict):
    dfbane       = pd.DataFrame()
    dfinterpolat = pd.DataFrame()
    if (switches['IBSswitch']==True) & (ibsdict['ibsMethod']=='baneApprox'):
        dfbane = pd.read_csv(ibsdict['bane'],delim_whitespace=True)
    if (switches['IBSswitch']==True) & (ibsdict['ibsMethod']=='interpolat'):
        dfinterpolat = pd.read_csv(ibsdict['int'],delim_whitespace=True)
    return dfbane, dfinterpolat

# benchmarked agains 'injection schemes' viewer for the colliding bunches
# beam1 : starts from 0 and bunch positions move down, moves clockwise
# beam2 : starts from 0 and bunch positions move down, moves counter-clockwise
# example ip2 : the bunch at b1 0 will move to 3564-445.5 while bunch at b2 891 will move to 445.5 and collide at ip2
# ip1 is at 0 for both beams
# ip2 is at 445.5 for beam2 and at 3564-445.5 for beam1
# ip5 is at 1782 for both beams
# ip8 is at 447 for beam 1 and 3117 for beam 2
def genlumikeys(b1,b2):
    tdc = OrderedDict()
    lumip1keys = []
    lumip2keys = []
    lumip5keys = []
    lumip8keys = []
    lumip1 = OrderedDict()
    lumip2 = OrderedDict()
    lumip5 = OrderedDict()
    lumip8 = OrderedDict()
    
    for i in range(3564*2):
        res = [None]*4
        bb1 = [(b-i/2.)%3564 for b in b1]
        bb2 = [(b-i/2.)%3564 for b in b2]
        for j in range(len(b1)):
            for k in range(len(b2)):
                if (bb1[j]==0) and (bb2[k]==0):
                    res[0]=[j,k]
                if (bb1[j]==3564-445.5) and (bb2[k]==445.5):
                    res[1]=[j,k]
                if (bb1[j]==1782) and (bb2[k]==1782):
                    res[2]=[j,k]
                if (bb1[j]==447) and (bb2[k]==3117):
                    res[3]=[j,k]
        if res != [None]*4:
            tdc[i/2.] = res 
            if res[0] is not None:
                lumip1keys.append('b1-'+str(b1[res[0][0]])+'-b2-'+str(b2[res[0][1]]))
            if res[1] is not None:
                lumip2keys.append('b1-'+str(b1[res[1][0]])+'-b2-'+str(b2[res[1][1]]))
            if res[2] is not None:
                lumip5keys.append('b1-'+str(b1[res[2][0]])+'-b2-'+str(b2[res[2][1]]))
            if res[3] is not None:
                lumip8keys.append('b1-'+str(b1[res[3][0]])+'-b2-'+str(b2[res[3][1]]))
    for k in lumip1keys:
        lumip1[k]= pd.DataFrame(columns= ['sim.turn','t(hours)', 'L(cm^-2 s^-1)','red factor','betas'])
    for k in lumip2keys:
        lumip2[k]= pd.DataFrame(columns= ['sim.turn','t(hours)', 'L(cm^-2 s^-1)','red factor','betas'])
    for k in lumip5keys:
        lumip5[k]= pd.DataFrame(columns= ['sim.turn','t(hours)', 'L(cm^-2 s^-1)','red factor','betas'])
    for k in lumip8keys:
        lumip8[k]= pd.DataFrame(columns= ['sim.turn','t(hours)', 'L(cm^-2 s^-1)','red factor','betas'])
    return tdc,lumip1,lumip2,lumip5,lumip8

# function that constructs dictionaries containing the bunch encounters
# given a list of bunchslots for beam1 and 2
# it also returns a list with the total number of collisions in each IP
# function benchmarked against 'Filling Scheme' viewer application

def getIPenc(b1bunches,b2bunches):
    B1 = OrderedDict()
    for b1 in b1bunches:
        B1[b1] = [None]*4
        if b1 in b2bunches:
            B1[b1][0] = b1
            B1[b1][2] = b1
        if b1 in (b2bunches-891)%3564:
            B1[b1][1] = (b1+891)%3564
        if b1 in (b2bunches-6234)%3564:
            B1[b1][3] = (b1+6234)%3564
        
    B2 = OrderedDict()
    for b2 in b2bunches:
        B2[b2] = [None]*4
        if b2 in b1bunches:
            B2[b2][0] = b2
            B2[b2][2] = b2
        if b2 in (b1bunches+891)%3564:
            B2[b2][1] = (b2-891)%3564
        if b2 in (b1bunches+6234)%3564:
            B2[b2][3] = (b2-6234)%3564
    coll = [sum(1 for k in B2.keys() if B2[k][0] is not None),
           sum(1 for k in B2.keys() if B2[k][1] is not None),
           sum(1 for k in B2.keys() if B2[k][2] is not None),
           sum(1 for k in B2.keys() if B2[k][3] is not None)]
    return B1,B2,coll

# function to build a list of bunches for beam1 and beam2 that are dependent on eachother through collisions
# input : 
# 1) a list of two lists of bunches, one for beam1 and one for beam2 of input bunches,
# in other words we are interested in these bunches but we need to track on which other bunches they depend 
# 2) two lists of bunch encounters = output of the getIPenc function from above
# returns a list of two sets of bunchslots, one for beam1 and one for beam2
# In other words their evolution is intertwined through collisions.

def updateenc(inlist,bunchencb1,bunchencb2):
    rellist = [set(inlist[0]),set(inlist[1])]
    currentlenb1 = len(rellist[0])
    currentlenb2 = len(rellist[1])
    for bunchone in rellist[0]:
        l1 = [b for b in bunchencb1[bunchone] if b is not None]
        rellist[1] = rellist[1] | set(l1)
        print rellist
    for bunchtwo in rellist[1]:
        #print bunchtwo
        if bunchtwo is not None:
            print bunchtwo,bunchencb2[bunchtwo]
        l2 = [b for b in bunchencb2[bunchtwo] if b is not None]
        rellist[0] = rellist[0] | set(l2)
        #print rellist
    if (len(rellist[0]) > currentlenb1) or (len(rellist[1]) > currentlenb2):
        rellist = updateenc(rellist,bunchencb1,bunchencb2)
        return rellist
    else:
        return rellist
    
def addparticles(df,tauhat,tcoeff,trev,v00,fnharm,fnharm2,vrf1,vrf2,gamma,blenprecision,nmacro,intmethod,cutoffamp,betax,betay):
    t,pt = faddpart.addlongitudinal(tauhat,tcoeff,2*np.pi/trev,v00,fnharm,fnharm2,vrf1,vrf2,
                                np.pi,5,.75,12763,gamma,df['rmsDelta'],df['blen'],blenprecision,const.c,nmacro,intmethod)
    x,px,y,py = faddpart.addtransverse(nmacro,cutoffamp,cutoffamp,betax,betay,df['ex'],df['ey'],12673)
    return t,pt,x,px,y,py

def writedistr(avglinedict,avglinexdict,avglineydict,exdict,eydict,nmacro,npart,qatom,nwrite,kturns,thib,nbins,betax,betay):
    for k in avglinedict.keys():
        dfwritemount = pd.DataFrame(columns = cwrite.mountainnrcols)
        dfwritemount = cwrite.mountainr(nmacro,npart,avglinedict[k],qatom,
                                         nwrite,kturns,thib,nbins,dfwritemount)
        dfwritemounttrans = pd.DataFrame(columns = cwrite.mounttransvcols)
        dfwritemounttrans = cwrite.mountTransv(avglinexdict[k],avglineydict[k],exdict[k],eydict[k],betax,betay,nbins,kturns,dfwritemounttrans)
    
