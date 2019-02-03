import numpy as np
import pandas as pd
import datetime
import time
from collections import OrderedDict
from itertools import groupby
from operator import itemgetter
import functools
import multiprocessing
from pandas import HDFStore
from pandas import read_hdf
import dask.dataframe as dd
import dask
from scipy import constants as const
from fnmatch import fnmatch, fnmatchcase
import glob
from random import randint

# naming - set for output file naming - see at end of program
fn = 'cte_run2_2016_pPb_6500_col1258_'

# set general parameters
ringparams                      = OrderedDict()
ringparams['aatomb1']           = 1         # atomic number of particle in beam 1
ringparams['aatomb2']           = 208         # atomic number of particle in beam 2

ringparams['nturns']            = 1000       # number of turns to simulate in the machine
ringparams['nMacro']            = 2000      # number of simulated particles
ringparams['timeRatio']         = 20000       # (=real turns/sim.turns)
ringparams['nwrite']            = 100         # step interval in turns to calculate output variables and write them to files

ringparams['vrfb1']             = 0.0
ringparams['nharmb1']           = 360
ringparams['vrf2b1']            = -12.0e6
ringparams['nharm2b1']          = 35640

ringparams['vrfb2']             = 0.0
ringparams['nharmb2']           = 360
ringparams['vrf2b2']            = -12.0e6
ringparams['nharm2b2']          = 35640
# Some of the parameters are determined directly from tfs files
# so it is very important to set the correct ones for BOTH beams !!!
ringparams['tfsb1']             = 'lhcb1-pPb-6500.tfs'
ringparams['tfsb2']             = 'lhcb2-pPb-6500.tfs'
# getting some of the settings from the tfs files
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

energy                          = float(dftfsb1head[dftfsb1head['NAME']=='ENERGY']['TWISS'].values[0])/float(dftfsb1head[dftfsb1head['NAME']=='CHARGE']['TWISS'].values[0])
ringparams['gammat']            = float(dftfsb1head[dftfsb1head['NAME']=='GAMMATR']['TWISS'].values[0])
ringparams['circ']              = float(dftfsb1head[dftfsb1head['NAME']=='LENGTH']['TWISS'].values[0])
ringparams['qatomb1']           = float(dftfsb1head[dftfsb1head['NAME']=='CHARGE']['TWISS'].values[0])
ringparams['qatomb2']           = float(dftfsb2head[dftfsb2head['NAME']=='CHARGE']['TWISS'].values[0])
ringparams['gamma1']            = float(dftfsb1head[dftfsb1head['NAME']=='GAMMA']['TWISS'].values[0])
ringparams['gamma2']            = float(dftfsb2head[dftfsb2head['NAME']=='GAMMA']['TWISS'].values[0])

# setting of initial parameters for the bunches to simulate
# note that the simulation takes a bunch length and RF voltage
# to calculate the energy spread for each bunch
#
# If collisions active :
# ----------------------
# bunch 1 is the main bunch to track so cannot be omitted
# bunch 2 is the bunch seen by the main bunch in IP1 and IP5
# bunch 3 is the bunch seen by the main bunch in IP2
# bunch 4 is the bunch seen by the main bunch in IP8
startingconditions                   = OrderedDict()
startingconditions['ex1']            = 3.5e-6/ringparams['gamma1'] # geometric emittance
startingconditions['ey1']            = 3.5e-6/ringparams['gamma1']
startingconditions['npart1']         = 1.0e11 # number of real particles in the bunch
startingconditions['blen1']          = 0.075 # bunch length
startingconditions['tauhat1']        = 1.25e-9  # half of bucket length

startingconditions['ex2']            = 1.5e-6/ringparams['gamma2']
startingconditions['ey2']            = 1.5e-6/ringparams['gamma2']
startingconditions['npart2']         = 1.0e8
startingconditions['blen2']          = 0.075
startingconditions['tauhat2']        = 1.25e-9

startingconditions['blen3']  = startingconditions['blen2'] # setting bunch length to 0.0 will make the simulation skip this bunch
startingconditions['ex3']    = startingconditions['ex2']
startingconditions['ey3']    = startingconditions['ey2']
startingconditions['npart3'] = int(startingconditions['npart2'])

startingconditions['blen4']  = startingconditions['blen2']
startingconditions['ex4']    = startingconditions['ex2']
startingconditions['ey4']    = startingconditions['ey2']
startingconditions['npart4'] = int(startingconditions['npart2'])

collisiondict                      = OrderedDict()
collisiondict['ips_leveling']      = {'IP1':False,'IP2':False,'IP5':False,'IP8':False} # leveling at IP1,IP2,IP5,IP8
collisiondict['ips_leveling_values'] = {'IP2':2.27e24}

# constants
protonmass                      = const.physical_constants['proton mass energy equivalent in MeV'][0]/1000 # GeV
electroncharge                  = const.physical_constants['elementary charge'][0]

# set physical processes to use in simulation
switches                        = OrderedDict()
switches['RFswitch']            = True    # RFswitch (set to 1 to activate synchrotron motion)
switches['betatronSwitch']      = True    # betatronSwitch (set to 1 to activate betatron motion)
switches['raddampSwitch']       = True    # raddampSwitch (set to 1 to activate radiation damping and quantum excitation)
switches['IBSswitch']           = True    # IBSswitch (set to 1 to activate ibs)
switches['collimationSwitch']   = False    # collimationSwitch (set to 1 to activate losses on aperture cuts)
switches['blowupSwitch']        = False    # blowupSwitch (set to 1 to activate artificial blowup - ADT loss maps etc)
switches['collisionSwitch']     = True    # collisionSwitch (set to 1 to activate collisions)
switches['writeAllCoordSwitch'] = False    # writeAllCoordSwitch
switches['writeMountSwitch']    = False    # writeMountSwitch

# set general parameters
ringparams['proton']            = protonmass
ringparams['ion']               = 193.7291748489224
ringparams['thib']              = 2.5e-9      # RF period
ringparams['dqmin1']            = 0.0     # (coupling between x-py, y-px)
ringparams['dqmin2']            = 0.0      # (coupling between x-py, y-px)
ringparams['k2L1']              = 1.1       # (thin sextupole strengths)
ringparams['k2Lskew1']          = 0.0        # (thin sextupole strengths)
ringparams['k2L2']              = 1.1       # (thin sextupole strengths)
ringparams['k2Lskew2']          = 0.0        # (thin sextupole strengths)

# setting the methods to generate the initial particle distributions
startingconditions['longcoor']       = 3                           # (0:parabolic with smokering, 1: read from file, 2: bi-Gaussian, 3: pseudo-Gaussian,exactly matched)
startingconditions['transcoor']      = 2                           # transvCoordMethod (1 or 2)
startingconditions['bunchLenPrecis'] = 0.01                        # (for longCoordMethod=3)
startingconditions['power']          = 0.75                        # (for longCoordMethod=0)
startingconditions['alint']          = 5                           # (for longCoordMethod=0)
startingconditions['coordfn1']       = 'startCoord.dat'
startingconditions['coordfn2']       = 'startCoord.dat'

# settings for radiation damping process
radiationdamping                     = OrderedDict()
radiationdamping['radMethod']        = 'lattic'          # adMethod: can be manual (input next line), approx (smooth latt. I4=0), or lattic (rad. int., twiss file required)
radiationdamping['tradlong']         = 23519.0           # in sec
radiationdamping['tradperp']         = 47072.0           # in sec
radiationdamping['siglong']          = 4.64e-12          # (eq. sigma from raddamp-exit.)
radiationdamping['sigperp']          = 7.43e-8           # (m) (only used with manual method)
radiationdamping['rho0']             = 2784.32           # (dipole bend. radius in m, used only with approx)

# settings for IBS calculations
ibsdict                              = OrderedDict()
ibsdict['ibsMethod']                 = 'nagaitsev'       # can be nagaitsev, piwiSmooth, piwLattice, modPiwLatt, baneApprox or interpolat
ibsdict['coupleIBS']                 = False                 # coupleIBS (0 gives separate growth rates, 1 gives same growth in x and y)
ibsdict['coulombLog']                = 20.0              # (used in nagaitsev method)
ibsdict['fracibstot']                = 1                 # factor of ibs strength
ibsdict['nbins']                     = 500               # number of bins
ibsdict['piwinski']                  = '/home/roderik/My_CERN_work/piwinski-ibs/ibs-rates-LHC-collision.dat'            
ibsdict['bane']                      = '/home/roderik/My_CERN_work/bane-ibs/gBaneTab.dat'              
ibsdict['intfile']                   = 'interpolate.dat'

# settings for collimation routine
colldict                     = OrderedDict()
colldict['refEmxy']          = 5.06158e-10
colldict['cutoffAmpl']       = 12.0 #5.7           # (ref. sigma at which init. distr. is cut)
colldict['collimAvgSwitch']  = 0            # 
colldict['emitMethod']       = 'stdev'       # (stdev,stcal,exfit)
colldict['nSigCutBeta']      = 12.0   #5.7           # (ref. sigma at which init. distr. is cut)
colldict['nSigCutMom']       = 10.0          # (momentum cut in beta-sigma)
colldict['betaxMom']         = 131.7        
colldict['dispxMom']         = 2.15         

# settings for blowup routine
blowupdict                 = OrderedDict()
blowupdict['pxKickFac']    = 0.02 
blowupdict['pyKickFac']    = 0.0
blowupdict['blowupMethod'] = 'unSum'         # (unifo,gauss,unSum)

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

collisiondict['longIntBins']       = 100


if (ringparams['aatomb'+str(1)] == 208):
        v001 = ringparams['ion']/  ringparams['qatomb'+str(1)] * 1e9
else:
        v001 = ringparams['proton']/  ringparams['qatomb'+str(1)] *1e9
        
if (ringparams['aatomb'+str(2)] == 208):
        v002 = ringparams['ion']/  ringparams['qatomb'+str(2)] * 1e9
else:
        v002 = ringparams['proton']/  ringparams['qatomb'+str(2)] *1e9

st = time.clock()
print 'start at:',time.ctime()
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

# function to get energyspread from bunchlength and RF voltage
def getsigpfromvrfandsigs(sigs,trev,h,vrf,eta,beta,eb):
    return 2*np.pi*sigs/trev*np.sqrt(h*vrf*eta/(2*np.pi*beta*eb))/eta/const.c

startingconditions['rmsDelta1'] = getsigpfromvrfandsigs(startingconditions['blen1'],trev1,ringparams['nharm2b1'],-1*ringparams['vrf2b1'],eta1,betar1,energy)
startingconditions['rmsDelta2'] = getsigpfromvrfandsigs(startingconditions['blen2'],trev2,ringparams['nharm2b2'],-1*ringparams['vrf2b2'],eta2,betar2,energy)
startingconditions['rmsDelta3'] = getsigpfromvrfandsigs(startingconditions['blen3'],trev2,ringparams['nharm2b2'],-1*ringparams['vrf2b2'],eta2,betar2,energy)
startingconditions['rmsDelta4'] = getsigpfromvrfandsigs(startingconditions['blen4'],trev2,ringparams['nharm2b2'],-1*ringparams['vrf2b2'],eta2,betar2,energy)

print startingconditions['blen1'],trev1,ringparams['nharm2b1'],-1*ringparams['vrf2b1'],eta1,betar1,energy
print startingconditions['rmsDelta1'],startingconditions['rmsDelta2'] 

betax1,betay1,betax2,betay2,pcoeff1,tcoeff1,rcoeff1,pcoeff2,tcoeff2,rcoeff2,dgamma_hat1,dgamma_hat2,tunes1,tunes2,tunesexact1,tunesexact2,dftfs1,dftfs2=setparam.set_parameters2(ringparams,startingconditions,omegarf1,trev1,eta1,betar1,omegarf2,trev2,eta2,betar2)

# reading bane or interpolat files for ibs if necessary
if (switches['IBSswitch']) & (ibsdict['ibsMethod']=='baneApprox'):
        dfbane = pd.read_csv(ibsdict['bane'],delim_whitespace=True)

if (switches['IBSswitch']) & (ibsdict['ibsMethod']=='interpolat'):
        dfinterpolat = pd.read_csv(ibsdict['int'],delim_whitespace=True)

# raddamping
if (switches['raddampSwitch'])&(radiationdamping['radMethod']!='manual'):
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


# initialisation of the variables to track evolutio
dflumip1 = pd.DataFrame(columns= ['sim.turn','t(hours)', 'L(cm^-2 s^-1)','red factor','betas'])
dflumip2 = pd.DataFrame(columns= ['sim.turn','t(hours)', 'L(cm^-2 s^-1)','red factor','betas'])
dflumip5 = pd.DataFrame(columns= ['sim.turn','t(hours)', 'L(cm^-2 s^-1)','red factor','betas'])
dflumip8 = pd.DataFrame(columns= ['sim.turn','t(hours)', 'L(cm^-2 s^-1)','red factor','betas'])

dfint1 = pd.DataFrame(columns=['sim.turn','t(hours)','N1_macro','N1_real','NlostLum1','SumL1','NlostDebunch1',
                               'SumD1','NLostBet1','SumB1','NlostMom1','SumM1'])
dfint2 = pd.DataFrame(columns=['sim.turn','t(hours)','N1_macro','N1_real','NlostLum1','SumL1','NlostDebunch1',
                               'SumD1','NLostBet1','SumB1','NlostMom1','SumM1'])
dfint3 = pd.DataFrame(columns=['sim.turn','t(hours)','N1_macro','N1_real','NlostLum1','SumL1','NlostDebunch1',
                               'SumD1','NLostBet1','SumB1','NlostMom1','SumM1'])
dfint4 = pd.DataFrame(columns=['sim.turn','t(hours)','N1_macro','N1_real','NlostLum1','SumL1','NlostDebunch1',
                               'SumD1','NLostBet1','SumB1','NlostMom1','SumM1'])

dfemit1 = pd.DataFrame(columns= ['sim.turn','t(hours)', 'ex1(m)','ey1(m)','el1(eV/s/charge)','sig1_T1','sig1_dP/P_2'])
dfemit2 = pd.DataFrame(columns= ['sim.turn','t(hours)', 'ex1(m)','ey1(m)','el1(eV/s/charge)','sig1_T1','sig1_dP/P_2'])
dfemit3 = pd.DataFrame(columns= ['sim.turn','t(hours)', 'ex1(m)','ey1(m)','el1(eV/s/charge)','sig1_T1','sig1_dP/P_2'])
dfemit4 = pd.DataFrame(columns= ['sim.turn','t(hours)', 'ex1(m)','ey1(m)','el1(eV/s/charge)','sig1_T1','sig1_dP/P_2'])

dfibs1  = pd.DataFrame(columns=['sim.turn','t(hours)','Tp(hours)','Tx(hours)','Ty(hours)'])
dfibs2  = pd.DataFrame(columns=['sim.turn','t(hours)','Tp(hours)','Tx(hours)','Ty(hours)'])
dfibs3  = pd.DataFrame(columns=['sim.turn','t(hours)','Tp(hours)','Tx(hours)','Ty(hours)'])
dfibs4  = pd.DataFrame(columns=['sim.turn','t(hours)','Tp(hours)','Tx(hours)','Ty(hours)'])

dfcolimpxbeta1 = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])
dfcolimpybeta1 = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])
dfcolimpxmom1  = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])
dfcolimpymom1  = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])

dfcolimpxbeta2 = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])
dfcolimpybeta2 = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])
dfcolimpxmom2  = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])
dfcolimpymom2  = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])

dfcolimpxbeta3 = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])
dfcolimpybeta3 = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])
dfcolimpxmom3  = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])
dfcolimpymom3  = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])

dfcolimpxbeta4 = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])
dfcolimpybeta4 = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])
dfcolimpxmom4  = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])
dfcolimpymom4  = pd.DataFrame(columns=['sim.turn','t(hours)','impact_par(m)','impact_par(sigma)'])


# adding particles - generating the distributions
# using the  startingconditions['blen1'] as test if this bunch needs to be simulated
if startingconditions['blen1']!= 0.0:
    bunch1 = True
else:
    print 'Bunch in beam 1 needed for simulation'
    quit()
    
if startingconditions['blen2']!= 0.0:
    bunch2 = True
else:
    bunch2 = False
    
if startingconditions['blen3']!= 0.0:
    bunch3 = True
else:
    bunch3 = False
    
if startingconditions['blen4']!= 0.0:
    bunch4 = True
else:
    bunch4 = False
    
t1,pt1 = faddpart.addlongitudinal(startingconditions['tauhat1'],
                                tcoeff1,2*np.pi/trev1,v001,fnharmb1,fnharm2b1,ringparams['vrfb1'],ringparams['vrf2b1'],
                                np.pi,5,.75,12763,ringparams['gamma1'],
                                startingconditions['rmsDelta1'],startingconditions['blen1'],startingconditions['bunchLenPrecis'],
                                const.c,ringparams['nMacro'] ,3)
x1,px1,y1,py1 = faddpart.addtransverse(ringparams['nMacro'],colldict['cutoffAmpl'],colldict['cutoffAmpl'],betax1,betay1,
                                   startingconditions['ex1'],startingconditions['ey1'],12673)


if bunch2:
    t2,pt2 = faddpart.addlongitudinal(startingconditions['tauhat2'],
                                tcoeff2,2*np.pi/trev2,v002,fnharmb2,fnharm2b2,ringparams['vrfb2'],ringparams['vrf2b2'],
                                np.pi,5,.75,12763,ringparams['gamma2'],
                                startingconditions['rmsDelta2'],startingconditions['blen2'],startingconditions['bunchLenPrecis'],
                                const.c,ringparams['nMacro'] ,3)
    x2,px2,y2,py2 = faddpart.addtransverse(ringparams['nMacro'],colldict['cutoffAmpl'],colldict['cutoffAmpl'],betax2,betay2,
                                   startingconditions['ex2'],startingconditions['ey2'],12673)


if bunch3:
    t3,pt3 = faddpart.addlongitudinal(startingconditions['tauhat2'],
                                tcoeff2,2*np.pi/trev2,v002,fnharmb2,fnharm2b2,ringparams['vrfb2'],ringparams['vrf2b2'],
                                np.pi,5,.75,12763,ringparams['gamma2'],
                                startingconditions['rmsDelta3'],startingconditions['blen3'],startingconditions['bunchLenPrecis'],
                                const.c,ringparams['nMacro'] ,3)
    x3,px3,y3,py3 = faddpart.addtransverse(ringparams['nMacro'],colldict['cutoffAmpl'],colldict['cutoffAmpl'],betax2,betay2,
                                   startingconditions['ex3'],startingconditions['ey3'],12673)


if bunch4:
    t4,pt4 = faddpart.addlongitudinal(startingconditions['tauhat2'],
                                tcoeff2,2*np.pi/trev2,v002,fnharmb2,fnharm2b2,ringparams['vrfb2'],ringparams['vrf2b2'],
                                np.pi,5,.75,12763,ringparams['gamma2'],
                                startingconditions['rmsDelta4'],startingconditions['blen4'],startingconditions['bunchLenPrecis'],
                                const.c,ringparams['nMacro'] ,3)
    x4,px4,y4,py4 = faddpart.addtransverse(ringparams['nMacro'],colldict['cutoffAmpl'],colldict['cutoffAmpl'],betax2,betay2,
                                   startingconditions['ex4'],startingconditions['ey4'],12673)


avgline1        = np.zeros(ibsdict['nbins'])
avglinex1       = np.zeros(ibsdict['nbins'])
avgliney1       = np.zeros(ibsdict['nbins'])
avgline2        = np.zeros(ibsdict['nbins'])
avglinex2       = np.zeros(ibsdict['nbins'])
avgliney2       = np.zeros(ibsdict['nbins'])
avgline3        = np.zeros(ibsdict['nbins'])
avglinex3       = np.zeros(ibsdict['nbins'])
avgliney3       = np.zeros(ibsdict['nbins'])
avgline4        = np.zeros(ibsdict['nbins'])
avglinex4       = np.zeros(ibsdict['nbins'])
avgliney4       = np.zeros(ibsdict['nbins'])
nLostLum1       =0
nLostLumSum1    =0
nLostLum2       =0
nLostLumSum2    =0
nLostDebunch1   =0
nLostDebunch2   =0
nLostDebunchSum1=0
nLostDebunchSum2=0
nLostMom1       =0
nLostMomSum1    =0
nLostBeta1      =0
nLostBetaSum1   =0
nLostMom2       =0
nLostMomSum2    =0
nLostBeta2      =0
nLostBetaSum2   =0
nLostLum3       =0
nLostLumSum3    =0
nLostDebunch3   =0
nLostDebunchSum3=0
nLostMom3       =0
nLostMomSum3    =0
nLostBeta3      =0
nLostBetaSum3   =0
nLostLum4       =0
nLostLumSum4    =0
nLostDebunch4   =0
nLostDebunchSum4=0
nLostMom4       =0
nLostMomSum4    =0
nLostBeta4      =0
nLostBetaSum4   =0
nLostLum1_ip1 = 0
nLostLum2_ip1 = 0
nLostLum1_ip2 = 0
nLostLum3_ip2 = 0
nLostLum1_ip5 = 0
nLostLum2_ip5 = 0
nLostLum1_ip8 = 0
nLostLum4_ip8 = 0

dfparticles1 = pd.DataFrame(columns=['x','px','y','py','t','pt'])
dfparticles2 = pd.DataFrame(columns=['x','px','y','py','t','pt'])
dfparticles3 = pd.DataFrame(columns=['x','px','y','py','t','pt'])
dfparticles4 = pd.DataFrame(columns=['x','px','y','py','t','pt'])

if bunch1:
    dfparticles1['x'] = x1
    dfparticles1['px']= px1
    dfparticles1['y'] = y1
    dfparticles1['py']= py1
    dfparticles1['t'] = t1
    dfparticles1['pt']= pt1
    pnumber1 = startingconditions['npart1']
else:
    pnubmer1 = 0
if bunch2:
    dfparticles2['x'] = x2
    dfparticles2['px']= px2
    dfparticles2['y'] = y2
    dfparticles2['py']= py2
    dfparticles2['t'] = t2
    dfparticles2['pt']= pt2
    pnumber2 = startingconditions['npart2']
else:
    pnumber2 = 0
if bunch3:
    dfparticles3['x'] = x3
    dfparticles3['px']= px3
    dfparticles3['y'] = y3
    dfparticles3['py']= py3
    dfparticles3['t'] = t3
    dfparticles3['pt']= pt3
    pnumber3 = startingconditions['npart3']
else:
    pnumber3 = 0
if bunch4:
    dfparticles4['x'] = x4
    dfparticles4['px']= px4
    dfparticles4['y'] = y4
    dfparticles4['py']= py4
    dfparticles4['t'] = t4
    dfparticles4['pt']= pt4
    pnumber4 = startingconditions['npart4']
else:
    pnumber4 = 0
    
n01 = ringparams['nMacro']
n02 = ringparams['nMacro']
n03 = ringparams['nMacro']
n04 = ringparams['nMacro']

if collisiondict['ips_leveling']['IP1']:
    lumimax1 = collisiondict['ips_leveling_values']['IP1']
else: 
    lumimax1 = 0.0
if collisiondict['ips_leveling']['IP2']:
    lumimax2 = collisiondict['ips_leveling_values']['IP2']
else: 
    lumimax2 = 0.0
if collisiondict['ips_leveling']['IP5']:
    lumimax5 = collisiondict['ips_leveling_values']['IP5']
else: 
    lumimax5 = 0.0
if collisiondict['ips_leveling']['IP8']:
    lumimax8 = collisiondict['ips_leveling_values']['IP8']
else: 
    lumimax8 = 0.0
print 'took :', time.clock()-st

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

# bane ibs
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

print 'start at:',time.ctime()
print'Starting main loop ...'

ex1, ey1 = fgetemit.getemittance(x1,y1,px1,py1,betax1,betay1,xcut1,ycut1,colldict['cutoffAmpl'],colldict['emitMethod'],colldict['refEmxy'])

if bunch2:
    ex2, ey2 = fgetemit.getemittance(x2,y2,px2,py2,betax2,betay2,xcut2,ycut2,colldict['cutoffAmpl'],colldict['emitMethod'],colldict['refEmxy'])
else:
    ex2 = 0.0
    ey2 = 0.0
if bunch3:
    ex3, ey3 = fgetemit.getemittance(x3,y3,px3,py3,betax2,betay2,xcut2,ycut2,colldict['cutoffAmpl'],colldict['emitMethod'],colldict['refEmxy'])
else:
    ex3 = 0.0
    ey3 = 0.0
if bunch4:
    ex4, ey4 = fgetemit.getemittance(x4,y4,px4,py4,betax2,betay2,xcut2,ycut2,colldict['cutoffAmpl'],colldict['emitMethod'],colldict['refEmxy'])
else:
    ex4 = 0.0
    ey4 = 0.0
print 'ex1={:.6e}, ex2={:.6e}, ex3={:.6e}, ex4={:.6e}'.format(ex1,ex2,ex3,ex4)
print 'ey1={:.6e}, ey2={:.6e}, ey3={:.6e}, ey4={:.6e}'.format(ey1,ey2,ey3,ey4)

#  write first turn
kturns=0
iwrite=1
eqtime  = ringparams['nturns'] * trev1 * ringparams['timeRatio']/3600.
if (switches['writeMountSwitch']):
    dfwritemount1 = pd.DataFrame(columns = cwrite.mountainnrcols)
    dfwritemount2 = pd.DataFrame(columns = cwrite.mountainnrcols)
    dfwritemount3 = pd.DataFrame(columns = cwrite.mountainnrcols)
    dfwritemount4 = pd.DataFrame(columns = cwrite.mountainnrcols)
    dfwritemount1 = cwrite.mountainr(ringparams['nMacro'],startingconditions['npart1'],avgline1,ringparams['qatomb1'],
                                         ringparams['nwrite'],kturns,ringparams['thib'],ibsdict['nbins'],dfwritemount1)
    if bunch2:
        dfwritemount2 = cwrite.mountainr(ringparams['nMacro'],startingconditions['npart2'],avgline2,ringparams['qatomb2'],
                                         ringparams['nwrite'],kturns,ringparams['thib'],ibsdict['nbins'],dfwritemount2)
    if bunch3:
        dfwritemount3 = cwrite.mountainr(ringparams['nMacro'],startingconditions['npart3'],avgline3,ringparams['qatomb2'],
                                         ringparams['nwrite'],kturns,ringparams['thib'],ibsdict['nbins'],dfwritemount3)
    if bunch4:
        dfwritemount4 = cwrite.mountainr(ringparams['nMacro'],startingconditions['npart4'],avgline4,ringparams['qatomb2'],
                                         ringparams['nwrite'],kturns,ringparams['thib'],ibsdict['nbins'],dfwritemount3)

    dfwritemounttrans1 = pd.DataFrame(columns = cwrite.mounttransvcols)
    dfwritemounttrans2 = pd.DataFrame(columns = cwrite.mounttransvcols)
    dfwritemounttrans3 = pd.DataFrame(columns = cwrite.mounttransvcols)
    dfwritemounttrans4 = pd.DataFrame(columns = cwrite.mounttransvcols)

    dfwritemounttrans1 = cwrite.mountTransv(avglinex1,avgliney1,ex1,ey1,betax1,betay1,ibsdict['nbins'],
                                                kturns,dfwritemounttrans1)
    if bunch2:
        dfwritemounttrans2 = cwrite.mountTransv(avglinex2,avgliney2,ex2,ey2,betax2,betay2,ibsdict['nbins'],
                                                kturns,dfwritemounttrans2)
    if bunch3:
        dfwritemounttrans3 = cwrite.mountTransv(avglinex3,avgliney3,ex3,ey3,betax2,betay2,ibsdict['nbins'],
                                                kturns,dfwritemounttrans3)
    if bunch4:
        dfwritemounttrans4 = cwrite.mountTransv(avglinex4,avgliney4,ex4,ey4,betax2,betay2,ibsdict['nbins'],
                                                kturns,dfwritemounttrans4)
    


dfemit1 = cwrite.writeEmi(ex1,ey1,t1,pt1,betar1,ringparams['gamma1'],
                              ringparams['aatomb1'],ringparams['qatomb1'],eqtime,kturns,ringparams['nturns'],dfemit1)
dfint1,nLostLumSum1,nLostDebunchSum1,nLostBetaSum1,nLostMomSum1 = cwrite.writeNb(dfparticles1,ringparams['nMacro'],startingconditions['npart1'], nLostLum1,
                       nLostLumSum1,nLostDebunch1,nLostDebunchSum1,nLostBeta1,nLostBetaSum1,nLostMom1,nLostMomSum1,
                       kturns,ringparams['nturns'],eqtime,dfint1)
    
if bunch2:    
    dfemit2 = cwrite.writeEmi(ex2,ey2,t2,pt2,betar2,ringparams['gamma2'],
                              ringparams['aatomb2'],ringparams['qatomb2'],eqtime,kturns,ringparams['nturns'],dfemit2)
    dfint2,nLostLumSum2,nLostDebunchSum2,nLostBetaSum2,nLostMomSum2 =cwrite.writeNb(dfparticles2,ringparams['nMacro'],startingconditions['npart2'], nLostLum2,
                      nLostLumSum2,nLostDebunch2,nLostDebunchSum2,nLostBeta2,nLostBetaSum2,nLostMom2,nLostMomSum2,
                       kturns,ringparams['nturns'],eqtime,dfint2)
    
if bunch3:
    dfemit3 = cwrite.writeEmi(ex3,ey3,t3,pt3,betar2,ringparams['gamma2'],
                              ringparams['aatomb2'],ringparams['qatomb2'],eqtime,kturns,ringparams['nturns'],dfemit3)
    dfint3,nLostLumSum3,nLostDebunchSum3,nLostBetaSum3,nLostMomSum3 = cwrite.writeNb(dfparticles3,ringparams['nMacro'],startingconditions['npart3'], nLostLum3,
                       nLostLumSum3,nLostDebunch3,nLostDebunchSum3,nLostBeta3,nLostBetaSum3,nLostMom3,nLostMomSum3,
                       kturns,ringparams['nturns'],eqtime,dfint3)

if bunch4:
    dfemit4 = cwrite.writeEmi(ex4,ey4,t4,pt4,betar2,ringparams['gamma2'],
                              ringparams['aatomb2'],ringparams['qatomb2'],eqtime,kturns,ringparams['nturns'],dfemit4)
    dfint4,nLostLumSum4,nLostDebunchSum4,nLostBetaSum4,nLostMomSum4 =cwrite.writeNb(dfparticles4,ringparams['nMacro'],startingconditions['npart4'], nLostLum4,
                       nLostLumSum4,nLostDebunch4,nLostDebunchSum4,nLostBeta4,nLostBetaSum4,nLostMom4,nLostMomSum4,
                       kturns,ringparams['nturns'],eqtime,dfint4)

if (switches['writeAllCoordSwitch']==1):   #     write all starting conditions to file
    dfparticles1.to_csv('all-coord-b1-'+str(kturn)+'.csv',index=False)
    if bunch2:
        dfparticles2.to_csv('all-coord-b15-'+str(kturn)+'.csv',index=False)
    if bunch3:
        dfparticles3.to_csv('all-coord-b2-'+str(kturn)+'.csv',index=False)
    if bunch4:
        dfparticles4.to_csv('all-coord-b8-'+str(kturn)+'.csv',index=False)

# calculation of parameters
psix1,psiy1,coeffchromx1,coeffchromy1,tpdqmin1 = setparam.settunechroma(1,ringparams)
psix2,psiy2,coeffchromx2,coeffchromy2,tpdqmin2 = setparam.settunechroma(2,ringparams)

nwrite = ringparams['nwrite']
nturns = ringparams['nturns']

print 'start at:',time.ctime(),' took: ',time.clock()-st

fmix1= 1
for kturns in  range(1,ringparams['nturns']+1):
    kcheck = kturns/nwrite
    kcheck = kcheck*nwrite
    
    if(kcheck==kturns):
        iwrite = 1
        print 'took :', time.clock()-st
    else:
        iwrite = 0
    
    if (kturns==nturns+1):
        print ' writing last turn'
        nwrite=-1
        iwrite=1
   
    ex1, ey1 = fgetemit.getemittance(x1,y1,px1,py1,betax1,betay1,xcut1,ycut1,colldict['cutoffAmpl'],colldict['emitMethod'],colldict['refEmxy'])
    if bunch2:
        ex2, ey2 = fgetemit.getemittance(x2,y2,px2,py2,betax2,betay2,xcut2,ycut2,colldict['cutoffAmpl'],colldict['emitMethod'],colldict['refEmxy'])
    if bunch3:
        ex3, ey3 = fgetemit.getemittance(x3,y3,px3,py3,betax2,betay2,xcut2,ycut2,colldict['cutoffAmpl'],colldict['emitMethod'],colldict['refEmxy'])
    if bunch4:
        ex4, ey4 = fgetemit.getemittance(x4,y4,px4,py4,betax2,betay2,xcut2,ycut2,colldict['cutoffAmpl'],colldict['emitMethod'],colldict['refEmxy'])


#   synchrotron motion
    if(switches['RFswitch']):
        if bunch1:
            np1,x1,px1,y1,py1,t1,pt1,nLostDebunch1,kept1 = frfupdate.rfupdate(x1,px1,y1,py1,t1,pt1,fmix1,nLostDebunch1,2*np.pi/trev1,v001,ringparams['thib'],
                                                            tcoeff1,fnharmb1,fnharm2b1,ringparams['vrfb1'],ringparams['vrf2b1'])
            x1  = x1[:kept1]
            px1 = px1[:kept1]
            y1  = y1[:kept1]
            py1 = py1[:kept1]
            t1  = t1[:kept1]
            pt1 = pt1[:kept1]
            if (iwrite==1):
                print kept1,nLostDebunch1 
        if bunch2:
            np2,x2,px2,y2,py2,t2,pt2,nLostDebunch2,kept2 = frfupdate.rfupdate(x2,px2,y2,py2,t2,pt2,fmix1,nLostDebunch2,2*np.pi/trev2,v002,ringparams['thib'],
                                                            tcoeff2,fnharmb2,fnharm2b2,ringparams['vrfb2'],ringparams['vrf2b2'])
            x2  = x2[:kept2]
            px2 = px2[:kept2]
            y2  = y2[:kept2]
            py2 = py2[:kept2]
            t2  = t2[:kept2]
            pt2 = pt2[:kept2]
        if bunch3:
            np3,x3,px3,y3,py3,t3,pt3,nLostDebunch3,kept3 = frfupdate.rfupdate(x3,px3,y3,py3,t3,pt3,fmix1,nLostDebunch3,2*np.pi/trev2,v002,ringparams['thib'],
                                                            tcoeff2,fnharmb2,fnharm2b2,ringparams['vrfb2'],ringparams['vrf2b2'])
            x3  = x3[:kept3]
            px3 = px3[:kept3]
            y3  = y3[:kept3]
            py3 = py3[:kept3]
            t3  = t3[:kept3]
            pt3 = pt3[:kept3]
        if bunch4:
            np4,x4,px4,y4,py4,t4,pt4,nLostDebunch4,kept4 = frfupdate.rfupdate(x4,px4,y4,py4,t4,pt4,fmix1,nLostDebunch4,2*np.pi/trev2,v002,ringparams['thib'],
                                                            tcoeff2,fnharmb2,fnharm2b2,ringparams['vrfb2'],ringparams['vrf2b2'])
            x4  = x4[:kept4]
            px4 = px4[:kept4]
            y4  = y4[:kept4]
            py4 = py4[:kept4]
            t4  = t4[:kept4]
            pt4 = pt4[:kept4]
            
#      betatron motion
    if(switches['betatronSwitch']):
        x1,px1,y1,py1,pt1   = fsptranschrom.sptranschrom(x1,px1,y1,py1,pt1,psix1,psiy1,coeffchromx1,coeffchromy1,tpdqmin1,ringparams['k2L1'],ringparams['k2Lskew1'])
        avglinex1,avgliney1 = fkeeptransprof.keeptransvprof(avglinex1,avgliney1,x1,y1,ex1,ey1,betax1,betay1)
            
        if bunch2:
            x2,px2,y2,py2,pt2   = fsptranschrom.sptranschrom(x2,px2,y2,py2,pt2,psix2,psiy2,coeffchromx2,coeffchromy2,tpdqmin2,ringparams['k2L2'],ringparams['k2Lskew2'])
            avglinex2,avgliney2 = fkeeptransprof.keeptransvprof(avglinex2,avgliney2,x2,y2,ex2,ey2,betax2,betay2)
            
        if bunch3:
            x3,px3,y3,py3,pt3   = fsptranschrom.sptranschrom(x3,px3,y3,py3,pt3,psix2,psiy2,coeffchromx2,coeffchromy2,tpdqmin2,ringparams['k2L2'],ringparams['k2Lskew2'])
            avglinex3,avgliney3 = fkeeptransprof.keeptransvprof(avglinex3,avgliney3,x3,y3,ex3,ey3,betax2,betay2)
            
        if bunch4:
            x4,px4,y4,py4,pt4   = fsptranschrom.sptranschrom(x4,px4,y4,py4,pt4,psix2,psiy2,coeffchromx2,coeffchromy2,tpdqmin2,ringparams['k2L2'],ringparams['k2Lskew2'])
            avglinex4,avgliney4 = fkeeptransprof.keeptransvprof(avglinex4,avgliney4,x4,y4,ex4,ey4,betax2,betay2)
            
#     radiation damping and quantum excitation
    if(switches['raddampSwitch']):
        x1,px1,y1,py1,pt1,iseed,np1 = frad.raddamp(x1,px1,y1,py1,pt1,tradlong1,tradperp1,trev1,siglong1,sigperp1,ringparams['timeRatio'],12763)
        if bunch2:
            x2,px2,y2,py2,pt2,iseed,np2 = frad.raddamp(x2,px2,y2,py2,pt2,tradlong2,tradperp2,trev2,siglong2,sigperp2,ringparams['timeRatio'],12763)
        if bunch3:
            x3,px3,y3,py3,pt3,iseed,np3 = frad.raddamp(x3,px3,y3,py3,pt3,tradlong2,tradperp2,trev2,siglong2,sigperp2,ringparams['timeRatio'],12763)
        if bunch4:
            x4,px4,y4,py4,pt4,iseed,np4 = frad.raddamp(x4,px4,y4,py4,pt4,tradlong2,tradperp2,trev2,siglong2,sigperp2,ringparams['timeRatio'],12763)  
# ibs  
    if (switches['IBSswitch']):
        px1,py1,t1,pt1,ibsrow1,iseed = fibs.ibslong(pnumber1,px1,py1,t1,pt1,ex1,ey1,ringparams['nMacro'],ringparams['circ'],
                                                        iwrite,ringparams['thib'],avgline1,betar1,ringparams['gamma1'],
                                                        ringparams['gammat'],betax1,betay1,ringparams['aatomb1'],ringparams['qatomb1'],
                                                        trev1,np.array(dftfs1['BETX'].values,order='F'),np.array(dftfs1['BETY'].values,order='F'),
                                                        np.array(dftfs1['DX'].values,order='F'),np.array(dftfs1['DY'].values,order='F'),
                                                        np.array(dftfs1['DPX'].values,order='F'),np.array(dftfs1['DPY'].values,order='F'),
                                                        np.array(dftfs1['ALFX'].values,order='F'),np.array(dftfs1['ALFY'].values,order='F'),
                                                        np.array(dftfs1['L'].values,order='F'),gtab,ibsdict['coulombLog'],xmin,xmax,ymin,ymax,zmin,zmax,ap,ax,ay,
                                                        pnuminterp,ibsdict['ibsMethod'],ibsdict['coupleIBS'],nturns,kturns,eqtime,
                                                        ibsdict['fracibstot'],const.c,np.pi,ringparams['timeRatio'] ,12673)
            
        dfibs1 = dfibs1.append(pd.DataFrame(ibsrow1.reshape((1,5)),columns=dfibs1.columns))
            
        if bunch2:
            px2,py2,t2,pt2,ibsrow2,iseed = fibs.ibslong(pnumber2,px2,py2,t2,pt2,ex2,ey2,ringparams['nMacro'],ringparams['circ'],
                                                        iwrite,ringparams['thib'],avgline2,betar2,ringparams['gamma2'],
                                                        ringparams['gammat'],betax2,betay2,ringparams['aatomb2'],ringparams['qatomb2'],
                                                        trev2,np.array(dftfs2['BETX'].values,order='F'),np.array(dftfs2['BETY'].values,order='F'),
                                                        np.array(dftfs2['DX'].values,order='F'),np.array(dftfs2['DY'].values,order='F'),
                                                        np.array(dftfs2['DPX'].values,order='F'),np.array(dftfs2['DPY'].values,order='F'),
                                                        np.array(dftfs2['ALFX'].values,order='F'),np.array(dftfs2['ALFY'].values,order='F'),
                                                        np.array(dftfs2['L'].values,order='F'),gtab,ibsdict['coulombLog'],xmin,xmax,ymin,ymax,zmin,zmax,ap,ax,ay,
                                                        pnuminterp,ibsdict['ibsMethod'],ibsdict['coupleIBS'],nturns,kturns,eqtime,
                                                        ibsdict['fracibstot'],const.c,np.pi,ringparams['timeRatio'] ,12673)
            dfibs2 = dfibs2.append(pd.DataFrame(ibsrow2.reshape((1,5)),columns=dfibs2.columns))
            
        if bunch3:
            px3,py3,t3,pt3,ibsrow3,iseed = fibs.ibslong(pnumber3,px3,py3,t3,pt3,ex3,ey3,ringparams['nMacro'],ringparams['circ'],
                                                        iwrite,ringparams['thib'],avgline3,betar2,ringparams['gamma2'],
                                                        ringparams['gammat'],betax2,betay2,ringparams['aatomb2'],ringparams['qatomb2'],
                                                        trev2,np.array(dftfs2['BETX'].values,order='F'),np.array(dftfs2['BETY'].values,order='F'),
                                                        np.array(dftfs2['DX'].values,order='F'),np.array(dftfs2['DY'].values,order='F'),
                                                        np.array(dftfs2['DPX'].values,order='F'),np.array(dftfs2['DPY'].values,order='F'),
                                                        np.array(dftfs2['ALFX'].values,order='F'),np.array(dftfs2['ALFY'].values,order='F'),
                                                        np.array(dftfs2['L'].values,order='F'),gtab,ibsdict['coulombLog'],xmin,xmax,ymin,ymax,zmin,zmax,ap,ax,ay,
                                                        pnuminterp,ibsdict['ibsMethod'],ibsdict['coupleIBS'],nturns,kturns,eqtime,
                                                        ibsdict['fracibstot'],const.c,np.pi,ringparams['timeRatio'] ,12673)
            dfibs3 = dfibs3.append(pd.DataFrame(ibsrow3.reshape((1,5)),columns=dfibs3.columns))
            
        if bunch4:
            px4,py4,t4,pt4,ibsrow4,iseed = fibs.ibslong(pnumber4,px4,py4,t4,pt4,ex4,ey4,ringparams['nMacro'],ringparams['circ'],
                                                        iwrite,ringparams['thib'],avgline4,betar2,ringparams['gamma2'],
                                                        ringparams['gammat'],betax2,betay2,ringparams['aatomb2'],ringparams['qatomb2'],
                                                        trev2,np.array(dftfs2['BETX'].values,order='F'),np.array(dftfs2['BETY'].values,order='F'),
                                                        np.array(dftfs2['DX'].values,order='F'),np.array(dftfs2['DY'].values,order='F'),
                                                        np.array(dftfs2['DPX'].values,order='F'),np.array(dftfs2['DPY'].values,order='F'),
                                                        np.array(dftfs2['ALFX'].values,order='F'),np.array(dftfs2['ALFY'].values,order='F'),
                                                        np.array(dftfs2['L'].values,order='F'),gtab,ibsdict['coulombLog'],xmin,xmax,ymin,ymax,zmin,zmax,ap,ax,ay,
                                                        pnuminterp,ibsdict['ibsMethod'],ibsdict['coupleIBS'],nturns,kturns,eqtime,
                                                        ibsdict['fracibstot'],const.c,np.pi,ringparams['timeRatio'] ,12673)
            dfibs4 = dfibs4.append(pd.DataFrame(ibsrow4.reshape((1,5)),columns=dfibs4.columns))
            
    if (switches['collimationSwitch'] ):
        if bunch1:
            x1,px1,y1,py1,t1,pt1,nLostMom1,nLostBeta1,collx1,colly1,momx1,kept1 = fcollim.collimation(x1,px1,y1,py1,t1,pt1,nLostMom1,nLostBeta1,
                    eqtime,betax1,betay1,colldict['refEmxy'],colldict['collimAvgSwitch'],colldict['nSigCutBeta'],
                    colldict['nSigCutMom'],colldict['betaxMom'],betar1,ringparams['gamma1'],colldict['dispxMom'],kturns,nturns)
            x1  = x1[:kept1]
            px1 = px1[:kept1]
            y1  = y1[:kept1]
            py1 = py1[:kept1]
            t1  = t1[:kept1]
            pt1 = pt1[:kept1]
            collx1 = collx1[~np.any(collx1==0, axis=1)]
            colly1 = colly1[~np.any(colly1==0, axis=1)]
            momx1  = momx1[~np.any(momx1==0, axis=1)]
            if (len(collx1) >0):
                dfcolimpxbeta1 = dfcolimpxbeta1.append(pd.DataFrame(collx1,columns=dfcolimpxbeta1.columns))
            if (len(colly1) >0):
                dfcolimpybeta1 = dfcolimpybeta1.append(pd.DataFrame(colly1,columns=dfcolimpybeta1.columns))
            if (len(momx1) >0):
                dfcolimpxmom1  = dfcolimpxmom1.append(pd.DataFrame(momx1,columns=dfcolimpxmom1.columns))
            
        if bunch2:
            x2,px2,y2,py2,t2,pt2,nLostMom2,nLostBeta2,collx2,colly2,momx2,kept2 = fcollim.collimation(x2,px2,y2,py2,t2,pt2,nLostMom2,nLostBeta2,
                    eqtime,betax2,betay2,colldict['refEmxy'],colldict['collimAvgSwitch'],colldict['nSigCutBeta'],
                    colldict['nSigCutMom'],colldict['betaxMom'],betar2,ringparams['gamma2'],colldict['dispxMom'],kturns,nturns)
            x2  = x2[:kept2]
            px2 = px2[:kept2]
            y2  = y2[:kept2]
            py2 = py2[:kept2]
            t2  = t2[:kept2]
            pt2 = pt2[:kept2]
            
            collx2 = collx2[~np.any(collx2==0, axis=1)]
            colly2 = colly2[~np.any(colly2==0, axis=1)]
            momx2  = momx2[~np.any(momx2==0, axis=1)]
            if (len(collx2)>0):
                dfcolimpxbeta2 = dfcolimpxbeta2.append(pd.DataFrame(collx2,columns=dfcolimpxbeta2.columns))
            if (len(colly2)>0):
                dfcolimpybeta2 = dfcolimpybeta2.append(pd.DataFrame(colly2,columns=dfcolimpybeta2.columns))
            if (len(momx2)>0):
                dfcolimpxmom2  = dfcolimpxmom2.append(pd.DataFrame(momx2,columns=dfcolimpxmom2.columns))
            
        if bunch3:
            x3,px3,y3,py3,t3,pt3,nLostMom3,nLostBeta3,collx3,colly3,momx3,kept3 = fcollim.collimation(x3,px3,y3,py3,t3,pt3,nLostMom3,nLostBeta3,
                    eqtime,betax2,betay2,colldict['refEmxy'],colldict['collimAvgSwitch'],colldict['nSigCutBeta'],
                    colldict['nSigCutMom'],colldict['betaxMom'],betar2,ringparams['gamma2'],colldict['dispxMom'],kturns,nturns)
            x3  = x3[:kept3]
            px3 = px3[:kept3]
            y3  = y3[:kept3]
            py3 = py3[:kept3]
            t3  = t3[:kept3]
            pt3 = pt3[:kept3]
            collx3 = collx3[~np.any(collx3==0, axis=1)]
            colly3 = colly3[~np.any(colly3==0, axis=1)]
            momx3  = momx3[~np.any(momx3==0, axis=1)]
            if (len(collx3)>0):
                dfcolimpxbeta3 = dfcolimpxbeta3.append(pd.DataFrame(collx3,columns=dfcolimpxbeta3.columns))
            if (len(colly3)>0):
                dfcolimpybeta3 = dfcolimpybeta3.append(pd.DataFrame(colly3,columns=dfcolimpybeta3.columns))
            if (len(momx3)>0):
                dfcolimpxmom3  = dfcolimpxmom3.append(pd.DataFrame(momx3,columns=dfcolimpxmom3.columns))
            
        if bunch4:
            x4,px4,y4,py4,t4,pt4,nLostMom4,nLostBeta4,collx4,colly4,momx4,kept4 = fcollim.collimation(x4,px4,y4,py4,t4,pt4,nLostMom4,nLostBeta4,
                    eqtime,betax2,betay2,colldict['refEmxy'],colldict['collimAvgSwitch'],colldict['nSigCutBeta'],
                    colldict['nSigCutMom'],colldict['betaxMom'],betar2,ringparams['gamma2'],colldict['dispxMom'],kturns,nturns)
            x4  = x4[:kept4]
            px4 = px4[:kept4]
            y4  = y4[:kept4]
            py4 = py4[:kept4]
            t4  = t4[:kept4]
            pt4 = pt4[:kept4]
            collx4 = collx4[~np.any(collx4==0, axis=1)]
            colly4 = colly4[~np.any(colly4==0, axis=1)]
            momx4  = momx4[~np.any(momx4==0, axis=1)]
            if (len(collx4)>0):
                dfcolimpxbeta4 = dfcolimpxbeta4.append(pd.DataFrame(collx4,columns=dfcolimpxbeta4.columns))
            if (len(colly4)>0):
                dfcolimpybeta4 = dfcolimpybeta4.append(pd.DataFrame(colly4,columns=dfcolimpybeta4.columns))
            if (len(momx4)>0):
                dfcolimpxmom4  = dfcolimpxmom4.append(pd.DataFrame(momx4,columns=dfcolimpxmom4.columns))
            
    if (switches['blowupSwitch']):
        px1,py1,iseed = fblowup.blowup(px1,py1,12673,blowupdict['blowupMethod'],betax1,betay1,colldict['refEmxy'],
                                blowupdict['pxKickFac'],blowupdict['pyKickFac'])
        if bunch2:
            px2,py2,iseed = fblowup.blowup(px2,py2,12673,blowupdict['blowupMethod'],betax2,betay2,colldict['refEmxy'],
                                blowupdict['pxKickFac'],blowupdict['pyKickFac'])
        if bunch3:
            px3,py3,iseed = fblowup.blowup(px3,py3,12673,blowupdict['blowupMethod'],betax2,betay2,colldict['refEmxy'],
                                blowupdict['pxKickFac'],blowupdict['pyKickFac'])
        if bunch4:
            px4,py4,iseed = fblowup.blowup(px4,py4,12673,blowupdict['blowupMethod'],betax2,betay2,colldict['refEmxy'],
                                blowupdict['pxKickFac'],blowupdict['pyKickFac'])
    if (switches['collisionSwitch']):
        if (collisiondict['collRoutine']=='6a'):
            # IP1
            if bunch2:
                hgold = hglassfacOld[0]
                hgfac,betass1,kept1,kept2,lumip1 = fcollision.collision6a(collisiondict['sigI'],ringparams['thib'],collisiondict['longIntBins'],
                                                                      ringparams['timeRatio'],ringparams['circ'],np.pi,const.c,12673,
                                                                      hgold,nLostLum1_ip1,nLostLum2_ip1,ringparams['nMacro'],ringparams['nMacro'],
                                                                      pnumber1,x1,px1,y1,py1,t1,pt1,ex1,ey1,betax1,betay1,bs1,lumimax1,theta1,
                                                                      pnumber2,x2,px2,y2,py2,t2,pt2,ex2,ey2,betax2,betay2)
                
                nLostLum1_ip1 = nLostLum1_ip1 + len(x1) - kept1
                nLostLum2_ip1 = nLostLum2_ip1 + len(x2) - kept2
                x1  = x1[:kept1]
                px1 = px1[:kept1]
                y1  = y1[:kept1]
                py1 = py1[:kept1]
                t1  = t1[:kept1]
                pt1 = pt1[:kept1]
                x2  = x2[:kept2]
                px2 = px2[:kept2]
                y2  = y2[:kept2]
                py2 = py2[:kept2]
                t2  = t2[:kept2]
                pt2 = pt2[:kept2]
                hglassfacOld[0] = hgfac
               
            # IP2
            if bunch3:
                hgold = hglassfacOld[1]
                hgfac,betass2,kept1,kept2,lumip2 = fcollision.collision6a(collisiondict['sigI'],ringparams['thib'],collisiondict['longIntBins'],
                                                                      ringparams['timeRatio'],ringparams['circ'],np.pi,const.c,12673,
                                                                      hgold,nLostLum1_ip2,nLostLum3_ip2,ringparams['nMacro'],ringparams['nMacro'],
                                                                      pnumber1,x1,px1,y1,py1,t1,pt1,ex1,ey1,betax1,betay1,bs2,lumimax2,theta2,
                                                                      pnumber3,x3,px3,y3,py3,t3,pt3,ex3,ey3,betax2,betay2)
                nLostLum1_ip2 = nLostLum1_ip2 + len(x1) - kept1
                nLostLum3_ip2 = nLostLum3_ip2 + len(x3) - kept2
                x1  = x1[:kept1]
                px1 = px1[:kept1]
                y1  = y1[:kept1]
                py1 = py1[:kept1]
                t1  = t1[:kept1]
                pt1 = pt1[:kept1]
                x3  = x3[:kept2]
                px3 = px3[:kept2]
                y3  = y3[:kept2]
                py3 = py3[:kept2]
                t3  = t3[:kept2]
                pt3 = pt3[:kept2]
                hglassfacOld[1] = hgfac
            # IP5 
            if bunch2:
                hgold = hglassfacOld[2]
                hgfac,betass5,kept1,kept2,lumip5 = fcollision.collision6a(collisiondict['sigI'],ringparams['thib'],collisiondict['longIntBins'],
                                                                    ringparams['timeRatio'],ringparams['circ'],np.pi,const.c,12673,
                                                                      hgold,nLostLum1_ip5,nLostLum2_ip5,ringparams['nMacro'],ringparams['nMacro'],
                                                                      pnumber1,x1,px1,y1,py1,t1,pt1,ex1,ey1,betax1,betay1,bs5,lumimax5,theta5,
                                                                      pnumber2,x2,px2,y2,py2,t2,pt2,ex2,ey2,betax2,betay2)
                nLostLum1_ip5 = nLostLum1_ip5 + len(x1) - kept1
                nLostLum2_ip5 = nLostLum2_ip5 + len(x2) - kept2
                x1  = x1[:kept1]
                px1 = px1[:kept1]
                y1  = y1[:kept1]
                py1 = py1[:kept1]
                t1  = t1[:kept1]
                pt1 = pt1[:kept1]
                x2  = x2[:kept2]
                px2 = px2[:kept2]
                y2  = y2[:kept2]
                py2 = py2[:kept2]
                t2  = t2[:kept2]
                pt2 = pt2[:kept2]
                hglassfacOld[2] = hgfac
            # IP8
            if bunch4:
                hgold = hglassfacOld[3]
                hgfac,betass8,kept1,kept2,lumip8 = fcollision.collision6a(collisiondict['sigI'],ringparams['thib'],collisiondict['longIntBins'],
                                                                      ringparams['timeRatio'],ringparams['circ'],np.pi,const.c,12673,
                                                                      hgold,nLostLum1_ip8,nLostLum4_ip8,ringparams['nMacro'],ringparams['nMacro'],
                                                                      pnumber1,x1,px1,y1,py1,t1,pt1,ex1,ey1,betax1,betay1,bs8,lumimax8,theta8,
                                                                      pnumber4,x4,px4,y4,py4,t4,pt4,ex4,ey4,betax2,betay2)
                nLostLum1_ip8 = nLostLum1_ip8 + len(x1) - kept1
                nLostLum4_ip8 = nLostLum4_ip8 + len(x4) - kept2
                x1  = x1[:kept1]
                px1 = px1[:kept1]
                y1  = y1[:kept1]
                py1 = py1[:kept1]
                t1  = t1[:kept1]
                pt1 = pt1[:kept1]
                x4  = x4[:kept2]
                px4 = px4[:kept2]
                y4  = y4[:kept2]
                py4 = py4[:kept2]
                t4  = t4[:kept2]
                pt4 = pt4[:kept2]
                hglassfacOld[3] = hgfac

            nLostLum1 = nLostLum1_ip1 + nLostLum1_ip5 + nLostLum1_ip2 + nLostLum1_ip8
            nLostLum2 = nLostLum2_ip1 + nLostLum2_ip5
            nLostLum3 = nLostLum3_ip2
            nLostLum4 = nLostLum4_ip8
            if (iwrite==1):
                print nLostLum1_ip1 , nLostLum1_ip5 , nLostLum1_ip2 , nLostLum1_ip8
                print nLostLum1,nLostLum2,nLostLum3,nLostLum4
        if (collisiondict['collRoutine']=='1d'):
            print 'under construction'
            quit()
    #first turn already written
    if ((kturns==1)&(switches['collisionSwitch'])):
        iwrite=0
        if bunch2:
            dflumip1 = cwrite.writeLumi(lumip1,hglassfacOld[0],dflumip1,kturns,nturns,eqtime,betass1)
        if bunch3:
            dflumip2 = cwrite.writeLumi(lumip2,hglassfacOld[1],dflumip2,kturns,nturns,eqtime,betass2)
        if bunch2:
            dflumip5 = cwrite.writeLumi(lumip5,hglassfacOld[2],dflumip5,kturns,nturns,eqtime,betass5)
        if bunch4:
            dflumip8 = cwrite.writeLumi(lumip8,hglassfacOld[3],dflumip1,kturns,nturns,eqtime,betass8)
    #     write output if desired
    if (iwrite==1):
        #call writemomentsshort(np1,y1,py1,t1,pt1,x1,px1,1)
        # call writemomentsshort(np2,y2,py2,t2,pt2,x2,px2,2)
        df1 = pd.DataFrame(columns=dfparticles1.columns)
        df2 = pd.DataFrame(columns=dfparticles2.columns)
        df3 = pd.DataFrame(columns=dfparticles3.columns)
        df4 = pd.DataFrame(columns=dfparticles4.columns)
        if bunch1:
                df1['x'] = x1
                df1['px'] = px1
                df1['y'] = y1
                df1['py'] = py1
                df1['t'] = t1
                df1['pt'] = pt1
        if bunch2:
                df2['x'] = x2
                df2['px'] = px2
                df2['y'] = y2
                df2['py'] = py2
                df2['t'] = t2
                df2['pt'] = pt2
        if bunch3:
                df3['x'] = x3
                df3['px'] = px3
                df3['y'] = y3
                df3['py'] = py3
                df3['t'] = t3
                df3['pt'] = pt3
        if bunch4:
                df4['x'] = x4
                df4['px'] = px4
                df4['y'] = y4
                df4['py'] = py4
                df4['t'] = t4
                df4['pt'] = pt4
        if switches['writeAllCoordSwitch']:   #     write all coordinates to file
            if bunch1:
                dfparticles1.to_csv('all-coord-b1-'+str(kturns)+'.csv',index=False)
            if bunch2:
                dfparticles2.to_csv('all-coord-b15-'+str(kturns)+'.csv',index=False)
            if bunch3:
                dfparticles3.to_csv('all-coord-b2-'+str(kturns)+'.csv',index=False)
            if bunch4:
                dfparticles4.to_csv('all-coord-b8-'+str(kturns)+'.csv',index=False)
           
        if switches['writeMountSwitch']:
            if bunch1:
                dfwritemount1 = cwrite.mountainr(ringparams['nMacro'],startingconditions['npart1'],avgline1,ringparams['qatomb1'],
                                         ringparams['nwrite'],kturns,ringparams['thib'],ibsdict['nbins'],dfwritemount1)
            if bunch2:
                dfwritemount2 = cwrite.mountainr(ringparams['nMacro'],startingconditions['npart2'],avgline2,ringparams['qatomb2'],
                                         ringparams['nwrite'],kturns,ringparams['thib'],ibsdict['nbins'],dfwritemount2)
            if bunch3:
                dfwritemount3 = cwrite.mountainr(ringparams['nMacro'],startingconditions['npart3'],avgline3,ringparams['qatomb2'],
                                         ringparams['nwrite'],kturns,ringparams['thib'],ibsdict['nbins'],dfwritemount3)
            if bunch4:
                dfwritemount4 = cwrite.mountainr(ringparams['nMacro'],startingconditions['npart4'],avgline4,ringparams['qatomb2'],
                                         ringparams['nwrite'],kturns,ringparams['thib'],ibsdict['nbins'],dfwritemount3)
              
            if bunch1:
                dfwritemounttrans1 = cwrite.mountTransv(avglinex1,avgliney1,ex1,ey1,betax1,betay1,ibsdict['nbins'],
                                                kturns,dfwritemounttrans1)
            if bunch2:
                dfwritemounttrans2 = cwrite.mountTransv(avglinex2,avgliney2,ex2,ey2,betax2,betay2,ibsdict['nbins'],
                                                kturns,dfwritemounttrans2)
            if bunch3:
                dfwritemounttrans3 = cwrite.mountTransv(avglinex3,avgliney3,ex3,ey3,betax2,betay2,ibsdict['nbins'],
                                                kturns,dfwritemounttrans3)
            if bunch4:
                dfwritemounttrans4 = cwrite.mountTransv(avglinex4,avgliney4,ex4,ey4,betax2,betay2,ibsdict['nbins'],
                                                kturns,dfwritemounttrans4)
           
        if bunch1:
            dfemit1 = cwrite.writeEmi(ex1,ey1,df1['t'],df1['pt'],betar1,ringparams['gamma1'],
                              ringparams['aatomb1'],ringparams['qatomb1'],eqtime,kturns,ringparams['nturns'],dfemit1)
            dfint1,nLostLumSum1,nLostDebunchSum1,nLostBetaSum1,nLostMomSum1  = cwrite.writeNb(df1,ringparams['nMacro'],startingconditions['npart1'], nLostLum1,
                       nLostLumSum1,nLostDebunch1,nLostDebunchSum1,nLostBeta1,nLostBetaSum1,nLostMom1,nLostMomSum1,
                       kturns,ringparams['nturns'],eqtime,dfint1)
    
        if bunch2:    
            dfemit2 = cwrite.writeEmi(ex2,ey2,df2['t'],df2['pt'],betar2,ringparams['gamma2'],
                              ringparams['aatomb2'],ringparams['qatomb2'],eqtime,kturns,ringparams['nturns'],dfemit2)
            dfint2,nLostLumSum2,nLostDebunchSum2,nLostBetaSum2,nLostMomSum2  = cwrite.writeNb(df2,ringparams['nMacro'],startingconditions['npart2'], nLostLum2,
                      nLostLumSum2,nLostDebunch2,nLostDebunchSum2,nLostBeta2,nLostBetaSum2,nLostMom2,nLostMomSum2,
                       kturns,ringparams['nturns'],eqtime,dfint2)
    
        if bunch3:
            dfemit3 = cwrite.writeEmi(ex3,ey3,df3['t'],df3['pt'],betar2,ringparams['gamma2'],
                              ringparams['aatomb2'],ringparams['qatomb2'],eqtime,kturns,ringparams['nturns'],dfemit3)
            dfint3,nLostLumSum3,nLostDebunchSum3,nLostBetaSum3,nLostMomSum3  = cwrite.writeNb(df3,ringparams['nMacro'],startingconditions['npart3'], nLostLum3,
                       nLostLumSum3,nLostDebunch3,nLostDebunchSum3,nLostBeta3,nLostBetaSum3,nLostMom3,nLostMomSum3,
                       kturns,ringparams['nturns'],eqtime,dfint3)


        if bunch4:
            dfemit4 = cwrite.writeEmi(ex4,ey4,df4['t'],df4['pt'],betar2,ringparams['gamma2'],
                              ringparams['aatomb2'],ringparams['qatomb2'],eqtime,kturns,ringparams['nturns'],dfemit4)
            dfint4,nLostLumSum4,nLostDebunchSum4,nLostBetaSum4,nLostMomSum4  = cwrite.writeNb(df4,ringparams['nMacro'],startingconditions['npart4'], nLostLum4,
                       nLostLumSum4,nLostDebunch4,nLostDebunchSum4,nLostBeta4,nLostBetaSum4,nLostMom4,nLostMomSum4,
                       kturns,ringparams['nturns'],eqtime,dfint4)
        

        if (switches['collisionSwitch']):
            if bunch2:
                dflumip1 = cwrite.writeLumi(lumip1,hglassfacOld[0],dflumip1,kturns,nturns,eqtime,betass1)
            if bunch3:
                dflumip2 = cwrite.writeLumi(lumip2,hglassfacOld[1],dflumip2,kturns,nturns,eqtime,betass2)
            if bunch2:
                dflumip5 = cwrite.writeLumi(lumip5,hglassfacOld[2],dflumip5,kturns,nturns,eqtime,betass5)
            if bunch4:
                dflumip8 = cwrite.writeLumi(lumip8,hglassfacOld[3],dflumip8,kturns,nturns,eqtime,betass8)
        print nLostLum1,nLostLum2,nLostLum3,nLostLum4
        nLostLum1       =0
        nLostLum2       =0
        nLostDebunch1   =0
        nLostDebunch2   =0
        nLostMom1       =0
        nLostBeta1      =0
        nLostMom2       =0
        nLostBeta2      =0
        nLostLum3       =0
        nLostDebunch3   =0
        nLostMom3       =0
        nLostBeta3      =0
        nLostLum4       =0
        nLostDebunch4   =0
        nLostMom4       =0
        nLostBeta4      =0
        nLostLum1_ip1 = 0
        nLostLum2_ip1 = 0
        nLostLum1_ip2 = 0
        nLostLum3_ip2 = 0
        nLostLum1_ip5 = 0
        nLostLum2_ip5 = 0
        nLostLum1_ip8 = 0
        nLostLum4_ip8 = 0
        
print 'took :', time.clock()-st
print 'end at:',time.ctime()

dflumip1.to_csv(fn+'lumip1.csv',index=False)
dflumip2.to_csv(fn+'lumip2.csv',index=False)
dflumip5.to_csv(fn+'lumip5.csv',index=False)
dflumip8.to_csv(fn+'lumip8.csv',index=False)

dfemit1.to_csv(fn+'emit1.csv',index=False)
dfint1.to_csv(fn+'int1.csv',index=False)
dfemit2.to_csv(fn+'emit2.csv',index=False)
dfint2.to_csv(fn+'int2.csv',index=False)
dfemit3.to_csv(fn+'emit3.csv',index=False)
dfint3.to_csv(fn+'int3.csv',index=False)
dfemit4.to_csv(fn+'emit4.csv',index=False)
dfint4.to_csv(fn+'int4.csv',index=False)

dfibs1.to_csv(fn+'ibs1.csv',index=False)
dfibs2.to_csv(fn+'ibs2.csv',index=False)
dfibs3.to_csv(fn+'ibs3.csv',index=False)
dfibs4.to_csv(fn+'ibs4.csv',index=False)
