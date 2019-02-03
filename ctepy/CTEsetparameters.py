import numpy as np
import pandas as pd
import datetime
import time
from collections import OrderedDict
import functools
from scipy import constants as const


def settunechroma(beam,ringparams):
    psix = 2 * const.pi * ringparams['tunex'+str(beam)]
    psiy = 2 * const.pi * ringparams['tuney'+str(beam)]
    coeffchromx = ringparams['chromx'+str(beam)] * 2 * const.pi/(ringparams['gamma'+str(beam)]-1./ringparams['gamma'+str(beam)])
    coeffchromy = ringparams['chromy'+str(beam)] * 2 * const.pi/(ringparams['gamma'+str(beam)]-1./ringparams['gamma'+str(beam)])
    tpdqmin = 2 * const.pi * ringparams['dqmin'+str(beam)]
    return psix,psiy,coeffchromx,coeffchromy,tpdqmin

def set_parameters1(ringparams,collisiondict):
    # hourglass effect 
    hglassfacOld = np.ones(10)
    # reading betastarr
    betaSMin     = [collisiondict['ips_betaSX'+str(i)] for i in range(1,
                        np.array([1 for c in collisiondict.keys() if c[:3]=='ips']).sum()/5+1)]
    # relativistic betas for each beam
    betar1       = np.sqrt(1-1/ringparams['gamma1']**2)
    betar2       = np.sqrt(1-1/ringparams['gamma2']**2)
    # velocities
    vrev1        = const.c*betar1
    vrev2        = const.c*betar2
    # revolution time or period
    trev1        = ringparams['circ']/vrev1      
    trev2        = ringparams['circ']/vrev2
    # rf perods 
    trf1         = trev1/ringparams['nharm2b1']
    trf2         = trev2/ringparams['nharm2b2']
    # rf frequencies (should be equal)
    frf1         = 1/trf1  # 400 MHz
    frf2         = 1/trf2  
    # angle omegas
    omegarf1     = 2*np.pi/trf1
    omegarf2     = 2*np.pi/trf2
    # slipfactors
    eta1         = 1/ringparams['gammat']**2 - 1/ringparams['gamma1']**2  # slipfactor see Lee (2.162)
    eta2         = 1/ringparams['gammat']**2 - 1/ringparams['gamma2']**2  #  slipfactor see Lee (2.162)
    print 'RF freq (b1,b2):',frf1,frf2

    eqtime1      = ringparams['nturns'] *trev1*ringparams['timeRatio']/3600.
    eqtime2      = ringparams['nturns'] *trev2*ringparams['timeRatio']/3600.
    print 'equivalent time (hours) b1 = ', eqtime1
    print 'equivalent time (hours) b2 = ', eqtime2

    fneqtime   = open('eqtime.out','w')
    fneqtime.write(str(eqtime1)+'\n')
    fneqtime.write(str(eqtime2)+'\n')
    fneqtime.close()
    return hglassfacOld,betaSMin,betar1,betar2,vrev1,vrev2,trev1,trev2,trf1,trf2,frf1,frf2,omegarf1,omegarf2,eta1,eta2,eqtime1,eqtime2

def longdyncoeff(ringparams,beam,omegarf,v,trev,eta,betar):
    assert(beam in [1,2])
    astr = 'aatomb'+str(beam)
    qstr = 'qatomb'+str(beam)
    gstr = 'gamma'+str(beam)
    # pcoeff - dp/dn
    if ringparams[astr]== 208:
        pcoeff = ringparams[qstr] * v * omegarf/(ringparams['ion']*1.0e9)
    else:
        pcoeff = ringparams[qstr] * v * omegarf/(ringparams['proton']*1e9)
    # tcoeff - dt/dn
    tcoeff = trev*eta/(betar**2*ringparams[gstr])
    rcoeff = pcoeff/tcoeff
    if(rcoeff >= 0):
        print 'wrong sign of voltage for sign of eta b1'
        quit()
    return pcoeff,tcoeff,rcoeff

def set_parameters2(ringparams,startingconditions,omegarf1,trev1,eta1,betar1,omegarf2,trev2,eta2,betar2):
    # main average betax and betay for both beams
    betax1=ringparams['circ']/(2*np.pi*ringparams['tunex1'] )
    betay1=ringparams['circ']/(2*np.pi*ringparams['tuney1'] )
    betax2=ringparams['circ']/(2*np.pi*ringparams['tunex2'] )
    betay2=ringparams['circ']/(2*np.pi*ringparams['tuney2'] )

    # part of code that I don't understand very well, have to check older versions with dual vrf system active
    if (ringparams['vrfb1']==0.0):
        v1 = ringparams['vrf2b1']
    else:
        v1 = ringparams['vrf2b1']
        
    if (ringparams['vrfb2']==0.0):
        v2 = ringparams['vrf2b2']
    else:
        v2 = ringparams['vrf2b2']

    pcoeff1,tcoeff1,rcoeff1 = longdyncoeff(ringparams,1,omegarf1,v1,trev1,eta1,betar1)
    pcoeff2,tcoeff2,rcoeff2 = longdyncoeff(ringparams,2,omegarf2,v2,trev2,eta2,betar2)
    
    dgamma_hat1 = startingconditions['tauhat1']*np.sqrt(-rcoeff1)
    dgamma_hat2 = startingconditions['tauhat2']*np.sqrt(-rcoeff2)
    
    # tauhat1=half bucket length (seconds) for beam 1 of the RF system with shortest wavelength
    print 'amp of de/e b1= ',dgamma_hat1/ringparams['gamma1']
    print 'amp of de/e b2= ',dgamma_hat2/ringparams['gamma2']
    print 'max change in tau per turn b1= ',dgamma_hat1*tcoeff1
    print 'max change in tau per turn b2= ',dgamma_hat2*tcoeff2

    # synchrotron tune
    tunes1 = np.sqrt(-pcoeff1*tcoeff1)/(2*np.pi)
    tunes2 = np.sqrt(-pcoeff2*tcoeff2)/(2*np.pi)

    tunesexact1 = np.arccos(1+pcoeff1*tcoeff1/2.)/(2*np.pi)
    tunesexact2 = np.arccos(1+pcoeff2*tcoeff2/2.)/(2*np.pi)
    print 'synchrotron frequency b1= ',tunes1/trev1
    print 'synchrotron frequency b2= ',tunes2/trev2
    print 'revolution period b1= ',trev1
    print 'revolution period b2= ',trev2
    print 'exact synchrotron tune b1 = ',tunesexact1
    print 'exact synchrotron tune b2 = ',tunesexact2

    # normalization for longitudinal profiles
    nwrite0 = ringparams['nwrite']

    dftfs1 = pd.read_csv(ringparams['tfsb1'],skiprows=45,nrows=1,delim_whitespace=True)
    dftfs1 = pd.read_csv(ringparams['tfsb1'],skiprows=47,delim_whitespace=True,names=dftfs1.columns[1:])
    dftfs2 = pd.read_csv(ringparams['tfsb2'],skiprows=45,nrows=1,delim_whitespace=True)
    dftfs2 = pd.read_csv(ringparams['tfsb2'],skiprows=47,delim_whitespace=True,names=dftfs2.columns[1:])
    return betax1,betay1,betax2,betay2,pcoeff1,tcoeff1,rcoeff1,pcoeff2,tcoeff2,rcoeff2,dgamma_hat1,dgamma_hat2, tunes1,tunes2,tunesexact1,tunesexact2,dftfs1,dftfs2
