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
from fnmatch import fnmatch,fnmatchcase
import tarfile
import glob
from dask import delayed
from dask.dot import dot_graph
from dask.diagnostics import ProgressBar
pbar=ProgressBar()
pbar.register()

electroncharge                  = const.physical_constants['elementary charge'][0]

def raddampingapprox(ringparams,radiationdamping,betax,betay):
    alfax = 0.
    alfay = 0.
    gamax = (1.+alfax**2)/betax
    gamay = (1.+alfay**2)/betay
    Dx    = ringparams['circ']/(2*np.pi*ringparams['gammat']**2)
    Dy    = 0.
    DxP   = 0.1             # should find an approximation formula. However not very important
    DyP   = 0.
    Hx    = (betax*DxP+2*alfax*Dx*DxP+gamax*Dx)
    Hy    = (betay*DyP+2*alfay*Dy*DyP+gamay*Dy)
#    define smooth approximation of radiation integrals
    I2    = 2*np.pi/radiationdamping['rho0']
    I3    = 2*np.pi/radiationdamping['rho0']**2
    I4x   = 0.0
    I4y   = 0.0
    I5x   = Hx*2*np.pi/radiationdamping['rho0']**2
    I5y   = Hy*2*np.pi/radiationdamping['rho0']**2
    return I2,I3,I4x,I4y,I5x,I5y

def raddampinglattic(ringparams,radiationdamping,dftfs):
    I2                = 0.
    I3                = 0.
    I4x               = 0.
    I4y               = 0.
    I5x               = 0.
    I5y               = 0.
    
    dftfscopy         = dftfs.copy()
    dftfscopy['rhoi'] = np.where(dftfscopy['ANGLE']!=0.0,dftfscopy['L']/dftfscopy['ANGLE'],0)
    dftfscopy['ki']   = np.where(dftfscopy['L']!=0.0,dftfscopy['K1L']/dftfscopy['L'],0)
    dftfscopy['I2']   = np.where(dftfscopy['rhoi']!=0.0,dftfscopy['L']/dftfscopy['rhoi']**2,0)
    dftfscopy['I3']   = np.where(dftfscopy['rhoi']!=0.0,dftfscopy['L']/dftfscopy['rhoi']**3,0)
    #  corrected to equations in accelerator handbook second edition p 220
    dftfscopy['I4x']  = np.where(dftfscopy['rhoi']!=0.0,
                                     ((dftfscopy['DY']/dftfscopy['rhoi']**3)*(1.+2.*dftfscopy['rhoi']**2*dftfscopy['ki']*dftfscopy['K1S']/dftfscopy['rhoi']))*dftfscopy['L'],
                                     0)
    dftfscopy['I4y']  = 0
    dftfscopy['gamax']= (1.+ dftfscopy['ALFX']**2)/dftfscopy['BETX']
    dftfscopy['gamay']= (1.+ dftfscopy['ALFY']**2)/dftfscopy['BETY']
    dftfscopy['Hx']   = dftfscopy['BETX']*dftfscopy['DPX']**2 +2.*dftfscopy['ALFX']*dftfscopy['DX']*dftfscopy['DPX']+dftfscopy['gamax']*dftfscopy['DX']**2
    dftfscopy['Hy']   = dftfscopy['BETY']*dftfscopy['DPY']**2 +2.*dftfscopy['ALFY']*dftfscopy['DY']*dftfscopy['DPY']+dftfscopy['gamay']*dftfscopy['DY']**2
    dftfscopy['I5x']  = dftfscopy['Hx']*2*np.pi/radiationdamping['rho0']**2*dftfscopy['L']
    dftfscopy['I5y']  = dftfscopy['Hy']*2*np.pi/radiationdamping['rho0']**2*dftfscopy['L']
    
    I2                = dftfscopy['I2'].sum()
    I3                = dftfscopy['I3'].sum()
    I4x               = dftfscopy['I4x'].sum()
    I4y               = dftfscopy['I4y'].sum()
    I5x               = dftfscopy['I5x'].sum()  
    I5y               = dftfscopy['I5y'].sum()
    return I2,I3,I4x,I4y,I5x,I5y

def raddamp(dfparticle,tradlong,tradperp,trev,siglong,sigperp,ringparams):
    # does radiation damping and quantum excitation once per turn
    # skip if transverse damping time is not positive
    if(tradperp<=0):
        return

    # damping time should scale down with the ratio of real turns/sim. turns
    # exact damping coefficient is exp(-trev/tradlong). expand to first order

    # get longitudinal radiation damping and excitation terms
    else:
        tratio          = ringparams['timeRatio']
        coeffdecaylong  = 1 - ((trev / tradlong) * tratio)
        
        # excitation uses a uniform deviate on [-1:1]
        coeffexcitelong = siglong * np.sqrt(3.) * np.sqrt(2 * (trev / tradlong) * tratio)
        
        # tradperp is the damping time for EMITTANCE, therefore need to multiply by 2
        # assume same damping in horizontal and vertical plane (I4x,I4y<<I2)
        coeffdecay      = 1 - ((trev /(2 * tradperp)) * tratio)
        
        # exact     coeffgrow= sigperp*sqrt(3.)*sqrt(1-coeffdecay**2)
        # but trev << tradperp so
        coeffgrow       = sigperp * np.sqrt(3.) * np.sqrt(2 * (trev /(2 * tradperp)) * tratio)

        df = dfparticle.copy()
        #print df.head()
        # not sure here, I think ran3(iseed) returns same number np.random not (no fixed seed)
        # longitudinal update
        df['pt'] = df['pt'] * coeffdecaylong + coeffexcitelong * np.random.uniform(-1,1,len(df))
        # transverse update
        df['x']  = coeffdecay * df['x']  + np.random.uniform(-1,1,len(df)) * coeffgrow 
        df['px'] = coeffdecay * df['px'] + np.random.uniform(-1,1,len(df)) * coeffgrow
        df['y']  = coeffdecay * df['y']  + np.random.uniform(-1,1,len(df)) * coeffgrow
        df['py'] = coeffdecay * df['py'] + np.random.uniform(-1,1,len(df)) * coeffgrow
        #print coeffdecay, coeffgrow
        #print df.head()
        return df


def raddampingconstfromI(I2,I3,I4x,I4y,I5x,I5y,beam,CalphaE3C,ringparams):
    assert(beam in [1,2])
    tradperp   = 1./(CalphaE3C*I2*(1.-I4x/I2))/2. # eq 9, but div by 2 to get rise time for emittance and not amplitude
    tradlong   = 1./(CalphaE3C*I2*(2.+(I4x+I4y)/I2))/2.
    if (ringparams['aatomb'+str(beam)] == 208):
        Cq    = 55./(32.*np.sqrt(3.0))*1.0546e-34*const.c / (ringparams['ion']*electroncharge*1e9)
      
    else:
        Cq    = 55./(32.*np.sqrt(3.0))*1.0546e-34*const.c / (ringparams['proton']*electroncharge*1e9)
    sigE0E0   = Cq *ringparams['gamma'+str(beam)]**2*I3/(2*I2+I4x+I4y)  # eq 18
    siglong   = ringparams['gamma'+str(beam)]*sigE0E0
    Jx        = 1.-I4x/I2
    Jy        = 1.-I4y/I2
    sigperp   = Cq*ringparams['gamma'+str(beam)]**2*I5x/(Jx*I2)
    return tradperp,tradlong,Cq,sigE0E0,siglong,Jx,Jy,sigperp