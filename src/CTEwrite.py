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

mountainnrcols         = ['kturns','kdt2','afvgline(k)']
mounttransvcols        = ['kturns','xc','yc','avglinex(k)','avgliney(k)']
writeemicols           = ['kturns','kturns/nturns*eqTime','epsx1','epsy1','epsl1','sigt1','dponp1',
                          'epsx2','epsy2','epsl2','sigt2','dponp2']
writenbcols            = ['kturns','kturns/nturns*eqTime',
                          'np1','np1/np10*pnumber1','nLostLum1','nLostLumSum1',
                          'nLostDebunch1','nLostDebunchSum1','nLostBeta1','nLostBetaSum1',
                          'nLostMom1','nLostMomSum1','np2','np2/np20*pnumber2', 'nLostLum2',
                          'nLostLumSum2','nLostDebunch2','nLostDebunchSum2','nLostBeta2',
                          'nLostBetaSum2','nLostMom2','nLostMomSum2']
writecoordcols         = ['x(k)','px(k)','y(k)','py(k)','t(k)','pt(k)']

def mountTransv(avglinex,avgliney,emix,emiy,betax,betay,nBins,kturns,dfoutxy):
    xmn  = -5. * np.sqrt(betax*emix)
    ymn  = -5. * np.sqrt(betay*emiy)
    dx   = -2.0 * xmn/float(nBins)
    dy   = -2.0 * ymn/float(nBins)
    ntot = 0.
    xc = []
    yc = []
    for k in range(nBins):
        xc.append(xmn + k * dx)
        yc.append(ymn + k * dy)
    dfdata  = [[kturns,xc[k],yc[k],avglinex[k],avgliney[k]] for k in range(nBins)]
    dfrow   = pd.DataFrame(dfdata,columns=mounttransvcols)
    dfoutxy = dfoutxy.append(dfrow,ignore_index=True) 
    return dfoutxy

def writeEmi(epsx1,epsy1,t1,pt1,relbeta1,gamma1,aatom1,qatom1,eqTime,kturns,nturns,dfout):
    sigt1  = 0.
    sigPt1 = 0.
    sigt2  = 0.
    sigPt2 = 0.

#      determine long. stand. dev.
#      beam 1
    sigt1  = np.sum([i**2 for i in t1])
    sigPt1 = np.sum([i**2 for i in pt1])
    sigt1  = np.sqrt(sigt1/float(len(t1)))
    sigPt1 = np.sqrt(sigPt1/float(len(t1)))
    dponp1 = sigPt1/(relbeta1 * relbeta1 * gamma1) # 1 sigma fractional momentum deviation
    epsl1  = sigt1*sigPt1*0.931494e9*aatom1*np.pi/float(qatom1) # long. emittance in eV s/charge at 1 sigma. not exact.

    dfdata  = [kturns,float(kturns)/float(nturns)*eqTime,epsx1,epsy1,epsl1,sigt1,dponp1]
    dfrow   = pd.Series(dfdata,index=dfout.columns)
    dfout = dfout.append(dfrow,ignore_index=True) 
    
    print '{:6d}{:15.5e}{:15.5e}{:15.5e}{:15.5e}{:15.5e}{:15.5e}'.format(kturns,float(kturns)/float(nturns)*eqTime,epsx1,epsy1,epsl1,sigt1,dponp1)
    return dfout

def writeNb(df,np10,pnumber1,nLostLum1,nLostLumSum1,nLostDebunch1,nLostDebunchSum1,
            nLostBeta1,nLostBetaSum1,nLostMom1,nLostMomSum1,kturns,nturns,eqTime,dfout):
    np1 = len(df)
    
    nLostLumSum1     = nLostLumSum1     + nLostLum1 # total n.o. lost particles over all turns from lumi
    nLostDebunchSum1 = nLostDebunchSum1 + nLostDebunch1 # total n.o. lost particles over all turns from debunching

    nLostBetaSum1    = nLostBetaSum1    + nLostBeta1
    nLostMomSum1     = nLostMomSum1     + nLostMom1

    dfdata = [kturns,float(kturns)/float(nturns)*eqTime,np1,float(np1)/float(np10)*pnumber1, 
              nLostLum1,nLostLumSum1,nLostDebunch1,nLostDebunchSum1,
              nLostBeta1,nLostBetaSum1,nLostMom1,nLostMomSum1]
    dfrow   = pd.Series(dfdata,index=dfout.columns)
    dfout = dfout.append(dfrow,ignore_index=True)
    # reset intermediate sums
    return dfout,nLostLumSum1,nLostDebunchSum1,nLostBetaSum1,nLostMomSum1

def mountainr(np0,pnumber,avgline,qatom,nwrite0,kturns,thib,nBins,dfout):
    dt2       = thib/nBins
    # put into units of current
    currcoeff = 1.602e-19*pnumber*qatom/(nwrite0*dt2)
    currcoeff = currcoeff/float(np0)
    peaki     = 0
    jdex      = 0
    
    for k in range(nBins):
        currk = avgline[k]*currcoeff
        if(currk>=peaki):
            peaki=currk
            jdex=k
        avgline[k]=currk
     
    # get the full width half max of the peak
    for k in range(jdex,nBins):
        pk = avgline[k]
        if (pk<=peaki/2.):
            break
    jhi=k
    
    for k in range(jdex,-1,-1):
        pk = avgline[k]
        if (pk <= peaki/2.0):
            break
    
    jlo       = k

    dfdata = pd.DataFrame([[kturns,k*dt2,avgline[k]] for k in range(nBins)],columns=mountainnrcols)
    dfout = dfout.append(dfdata,ignore_index=True)
    return dfout

def writeLumi(lumi,redfac,dfout,kturns,nturns,eqTime,betas):
    print 'luminosity=',lumi,'cm^-2 s^-1, ','geom. reduction factor=',redfac
    dfdata = [kturns,float(kturns)/float(nturns)*eqTime,lumi,redfac,betas]
    dfrow   = pd.Series(dfdata,index=dfout.columns)
    dfout = dfout.append(dfrow,ignore_index=True)
    return dfout

