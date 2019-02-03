! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       IBS ROUTINE
!
!   DEPENDENCIES:
!        IBSLONG ->
!                  ibspiwismooth.f95
!                  ibspiwilattice.f95
!                  ibsmodpiwilattice.f95
!                  ibsbane.f95
!                  ibsinterpolat.f95
!                  ibsnagaitsev.f95
! ---------------------------------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------------------------
! pnumber   : total number of particles in the bunch
! np0       : number of macro particles in the bunch
! circ      : accelerator circumference
! iwrite    : write distributions to file
! nBins     : number of bins to bin the particles in
! thib      : height of bunch longitudinally
! betar     : relativistic beta
! gamma0    : relativistic gamma
! betax     : ring average horizontal beta
! betay     : ring average vertical beta
! aatom     : atomic number
! qatom     : charge
! ibsrow    : output ibs growthrates to use in output files in main program
! trev      : revolution time
! betx      : horizontal beta madx
! bety      : vertical beta madx
! dispx     : DX madx
! dispy     : DY madx
! dxp       : DPX madx
! dyp       : DPY madx
! alfx      : ALFX madx
! leng      : L madx
! gtab      : Bane input table
! coulomb   : coulomb constant, input nagaitsev
! xmin      : x min in gtab
! ymin      : y min in gtab
! zmin      : z min in gtab
! xmax      : x max in gtab
! ymax      : y max in gtab
! zmax      : z max in gtab
! Ap        : interpolation table
! Ax        : interpolation table
! Ay        : interpolation table
! pnumInterp: scaling number
! coupleibs : coupling between x-y planes
! nturns    : number of turns to perform
! kturns    : turn number
! eqTime    : sim vs real time
! fracibstot: scaling factor that can be set to increase or decrease ibs
! clight    : speed of light
! pi        : constant
! tratio    : ratio real vs sim time
! iseed     : random number, needs to be inout as the value is updated by the random generator
! np        : number of macro particles in bunch
! nx        : table size
! ny        : table size
! xBins     : number of bins horizontal
! yBins     : number of bins vertical
! zBins     : number of bins longitudinal
! nElem     : number of elements (madx)
! ----------------------------------------------------------------------------------------------------------

subroutine ibslong(pnumber,px,py,t,pt,epsx,epsy,np0,&
circ,iwrite,nBins,thib,avgline,betar,gamma0,gammat,betax,betay,aatom,qatom,ibsrow,trev,&
betx,bety,dispx,dispy,dxp,dyp,alfx,alfy,leng,&
gtab,coulomb,xmin,xmax,ymin,ymax,zmin,zmax,Ap,Ax,Ay,pnumInterp,&
ibsMethod,coupleIBS,nturns,kturns,eqTime,fracibstot,clight,pi,tratio,iseed,np,nx,ny,xBins,yBins,zBins,nElem)

!     pnumber is the total number of particles in the bunch
!     pnnorm is the bunch population that IBS should be normalized to.
!     In the original code, the IBS strength is a factor pnumber/np0 higher,
!     but we need the same normalization factor for BOTH bunches, otherwise
!     their equivalent time in the machine will not match!

character(10), intent(in) :: ibsMethod
!f2py intent(in) :: ibsMethod
integer, intent(in) :: np,iwrite,nBins,nturns,kturns,nElem,np0,nx,ny
!f2py intent(in) :: np,iwrite,nBins,nturns,kturns,iseed,nElem
double precision, dimension(nx,ny), intent(in) :: gtab
!f2py intent(in) :: gtab
double precision, intent(inout) :: iseed
!f2py intent(in,out) :: iseed
integer, intent(in) :: xBins,yBins,zBins,pnumInterp
!f2py intent(in) :: xBins,yBins,zBins,pnumInterp
double precision, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
!f2py intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax

double precision, dimension(xBins,yBins,zBins), intent(in) :: Ap,Ax,Ay
!f2py intent(in) :: Ap,Ax,Ay

double precision, dimension(nElem), intent(in) :: betx,bety,dispx,dispy,dxp,dyp,alfx,alfy,leng
!f2py intent(in) :: betx,bety,dispx,dispy,dxp,dyp,alfx,alfy,leng

double precision, dimension(np), intent(inout) :: px,py,t,pt
!f2py intent(out) :: px,py,t,pt

double precision, intent(in) :: betar,gamma0,betax,betay,aatom,qatom,gammat,pnumber,coupleIBS,epsx,epsy
!f2py intent(in) :: betar,gamma0,betax,betay,aatom,qatom,gammat,pnumber,coupleIBS,epsx,epsy

double precision, intent(in) :: thib,coulomb,trev,eqTime,fracibstot,clight,tratio,pi,circ
!f2py intent(in) :: thib,coulomb,trev,eqTime,fracibstot,clight,np0,tratio,pi

double precision, dimension(5), intent(out) :: ibsrow
!f2py intent(out) :: ibsrow

!double precision,external :: ran3
double precision :: alfap,alfax,alfay,dx,dy,dg,sigs,dponp
double precision :: rmsg,rmst,tk,gk,dtsamp2,thib2,coeffs,coeffx,coeffy,rmsx,rmsy
double precision :: grv1,grv2,alfaAvg,coeffmult,denlonn,denlon2k,alfap0,alfax0,alfay0
double precision, dimension(nBins) :: denlon2,avgline

integer :: k, nk
!     this routine copied as is from Mike, except small changes:
!     --using timeRatio parameter to scale the strength instead of pnumber1/np10
!     --changed notation x2->x and x->y for consistency with the rest of the code
!     --implemented option of coupling that averages the growth rates over x and y

  rmsg=0
  rmst=0
! zero out the array for line density
  do k=1,nBins
     denlon2(k)=0
  enddo
  
  dtsamp2 = thib/float(nBins)
  thib2   = thib/2.
  do  k=1,np
!     t=0 is stable fixed point
         tk = t(k)
         gk = pt(k)
         rmst = rmst + tk*tk
         rmsg = rmsg + gk*gk
         tk = (tk+thib2)/dtsamp2
! just use nearest integer
         nk = nint(tk)
         if(nk.le.0)nk=1
         if(nk.gt.nBins)nk=nBins
         denlon2(nk)=denlon2(nk)+1
      enddo
! add into avgline(k)
      do k=1,nBins
         avgline(k) = avgline(k)+denlon2(k)
      enddo

! nothing for no ibs
      if(fracibstot.le.0)return

      rmst = sqrt(rmst/np)
      rmsg = sqrt(rmsg/np)
      dponp = rmsg/(betar*betar*gamma0)
      sigs = rmst*clight*betar

      rmsx = sqrt(epsx*betax)
      rmsy = sqrt(epsy*betay)

!      write(6,*) ibsMethod
      if (ibsMethod.eq.'piwiSmooth') then
         call ibsPiwSmooth(pnumber,epsx,epsy,sigs,pi,gamma0,gammat,circ,betar,&
betax,betay,clight,qatom,aatom,dponp,alfap0,alfax0,alfay0)
      elseif (ibsMethod.eq.'piwLattice') then
         call ibsPiwLattice(pnumber,epsx,epsy,sigs,dponp,pi,circ,&
clight,qatom,aatom,betar,betx,bety,dispx,leng,nElem,gamma0,alfap0,alfax0,alfay0)
      elseif (ibsMethod.eq.'modPiwLatt') then
         call modPiwLattice(pnumber,epsx,epsy,sigs,dponp,pi,circ,&
clight,qatom,aatom,betar,betx,bety,dispx,dxp,leng,alfx,nElem,gamma0,alfap0,alfax0,alfay0)
      elseif (ibsMethod.eq.'baneApprox') then
         call ibsBane(pnumber,epsx,epsy,sigs,dponp,qatom,aatom,gamma0,circ,clight,&
betx,bety,dispx,dispy,dxp,dyp,alfx,alfy,leng,alfap0,alfax0,alfay0,nElem,gtab,nx,ny)
      elseif (ibsMethod.eq.'interpolat') then
         call ibsInterpolat(pnumber,epsx,epsy,sigs,dponp,xmin,&
xmax,ymin,ymax,zmin,zmax,xBins,yBins,zBins,Ap,Ax,Ay,pnumInterp,alfap0,alfax0,alfay0)
      elseif (ibsMethod.eq.'nagaitsev') then
         call nagaitsevIBSlattice(pnumber,epsx,epsy,dponp,sigs,qatom,aatom,nElem,dxp,alfx,betx,dispx,bety,&
gamma0,circ,clight,pi,coulomb,alfap0,alfax0,alfay0,leng,betar)
      else
         write(*,*) 'Stop - unknown ibs method: ',ibsMethod
         stop
      endif

!      write(*,*) alfap0,alfax0,alfay0
      if (iwrite.eq.1) then
         ibsrow(1) = kturns 
         ibsrow(2) = real(kturns)/real(nturns)*eqTime
         ibsrow(3) = 1/alfap0/3600/2
         ibsrow(4) = 1/alfax0/3600/2
         ibsrow(5) = 1/alfay0/3600/2          
      endif

! alphas are amplitude growth rates for pnumber particles
! and no losses
! correct for small number of macro particles, normalize to pnnorm
! mult by 2 so it a kick for emittance
! reduce by input fraction
      alfap = 2*fracibstot*float(np)*tratio*alfap0/float(np0)
      alfay = 2*fracibstot*float(np)*tratio*alfay0/float(np0)
      alfax = 2*fracibstot*float(np)*tratio*alfax0/float(np0)
! with rms growth rate for 1 dim of shm need to kick
! sqrt(2.) times harder
      if(alfap.gt.0)then
          coeffs = sqrt(6*alfap*trev)*rmsg
      else
          coeffs=0
      endif
      if(coupleIBS.eq.0) then   ! use separate growth rates in x and y. Original version from Mike
         if(alfay.gt.0)then
            coeffy = sqrt(6*alfay*trev)*rmsy
         else
            coeffy = 0
         endif
         if(alfax.gt.0)then
            coeffx = sqrt(6*alfax*trev)*rmsx
         else
            coeffx = 0
         endif
      else                      ! use average alfa in x and y. corresponds to COOL07 paper and early LHC simulations
         alfaAvg = (alfay+alfax)/2
         if (alfaAvg.gt.0) then
            coeffy = sqrt(6*alfaAvg*trev)*rmsy
            coeffx = sqrt(6*alfaAvg*trev)*rmsx
         else
            coeffy=0
            coeffx=0
         endif
      endif

      coeffmult = 2*sqrt(pi)/(np*dtsamp2*clight) * sigs
      denlonn   = 0
      do k=1,np
          tk = (t(k)+thib2)/dtsamp2
! just use nearest integer
          nk = nint(tk)
          if(nk.le.0)nk=1
          if(nk.gt.nBins)nk=nBins

          denlon2k = denlon2(nk)*coeffmult
          denlonn = denlonn+denlon2k
          denlon2k = sqrt(denlon2k)
! that easy??
!         dy = denlon2k*coeffy*(2*ran3(iseed)-1)
!         dg = denlon2k*coeffs*(2*ran3(iseed)-1)
         call get2gaussrv(iseed,grv1,grv2)
         dy = denlon2k*coeffy*grv1
         dg = denlon2k*coeffs*grv2
         call get2gaussrv(iseed,grv1,grv2)
         dx = denlon2k*coeffx*grv1

         py(k)=py(k)+dy
         px(k)=px(k)+dx
         pt(k)=pt(k)+dg
      enddo
      return
end

subroutine get2gaussrv(iseed,grv1,grv2)

double precision, intent(inout) :: iseed
!f2py intent(in,out) :: iseed
double precision, intent(out):: grv1,grv2
!f2py intent(out) :: grv1,grv2

double precision :: r1,r2,facc,amp

!double precision, external :: ran3

44 continue
   call random_number(iseed)
   r1 = 2*iseed-1
   call random_number(iseed)
   r2 = 2*iseed-1
   amp = r1**2 + r2**2
 
   if ((amp.ge.1).or.(amp.lt. 1.e-8)) go to 44
   facc = sqrt(-2.*log(amp)/amp)
   grv1 = r1*facc/sqrt(3.)
   grv2 = r2*facc/sqrt(3.)
   return
end subroutine get2gaussrv

!include 'ran3.f95'
include 'ibspiwismooth.f95'
include 'ibspiwilattice.f95'
include 'ibsmodpiwilattice.f95'
include 'ibsbane.f95'
include 'ibsinterpolat.f95'
include 'ibsnagaitsev.f95'
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
