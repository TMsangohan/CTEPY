! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       IBS ROUTINE PIWINSKI USING SMOOTH PARAMETERS
! ---------------------------------------------------------------------------------------------------------------

! ===============================================================================================================
! pnumber  : number of real particles in the bunch
! epsx     : hor emit
! epsy     : ver emit
! sigs     : bunch len
! dponp    : momentum spread
! pi       : constant
! circ     : accelerator circumference
! clight   : speed of light
! qatom    : charge
! aatom    : atomic number
! betar    : relativistic beta
! betx     : BETX madx
! bety     : BETY madx
! dispx    : DX madx
! dispxp   : DPX madx
! leng     : L madx
! alfx     : ALFX madx
! nElem    : number of elements
! gamma0   : relativistic gamma
! alfap0   : longitudinal growth rate
! alfax0   : hor growth rate
! alfay0   : ver growth rate
! ===============================================================================================================

subroutine ibsPiwSmooth(pnumber,epsx,epsy,sigs,pi,gamma0,gammat,circ,betar,&
betax,betay,clight,qatom,aatom,dponp,alfap0,alfax0,alfay0)
! uses Piwinski's formulae for ibs, pg 126 in handbook
! smooth lattice approximation

double precision, intent(in) :: circ,pi,gammat, epsx, epsy, sigs, dponp, pnumber,betax,betay
double precision, intent(in) :: qatom,aatom, clight, betar,gamma0
double precision, intent(out) :: alfap0,alfax0,alfay0

double precision, external :: fmohl

double precision :: xdisp, sigh2inv, sigh
double precision :: r0,atop,abot,rmsx,rmsy,ca,a,b,d,q,fmohlp,fmohlx,fmohly
integer, PARAMETER :: npp=1000

! dispersion in smoooth approximation
      xdisp = circ/(2*pi*gammat**2)

      rmsx=sqrt(epsx*betax)
      rmsy=sqrt(epsy*betay)

      sigh2inv =  1/dponp**2  + (xdisp/rmsx)**2
      sigh = 1/sqrt(sigh2inv)

! get capitol A
!      write(6,*)sigx,sigy,dponp,beta,sigh
! classical radius
      r0 = (qatom**2/aatom)*1.54e-18
      atop = r0*r0*clight*pnumber
      abot = 64*pi*pi*betar**3*gamma0**4*epsy*epsx*sigs*dponp
! ca is Piwinski's A
      ca = atop/abot
! mohl's a,b,and q
! a is horizontal, x
      a = sigh*betax/(gamma0*rmsx)
! b is vertical
      b = sigh*betay/(gamma0*rmsy)
! log is good enough
      if (rmsx.le.rmsy) then
         d=rmsx
      else
         d=rmsy
      endif
      q = sigh*betar*sqrt(2*d/r0)
! calculate fmohl(a,b,q) with 1000 points
      fmohlp = fmohl(a,b,q,npp)
      fmohlx = fmohl(1/a,b/a,q/a,npp)
      fmohly = fmohl(1/b,a/b,q/b,npp)
      alfap0 =ca*fmohlp*(sigh/dponp)**2
      alfax0 =ca*(fmohlx+fmohlp*(xdisp*sigh/rmsx)**2)
      alfay0 =ca*fmohly
      return
      end

include 'fmohl.f95'
