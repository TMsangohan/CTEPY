! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       IBS ROUTINE PIWINSKI USING LATTICE
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

subroutine ibsPiwLattice(pnumber,epsx,epsy,sigs,dponp,pi,circ,&
clight,qatom,aatom,betar,betx,bety,dispx,leng,nElem,gamma0,alfap0,alfax0,alfay0)
! uses Piwinski's formulae for ibs, pg 126 in handbook
! taking variations of optical functions in lattice into account

integer, intent(in) :: nElem
double precision, dimension(nElem) :: betx,bety,dispx,leng
double precision, intent(out) :: alfap0,alfax0,alfay0
double precision, intent(in) :: pnumber,epsx,epsy,sigs,dponp
double precision, intent(in) :: qatom,aatom,betar,gamma0,pi,circ,clight
double precision, external :: fmohl
double precision :: r0,atop,abot,ca,rmsx,rmsy,d,sigh2inv,sigh
double precision :: a,b,q,fmohlx,fmohly,fmohlp,bx,by,xdisp
integer, PARAMETER :: npp=1000

!     classical radius
      r0 = (qatom**2/aatom)*1.54e-18

!     get capitol A
      atop = r0*r0*clight*pnumber
      abot = 64*pi*pi*betar**3*gamma0**4*epsy*epsx*sigs*dponp

! ca is Piwinski's A
      ca = atop/abot

! zero out average quantities
      alfax0=0
      alfay0=0
      alfap0=0

      do i=1,nElem
         bx = betx(i)
         by = bety(i)
         xdisp = dispx(i)
         rmsx=sqrt(epsx*bx)
         rmsy=sqrt(epsy*by)
         if (rmsx.le.rmsy) then
            d=rmsx
         else
            d=rmsy
         endif

         sigh2inv =  1./dponp**2  + (xdisp/rmsx)**2
         sigh = 1./sqrt(sigh2inv)

! mohl's a,b,and q
! a is horizontal
         a = sigh*bx/(gamma0*rmsx)
! b is vertical
         b = sigh*by/(gamma0*rmsy)
         q = sigh*betar*sqrt(2*d/r0)
! calculate fmohl(a,b,q) with 1000 points
         fmohlp = fmohl(a,b,q,npp)
         fmohlx = fmohl(1/a,b/a,q/a,npp)
         fmohly = fmohl(1/b,a/b,q/b,npp)
         alfap0 =ca*fmohlp*(sigh/dponp)**2*leng(i)+alfap0
         alfax0 =ca*(fmohlx+fmohlp*(xdisp*sigh/rmsx)**2)*leng(i)+alfax0
         alfay0 =ca*fmohly*leng(i)+alfay0
      enddo
      alfap0=alfap0/circ
      alfax0=alfax0/circ
      alfay0=alfay0/circ
      return
      end

! already included in ibslong no need to reload
! include 'fmohl.f95'
