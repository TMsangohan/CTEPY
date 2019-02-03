! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       IBS from bane method
! ---------------------------------------------------------------------------------------------------------------

! ==============================================================================================================
! nparti    :
! ex        : horizontal emittance
! ey        : vertical emittance
! sigs      : bunch length
! sigp      : momentum spread
! qatom     : particle charge
! aatom     : atomic number 
! gamma0    : relativistic gamma
! circ      : accelerator circumference
! clight    : speed of light
! betax     : horizontal beta
! betay     : vertical  beta
! dispx     : horizontal dispersion
! dispy     : vertical dispersion
! dxp       : madx DPX
! dyp       : madx DPY
! alfx      : madx ALFX
! alfy      : madx ALFY
! leng      : element length
! alfap0    : growth rate longitudinal
! alfax0    : growth rate horizontal
! alfay0    : growth rate vertical
! nElem     : number of element (madx)
! gtab      : g number table needed for the Bane algorithm
! nx        : dimension for gtab
! ny        : dimension for gtab
! ==============================================================================================================
subroutine ibsBane(nparti,ex,ey,sigs,sigp,qatom,aatom,gamma0,circ,clight,&
betx,bety,dispx,dispy,dxp,dyp,alfx,alfy,leng,alfap0,alfax0,alfay0,nElem,gtab,nx,ny)
!  high-energy IBS approximation. Reference: K.L.F Bane, A Simplified Model of Intrabeam Scattering, SLAC-PUB-9226 (2002)

integer, intent(in) :: nElem,nx,ny
double precision, dimension(nx,ny), intent(in) :: gtab
double precision, dimension(nElem) :: betx,bety,dispx,dispy,alfx,alfy,leng,dxp,dyp
double precision, intent(in) :: nparti,ex,ey,sigs,sigp
double precision, intent(out) :: alfap0,alfax0,alfay0
double precision, intent(in) :: qatom,aatom,gamma0,circ,clight
double precision, external :: gbane
double precision :: r0,d,sigh
double precision :: a,b,bx,by,dxx,dyy,Hx,Hy
double precision :: brAVG,HxAVG,HyAVG,betxp,betyp,sigx,sigy,coulLog,Tp,Tx,Ty
integer, PARAMETER :: npp=1000

!     classical radius
      r0 = (qatom**2/aatom)*1.54e-18

!     zero out average quantities
      brAVG=0  ! br is the quantity in brackets < > to be averaged in eq. 11
      HxAVG=0
      HyAVG=0

!     loop over all elements in lattice, calculate average quantities required by eqs 11,14 in slacpub9226

      do i=1,nElem
         bx=betx(i)
         by=bety(i)
         dxx=dispx(i)
         dyy=dispy(i)
         betxp=-2*alfx(i)
         betyp=-2*alfy(i)

         sigx=sqrt(ex*bx+sigp**2*dxx)
         sigy=sqrt(ey*by+sigp**2*dyy)

         Hx=(dxx**2+(bx*dxp(i)-0.5*betxp*dxx)**2)/bx
         Hy=(dyy**2+(by*dyp(i)-0.5*betyp*dyy)**2)/by

         sigh=1./sqrt(1./sigp**2+Hx/ex+Hy/ey)
         a=sigh/gamma0*sqrt(bx/ex)
         b=sigh/gamma0*sqrt(by/ey)

         if (sigx<sigy) then
            d=sigx
         else
            d=sigy
         endif
         coulLog=log(d*sigh**2/(4*r0*a**2))

!     calcualte average over lattice of quantities required
         brAVG=brAVG+sigh*gBane(a/b,gtab,nx,ny)*(bx*by)**(-0.25)*coulLog*leng(i)  ! the factor to be averaged over in Eq.11, slacpub9226
         HxAVG=HxAVG+Hx*leng(i)
         HyAVG=HyAVG+Hy*leng(i)
      enddo

      brAVG=brAVG/circ
      HxAVG=HxAVG/circ
      HyAVG=HyAVG/circ

      Tp=1/(r0**2*clight*nparti/(16*gamma0**3*(ex*ey)**0.75*sigs*sigp**3)*brAVG)

      Tx=Tp*ex/(HxAvg*sigp**2)
      Ty=Tp*ey/(HyAvg*sigp**2)

      alfap0=1/Tp
      alfax0=1/Tx
      alfay0=1/Ty

      return
      end

!*********************************************************************************
double precision  function gBane(alpha,gTab,nx,ny)
!     calculates value of Bane's g-function through linear interpolation of tabulated values in file
integer :: nLow
integer, intent(in) :: nx,ny
double precision, intent(in) :: alpha
double precision, dimension(nx,ny),intent(in) :: gTab

      if((alpha<0.02).or.(alpha>=15.0)) then
         write(*,*) 'gBane: called with argument ',alpha,' outside range. allowed range is 0.02<arg<10'//CHAR(13)//CHAR(10)
         stop
      endif
!     find the closest lower tabulated value for alpha--nLow is the bin in the list where interpolation starts
      nLow=int(100*alpha+1)
      gBane=(gTab(nLow+1,2)-gTab(nLow,2))*(alpha-gTab(nLow,1))/(gTab(nLow+1,1)-gTab(nLow,1))+gTab(nLow,2)
      return
      end
