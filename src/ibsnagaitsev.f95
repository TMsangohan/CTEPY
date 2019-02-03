! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       IBS ROUTINE NAGAITSEV
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
! dx       : DX madx
! dxp      : DPX madx
! l        : L madx
! alfx     : ALFX madx
! nelem    : number of elements
! gamma0   : relativistic gamma
! alfap0   : longitudinal growth rate
! alfax0   : hor growth rate
! alfay0   : ver growth rate
! ===============================================================================================================

subroutine nagaitsevIBSlattice(pnumber,epsx,epsy,dponp,sigs,qatom,aatom,nelem,dxp,alfx,betx,dx,bety,&
gamma0,circ,clight,pi,coulomb,alfap0,alfax0,alfay0,l,betar)
!     calcualtes IBS growth rates using the Nagaitsev formulation of Bjorken-Mtingwa. Reference: PRSTAB 8, 064403 (2005)

integer :: i
double precision :: rIon

double precision, external :: rds

integer, intent(in) :: nelem
double precision, intent(in) ::  pnumber,epsx,epsy,sigs,dponp,qatom, aatom,gamma0,circ,clight,pi,coulomb,betar
double precision, intent(in), dimension(nelem) :: dxp, alfx,betx,dx,bety,l
double precision, intent(out) :: alfap0,alfax0,alfay0
double precision :: phi,axx,ayy,sigmax,sigmay,as,a1,a2,lambda1,lambda2,lambda3 
double precision :: R1,R2,R3,sp,sx,sxp,alfapp,alfaxx,alfayy,b1

rIon   = (qatom**2/aatom)*1.54e-18
alfap0 = 0.
alfax0 = 0.
alfay0 = 0.
do i=1,nelem
   phi = dxp(i)+(alfx(i)*(dx(i)/betx(i)))
   axx = betx(i)/epsx
   ayy = bety(i)/epsy

   sigmax = sqrt(dx(i)**2*dponp**2+epsx*betx(i))
   sigmay = sqrt(epsy*bety(i))

   as = axx*(dx(i)**2/betx(i)**2+phi**2)+(1/(dponp**2))
   a1 = 0.5d0*(axx+gamma0**2*as)
   a2 = 0.5d0*(axx-gamma0**2*as)
   b1 = sqrt(a2**2 + gamma0**2*axx**2*phi**2)
   
   lambda1=ayy
   lambda2 = a1 + b1
   lambda3 = a1 - b1
 
   R1 = (1/lambda1)*rds((1./lambda2),(1./lambda3),(1./lambda1))
   R2 = (1/lambda2)*rds((1./lambda3),(1./lambda1),(1./lambda2))
   R3 = 3*sqrt((lambda1*lambda2)/lambda3)-(lambda1/lambda3)*R1-(lambda2/lambda3)*R2
   
   sp  = (gamma0**2/2.) * ( 2.*R1 - R2*( 1. - 3.*a2/b1 ) - R3*( 1. + 3.*a2/b1 ))
   sx  = 0.5 * (2.*R1 - R2*(1. + 3.*a2/b1) -R3*(1. - 3.*a2/b1))
   sxp = (3.*gamma0**2*phi**2*axx)/b1*(R3-R2)
   
   alfapp = sp/(sigmax*sigmay)
   alfaxx = (betx(i)/(sigmax*sigmay)) * (sx+sxp+sp*(dx(i)**2/betx(i)**2+phi**2))
   alfayy = (bety(i)/(sigmax*sigmay)) * (-2.d0*R1+R2+R3)
   
   alfap0 = (alfapp*l(i)/circ) / dponp**2 * (pnumber*rIon**2*clight*coulomb)/(12.d0*pi*betar**3*gamma0**5*sigs)/2 + alfap0
   alfax0 = (alfaxx*l(i)/circ) / epsx *     (pnumber*rIon**2*clight*coulomb)/(12.d0*pi*betar**3*gamma0**5*sigs)/2 + alfax0
   alfay0 = (alfayy*l(i)/circ) / epsy *     (pnumber*rIon**2*clight*coulomb)/(12.d0*pi*betar**3*gamma0**5*sigs)/2 + alfay0
enddo
! extra division by 2 because later in program is multiplied with, comes from other models
! that simulate sigmas and not emittances
return
end subroutine nagaitsevIBSlattice

include 'rds.f95' 
