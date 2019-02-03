! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       FUNCTIONS TO MATCH DISTRIBUTION LONGITUDINALLY WITH HAMILTONIAN
!   DEPENDENCY : 
!       -> ham.f95
! ---------------------------------------------------------------------------------------------------------------

subroutine sampleLongMatched(np,t,pt,fnharm,fnharm2,ham1sig,sigs,&
tauhat,ptmax,v00,omega0,hammax,clight,tcoeff,vrf1,vrf2,iseed)
implicit none

integer, intent(in) :: np
double precision, intent(inout) :: iseed
double precision, dimension(np), intent(out) :: t,pt
double precision, intent(out) :: sigs
double precision, intent(in) :: ham1sig,tauhat,ptmax,v00,omega0,hammax
double precision, intent(in) :: fnharm,fnharm2,clight
double precision, intent(in) :: tcoeff,vrf1,vrf2

double precision :: tk,ptk,p1,p2,test,hamm,prob
double precision, external :: ham
integer :: k

 sigs=0.

 do k=1,np
55  continue
    call random_number(iseed)
    tk  = tauhat*(2*iseed-1)
    call random_number(iseed)
    ptk = ptmax*(2*iseed-1)

    p1 = fnharm*tk*omega0   
    p2 = fnharm2*tk*omega0

    ! exact hamiltonian:
    hamm = ham(tk,ptk,tcoeff,omega0,v00,fnharm,fnharm2,vrf1,vrf2)     

    if(hamm>hammax)go to 55 ! restart sampling if outside bucket

    prob = exp(-hamm/ham1sig)
    call random_number(iseed)
    test = iseed
    
    if(prob.lt.test) go to 55
    
    sigs  = sigs + tk**2
    t(k)  = tk
    pt(k) = ptk
 enddo
  
 sigs = clight * sqrt(sigs/np)
end subroutine sampleLongMatched

include 'ham.f95'
