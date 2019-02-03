! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       FUNCTIONS TO KEEP VERTICAL DATA
! ---------------------------------------------------------------------------------------------------------------
subroutine keepturn(np,nb,y,py,ykeep)
implicit none    

integer :: k
integer, intent(in) ::np,nb
!f2py intent(in) :: np,nb
double precision :: yk,pk
double precision, intent(in), dimension(nb) :: y,py
!f2py intent(in) :: y,py
double precision, intent(out), dimension(7,nb):: ykeep
!f2py intent(out) :: ykeep

! accumulates the running sums
 do k=1,np
    yk = y(k)
    pk = py(k)
    ykeep(1,k) = ykeep(1,k) + yk*yk + pk*pk
 enddo
!      write(60,*)ykeep(1,1),ykeep(1,2)
return
end

