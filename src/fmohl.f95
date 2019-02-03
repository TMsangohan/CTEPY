! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       FUNCTION TO CALCULATE FMOHL A HELPER FUNCTION IN THE NAGAITSEV IBS MODULE
! ---------------------------------------------------------------------------------------------------------------

double precision function fmohl(a,b,q,np)
integer, intent(in) ::np
double precision, intent(in) :: a,b,q
double precision :: du,u,cp,cq,dsum,sum

! direct crib from page 126 in handbook
      ceuler = 0.577215
      pi = 3.14159265
      sum = 0
      du = 1/float(np)
      do k=0,np
         u = k*du
         cp = sqrt(a*a+(1-a*a)*u*u)
         cq = sqrt(b*b+(1-b*b)*u*u)
         dsum = 2*log(q*(1/cp+1/cq)/2) - ceuler
         dsum = dsum*(1-3*u*u)/(cp*cq)
         if(k.eq.0)dsum=dsum/2
         if(k.eq.np)dsum=dsum/2
         sum = sum + dsum
      enddo
      fmohl = 8*pi*du*sum
!      write(6,*)'a,b,q,fmohl',a,b,q,fmohl
      return
end

