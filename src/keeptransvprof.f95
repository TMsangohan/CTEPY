! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       FUNCTIONS TO DETERMINE TRANSVERSE PROFILE AND STORE
! ---------------------------------------------------------------------------------------------------------------

subroutine keepTransvProf(avglinex,avgliney,x,y,emix,emiy,betax,betay,nBins,np)
! bin the transverse coord., store the  profile in arrays
implicit none
integer :: k,nk
double precision :: xmn,dx,ymn,dy

integer, intent(in) :: np,nBins
double precision, intent(in), dimension(np) :: x,y
!f2py intent(in) :: x,y
double precision, intent(in) :: emix,emiy,betax,betay
double precision, intent(inout), dimension(nBins) :: avglinex,avgliney
!f2py intent(in,out) :: avglinex,avgliney
      xmn=-5*sqrt(betax*emix)
      ymn=-5*sqrt(betay*emiy)
      dx =-2*xmn/nBins
      dy =-2*ymn/nBins

      do k=1,np
! use nearest integer to find bin in x
         nk = nint((x(k)-xmn)/dx)
         if(nk.le.1)nk=1
         if(nk.gt.nBins)nk=nBins
         avglinex(nk)=avglinex(nk)+1

! same in y
         nk = nint((y(k)-ymn)/dy) ! using DIFFERENT scale in x and y
         if(nk.le.1)nk=1
         if(nk.gt.nBins)nk=nBins
         avgliney(nk)=avgliney(nk)+1
      enddo

      return
end
