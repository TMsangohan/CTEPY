! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       ROUTINE TO TAKE INTO ACCOUNT CHROMATICITY 
! ---------------------------------------------------------------------------------------------------------------

! ==============================================================================================================
! psix        : MUX madx
! psiy        : MUY madx
! coeffchromx : horizontal chromaticity (input)
! coeffchromy : vertical chromaticity (input)
! tpdqmin     :
! ==============================================================================================================

subroutine sptransChrom(np,x,px,y,py,pt,psix,psiy,coeffchromx,coeffchromy,tpdqmin,k2L,k2Lskew)
! full turn update for x, horizontal
implicit none
integer, intent(in) :: np
double precision, intent(inout), dimension(np) :: x,px,y,py,pt
!f2py intent(in,out) :: x,px,y,py,pt
double precision, intent(in) :: psix,psiy,coeffchromx,coeffchromy,k2L,k2Lskew,tpdqmin
double precision :: xk,pkx,xk1,pxk1,yk,pky,yk1,pyk1,psi1,ptk,a11,a12
integer :: k

      do k = 1,np

! rotate in y-py plane
         ptk  = pt(k)
         psi1 = psiy + ptk * coeffchromy
         a11  = cos(psi1)
         a12  = sin(psi1)
         yk   = y(k)
         pky  = py(k)
         yk1  = yk * a11 + pky * a12
         pyk1 = pky * a11 - yk * a12
         y(k) = yk1

! rotate in x-px plane
         psi1 = psix + ptk*coeffchromx
         a11  = cos(psi1)
         a12  = sin(psi1)
         xk   = x(k)
         pkx  = px(k)
         xk1  = xk * a11 + pkx * a12
         pxk1 = pkx * a11 - xk * a12
         x(k) = xk1

! now have dqmin part - coupling between x and y
         px(k) = pxk1 + tpdqmin * yk1
         py(k) = pyk1 + tpdqmin * xk1

! thin sextupole kick 
         px(k)= px(k) + 0.5 * k2L * (xk1**2 - yk1**2) - k2Lskew * (xk1 * yk1)
         py(k)= py(k) + k2L * (xk1 * yk1) + 0.5 * k2Lskew * (xk1**2 - yk1**2)

      enddo
!      write(44,*)psi1,ptk*coeffchrom,ptk
      return
end
