! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       FUNCTIONS TO GET TRANSVERSE EMITTANCES FROM DISTRIBUTIONS
! ---------------------------------------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------------------------------------
! RETURNS TRANSVERSE EMITTANCES
! -----------------------------
! x1         = array containing the horizontal position distribution
! px1        = array containing the horizontal momentum distribution
! y1         = array containing the vertical position distribution
! py1        = array containing the vertical momentum distribution
! np1        = size of these arrays
! ex1        = returned horizontal emittance
! ey1        = returned vertical emittance
! betax      = average horizontal optical beta of the ring
! betay      = average vertical optical beta of the ring
! xcut       = cutting the horizontal distribution off at a set number of sigma
! ycut       = cutting the vertical distribution off at a set number of sigma 
! cutoffAmpl = the number of sigma to cut at (xcut = cutoffamp * sqrt(betax * refex))
! emitMethod = method to use to calculate the emittances
! refEmxy    = reference emittance to use for beam size calculation for the cut
! ---------------------------------------------------------------------------------------------------------------
!*********************************************************************************
! schema
! ------
! getemittance <- emitStDev
!              <- calcEmitWCut
!              <- fitemittance <- getNewtonSigma
!*********************************************************************************

subroutine getemittance(x1,y1,px1,py1,np1,ex1,ey1,betax,betay,xcut,ycut,cutoffAmpl,emitMethod,refEmxy)
implicit none

character(len=5) :: emitMethod
integer, intent(in) :: np1
double precision, intent(in) :: betax, betay, xcut,ycut,cutoffAmpl,refEmxy
double precision, intent(in), dimension(np1) :: x1,px1,y1,py1
double precision, intent(out) :: ex1,ey1

  if (emitMethod.eq."stdev") then
     call emitStDev(x1,px1,np1,betax,ex1)
     call emitStDev(y1,py1,np1,betay,ey1)
  elseif (emitMethod.eq."stcal") then
     call calcEmitWCut(x1,px1,np1,xcut,betax,ex1)
     call calcEmitWCut(y1,py1,np1,ycut,betay,ey1)
  elseif (emitMethod.eq."exfit") then
     call fitemittance(x1,px1,np1,betax,cutoffAmpl,ex1,refEmxy)
     call fitemittance(y1,py1,np1,betay,cutoffAmpl,ey1,refEmxy)
  else
     write(6,*) 'unknown emittance routine. Stop.'
     stop
  endif
  ! write(6,*) ex1
  ! write(6,*) 'emittance routine = ',emitMethod
  return
end subroutine getemittance

!*********************************************************************************
!*********************************************************************************

subroutine emitStDev(x,px,np,bx,emit)
integer :: i
integer, intent(in) :: np
double precision, intent(in) :: bx
double precision, intent(in), dimension(np) :: x, px
double precision, intent(out) :: emit
       
  emit=0.0
  do i=1,np
     emit=emit+x(i)**2+px(i)**2
  enddo
  emit=emit/bx/np/2.
  return
end 

!*********************************************************************************
!*********************************************************************************

subroutine fitEmittance(x,px,np,bx,Jmax,emit,refEmxy)
implicit none

integer :: i,bin,maxbin
integer, intent(in) :: np
double precision, intent(in), dimension(np) :: x, px
double precision, intent(in) :: Jmax, bx, refEmxy
double precision, intent(out) :: emit
double precision, dimension(1000) :: binVals
double precision :: Jx,sigma1,sigma0,getNewtonSigma,bin1,binwidth
! calculate the betatron action for each particle. 
! Then bin the action and fit the resulting bin values to an exponential. 
! Not as robust as calculating backwards from standard deviation as in StdDev routine

  maxbin=1000
  binwidth=0.025 ! fraction of a sigma used as bin width

  do i=1,maxbin
     binVals(i)=0.0
  enddo

  do i=1,np
     Jx=(x(i)**2+px(i)**2)/(2*bx*refEmxy)      ! calculate the normalized action.
     bin=int(Jx/binwidth)+1
     if (bin<=maxbin) binVals(bin)=binVals(bin)+1.
  enddo

  bin1=binVals(1)

  do i=1,maxbin ! for better numerical accuracy, normalize everything to first bin
     binVals(i)=binVals(i)/bin1
  enddo

! solve f'=0 with Newton's method -- see function getNewtonSigma for details
  sigma0 = 1.

  do i=1,100
     sigma1=getNewtonSigma(sigma0,binVals,Jmax,binwidth)
!     write(6,*) sigma0,sigma1
     if (abs(sigma1-sigma0).le.0.01) exit
     sigma0=sigma1
  enddo
! the thing that is returned by this mess
  emit = sigma1 * refEmxy

  return
end 

!*********************************************************************************
!*********************************************************************************

double precision function getNewtonSigma(sigmaOld,binVals,Jmax,binwidth)
implicit none
! for fitting sigma in the function g(x)=binVals(1)*Exp[-x/sigma] with Newton's method to the binned betatron actions, in order to calculate the emittance with a least square fit
! The sum of the squares is f(sigma)=sum_i[ (g(xi)-binCount(i))^2 ]. Find sigma that minimizes f
! a minimum is found for f'(sigma)=0
! with newton's method, we thus iterate sigma_(n+1) = sigma_n - f'(sigma_n)/f''(sigma_n)
! the fit only takes bins up to Jmax into account. Therefore we will get the correct emittance even if the tail of the Gaussian is cut.
integer :: maxBin,i
double precision, dimension(500) :: binVals
double precision :: fprime,fbis,A,AExpXiSig,xi,yi
double precision, intent(in)  :: sigmaOld,Jmax,binwidth
  maxbin=nint(Jmax/binwidth)
  fprime=0.
  fbis=0.
  A=binVals(1) ! amplitude of f

  do i=1,maxBin ! sum over all bins to calculate f'(sigmaOld) and f''(sigmaOld)
     xi = DBLE(i)*binwidth-binwidth/2. ! get the center of the bin
     yi = binVals(i)
     AExpXiSig = A*Exp(-xi/sigmaOld)
     fprime=fprime + 2*AExpXiSig*xi*( AExpXiSig-yi )/sigmaOld**2
     fbis=fbis + (2*A*xi*AExpXiSig**2*(-(1/AExpXiSig*yi*(xi - 2*sigmaOld)) + 2*A*(xi - sigmaOld)))/sigmaOld**4
  enddo

! output
  getNewtonSigma = sigmaOld - fprime/fbis
  return
end 

!*********************************************************************************
!*********************************************************************************

subroutine calcEmitWCut(x,px,np,cutx,bx,emit)
! routine to calculate the emittance if the distribution is Gaussian with cut tails at cutx
! just taking the standard deviation then gives a too small value
! calculate backwards from the standard deviation with a cut (see cut_Gaussian_tails_fit_emittance.nb for derivation)
integer, intent(in) :: np
integer :: i
double precision, intent(in), dimension(np) :: x, px
double precision, intent(in) :: bx, cutx
double precision, intent(out) :: emit
double precision :: sigObtained,sigma0,sigma1,func,funcP

  sigObtained=0.0

  do i=1,np
     sigObtained=sigObtained+x(i)**2+px(i)**2
  enddo

  sigObtained=sqrt(sigObtained/np/2.0) ! divide by 2 since we are including both x and px in the sum. this is now raw sigma^2

! now solve for original sigma without tails with Newton's method
! use obtained sigma as starting value

  sigma0=sigObtained
  sigma1=-1.0

  do i=1,100 ! max 100 iterations
     func = cutx**2/(2-2*Exp(cutx**2/(2*sigma0**2))) + sigma0**2 - sigObtained**2
     funcP = 2*sigma0 - cutx**4 / ( 8*sigma0**3 *sinh(cutx**2/(4*sigma0**2)) )
     sigma1=sigma0-func/funcP
     if (abs((sigma1-sigma0)/sigma0).le.0.005) exit
     sigma0=sigma1
  enddo

  if (i.eq.100) then
     write(6,*) 'emittance calculation not converging - STOP'
     stop
  endif

  emit=sigma1**2/bx

  ! write(6,*) 'emit = ',emit
  ! write(6,*) '**********************'
  ! write(6,*)
  return
end
