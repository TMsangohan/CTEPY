! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : RODERIK BRUCE
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       FUNCTIONS TO GENERATE THE EFFECT OF THE COLLIMATION SYSTEM - BETATRON AND MOMENTUM CLEANING
!       RETURNS IMPACT DISTRIBUTIONS ON COLLIMATORS
! ---------------------------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------
! np               : size of x,px,y,py,t,pt or the current number of macro particles
! x,px,y,py,t,pt   : 6D distribution of macro particles
! nLostMom         : variable that keeps track of the number of particles lost due to momentum cleaning
! NLostBeta        : variable that keeps track of the number of particles lost due to betatron cleaning
! eqTime           : time equivalent in hours of a simulation turn
! betax            : average betax in the ring
! betay            : average betay in the ring
! refEmxy          : reference emittance
! collimAvgSwitch  : switch for writing impact distributions
! collx            : betatron x impact distribution
! colly            : betatron y impact distribution
! momx             : momentum x impact distribution
! kept             : number of particle remaing after betatron and momentum cleaning
! nSigCutBeta      : number of beam sigma collimation cuts
! nSigCutMom       : number of beam sigma collimation cuts
! betaxMom         : betax at momentum cleaning
! betar            : relativistic beta
! gamma0           : relativistic gamma
! dispxMom         : dispersion at momentum cleaning
! kturns           : current turn number
! nturns           : total number of simulation turns
! ---------------------------------------------------------------------------------------------------------------
subroutine collimation(np,x,px,y,py,t,pt,nLostMom,nLostBeta,&
eqTime,betax,betay,refEmxy,collimAvgSwitch,collx,colly,momx,&
kept,nSigCutBeta,nSigCutMom,betaxMom,betar,gamma0,dispxMom,kturns,nturns)
implicit none

integer, intent(in) :: np,kturns,nturns
integer, intent(inout) :: nLostMom,nLostBeta
double precision, dimension(np,4), intent(out) :: collx,colly,momx
double precision, intent(in) :: nSigCutBeta,nSigCutMom,betaxMom,betar,gamma0,dispxMom
double precision, intent(in):: eqTime,betax,betay,refEmxy
double precision, dimension(np), intent(inout) :: x,px,y,py,t,pt
integer, intent(in) :: collimAvgSwitch
double precision :: xCutBeta,yCutBeta,xCutBeta2,yCutBeta2,deltaPoP,xCutMom
double precision :: xDx,xmax2,ymax2,betaRatio

integer :: koo,i,lost
integer, intent(out) :: kept

!      open(unit=80, file='collimator_impacts_x_betatron.out',access='append')
!      open(unit=81, file='collimator_impacts_y_betatron.out',access='append')
!      open(unit=82, file='collimator_impacts_x_momentum.out',access='append')

      ! loop through particles. For each particle, calculate betatron amplitude and momentum amplitude. loop through aperture cuts and check if outside.
      ! amplitude can be calcualted either taking into account both amplitude and phase, or amplitude only.

koo=0

! cutoff values for betatron collimation
xCutBeta  = nSigCutBeta*sqrt(betax*refEmxy) ! sufficient to use beta function of ring
xCutBeta2 = xCutBeta**2
yCutBeta  = nSigCutBeta*sqrt(betay*refEmxy)
yCutBeta2 = yCutBeta**2
betaRatio = sqrt(betaxMom/betax)

! cutoff values for momentum collimation.
xCutMom = nSigCutMom * sqrt(refEmxy*betaxMom)

do i=1,np
   lost=0
   deltaPoP = pt(i)/(betar*betar*gamma0)
   if (collimAvgSwitch.eq.1) then ! max amplitude, disregarding phase (if 1 sim turn is many machine turns)
      xmax2 = x(i)**2 + px(i)**2
      ymax2 = y(i)**2 + py(i)**2
      xDx = sqrt(xmax2)*betaRatio + abs(dispxMom*deltaPoP)
   else ! take into account phase
      xmax2 = x(i)**2
      ymax2 = y(i)**2
      xDx = abs( x(i)*betaRatio + dispxMom*deltaPoP )
   endif

   ! check hor. betatron cut (checking squares to save sqrt operation)
   if (xmax2.ge.xCutBeta2) then
      lost = 1
      nLostBeta = nLostBeta+1
      if (collimAvgSwitch.eq.0) then
         collx(i,1) = kturns
         collx(i,2) = real(kturns)/real(nturns)*eqTime
         collx(i,3) = abs(x(i))-xCutBeta
         collx(i,4) = (abs(x(i))-xCutBeta)/sqrt(refEmxy*betax)
      endif
   endif
   ! check ver. betatron cut
   if ((ymax2.ge.yCutBeta2).and.(lost.eq.0)) then ! don't double count the particle if lost in both x and y
      lost = 1
      nLostBeta = nLostBeta+1
      if (collimAvgSwitch.eq.0) then
         colly(i,1) = kturns
         colly(i,2) = real(kturns)/real(nturns)*eqTime
         colly(i,3) = abs(y(i))-yCutBeta
         colly(i,4) = (abs(y(i))-yCutBeta)/sqrt(refEmxy*betay)
      endif
   endif
   ! check momentum collimator
   if ((xDx.ge.xCutMom).and.(lost.eq.0)) then
      lost=1
      nLostMom = nLostMom+1
      if (collimAvgSwitch.eq.0) then
         momx(i,1) = kturns
         momx(i,2) = real(kturns)/real(nturns)*eqTime
         momx(i,3) = xDx-xCutMom
         momx(i,4) = (xDx-xCutMom)/sqrt(refEmxy*betaxMom)
      endif
   endif
   if (lost.eq.0) then ! keep particle
            koo = koo + 1
            t(koo)=t(i)
            pt(koo)=pt(i)
            x(koo)=x(i)
            y(koo)=y(i)
            px(koo)=px(i)
            py(koo)=py(i)
   endif

enddo
kept = koo+1 ! plus 1 because python arrays start from zero compared to fortran arrays that start at one, makes it easier in the python main program
return
end subroutine collimation
