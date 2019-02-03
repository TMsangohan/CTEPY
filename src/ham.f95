! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       FUNCTIONS TO GENERATE PARTICLE DISTRIBUTIONS FOR THE SIMULATION
! ---------------------------------------------------------------------------------------------------------------

! ================================================================================================================
! t        : UNIFORMLY DISTRIBUTED WITH SIGMA 2 TAUHAT
! pt       : UNIFORMLY DISTRIBUTED WITH SIGMA 2 PTMAX
! tauhat   : HALF LENGTH OF THE BUCKET
! p1       : UNIFORMLY DISTRIBUTED PHASE RF SYSTEM 1
! p2       : UNIFORMLY DISTRIBUTED PHASE RF SYSTEM 2
! v00      : VOLTAGE AT WHICH THE RF SYSTEM CROSSES ZERO
! omega0   : INPUT ANGULAR FREQUENCY
! ham      : THE CALCULATED HAMILTONIAN
! tcoeff   : TIME COEFFICIENT tcoeff = trev*eta/(betar**2*gamma)
! fnharm   : HARMONIC NUMBER OF THE FIRST RF SYSTEM
! fnharm2  : HARMONIC NUMBER OF THE SECOND RF SYSTEM
! vrf1     : RF VOLTAGE OF FIRST SYSTEM
! vrf2     : RF VOLTAGE OF SECOND SYSTEM
! ------------------------------------------------------------------------------------------------------------
!   DESCRIPTION :
!         CALCULATES THE LONGITUDINAL HAMILTONIAN FOR A PARTICLE WITH GIVEN INPUT PARAMETERS IN A DOUBLE RF
!         SYSTEM.
! ==============================================================================================================
double precision function ham(t,pt,tcoeff,omega0,v00,fnharm,fnharm2,vrf1,vrf2)
   implicit none

   double precision, intent(in) :: t,pt,tcoeff,omega0,v00,fnharm,fnharm2,vrf1,vrf2
   double precision :: p1,p2

   p1 = fnharm*t*omega0
   p2 = fnharm2*t*omega0
   ham = 0.5*pt*pt*tcoeff+(cos(p1)-1)*vrf1/(fnharm*v00*omega0) + (cos(p2)-1)*vrf2/(fnharm2*v00*omega0)
   return
end function ham


! ==============================================================================================================
! tauhat   : HALF LENGTH IN TIME OF THE BUCKET
! v00      : VOLTAGE AT WHICH THE RF SYSTEM CROSSES ZERO
! omega0   : INPUT ANGULAR FREQUENCY
! hammax   : THE CALCULATED MAX HAMILTONIAN
! fnharm   : HARMONIC NUMBER OF THE FIRST RF SYSTEM
! fnharm2  : HARMONIC NUMBER OF THE SECOND RF SYSTEM
! vrf1     : RF VOLTAGE OF FIRST SYSTEM
! vrf2     : RF VOLTAGE OF SECOND SYSTEM
! ------------------------------------------------------------------------------------------------------------
!   DESCRIPTION :
!         CALCULATES THE MAX LONGITUDINAL HAMILTONIAN FOR A PARTICLE WITH GIVEN INPUT PARAMETERS IN A DOUBLE RF
!         SYSTEM.
! ==============================================================================================================
double precision function hammax(tauhat,omega0,v00,fnharm,fnharm2,vrf1,vrf2)
   implicit none

   double precision, intent(in) :: tauhat,omega0,v00,fnharm,fnharm2,vrf1,vrf2
   double precision :: p1,p2

   p1 = fnharm*tauhat*omega0
   p2 = fnharm2*tauhat*omega0
   hammax  = (cos(p1)-1)*vrf1/(fnharm*v00*omega0) + (cos(p2)-1)*vrf2/(fnharm2*v00*omega0)
   return
end function hammax

! ==============================================================================================================
! pi       : THE REAL NUMBER PI
! tauhat   : HALF LENGTH IN TIME OF THE BUCKET
! v00      : VOLTAGE AT WHICH THE RF SYSTEM CROSSES ZERO
! omega0   : INPUT ANGULAR FREQUENCY
! hammax   : THE CALCULATED MAX HAMILTONIAN
! fnharm   : HARMONIC NUMBER OF THE FIRST RF SYSTEM
! fnharm2  : HARMONIC NUMBER OF THE SECOND RF SYSTEM
! vrf1     : RF VOLTAGE OF FIRST SYSTEM
! vrf2     : RF VOLTAGE OF SECOND SYSTEM
! ------------------------------------------------------------------------------------------------------------
!   DESCRIPTION :
!         CALCULATES THE LONGITUDINAL HAMILTONIAN FOR A PARTICLE WITH GIVEN INPUT PARAMETERS IN A DOUBLE RF
!         SYSTEM WHEN RF SYSTEM ARE AT ZERO CROSSING PHASE
! ==============================================================================================================
double precision function hamzero(phizero,omega0,v00,fnharm,fnharm2,vrf1,vrf2)
   implicit none

   double precision, intent(in) :: phizero,omega0,v00,fnharm,fnharm2,vrf1,vrf2
   double precision :: p2

   p2      = phizero*fnharm2/fnharm
   hamzero = (cos(phizero)-1)*vrf1/(fnharm*v00*omega0) + (cos(p2)-1)*vrf2/(fnharm2*v00*omega0)
   return
end function hamzero

