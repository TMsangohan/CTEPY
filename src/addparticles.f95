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
!   DEPENDENCY :
!           longmatch.f95 -> ham.f95
! ---------------------------------------------------------------------------------------------------------------
! GENERATE LONGITUDINAL DISTRIBUTION 
! REMARK : SPECIFICATION OF INTENT IS VERY IMPORTANT TO BE ABLE TO USE THE FORTRAN ROUTINES FROM PYTHON
! ===============================================================================================================
! tauhat         : half bucket length in nano seconds
! tcoeff         : time coefficient trev * eta / (betarel**2 gamma0)
! omega0         : 2 * pi / trev
! v00            : particle mass / charge * 1e9  
! fnharm         : harmonic number of first RF system
! fnharm2        : harmonic number of second RF system
! vrf1           : peak RF voltage of first RF system
! vrf2           : peak RF voltage of second RF system
! power          : parameter that can be set to change longitudinal distribution in smoke ring distribution method
! alint          : cfr power
! iseed          : input number of the random generator, note that this needs to be intent(inout) as its value is updated by the random routine call
! gamma0         : relativistic gamma of the particle
! rmsdelta       : delta p over p
! rmsbunchlen    : desired bunchlength of the distribution
! bunchlenprecis : desired precision for the distribution bunchlength to match with the desired one
! clight         : light speed - used as input to avoid as global variable which allows to split fortran code in multiple files
! nmacro         : number of macro particles
! t              : returned array containing the distribution in t
! pt             : returned array containing the distribution in pt
! method         : integer in set(1,2,3) selecting the method to generate the longitudinal distribution
! ===============================================================================================================
subroutine addlongitudinal(tauhat,tcoeff,omega0,v00,fnharm,fnharm2,vrf1,vrf2,pi,power,alint,iseed, &
gamma0,rmsdelta,rmsbunchlen,bunchlenprecis,clight,nmacro,t,pt,method)
implicit none
integer, intent(in) :: nmacro,method
double precision, intent(inout) :: iseed
double precision, intent(in) :: tauhat,tcoeff,omega0,v00,fnharm
double precision, intent(in) :: fnharm2,vrf1,vrf2,power,alint,bunchlenprecis
double precision, intent(in) :: gamma0,rmsdelta,rmsbunchlen,clight,pi
double precision, dimension(nmacro), intent(out) :: t, pt
double precision :: phizero,vzero,phik,vk,phik0,ptmax
double precision :: hammaxx,p1,p2,tk,ptk,hamm,prob,test,ampt,ampPt,amp
double precision, external :: hammax,ham,hamzero
double precision :: ham1sig0,sigs0,sigs1,ham1sig1,fprime,r1,r2,facc,hamzeroo

integer :: k,np
   phizero = fnharm/fnharm2*3*pi /4.
   vzero   = abs(vrf1)
   ! looking for zero crossing using binning for phizero
   do k=-1000,1000
     phik = fnharm/fnharm2*(pi+k*pi/2000)
     vk   = vrf1*sin(phik) + vrf2*sin(fnharm2*phik/fnharm)
     if (abs(vk).lt.vzero) then
       vzero = abs(vk)
       phizero = phik
     endif
   enddo

   phik0=phizero

   do k=-1000,1000
     phik = phik0+ k/100000.
     vk = vrf1*sin(phik)+vrf2*sin(fnharm2*phik/fnharm)
     if (abs(vk).lt.vzero) then
       vzero = abs(vk)
       phizero=phik
     endif
   enddo

   
   hammaxx  = hammax(tauhat,omega0,v00,fnharm,fnharm2,vrf1,vrf2)
   ptmax   = sqrt(2 * hammaxx/tcoeff)
   hamzeroo = hamzero(phizero,omega0,v00,fnharm,fnharm2,vrf1,vrf2)
   if (method==0) then ! smoke ring dist.
          do np = 1,nmacro
 43          continue
             call random_number(iseed)
             tk  = tauhat*(2*iseed-1)
             call random_number(iseed)
             ptk = ptmax*(2*iseed-1)

             p1 = fnharm*tk*omega0
             p2 = fnharm2*tk*omega0
             hamm = ham(tk,ptk,tcoeff,omega0,v00,fnharm,fnharm2,vrf1,vrf2)
             if (hamm>=hammaxx) go to 43
             prob = 1-(hamm/hammaxx)
             prob = prob**power
             if((abs(p1)>phizero).and.(hamm<=hamzeroo))then
                prob = prob*(hamm/hamzeroo)**alint
             endif
             
             call random_number(iseed)
             test = iseed
             if(prob<test)go to 43
             pt(np) = ptk
             t(np)=tk
          enddo
   ! read from file method is done in main python program
   elseif (method==2) then 
      ! generate bi-gaussian (possibly unmatched, can only be matched in small angle approximation)
      do np=1,nmacro
56      continue
        ampPt = gamma0 * rmsdelta
        ampt  = rmsbunchlen/clight
        ! r1 and r2 are uniform on (-1,1)
        call random_number(iseed)
        r1 = 2*iseed-1
        call random_number(iseed)
        r2 = 2*iseed-1
        amp = r1*r1 + r2*r2
      
        if(amp.ge.1) go to 56
      
        facc = sqrt(-2.*log(amp)/amp) ! amp is gaussian
        tk   = ampt*r1*facc
        ptk  = ampPt*r2*facc
        p1   = fnharm*tk*omega0
        p2   = fnharm2*tk*omega0
        hamm = ham(tk,ptk,tcoeff,omega0,v00,fnharm,fnharm2,vrf1,vrf2)

        if(hamm>=hammaxx)go to 56 ! restart sampling if outside bucket
        if (abs(tk).ge.tauhat) go to 56 ! restart sampling if t is too large
      
        t(np)=tk
        pt(np)=ptk
    enddo

   elseif (method==3) then ! generate matched 'pseudo-Gaussian' phase space
      ! Assuming distribution function for hamiltonian rho(h)=1/H*exp(-h/H), so rho(pt,t)~exp(-C1*pt^2-C2*(1-cos(C3*t)). at 2nd order in t this is a bi-gaussian
      ! average H can not be calculated analytically for the exact hamiltonian, only for small angle approximation.
      ! Use first this value to calculate bunch length, then use Newton's method to calculate new H and iterate until convergence
      ! get average H in small osc. approx. see derivation in check_phase_space_matching.nb
      ham1sig0 = -fnharm2*omega0*vrf2/(clight**2 * v00) * rmsbunchlen**2
      
      call samplelongmatched(nmacro,t,pt,fnharm,fnharm2,ham1sig0,sigs0,&
tauhat,ptmax,v00,omega0,hammaxx,clight,tcoeff,vrf1,vrf2,iseed)
      
      ham1sig1 = ham1sig0 * (rmsbunchlen/sigs0)**2

      call samplelongmatched(nmacro,t,pt,fnharm,fnharm2,ham1sig1,sigs1,&
tauhat,ptmax,v00,omega0,hammaxx,clight,tcoeff,vrf1,vrf2,iseed) 
      do while (abs(sigs1/rmsbunchlen-1.)>=bunchlenprecis)
             fprime   = (sigs1-sigs0) / (ham1sig1-ham1sig0)
             ham1sig0 = ham1sig1
             sigs0    = sigs1
             ham1sig1 = ham1sig0 - (sigs0-rmsBunchLen)/fPrime
             call samplelongmatched(nmacro,t,pt,fnharm,fnharm2,ham1sig1,sigs1,&
tauhat,ptmax,v00,omega0,hammaxx,clight,tcoeff,vrf1,vrf2,iseed)
      enddo
   endif     
   !sigv = sigv + pt(np)**2
end subroutine addlongitudinal

include 'longmatch.f95'

! GENERATE TRANSVERSE DISTRIBUTION 
! REMARK : SPECIFICATION OF INTENT IS VERY IMPORTANT TO BE ABLE TO USE THE FORTRAN ROUTINES FROM PYTHON
! ================================================================================================================
! iseed          : input number of the random generator, note that this needs to be intent(inout) as its value is updated by the random routine call
! nmacro         : number of macro particles
! x              : returned array containing the distribution in x
! px             : returned array containing the distribution in px
! y              : returned array containing the distribution in y
! py             : returned array containing the distribution in py
! xcut           : cut off amplitude in the horizontal plane from collimation
! ycut           : see xcut but in vertical plane
! betax          : average betax in the ring 
! betay          : average betay in the ring
! emix           : horizontal geometric emittance
! emiy           : vertical geometric emittance
! ===============================================================================================================
subroutine addtransverse(nmacro,x,px,y,py,xcut,ycut,betax,betay,emix,emiy,iseed)
implicit none

integer, intent(in) :: nmacro
double precision, intent(inout) :: iseed
double precision, intent(in) :: xcut,ycut
double precision, dimension(nmacro), intent(out) :: x,px,y,py
double precision, intent(in) :: betax,betay,emix,emiy

double precision :: ampx,ampy,r1,r2,facc,amp
integer :: k

  ampx = sqrt(betax*emix)
  ampy = sqrt(betay*emiy)
  
  do k = 1,nMacro
44   continue
     ! r1 and r2 are uniform on (-1,1)
     call random_number(iseed)
     r1 = 2*iseed-1
     call random_number(iseed)
     r2 = 2*iseed-1
     amp = r1*r1 + r2*r2
     if((amp>=1).or.(amp<= 3.e-6))go to 44
     facc = sqrt(-2.*log(amp)/amp)
     x(k) = ampx*r1*facc
     ! px has same amplitude as x ie px = beta_L*x'
     px(k) = ampx*r2*facc
     ! reject if x>initial cut off (given in sigmas with the reference emittance)
     if (sqrt(x(k)**2+px(k)**2)>=xcut) go to 44
     ! other transverse variable
45   continue
     call random_number(iseed)
     r1 = 2*iseed-1
     call random_number(iseed)
     r2 = 2*iseed-1
     amp = r1*r1 + r2*r2
     if((amp>=1).or.(amp<= 3.e-6))go to 45
     facc = sqrt(-2.*log(amp)/amp)
     y(k)= ampy*r1*facc
     py(k)= ampy*r2*facc
     ! reject if x>initial cut off (given in sigmas with the reference emittance)
     if (sqrt(y(k)**2+py(k)**2)>=ycut) go to 45
  enddo
end subroutine addtransverse
