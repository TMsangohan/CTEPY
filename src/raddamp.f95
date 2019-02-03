! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       ROUTINE FOR RADIATION DAMPING
! ---------------------------------------------------------------------------------------------------------------

! ==============================================================================================================
! tradlong : longitudinal radiation damping lifetime
! tradperp : transverse radiation damping lifetime
! trev     : revolution time
! siglong  : bunch len
! sigperp  : transverse beam size
! tratio   : ratio between simulation and real time
! iseed    : random number, needs to be intent(inout) because its value is updated by the random generator
! np       : number of macro particles
! ==============================================================================================================
subroutine raddamp(x,px,y,py,pt,tradlong,tradperp,trev,siglong,sigperp,tratio,iseed,np)
! does radiation damping and quantum excitation once per turn

implicit none
integer, intent(inout) :: np
!f2py intent(in,out) :: np
double precision, intent(inout):: iseed
!f2py intent(in,out) :: iseed
integer ::k
double precision, intent(in) :: tratio,trev,tradlong,tradperp,siglong,sigperp
!f2py intent(in) :: tratio,trev,tradlong,siglong,sigperp
double precision :: coeffdecaylong,coeffexcitelong,coeffgrow,coeffdecay
double precision, intent(inout), dimension(np) ::x,px,y,py,pt
!f2py intent(in,out) :: x,px,y,py,pt
double precision,external :: ran
      coeffdecaylong  = 1 - ((trev / tradlong) * tratio)
        
      ! excitation uses a uniform deviate on [-1:1]
      coeffexcitelong = siglong * sqrt(3.) * sqrt(2 * (trev / tradlong) * tratio)
        
      ! tradperp is the damping time for EMITTANCE, therefore need to multiply by 2
      ! assume same damping in horizontal and vertical plane (I4x,I4y<<I2)
      coeffdecay      = 1 - ((trev /(2 * tradperp)) * tratio)
        
      ! exact     coeffgrow= sigperp*sqrt(3.)*sqrt(1-coeffdecay**2)
      ! but trev << tradperp so
      coeffgrow       = sigperp * sqrt(3.) * sqrt(2 * (trev /(2 * tradperp)) * tratio)


! skip if transverse damping time is not positive
      if(tradperp.le.0)return

      do k=1,np
! longitudinal
         call random_number(iseed)
         pt(k) = pt(k)*coeffdecaylong +coeffexcitelong*(2*iseed-1)
! transverse
         x(k)  = coeffdecay*x(k) + (2*iseed-1)*coeffgrow
         px(k) = coeffdecay*px(k)+(2*iseed-1)*coeffgrow
         y(k)  = coeffdecay*y(k) + (2*iseed-1)*coeffgrow
         py(k) = coeffdecay*py(k)+(2*iseed-1)*coeffgrow
      enddo
      return
end

