! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       FUNCTIONS TO GENERATE ARTIFICIAL BLOWUP OF THE BEAM
! ---------------------------------------------------------------------------------------------------------------
! DEPENDENCIES :  
!     -> getgauss.f95
! ===============================================================================================================
! np           : size of px and py array, the current number of macro particles
! px           : array containing the px distribution, which is updated in this routine
! py           : array containing the px distribution, which is updated in this routine
! iseed        : number for random generator, needs to be intent(inout) as it is updated during call of random generator
! blowupMethod : method used to blow up the beam, in set("gauss","unifo","unSum")
! betax        : average betax of the ring
! betay        : average betay of the ring
! refEmxy      : reference emittance
! pxKickFac    : strength factor of the applied kick in the horizontal plane
! pyKickFac    : strength factor of the applied kick in the vertical plane
! ===============================================================================================================

subroutine blowup(np,px,py,iseed,blowupMethod,betax,betay,refEmxy,pxKickFac,pyKickFac)
integer, intent(in):: np
double precision, intent(inout) :: iseed
character(5), intent(in) :: blowupMethod
double precision, dimension(np), intent(inout) :: px,py
double precision, intent(in) :: pxKickFac,pyKickFac
double precision :: pxKick,pyKick
        ! give random kicks in angle - Gaussian with standard deviation pxKickFac or
        ! uniform in intervall {-pxKickFac,+pxKixkFac}.
        ! unit of pxKickFac: sigma with collimaiton emittance

        if (blowupMethod.eq."gauss") then
           do i=1,np
              call getgaussrv(iseed,pxKick,pyKick) ! Gaussian random with sigma=1/sqrt(3)
              px(i)=px(i)+pxKickFac*1.73205*pxKick*sqrt(betax*refEmxy) ! mult. by sqrt(3)=1.73205, since get2gaussrv returns gaussian with sigma=1/sqrt(3). Final standard dev. of kick is thus 
              py(i)=py(i)+pyKickFac*1.73205*pyKick*sqrt(betay*refEmxy)
            enddo
        elseif(blowupMethod.eq."unifo") then ! uniform distribution - 1 kick
           do i=1,np
              call random_number(iseed)
              pxKick = pxKickFac*sqrt(betax*refEmxy) * (2*iseed-1)   ! uniform random number between -1 and 1, scaled by sigma * pxKickFac. Thus: pxKickFac is max amplitude of total kick
              pyKick = pyKickFac*sqrt(betay*refEmxy) * (2*iseed-1)
              px(i) = px(i) + pxKick
              py(i) = py(i) + pyKick
           enddo
           elseif(blowupMethod.eq."unSum") then
              do i=1,np
                 pxKick = pxKickFac*sqrt(betax*refEmxy) * (2*(iseed+iseed+iseed+iseed)-4)   ! sum of four uniform random number between -1 and 1, scaled by sigma * pxKickFac. pxKickfac is thus amplitude in sigma for EACH ONE of the four ADTs
                 pyKick = pyKickFac*sqrt(betay*refEmxy) * (2*(iseed+iseed+iseed+iseed)-4)
                 px(i) = px(i) + pxKick
                 py(i) = py(i) + pyKick
              enddo
           else
              write(6,*) "Unknown blowup method - stop"
              stop
        endif
end subroutine blowup

include 'getgauss.f95'
