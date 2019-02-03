! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       FUNCTIONS TO CALCULATE NAGAITSEV FUNCTION  RDS
! ---------------------------------------------------------------------------------------------------------------


double precision function rds(x,y,z)
!implicit none

integer :: iter
!double precision, intent(out) :: rds
double precision :: errtol,tiny,big,c1,c2,c3,c4,c5,c6
double precision :: xt,yt,zt,sum,fac
double precision :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,sqrtx,sqrty,sqrtz
double precision, intent(in) :: x,y,z
!f2py intent(in) ::  x, y,z

! init
errtol = 0.05
tiny   = 1.0e-25
big    = 4.5e21 
c1     = 3.0/14.0
c2     = 1.0/6.0
c3     = 9.0/22.0
c4     = 3.0/26.0
c5     = 0.25*c3
c6     = 1.5*c4

xt  = x
yt  = y
zt  = z
sum = 0.0
fac = 1.0
iter= 0

do
   iter=iter+1
   sqrtx=sqrt(xt)
   sqrty=sqrt(yt)
   sqrtz=sqrt(zt)
   alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
   sum=sum+fac/(sqrtz*(zt+alamb))
   fac=0.25*fac
   xt=0.25*(xt+alamb)
   yt=0.25*(yt+alamb)
   zt=0.25*(zt+alamb)
   ave=0.2*(xt+yt+3.0*zt)
   delx=(ave-xt)/ave
   dely=(ave-yt)/ave
   delz=(ave-zt)/ave
   if (max(abs(delx),abs(dely),abs(delz)) <= ERRTOL) exit
end do
ea=delx*dely
eb=delz*delz
ec=ea-eb
ed=ea-6.0*eb
ee=ed+ec+ec
rds=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
return
end
