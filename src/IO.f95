! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       FUNCTIONS TO CALCULATE MODIFIED BESSELFUNCTION IO(X)
! ---------------------------------------------------------------------------------------------------------------

double precision function IO(iox)
!
!       =========================================================
!       Compute modified Bessel function I0(x)
!       modified to use real*4 instead of real*8
!       using real*8 causes many strange numerical instabilities

double precision, intent(in) :: iox
double precision :: EL,BI0,X2,R,CA,XR,A(12),B(12),pi
INTEGER :: K,K0
pi = 3.1415926535
!        DIMENSION A(12),B(12),A1(8)
!        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        X2=iox*iox
        IF (iox.EQ.0.0) THEN
           BI0=1.0
           IO = BI0
           RETURN
        ELSE IF (iox.LE.18.0) THEN
           BI0=1.0
           R=1.0
           DO 15 K=1,50
              R=0.25*R*X2/(K*K)
              BI0=BI0+R
              IF (ABS(R/BI0).LT.1.0E-15) GO TO 20
15         CONTINUE
20         CONTINUE
           R=1.0

        ELSE
           DATA A/0.125,7.03125E-2,&
                 7.32421875E-2,1.1215209960938E-1,&
                 2.2710800170898E-1,5.7250142097473E-1,&
                 1.7277275025845E0,6.0740420012735E0,&
                 2.4380529699556E01,1.1001714026925E02,&
                 5.5133589612202E02,3.0380905109224E03/
           DATA B/-0.375E0,-1.171875E-1,&
                 -1.025390625E-1,-1.4419555664063E-1,&
                 -2.7757644653320E-1,-6.7659258842468E-1,&
                 -1.9935317337513E0,-6.8839142681099E0,&
                 -2.7248827311269E01,-1.2159789187654E02,&
                 -6.0384407670507E02,-3.3022722944809E03/
           K0=12
           IF (iox.GE.35.0) K0=9
           IF (iox.GE.50.0) K0=7
           CA=EXP(iox)/SQRT(2.0E0*PI*iox)
           BI0=1.0E0
           XR=1.0E0/iox
           DO 35 K=1,K0
35            BI0=BI0+A(K)*XR**K
           BI0=CA*BI0
        ENDIF
IO=BI0
RETURN
END function IO
