! ---------------------------------------------------------------------------------------------------------------
! ORIGINAL AUTHORS : MIKE BLASKIEWISC, RODERIK BRUCE, MICHAELA SCHAUMANN, TOM MERTENS
! ---------------------------------------------------------------------------------------------------------------
! VERSION 2.0 : UPDATE TO TRACK MULTIPLE BUNCHES AND HAVE DIFFERENT PARTICLE TYPES (P-PB) 
!   AUTHOR    : TOM MERTENS
!   DATE      : 19/12/2016
!   COPYRIGHT : CERN
!   
!   DESCRIPTION : 
!       IBS from interpolation, input table needed to interpolate over
! ---------------------------------------------------------------------------------------------------------------

! ===============================================================================================================
! pnumber   : number of real particles in bunch
! epsx      : horizontal emittance
! epsy      : vertical emittance
! sigs      : bunch length
! dponp     : momentum spread
! xmin      : x min in table
! ymin      : y min in table
! zmin      : z min in table
! xmax      : x max in table
! ymax      : y max in table
! zmax      : z max in table
! xBins     : number of bins horizontal
! yBins     : number of bins vertical
! zBins     : number of bins longitudinal
! Ap        : growth rates table longitudinal
! Ax        : growth rates table horizontal
! Ay        : growth rates table vertical
! pnumInterp: rescaling number 
! alfap0    : growth rate longitudinal
! alfax0    : growth rate horizontal
! alfay0    : growth rate vertical
! ===============================================================================================================
subroutine ibsInterpolat(pnumber,epsx,epsy,sigs,dponp,xmin,&
xmax,ymin,ymax,zmin,zmax,xBins,yBins,zBins,Ap,Ax,Ay,pnumInterp,alfap0,alfax0,alfay0)
!     trilinear interpolation of ibs growth rates Ai(x,y,z) , for i=x,y,z, from tabulated values in external file
!     reading A as a 3D array A(i,j,k)
!     reference: see for example http://en.wikipedia.org/wiki/Trilinear_interpolation

integer, intent(in) :: xBins,yBins,zBins,pnumInterp
double precision, intent(in) :: pnumber,epsx,epsy,sigs,dponp
double precision, intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
double precision, intent(inout) :: alfap0,alfax0,alfay0
double precision, dimension(xBins,yBins,zBins), intent(in) :: Ax,Ay,Ap
integer :: i,j,k
double precision :: a000,a001,a010,a100,a110,a101,a011,a111,a00,a01,a10,a11,a0,a1,invDeltaX,invDeltaZ,invDeltaY
double precision :: DeltaX,DeltaY,DeltaZ,x0,y0,z0
double precision :: x,y,z ! point where function value should be interpolated

!     rewrite input arguments as x,y,z for use in trilinear program
      x=(epsx+epsy)/2. ! average transverse emittance
      y=sigs
      z=dponp

!     if arguments outside tabulated range, exit
      if ((x>xmax).or.(x<xmin).or.(y>ymax).or.(y<ymin).or.(z>zmax).or.(z<zmin)) then
         write(*,*) 'input argument out of bounds'
         write(*,*) 'x=',x,' xmin= ',xmin,' xmax=',xmax
         write(*,*) 'y=',y,' ymin= ',ymin,' ymax=',ymax
         write(*,*) 'z=',z,' zmin= ',zmin,' zmax=',zmax
         stop
      endif

!     find correct bins in 000 corner
      i=int(1+(x-xMin)/(xMax-xMin)*(xBins-1))
      j=int(1+(y-ymin)/(ymax-ymin)*(yBins-1))
      k=int(1+(z-zmin)/(zmax-zmin)*(zBins-1))

!     calculate bin width in all 3 directions.
      DeltaX=(xMax-xMin)/(xBins-1)
      DeltaY=(yMax-yMin)/(yBins-1)
      DeltaZ=(zMax-zMin)/(zBins-1)

      invDeltaX=1/DeltaX
      invDeltaY=1/DeltaY
      invDeltaZ=1/DeltaZ

!     calculate the lower values of (x,y,z) in the corner of the cube
      x0=xmin+DeltaX*(i-1)
      y0=ymin+DeltaY*(j-1)
      z0=zmin+DeltaZ*(k-1)

!     interpolate alfap
!------------------------------------------------------
!     get tabulated function values in corners of cube
      a000=Ap(i,j,k)
      a100=Ap(i+1,j,k)
      a010=Ap(i,j+1,k)
      a001=Ap(i,j,k+1)
      a110=Ap(i+1,j+1,k)
      a101=Ap(i+1,j,k+1)
      a011=Ap(i,j+1,k+1)
      a111=Ap(i+1,j+1,k+1)

!     find 4 new T-values by interpolation along x (middle of each side along x in cube)
      a00=a000+(a100-a000)*invDeltaX*(x-x0)
      a10=a010+(a110-a010)*invDeltaX*(x-x0)
      a01=a001+(a101-a001)*invDeltaX*(x-x0)
      a11=a011+(a111-a011)*invDeltaX*(x-x0)

!     find 2 new T-values by interpolation along y (middle points of the two planes with different y-values)
      a0=a00+(a10-a00)*invDeltaY*(y-y0)
      a1=a01+(a11-a01)*invDeltaY*(y-y0)

!     find final T-value by interpolation along z to the point inside cube
      alfap0=a0+(a1-a0)*invDeltaZ*(z-z0)

!     interpolate alfax
!------------------------------------------------------
!     get tabulated function values in corners of cube
      a000=Ax(i,j,k)
      a100=Ax(i+1,j,k)
      a010=Ax(i,j+1,k)
      a001=Ax(i,j,k+1)
      a110=Ax(i+1,j+1,k)
      a101=Ax(i+1,j,k+1)
      a011=Ax(i,j+1,k+1)
      a111=Ax(i+1,j+1,k+1)

!     find 4 new T-values by interpolation along x (middle of each side along x in cube)
      a00=a000+(a100-a000)*invDeltaX*(x-x0)
      a10=a010+(a110-a010)*invDeltaX*(x-x0)
      a01=a001+(a101-a001)*invDeltaX*(x-x0)
      a11=a011+(a111-a011)*invDeltaX*(x-x0)

!     find 2 new T-values by interpolation along y (middle points of the two planes with different y-values)
      a0=a00+(a10-a00)*invDeltaY*(y-y0)
      a1=a01+(a11-a01)*invDeltaY*(y-y0)

!     find final T-value by interpolation along z to the point inside cube
      alfax0=a0+(a1-a0)*invDeltaZ*(z-z0)

!     interpolate alfay
!------------------------------------------------------
!     get tabulated function values in corners of cube
      a000=Ay(i,j,k)
      a100=Ay(i+1,j,k)
      a010=Ay(i,j+1,k)
      a001=Ay(i,j,k+1)
      a110=Ay(i+1,j+1,k)
      a101=Ay(i+1,j,k+1)
      a011=Ay(i,j+1,k+1)
      a111=Ay(i+1,j+1,k+1)

!     find 4 new T-values by interpolation along x (middle of each side along x in cube)
      a00=a000+(a100-a000)*invDeltaX*(x-x0)
      a10=a010+(a110-a010)*invDeltaX*(x-x0)
      a01=a001+(a101-a001)*invDeltaX*(x-x0)
      a11=a011+(a111-a011)*invDeltaX*(x-x0)

!     find 2 new T-values by interpolation along y (middle points of the two planes with different y-values)
      a0=a00+(a10-a00)*invDeltaY*(y-y0)
      a1=a01+(a11-a01)*invDeltaY*(y-y0)

!     find final T-value by interpolation along z to the point inside cube
      alfay0=a0+(a1-a0)*invDeltaZ*(z-z0)

!     rescale the bunch population from the fixed value used in tabulation
      alfap0=alfap0*pnumber/pnumInterp
      alfax0=alfax0*pnumber/pnumInterp
      alfay0=alfay0*pnumber/pnumInterp

      return
      end

