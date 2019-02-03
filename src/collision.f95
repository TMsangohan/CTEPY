subroutine updatelumi(x,px,y,py,t,pt,xx,pxx,yy,pyy,tt,ptt,&
betax1,betay1,betax2,betay2,ex1,ey1,ex2,ey2,&
lumprefac,pnumber1,pnumber2,np10,np20,nlostb1,nlostb2,np1,np2,&
pi,iseed,kept1,kept2)
implicit none

integer, intent(in) :: np1,np2
!f2py intent(in) :: np1,np2
integer, intent(in) :: np10,np20
!f2py intent(in) :: np10,np20
integer, intent(out) :: kept1,kept2
!f2py intent(out) :: kept1,kept2

integer, intent(inout) :: nlostb1,nlostb2
!f2py intent(in,out) :: nlostb1,nlostb2

double precision, intent(inout) :: iseed
!f2py intent(in,out) :: iseed
double precision, intent(in) :: ex1,ey1,ex2,ey2,lumprefac,pi,pnumber1,pnumber2
!f2py intent(in) :: ex1,ey1,ex2,ey2,lumiprefac,pi,pnumber1,pnumber2

double precision, intent(in) :: betax1,betay1,betax2,betay2
!f2py intent(in) :: betas

double precision, dimension(np1), intent(inout) :: x,px,y,py,t,pt
!f2py intent(in,out) :: x,px,y,py,t,pt

double precision, dimension(np2), intent(inout) :: xx,pxx,yy,pyy,tt,ptt
!f2py intent(in,out) :: xx,pxx,yy,pyy,tt,ptt

double precision, external :: IO

integer :: k,ko,koo
double precision :: jx,jy,argx,argy,prob,sigx1,sigy1,sigx2,sigy2

!probavg1=0
!probavg2=0

! average beam size
sigx1 = sqrt(ex1*betax1)
sigy1 = sqrt(ey1*betay1)
sigx2 = sqrt(ex2*betax2)
sigy2 = sqrt(ey2*betay2)

ko = 0
koo= 0
do k=1,np1
   jx = (x(k)**2+px(k)**2)/(2*betax1) ! action variable. x1 given in ring so need to use beta of the ring
   jy = (y(k)**2+py(k)**2)/(2*betay1)
   argx = jx/(2*ex2)      ! argument for the density function
   argy = jy/(2*ey2)
   prob = lumprefac * np2 * exp(-argx)*IO(argx)*exp(-argy)*IO(argy)* & !interaction probability averaged over phase
              (pnumber2/np20)/ & ! scaling factor for the number of real particles in beam 2
              (2*pi*sigx2*sigy2)
!   probavg1 = probavg1 + prob
   call random_number(iseed)
   if (iseed.ge.prob) then
     ko     = ko + 1
     t(ko)  = t(k)
     pt(ko) = pt(k)
     x(ko)  = x(k)
     px(ko) = px(k)
     y(ko)  = y(k)
     py(ko) = py(k)
   endif
enddo
kept1 = ko
nlostb1 = nlostb1 + np1 - ko
do k=1,np2
   jx = (xx(k)**2+pxx(k)**2)/(2*betax2) ! action variable. x1 given in ring so need to use beta of the ring
   jy = (yy(k)**2+pyy(k)**2)/(2*betay2)
   argx = jx/(2*ex1)      ! argument for the density function
   argy = jy/(2*ey1)
   prob = lumprefac * np1 * exp(-argx)*IO(argx)*exp(-argy)*IO(argy)* & !interaction probability averaged over phase
              (pnumber1/np10)/ & ! scaling factor for the number of real particles in beam 2
              (2*pi*sigx1*sigy1)
!  probavg2 = probavg2 + prob
   call random_number(iseed)
   if (iseed.ge.prob) then
     koo     = koo + 1
     tt(koo)  = tt(k)
     ptt(koo) = ptt(k)
     xx(koo)  = xx(k)
     pxx(koo) = pxx(k)
     yy(koo)  = yy(k)
     pyy(koo) = pyy(k)
   endif
enddo
kept2 = koo  
nlostb2 = nlostb2 + np2 - koo
end


double precision function betalevel(lumimax,npp1,npp10,npp2,npp20,pnum1,pnum2,betas, &
hglassold,circ,ex1,ey1,ex2,ey2,clight,pi)
implicit none
integer, intent(in) :: npp1,npp10,npp2,npp20
!f2py intent(in) :: npp1,npp10,npp2,npp20
double precision, intent(in) :: betas,circ,ex1,ex2,ey1,ey2,clight,hglassold,pi,lumimax,pnum1,pnum2
!f2py intent(in) :: betas,circ,ex1,ex2,ey1,ey2,clight,hglassold,pi,lumimax,pnum1,pnum2
double precision :: H,N1,N2,f,lumitest

H  = hglassold
N1 = (npp1*pnum1/npp10)
N2 = (npp2*pnum2/npp20)
f  = (clight/circ)

lumitest = (H*N1*N2*f)/(2*pi*sqrt(betas**2 * (ex1 + ex2) * &
           (ey1 + ey2)))*1.e-4 !scaling factor 1.e-4 to get units cm^-2 s^-1
if (lumitest>=lumimax) then
    betalevel = (H*N1*N2*f)/(2*pi*(lumimax)*sqrt((ex1+ex2)*(ey1+ey2)))* 1.0e-4 
else
    betalevel = betas
endif
return
end function betalevel


subroutine binning(t,denlon,thib2,dtsamp,longintbins,np)
implicit none
integer, intent(in) :: np,longintbins
!f2py intent(in) :: np,longintbins
double precision, intent(in) :: thib2,dtsamp
!f2py intent(in) :: thib2,dtsamp
double precision, dimension(np), intent(in) :: t
!f2py intent(in) :: t,pt
double precision, dimension(np), intent(inout) :: denlon
!f2py intent(in,out) :: denlon
integer :: k,nkt

double precision tk

do k=1,np
   tk = (t(k) + thib2)/dtsamp
   nkt = nint(tk)
   if (nkt.le.0) nkt=1
   if (nkt.gt.longIntBins) nkt=longIntBins
      denlon(nkt)=denlon(nkt)+1 ! binning in t, but same binning valid in z with dzsamp=clight*dtsamp
   enddo
end subroutine binning

! same routine as before but for 1d collision calculations
subroutine binning3d(x,y,t,denlonx,denlony,denlon,ampx,ampy,thib2,ddx,ddy,dtsamp,np,longintbins)

integer, intent(in) :: np,longintbins
!f2py intent(in) :: np,longintbins
double precision, intent(in) :: thib2,dtsamp,ampx,ampy,ddx,ddy
!f2py intent(in) :: thib2,dtsamp,ampx,ampy,ddx,ddy
double precision, dimension(np), intent(in) :: x,y,t
!f2py intent(in) :: x,px,y,py, t,pt

double precision, dimension(np), intent(inout) :: denlonx,denlony,denlon
!f2py intent(in,out) :: denlonx,denlony,denlon

integer :: k,nkx,nky,nkt
double precision :: xk,yk,tk

do k=1,np
         xk  = (x(k)+5*ampx)/ddx
         nkx = nint(xk)  ! find nearest integer = bin number
         if(nkx.le.0)  nkx=1
         if(nkx.gt.longIntBins) nkx=longIntBins
         denlonx(nkx)=denlonx(nkx)+1

         yk = (y(k)+5*ampy)/ddy
         nky = nint(yk)  ! find nearest integer = bin number
         if(nky.le.0)  nky=1
         if(nky.gt.longIntBins) nky=longIntBins
         denlony(nky)=denlony(nky)+1

         tk = (t(k)+thib2)/dtsamp
         nkt = nint(tk)
         if(nkt.le.0) nkt=1
         if(nkt.gt.longIntBins) nkt=longIntBins
         denlon(nkt)=denlon(nkt)+1 ! binning in t, but same binning valid in z with dzsamp=clight*dtsamp
enddo
end subroutine binning3d

subroutine collision6a(sigI,thib,longintbins,timeRatio,circ,pi,clight,iseed,hglassold, &
nLostLum1,nLostLum2,nmacro1,nmacro2, &
pnumber1,np1,x1,px1,y1,py1,t1,pt1,ex1,ey1,bx1,by1,bs1,lumimax1,theta1,&
pnumber2,np2,x2,px2,y2,py2,t2,pt2,ex2,ey2,bx2,by2,&
hglassfac,betass,kept1,kept2,lumip1)
implicit none

integer, intent(in) :: np1,np2
integer, intent(in) :: nmacro1,nmacro2,longintbins
integer, intent(inout) :: nLostLum1,nLostLum2
integer, intent(out) :: kept1,kept2
double precision, intent(in) :: ex1,ey1,bx1,by1,bs1,lumimax1,theta1
double precision, intent(in) :: ex2,ey2,bx2,by2
double precision, intent(in) :: sigI,thib,timeRatio,circ,pi,clight,pnumber1,pnumber2
double precision, dimension(np1), intent(inout) :: x1,px1,y1,py1,t1,pt1 
double precision, dimension(np2), intent(inout) :: x2,px2,y2,py2,t2,pt2
double precision, intent(out) :: lumip1
double precision, intent(in) :: hglassold
double precision, intent(out) :: hglassfac,betass
double precision, intent(inout) :: iseed

integer :: k,m
double precision :: sigm2,dtsamp,thib2,f,H1,N11,N12
double precision, dimension(longintbins) :: denlon1,denlon2
double precision :: betas1,z1,z2
double precision :: thetas,nps,lumprefac
double precision, dimension(1,longintbins) :: denlons
double precision,external :: betalevel

sigm2  = sigI * 1.e-28   ! convert cross section from barn to m^2
dtsamp = thib/float(longIntBins) ! binsize longitudinal
thib2  = thib/2. 

! init longitudinal binning arrays
do k=1,longintbins
  denlon1(k) = 0.0   
  denlon2(k) = 0.0
enddo

! fill binning arrays
call binning(t1,denlon1,thib2,dtsamp,longintbins,np1)
call binning(t2,denlon2,thib2,dtsamp,longintbins,np2)

! if luminosity levelling is active, calculate the corresplonding beta to level
betas1 = bs1
    
if (lumimax1.ne.0) then
  betas1 = betalevel(lumimax1,np1,nmacro1,np2,nmacro2,pnumber1,pnumber2,betas1, &
                     hglassold,circ,ex1,ey1,ex2,ey2,clight,pi)
endif

! calculating hourglass factors
hglassfac = 0.0
betass = betas1
thetas = theta1

do k =1,longintbins
  denlons(1,k) = denlon2(k)
enddo

nps = np2

do k=1,longintbins
  z1 = ((real(k)/real(longintbins))*thib-thib2)*clight
     do m=1,longIntBins
        z2 = -((real(m)/real(longintbins))*thib-thib2)*clight
                hglassfac = hglassfac +  &
                    exp(-(2*(z1+z2)**2*betass*sin(thetas)**2)/ & ! theta is half of the crossing angle
                    ((ex1+ex2)*((z1+z2)**2+2*betass**2*(1+cos(2*thetas)))))* &
                    denlon1(k)/np1*denlons(1,m)/nps/(1+((z1+z2)/(2*betass*cos(thetas)))**2)
     enddo
enddo

lumprefac = hglassfac*sigm2*timeRatio / betass * sqrt(bx1*by1)

f = clight/circ
! main bunch colliding in IP 1
H1  = hglassfac
N11 = (np1*pnumber1/nmacro1)
N12 = (np2*pnumber2/nmacro2)
lumip1 = (H1*N11*N12*f)/(2* pi* betass *sqrt((ex1+ex2)*(ey1+ey2)))* 1.e-4
call updatelumi(x1,px1,y1,py1,t1,pt1,x2,px2,y2,py2,t2,pt2,&
        bx1,by1,bx2,by2,ex1,ey1,ex2,ey2,&
        lumprefac,pnumber1,pnumber2,nmacro1,nmacro2,nLostLum1,nLostLum2,np1,np2,&
        pi,iseed,kept1,kept2)

nLostLum1 = nLostLum1 + np1 - kept1
nLostLum2 = nLostLum2 + np2 - kept2

end subroutine collision6a

include 'IO.f95'
