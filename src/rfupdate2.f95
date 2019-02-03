subroutine rfupdate(np,x,px,y,py,t,pt,fmix,nLostDebunch,omega0,v00,thib,tcoeff,nharm1,nharm2,vrf1,vrf2,kept)
implicit none

integer :: koo,k
integer, intent(out) :: kept
!f2py intent(out) :: kept
integer, intent(inout) :: np
!f2py intent(in,out) ::np
double precision, intent(in) :: fmix,omega0,v00,thib,tcoeff,nharm1,nharm2,vrf1,vrf2
double precision :: volt,ptk,tk,p1,p2,thib2 !ptot,ham
double precision, intent(inout) :: nLostDebunch
!f2py intent(in,out) :: nLostDebunch
double precision, intent(inout), dimension(np) :: x, px,t,pt,y,py
!f2py intent(in,out) ::  x, px,t,pt,y,py


! integer, intent(in) :: iwrite
! iwrite is usually used to run checks when debugging the code
      koo   = 0
!      pttot = 0.0
      thib2 = thib/2.0
      do k=1,np
         tk = t(k)
         ptk = pt(k)
! update the time
         tk = tk + fmix*tcoeff*ptk
         p1 = tk*nharm1*omega0
         p2 = tk*nharm2*omega0
         volt = vrf1*sin(p1)+vrf2*sin(p2)
         ptk = ptk + fmix*volt/v00
!         pttot=pttot+ptk
! particles with abs(tk).gt.thib are lost. if thib==0, an infinite boundary is assumed and all particles are kept.
         if((abs(tk).le.thib2).or.(thib.eq.0.0)) then
            koo = koo + 1
            t(koo)=tk
            pt(koo)=ptk
            x(koo)=x(k)
            y(koo)=y(k)
            px(koo)=px(k)
            py(koo)=py(k)
!            ham=0.5*ptk*ptk*tcoeff+(cos(p1)-1)*vrf1/(nharm1*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0)
! longitudinal hamiltonian
!            xpavg(koo)=ham
         endif
      enddo
      nLostDebunch=nLostDebunch+ np - koo
      kept = koo
! put the normaization in np+1
      p1 = 0.5*thib*nharm1*omega0
      p2 = 0.5*thib*nharm2*omega0
!      xpavg(np+1)= (cos(p1)-1)*vrf1/(nharm1*v00*omega0) + (cos(p2)-1)*vrf2/(nharm2*v00*omega0)
! write the average long. momentum as a check. radiation and RF should cancel, mean should be 0
!      if(iwrite.eq.1) write(16,*) 'mean pt=',pttot/real(np)
      return
end
