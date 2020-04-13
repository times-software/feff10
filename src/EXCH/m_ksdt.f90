! Copyright (C) 2014 Orbital-free DFT group at University of Florida, USA
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!  
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!  
! You should have received a copy of the GNU Lesser General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
! 


module ksdt_mod
implicit none
private
public :: ksdt_vxc, fxc_ksdt_01, exc_ksdt_01
contains

subroutine ksdt_vxc( rs, t, vxc)
! vxc potential according to KSDT
! 
! INPUT: rs, the Wigner-Seitz radius [Bohr],
!        t,  k_b * T  [Ha]
! OUTPUT: vxc, xc correlation potential [Ha]
! 
! Lazaro Calderin, Mar/2015
!     
   implicit none
   real(8), intent(in)  :: rs, t
   real(8), intent(out) :: vxc

   real(8) :: fxc, tr
   integer :: iz

   real(8), parameter :: pi = 4.0d0*atan(1.0d0)

   iz=0 ! No spin 

   ! Reduced temperature
   tr = 2.0d0*t*rs*rs /( 9.0d0*pi/4.0d0 )**(2.0d0/3.0d0)

   call fxc_ksdt_01(fxc,vxc,rs,tr,iz)

end subroutine ksdt_vxc

!-------------------------------------------------------------------------------
subroutine fxc_ksdt_01(fxc,vxc,rs,t,iz)
  !-----------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !   LDA XC free-energy parameterization from Monte Carlo data (unpol/pol) 
  !
  !    rs: Wigner-Seitz radius (bohr)
  !     t: reduced temperature
  !    iz: 0 - spin-unpolarized, 1 - fully polarized
  !
  !   fxc: XC free-energy per particle (hartree)
  !   vxc: xc potential (hartree)
  !
  ! REFERENCES:
  ! 

  !
  !-----------------------------------------------------------------------------
  ! REVISION LOG:
  !  21-NOV-2013 Subroutine created (V.V. Karasiev)
  !-----------------------------------------------------------------------------
  !
  implicit none
  real(8), intent(in) :: rs,t
  integer, intent(in) :: iz
  real(8), intent(out) :: fxc,vxc
  !
  real(8), parameter :: &
       onethird=1.D0/3.D0, &
       threehalf=1.5d0, &
       pi=4.d0*atan(1.d0), &
       lambda=(4.d0/9.d0/pi)**onethird, &
       a0=1.d0/(pi*lambda)
  real(8) :: aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee
  real(8) :: dtdn,tanht,dtanht,tanhsqrt,dtanhsqrt
  real(8) :: f1,dfxcdt,dfxcdrs
  real(8) :: num,dnum,den,dden,dnumdrs,ddendrs,n,drsdn
  !
  real(8) :: omega,a(6),b(0:1,4),c(0:1,3),d(0:1,5),e(0:1,5)
  !
  data a/0.75d0,3.04363d0,-0.092270d0,1.70350d0,8.31051d0,5.1105d0/
  !
  data b(0,:)/0.283997d0,48.932154d0,0.370919d0,61.095357d0/
  data c(0,:)/0.870089d0,0.193077d0,2.414644d0/
  data d(0,:)/0.579824d0,94.537454d0,97.839603d0,59.939999d0,24.388037d0/
  data e(0,:)/0.212036d0,16.731249d0,28.485792d0,34.028876d0,17.235515d0/
  !
  data b(1,:)/0.329001d0,111.598308d0,0.537053d0,105.086663d0/
  data c(1,:)/0.848930d0,0.167952d0,0.088820d0/
  data d(1,:)/0.551330d0,180.213159d0,134.486231d0,103.861695d0,17.750710d0/
  data e(1,:)/0.153124d0,19.543945d0,43.400337d0,120.255145d0,15.662836d0/
  !
  if(iz==0) then
    omega=1.d0
  elseif(iz==1) then
    omega=2.d0**onethird
  endif
!
  if(t==0.d0) then
!fxc
    f1=-1.d0/rs
    num=omega*a0*a(1)+b(iz,1)*sqrt(rs)+c(iz,1)*e(iz,1)*rs
    den=1.d0+d(iz,1)*sqrt(rs)+e(iz,1)*rs
    fxc=f1*num/den
!
    dnumdrs=b(iz,1)/sqrt(rs)/2.d0+c(iz,1)*e(iz,1)
    ddendrs=d(iz,1)/sqrt(rs)/2.d0+e(iz,1)

    n=3.d0/(4.d0*pi*rs**3) ! density
    drsdn=-onethird*rs/n ! (drs/dn)
    !Fxc=n*fxc=n*f1*num/den=n*A*B*C
    !dFxc/dn=fxc+n*(dA/dn)*B*C+n*a*(dB/dn)*C+n*a*B*(dC/dn)
    vxc=fxc+(onethird*f1)*num/den &   ! fxc+n*(dA/dn)*B*C
           + n*f1*(dnumdrs*drsdn)/den &  ! n*a*(dB/dn)*C
           - n*f1*num*(ddendrs*drsdn)/den**2 ! n*a*B*(dC/dn)   
  else
    tanht=tanh(1.d0/t)
    tanhsqrt=tanh(1.d0/sqrt(t))
    dtanht=(tanht**2-1.d0)/t**2 !d/dt tanh(1/t)
    dtanhsqrt=(tanhsqrt**2-1.d0)/t**threehalf/2.d0 !d/dt tanh(1/sqrt(t))
!
! a(t)
    num=a(1)+a(2)*t**2+a(3)*t**3+a(4)*t**4
    den=1.d0+a(5)*t**2+a(6)*t**4
!
    dnum=a(2)*2.d0*t+a(3)*3.d0*t**2+a(4)*4.d0*t**3
    dden=a(5)*2.d0*t+a(6)*4.d0*t**3
! 
    aa=a0*tanht*num/den
    daa=a0*(dtanht*num/den+tanht*dnum/den-tanht*num*dden/den**2)
! 
! b(t)
    num=b(iz,1)+b(iz,2)*t**2+b(iz,3)*t**4
    den=1.d0+b(iz,4)*t**2+omega*sqrt(3.d0)*b(iz,3)/sqrt(2.d0*lambda**2)*t**4
! 
    dnum=b(iz,2)*2.d0*t+b(iz,3)*4.d0*t**3
    dden=b(iz,4)*2.d0*t+omega*sqrt(3.d0)*b(iz,3)/sqrt(2.d0*lambda**2)*4.d0*t**3
!
    bb=tanhsqrt*num/den
    dbb=dtanhsqrt*num/den+tanhsqrt*dnum/den-tanhsqrt*num*dden/den**2
!
! d(t)
    num=d(iz,1)+d(iz,2)*t**2+d(iz,3)*t**4
    den=1.d0+d(iz,4)*t**2+d(iz,5)*t**4
! 
    dnum=d(iz,2)*2.d0*t+d(iz,3)*4.d0*t**3
    dden=d(iz,4)*2.d0*t+d(iz,5)*4.d0*t**3
!
    dd=tanhsqrt*num/den
    ddd=dtanhsqrt*num/den+tanhsqrt*dnum/den-tanhsqrt*num*dden/den**2
!
! e(t)
    num=e(iz,1)+e(iz,2)*t**2+e(iz,3)*t**4
    den=1.d0+e(iz,4)*t**2+e(iz,5)*t**4
! 
    dnum=e(iz,2)*2.d0*t+e(iz,3)*4.d0*t**3
    dden=e(iz,4)*2.d0*t+e(iz,5)*4.d0*t**3
!
    ee=tanht*num/den
    dee=dtanht*num/den+tanht*dnum/den-tanht*num*dden/den**2
!
! c(t)
    num=c(iz,1)+c(iz,2)*exp(-c(iz,3)/t)
    dnum=c(iz,2)*c(iz,3)*exp(-c(iz,3)/t)/t**2
    cc=num*ee
    dcc=dnum*ee+num*dee
!
!fxc
    f1=-1.d0/rs
    num=omega*aa+bb*sqrt(rs)+cc*rs
    den=1.d0+dd*sqrt(rs)+ee*rs
    fxc=f1*num/den
!
    dnum=omega*daa+dbb*sqrt(rs)+dcc*rs
    dnumdrs=bb/sqrt(rs)/2.d0+cc
    dden=ddd*sqrt(rs)+dee*rs
    ddendrs=dd/sqrt(rs)/2.d0+ee

    n=3.d0/(4.d0*pi*rs**3) ! density
    dtdn = -2.d0/3.d0*t/n ! (dt/dn)
    drsdn=-onethird*rs/n ! (drs/dn)
    !Fxc=n*fxc=n*f1*num/den=n*A*B*C
    !dFxc/dn=fxc+n*(dA/dn)*B*C+n*a*(dB/dn)*C+n*a*B*(dC/dn)
    vxc=fxc+(onethird*f1)*num/den &   ! fxc+n*(dA/dn)*B*C
      + n*f1*(dnum*dtdn+dnumdrs*drsdn)/den &  ! n*a*(dB/dn)*C
      - n*f1*num*(dden*dtdn+ddendrs*drsdn)/den**2 ! n*a*B*(dC/dn)
  endif
  return
end subroutine fxc_ksdt_01
!
!-------------------------------------------------------------------------------
subroutine exc_ksdt_01(exc,rs,t,iz)
  !-----------------------------------------------------------------------------
  !
  ! DESCRIPTION:
  !   LDA XC internal-energy from parameterization of Monte Carlo data (unpol/pol) 
  !
  !    rs: Wigner-Seitz radius (bohr)
  !     t: reduced temperature
  !    iz: 0 - spin-unpolarized, 1 - fully polarized
  !
  !   exc: XC internal-energy per particle (hartree)
  !
  ! REFERENCES:
  ! 

  !
  !-----------------------------------------------------------------------------
  ! REVISION LOG:
  !  24-NOV-2013 Subroutine created (V.V. Karasiev)
  !-----------------------------------------------------------------------------
  !
  implicit none
  real(8), intent(in) :: rs,t
  integer, intent(in) :: iz
  real(8), intent(out) :: exc
  real(8)              :: fxc,vxc,sxc,tempF
  !
  real(8), parameter :: &
       onethird=1.D0/3.D0, &
       threehalf=1.5d0, &
       pi=4.d0*atan(1.d0), &
       lambda=(4.d0/9.d0/pi)**onethird, &
       a0=1.d0/(pi*lambda)
  real(8) :: aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee
  real(8) :: dtdn,tanht,dtanht,tanhsqrt,dtanhsqrt
  real(8) :: f1,dfxcdt,dfxcdrs
  real(8) :: num,dnum,den,dden,dnumdrs,ddendrs,n,drsdn
  !
  real(8) :: omega,a(6),b(0:1,4),c(0:1,3),d(0:1,5),e(0:1,5)
  !
  data a/0.75d0,3.04363d0,-0.092270d0,1.70350d0,8.31051d0,5.1105d0/
  !
  data b(0,:)/0.283997d0,48.932154d0,0.370919d0,61.095357d0/
  data c(0,:)/0.870089d0,0.193077d0,2.414644d0/
  data d(0,:)/0.579824d0,94.537454d0,97.839603d0,59.939999d0,24.388037d0/
  data e(0,:)/0.212036d0,16.731249d0,28.485792d0,34.028876d0,17.235515d0/
  !
  data b(1,:)/0.329001d0,111.598308d0,0.537053d0,105.086663d0/
  data c(1,:)/0.848930d0,0.167952d0,0.088820d0/
  data d(1,:)/0.551330d0,180.213159d0,134.486231d0,103.861695d0,17.750710d0/
  data e(1,:)/0.153124d0,19.543945d0,43.400337d0,120.255145d0,15.662836d0/
  !
  if(iz==0) then
    omega=1.d0
  elseif(iz==1) then
    omega=2.d0**onethird
  endif
!
  if(t==0.d0) then
    call fxc_ksdt_01(fxc,vxc,rs,t,iz)
    exc=fxc
  else
    tanht=tanh(1.d0/t)
    tanhsqrt=tanh(1.d0/sqrt(t))
    dtanht=(tanht**2-1.d0)/t**2 !d/dt tanh(1/t)
    dtanhsqrt=(tanhsqrt**2-1.d0)/t**threehalf/2.d0 !d/dt tanh(1/sqrt(t))
!
! a(t)
    num=a(1)+a(2)*t**2+a(3)*t**3+a(4)*t**4
    den=1.d0+a(5)*t**2+a(6)*t**4
!
    dnum=a(2)*2.d0*t+a(3)*3.d0*t**2+a(4)*4.d0*t**3
    dden=a(5)*2.d0*t+a(6)*4.d0*t**3
! 
    aa=a0*tanht*num/den
    daa=a0*(dtanht*num/den+tanht*dnum/den-tanht*num*dden/den**2)
! 
! b(t)
    num=b(iz,1)+b(iz,2)*t**2+b(iz,3)*t**4
    den=1.d0+b(iz,4)*t**2+omega*sqrt(3.d0)*b(iz,3)/sqrt(2.d0*lambda**2)*t**4
! 
    dnum=b(iz,2)*2.d0*t+b(iz,3)*4.d0*t**3
    dden=b(iz,4)*2.d0*t+omega*sqrt(3.d0)*b(iz,3)/sqrt(2.d0*lambda**2)*4.d0*t**3
!
    bb=tanhsqrt*num/den
    dbb=dtanhsqrt*num/den+tanhsqrt*dnum/den-tanhsqrt*num*dden/den**2
!
! d(t)
    num=d(iz,1)+d(iz,2)*t**2+d(iz,3)*t**4
    den=1.d0+d(iz,4)*t**2+d(iz,5)*t**4
! 
    dnum=d(iz,2)*2.d0*t+d(iz,3)*4.d0*t**3
    dden=d(iz,4)*2.d0*t+d(iz,5)*4.d0*t**3
!
    dd=tanhsqrt*num/den
    ddd=dtanhsqrt*num/den+tanhsqrt*dnum/den-tanhsqrt*num*dden/den**2
!
! e(t)
    num=e(iz,1)+e(iz,2)*t**2+e(iz,3)*t**4
    den=1.d0+e(iz,4)*t**2+e(iz,5)*t**4
! 
    dnum=e(iz,2)*2.d0*t+e(iz,3)*4.d0*t**3
    dden=e(iz,4)*2.d0*t+e(iz,5)*4.d0*t**3
!
    ee=tanht*num/den
    dee=dtanht*num/den+tanht*dnum/den-tanht*num*dden/den**2
!
! c(t)
    num=c(iz,1)+c(iz,2)*exp(-c(iz,3)/t)
    dnum=c(iz,2)*c(iz,3)*exp(-c(iz,3)/t)/t**2
    cc=num*ee
    dcc=dnum*ee+num*dee
!
!fxc
    f1=-1.d0/rs
    num=omega*aa+bb*sqrt(rs)+cc*rs
    den=1.d0+dd*sqrt(rs)+ee*rs
    fxc=f1*num/den
!
    dnum=omega*daa+dbb*sqrt(rs)+dcc*rs
    dden=ddd*sqrt(rs)+dee*rs

    n=3.d0/(4.d0*pi*rs**3) ! density
    tempF = (3.d0*PI**2*n)**(2.d0/3.d0)/2.d0*omega**2
    !dfxc/dt=A*(dB/dt)*C+A*B*(dC/dt)
    sxc= f1*(dnum)/den &  ! A*(dB/dn)*C
       - f1*num*(dden)/den**2 ! A*B*(dC/dn)
    sxc=-sxc/tempF !sxc=-(t/T)*dfxc/dt=-(1/tempF)*dfxc/dt
    exc=fxc+t*tempF*sxc !exc=fxc+T*sxc=fxc+(t*tempF)*sxc
  endif
  return
end subroutine exc_ksdt_01
end module ksdt_mod
