!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: ppole.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function wdisp(q)
! dispertion relation
! input: q - momentum or wavenumber
! input from common blocks
!       ompl - zero q energy of mode
!       adisp - dispersion parameter for dispersion relation,
      implicit none
      double precision q
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      wdisp=sqrt(ompl**2+adisp*q**2+q**4/4.d0)
      return
      end

      double precision function dwdq(q)
! the derivative of the dispertion relation 
! with respect to the plasmon momentum q
! input: q - momentum or wavenumber
! input from common blocks
!       ompl - zero q energy of mode
!       adisp - dispersion parameter for dispersion relation,
      implicit none
      double precision q,wdisp
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      external wdisp
      dwdq=(q**3+2.d0*adisp*q)/(2.d0*wdisp(q))
      return
      end

      double precision function d2wdq2(q)
! the second derivative of the dispertion relation 
! with respect to the plasmon momentum q
! input: q - momentum or wavenumber
! input from common blocks
!       omp - plasma frequency omega_p
!       adisp - dispersion parameter for dispersion relation,
      implicit none
      double precision q,wdisp,dwdq
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      external wdisp,dwdq
      d2wdq2=(3.d0*q**2+2.d0*adisp)*wdisp(q)                            &
     &       -(q**3+2.d0*adisp*q)*dwdq(q)
      d2wdq2=d2wdq2/(2.d0*wdisp(q)**2)
      return
      end

      double precision function qdisp(w)
! The inverse dispertion relation
! input: w - energy (omega)
! input from common blocks
!       ompl - zero q energy of mode
!       adisp - dispersion parameter for dispersion relation,
!               w(q)**2=ompl+adisp*q**2+q**4/4
      implicit none
      double precision w,x,y
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      x=adisp**2+w**2-ompl**2
      if (x.ge.0.d0) then
        y=-2.d0*adisp+2.d0*dsqrt(x)
        if (y.ge.0.d0) then
          qdisp=dsqrt(y)
          return
        endif
      endif
      qdisp=0.d0
      return
      end

      double precision function qthresh(ompl,AA,ef,qf)
! Find the photoelectron momentum corresponding to the onset of plasmon
! losses.
! input: ompl - zero q energy of mode
!        AA - dispersion parameter omega(q)**2=ompl**2+AA*q**2+q**4/4
!        ef - fermi energy
!        qf - fermi momentum
      implicit none
      integer i,nrts,nrts2,nflip,nqa,nqb
      double precision AA,ompl,ef,qf,a,b,c,d,qthresh1,qthresh2,q01,     &
     &                 q02,xfact,test1,test2,test3,qh,ek,q1a,q2a,       &
     &                 q3a,q0a,q1b,q2b,q3b,q0b,w,x,y
      complex*16 rt1,rt2,rt3,rtt1,rtt2,rtt3,rt(3),rtt(3),qst
      a=1.d0
      b=-3.d0*AA
      c=3.d0*AA**2-27.d0*ompl**2/4.d0
      d=-AA**3
      call croots(a,b,c,d,rt1,rt2,rt3,nrts)
      rt(1)=rt1
      rt(2)=rt2
      rt(3)=rt3
      if (nrts.eq.1) then
 10     continue
          nflip=0
          do i=1,2
            if (dimag(rt(i)).lt.dimag(rt(i+1))) then
              qst=rt(i)
              rt(i)=rt(i+1)
              rt(i+1)=qst
              nflip=nflip+1
            endif
          enddo
        if (nflip.ne.0) goto 10
        qthresh1=dble(rt(2))
      else
        qthresh1=max(dble(rt1),dble(rt2),dble(rt3))
      endif
!      write(6,*) rt1,rt2,rt3
!      write(6,*) qthresh1
      if (qthresh1.gt.0.d0) then
        qthresh1=dsqrt(qthresh1)
      else
        qthresh1=0.d0
      endif
!      write(6,*) qthresh1
      a=1.d0
      b=1.5d0*qf+AA/qf
      c=qf**2+2.d0*AA
      d=qf**3/4.d0+AA*qf+ompl**2/qf
      call croots(a,b,c,d,rt1,rt2,rt3,nrts)
      rt(1)=rt1
      rt(2)=rt2
      rt(3)=rt3
      b=-b
      d=-d
      call croots(a,b,c,d,rtt1,rtt2,rtt3,nrts2)
      rtt(1)=rtt1
      rtt(2)=rtt2
      rtt(3)=rtt3
      if (nrts.eq.1) then
 11     continue
          nflip=0
          do i=1,2
            if (dimag(rt(i)).lt.dimag(rt(i+1))) then
              qst=rt(i)
              rt(i)=rt(i+1)
              rt(i+1)=qst
              nflip=nflip+1
            endif
          enddo
        if (nflip.ne.0) goto 11
        q01=dble(rt(2))
      else
        xfact=sqrt(AA**2+(dble(rt(1))**2/2.d0)**2-ompl**2)-AA
        test1=dble(rt(1))-qf-sqrt(2.d0*xfact)
        xfact=sqrt(AA**2+(dble(rt(2))**2/2.d0)**2-ompl**2)-AA
        test2=dble(rt(2))-qf-sqrt(2.d0*xfact)
        xfact=sqrt(AA**2+(dble(rt(3))**2/2.d0)**2-ompl**2)-AA
        test3=dble(rt(3))-qf-sqrt(2.d0*xfact)
        if(test1.lt.test2.and.test1.lt.test3) then
          q01=dble(rt(1))
        elseif(test2.lt.test3) then
          q01=dble(rt(2))
        else
          q01=dble(rt(3))
        endif
      endif
      if (nrts2.eq.1) then
 12     continue
          nflip=0
          do i=1,2
            if (dimag(rtt(i)).lt.dimag(rtt(i+1))) then
              qst=rtt(i)
              rtt(i)=rtt(i+1)
              rtt(i+1)=qst
              nflip=nflip+1
            endif
          enddo
        if (nflip.ne.0) goto 12
        q02=dble(rtt(2))
      else
        xfact=sqrt(AA**2+(dble(rtt(1))**2/2.d0)**2-ompl**2)-AA
        test1=dble(rtt(1))+qf-sqrt(2.d0*xfact)
        xfact=sqrt(AA**2+(dble(rtt(2))**2/2.d0)**2-ompl**2)-AA
        test2=dble(rtt(2))+qf-sqrt(2.d0*xfact)
        xfact=sqrt(AA**2+(dble(rtt(3))**2/2.d0)**2-ompl**2)-AA
        test3=dble(rtt(3))+qf-sqrt(2.d0*xfact)
        if(test1.lt.test2.and.test1.lt.test3) then
          q02=dble(rt(1))
        elseif(test2.lt.test3) then
          q02=dble(rt(2))
        else
          q02=dble(rt(3))
        endif
      endif
!      write(6,*) q01,q02
      qthresh2=min(abs(q01),abs(q02))
!      write(6,*) qthresh1,qthresh2

      qh=1000.d0*qf
      ek=qthresh1**2/2.d0
      call qlimits(ek,qthresh1,ompl,AA,qh,nqa,q1a,q2a,q3a)
      q0a=0.d0
      w=ek-ef
      x=AA**2+w**2-ompl**2
      if (x.ge.0.d0) then
        y=-2.d0*AA+2.d0*dsqrt(x)
        if (y.ge.0.d0) then
          q0a=dsqrt(y)
        endif
      endif
      ek=qthresh2**2/2.d0
      call qlimits(ek,qthresh2,ompl,AA,qh,nqb,q1b,q2b,q3b)
      q0b=0.d0
      w=ek-ef
      x=AA**2+w**2-ompl**2
      if (x.ge.0.d0) then
        y=-2.d0*AA+2.d0*dsqrt(x)
        if (y.ge.0.d0) then
          q0b=dsqrt(y)
        endif
      endif

      if (nqa.eq.0) then
        qthresh=qthresh1
      elseif(abs(q1a-q2a).lt.abs(q1b-q0b)) then
        qthresh=qthresh1
      else
        qthresh=qthresh2
      endif
!      write(6,*) qthresh
      return
      end

      double precision function vpp2(q)
! the square of the coupling potential
! input: q - momentum or wavenumber
! input from common blocks
!       pi - ratio of circumference to diameter of a circle in
!            euclidian geometry
!       ompl - zero q energy of mode
      implicit none
      double precision q,wdisp
      double precision pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      common /convsf/ pi,ef,xmu,qf,omp,ompl,wt,ekp,ek,pk,acc,brd,adisp
      external wdisp
      vpp2=2*pi*omp**2/(q**2*wdisp(q))
      return
      end

