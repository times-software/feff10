!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: cwig3j.f90,v $:
! $Revision: 1.3 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function cwig3j (j1,j2,j3,m1,m2,ient)
!     wigner 3j coefficient for integers  (ient=1)
!                         or semiintegers (ient=2)
!     other arguments should be multiplied by ient
 
      implicit double precision (a-h,o-z)
      integer,parameter :: idim = 58 !KJ 7-09 for NRIXS - used to be 58)
      character*512 slog
!     dimensions  modified for larger arguments by ala 12.12.94
      dimension al(idim+1),m(12)
      save ini, al
      data ini/1/
!     idim-1 is the largest argument of factorial to calculate

      m3=-m1-m2
      if (ini) 1,21,1
!        initialisation of the log's of the factorials
 1    ini=0
      al(1)=0.0d00
      do i=1,idim
         b=i
         al(i+1)=al(i)+ log(b)
      enddo
 21   cwig3j=0.0d00
      if (((ient-1)*(ient-2)).ne.0) go to 101
      ii=ient+ient
!        test triangular inequalities, parity and maximum values of m
      if (( abs(m1)+ abs(m2)).eq.0.and.mod(j1+j2+j3,ii).ne.0) go to 99
      m(1)=j1+j2-j3
      m(2)=j2+j3-j1
      m(3)=j3+j1-j2
      m(4)=j1+m1
      m(5)=j1-m1
      m(6)=j2+m2
      m(7)=j2-m2
      m(8)=j3+m3
      m(9)=j3-m3
      m(10)=j1+j2+j3+ient
      m(11)=j2-j3-m1
      m(12)=j1-j3+m2
      do 41 i=1,12
!	     write(*,*) 'i,m(i)',i,m(i)
!	     write(*,*) 'error1?'
         if (i.gt.10) go to 31
         if (m(i).lt.0) go to 99
!        write(*,*) 'error2?', mod(m(i),ient)
31       if (mod(m(i),ient).ne.0) go to 101
         m(i)=m(i)/ient
!		 write(*,*) 'error3?',idim
         if (m(i).gt.idim) go to 101
 41   continue

!        calculate 3j coefficient
      max0= max(m(11),m(12),0)+1
      min0= min(m(1),m(5),m(6))+1
      isig=1
      if (mod(max0-1,2).ne.0) isig=-isig
      c=-al(m(10)+1)
      do i=1,9
         c=c+al(m(i)+1)
      enddo
      c=c/2.0d00
      do i=max0,min0
         j=2-i
         b=al(i)+al(j+m(1))+al(j+m(5))+al(j+m(6))+al(i-m(11))+al(i-m(12))
         cwig3j=cwig3j+isig* exp(c-b)
         isig=-isig
      enddo
      if (mod(j1-j2-m3,ii).ne.0) cwig3j=-cwig3j
 99   return

 101  write(slog,'(a,6i5)') 'error in cwig3j ',j1,j2,j3,m1,m2,ient
      call wlog(slog)
      stop
      end
