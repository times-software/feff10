!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: makerotations.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2010/05/27 22:42:30 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine makerotations
      use boundaries,only : nkmmax,nl=>nlmax
        use trafo,only : mrotr
        use kklist,only : drot
!        use struct
      implicit none

        real*8 eulangcub(3,24),alfa,beta,gamma,w
      real*8     fact(0:100)
      character*4 symsymcub(48)
        integer i,l,m,isym,isymp,nkm,k,m1,m2,kap,lmax,kmin,kmax,ipos
        integer nmu
        real*8 mu1,mu2,j,cb,sb,db
        complex*16 eima,eimg
        integer,parameter :: nsymh=24,iinv=25
      complex*16 C0,C1,CI
      complex*16 CS,USC(3,3),W3X3(3,3)
      parameter (C0=(0.0D0,0.0D0),C1=(1.0D0,0.0D0),CI=(0.0D0,1.0D0))
      REAL*8 PI
      PARAMETER ( PI = 3.141592653589793238462643D0 )   
	  logical,parameter :: debug=.false.

!=======================================================================
! the following tables give the Euler angles corresponding to the 
! proper rotations of the cubic group Oh and the hexagonal group D6h
! taken from Bradley & Cracknell's table 1.4. 
! Note: These authors use active fixed rotation axes
!       Rose's convention used here is to use ACTIVE TEMPORARY AXES
!       accordingly alpha and gamme had to be interchanged
!
      DATA EULANGCUB /                                                  &
     &      0.0D0,    0.0D0,    0.0D0,                                  &
     &      0.0D0,  180.0D0,  180.0D0,                                  &
     &      0.0D0,  180.0D0,    0.0D0,                                  &
     &      0.0D0,    0.0D0,  180.0D0,                                  &
     &      0.0D0,   90.0D0,   90.0D0,                                  &
     &    180.0D0,   90.0D0,  270.0D0,                                  &
     &    180.0D0,   90.0D0,   90.0D0,                                  &
     &      0.0D0,   90.0D0,  270.0D0,                                  &
     &     90.0D0,   90.0D0,  180.0D0,                                  &
     &    270.0D0,   90.0D0,    0.0D0,                                  &
     &     90.0D0,   90.0D0,    0.0D0,                                  &
     &    270.0D0,   90.0D0,  180.0D0,                                  &
     &    270.0D0,   90.0D0,   90.0D0,                                  &
     &      0.0D0,   90.0D0,    0.0D0,                                  &
     &      0.0D0,    0.0D0,   90.0D0,                                  &
     &     90.0D0,   90.0D0,  270.0D0,                                  &
     &    180.0D0,   90.0D0,  180.0D0,                                  &
     &      0.0D0,    0.0D0,  270.0D0,                                  &
     &      0.0D0,  180.0D0,   90.0D0,                                  &
     &      0.0D0,  180.0D0,  270.0D0,                                  &
     &      0.0D0,   90.0D0,  180.0D0,                                  &
     &     90.0D0,   90.0D0,   90.0D0,                                  &
     &    180.0D0,   90.0D0,    0.0D0,                                  &
     &    270.0D0,   90.0D0,  270.0D0/          
!
      DATA SYMSYMCUB / 'E   ','C2x ','C2y ','C2z ','C+31','C+32',       &
     &                 'C+33','C+34','C-31','C-32','C-33','C-34',       &
     &                 'C+4x','C+4y','C+4z','C-4x','C-4y','C-4z',       &
     &                 'C2a ','C2b ','C2c ','C2d ','C2e ','C2f ',       &
     &                 'I   ','sx  ','sy  ','sz  ','S-61','S-62',       &
     &                 'S-63','S-64','S+61','S+62','S+63','S+64',       &
     &                 'S-4x','S-4y','S-4z','S+4x','S+4y','S+4z',       &
     &                 'sda ','sdb ','sdc ','sdd ','sde ','sdf '/


      fact(0) = dble(1)
      do i = 1,100
         fact(i) = fact(i-1)*dble(i)
      enddo

      lmax=nl-1
        drot=dcmplx(0)
      mrotr=dble(0)

!----------------
!-------------------------------------------------------
!  create transformation matrix   U  cartesian/sperical ccordinates
!-----------------------------------------------------------------------
!  RC,RCP  vectors in cartesian coordinates
!  RS,RSP  vectors in spherical coordinates
!         RS  = USC * RC                                 (4.40)
!         RSP = MS  * RS                                 (4.37)
!     MS(i,j) = D(j,i)                                   (4.42)
!     D  rotation matrix for complex spherical harmonics
!
      W = 1.0D0/SQRT(2.0D0)
!
! ordering of: m=-1,0,+1 >>> row 1 and 3 interchanged compared to (4.44)
      USC(1,1) = W
      USC(1,2) = -CI*W
      USC(1,3) = 0.0D0
      USC(2,1) = 0.0D0
      USC(2,2) = 0.0D0
      USC(2,3) = 1.0D0
      USC(3,1) = -W
      USC(3,2) = -CI*W
      USC(3,3) = 0.0D0


      nkm=nkmmax

      do i=1,nsymh

           alfa=eulangcub(1,i)*pi/dble(180)  ! degrees to radians
           beta=eulangcub(2,i)*pi/dble(180)
           gamma=eulangcub(3,i)*pi/dble(180)
         cb=dcos(beta/dble(2))
           sb=-dsin(beta/dble(2))

!  SET UP NONRELATIVISTIC ROTATION MATRIX
         do l=0,lmax
           do m1=-l,l
           do m2=-l,l

            eima=cdexp(-ci*m2*alfa)
              eimg=cdexp(-ci*m1*gamma)
              db=dble(0)
              kmin=max(0,m1-m2)
              kmax=min(l-m2,l+m1)
              do k=kmin,kmax
                 db=db+(-1)**k * cb**(2*l+m1-m2-2*k) * sb**(m2-m1+2*k)  &
     &          / (fact(l-m2-k)*fact(l+m1-k)*fact(k+m2-m1)*fact(k))
            enddo
              db=db*dsqrt(fact(l+m1)*fact(l-m1)*fact(l+m2)*fact(l-m2))

            drot(l*l+l+1+m2,l*l+l+1+m1,i,2)=eima*eimg*db

         enddo
           enddo
           enddo

!  SET UP RELATIVISTIC ROTATION MATRIX
         ipos=1
         do kap=1,2*lmax+1
           l=kap/2
           j=l+dble(0.5)
         if(2*l.eq.kap) j=j-dble(1)
           nmu=nint(2*j+1)

           do m1=0,nmu-1
           do m2=0,nmu-1
         mu1=-j+m1
           mu2=-j+m2

            eima=cdexp(-ci*mu2*alfa)
              eimg=cdexp(-ci*mu1*gamma)
              db=dble(0)
              kmin=max(0,nint(mu1-mu2))
              kmax=min(nint(j-mu2),nint(j+mu1))
              do k=kmin,kmax
                 db=db+(-1)**k * cb**(nint(2*j+mu1-mu2)-2*k) *          &
     &		 sb**(nint(mu2-mu1)+2*k) / (fact(nint(j-mu2)-k)                 &
     &         *fact(nint(j+mu1)-k)*fact(k+nint(mu2-mu1))*fact(k))
            enddo
              db=db*dsqrt(fact(nint(j+mu1))*fact(nint(j-mu1))*          &
     &               fact(nint(j+mu2))*fact(nint(j-mu2)))

            drot(ipos+m2,ipos+m1,i,1)=eima*eimg*db

         enddo
           enddo
           ipos=ipos+nmu
           enddo

      enddo



!-----------------------------------------------------------------------
! create the rotation matrix MROTR for vectors in cartesian coordinates
! NOTE:  U^+ D^T U gives the inverse of the real matrix  M
!        for that reason  the transposed matrix is stored as MROTR(J,I)
!-----------------------------------------------------------------------
      do ISYM = 1,NSYMH
!
         do I = 1,3
            do L = 1,3
               CS = 0.0D0
               do K = 1,3
                  CS = CS + DROT(K+1,I+1,ISYM,2)*USC(K,L)
               enddo
               W3X3(I,L) = CS
            enddo
         enddo
!
         do I = 1,3
            do L = 1,3
               CS = 0.0D0
               do K = 1,3
                  CS = CS + dconjg(USC(K,I))*W3X3(K,L)
               enddo
               if ( dimag(CS).GT.1D-8 ) write (*,*) 'ISYM=',ISYM,       &
     &              ' MROT',I,L,CS,' ???????????'
! see above >> MROTR(I,L,ISYM) = Dreal(CS)
               MROTR(L,I,ISYM) = Dreal(CS)
            enddo
         enddo
!
      enddo



!                     create matrix for inversion
!-----------------------------------------------------------------------
      mrotr(:,:,iinv)=dble(0)
        do i=1,3
           mrotr(i,i,iinv)=dble(-1)
        enddo

      drot(:,:,iinv,:)=dcmplx(0)
        i=0
        k=0
        do l=0,lmax
           do m=1,2*l+1
              i=i+1
              drot(i,i,iinv,2)=dcmplx(-1,0)**l
           enddo
           do m=1,2*(2*l+1)
              k=k+1
              drot(k,k,iinv,1)=dcmplx(-1,0)**l
           enddo
        enddo


!-----------------------------------------------------------------------
!                         include inversion
!-----------------------------------------------------------------------
      do ISYM = 2,NSYMH
         ISYMP = NSYMH + ISYM
        
         call ZGEMM('N','N',NKM,NKM,NKM,C1,DROT(1,1,ISYM,2),NKMMAX,     &
     &              DROT(1,1,IINV,2),NKMMAX,C0,DROT(1,1,ISYMP,2),NKMMAX)
         call ZGEMM('N','N',NKM,NKM,NKM,C1,DROT(1,1,ISYM,1),NKMMAX,     &
     &              DROT(1,1,IINV,1),NKMMAX,C0,DROT(1,1,ISYMP,1),NKMMAX)

         call DGEMM('N','N',3,3,3,1D0,MROTR(1,1,ISYM),3,MROTR(1,1,IINV),&
     &              3,0D0,MROTR(1,1,ISYMP),3)
      enddo

!=======================================================================
!            set up of transformation matrices completed
!=======================================================================

      if(debug)then
        open(77,file='mrotr.txt')
        do i=1,48
        if(i.le.nsymh) then
           write(77,*) 'SYMMETRY OPERATION ',i,symsymcub(i)
        else
           write(77,*) 'SYMMETRY OPERATION ',i,symsymcub(i-nsymh),'+I'
        endif
        do l=1,3
        do k=1,3
        if(dabs(mrotr(l,k,i)).lt.0.000001) mrotr(l,k,i)=dble(0)
        enddo
        write(77,'(3f14.6)') mrotr(l,:,i)
        enddo
        enddo
        close(77)
      endif


        return
        end

