!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: bcoefjas.f90,v $:
! $Revision: 1.5 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Assumes 
!     1. no spin trace
!     2. no trace over mg (Green's functions m_l)
!     3. only one initial state m_j

!KJ      subroutine bcoefjas(kinit,minit,kfinmax,jmax,ljmax, ispin, kind,lind,lgind,ljind,lx,hbmat,indmax,nspx)
      subroutine bcoefjas(kinit,minit,lx,hbmat)

      use nrixs_inp,only: lgind,ljind,indmax,kfinmax,jmax,ljmax
      implicit none
      complex*16, parameter :: coni=(0.d0,1.d0)
      integer,intent(in) :: kinit,minit,lx
      real*8,intent(out) :: hbmat(0:1,kfinmax)

!     local stuff
      integer jind(kfinmax),mind(kfinmax)
      integer jinit,ainit,jspin,ient,jfin,lfin,afin
      real*8 simp3j(kfinmax),lstoj(0:1,kfinmax)
      integer i1,j1,lj,ind,mfin,msg,mg,lg,lt
      real*8,external :: cwig3j


!     get the factor of 2 into play
      ient=2
      jspin=1
      jinit=2*abs(kinit)-1
      if (kinit.gt.0) then
         ainit=-1
      else
         ainit=1
      end if

!     Now set up some local working arrays.  (The others have moved to COMMON/m_nrixs.f90.)
!     (notice again lj=2*lj, and ljmax=2*ljmax)
      ind=0
      do lj=0,abs(ljmax)
!     triangle
         do jfin=max(abs(2*lj-jinit),1), min(2*lj+jinit,jmax), 2
!     use parity as an way to check/get afin  
            if (mod(jinit+jfin+2*lj,4).eq.0) then
               afin=-ainit
            else
               afin=ainit
            end if
            if (afin.gt.0) then 
               lfin=jfin-1
            else
               lfin=jfin+1
            end if
            if (lfin.le.2*lx) then 
               ind=ind+1
               jind(ind)=jfin
               mind(ind)=minit
            end if
         end do                 ! jfin
      end do                    ! lj


!     start filling up the matrices
!     since same jfin can show up many times 
!     only do one run
!     this is 
!     jfin lj linit
!     -minit 0 minit
      do j1=1,indmax
         jfin=jind(j1)
         lj=ljind(j1) !KJ ljind does not need to be *2
         lg=lgind(j1)*2  !KJ have to add *2 now 7-09
         mfin=mind(j1)
         lt=2*lj
         if (abs(minit).le.jfin) then 
            simp3j(j1)=cwig3j(jinit,lt,jfin,-minit,0,ient)
         else
            simp3j(j1)=0.0d0
         end if
         if (mod(minit+1,4).ne.0) simp3j(j1)=-simp3j(j1)
         
         do i1=0,1
            lstoj(i1,j1)=0.0d0
            if (abs(mfin).le.jfin .and. abs(mfin)-1.le.2*lx) then
               msg=2*i1-1
               mg=mfin-msg
!     This assumes that initial spin is the same as the final 
!     and since mfin=minit -> mg=ml_init
               lstoj(i1,j1)= cwig3j(lg,jspin,jfin,mg,msg,ient)
                  if (mod(lg-jspin+mfin,4).ne.0) then 
                     lstoj(i1,j1)=-lstoj(i1,j1)
                  end if  
                  lstoj(i1,j1)=sqrt(dble(jfin+1)) *lstoj(i1,j1)
            end if
         end do
      end do

!     start doing multiplications
	  hbmat(0:1,1:indmax)=dcmplx(0.d0)
      
!     do the multiplication and 
!     get rid of the factor of 2, in lgind,lind,and kind
      do j1=1,indmax
         do i1=0,1
            hbmat(i1,j1)=lstoj(i1,j1)*simp3j(j1)
         end do 
      end do 
         

      return
      end






