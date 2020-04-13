!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: specupd.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2011/11/29 00:04:07 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine specupd(icalc,isp,kfinmax,indmax,indcalc,indmap, &
           ljind,jinit,hbmat,ljmax,xirflj, xsec,xnorm,imode)  !KJ ,qw)
	  use global_inp,only: imdff,mixdff,qw,nq,cosmdff !KJ cosleg=>cosmdff
      implicit none
      integer icalc,isp,kfinmax,indmax
      integer indcalc(kfinmax,3),indmap(kfinmax)
      integer ljind(kfinmax)
      integer jinit
      double precision hbmat(0:1,kfinmax,-jinit:jinit)
      integer ljmax
      complex*16 xirflj(0:ljmax)
      double precision xnorm,pleg(ljmax+1) !KJ ,qw
      complex*16 xsec(0:ljmax)
      integer imode
! c APS added cosleg and pleg the pleg is calculated with cpl0 MATH subroutine. cosleg is the angle between q and q'
      real*8 cosleg
	  real*8 plegqq(nq,nq,ljmax+1) !KJ legendre polynomials for all q,q' pairs
	  complex*16 qweights(nq)
	  integer iq,iqq,iqmin,iqmax,iqqmin,iqqmax
	  	  
!     Everyhing else is clear from the calling program
!     imode=1 is for the regular part
!     imode=2 is for the irregular
!

!
!  This part is for the trace stuff   
!
      double precision traceres
      integer ii
      integer mjinit
      integer lj
      complex*16 aa
      double complex xirnorm

!   KJ : put Adam's code in a loop
    do iq=1,nq
	do iqq=1,nq
	   cosleg=cosmdff(iq,iqq)
!      APS get legendre polynomials for MDFF:
       call cpl0 (cosleg,pleg,ljmax+1)
!      now I have legendres from 0 to ljmax indexed from 1 to ljmax+1        plegqq(iq,iqq,:)=pleg(:)
       plegqq(iq,iqq,:)=pleg(:)
    enddo
    enddo


	 if(mixdff.and.(imdff.eq.1)) then
	    iqmin=1
		iqmax=nq
		iqqmin=1
		iqqmax=nq
	 elseif(mixdff.and.(imdff.eq.2)) then
	    iqmin=1
		iqmax=1
		iqqmin=2
		iqqmax=2
	 elseif(.not.mixdff) then
	    iqmin=1
		iqmax=nq
		!set the others inside the loop
	 else 
	    call par_stop('What is this - invalid MDFF option in getgtrjas')
	 endif

!     In the old code, SUM{iq}      weight{iq}  spectrum{iq}
!     In the new code, SUM{iq,iqq}  weight{iq}  weight{iqq}  spectrum{iq,iqq}
      if(mixdff) then
	     !In this case, user input is ~ a beam electron wave function coefficient
	     qweights=qw 
	  else
	     !In this case, user input is ~ a beam electron density coefficient (square of wave function coefficient)
	     do iq=1,nq
		    qweights(iq)=cdsqrt(qw(iq))
		 enddo
	  endif




      do ii=1,indmax
         if (abs(indmap(ii)).eq. icalc) then
!
!     do trace over the final m_l(Greens function)
!     however, this is non zero only for m_l=mjinit-m_s 
!     so do trace over minit
!
            traceres=0.0d0
            do mjinit=-jinit,jinit,2
               traceres=traceres+hbmat(isp,ii,mjinit) *hbmat(isp,ii,mjinit)
            end do
!
!     i**lj*(-i)**lj disappears 
!


	       do iq=iqmin,iqmax 
		      if(.not.mixdff) then
		         iqqmin=iq; iqqmax=iq !diagonal element  - !KJ
		      endif
		      do iqq=iqqmin,iqqmax
	 
            lj=ljind(ii)
            if (imode.eq.1) then 
               aa=-dcmplx(0.0d0,1.0d0)*xirflj(lj)*xirflj(lj) *   plegqq(iq,iqq,ljind(ii)+1)  !KJ pleg(lj+1) !KJ pleg turns it into MDFF
            else if (imode .eq.2) then
               aa=xirflj(lj)*  plegqq(iq,iqq,ljind(ii)+1)  !KJ   pleg(lj+1) !KJ again
            else
               write(6,*) "Something wrong in specupd"
               stop
            end if

            xsec(lj)=xsec(lj)-aa*traceres *qweights(iq) *qweights(iqq) !KJ qw !KJ qw for list of q-vectors
            if (imode.eq.1) then
               xnorm=xnorm+dble(xirflj(lj)*dconjg(xirflj(lj))) /dble(2*lj+1) *qweights(iq) *qweights(iqq) !KJ *qw !KJ again             
            end if
			
			   enddo
			enddo
			
         end if ! abs(indmap(ii)) == icalc
      end do ! ii

      return 
      end
