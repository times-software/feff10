!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_nrixs.f90,v $:
! $Revision: 1.10 $
! $Author: jorissen $
! $Date: 2012/03/22 18:59:18 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		MODULE nrixs_inp
		!KJ a module containing a few NRIXS working variables.  I may yet move this baby around.
		!Note: all INPUT is in m_inpmodules.f90, as it should be.  This nrixs_inp module is a WORKSPACE.

		use global_inp,only: le2,ldecmx,xivnorm,do_nrixs
		implicit none
		double precision qtrans
		integer jmax,kfinmax,jinit,indmax,ljmax,ldecs
        integer,allocatable:: lgind(:),kind(:),lind(:),ljind(:)
		

		CONTAINS

		subroutine nrixs_init
			! Aleksi's 		  call getsizes(le2,jinit,jmax,kfinmax)
			! getsizes is inserted here, it's not much anyway
			!KJ NOTE THAT WE DO NOT SAVE LINIT.  Aleksi has a bunch of routines where he uses the double of certain indices,
			! including lj and linit, presumably to make some arithmetic easier, I'm not sure.  But makes it a little dangerous
			! to make the variable global by saving it in the module ...
			use dimsmod,only: ltot
			integer ihole,kinit,linit,i1,i2,i3,i4,jfin,lfin,afin,ainit
			integer lj !lj is just a loop variable here
            integer ll,jspin,ient,ind,indt,kfin


            ljmax=abs(le2)  !KJ passed in a messy way, copied from Aleksi but would like to make cleaner ... - corresponds to "LJMAX lj" card
			ldecs=iabs(ldecmx) !KJ similar 12-2011
!          START of getsizes routine
	!c		get ihole
			open (unit=99, file='pot.bin', status='old')
			read(99,'(5(1x,i4))') i1, i2, i3, i4, ihole
			close(99)
	!c	    get kinit
			call setkap(ihole, kinit, linit)
			jinit=2*abs(kinit)-1
			jmax=abs(2*le2)+jinit
			kfinmax=2*(2*abs(kinit)+1)*(abs(le2)+1)*(jinit+1)
!          END of getsizes routine

            allocate(lgind(kfinmax),kind(kfinmax),lind(kfinmax),ljind(kfinmax))
			lgind=0 ; kind(:)=0 ; lind=0 ; ljind=0
			!!! next line taken from rexsphjas.f - not from getsizes
			qtrans = xivnorm
			if(do_nrixs.eq.0) kfinmax=8  !KJ THIS VERY IMPORTANT STATEMENT makes sure the dimension of arrays is the good old-fashioned 8 if nrixs is not used.

!           Next section of code taken from bcoefjas.f90 and identically calclbcoef.f90.
!           Need to put it here since we don't call bcoefjas in FMS anymore (now in MKGTR) but still need some of the arrays.
!           Note that "mind" and "jind" are working arrays local to bcoefjas/calclbcoef and therefore stay there.

!			get the factor of 2 into play
			ient=2
			jspin=1

			if (kinit.gt.0) then
				linit=jinit+1
				ainit=-1
			else
				linit=jinit-1
				ainit=1
			end if

		!     Now start figuring out the possible kfin. (Notice again lj=2*lj, and ljmax=2*ljmax)
			  ind=0
			  do lj=0,abs(ljmax) !KJ added abs
		!     triangle
				 do jfin= max(abs(2*lj-jinit),1), min(2*lj+jinit,jmax), 2
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
					kfin=-(jfin+1)*afin
					indt=0
					if (lfin.le.2*ltot) then 
					   ind=ind+1
					   indt=indt+1
					   kind(ind)=kfin
					   lgind(ind)=lfin
					   ljind(ind)=lj
					   lind(ind)=lfin
					end if
				 end do                 ! jfin
			  end do                    ! lj
			  indmax=ind
			  if (indmax.gt.kfinmax .and. do_nrixs .eq. 1) stop "Error : indmax > kfinmax"

		!     get rid of the factor of 2, in lgind,lind,and kind
			  do ll=1,indmax
				 lgind(ll)=lgind(ll)/2
				 lind(ll)=lind(ll)/2
				 kind(ll)=kind(ll)/2
			  end do 

            return
		end subroutine nrixs_init

		end module

