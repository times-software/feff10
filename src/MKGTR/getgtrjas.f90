!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: getgtrjas.f90,v $:
! $Revision: 1.9 $
! $Author: jorissen $
! $Date: 2011/11/29 00:04:07 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getgtrjas

!     Calculates FMS contribution to absorption
!     uses Bruce Ravel subroutine to do FMS in self-consistency loop
!     notice that it can do FMS with polarization dependence and
!     always include l-->l-1 channel.
!     written by alexei ankudinov 06.1997

      use dimsmod, only: nphx=>nphu, ltot, nspx=>nspu, nex, lx
	  use constants
      use IOMod
	  use par
      use eels_inp,only: ipmin,ipmax,ipstep
      use fms_inp,rclust=>rfms2,sigma2=>sig2g
      use global_inp,only: ipol,ispin,le2,angks,ptz,ldecmx,xivnorm,elpty,cosmdff,qw,qtrig,qaverage,nq,imdff,mixdff
      use atoms_inp,only: nat,iphat,ratdbl=>rat
      use nrixs_inp,only: jmax,kfinmax,jinit,indmax,kiind=>kind,lgind,ljind,lind,ljmax

      implicit none
      integer,parameter :: npadx = 8, iblock = 1

      integer ne, ne1, ne3, ihole  ! ,nph  !KJ 7-09 nph is now a module variable.  Will be reset in rexsph below
      integer iz(0:nphx)

!     work space
      complex*16, allocatable :: ph(:,:,:,:)
      integer, allocatable    :: lmax(:,:)
!     complex energy grid emg is decomposed into em and eref to have the same structure in phase.bin
      complex*16, allocatable  :: em(:), eref(:,:)
      character*6, allocatable  :: potlbl(:)
!     fms staff
      integer, allocatable :: indmap(:),map(:)
      complex, allocatable :: gtr(:),gtrloc(:),gtrl(:,:,:),gtrlloc(:,:,:)
      real*8, allocatable :: clbcoef(:,:,:,:),hbmat(:,:,:),res(:,:)
      complex*16, allocatable :: rkk(:,:,:,:),dum(:),gdummy(:,:)
      complex, allocatable :: gg(:,:,:,:)
      complex*16 hbrkk1,hbrkk2
      integer ix1,ix2,mj1,mj2,ig1,ig2,lg1,lg2
      integer ims1,ims2,jfin1,jfin2
      integer  nsp, ie
      character*25 innerform
      complex*16 itolj1,itolj2,itolg1,itolg2
      integer ilm1,ilm2,jlmax,mjlmax,j1ind,m1ind
	  real*8 rnrmav,xmu,edge,aa,bb
	  integer ik0,lmaxp1,kinit,linit,minit,nphtmp,iph,l1,l2,ix,is1,m1,isize,is2,m2,ios
	  integer iq,iqq,iqmin,iqmax,iqqmin,iqqmax
!      save  !KJ wtf??
      real*8 pleg(ljmax+1)  !APS legendre polynomials at argument=cosleg from l=0 to ljmax
	  real*8 qcst(nq),qsnt(nq),qcsf(nq),qsnf(nq),beta(nq),cosleg
	  complex pha(nq),ggrot(nspx*(lx+1)**2,nspx*(lx+1)**2)
	  real*8 plegqq(nq,nq,ljmax+1) !KJ legendre polynomials for all q,q' pairs
	  complex*16 qweights(nq)


  ! Allocate local variables
    allocate(ph(nex, -ltot:ltot, nspx, 0:nphx))
    allocate(lmax(nex, 0:nphx))
    allocate(em(nex), eref(nex, nspx))
    allocate(potlbl(0:nphx))
    allocate(rkk(nex,nq,kfinmax,nspx))
    allocate(dum(nex))
    allocate(gg(nspx*(lx+1)**2, nspx*(lx+1)**2, 0:nphx, nex))
    allocate(map(nex),gtr(nex),gtrloc(nex),gtrl(0:(iabs(ldecmx)),0:(iabs(ldecmx)),nex))
	allocate(gtrlloc(0:(iabs(ldecmx)),0:(iabs(ldecmx)),nex),gdummy(0:(iabs(ldecmx)),0:(iabs(ldecmx))))
	allocate(clbcoef(1:(2*lx+2),1:(lx+2),0:1,0:lx),hbmat(0:1,kfinmax,-jinit:jinit))
    allocate(indmap(kfinmax))
    allocate(res(0:1,-jinit:jinit))

    jlmax=lx+2
    mjlmax=2*jlmax-2

    gtr(:) = 0
	res(:,:)=0.0d0
	gtrl(:,:,:)=cmplx(0.0d0,0.0d0)

!   need less data than rphbin.f provides, also dimensions of ph array are different.
    call rdxsphjas (ne, ne1, ne3, nph, ihole, rnrmav, xmu, edge, ik0, em, eref, iz, potlbl, ph, rkk, lmax, lmaxp1)  !KJ 7-09 moved kfinmax,indmax to module nrixs_inp
    call setkap (ihole, kinit, linit)

	isize = 0
    nsp = 1
    if (abs(ispin).eq.1 ) nsp = nspx
!    check also nsp
     if (nsp.gt.2) stop "nsp too large"


!   KJ : put Adam's code in a loop
    do iq=1,nq
	do iqq=1,nq
	   cosleg=cosmdff(iq,iqq)
!      APS get legendre polynomials for MDFF:
       call cpl0 (cosleg,pleg,ljmax+1)
       plegqq(iq,iqq,:)=pleg(:)
    enddo
    enddo
	

!     read in the angles of the rotations for the different q-vectors
!     KJ : this used to be in a subroutine "angread" ; I've substituted the code below.  3-2011
      if(.not. qaverage ) then !.and. nq.gt.1) then   !KJ if (nq.gt.0) then
	     qcst=qtrig(:,1)
		 qsnt=qtrig(:,2)
		 qcsf=qtrig(:,3)
		 qsnf=qtrig(:,4)
!         call angread(nqloc,pha,beta,qweights)
         do iq=1,nq
		    pha(iq)=cmplx(qcsf(iq),qsnf(iq))
			pha(iq)=conjg(pha(iq))
			beta(iq)=atan2(qsnt(iq),qcst(iq))
!			if(beta(iq).lt.0.0d0 .or. beta(iq).gt.pi) stop 'something wrong with angles in getgtrjas'
		 enddo
      else
         pha(:)=cmplx(1.0d0,0.0d0)
         beta(:)=0.0d0
         qweights(:)=1.0d0
      end if


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


!   get hbmat
    do minit=-jinit,jinit,2
       if (elpty.ge.0) then
!	      this is the normal case
          call bcoefjas(kinit,minit,ltot,hbmat(0,1,minit))
       else
!	      this is the spherical averaging case 
          call calclbcoef(lx,jlmax,mjlmax,clbcoef)
       end if
    end do
        
  ! Read gg
    DO ie = 1, ne
		IF(.FALSE.) THEN
			nphtmp = nph
		ELSE
			nphtmp = 0
		END IF
		DO iph = 0, nphtmp
			CALL Read2D('gg.bin', gg(1:nspx*(lmaxph(iph)+1)**2,1:nspx*(lmaxph(iph)+1)**2,iph,ie), L1, L2)
			! Check that bounds are correct.
			IF(L1.ne.nspx*(lmaxph(iph)+1)**2.or.L2.ne.nspx*(lmaxph(iph)+1)**2) CALL Error('Error when reading gg.bin')
		END DO
	END DO
	CALL CloseFl('gg.bin')

!	 figure out what calculations are needed
     do ix=1,indmax
        if (lgind(ix).le.lx) then 
           indmap(ix)=1
        else
           indmap(ix)=0
        end if
     end do
	 
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



	 do ie =  1,ne 
	 
	 
	 do iq=iqmin,iqmax 
	 
	    !KJ rotate the green's function matrix:
!		if(nq.eq.1)then
!	    call rotgmatrix(0,elpty,pha(iq),beta(iq),qweights(iq),nsp,nspx,lx,gg(:,:,0,ie),ggrot)  !KJ qweights not used in current version		
!		else
	    call rotgmatrix(nq,elpty,pha(iq),beta(iq),qweights(iq),nsp,nspx,lx,gg(:,:,0,ie),ggrot)  !KJ qweights not used in current version
!	    endif
	 
	 
          do ix1=1,indmax
             jfin1=2*abs(kiind(ix1))-1
             lg1=lgind(ix1)
             itolj1=dcmplx(0.0d0,1.0d0)**ljind(ix1)
             itolg1=dcmplx(0.0d0,1.0d0)**lg1
             if(indmap(ix1).gt.0 ) then 
                if (elpty.ge.0.0d0) then 

!				not spherical average
                do mj1=-jinit,jinit,2
                   do ims1=1,nsp
                      is1=ims1-1
                      m1=(mj1-(2*is1-1))/2
                      ilm1=lg1*(lg1+1)+m1+1
                      if (abs(m1).le.lg1 .and. abs(mj1).le.jfin1) then
                         ig1=nsp*(lg1*lg1+lg1)+nsp*m1+ims1
                         hbrkk1=hbmat(is1,ix1,mj1)*rkk(ie,iq,ix1,ims1) *qweights(iq) !KJ added weight
						 if(.not.mixdff) then
						    iqqmin=iq; iqqmax=iq !diagonal element  - !KJ
						 endif
						 do iqq=iqqmin,iqqmax
                         do ix2=1,indmax
                            lg2=lgind(ix2)
                            itolj2=dcmplx(0.0d0,1.0d0)**ljind(ix2)
                            itolg2=dcmplx(0.0d0,1.0d0)**lg2
                            jfin2=2*abs(kiind(ix2))-1
                            if (indmap(ix2).gt.0) then
                               do ims2=1,nsp
									is2=ims2-1
									mj2=mj1
									m2=(mj2-(2*is2-1))/2
									ilm2=lg2*(lg2+1)+m2+1
									if (abs(m2).le.lg2 .and. abs(mj2).le.jfin2) then
							!		write(*,*) 'ix1,ix2',ix1,ix2,'hbrkk1,hbrkk1,ggrot,gg',hbrkk1,hbrkk2,ggrot(ig2,ig1),gg(ig2,ig1,0,ie)
										ig2=nsp*(lg2*lg2+lg2)+nsp*m2+ims2
										hbrkk2=hbmat(is2,ix2,mj2)*rkk(ie,iqq,ix2,ims2)  *qweights(iqq) !KJ added weight !KJ iq -> iqq
!										gtr(ie)=gtr(ie)+gg(ig2,ig1,0,ie)*hbrkk1*hbrkk2*itolj2*dconjg(itolj1)*dconjg(itolg1)*itolg2
										gtr(ie)=gtr(ie)+ggrot(ig2,ig1)*hbrkk1*hbrkk2*itolj2*dconjg(itolj1)*dconjg(itolg1)*itolg2*plegqq(iq,iqq,ljind(ix1)+1) !KJ
										if (lg1.le.ldecmx .and. lg2.le.ldecmx) &
											gtrl(lg2,lg1,ie)=gtrl(lg2,lg1,ie)+ggrot(ig2,ig1)*hbrkk1*hbrkk2* &
        itolj2*dconjg(itolj1)*dconjg(itolg1)*itolg2*plegqq(iq,iqq,ljind(ix1)+1) !KJ
!											gtrl(lg2,lg1,ie)=gtrl(lg2,lg1,ie)+gg(ig2,ig1,0,ie)*hbrkk1*hbrkk2*itolj2*dconjg(itolj1)*dconjg(itolg1)*itolg2
                                     end if !if (abs(mj2))..
                               end do !ims2
                            end if ! if (indmap(ix2)... 
                         enddo ! ix2
						 enddo !iqq
                      end if ! if (abs(mj1)...)
                   end do ! ims1
                end do ! mj1

                else
!				this is the spherically averaged case

                   j1ind=(jfin1+1)/2
                   do mj1=-jfin1,jfin1,2
                      do ims1=1,nsp
						is1=ims1-1
						m1=(mj1-(2*is1-1))/2
						if (abs(m1).le.lg1) then 
							ig1=nsp*(lg1*lg1+lg1)+nsp*m1+ims1
							m1ind=(mj1+jfin1)/2+1
							do ims2=1,nsp
								is2=ims2-1
								m2=(mj1-(2*is2-1))/2
								ig2=nsp*(lg1*lg1+lg1)+nsp*m2+ims2
						 if(.not.mixdff) then
						    iqqmin=iq; iqqmax=iq !diagonal element  - !KJ
						 endif
						 do iqq=iqqmin,iqqmax
								
								hbrkk1=rkk(ie,iq,ix1,ims1)*rkk(ie,iqq,ix1,ims2) *qweights(iq)*qweights(iqq) !KJ added weights !KJ iq->iqq
								hbrkk1= hbrkk1*clbcoef(m1ind,j1ind,is1,lg1)
								hbrkk1= hbrkk1*clbcoef(m1ind,j1ind,is2,lg1)
								gtr(ie)=gtr(ie)+ggrot(ig1,ig2)*hbrkk1/dble(2*ljind(ix1)+1) *plegqq(iq,iqq,ljind(ix1)+1) !KJ
								if (lg1.le.ldecmx) then
									gtrl(lg1,lg1,ie)=gtrl(lg1,lg1,ie)+ggrot(ig1,ig2)*hbrkk1/dble(2*ljind(ix1)+1) *plegqq(iq,iqq,ljind(ix1)+1) !KJ
								end if
						  enddo !iqq
							end do
						end if
                      end do
                   end do
                end if             ! if elpty.ge. ...
             end if                ! if (indmap(x1)...
          end do ! ix1
		  enddo ! iq
     enddo  ! do ie=1,ne


!KJ NRIXS OUTPUT :
!KJ   - fms.bin same as for regular feff, contains gtr
!KJ   - fmsl.bin new file, contains gtrl
!KJ   - gtrl.dat new file, contains gtrl in readable format
!     write fms.bin
      if (master) then
         open (unit=1, file='fms.bin', status='unknown', iostat=ios)
!        write title line
         write(1,105) rclust*bohr
 105     format('FMS rfms=', f7.4)
         write(1, 110) ne, ne1, ne3,  nph, npadx
 110     format(5(1x,i3))
         do 120 ie = 1, ne
           aa = dble ( real ( gtr(ie) ) )
           bb = dble (aimag ( gtr(ie) ) )
           dum(ie) = dcmplx (aa,bb)
 120     continue
         call wrpadx(1, npadx, dum(1), ne)
         close (unit=1)

         open(99,file='gtr.dat',form='formatted',status='unknown')
         do ie=1,ne
            write(99,'(10f13.6)') em(ie),gtr(ie)
         enddo
         close(99)

         if (ldecmx.ge.0) then 
            !write(innerform,'(''(I5,'',I2,''e18.8)'')') 19
            innerform='(I5,100e18.8)'
            open(unit=66,file='gtrl.dat',form='formatted',status='unknown')
            rewind(66)
            open(unit=7,file='fmsl.bin',status='unknown',iostat=ios)
            call chopen(ios,'fmsll.bin','fmstotjas')
            do ie=1,ne 
               do ilm1=0,ldecmx
                  do ilm2=0,ldecmx
                     aa = dble ( real ( gtrl(ilm2,ilm1,ie) ) )
                     bb = dble (aimag ( gtrl(ilm2,ilm1,ie) ) )    
                     gdummy(ilm2,ilm1) = dcmplx(aa,bb)
                  end do
               end do
               call wrpadx(7,npadx,gdummy,(ldecmx+1)*(ldecmx+1)) 
               write(66,innerform) ie,dble(em(ie)), &
                    ((real(gtrl(ilm2,ilm1,ie)),ilm2=0,2),ilm1=0,ldecmx), &
                    ((imag(gtrl(ilm2,ilm1,ie)),ilm2=0,2),ilm1=0,ldecmx)
            end do
            close(7)
            close(66)
         end if
      endif

  !   Deallocate local variables
      deallocate(ph,lmax,em,eref,potlbl,gg,rkk,dum)
      deallocate(map,gtr,gtrloc,gtrl,gtrlloc,gdummy,res)
	  deallocate(clbcoef,hbmat,indmap)


      return
      end
