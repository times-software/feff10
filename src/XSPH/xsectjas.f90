!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: xsectjas.f90,v $:
! $Revision: 1.18 $
! $Author: jorissen $
! $Date: 2012/02/03 00:45:55 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine xsectjas (ipr2, dx, x0, ri, ne, ne1, ik0, em, edge, &
          ihole, emu, corr, dgc0, dpc0, jnew, dgcx0, dpcx0, &
          ixc, lreal, rmt, rnrm, xmu, &
          vi0, iPl, NPoles, Eps0, EGap, & !iPl new JK, NPoles and Eps0 new JJK 3/9/2010
		  gamach, &
          vtot, vvalgs, edens, dmag, edenvl, &
          dgcn, dpcn, adgc, adpc, xsec, xsnorm, rkk, &
          iz, xion, iunf, xnval, iorb, l2lp, &
          ipol, ispin, ljmaxdummy, angks, ptz, iph) !KJ iph      !, &
          !qtrans,kfinmax,jmax,jinit,indmax, &
          !ldecmx)

!     right now the same self-energy is used for calculation
!     of the central atom part (xsec) and dipole m.e. for
!     scattering (rkk). You may want to run xsect separately
!     for xsec and for rkk, if you want to use different self-energy
!     for central and scattering parts.  ala. fix later


!     INPUT
!     dx, x0, ri(nr)
!     Loucks r-grid, ri=exp((i-1)*dx-x0)
!     ne, em(ne)   number of energy points, real energy grid
!     edge         chemical potential (energy for k=0)
!     ihole        hole code
!     emu          position of chemical potential in absorption specrum
!     dgc0(nr)     dirac upper component, ground state hole orbital
!     dpc0(nr)     dirac lower component, ground state hole orbital
!     ixc          0  Hedin-Lunqist + const real & imag part
!     1  Dirac-Hara + const real & imag part
!     2  ground state + const real & imag part
!     3  Dirac-Hara + HL imag part + const real & imag part
!     5  Dirac-Fock exchange with core electrons +
!     ixc=0 for valence electron density
!     lreal        logical, true for real phase shifts only
!     rmt          r muffin tin
!     xmu          fermi level
!     vi0          const imag part to add to complex potential
!     gamach       core hole lifetime
!     vtot(nr)     total potential, including gsxc, final state
!     edens(nr)    density, hole orbital, final state
!     dmag(251)     density magnetization
!     edenvl      valence charge density
!     dgcn(dpcn)   large (small) dirac components for central atom
!     adgc(adpc)   their development coefficients
!     
!     OUTPUT
!     xsec(ne)    atomic absorption cross section to multiply \chi
!     (atomic background for XMCD)
!     xsnorm(ne)  atomic  absorption cross section (norm for XMCD)
!     rkk(ne, 8)  normalized reduced matrix elements for construction
!     of termination matrix in genfmt.
      use dimsmod, only: nex, nrptx, ltot, MxPole, nspx=>nspu, nphx=>nphu
	  use IOMOD
	  use SelfEnergyMod
      use constants
	  use global_inp,only: nq,qweights=>qw,qnorms=>qn,qaverage 
	  use nrixs_inp,kiind=>kind
      implicit none !double precision (a-h, o-z)
      integer, parameter :: npadx=8

!      integer ljmax,jinit,ipr2,iunf,ipol,ne1,ik0,ios,ihole,kinit,linit,imt,jri,jri1,inrm,jnrm,i,lreal,index,ixc,minit,ispin,ilast,jnew,isp,ie,k1
      integer ipr2,iunf,ipol,ne1,ik0,ios,ihole,kinit,linit,imt,jri,jri1,inrm,jnrm,i,lreal,index,ixc,minit,ispin,ilast,jnew,isp,ie,k1
	  integer ifirst,ne,iph,ncycle,ikap,lfin,iold,ic3,irr,iz,ilp,ll,ind,kdif,l2lp
	  real*8 edge,emu,gamach,rmt,x0,dx,rnrm,xinorm,del,xmu,omega,xk0,xion,sign,prefac,xnorm,corr,vi0,angks


!      double precision qtrans
!      integer indmax,kfinmax,jmax
      complex*16 ptz(-1:1, -1:1)
      complex*16 em(nex)
	  !KJ 11-2011 In an attempt to make the NRIXS routines less vulnerable to either stack overflow or compiler optimization stack issues,
	  ! I've turned the following arrays into allocatable arrays (so they won't go on the stack).	  
      real*8 ri(nrptx), vtot(nrptx), edens(nrptx),dmag(nrptx)
      real*8 dgc0(nrptx), dpc0(nrptx), vvalgs(nrptx), edenvl(nrptx)
      real*8 dgcx0(nrptx), dpcx0(nrptx)
      real*8 dgcn(nrptx,30), dpcn(nrptx,30)
!KJ
	  
      real*8 adgc(10,30), adpc(10,30), xnval(30)
	  integer iorb(-4:3)
      complex*16 rkk(nex, nq, kfinmax), xsec(nex),xsecl(0:ljmax) !KJ added nq
      complex*16 atomxsec(kfinmax,nex)

      complex*16 xseclg(0:ljmax),lgsum, ortcor(0:ljmax,nq)
!     double precision bmat(-lx:lx,0:1,8, -lx:lx,0:1,8)
!     additional data needed for relativistic version

!     
!     Aleksi changed
!     
!        call xsectjas (ipr2, dxnew, x0, ri, ne, ne1, ik0, em, edge, &
!             ihole, emu, corr, dgcx, dpcx, jnew,dgcx0,dpcx0, & !dgcx0,dpcx0 are new for jas
!             ixc0, lreal, rmt(0), rnrm(0), xmu, vi0, &
!             gamach, vtotph, vvalph, rhoph, dmagx, rhphvl, & 
!             dgcn, dpcn, adgc(1,1,iph), adpc(1,1,iph), xsec(1,isp), &
!             xsnorm(1,isp), rkk(1,1,isp),iz(0), xion(0), iunf, &
!             xnval(1,iph), iorb(-4,iph), l2lp, &
!             ipol, ispinp, abs(le2), angks,ptz, & ! abs is new
!             qtrans,kfinmax,jmax,jinit,indmax,ldecmx) ! nrixs vars

!     double precision hbmat(0:1,kfinmax)
!     USE THESE WHEN THE MINIT IS NOT SET ARTIFICIALLY
!     
      real*8 hbmat(0:1,kfinmax,-jinit:jinit)
!     
!     
!     if want all the minit

      !integer kiind(kfinmax)!, lgind(kfinmax),ljind(kfinmax)  !KJ renamed kind to kiind to avoid conflicts with reserved name as used in modules
      !integer lind(kfinmax)
      integer indmap(kfinmax),indcalc(kfinmax,3)
      integer ljneeded(0:ljmax)
      real*8 qjbess(nrptx,0:ljmax,nq)
      integer ljj,ncalc,icalc,know,ljnow,lnow,ii,ic3no
      character*25 innerform
    
      real*8 xsnorm(nex)

      real*8 xp(nrptx), xq(nrptx)

!     work space for xcpot
      real*8 vxcrmu(nrptx), vxcimu(nrptx), gsrel(nrptx)
      real*8 vvxcrm(nrptx), vvxcim(nrptx)
      DOUBLE PRECISION WpCorr(MxPole), Gamma(MxPole), AmpFac(MxPole)
      integer iPl,NPoles,NData
	  real*8,allocatable :: Energy(:), Loss(:)
	  real*8 EGap,Eps0
      character(LEN=10) ColumnLabels(20)
	  
!     work space for fovrg
      complex*16 p(nrptx), q(nrptx), pn(nrptx), qn(nrptx)
!     storage for calculation of cross term (SPIN 1 only)
      complex*16 xrcold(nrptx,0:1) , xncold(nrptx,0:1)

      complex*16  p2, ck, xkmt, xkmtp
      complex*16  pu, qu, dum1, factor
      complex*16  xfnorm, xirf
      complex*16  temp, aa, bb, cc, phold,ph0
      complex*16  rkk1(0:ljmax),rkk0(0:ljmax)
      complex*16  phx(kfinmax)
      complex*16  xirflj(0:ljmax,nq)
      complex*16  eref, xm1, xm2, xm3, xm4

      complex*16  jl,jlp1,nl,nlp1
      complex*16 jlall(ltot+2), nlall(ltot+2)
!KJ 11-2011 In an attempt to reduce the stack footprint of the NRIXS routines, the following big guys are now dynamically allocated:
      complex*16,allocatable :: v(:), vval(:), xrc(:,:,:), xnc(:,:,:)	  
!KJ      complex*16  v(nrptx), vval(nrptx)
!KJ      complex*16  xrc(nrptx,0:ljmax,nq), xnc(nrptx,0:ljmax,nq)
      character*512 slog
      logical ltrace
!     nesvi:  
      complex*16 xrhoce(nex), xrhopr(nex), chia(nex), cchi(nex)
      real*8 omega1(nex), bf(0:2, nrptx)

      real*8 pat(nrptx),qat(nrptx)
      complex*16 intr(nrptx),var(nrptx) 

      external besjn
      real*8 w1,theta,w2,qnew,xsdum
      integer iorthg,iq,llim,lfin0,ljmaxdummy




!KJ 11-2011 Now allocate arrays that used to be dimensioned at compile time (making them end up on the stack) :
      allocate (v(nrptx), vval(nrptx), xrc(nrptx,0:ljmax,nq), xnc(nrptx,0:ljmax,nq))



!     set imt and jri (use general Loucks grid)
!     rmt is between imt and jri (see function ii(r) in file xx.f)
      imt = (log(rmt) + x0) / dx  +  1
      jri = imt+1
      jri1 = jri+1
      if (jri1 .gt. nrptx)  call par_stop('jri .gt. nrptx in phase')

      WpCorr(:) = 0.d0
      Gamma(:) = 0.d0
      AmpFac(:) = 0.d0
!KJ 3-2011 : I blindly copied the following section from xsect ; hope it works!
!     Josh - if PLASMON card is set, and using HL exc,
!          - read loss function from loss.dat and find poles etc.
      ! Josh - Changing to read directly from loss.dat : 3/9/2010
      IF ( (iPl.gt.0).and.(ixc.eq.0) ) THEN
	     write(6,*) 'doing the self-energy'
		 write(6,*) 'ipl,ixc,ndata',ipl,ixc,ndata
         WpCorr(:) = -1.d30
         CALL OpenFl('loss.dat', FileStatus = 'OLD')
         NData = NumberOfLines('loss.dat')

         ALLOCATE(Energy(NData), Loss(NData))

         CALL ReadArrayData('loss.dat', Double1 = Energy, Double2 = Loss)
  
         !NPoles = 100 ! 100 poles should be enough for anything. Can add a card for this later.
         !Eps0 = -2.d0 ! This will not do anything. Will add a card to set Eps0 later.
         CALL MkExc(Energy, Loss, Eps0, WpCorr, AmpFac, NPoles)

         Gamma(:)  = Gamma(:)/hart
         WpCorr(:) = (WpCorr(:)/hart) / SQRT(3.d0/((3 / (4*pi*edens(jri+1))) ** third)**3)
         CALL CloseFl('loss.dat')
         DEALLOCATE(Energy,Loss)
      END IF

      ! Write wp as calculated from density to sigma.dat
      ColumnLabels(:) = '          '
      ColumnLabels(1) = 'Rs        '
      ColumnLabels(2) = 'wp (eV)   '
      CALL WriteData('mpse.dat',                                                       &
           &  Double1 = dble((3 / (4*pi*edens(jri+1))) ** third),                          &
           &  Double2 = dble(SQRT(3.d0/((3 / (4*pi*edens(jri+1))) ** third)**3)*hart),     &
           &  ColumnLabels = ColumnLabels, WriteDataInHeader = .TRUE.,               &
           &  Headers = (/ 'This file contains information about the self-energy.' /))

!     Josh END
!!KJ




      lfin0=20
!      w1=29300.0d0/27.21140d0
       w2=9890.00d0/27.21140d0
       theta=165.0d0/180.0d0*4.0d0*atan(1.0d0)
!      theta=148.0d0/180.0d0*4.0d0*atan(1.0d0)
!
!     innerform
!
      ii=ljmax+1
      ii=2*ii+3
      if (ii.lt. 10) then 
         write(innerform,'(''('',I1,''e18.8)'')') ii
      else if (ii.lt.100) then
         write(innerform,'(''('',I2,''e18.8)'')') ii
      else
         write(6,*) "xsecl write fails in xsectjas"
         stop
      end if
      open(unit=66,file='xsecl.dat',form='formatted',status='unknown')
      rewind(66)
      write(66,*) ne1,ik0,edge,emu,gamach
      open(unit=77,file='xsecl2.dat',form='formatted',status='unknown')
      rewind(77)
      write(77,*) ne1,ik0,edge,emu,gamach
      if (ldecmx.ge.0) then 
         open (unit=7,file='xsecl.bin',status='unknown',iostat = ios)
         call chopen(ios,'xsecl.bin','xsectjas')
      end if


!     No spin-dependent calculations
!     
      ic3no=1



      call setkap(ihole, kinit, linit)


!     set imt and jri (use general Loucks grid)
!     rmt is between imt and jri (see function ii(r) in file xx.f)
      imt = (log(rmt) + x0) / dx  +  1
      jri = imt+1
      jri1 = jri+1
      if (jri1 .gt. nrptx)  stop 'jri .gt. nrptx in phase'

!     nesvi: define jnrm
      inrm = (log(rnrm) + x0) / dx + 1
      jnrm = inrm + 1

!     We'll need <i|i> later to normalize dipole matrix elements
!     <i|r|f>.  NB, dgc and dpc are r*wave_fn, so use '0' in somm to
!     get integral  psi**2 r**2 dr.
!     Square the dgc0 and dpc0 arrays before integrating.
!     <i|i> == xinorm.
!     dgc and dpc should be normalized <i|i>=1, check this here
      do i = 1, nrptx
         xp(i) = dpc0(i)**2
         xq(i) = dgc0(i)**2
      end do
!     nb, xinorm is used for exponent on input to somm
      xinorm = 2*linit + 2
      call somm (ri, xp, xq, dx, xinorm, 0, jnrm)
      del = abs (abs(xinorm) - 1.0d0)
      if (del .gt. 1.e-2) then
         write(slog,'(a,i8,1p2e13.5)') ' ihole, xinorm ', ihole , xinorm
         call wlog(slog)
!     if using real phase shifts, don't expect great results
         if (lreal.lt.2)  then
            call wlog(' There may be convergence problems.')
            call wlog(' Xinorm should be 1. If you set the RGRID, minor interpolation errors ')
            call wlog(' that will not affect final results may occur')
         end if
      end if

!     use ixc for testing
      index = ixc
!     Always use ground state self energy for xsection, quick fix
!     JJR, Jan 93
!     change for testing broadened plasmon pole 6/93
!     index = 2
!     ALA found that it is better to use index=ixc and real part of 
!     self-energy for atomic xsection. 12/96
!     ltrace = .true.
!     ltrunc = 2
!     call bcoefp(kinit, ltrunc, ltrace, ispin, angks, 
!     1           kiind, lind, bmat)
!     set spin index to use bmat

      do minit=-jinit,jinit,2
         call bcoefjas(kinit,minit,ltot,hbmat(0,1,minit))
      end do
      

!     Get all the k-l_{bessel} pairs ordered for 
!     optimal computation
!     
!     indcalc(i,1) gives kappa for a calculation
!     indcalc(i,2) gives the maximum l- for the bessel 
!     function in the radial matrix element 
!     

      call mincalc(kfinmax,indmax,kiind,lind,ljind,indcalc,indmap,ljj,ncalc)
!     
!     get all the bessel functions that are needed
!     

	  if (nq.le.0) then  !KJ fix later  !KJ later : this shouldn't happen anymore ...  qtrans=qnorms(1)
         call qbesselget(qtrans,nrptx,ljmax,ri,qjbess(1,0,1))
      else
	     do iq=1,nq
		    call qbesselget(qnorms(iq),nrptx,ljmax,ri,qjbess(1,0,iq))
	     enddo
	  endif

!
! get the correction to be used in radjas
!

      ilast = jnrm+6
      if (ilast.lt.jnew) ilast = jnew
      if (ilast.gt.nrptx) ilast = nrptx
      
      call getcorrection(nrptx,ri,dx,jinit,linit,dpc0,dgc0,ljmax,qjbess,ilast,ortcor)


      isp = 0
      if (ispin.eq.1) isp = nspx - 1
!     zero rkk and phx
	  rkk=dcmplx(0) !KJ
	  atomxsec=dcmplx(0)
	  phx=dcmplx(0)

      ifirst = 0
      do ie = 1, ne
	  
         if (ie.gt.ne1) then
            iorthg=1
         else
            iorthg=0
         end if
         iorthg=2
         iph = 0

         EGap=dble(0)
         call xcpot (iph, ie, index, lreal, ifirst, jri,em(ie), xmu,vtot, vvalgs, edens, dmag, edenvl, &
             eref, v, vval, iPl, WpCorr, Gamma, AmpFac, EGap, & ! for MPSE - I've just hacked with zeroes now FIX LATER
         &    vxcrmu, vxcimu, gsrel, vvxcrm, vvxcim)
!		set the method to calculate atomic cross section
!		p2 is (complex momentum)**2 referenced to energy dep xc
         p2 = em(ie) - eref
         ck = sqrt (2*p2 + (p2*alphfs)**2)
         xkmt = rmt * ck
!
!		Get the bessels here
!        
!         write(6,*) "Here we call it ",xkmt,ie, em(ie),ck,p2 
         call besjn (xkmt, jlall, nlall)
 
         if (mod(index,10) .lt. 5) then
            ncycle = 0
         else
!		fix later . may be ncycle can be less
            ncycle = 3
         end if
!
!		JAS is leaving these here. Basically they are for the "real"
!		experimental setup for XRS. The momentum transfer is function
!		of the scattering angle (theta), incoming energy (w1) and outgoing
!		energy (w2)  
!
		!if (nq.le.0) then  !KJ fix later
		if(qaverage) then  !KJ my fix
		
!			dble(em(ie))-edge+emu,
!			w2=w1-(dble(em(ie))-edge+emu)
!				-dble(em(ie))
!			w1=w2+dble(em(ie))-edge+emu
!				+dble(em(ie))
!
!			JAS: if one is doing "scanning" calculation uncomment
!			line below and comment the one below it.  
!			qnew=sqrt(w1*w1+w2*w2-2.0d0*w1*w2*cos(theta))*alphfs
			qnew=qtrans
			if (ie.gt.ne1) then 
!				For energies on the imaginary axis use the qtrans (this assumed to be the momentum transfer at the edge).
				call qbesselget(qtrans,nrptx,ljmax,ri,qjbess(1,0,1))
			else
!				For "real" energies use the new momentum transfer
				call qbesselget(qnew,nrptx,ljmax,ri,qjbess(1,0,1))
			end if
		endif !KJ fix later

         omega = (dble(em(ie)) - edge) + emu
         omega = max (omega, 0.1d0 / hart)
!		nesvi: add omega1(ie)- need it later
         omega1(ie) = omega

!		remember the bessel functions for multipole matrix elements
         xk0 = omega * alphfs
         ilast = jnrm+6
         if (ilast.lt.jnew) ilast = jnew
         if (ilast.gt.nrptx) ilast = nrptx
         xsnorm(ie) = 0 
         xsec(ie) = dcmplx(0.0d0,0.0d0)
         do ii=0,ljmax
            xsecl(ii)=dcmplx(0.0d0,0.0d0)
            xseclg(ii)=dcmplx(0.0d0,0.0d0)
         end do
         if (dble(em(ie)).ge.-10.d0) then 
            if (dimag(p2).gt.0.d0 .or. dble(p2).gt.0.d0) then 
               do icalc=1,ncalc
                  know=indcalc(icalc,1)
                  ikap=know
                  ljnow=indcalc(icalc,2)
                  call ljneeded0(ljmax,ljneeded,kfinmax,indmax,ljind,icalc,indmap,indcalc)
                  lnow = indcalc(icalc,3)
                  lfin = lnow
                  iold = 0
                  ic3=0
                  irr = -1
                  ilast = jnrm + 6
                  if (ilast.lt.jnew) ilast = jnew
                  call dfovrg ( ncycle, ikap, rmt, ilast, jri, p2, dx,ri, v, vval, dgcn, dpcn, adgc, adpc, &
                      xnval, pu, qu, p, q,iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph
                  

                  ilp = lfin - 1
                  if (ikap .lt. 0) ilp = lfin + 1
                  jl=jlall(lfin+1)
                  nl=nlall(lfin+1)
                  jlp1=jlall(ilp+1)
                  nlp1=nlall(ilp+1)
                  call phamp(rmt,pu,qu, ck, jl,nl,jlp1,nlp1, ikap, ph0,temp)
                  
                  sign = -1.0
                  if (ikap.gt.0) sign = 1.0
                  factor = ck*alphfs 
                  factor = sign * factor/(1+sqrt(1+factor**2))
                  dum1 = 1/ sqrt(1+factor**2)
                  xfnorm = 1 / temp *dum1
!					normalization factor
!					xfnorm = dum1*rmt*(jl*cos(delta) - nl*sin(delta))/ Rl(rmt)
!					dum1 is relativistic correction to normalization
                  
!					normalize regular solution
                  do i = 1,ilast
                     p(i)=p(i)*xfnorm
                     q(i)=q(i)*xfnorm
                  end do
				  do iq=1,nq
                  call radjas(1, kinit, dgc0, dpc0,dgcx0,dpcx0, ikap, p, q, &
                    pn, qn,ri, dx, ilast, iold, xrc(1,0,iq), xnc(1,0,iq), &
                    xrcold, xncold, xirflj(0,iq),ljmax,ljneeded, &
                    qjbess(1,0,iq),iorthg,ortcor(0,iq))
			      enddo
                  if (ic3.eq.0) then
!     
!					map the calculated results to the correct 
!					matrix elements and phases 
!     
                     
                     do ii=1,indmax
                        if (abs(indmap(ii)).eq.icalc) then 
                           phx(ii) = ph0
                           ll=ljind(ii)
                           if (lfin.le.5 .or. (ie.gt.ik0 .and. ie.le.ne1)) then 
                              rkk(ie,1:nq,ii)=xirflj(ll,1:nq) 
                           else
                              rkk(ie,:,ii)=dcmplx(0.0d0,0.0d0)
                              xirflj(ll,:)=dcmplx(0.0d0,0.0d0)
                           end if
                        end if
                     end do
!     
!					Update xsec and xnorm
!             
                     do iq=1,nq
                     call specupd(icalc,isp,kfinmax,indmax,indcalc,indmap,ljind,jinit,hbmat,ljmax,xirflj(0,iq),xsecl,xsnorm(ie),1) !KJ 2012 ,qweights(iq))
                     xsdum=0.d0
					 !xsdum=xsnorm(ie)  !KJ test
                     call specupdatom(icalc,isp,kfinmax,indmax,indcalc, &
                       indmap,ljind,jinit,hbmat,ljmax,xirflj(0,iq), &
                       atomxsec(1,ie),xsdum,1) !KJ 2012 ,qweights(iq))
                     !xsnorm(ie)=xsdum  !KJ test

                     
                     call specupdlg(icalc,isp,kfinmax,indmax,indcalc,indmap,lind,ljind,jinit,hbmat,ljmax,xirflj(0,iq),xseclg,1) !KJ 2012 ,qweights(iq))
                     enddo


                  elseif (iold.eq.1 .and. ic3no.eq.0) then
				     stop 'should not be here'
!                     do ll=0,jmax
!                        rkk1(ll)=xirflj(ll)
!                     end do
!                     phold = ph0
                  elseif (iold.eq.2) then
				     stop 'should not be here'
!                     do ll=0,jmax
!                        rkk0(ll)=xirflj(ll)
!                     end do
                  end if

!					get irregular solution and atomic cross-section xsec
!					find irregular solution

                  if(dimag(em(ie)).gt.0.d0) then
                     irr = 1
!					set pu, qu - initial condition for irregular solution 
                     pu = (nl*cos(ph0)+jl*sin(ph0)) *rmt * dum1
                     qu=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
                     
!					test on bessel functions
!					if (ikap.gt.0) print*,'test1',xkmt**2*(jl*nlp1-nl*jlp1)
                     
                     do i = 1, ilast
                        pn(i) = 0
                        qn(i) = 0
                     end do
                     call dfovrg (ncycle, ikap, rmt, ilast, jri, p2, dx, ri, v,vval, dgcn, dpcn, adgc, adpc, &
                         xnval, pu, qu, pn, qn,iz, ihole, xion, iunf, irr, ic3, iph) !KJ iph
!c feff82D changes after irr = 1 and call dfovrg
!c            set N- irregular solution , which is outside
!c            N=(nlp1*cos(ph0)+jlp1*sin(ph0))*factor *rmt * dum1
!c  
                     temp = exp(coni*ph0)
                     do i = 1, ilast
                        pn(i) = coni * p(i) - temp * pn(i)
                        qn(i) = coni * q(i) - temp * qn(i)
                     enddo
                  else
                     do i = 1, ilast
                        pn(i) = 0
                        qn(i) = 0
                     end do
                  end if

!				combine regular and irregular solution into the
!				central atom absorption coefficient xsec (mu = dimag(xsec))
!				thus for real energy dimag(xsec)=xsnorm
                  iorthg=2
				  do iq=1,nq
                  call radjas(2, kinit, dgc0, dpc0, dgcx0,dpcx0,ikap, p, q, &
                    pn, qn,ri, dx, ilast, iold, xrc(1,0,iq), xnc(1,0,iq), &
                    xrcold, xncold, xirflj(0,iq),ljmax,ljneeded, &
                    qjbess(1,0,iq),iorthg,ortcor(0,iq))
				  enddo
                  if (ic3.eq.0.and. (lfin.le.5 .or. (ie.gt.ik0 .and. ie.le.ne1)) .and. lfin.le.lfin0) then
				     do iq=1,nq
                     call specupd(icalc,isp,kfinmax,indmax,indcalc,indmap,ljind,jinit,hbmat,ljmax,xirflj(0,iq),xsecl,xsnorm(ie),2) !KJ 2012 ,qweights(iq))

                     call specupdatom(icalc,isp,kfinmax,indmax,indcalc,indmap,ljind,jinit, &
                          &           hbmat,ljmax,xirflj(0,iq),atomxsec(1,ie),xsdum,2) !KJ 2012 ,qweights(iq))
					enddo
!     
!     
!				Do not understand why xirf and not xirf*xirf (or xirf*dconjg(xirf))
!				correct this
!				xsec(ie) = xsec(ie) - xirf * bmat(0,isp,ind, 0,isp,ind)
!				also why is the xnorm not updated 
                  end if
                  
                  if (iold.gt.0 .and. ic3no.eq.0) then
                     write(6,*) "This part is not in use"
                     stop "Do not go here"
!
!THIS IS FAKE
!
                     ind = 0
!     calculate cross term contribution to XMCD
!     in both cases coupling between neighbors 
!     need to remove SO interaction (ic3=1) in order
!     to avoid unphysical peak in Gd XMCD. a.l. ankudinov
                     k1 = ind - 1
                     if (k1.ge.1 .and.k1.le.8) then
                        if (lgind(k1).eq.lgind(ind) .and. lgind(k1).gt.0) then
!     aa = exp( coni*(ph0 - phold))
!     bb = 1/aa
!     cc = - ( bmat(0,isp,k1, 0,isp,ind) +
!     1                    bmat(0,isp,ind, 0,isp,k1) ) / 2.d0
!     xsec(ie) = xsec(ie) - coni * rkk1 * rkk0 
!     1                    * (bb+aa) * cc
!     c              combine regular and irregular solution into the
!     c              central atom absorption coefficient (mu=dimag(xsec))
!     c              thus for real energy dimag(xsec)=xsnorm
!     call radjas (3, lls, bf, kinit, dgc0, dpc0, ikap, p, 
!     1                    q, pn, qn, ri, dx, ilast, iold, xrc, xnc, 
!     2                    xrcold, xncold, xirf)
!     xsec(ie) = xsec(ie) + xirf * cc * bb
!     
!     call radjas (4, lls, bf, kinit, dgc0, dpc0, ikap, p, 
!     1                    q, pn, qn, ri, dx, ilast, iold, xrc, xnc, 
!     2                    xrcold, xncold,xirf)
!     correct this
!     xsec(ie) = xsec(ie) + xirf * cc * aa
                        end if
                     end if
                  end if
!     c          end of |ispin=1| cross term calculations
                  
!     prepare for ic3=1 cross term calculations if needed
                  if (ic3.eq.0 .and. abs(ispin).eq.1 .and. ic3no.eq.0) then
                     iold = 0
                     if (ind.lt.8 .and. lgind(ind).gt.0) then
                        k1 = ind + 1
                        if (lgind(k1).eq.lgind(ind)) iold = 1
                     end if
                     if (ind.gt.1 .and. lgind(ind).gt.0) then
                        k1 = ind - 1
                        if (lgind(k1).eq.lgind(ind)) iold = 2
                     end if
!     need to remove SO interaction to calculate cross term
!     big effect for Gd XMCD calculations
                     if (iold.gt.0) then
!     repeat calculation for current kdif with SO turned off
                        ic3 = 1
                        write(6,*) "Do not want to go anywhere"
                        stop
!                        goto 100
                     end if
                  end if
                  
!     300        continue
               end do           ! icalc
            end if
         end if
!     350   continue

         if (omega.gt.0.0) then
!     prefac = (8 * pi / 3)  * alphfs * omega  -- nonrelativistic
!     relativistic is (for alpha form)
!     prefac = 4 * pi * alpinv / omega * bohr**2
!     
!     NRIXS prefac the (4*pi)^2 comes from 
!     bessel function expansion of e^iqr and 1/pi from 
!     writing the spectra in a im(1/(something+i\delta)) form
!
!            prefac=16.0d0*pi
!     ofcourse this time we do not have (4pi)^2
!
!     
!     Trying to make spectra normalized to "one electron"
!     as opposed to total contribution from a sub-level for a spin
!     Ofcourse this is easy to change back. 
!
!            prefac=1.0d0/pi
            prefac=2.0d0/pi/dble(jinit+1)

            do ii=0,ljmax
               xsec(ie)=xsec(ie)+xsecl(ii)
               xsecl(ii)=xsecl(ii)* prefac* 2.0d0*ck
               xseclg(ii)=xseclg(ii)* prefac* 2.0d0*ck
            end do
            do ii=1,kfinmax !KJ 7-09 I presume this was a bug - used to be : 0,kfinmax  !Aleksi's version had the same problem, by the way.
               atomxsec(ii,ie)=atomxsec(ii,ie)*prefac*2.0d0*ck
            end do
!
! making the xsnorm = 1.0*prefac*2.0d0*abs(ck)
! instead of xsnorm(ie) * prefac * 2.0d0*abs(ck) 
! where xsnorm depends on the qtrans
!

            xsnorm(ie) =  xsnorm(ie)*prefac* 2.0d0*ck
            xnorm= sqrt( xsnorm(ie) )
            xsec(ie) = xsec(ie) * prefac* 2.0d0*ck
!     put complex sqrt(prefactor) into reduced matrix elements rkk
            ck = sqrt ( prefac * (2.0d0*ck))
!     guarantee that we have the right root
            if (dimag(ck) .lt. 0) ck = -ck
!     add central atom phase shift here. 
            do kdif = 1 , indmax
			do iq=1,nq
               rkk(ie,iq,kdif)= rkk(ie,iq,kdif) * ck * exp(coni*phx(kdif))/xnorm
	        enddo
            enddo
         endif
!
!     here write out the xsecl
!
         lgsum=dcmplx(0.0d0,0.0d0)
         do ii=0,ljmax
            lgsum=lgsum+xseclg(ii)
         end do
         write(66,innerform) dble(em(ie))-edge+emu,(dble(xseclg(ii)),ii=0,ljmax), &
             (imag(xseclg(ii)),ii=0,ljmax),dble(lgsum),dimag(lgsum)
         lgsum=dcmplx(0.0d0,0.0d0)
         do ii=0,ljmax
            lgsum=lgsum+xsecl(ii)
         end do
         write(77,innerform) dble(em(ie))-edge+emu,(dble(xsecl(ii)),ii=0,ljmax), &
             (imag(xsecl(ii)),ii=0,ljmax),dble(lgsum),dimag(lgsum)
      end do                    ! ie
!     end of energy cycle

!
!     Now output the spectra in padlib form to unit 7
!
!     add all kinds of other useful information
!
!
      if (ldecmx.ge.0) then 
         write(7,'(3I5)') kfinmax,indmax,jinit
         do ii=1,indmax
            write(7,'(4I5)') kiind(ii),lgind(ii),ljind(ii),lind(ii)
         end do
         do ie=1,nex
            call wrpadx (7, npadx, atomxsec(1,ie),kfinmax )
         end do
         close(7)
      end if
      close(66)
	  
	  
	  
      return
      end
