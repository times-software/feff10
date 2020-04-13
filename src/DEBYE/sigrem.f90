!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: sigrem.f90,v $:
! $Revision: 1.12 $
! $Author: jorissen $
! $Date: 2011/11/29 05:48:20 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
!     program sigem
!
!     calculate the Debye-Waller factors for each MS path
!     using the equation-of-motion methods
!
!     input files:  feff.inp and spring.inp
!
!     version 2  ( January 99)
!
!     coded by  A. Poiarkova
!
!---------------------------------------------------------------------
!  References:  
!             for the EM method: Phys. Rev. B , 59, p.948, 1999
!     also see dissertation 
!        "X-ray Absorption Fine Structure Debye-Waller Factors"
!         by Anna V. Poiarkova
!
!---------------------------------------------------------------------
!         tk temperature in degrees K
!         nleg  nlegs in path
!         rat   positions of each atom in path
!         NB Units of distance in this routine
!            are angstroms, including sig2. 
!         sig2 is output DW factor
!
!---------------------------------------------------------------------
      subroutine sigem (sig2mx, sig2x, iem, tk, ipath, nleg, rat, sig2)
      use par
	  use dimsmod,only: nphx=>nphu,natx,legtot  !KJ 7-09 replaces "includ 'dwpar.h'"
      implicit double precision (a-h, o-z)

!KJ      include 'dwpar.h'
!      parameter (nphx = nphx1)
!      parameter (natx = natxdw)

! local parameters:
      parameter (amu0  = 1.66054)
      parameter (pi = 3.1415926535897932384626433d0)
      parameter (nwx = 700)

      double precision sig2mx, sig2x(0:nphx,0:nphx)
      dimension iphat(natx), izph(0:nphx)

! variables shared with rdspr.f:
      dimension rat1(3,natx), iz(natx)
	  !KJ 11-2011: made following allocatable to pull them out of (limited - 32b) static memory space:
	  real*8, allocatable :: dm (:,:,:,:),rnn(:,:,:)
	  integer,allocatable :: nnl(:,:)
      !dimension dm(3,3,natx,natx)
      !dimension rnn(3,natx,natx)
      !dimension nnl(natx,natx)

! local variables:
      dimension rat(3,0:legtot) !KJ 7-09 replaced nlegx1 with legtot
      dimension nconv(0:legtot) !KJ 7-09 replaced nlegx1 with legtot
      dimension q0(3,natx)
      dimension gr(nwx), w(nwx)
      dimension nq0(0:legtot) !KJ 7-09 replaced nlegx1 with legtot 
      dimension uu(3,natx), up(3,natx), ff(3,natx)

      character*30  fname
      parameter (ntitx1 = 10)
      character*71  title(ntitx1)
      dimension ltit(ntitx1)
!     character*80  titlep(ntitx1)
integer nat !KJ

      character*512 slog
      logical iem_open
      integer istrln
      external istrln

      save 
      data nsigc /0/
!-------------------------------------------------------------
      inquire(unit=iem,opened=iem_open)
      if (nsigc.eq.0) then


! Read coordinates and potentials from feff.inp
      call dwrdin (rat1, iphat, izph, natom, ntitle, title, ltit)
      if (natom.gt.natx) natom=natx
      do 5 iat=1, natom
         iz(iat) = izph(iphat(iat))
         if (iphat(iat).eq.0) i0=iat
5     continue

            allocate( dm(3,3,natom,natom),rnn(3,natom,natom),nnl(natom,natom)) !KJ 02-14
!      allocate( dm(3,3,natx,natx),rnn(3,natx,natx),nnl(natx,natx)) !KJ 11-11

      nat=natom
! Read spring.inp and build dynamical matrix
      call rdspr(rat1, iz, natom, i0,                                   &
     &           dm, rnn,                                               &
     &           acut, res, wmax, dosfit, zshell, w0,                   &
     &           rintr, iprdos, nnl,nat )

            write(slog,7)
   7        format(2x,'Calculating Debye-Waller factors via EMM...')
            call wlog(slog)
            write(slog,9)
   9        format(2x,'This might take a while.')
            call wlog(slog)

      if (ipath.ne.0.and.iem_open) then
!           Echo title cards to s2_em.dat
            do 10  i = 1, ntitle
               write(iem,12)  title(i)(1:ltit(i))
  10        continue
  12        format (1x, a)
            write(iem,17) tk, natom
  17        format(1x,'temperature =',f7.2,2x,'N_at =',i4)
            write(iem,19)
  19        format (1x, 71('-'))
            write(iem,25)
            write(slog,25)
            call wlog(slog)
  25        format(3x,'ipath',4x,'nleg',3x,'sig2',5x,                   &
     &            'mu_ipath',2x,'check0(%)')
      endif

! Integration parameters:
      wmaxx=sqrt(zshell)
      dt=2.*pi/wmaxx/15.
! top limit in t integration:
      cutoff=2.*sqrt(2.*acut)/res/wmaxx 
      nstep=cutoff/dt
      xlam=acut/(cutoff)**2
      wl=0.0000001
! top limit in w integration:
      wm=wmax*wmaxx 
      dw=0.01 
      nw=(wm-wl)/dw + 1
      if (nw .gt. nwx) then
          nw = nwx
          dw = (wm-wl)/(nw -1)
      endif
      nfit = dosfit*nw/20.

      endif
!------------------------------------

!c    Open path input file (unit in) and read title.  Use unit 2.
!     ntitle2 = 5
!     open(unit=2,file='paths.dat',status='old', iostat=ios)
!     call chopen (ios, 'paths.dat', 'sigem')
!     call rdhead (2, ntitle2, titlep, ltit)
!c    if (ntitle2 .le. 0)  then
!c       titlep(1) = ' '
!c    endif

! 84  continue
!     read(2,*,end=1010) ipath, nleg
!     skip label (x y z ipot rleg beta eta) and read the path
!     read(2,*)
      do 78 ileg=0,nleg
!        read(2,*,end=1010) (rat(j,ileg),j=1,3)
         nconv(ileg)=0
  78  continue

      do 88 n=1,3
         aa = rat(n,nleg)
         do 87 i=0,(nleg-2)
            j=nleg-i
            rat(n,j)=rat(n,j-1)
  87     continue
         rat(n,1)=aa
  88  continue
      do 89 i=1,nleg
         nq0(i)=0.
  89  continue

! nconv converts # of an atom in the nleg list of coordinates (rat) to
! its # in the full list of all atomic coordinates (rat1)
      do 94 i=1,natom
         do 91 n=1,3
   91    q0(n,i)=0.
         do 95 jl=1,nleg
            m=0
            do 93 n=1,3
               l=nint(100.*rat(n,jl))
               l1=nint(100.*rat1(n,i))
               if (abs(l-l1).le.1) m=m+1
   93       continue
            if (m.eq.3) then
              nconv(jl)=i
            endif
   95    continue
   94 continue
!     check that all path atoms are found
      do 96 jl=1,nleg
        if (nconv(jl).eq.0) then
           print*,' did not find atom jl=', jl
           print*, rat(1,jl),rat(2,jl),rat(3,jl)
           call par_stop('SIGREM-1')
        endif
  96  continue

      atmu=0.
      iq0=0
      nconv(0)=nconv(nleg)
      do 100 il=1,nleg
         l=nconv(il)
         do 101 jq=1,iq0
  101    if(nq0(jq).eq.l) go to 102
         iq0=iq0+1
         nq0(iq0)=l
  102    continue
         nq0x=iq0
         i=nconv(il)
         im=nconv(il-1)
         ip=nconv(il+1)
!        if (il.eq.1) im=nconv(nleg)
         if (il.eq.nleg) ip=nconv(1)
         atmass=atwtd(iz(i))
      do 100 n=1,3
         atmu=atmu + 0.25*( rnn(n,i,im)+rnn(n,i,ip) )**2 /atmass
  100 continue
      atmu=1./atmu
      icount=1
  108 continue
      icount= icount+1
      if (icount.gt.10) call par_stop('SIGREM-2')

      do 115 i=1,natom
      do 115 n=1,3
  115 q0(n,i)=0.

! Build initial state vector |Q_j(0)> for the current path
      do 116 n=1,3
      do 116 il=1,nleg
         i=nconv(il)
         im=nconv(il-1)
         ip=nconv(il+1)
         if (il.eq.1) im=nconv(nleg)
         if (il.eq.nleg) ip=nconv(1)
         atmass=atwtd(iz(i))
         q0(n,i)=q0(n,i)+sqrt(atmu/atmass)*(rnn(n,im,i)-rnn(n,i,ip))/2.
  116 continue

! make sure it's normalized <Q_j(0)|Q_j(0)>=1
      q0q0=0.
      do 120 iq0=1,nq0x
      i=nq0(iq0)
      do 120 n=1,3
         q0q0=q0q0+q0(n,i)*q0(n,i)
 120  continue
      p00=nint(q0q0*1000.d0)/1000.d0
      if (abs(p00-1.d0).gt.5.d-4) then
         atmu=atmu/q0q0
         go to 108
      endif

!     to get THz units:
      wnorm=100.*w0/sqrt(amu0*10.) 
!*** moments
      a0=0.
      do 132 il=1,nq0x
      do 132 im=1,nq0x
         l=nq0(il)
         m=nq0(im)
      do 132 n1=1,3
      do 132 n2=1,3
         a0 = a0 + q0(n1,l)*dm(n1,n2,l,m)*q0(n2,m)/w0/w0
  132 continue
      a0=wnorm*sqrt(a0)

      do 125 kw=1, nwx
         gr(kw) = 0.
  125 w(kw) = (wl+(kw-1)*dw)

!  make file prdennnnn.dat
      if (master.and.ipath.ne.0.and.ipath.le.iprdos) then 
         write(fname,130)  ipath
  130    format('prden', i4.4, '.dat')
         open (unit=25, file=fname, status='unknown',iostat=ios)
         call chopen (ios, fname, 'sigem')
         do 134  i = 1, ntitle
            write(25,136)  title(i)(1:ltit(i))
  134    continue
         write(25,135) natom
  135    format('#',1x,'N_at =', i4)
  136    format ('#',1x, a)
         write(25,138)
  138    format ('#',1x, 71('-'))
      endif


!  set initial conditions
      do 150 i=1, natom
      do 150 n=1,3
         uu(n,i)=q0(n,i)
         up(n,i)=uu(n,i)
  150 continue

! Solve 3*natom equations of motion and find projected VDOS (gr)
      dt2=dt*dt
      t=dt/2.
      do 200 kstep = 1, nstep
!        damping factor:
         e1=exp(-xlam*t*t) 
         xat=0.
         do 167 i=1, natom
         do 167 n=1,3
  167    xat = xat + uu(n,i)*q0(n,i)
         xat=xat*e1
         do 170 kw=1, nw
  170    gr(kw) = gr(kw) + xat*cos(w(kw)*t)*dt
         if(kstep.eq.nstep) go to 200

         do 175 i=1,natom
         do 175 n=1,3
  175    ff(n,i)=0.

         do 180 i=1,natom
            jn=1
  185       if (nnl(i,jn).ne.0) then
               j=nnl(i,jn)
               am=w0*w0
               do 187 n1=1,3
               do 187 n2=1,3
                  ff(n1,i)=ff(n1,i)-dm(n1,n2,i,j)*uu(n2,j)/am
                  if(i.ne.j) ff(n1,j)=ff(n1,j)-dm(n1,n2,j,i)*uu(n2,i)/am
  187          continue
               jn = jn + 1
               go to 185
            endif
  180    continue

         do 199 i=1,natom
         do 199 n=1,3
            put=2.*uu(n,i)-up(n,i)+dt2*ff(n,i)
            up(n,i)=uu(n,i)
            uu(n,i)=put
  199    continue

  200 t=t+dt

      afit = 0.
      if (nfit.ne.0) then
         if (w(nfit).ne.0.) afit=gr(nfit)/(w(nfit)**4)
      endif

! fit vibr.density to A*w^4, for low w
      do 225 kw=1, nfit
         gr(kw)=afit*w(kw)**4
  225 continue

! Normalization of the pr.density of modes 
! (it's the 2/pi factor which was left out till now with,
! perhaps, a small diffrence due to the fit) 
      gr(nw)=0.
      if (gr(1).lt.0.) gr(1)=0.
      xx=(gr(1)+gr(nw))*dw/2.
      do 247 kw=2, (nw-1)
         if (gr(kw).lt.0.) gr(kw)=0.
  247 xx = xx + gr(kw)*dw
      cn1=1./xx

      if (master.and.ipath.ne.0.and.ipath.le.iprdos) then
! to get THz units:
         wnorm=100.*w0/sqrt(amu0*10.) 
         write(25,349) ipath, nleg
  349    format('#',2x,'ipath =',i3,2x,'nleg =',i2)
         write(25,350)
! 350    format(1h#,6x,5hw,THz,18x,6hrho(w))
  350    format(1h#,6x,'cm^-1',18x,6hrho(w))
         do 370 kw=1,nw
            write(25,360) w(kw)*wnorm*100./6./pi, gr(kw)*cn1/wnorm
!           write(25,360) w(kw)*wnorm, gr(kw)*cn1/wnorm
  360       format(2x,f10.3,15x,f10.7)
  370    continue
         close (unit=25)
      endif

      wt=tk/187.64/w0
      ccc=cn1
      check0 = abs((2./pi - cn1)/(2./pi))
      check0=check0*100.
      coef = ccc*0.5*0.2587926/atmu/w0
! integrate over w to get sig2
      cth=0.
      s2=0.
!     gr(1)=0.
      do 400 kw=2, (nw-1)
         cth = 1./tanh( w(kw)/(2.*wt) )
         s2 = s2 + coef*gr(kw)*cth*dw/w(kw)
  400 continue
      sig2 = s2

      if (ipath.ne.0.and.iem_open) then
         write(iem,473) ipath, nleg, sig2,atmu,check0
         write(slog,473) ipath, nleg, sig2,atmu,check0
         call wlog(slog)
  473    format(4x,i3,4x,i3,4x,f7.5,3x,f7.3,4x,f5.2)
      endif

      nsigc = nsigc + 1
!1000 go to 84
!1010 continue
!     close (unit=2)
      if (sig2.gt.1.0) then
        sig2 = 1.0d0
        call wlog (' WARNING: Found sig**2>1. Set sig2=1. ')
        write (slog,1011) nconv(1), nconv(2)
 1011   format('          Possible zero frequency modes with atoms', i4,&
     &   ' or', i4)
         call wlog(slog)
         call wlog('          Check springs.inp')
      endif
      if (check0.gt.5.0) then
        write (slog,*) ' WARNING: Failed check0 test:missing VDOS.',    &
     &  ' Reduce dosfit and/or increase wmax.'
        call wlog(slog)
      endif

!     update maximum DW factors
      if (sig2.gt.sig2mx) sig2mx=sig2
      if (sig2.gt.sig2x( iphat(nconv(1)),  iphat(nconv(2)) )) then
         sig2x( iphat(nconv(1)),  iphat(nconv(2)) ) = sig2
         sig2x( iphat(nconv(2)),  iphat(nconv(1)) ) = sig2
      endif

      return
      end
!----------------------------------------------------
      subroutine dwrdin (rat, iphat, izph, nat,                         &
     &            ntitle, title, ltit)

!     Read feff.inp for sigem.f
!     (here we need only coordinates and potentials)
	  use dimsmod,only: nphx=>nphu,natx,legtot,nwordx  !KJ 7-09 replaces "includ 'dwpar.h'"
         implicit double precision (a-h, o-z)

!      include 'dwpar.h'
!      parameter (nphx = nphx1)
!      parameter (natx = natxdw)

      dimension iphat(natx)
      dimension rat(3,natx)
      dimension iatph(0:nphx)

      character*6  potlbl(0:nphx)

!     Local stuff
      character*150  line
!      parameter (nwordx = 20)
      character*20 words(nwordx)

      parameter (ntitx = 10)
      character*71  title(ntitx)
      dimension ltit(ntitx)
      dimension izph(0:nphx)
      logical iscomm
      parameter (nssx = 16)

      parameter (big = 1.0e5)
      character*512 slog
      external itoken

   10 format (a)
   20 format (bn, i15)
   30 format (bn, f15.0)

!c    initialize things

      ntitle = 0

      nat = 0
      do 100  iat = 1, natx
         iphat(iat) = -1
  100 continue

      nph = 0
      do 110  iph = 0, nphx
         iatph(iph) = 0
         izph(iph) = 0
         potlbl(iph) = ' '
  110 continue

!     Open feff.inp, the input file we're going to read
      open (unit=1, file='feff.inp', status='old', iostat=ios)
      call chopen (ios, 'feff.inp', 'rdinp')

!     tokens  0 if not a token
!             1 if ATOM (ATOMS)
!             7 if TITL (TITLE)
!            10 if DEBY (DEBYE)
!            13 if PRIN (PRINT)
!            14 if POTE (POTENTIALS)
!            -1 if END  (end)
!     mode flag  0 ready to read a keyword card
!                1 reading atom positions
!                2 reading overlap instructions for unique pot
!                3 reading unique potential definitions

      mode = 0
  200 read(1,10,iostat=ios)  line
         if (ios .lt. 0)  line='END'
         call triml (line)
         if (iscomm(line))  goto 200
         nwords = nwordx
         call bwords (line, nwords, words)
         itok = itoken (words(1),'feff.inp')

!        process the card using current mode
  210    continue

         if (mode .eq. 0)  then
            if (itok .eq. 1)  then
!              ATOM
!              Following lines are atom postions, one per line
               mode = 1

            elseif (itok .eq. 7)  then
!              TITLE title...
               ntitle = ntitle + 1
               if (ntitle .le. ntitx)  then
                  title(ntitle) = line(6:)
                  call triml (title(ntitle))
               else
                  call wlog(' Too many title lines, title ignored')
                  call wlog(' ' // line(1:71))
               endif
               mode = 0

!           elseif (itok .eq. 10)  then
!c             DEBYE  temp debye-temp
!c                  temps in kelvin
!c                  These add to any sig2 from SIG2 card or files.dat
!              read(words(2),30,err=900)  tk
!              read(words(3),30,err=900)  thetad
!              idwopt=0
!              read(words(4),20,err=900)  idwopt
!              mode = 0
            elseif (itok .eq. 14)  then
!              POTENTIALS
!              Following lines are unique potential defs, 1 per line
               mode = 3

            elseif (itok .eq. -1)  then
!c             END
               goto 220
            else
               mode = 0
! *            write(slog,'(1x,a)') line(1:70)
! *            call wlog(slog)
! *            write(slog,'(1x,a)') words(1)
! *            call wlog(slog)
! *            write(slog,'(a,i8)') ' Token ', itok
! *            call wlog(slog)
! *            call wlog(' Keyword unrecognized.')
! *            call wlog(' See FEFF document -- some old features')
! *            call wlog(' are no longer available.')
! *            call par_stop('DWRDIN-1')
            endif
         elseif (mode .eq. 1)  then
            if (itok .ne. 0)  then
!c             We're done reading atoms.
!c             Change mode and process current card.
               mode = 0
               goto 210
            endif
            nat = nat+1
            if (nat .gt. natx)  then
               write(slog,'(1x,a,i5)') 'Too many atoms, max is ', natx
               call wlog(slog)
               write(slog,'(1x,a,i5,a)') 'Only', natx,                  &
     &       ' atoms will be considered in the DW factor calculations.'
               call wlog(slog)
               nat = nat-1
               mode = 0
               goto 210
!              call par_stop('DWRDIN-2')
            endif
            if (nat.le.natx) then
               read(words(1),30,err=900)  rat(1,nat)
               read(words(2),30,err=900)  rat(2,nat)
               read(words(3),30,err=900)  rat(3,nat)
               read(words(4),20,err=900)  iphat(nat)
               if (iphat(nat).eq.0) iat0 = nat
            else
               mode = 0
               goto 210
            endif
         elseif (mode .eq. 3)  then
            if (itok .ne. 0)  then
!c             We're done reading unique potential definitions
!c             Change mode and process current card.
               mode = 0
               goto 210
            endif
            read(words(1),20,err=900)  iph
            if (iph .lt. 0  .or.  iph .gt. nphx)  then
               write(slog,'(a,i8)')                                     &
     &             'Unique potentials must be between 0 and ',          &
     &             nphx
               call wlog(slog)
               write(slog,'(i8,a)') iph, ' not allowed'
               call wlog(slog)
               write(slog,'(1x,a)') line(1:71)
               call wlog(slog)
               call par_stop('DWRDIN-3')
            endif
            read(words(2),20,err=900)  izph(iph)
!c          No potential label if user didn't give us one
!c          Default set above is potlbl=' '
            if (nwords .ge. 3)  potlbl(iph) = words(3)
         else
            write(slog,'(a,i8)')                                        &
     &        'DWRDIN-4: Mode unrecognized, mode ', mode
!           call wlog(slog)
            call par_stop(slog)
         endif
      goto 200
  220 continue

!c    We're done reading the input file, close it.
      close (unit=1)
            if (nat .gt. natx)  then
               write(slog,'(a,i8)') 'Too many atoms for DW calculations,&
     &         max is ', natx
               call wlog(slog)
               write(slog,'(a,i8,a)') 'Only atoms up to #',natx,        &
     &         '  will be considered'
               call wlog(slog)
            endif

      do 250 iat = 1, nat
      do 250 i = 1,3
        if (iat.ne. iat0) rat(i,iat) = rat(i,iat) - rat(i,iat0)
 250  continue
      do 251 i = 1,3
 251  rat(i,iat0) = 0.d0

!c    Find out how many unique potentials we have
      nph = 0
      do 300  iph = nphx, 0, -1
         if (izph(iph) .gt. 0)  then
            nph = iph
            goto 301
         endif
  300 continue
  301 continue
!c    Must have central atom
      if (izph(0) .le. 0)  then
         call wlog(' No absorbing atom (unique pot 0) was defined.')
         call par_stop('DWRDIN-5')
      endif
!c    Find central atom (only 1 permitted)
      iatabs = -1
      do 400  iat = 1, nat
         if (iphat(iat) .eq. 0)  then
            if (iatabs .lt. 0)  then
               iatabs = iat
            else
               call wlog(' More than one absorbing atom (potential 0)')
               call wlog(' Only one absorbing atom allowed')
               call par_stop('DWRDIN-6')
            endif
         endif
  400 continue

!c    Then find model atoms for unique pots that have them
!c    Use atom closest to absorber for model
      do 330  iph = 0, nphx
         rabs = big
         do 320  iat = 1, nat
            if (iph .eq. iphat(iat))  then
               tmp = dist (rat(1,iat), rat(1,iatabs))
               if (tmp .lt. rabs)  then
!c                this is the closest so far
                  rabs = tmp
                  iatph(iph) = iat
               endif
            endif
  320    continue
  330 continue
!c    if iatph > 0, a model atom has been found.

      if (ntitle .le. 0)  then
         ntitle = 1
         title(1) = 'Null title'
      endif
      do 490  i = 1, ntitle
         ltit(i) = istrln (title(i))
  490 continue

      return

  900 continue
      call wlog(' Error reading input, bad line follows:')
      write(slog,'(1x,a)') line(1:71)
      call wlog(slog)
      call par_stop('DWRDIN-7 fatal error.')

!      return
      end

!----------------------------------------------------------
      subroutine rdspr(rat1, iz, natom, i0,                             &
     &           dm, rnn,                                               &
     &           acut, res, wmax, dosfit, zshell, w0,                   &
     &           rintr, iprdos, nnl, nat)

!     Read spring.inp for multiple scattering feff and
!     build dynamical matrix.
      use par
	  use dimsmod,only: nphx=>nphu,natx,legtot,nwordx  !KJ 7-09 replaces "includ 'dwpar.h'"
         implicit double precision (a-h, o-z)
!      include 'dwpar.h'
!!     parameter (nphx = nphx1)
!      parameter (natx = natxdw)
       integer,intent(in)::nat 

! new local parameters:
      parameter (nangx = 7*natx)
      parameter (nsprx = 40)
      parameter (nshx = 100)

! variables shared with sigem.f:
      dimension rat1(3,natx), iz(natx)
      dimension dm(3,3,nat,nat) !KJ natx=>nat
      dimension rnn(3,nat,nat) !KJ
      dimension nnl(nat,nat) !KJ

! local variables:
    real*8, allocatable :: rshell(:,:), ang(:), dang(:)
    integer, allocatable :: nspr(:,:), nang(:,:)
      dimension si(3), sj(3), sk(3)
	  !KJ :
      real*8,allocatable :: dmstr(:,:,:,:),dma(:,:,:,:), drij(:,:), str(:,:)
      !dimension dmstr(3,3,natx,natx),dma(3,3,natx,natx), drij(natx,natx)
      !dimension str(natx,natx)

      character*150  line
!      parameter (nwordx = 20)
      character*20  words(nwordx)
      character*512 slog

      logical iscomm
      integer nnang

   10 format (a)
   20 format (bn, i15)
   30 format (bn, f15.0)

      nnang=7*nat
      
      allocate (dmstr(3,3,nat,nat),dma(3,3,nat,nat), drij(nat,nat), str(nat,nat))
!      allocate(rshell(natx,0:nshx),nspr(2,nsprx),ang(nangx), dang(nangx),nang(3,nangx))
      allocate(rshell(natx,0:nshx),nspr(2,nsprx),ang(nnang), dang(nnang),nang(3,nnang))

! initialize things


      str(:,:)=0.
      drij(:,:)=0.02
      ang(:)=0.
      dang(:)=0.
      nang(:,:)=0
      nspr(:,:)=0

      acut=3.
      res=0.05
      dosfit=0.
      wmax=1.
      na=1
      nintr=0
      strx=10000.
      ispr=1
      iprdos = 0
      ddrij=0.02
      ddang=0.02

      open(unit=1,file='spring.inp',status='old', iostat=ios)
      call chopen (ios, 'spring.inp', 'rdspr')

!     tokens  0 if not a token
!             1 if STRE (STRETCHES)
!             2 if ANGL (ANGLES)
!             3 if VDOS
!             4 if PRDOS 
!            -1 if END  (end)
!     mode flag  0 ready to read a keyword card
!                1 reading stretches
!                2 reading angle-bends
 
      mode = 0
  200 read(1,10,iostat=ios)  line
         if (ios .lt. 0)  line='END'
         call triml (line)
         if (iscomm(line))  goto 200
         nwords = nwordx
         call bwords (line, nwords, words)
         itok = itoken (words(1),'spring.inp')

!        process the card using current mode
  210    continue

         if (mode .eq. 0)  then
            if (itok .eq. 1)  then
!              STRE
!              Following lines are stretches, one per line
!              read(words(2),20,err=900)  nintr
               mode = 1
            elseif (itok .eq. 2)  then
!              ANGL
!              Following are ...
               mode = 2
            elseif (itok .eq. 3)  then
!              VDOS
!              VDOS  resolution, a_cut, wmax, dosfit
!               0 - do not run modules, 1 - run module
               read(words(2),30,err=900)  res
               read(words(3),30,err=900)  wmax
               read(words(4),30,err=900)  dosfit
               if (nwords.gt.4) then
                   read(words(5),30,err=900)  acut
               endif
               mode = 0
            elseif (itok .eq. 4)  then
!              PRINT  iprdos
!              to print or not to print prdennnnn.dat files;
!              if the card is present, these files will be
!              printed for paths 1 through iprdos
               iprdos = 1 
               read(words(2),20,err=900)  iprdos
               mode = 0
            elseif (itok .eq. -1)  then
!              END
               goto 220
            else
               write(slog,'(1x,a)') line(1:71)
               call wlog(slog)
               write(slog,'(1x,a)') words(1)
               call wlog(slog)
               write(slog,'(a,i8)') ' Token ', itok
             call wlog(slog)
               call wlog(' Keyword unrecognized.')
               call wlog(' See FEFF document -- some old features')
               call wlog(' are no longer available.')
             call par_stop('RDSPR-1')
            endif
         elseif (mode .eq. 1)  then
            if (itok .ne. 0)  then
!              We're done reading stretches
!              Change mode and process current card.
               mode = 0
               goto 210
            endif
            read(words(1),20,err=900) ii
            i=ii+1
            call chekin (ii, natom, line)
            read(words(2),20,err=900) jj
            j=jj+1
            call chekin (jj, natom, line)
            read(words(3),30,err=900) str(i,j)
            if (str(i,j).lt.strx) then
               strx=str(i,j)
               ix=i
               jx=j
            endif
            nspr(1,ispr)=i
            nspr(2,ispr)=j
            ispr=ispr+1
            read(words(4),30,err=900) ddrij
            drij(i,j) = abs(ddrij)/100.
            drij(j,i) = abs(ddrij)/100.
         elseif (mode .eq. 2)  then
            if (itok .ne. 0)  then
!              We're done reading angle-bends
!              Change mode and process current card.
               mode = 0
               goto 210
            endif
            read(words(1),20,err=900) ii
            i=ii+1
            call chekin (ii, natom, line)
            read(words(2),20,err=900) jj
            j=jj+1
            call chekin (jj, natom, line)
            read(words(3),20,err=900) kk
            k=kk+1
            call chekin (kk, natom, line)
            read(words(4),30,err=900) ang(na)
            nang(1,na)=i
            nang(2,na)=j
            nang(3,na)=k
            read(words(5),30,err=900) ddang
            dang(na) = abs(ddang)/100.
            na=na+1
         else
            write(slog,'(a,i8)') 'Mode unrecognized, mode ', mode
            call wlog(slog)
            call par_stop('RDSPR-2')
         endif
      goto 200
  220 continue

!     We're done reading the input file, close it.
      close (unit=1)
      nax=na-1
      call wlog('done reading spring.inp')

!     write statistics on found bonds and angles into spring.dat 
      if (master) then
        open (unit=2, file='spring.dat', status='unknown',iostat=ios)
        call chopen (ios, 'spring.dat', 'spring')
        write(2,*) ' Statistics on spring constants in spring.inp.'
        write(2,*) '   STRETCHES  i  j  aa   found_number'  
      endif

! find all stretching bonds
      do 321 jspr=1, (ispr-1)
         icnt=0
         i=nspr(1,jspr)
         j=nspr(2,jspr)
         aa = str(i,j)
         ddrij=drij(i,j)
         ip=iz(i)
         jp=iz(j)
         rij = dist (rat1(1,i), rat1(1,j))
         if (aa.eq.0.) go to 321
         do 320 k=1, natom
            do 320 l=k+1, natom
               kp=iz(k)
               lp=iz(l)
               rkl = dist (rat1(1,k), rat1(1,l))
               comp = abs(rij/rkl - 1.)
               if (comp.gt.ddrij) go to 320
               if (ip.ne.kp.or.jp.ne.lp) then
                  if (ip.ne.lp.or.jp.ne.kp) go to 320
               endif
               str(k,l) = aa
               str(l,k) = aa
!alex       to check the bonds, that were found
!alex        print*, k,l,aa
                icnt = icnt+1
  320    continue
         str(j,i) = aa
         if (master) write (2,*) i-1, j-1, aa, icnt
  321 continue
      if (master) write(2,*) '   BENDS   i  j  k   aa   found_number'  

! find all bending angles
      naxx=nax
      do 323 na=1,nax
         icnt=1
         i=nang(1,na)
         j=nang(2,na)
         k=nang(3,na)
         ddrij=drij(i,j)
         ddrkj=drij(k,j)
         ip=iz(i)
         jp=iz(j)
         kp=iz(k)
         call coss(rat1(1,i),rat1(1,j),rat1(1,k),cosijk)
         rij = dist (rat1(1,i), rat1(1,j))
         rkj = dist (rat1(1,k), rat1(1,j))
         aa=ang(na)
!        print*, na, i,j,k, aa
         do 326 ii=1, natom
         do 326 jj=1, natom
            if (ii.eq.jj) go to 326
            rrij=dist (rat1(1,ii), rat1(1,jj))
            do 322 kk=ii+1, natom
               if (kk.eq.jj) go to 322
               rrkj=dist (rat1(1,kk), rat1(1,jj))
               comp1 = abs(rrij/rij - 1.)
               comp2 = abs(rrkj/rkj - 1.)
               if (comp1.gt.ddrij.or.comp2.gt.ddrkj) then
                  comp1 = abs(rrkj/rij - 1.)
                  comp2 = abs(rrij/rkj - 1.)
                  if (comp1.gt.ddrij.or.comp2.gt.ddrkj) go to 322
               endif
               iip=iz(ii)
               jjp=iz(jj)
               kkp=iz(kk)
            if (iip.ne.ip.or.jjp.ne.jp.or.kkp.ne.kp) then
               if (kkp.ne.ip.or.jjp.ne.jp.or.iip.ne.kp) go to 322
            endif
               call coss(rat1(1,ii),rat1(1,jj),rat1(1,kk),cssijk)
               if (dacos(cosijk).eq.0.) go to 322
                  comp = abs( dacos(cssijk)/dacos(cosijk) -1.)
               if (comp.ge.dang(na)) go to 322
               do 324 na1=1,naxx
                  ii1=nang(1,na1)
                  jj1=nang(2,na1)
                  kk1=nang(3,na1)
               if (ii.eq.ii1.and.jj.eq.jj1.and.kk.eq.kk1) go to 322
               if (kk.eq.ii1.and.jj.eq.jj1.and.ii.eq.kk1) go to 322
 324           continue
               naxx=naxx+1
               ang(naxx)=aa
               nang(1,naxx)=ii
               nang(2,naxx)=jj
               nang(3,naxx)=kk
!alex          to check the bends, that were found
!              print*, naxx, ii,jj,kk,aa
               icnt = icnt + 1
               if (naxx.eq.nangx) goto 333
 322        continue
 326     continue
         if (master) write (2,*) i-1, j-1, k-1, aa, icnt
 323  continue
 333  continue

      if (master) close (unit=2)


         rshell(:,:)=0.


! find shells
      rintr=0.
      nintr=1
      do 330 i=1, natom
      nshell=0
      do 335 j=1, natom
         if (j.eq.i) go to 332
         if (nshell.gt.nshx) go to 332
         rij = dist (rat1(1,i), rat1(1,j))
         ddrij=drij(i,j)
         ncount=0
         do 331 ish=0, nshell
            b = real(rshell(i,ish))
            dif=1.
            if (b.ne.0.) dif = abs(rij -b)/b
            if (dif.le.ddrij) ncount=ncount+1 
 331     continue
         if (ncount.eq.0) then
            nshell = nshell + 1
            if (str(i,j).ne.0.and.rij.gt.rintr) rintr=rij
            rshell(i,nshell) = rij
         endif
 332     do 335 n=1,3
         rnn(n,i,j)=0.
         do 335 m=1,3
            dmstr(n,m,i,j)=0.
            dma(n,m,i,j)=0.
            dm(n,m,i,j)=0.
 335  continue
! sort rshell into ascending numerical order
! and find maximum order of interacting neighbor nintr
      do 342 jsh=2,nshell 
         aa = rshell(i,jsh)
         do 341 ish=jsh-1,1,-1
            if(rshell(i,ish).le.aa) go to 340
            rshell(i,ish+1)=rshell(i,ish)
 341     continue
         ish=0
 340     rshell(i,ish+1) = aa
         if (aa.le.rintr.and.(ish+1).gt.nintr) nintr = ish+1
 342  continue
 330  continue
call wlog('rshells sorted')
      zshell=0.
      i1=0
      do 350 i=1,natom
         nnl(i,:)=0
         do 350 j=i+1,natom
            dx=rat1(1,j)-rat1(1,i)
            dy=rat1(2,j)-rat1(2,i)
            dz=rat1(3,j)-rat1(3,i)
            dr=sqrt(dx*dx+dy*dy+dz*dz)
            rnn(1,i,j)=dx/dr
            rnn(2,i,j)=dy/dr
            rnn(3,i,j)=dz/dr
            rnn(1,j,i)=-rnn(1,i,j)
            rnn(2,j,i)=-rnn(2,i,j)
            rnn(3,j,i)=-rnn(3,i,j)
            rrij = abs( dr/rshell(1,1) -1.)
!           if (i.eq.1.and.rrij.le.drij(i,j)) zshell=zshell+1
            if (i.eq.i0.and.rrij.le.drij(i,j)) then
               zshell=zshell+1
               if (i1.eq.0.and.str(i,j).ne.0.) i1 = j
            endif
  350 continue

! Build dynm. matrix for angle bends 
      nan=0
      do 355 na=1,naxx
!        print*,na, naxx, nangx
         i=nang(1,na)
         j=nang(2,na)
         k=nang(3,na)
         if (i.eq.j.or.j.eq.k) go to 355
         if(ang(na).eq.0.) go to 355
         nan=nan+1
         rij = dist(rat1(1,i),rat1(1,j))
         rkj = dist(rat1(1,k),rat1(1,j))
         if (rij.gt.rintr.or.rkj.gt.rintr) go to 355
         call sang (i, j, k, rat1, si, sj, sk)
         do 357 n1=1,3
         do 357 n2=1,3
            dma(n1,n2,i,j)=dma(n1,n2,i,j)+ang(na)*si(n1)*sj(n2)
            dma(n1,n2,j,k)=dma(n1,n2,j,k)+ang(na)*sj(n1)*sk(n2)
            dma(n1,n2,i,k)=dma(n1,n2,i,k)+ang(na)*si(n1)*sk(n2)

            dma(n1,n2,i,i)=dma(n1,n2,i,i)+ang(na)*si(n1)*si(n2)
            dma(n1,n2,j,j)=dma(n1,n2,j,j)+ang(na)*sj(n1)*sj(n2)
            dma(n1,n2,k,k)=dma(n1,n2,k,k)+ang(na)*sk(n1)*sk(n2)

            dma(n2,n1,j,i)=dma(n1,n2,i,j)
            dma(n2,n1,k,j)=dma(n1,n2,j,k)
            dma(n2,n1,k,i)=dma(n1,n2,i,k)
  357    continue
  355 continue

! Build dynm. matrix for stretches 
      do 375 l=1,natom
      do 375 m=l,natom
         do 373 n1=1,3
            x2=str(l,m)*rnn(n1,l,m)
         do 373 n2=1,3
         dmi=0.
         if (l.eq.m) then
            do 377 i=1,natom
               if (real(str(i,m)).eq.0.) go to 377
               dmi = dmi + str(i,m)*rnn(n1,i,m)*rnn(n2,i,m)
 377        continue
         endif
         dmstr(n1,n2,l,m) = dmi - x2*rnn(n2,l,m)
         dmstr(n2,n1,m,l)=dmstr(n1,n2,l,m)
 373  continue
 375  continue

! Add two dynm. matrices D_str+D_ang
      lnx=0
      do 380 i=1,natom
      ami=sqrt(atwtd(iz(i)))
      in=0
      do 380 j=1,natom
      amj=sqrt(atwtd(iz(j)))
      sumdm=0.
      do 381 n1=1,3
      do 381 n2=1,3
         dmdm=dma(n1,n2,i,j)+dmstr(n1,n2,i,j)
         sumdm=sumdm+abs(dmdm)
         dm(n1,n2,i,j)=dmdm/ami/amj
 381  continue
      if (real(sumdm).ne.0.and.i.le.j) then
         in=in+1
         nnl(i,in)=j
      endif
      if (in.ge.lnx) lnx=in
 380  continue

      atmu = 1./(1./atwtd(iz(i0)) + 1./atwtd(iz(i1)))
      a0=0.
      do 450 i=1,2
      do 450 j=1,2
      if (i.eq.1) l=i0
      if (i.eq.2) l=i1
      if (j.eq.1) m=i0
      if (j.eq.2) m=i1
      do 450 n1=1,3
      do 450 n2=1,3
         fact = (-1)**(i+j)*atmu
         atmass = 1./atwtd(iz(l))/atwtd(iz(m))
         a0 = a0 + fact*sqrt(atmass)*rnn(n1,i0,i1)*                     &
     &             dm(n1,n2,l,m)*rnn(n2,i0,i1)
  450 continue
!     effective freq. for the 1st shell:
      w0=sqrt(a0) 
      if (w0.eq.0.) then
         atmux = 1./(1./atwtd(iz(ix)) + 1./atwtd(iz(jx)))
         w0=sqrt(strx/atmux)
      endif

      deallocate(rshell,nspr,ang,dang,nang)

      
      return

  900 continue
      call wlog(' Error reading input, bad line follows:')
      write(slog,'(1x,a)') line(1:71)
      call wlog(slog)
      call par_stop('RDSPR fatal error.')

!c      return
      end

!---------------------------------------------------

      subroutine sang (i, j, k, rat1, si, sj, sk)

!*calculates coeficients s  (Sint=sum {si*ui}) connecting internal
!*coordinate delta_phi(ijk)=(valence ijk angle bend) with ui (atomic
!*displacements)

	  use dimsmod,only: nphx=>nphu,natx,legtot  !KJ 7-09 replaces "includ 'dwpar.h'"
      implicit double precision (a-h, o-z)
!      include 'dwpar.h'
!      parameter (natx = natxdw)

      dimension rat1(3,natx), rji(3), rjk(3)
      dimension eji(3), ejk(3), ej(3)
      dimension si(3), sk(3), sj(3)

      dji=0.
      djk=0.
      dik=0.
      do 905 m = 1, 3
         rji(m) = rat1(m,i) - rat1(m,j)
         rjk(m) = rat1(m,k) - rat1(m,j)
         dji = dji + rji(m)**2
         djk = djk + rjk(m)**2
         dik = dik + ( rat1(m,k) - rat1(m,i) )**2
         si(m) = 0.
         sj(m) = 0.
         sk(m) = 0.
  905 continue
      dji = sqrt(dji)
      djk = sqrt(djk)
      dik = sqrt(dik)

      dotj=0.
      do 910 m = 1, 3
         eji(m) = rji(m)/dji
         ejk(m) = rjk(m)/djk
         dotj = dotj + eji(m) * ejk(m)
  910 continue
!     ri = dji
!     rk = djk
!     rj = sqrt(dji*djk)
!     rj = dik
      rj=1.
      call vect (eji, ejk, ej, sinj)
      do 920 m = 1, 3
         si(m) = rj*(dotj * eji(m) - ejk(m))/dji/sinj
         sk(m) = rj*(dotj * ejk(m) - eji(m))/djk/sinj
         sj(m) = rj*((dji - djk * dotj)*eji(m) +                        &
     &              (djk - dji * dotj)*ejk(m))/dji/djk/sinj
  920 continue
      return
      end

!-----------------------------------------------------------
      subroutine vect (v1, v2, v3, sin12)

!*calculates vector product v3 = [v1 x v2] and sin of the angle

      implicit double precision (a-h, o-z)

      dimension v1(3), v2(3), v3(3)

      v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
      v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
      d1 = 0.
      d2 = 0.
      d3 = 0.
      do 990 m = 1, 3
         d1 = d1 + v1(m)**2
         d2 = d2 + v2(m)**2
         d3 = d3 + v3(m)**2
  990 continue
      sin12 = sqrt(d3/d1/d2)
      return
      end

!-----------------------------------------------------------
      subroutine coss (v1,v2,v3,cos12)
!* calculates cos between two vectors v1-v2 and v3-v2
      implicit double precision (a-h, o-z)
      dimension v1(3), v2(3), v3(3)

      vv1=0.
      vv2=0.
      scal=0.
      do 995 m=1,3
         vv1=vv1+(v1(m)-v2(m))**2
         vv2=vv1+(v3(m)-v2(m))**2
         scal=scal+(v1(m)-v2(m))*(v3(m)-v2(m))
  995 continue
      cos12=scal/vv1/vv2
      return
      end

!-----------------------------------------------------------
      subroutine chekin (i, natom, line)
      character*150  line
      character*512  slog

            if (i .gt. (natom-1) .or. i .lt. 0) then
               write(slog,'(a,i8)')                                     &
     &             'the atomic indexes must be between 0 and ',         &
     &             (natom - 1)
               call wlog(slog)
               write(slog,'(i8,a)') i, ' not allowed'
               call wlog(slog)
               write(slog,'(1x,a)') line(1:71)
               call wlog(slog)
               call par_stop('RDSPR')
            endif
        return
        end

!---------------------------------------------------------------------
!     program sigrm
!
!     calculate the Debye-Waller factors for each MS path
!     using the recursion method
!
!     input files:  feff.inp and spring.inp
!
!     version 2  ( January 99)
!
!     coded by  A. Poiarkova
!
!---------------------------------------------------------------------
!  References:
!             for the RM: J. Synchrotron Rad., 1999 (to bu published)
!     also see dissertation
!        "X-ray Absorption Fine Structure Debye-Waller Factors"
!         by Anna V. Poiarkova
!
!---------------------------------------------------------------------
!         tk temperature in degrees K
!         nleg  nlegs in path
!         rat   positions of each atom in path
!         NB Units of distance in this routine
!            are angstroms, including sig2.
!         sig2 is output DW factor
!
!---------------------------------------------------------------------
      subroutine sigrm (sig2mx, sig2x,ir1, ir2, tk,ipath,nleg,rat,sig2)
	  use dimsmod,only: nphx=>nphu,natx,legtot  !KJ 7-09 replaces "includ 'dwpar.h'"
      implicit double precision (a-h, o-z)
!         include 'dwpar.h'
!      parameter (nphx = nphx1)
!      parameter (natx = natxdw)

! local parameters:
      parameter (amu0  = 1.66054)
      double precision sig2mx, sig2x(0:nphx,0:nphx)
      dimension iphat(natx), izph(0:nphx)

! variables shared with rdspr.f:
      dimension rat1(3,natx), iz(natx)
      !dimension dm(3,3,natx,natx)
      !dimension rnn(3,natx,natx)
      !dimension nnl(natx,natx) !KJ:
	  real*8, allocatable :: dm(:,:,:,:),rnn(:,:,:)
	  integer,allocatable :: nnl(:,:)

! local variables:
      dimension rat(3,0:legtot) !KJ 7-09 replaced nlegx1 with legtot
      dimension nconv(0:legtot) !KJ 7-09 replaced nlegx1 with legtot
      dimension q0(3,natx)
!   list of atoms |0>=|Q>:
      dimension nq0(0:legtot) !KJ 7-09 replaced nlegx1 with legtot  
!   state |1>=D|Q>:
      dimension q1(3,natx)  
!   list of atoms in |1>:
      dimension nq1(natx)   

!     character*30  fname
      parameter (ntitx1 = 10)
      character*71  title(ntitx1)
      dimension ltit(ntitx1)
!     character*80  titlep(ntitx1)
integer nat !KJ
      character*512 slog

      logical ir1_open, ir2_open

      save 
      data nsigc /0/
!-------------------------------------------------------------

      inquire(unit=ir1,opened=ir1_open)
      inquire(unit=ir2,opened=ir2_open)



      if (nsigc.eq.0) then
!        Read coordinates and potentials from feff.inp
         call dwrdin (rat1, iphat, izph, natom, ntitle, title, ltit)

         if (natom.gt.natx) natom=natx
         do 5 iat=1, natom
            iz(iat) = izph(iphat(iat))
            if (iphat(iat).eq.0) i0=iat
  5      continue

         write(slog,7)
   7     format(2x,'Calculating Debye-Waller factors via RM...')
         call wlog(slog)
         write(slog,9)
   9     format(2x,'This might take a while.')
         call wlog(slog)
  
  nat=natx !KJ 2011
  nat=natom !KJ 2014
  !KJ   I added an allocate statement in 2009.  Now (2011) I am getting errors because the arrays are already allocated.  Perhaps I just copied this from the EM
!KJ   code or whatever, and nobody ever used this method, hence the error went undetected.  Or maybe the allocation is needed for FMS+DEBYE calculations, but not
!KJ   for EXAFS+DEBYE calculations.  I have no idea, but this will work either way:
!      if(.not.allocated(dm)) allocate(dm(3,3,natx,natx))
!	  if(.not.allocated(rnn)) allocate(rnn(3,natx,natx))
!	  if(.not.allocated(nnl)) allocate(nnl(natx,natx))
      if(.not.allocated(dm)) allocate(dm(3,3,nat,nat)) !KJ 2014 code
	  if(.not.allocated(rnn)) allocate(rnn(3,nat,nat))
	  if(.not.allocated(nnl)) allocate(nnl(nat,nat))
      
!        Read spring.inp and build dynamical matrix
         call rdspr(rat1, iz, natom, i0,                                &
     &           dm, rnn,                                               &
     &           acut, res, wmax, dosfit, zshell, w0,                   &
     &           rintr, iprdos, nnl, nat ) !KJ nat

         if (ipath.ne.0) then
            if(ir1_open) then
!             Echo title cards to s2_rm2.dat
              do 10  i = 1, ntitle
                write(ir1,12)  title(i)(1:ltit(i))
  10          continue
  12          format (1x, a)
              write(ir1,17) tk, natom
  17          format(1x,'temperature =',f7.2,2x,'N_at =',i4)
              write(ir1,19)
  19          format (1x, 71('-'))
              write(ir1,25)
              write(slog,25)
  25          format(1x,'ipath',2x,'nleg',4x,'sig2',3x,'mu_ipath',4x,   &
     &          'w_1',6x,'w_2',7x,'A1',5x,'A2')
              call wlog(slog)
            endif
            if (iprdos.ge.0.and.ir2_open) then
!              Echo title cards to s2_rm1.dat
               do 30  i = 1, ntitle
                  write(ir2,12)  title(i)(1:ltit(i))
  30           continue
               write(ir2,17) tk, natom
               write(ir2,19)
               write(ir2,35)
  35           format(1x,'ipath',2x,'nleg',4x,'sig2',3x,'mu_ipath',     &
     &              4x,'w_e')
            endif
         endif
      endif
      nsigc = nsigc + 1
!---- end of first time reading -------

!c    Open path input file (unit in) and read title.  Use unit 2.
!     ntitle2 = 5
!     open(unit=2,file='paths.dat',status='old', iostat=ios)
!     call chopen (ios, 'paths.dat', 'sigrm')
!     call rdhead (2, ntitle2, titlep, ltit)
!c    if (ntitle2 .le. 0)  then
!c       titlep(1) = ' '
!c    endif

! 84  continue
!     read(2,*,end=1010) ipath, nleg
!     skip label (x y z ipot rleg beta eta) and read the path
!     read(2,*)
      do 78 ileg=0,nleg
!        read(2,*,end=1010) (rat(j,ileg),j=1,3)
         nconv(ileg)=0
  78  continue

      do 88 n=1,3
         aa = rat(n,nleg)
         do 87 i=0,(nleg-2)
            j=nleg-i
            rat(n,j)=rat(n,j-1)
  87     continue
         rat(n,1)=aa
  88  continue
      do 89 i=1,nleg
         nq0(i)=0.
  89  continue

! nconv converts # of an atom in the nleg list of coordinates (rat) to
! its # in the full list of all atomic coordinates (rat1)
      do 94 i=1,natom
         do 91 n=1,3
            q1(n,i)=0.
   91    q0(n,i)=0.
            do 95 jl=1,nleg
            m=0
            do 93 n=1,3
               l=nint(100.*rat(n,jl))
               l1=nint(100.*rat1(n,i))
               if (abs(l-l1).le.1) m=m+1
   93       continue
            if (m.eq.3) then
              nconv(jl)=i
              go to 95
            endif
   95    continue
   94 continue

      atmu=0.
      inn=1
      nq1(inn)=1
      iq0=0
      nconv(0)=nconv(nleg)
      do 100 il=1,nleg
         l=nconv(il)
         do 101 jq=1,iq0
  101    if(nq0(jq).eq.l) go to 102
         iq0=iq0+1
         nq0(iq0)=l
  102    continue
         do 105 ii=1,natom
            a = dist(rat1(1,ii),rat1(1,l))
            if (a.le.rintr) then
               do 103 jn=1,inn
  103          if(nq1(jn).eq.ii) go to 105
               inn=inn+1
               nq1(inn)=ii
            endif
  105    continue
         nq1x=inn
         nq0x=iq0
         i=nconv(il)
         im=nconv(il-1)
         ip=nconv(il+1)
!        if (il.eq.1) im=nconv(nleg)
         if (il.eq.nleg) ip=nconv(1)
         atmass=atwtd(iz(i))
      do 100 n=1,3
         atmu=atmu + 0.25*( rnn(n,i,im)+rnn(n,i,ip) )**2 /atmass
  100 continue
      atmu=1./atmu
  108 continue

      do 115 i=1,natom
      do 115 n=1,3
  115 q0(n,i)=0.

! Build initial state vector |Q_j(0)> for the current path
      do 116 n=1,3
      do 116 il=1,nleg
         i=nconv(il)
         im=nconv(il-1)
         ip=nconv(il+1)
         if (il.eq.1) im=nconv(nleg)
         if (il.eq.nleg) ip=nconv(1)
         atmass=atwtd(iz(i))
         q0(n,i)=q0(n,i)+sqrt(atmu/atmass)*(rnn(n,im,i)-rnn(n,i,ip))/2.
  116 continue

! make sure it's normalized <Q_j(0)|Q_j(0)>=1
      q0q0=0.
      do 120 iq0=1,nq0x
      i=nq0(iq0)
      do 120 n=1,3
         q0q0=q0q0+q0(n,i)*q0(n,i)
 120  continue
      p00=nint(q0q0*1000.)/1000.
      if (abs(p00-1.d0).gt.5.d-4) then
         atmu=atmu/q0q0
         go to 108
      endif

!     to get THz units:
      wnorm=100.*w0/sqrt(amu0*10.) 
!*** moments
      a0=0.
      do 132 il=1,nq0x
      do 132 im=1,nq0x
         l=nq0(il)
         m=nq0(im)
      do 132 n1=1,3
      do 132 n2=1,3
         a0 = a0 + q0(n1,l)*dm(n1,n2,l,m)*q0(n2,m)/w0/w0
  132 continue
      we=wnorm*sqrt(a0)
      if (we.lt.1) then
!        recursion method is inapplicable, use statistics to set sig2
         sig2 = sig2x ( iphat(nconv(1)),  iphat(nconv(1)) )
         if (sig2.lt.1.d-6) sig2 = sig2mx
         return
      endif

      do 137 iset=1,nq1x
         i=nq1(iset)
      do 137 n1=1,3
            q1i=0.
            do 138 im=1,nq0x
               m=nq0(im)
            do 138 n2=1,3
               q1i=q1i+dm(n1,n2,i,m)*q0(n2,m)/w0/w0
  138       continue
         q1(n1,i) = q1i - a0*q0(n1,i)
  137 continue

      b0=0.
      do 139 i=1,natom
      do 139 n1=1,3
         b0=b0+q1(n1,i)*q1(n1,i)
  139 continue

      a1=0.
      do 150 iset=1,nq1x
         i=nq1(iset)
         do 150 n1=1,3
         q2=0.
         do 151 jset=1,nq1x
         j=nq1(jset)
            do 151 n2=1,3
               q2 = q2 + dm(n1,n2,i,j)*q1(n2,j)/w0/w0
  151       continue
            a1 = a1 + q1(n1,i)*q2
  150 continue

      a0=a0*wnorm**2
      a1=a1/b0
      a1=a1*wnorm**2
      b0=b0*wnorm**4

!** recursion sigma^2
      dd = (a0+a1)**2 - 4.*(a0*a1-b0)
      x1 = (a0+a1+sqrt(dd))/2.
      x2 = (a0+a1-sqrt(dd))/2.
      aa2 = (a1-x2)/(x1-x2)
!     aa2 = (a1-x2)/(x1-x2)*9./8.
      aa1 = (x1-a1)/(x1-x2)
      w1 = sqrt(x1)
      w2 = sqrt(x2)
      s1 = 3.1746/(atmu*w1*tanh(w1*7.6383/2./tk))
      s2 = 3.1746/(atmu*w2*tanh(w2*7.6383/2./tk))
      sigma2 = aa1*s1+aa2*s2
      sig2e = 3.1746/(atmu*we*tanh(we*7.6383/2./tk))

      if (ipath.ne.0) then
         write(slog,250) ipath,nleg,sigma2,atmu,w1,w2,aa1,aa2
         call wlog(slog)
         if (ir1_open) then
           write(ir1,250) ipath,nleg,sigma2,atmu,w1,w2,aa1,aa2
  250      format(1x,i3,4x,i1,3x,f9.5,2x,f7.3,2x,f7.2,2x,f7.2,          &
     &       4x,f5.3,2x,f5.3)
         endif
         if (iprdos.ge.0.and.ir2_open) then
            write(ir2,260) ipath,nleg,sig2e,atmu,we
  260       format(1x,i3,4x,i3,3x,f7.5,2x,f7.3,2x,f7.2)
         endif
      endif
      sig2 = sigma2

!     update maximum DW factors
      if (sig2.gt.sig2mx) sig2mx=sig2
      if (sig2.gt.sig2x( iphat(nconv(1)),  iphat(nconv(2)) )) then
         sig2x( iphat(nconv(1)),  iphat(nconv(2)) ) = sig2
         sig2x( iphat(nconv(2)),  iphat(nconv(1)) ) = sig2
      endif

      return
      end
!----------------------------------------------------
