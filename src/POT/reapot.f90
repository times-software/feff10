!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: reapot.f90,v $:
! $Revision: 1.22 $
! $Author: jorissen $
! $Date: 2012/05/30 00:55:55 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine reapot(run_kpreppot) !KJ put everything in modules 7-09

        use controls
        use struct,nphkevin=>nph
        use kklist
        use strfacs
        use constants
		use potential_inp
		use atoms_inp
		use reciprocal_inp
        implicit none

!     Local stuff
      logical,intent(in) :: run_kpreppot
      character*32 s1, s2, s3
	  real*8 celvin !KJ scaling factor for volume of unit cell - hidden option
	  integer i,iat,istart,iend,ilen,iph,iovr
	  integer,external :: istrln


      call atoms_read !read geom.dat
      call potential_read !read pot.inp
	  
	  if(mpot.le.0) return  !KJ 2-2012

      !KJ Next section added for k-space calculations
      nphkevin=nph !KJ temp fix here - need to have same value in both modules!!
      call init_controls
      call reciprocal_read(celvin)  ! read reciprocal.inp
	  if (nohole.ge.0) corehole=.false. !KJ if NOHOLE 2, corehole needs to be used in fms but *not* in pot !!

      if(ispace.eq.0) then
		   !KJ next lines : initialize nsp in the struct module (old routines will use same value reinitialized by fmstot).
           nsp = 1
           lpot(0:nph)= lmaxsc(0:nph)  ! copy lmaxsc into lpot
      endif

      makekmeshnow=.true.
      if(ispace.eq.0 .and. run_kpreppot )  then !KJ prepare the k-mesh
	     if(ktype.eq.2 .or. ktype.eq.3) then
		    !KJ 1-2012 potentials optionally run with k-mesh that's 5 times smaller than k-mesh for fms/ldos.
			! not possible to specify x/y/z k-mesh.
		    nkp=nkp/5
			nkx=nkp
			nky=0
			nkz=0
		 endif
         a1=a1/bohr  ! lattice constants in bohr
         a2=a2/bohr
         a3=a3/bohr
         celvin=celvin/(bohr**3)
         call crystalstructure(celvin)

		 if(nscmt.gt.0) then  !If no SCF we don't need to set up structure constants - saves loads of time for large unit cells  KJ 2-2012
             if(makekmeshnow) then
                call kmesh

             else
    !           klist.inp  !KJ added 8/06
                open(3,file='klist.inp',form='formatted',status='old')
                read(3,*)
                read(3,*) nkp,usesym,ktype
                call init_kklist(nkp,nsym) !KJ 6-09
                read(3,*)
                do i=1,nkp
                   read(3,*) bk(:,i),weight(i)
                enddo
                close(3)
    !           k-mesh in fractional units 0-1 for strfac
                call wlog('assuming a klist.inp in Angstrom')
                call wlog('update these instructions !!')
                do i=1,3
                   bk(i,:)=bk(i,:) * alat(i)/(dble(2)*pi) *bohr  !remove *bohr if klist.inp is in a.u.
                enddo
             endif

             !KJ kprep called here in mod1 (not in mod2 and mod3) :
             if (run_kpreppot) call kpreppot   ! executed for "pot" but not for "atomic"
		 endif ! nscmt > 0 (i.e. do SCF)
      endif


!     transform to code units (bohrs and hartrees - atomic units)
      rfms1 = rfms1 / bohr
      rfms1_start = rfms1_start/bohr
      gamach = gamach / hart
      ecv   = ecv   / hart
      totvol = totvol / bohr**3
      do 210 iat = 1, nat
      do 210 i = 1,3
        rat(i,iat) = rat (i, iat) / bohr
  210 continue
      do 220 iph = 0, nph
      do 220 iovr = 1, novr(iph)
         rovr(iovr,iph) = rovr(iovr,iph) / bohr
  220 continue

!     add lines to the title
      if (mpot.eq.1) then
         ntitle = ntitle + 1
         if (nat.gt.1) then
           if (rfms1.lt.0) rfms1 = 0
           if (nscmt.gt.0) then
             write(s1, 230) nscmt, rfms1*bohr, lfms1
  230        format(' POT  SCF', i4, f8.4, i4)
           else
             write(s1, 235) 
  235        format(' POT  Non-SCF' )
           endif
         else
           write(s1, 240) 
  240      format(' POT  used OVERLAP geometry,')
         endif
         if (nohole.eq.0) then
           write(s2, 310) 
  310      format(', NO core-hole,')
         elseif (nohole.eq.2) then
           write(s2, 315) 
  315      format(', screened core-hole,')
         else
           write(s2, 320) 
  320      format(', core-hole,')
         endif
         if (iafolp.lt.0) then
           write(s3, 330) folp(0)
  330      format(' FOLP (folp(0)=', f6.3, ')' )
         else
           write(s3, 340) folp(0)
  340      format(' AFOLP (folp(0)=', f6.3, ')' )
         endif
!        concatenate 3 strings into 1
         title(ntitle) = ' '
         ilen = istrln(s1)
         istart = 1
         iend = ilen
         title(ntitle)(istart:iend) = s1(1:ilen)
         ilen = istrln(s2)
         istart = iend + 1
         iend = iend + ilen
         title(ntitle)(istart:iend) = s2(1:ilen)
         ilen = istrln(s3)
         istart = iend + 1
         iend = iend + ilen
         title(ntitle)(istart:iend) = s3(1:ilen)
      endif

      return
      end
