!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rexsph.f90,v $:
! $Revision: 1.21 $
! $Author: jorissen $
! $Date: 2012/02/10 05:18:49 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rexsph


!KJ added modules
      use controls
      use struct,only: nphkevin=>nph
      use kklist
      use constants,only: pi,bohr, hart
      use eels_inp
      use reciprocal_inp
      use global_inp
      use xsph_inp
      use atoms_inp
      use nrixs_inp
      use dimsmod, only: nspx=>nspu
      use hubbard_inp
!KJ

      implicit none
!     Local stuff
      integer i,iat
      real*8 celvin  !optional scaling factor for unit cell volume - "hidden option"

      call atoms_read   ! read geom.dat
      call global_read  ! read global.inp
      call potential_read
      call xsph_read    ! read xsph.inp (formerly mod2.inp)
      call eels_read    ! read eels.inp
      call hubbard_read
          ! Josh - only call nrixs_init if xsph is set to run
          IF(mphase.ne.0) call nrixs_init ! initialize some nrixs settings
!KJ next section added for ELNES calculations 1-06
      if(eels.eq.1.and.mphase.eq.1) then
        !call wlog(':INFO : rexsph reduces your polarization tensor to the unit matrix for EELS.')
        ptz(:,:)=dcmplx(0,0)
        do i=-1,1
           ptz(i,i)=dble(1)/dble(3)
        enddo
      endif

! !KJ Next section added for k-space calculations
      nphkevin=nph !KJ temp fix here - need to have same value in both modules!!
      call init_controls
      call reciprocal_read(celvin)  ! read reciprocal.inp
!      call init_struct(nph)

      if(ispace.eq.0) then
           !KJ next lines : initialize nsp in the struct module (old routines will use same value reinitialized by fmstot).
           nsp = 1
           if (abs(ispin).eq.1 ) nsp = nspx
           lpot(0:nph)=lmaxph(0:nph)
           
           ! We only really need to calculate anything in k-space (FMS) if szlz is going to be called in xsphsub.  (2-2012)
           ! For large calculations, it saves time to skip k-space set-up (structure factors etc.)
           makekmeshnow=(ipr2.ge.3)
      endif

      if(ispace.eq.0)  then !KJ set the k-mesh!
         a1=a1/bohr  ! lattice constants in bohr
         a2=a2/bohr
         a3=a3/bohr
         celvin=celvin/(bohr**3)
         call crystalstructure(celvin)

         if(makekmeshnow) then
            call kmesh
         endif

      endif
! !KJ end my changes

!     transform to code units (bohrs and hartrees - atomic units)
      rfms2 = rfms2 / bohr
      vr0   = vr0 / hart
      vi0   = vi0 / hart
      gamach = gamach / hart
      vixan = vixan / hart
      xkstep = xkstep * bohr
      xkmax  = xkmax  * bohr
      EGap = EGap/hart
      do i = 1,3
      do iat = 1, nat
        rat(i,iat) = rat(i,iat) / bohr
      enddo
      enddo

      return
      end
