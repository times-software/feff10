!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: getorb.f90,v $:
! $Revision: 1.8 $
! $Author: jorissen $
! $Date: 2013/01/18 02:19:08 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine getorb (iz, ihole, xion, iunf, norb, norbco, iorb, iholep, nqn, nk, xnel, xnval, xmag, iph)  !KJ 12-2010 added iph
!     Gets orbital data for chosen element.  Input is:
!       iz - atomic number of desired element,
!       ihole - index of core-hole orbital
!       xion  - ionicity (usually zero)
!     other arguments are output.
!       norb - total number of orbitals
!       norbco - number of core orbitals
!       iorb - index of orbital for making projections (last occupied)
!       iholep - index of core hole orbital in compacted list
!       nqn - principal quantum number for each orbital
!       nk - quantum number kappa for each orbital
!       xnel - occupation for each orbital
!       xnval - valence occupation for each orbital
!       xmag - spin magnetization for each orbital
      use config,izdum=>iz  ! replaces old data statements
	  
	  ! KJ bugfix Sept 2012 : added initialization statement setting xnel,xmag,xnval to 0 since getorb is called many times ...  Yikes!
	  
	  implicit none
	  
      integer,intent(in)  :: iz, ihole, iph, iunf
	  real*8,intent(in)   :: xion
	  integer,intent(out) :: norb, norbco, iorb(-4:3), iholep, nqn(30), nk(30)
	  real*8,intent(out)  :: xnel(30), xnval(30), xmag(30)

!     Written by Steven Zabinsky, July 1989
!     modified (20 aug 1989)  table increased to at no 100
!     Recipe for final state configuration is changed. Valence electron occupations are added. ala 17.1.1996
!     Structure changed - using modules, allowing user input.  KJ 12.2010

      character*512 slog
	  integer ion,index,ilast,iscf,iion,i,index1,iscr,iphl
	  real*8 delion



!  How an ION is treated :
!  Say we have ION 1 0.6 in feff.inp, where potential "1" is for Fe (Z=26), i.e. we have a Fe^+0.6 atom
!  We need to remove 0.6 e charge from this Fe atom.
!  We will start from the default Z=25 configuration of Mn.  
!  Then, we will find the orbital for which the occupation differs between the default Mn and Fe occupations.
!  In this orbital, we will add 0.4 e.  (Unless there is a corehole and the screening electron would already fill all available
!  space in this orbital; in that case, we choose the highest orbital with any occupation; or, if that one too lacks sufficient space, the first empty orbital.)
!  The parameters used below will evaluate to (for a NOHOLE calculation) :
!  (iz=26; index=25; index1=26; xion=0.6; ion=1; delion=-0.4; iion=9 [3d]; ilast=11; iscr=9 )

!  If instead we had ION 1 -0.6, then we would start from the default Z=26 configuration for Fe.
!  We would find the highest orbital with sufficient occupation to remove 0.6 e.
!  And we'd just take them out :).


!     write(*,*) '########## START GETORB ##########'
!     write(*,*) 'iz,ihole,xion',iz,ihole,xion
!	  write(*,*) 'iunf,norb,norbco',iunf,norb,norbco
!	  write(*,*) 'iorb',iorb
!	  write(*,*) 'iholep,nqn',iholep,nqn
!	  write(*,*) 'nk',nk
!	  write(*,*) 'xnel',xnel
!	  write(*,*) 'xnval',xnval
!	  write(*,*) 'xmag,iph',xmag,iph
!	  write(*,*) '       ######'


      call InitConfig   !Initialize default electronic configurations, taking user input into account.
      if (iz .lt. 1  .or.  iz .gt. 100)  then
         write(slog,'(" Atomic number ",i5," not available.")')  iz
         call wlog(slog)
         call par_stop('GETORB-0')
      endif
!	  if (iph.eq.1) then
!	  write(*,*) "O configuration iocc:"
!	  write(*,'(100f7.2)') iocc(1,15:24)
!	  endif


! ###### The following section figures out which orbitals will accommodate the core hole, screening electron, and ionization :
!        (i.e., find iscr and iion)

      ion = nint(xion)
      delion=xion-ion  !we will subtract delion electrons from the base configuration Z=index
      index = iz - ion
      ilast = 0
      iscr = 0
      iion = 0
      iholep = ihole
      if(index.eq.iz) then !KJ
         iphl=iph  !=> use potential index to retrieve occupation numbers
      else
         iphl=-743 !negative number  => use modified (ionized) atomic number to retrieve atomic numbers
      endif


!     find last occupied orbital (ilast) and iion for delion.ge.0
!     (if delion<0, this produces nonsense, but it will be fixed below, once we know iscr.)
      do i=29,1,-1
         if (iion.eq.0  .and. f_iocc(index,i,iphl).gt.delion) iion=i
         if (ilast.eq.0 .and. f_iocc(index,i,iphl).gt.0) ilast=i
      enddo
!     write(*,*) 'iion,ilast',iion,ilast
      !check that the core hole orbital contains sufficient charge to create a hole:
      if (ihole.gt.0) then
         if ( f_iocc(index,ihole,iphl) .lt. 1 ) then
           call wlog(' Cannot remove an electron from this level')
           call par_stop('GETORB-1')
         endif
      endif
	  !even if we are also ionizing that same orbital:
      if (ihole.eq.iion .and. delion.gt.0) then   !KJ ilast->iion 5-2012; if delion<0 we will add electrons, so there's no need to worry here, but skip the test, since iion will be set to a nonsense value.
         if ( f_iocc(index,ihole,iphl)-delion.lt.1) then
           call wlog(' Cannot remove an electron from this level')
           call par_stop('GETORB-1')
        endif
      endif

!        the recipe for final state atomic configuration is changed from iz+1 prescription, since sometimes it changed occupation
!        numbers in more than two orbitals. This could be consistent only with s02=0.0. New recipe remedies this deficiency.

!     find where to put screening electron : put it in the orbital where the "Z+1" atom has its extra electron
      index1 = index + 1
      do i = 1, 29
         if (iscr.eq.0 .and. (f_iocc(index1,i,-1)-f_iocc(index,i,iphl)).gt.0.5) iscr=i  
      enddo

!     If core-hole orbital only has one electron, set iscr to ihole - Josh Kas
      if(ihole.gt.0) then
	     if (f_iocc(index,ihole,iphl).lt.1.5) iscr = ihole !KJ array bounds get exceeded here when ihole=0 - fixed!
	  endif
!     special case of hydrogen like ion
!     if (index.eq.1) iscr=2

!     find where to add or subtract charge delion (iion).
!     if (delion .ge. 0) then
!        removal of electron charge
!        iion is already found
      if (delion .lt. 0) then
!        addition of electron charge in the amount of -delion to the same orbital where we'd put the screening electron
         iion = iscr
!        except special cases where there's not enough space for the ionization charge
         if (ihole.ne.0 .and. f_iocc(index,iscr,iphl)+1-delion.gt.2*abs(kappa(iscr))) then
             iion = ilast
             if (ilast.eq.iscr .or. f_iocc(index,ilast,iphl)-delion.gt. 2*abs(kappa(ilast)) ) iion = ilast + 1
         endif
      endif

!      if (iph.eq.1)      write(*,*) 'iscr,iion,index,ihole,delion',iscr,iion,index,ihole,delion
!      write(*,*) 'f_iocc(index,i,iphl',(f_iocc(index,i,iphl),i=1,29)
! 	   write(*,*) 'v_ival(index,i,iphl',(f_ival(index,i,iphl),i=1,29)
	  
	  
! ######## Now we know everything.  Start filling up the occupation arrays.	 

! Note that xnel,xnval,xmag do not track all orbitals, but only the ones containing charge.  e.g. "2 2 4 1 0 0 0 2" --> "2 2 4 1 2" 
      norb = 0
      iorb(-4:3) = 0
	  xnel(:)=0.d0
	  xnval(:)=0.d0
	  xmag(:)=0.d0
      do i = 1, 29
! Modified by FDV
! Split line to avoid error in Solaris Studio
         if (f_iocc(index,i,iphl).gt.0 .or. (i.eq.iscr .and. ihole.gt.0) &
             .or. (i.eq.iion .and. f_iocc(index,i,iphl)-delion.gt.0) )  then
		     !the template has e here        !put the screening e here         !we put ionization charge here
            if (i.ne.ihole .or. f_iocc(index,i,iphl).ge.1) then
               norb = norb + 1
               nqn(norb) = nnum(i) !principal quantum number
               nk(norb)  = kappa(i) !relativistic quantum number
               xnel(norb) = f_iocc(index,i,iphl) !occupation of the template
               if (i.eq.ihole) then   !create the core hole
                  xnel(norb) = xnel(norb) - 1
                  iholep = norb
               endif
               if (i.eq.iscr .and. ihole.gt.0)  xnel(norb)=xnel(norb)+1  !add the screening electron
               xnval(norb)= f_ival(index,i,iphl) !valence occupation
               if ((kappa(i).eq.-4 .or. kappa(i).eq.3) .and. iunf.eq.0)  xnval(norb) = 0 !put f-electron in the core, i.e. "freeze" them
               xmag(norb) = f_ispn(index,i,iphl)
               iorb(nk(norb)) = i
               if (i.eq.ihole .and. xnval(norb).ge.1)   xnval(norb) = xnval(norb) - 1 !adjust valence occupation for core hole
               if (i.eq.iscr .and. ihole.gt.0) xnval(norb) = xnval(norb) + 1 !adjust valence occupation for screening electron
               if (i.eq.iion)  xnel(norb) = xnel(norb) - delion !adjust occupation for ionization  !KJ 5-2012 iscr-> iion bugfix
               if (i.eq.iion)  xnval(norb) = xnval(norb) - delion !adjust valence occupation for ionization
            endif
         endif
      enddo
      norbco = norb

!     check that all occupation numbers are within limits
      do i = 1, norb
         if ( xnel(i).lt.0 .or.  xnel(i).gt.2*abs(nk(i)) .or. xnval(i).lt.0 .or. xnval(i).gt.2*abs(nk(i)) ) then
            write (slog,55) i
   55       format(' error in getorb.f. Check occupation number for ',i3, '-th orbital. May be a problem with ionicity.')
            call wlog(slog)
            call par_stop('GETORB-99')
         endif
      enddo
	  
!	  if (iph.eq.1) write(*,'(100f7.2)') xnel(15:24)
!     if (iph.eq.1) stop

!	  write(*,*) '       ######'
!     write(*,*) 'iz,ihole,xion',iz,ihole,xion
!	  write(*,*) 'iunf,norb,norbco',iunf,norb,norbco
!	  write(*,*) 'iorb',iorb
!	  write(*,*) 'iholep,nqn',iholep,nqn
!	  write(*,*) 'nk',nk
!	  write(*,*) 'xnel',xnel
!	  write(*,*) 'xnval',xnval
!	  write(*,*) 'xmag,iph',xmag,iph
!	  write(*,*) '       ######'
!     write(*,*) '########## END GETORB ##########'

            
      return
      end  ! subroutine getorb

      ! JK - added getspin to use for XMCD when no spin occupations
      !      are defined in potentials input. In this case it will 
      !      just use the sum of spin occupation as defined in the
      !      atomic configuration.
      real*8 function getspin (iz, ihole, xion, iph)
          IMPLICIT NONE
          integer,intent(in)  :: iz, ihole
	  real*8,intent(in)   :: xion
          integer iph, iunf, i
	  integer :: norb, norbco, iorb(-4:3), iholep, nqn(30), nk(30)
	  real*8  :: xnel(30), xnval(30), xmag(30)
          iunf = 0
          xmag = 0.d0
          call getorb (iz, ihole, xion, iunf, norb, norbco, iorb, iholep, nqn, nk, xnel, xnval, xmag, iph)  
          getspin=0.0
          DO i = 1, 30
             getspin = getspin + xmag(i)
          END DO
      end function getspin
