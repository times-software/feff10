!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: consistency.f90,v $:
! $Revision: 1.14 $
! $Author: jorissen $
! $Date: 2012/03/27 18:15:07 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine consistency_checker(c)
! Added KJ 7-09
! This routine knows which cards have been read from feff.inp by looking at the c array.
! It checks for invalid combinations, eg. two incompatible cards specified, or a missing required card.

! Since we now have so many cards, it would be good to have all these "rules" in one place.
! Also, this setup should be easy to port into the GUI.
! More sophisticated checks (eg. dependent on value of a given card option) should still be done in rdinp.f90


      implicit none
	  logical,intent(in) :: c(150)
!     "c" keeps track of whether a card was listed in feff.inp or no.
!     Check that no incompatible combinations exist.
      integer i(150)
	  integer j
      CHARACTER message(300)

!     set up equivalent array i - integers easier to check groups of exclusive cards
      i=0
      do j=1,150
	     if (c(j)) i(j)=1
	  enddo


!  1/ Not more than one spectroscopy selected (only 1 choice of exafs,exelfs,xanes,elnes,xes,danes,fprime)
      if(i(21)+i(24)+i(42)+i(43)+i(44)+i(56)+i(57).gt.1) stop 'ERROR more than one type of spectroscopy selected'
!  2/ Nrixs must be combined with XANES or EXAFS, but no other card
      if(c(78).and.(i(21)+i(24).ne.1)) stop 'NRIXS must be combined with XANES or EXAFS'
	  if(c(78).and.(i(42)+i(43)+i(44)+i(56)+i(57).gt.0)) stop 'NRIXS combined with incompatible spectroscopy card'	        
!  3/ NRIXS check 2 : This may be overly conservative but I'm worried about overlap and reuse of le2 MULTIPOLE <-> LJMAX
      if((c(79).or.c(80)).and.(.not.c(78))) stop 'LDEC and LJMAX cards only allowed with NRIXS'
      if(c(78).and.c(47)) stop 'you cannot combine NRIXS and MULTIPOLE'
!  4/ NRIXS check 3 :
      if(c(78).and.(i(25)+i(26)+i(29)+i(34)+i(40)+i(46)+i(28)+i(42)+i(49)+i(50)+i(104).gt.0)) then
         
         message = 'The following cards explicitly forbidden for NRIXS : ' // &
      &       'ELLIP,POLARIZATION,NSTAR,SPIN,CFAVERAGE,XNCD,XMCD,RPHASES,TDLDA,XES,PMBSE,HUBBARD'
         stop
      end if
!  5/ k-space needs lattice vectors and k-mesh
      if(c(62)) then
	     if (((i(65)+i(71)).ne.2)) stop 'KMESH and TARGET are required for RECIPROCAL card'
		 if ((i(64)+i(92)).ne.1) stop 'use either LATTICE or CIF with RECIPROCAL card'
	  endif
!  6/ NOHOLE card and COREHOLE card do the same things. JK 08/09
      if(c(30).and.c(68)) stop 'Please use only one of the NOHOLE and COREHOLE cards. They are redundant.'

!!  7/ MDFF needs ELNES or EXELFS
!      if(c(88).and.(.not.(c(56).or.c(57)))) stop 'MDFF must be used with ELNES or EXELFS.'


!  8/ No COMPTON options if compton not enabled
	  if((.not.(c(94).or.c(95))) .and. c(96)) stop 'Cannot use CGRID without COMPTON or RHOZZP.  Exiting.'

!  9/ HUBBARD not compatible with KSPACE
      if( c(104) .and. c(62)) stop 'Cannot use RECIPROCAL with HUBBARD.'

!!! Everybody please add their own checks!




      return
	  end


      logical function have_card(card,c)
      implicit none
	  logical,intent(in) :: c(150)
!     "c" keeps track of whether a card was listed in feff.inp or no.
      character*(*),intent(in) :: card
      integer i
      integer,external :: itoken
      character*4 default_string

      have_card=.false.
      default_string='    '
      !itoken expects a string of at least 4 characters.
      if(len(card).lt.4) then
          default_string=card
          i=itoken(default_string,"feff.inp")
      else
          i=itoken(card,"feff.inp")
      endif
      if (i.gt.0 .and. i.le.150)  have_card=c(i)
      !Note that itoken=0 corresponds to the END card, which is not only useless in the current context, but also doesn't have a corresponding field in "c"

      return
      end





      character*10 function name_spectroscopy(c)
      implicit none
	  logical,intent(in) :: c(150)
      character*10 s
      logical,external :: have_card

      s='          '
      if(have_card('XANES',c)) s='XANES'
      if(have_card('EXAFS',c)) s='EXAFS'
      if(have_card('ELNES',c)) s='ELNES'
      if(have_card('EXELFS',c)) s='EXELFS'
      if(have_card('COMPTON',c)) s='COMPTON'
      if(have_card('NRIXS',c)) s='NRIXS'
      if(have_card('RIXS',c)) s='RIXS'
      if(have_card('XES',c)) s='XES'
      if(have_card('FPRIME',c)) s='FPRIME'
      if(have_card('XMCD',c)) s='XMCD'
!      if(have_card('XNCD',c)) s='XNCD'
!     XNCD and XMCD have the same "itoken" value and are the same card as far as FEFF is concerned.
!     Always return "XMCD" since this is the more common option.
      if(have_card('DANES',c)) s='DANES'

      ! Annoyingly, FEFF allows a default of EXAFS:
      if(s.eq.'          ') s='EXAFS'

      name_spectroscopy=s
      return
      end





      character*10 function name_corehole(i)
      implicit none
      integer,intent(in) :: i

      name_corehole='FSR       '  !default
      if(i.eq.2) name_corehole='RPA'
      if(i.eq.0) name_corehole='no'
      return
      end





      character*3000 function list_cards(c)
      implicit none
	  logical,intent(in) :: c(150)
      character*3000 cardlist
      character*20,external :: itoken_reverse
      integer i

      cardlist=' '
      do i=1,150
        if(c(i)) cardlist=trim(cardlist)//' '//trim(itoken_reverse("feff.inp",i))//' '
      enddo
      list_cards=cardlist
      return
      end

      




      character*3000 function list_features(c)
      implicit none
	  logical,intent(in) :: c(150)
      character*3000 featurelist
      character*20,external :: itoken_reverse
      character*100 feature_card,feature_prose
      dimension feature_card(1:11),feature_prose(1:11)
      logical,external :: have_card
      integer i
      data feature_card /'DEBYE','MPSE','TDLDA','SPIN','SCF','UNFREEZEF','PMBSE','S02CONV','EXTPOT','CONFIG','TEMP'/
      data feature_prose &
        /'Debye-Waller factors','Many-Pole Self-Energy', &
         'Time-Dependent Density-Functional-Theory','Spin Polarization', &
         'Self-Consistent Field potentials','SCF-converged f-states', &
         'approximated Bethe-Salpeter cross-section', &
         'Satellite Spectral Function','External Potentials', &
         'Custom electron configuration','Finite Temperature Fermi Distribution'/

      featurelist=' '
      do i=1,11
        if(have_card(feature_card(i),c)) featurelist=trim(featurelist)//'   * '//trim(feature_prose(i))//' '
      enddo
      list_features=featurelist
      return
      end





!     Copied from COMMON/itoken.f90 on 7-09 :
! ATOM = 1
! HOLE = 2
! OVER = 3
! CONT = 4
! EXCH = 5
! ION  = 6
! TITL = 7
! FOLP = 8
! RPAT/RMAX = 9
! DEBY = 10
! RMUL = 11
! SS =   12
! PRIN = 13
! POTE = 14
! NLEG = 15
! CRIT = 16
! NOGE = 17
! IORD = 18
! PCRI = 19
! SIG2 = 20
! XANE = 21
! CORR = 22
! AFOL = 23
! EXAF = 24
! POLA = 25
! ELLI = 26
! RGRI = 27
! RPHA = 28
! NSTA = 29
! NOHO = 30
! SIG3 = 31
! JUMP = 32
! MBCO = 33
! SPIN = 34
! EDGE = 35
! SCF  = 36
! FMS  = 37
! LDOS = 38
! INTE = 39
! CFAV = 40
! S02  = 41
! XES  = 42
! DANE = 43
! FPRI = 44
! RSIG = 45
! XNCD = 46
! XMCD = 46
! MULT = 47
! UNFR = 48
! TDLD = 49
! PMBS = 50
! PLAS = 51
! SO2C = 52
! SELF = 53
! SFSE = 54
! RCONV = 55
! ELNE = 56
! EXEL = 57
! MAGI = 58
! ABSO = 59
! SYMM = 60
! REAL = 61
! RECI = 62
! SGRO = 63
! LATT = 64
! KMES = 65
! STRF = 66
! BAND = 67
! TARG = 71
! EGRI = 72
! COOR = 73	
! EXTP = 74
! CHBR = 75
! CHSH = 76
! DIMS = 77
! NRIX = 78
! LJMA = 79
! LDEC = 80
