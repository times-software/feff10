!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: itoken.f90,v $:
! $Revision: 1.29 $
! $Author: jorissen $
! $Date: 2013/01/07 15:19:46 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function itoken (word,flname)
!     chars in word assumed upper case, left justified
!     returns 0 if not a token, otherwise returns token

      character*(*) word
      character*4   w
      character*(*) flname
      integer itoken

      w = word(1:4)
      call upper(w)

!     Tokens for feff.inp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (flname(1:8).eq.'feff.inp') then
         if     (w .eq. 'ATOM')  then
            itoken = 1
         elseif (w .eq. 'HOLE')  then
            itoken = 2
         elseif (w .eq. 'OVER')  then
            itoken = 3
         elseif (w .eq. 'CONT')  then
            itoken = 4
         elseif (w .eq. 'EXCH')  then
            itoken = 5
         elseif (w .eq. 'ION ')  then
            itoken = 6
         elseif (w .eq. 'TITL')  then
            itoken = 7
         elseif (w .eq. 'FOLP')  then
            itoken = 8
         elseif (w .eq. 'RPAT' .or. w .eq. 'RMAX')  then
            itoken = 9
         elseif (w .eq. 'DEBY')  then
            itoken = 10
         elseif (w .eq. 'RMUL')  then
            itoken = 11
         elseif (w .eq. 'SS  ')  then
            itoken = 12
         elseif (w .eq. 'PRIN')  then
            itoken = 13
         elseif (w .eq. 'POTE')  then
            itoken = 14
         elseif (w .eq. 'NLEG')  then
            itoken = 15
         elseif (w .eq. 'CRIT')  then
            itoken = 16
         elseif (w .eq. 'NOGE')  then
            itoken = 17
         elseif (w .eq. 'IORD')  then
            itoken = 18
         elseif (w .eq. 'PCRI')  then
            itoken = 19
         elseif (w .eq. 'SIG2')  then
            itoken = 20
         elseif (w .eq. 'XANE')  then
            itoken = 21
         elseif (w .eq. 'CORR')  then
            itoken = 22
         elseif (w .eq. 'AFOL')  then
            itoken = 23
         elseif (w .eq. 'EXAF')  then
            itoken = 24
         elseif (w .eq. 'POLA')  then
            itoken = 25
         elseif (w .eq. 'ELLI')  then
            itoken = 26
         elseif (w .eq. 'RGRI')  then
            itoken = 27
         elseif (w .eq. 'RPHA')  then
            itoken = 28
         elseif (w .eq. 'NSTA')  then
            itoken = 29
         elseif (w .eq. 'NOHO')  then
            itoken = 30
         elseif (w .eq. 'SIG3')  then
            itoken = 31
         elseif (w .eq. 'JUMP')  then
            itoken = 32
         elseif (w .eq. 'MBCO')  then
            itoken = 33
         elseif (w .eq. 'SPIN')  then
            itoken = 34
         elseif (w .eq. 'EDGE')  then
            itoken = 35
         elseif (w .eq. 'SCF ')  then
            itoken = 36
         elseif (w .eq. 'FMS ')  then
            itoken = 37
         elseif (w .eq. 'LDOS')  then
            itoken = 38
         elseif (w .eq. 'INTE')  then
            itoken = 39
         elseif (w .eq. 'CFAV')  then
            itoken = 40
         elseif (w .eq. 'S02 ')  then
            itoken = 41
         elseif (w .eq. 'XES ')  then
            itoken = 42
         elseif (w .eq. 'DANE')  then
            itoken = 43
         elseif (w .eq. 'FPRI')  then
            itoken = 44
         elseif (w .eq. 'RSIG')  then
            itoken = 45
         elseif (w .eq. 'XNCD')  then
            itoken = 46
         elseif (w .eq. 'XMCD')  then
            itoken = 46
         elseif (w .eq. 'MULT')  then
            itoken = 47
         elseif (w .eq. 'UNFR')  then
            itoken = 48
         elseif (w .eq. 'TDLD')  then
            itoken = 49
         elseif (w .eq. 'PMBS')  then
            itoken = 50
         elseif (w .eq. 'PLAS' .or. w .eq. 'MPSE')  then
            itoken = 51
         elseif (w .eq. 'SO2C' .or. w .eq. 'SFCO')  then
            itoken = 52
         elseif (w .eq. 'SELF')  then
            itoken = 53
         elseif (w .eq. 'SFSE')  then
            itoken = 54
         elseif (w .eq. 'RCONV') then
            itoken = 55
         elseif (w .eq. 'ELNE') then !KJ new card for EELS 1-06
            itoken = 56
         elseif (w .eq. 'EXEL') then !KJ new card for EELS 1-06
            itoken = 57
         elseif (w .eq. 'MAGI') then !KJ new card for EELS 1-06
            itoken = 58
         elseif (w .eq. 'ABSO') then !KJ new card 3-06
            itoken = 59  
         elseif (w .eq. 'SYMM') then !KJ new card 6-06
            itoken = 60  
         elseif (w .eq. 'REAL') then  !KJ 8/06
              itoken = 61
         elseif (w .eq. 'RECI') then  !KJ 8/06
              itoken = 62
         elseif (w .eq. 'SGRO') then  !KJ 8/06
              itoken = 63
         elseif (w .eq. 'LATT') then  !KJ 8/06
              itoken = 64
         elseif (w .eq. 'KMES') then  !KJ 8/06
              itoken = 65
         elseif (w .eq. 'STRF') then  !KJ 8/06
              itoken = 66
         elseif (w .eq. 'BAND') then  !KJ 8/06
              itoken = 67
         elseif (w .eq. 'CORE') then  !JK 8/09
              itoken = 68
         elseif (w .eq. 'MARK' .or. w .eq. 'TARG') then  !KJ 8/06 !KJ 12-2010
            itoken = 71
         elseif (w .eq. 'EGRI') then  !KJ 1/07
            itoken = 72
         elseif (w .eq. 'COOR') then   !KJ 4/07
            itoken = 73	    
         elseif (w .eq. 'EXTP') then  !JK 4/19/08
            itoken = 74
         elseif (w .eq. 'CHBR') then !JK 4/19/08
            itoken = 75
         elseif (w .eq. 'CHSH') then ! Added by Fer
            itoken = 76
         elseif (w .eq. 'DIMS') then !JPR 4/20/09 renumbered KJ
            itoken = 77
		 elseif (w .eq. 'NRIX') then !KJ 7/09
		    itoken = 78
		 elseif (w .eq. 'LJMA') then !KJ 7/09
		    itoken = 79
		 elseif (w .eq. 'LDEC') then !KJ 7/09
		    itoken = 80
         elseif (w .eq. 'SETE') then ! JJK 1/2010
            itoken = 81
         elseif (w .eq. 'EPS0') then ! JJK 3/2010
            itoken = 82
         elseif (w .eq. 'OPCO') then ! JJK 3/2010
            itoken = 83
         elseif (w .eq. 'NUMD') then ! JJK 3/2010
            itoken = 84
         elseif (w .eq. 'PREP') then ! JJK 3/2010
            itoken = 85 
         elseif (w .eq. 'EGAP') then ! JJK 4/2010
            itoken = 86
         elseif (w .eq. 'CHWI') then ! KJ 6/2010
            itoken = 87 
		 elseif (w .eq. 'MDFF') then ! KJ 11/2010
		    itoken = 88
		 elseif (w .eq. 'REST') then !KJ 12/2010
		    itoken = 89
         elseif (w .eq. 'CONF') then !KJ 12/2010
		    itoken = 90
		 elseif (w .eq. 'SCRE') then !KJ 7/2011
		    itoken = 91
		 elseif (w .eq. 'CIF ') then !KJ 10/2011
		    itoken = 92
		 elseif (w .eq. 'EQUI') then !KJ 1/2012
		    itoken = 93
		 elseif (w .eq. 'COMP') then !BAM 2/2012
		    itoken = 94
		 elseif (w .eq. 'RHOZ') then !BAM 2/2012
		    itoken = 95
		 elseif (w .eq. 'CGRI') then !BAM 2/2012
		    itoken = 96
		 elseif (w .eq. 'CORV') then !KJ 10/2012
		    itoken = 97
		 elseif (w .eq. 'SIGG') then !KJ 01/2013
		    itoken = 98
		 elseif (w .eq. 'TEMP') then !BAM 2/2013
		    itoken = 99
		 elseif (w .eq. 'DENS') then !BAM 2/2013
		    itoken = 100
		 elseif (w .eq. 'RIXS') then !JJK 5/8/2014
		    itoken = 101
		 elseif (w .eq. 'RLPR') then !JJK 5/8/2014
		    itoken = 102
		 elseif (w .eq. 'ICOR') then !JJK 5/8/2014
		    itoken = 103
		 elseif (w .eq. 'HUBB') then !CV 11/2013
		    itoken = 104
		 elseif (w .eq. 'CRPA') then !CV 04/2014
                    itoken = 105
         elseif (w .eq. 'FULL') then !KJ 09/2014
            itoken = 106
         elseif (w .eq. 'SCXC') then !LC 03/2015
            itoken = 107
         elseif (w .eq. 'HIGH') then !JK 9/2020
            itoken = 108
         elseif (w .eq. 'SCFT') then !TS 07/2020
            itoken = 109
         elseif (w .eq. 'WARN') then !JK 5/2025
            itoken = 110
         elseif (w .eq. 'END ') then
            itoken = -1            
         else
            itoken = 0
         endif
      elseif (flname(1:10).eq.'spring.inp') then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     These tokens are for spring.inp (input for eq of motion method)
         if (w .eq. 'STRE')  then
            itoken = 1
         elseif (w .eq. 'ANGL')  then
            itoken = 2
         elseif (w .eq. 'VDOS')  then
            itoken = 3
         elseif (w .eq. 'PRDO' .or. w .eq. 'PRIN') then
            itoken = 4
         elseif (w .eq. 'END ')  then
            itoken = -1            
         else
            itoken = 0
         endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      endif
      
      
      return
      end


     character*20 function itoken_reverse (flname,itoken)
      character*20   w
      character*(*),intent(in) :: flname
      integer,intent(in) :: itoken

      w = '                    '

!     Tokens for feff.inp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (flname(1:8).eq.'feff.inp') then
            if(itoken.eq.1) w='ATOMS'
            if(itoken.eq.2) w='HOLE'
            if(itoken.eq.3) w='OVERLAP'
            if(itoken.eq.4) w='CONTROL'
            if(itoken.eq.5) w='EXCHANGE'
            if(itoken.eq.6) w='ION'
            if(itoken.eq.7) w='TITLE'
            if(itoken.eq.8) w='FOLP'
            if(itoken.eq.9) w='RPATH'
            if(itoken.eq.10) w='DEBYE'
            if(itoken.eq.11) w='RMULT'
            if(itoken.eq.12) w='SS'
            if(itoken.eq.13) w='PRINT'
            if(itoken.eq.14) w='POTENTIALS'
            if(itoken.eq.15) w='NLEG'
            if(itoken.eq.16) w='CRITERIA'
            if(itoken.eq.17) w='NOGEOM'
            if(itoken.eq.18) w='IORD'
            if(itoken.eq.19) w='PCRITERIA'
            if(itoken.eq.20) w='SIG2'
            if(itoken.eq.21) w='XANES'
            if(itoken.eq.22) w='CORRECTIONS'
            if(itoken.eq.23) w='AFOLP'
            if(itoken.eq.24) w='EXAFS'
            if(itoken.eq.25) w='POLARIZATION'
            if(itoken.eq.26) w='ELLIPTICITY'
            if(itoken.eq.27) w='RGRID'
            if(itoken.eq.28) w='RPHASES'
            if(itoken.eq.29) w='NSTAR'
            if(itoken.eq.30) w='NOHOLE'
            if(itoken.eq.31) w='SIG3'
            if(itoken.eq.32) w='JUMPRM'
            if(itoken.eq.33) w='MBCONV'
            if(itoken.eq.34) w='SPIN'
            if(itoken.eq.35) w='EDGE'
            if(itoken.eq.36) w='SCF'
            if(itoken.eq.37) w='FMS'
            if(itoken.eq.38) w='LDOS'
            if(itoken.eq.39) w='INTERSTITIAL'
            if(itoken.eq.40) w='CFAVERAGE'
            if(itoken.eq.41) w='S02'
            if(itoken.eq.42) w='XES'
            if(itoken.eq.43) w='DANES'
            if(itoken.eq.44) w='FPRIME'
            if(itoken.eq.45) w='RSIGMA'
            if(itoken.eq.46) w='XMCD'
            if(itoken.eq.47) w='MULT'
            if(itoken.eq.48) w='UNFREEZEF'
            if(itoken.eq.49) w='TDLDA'
            if(itoken.eq.50) w='PMBSE'
            if(itoken.eq.51) w='MPSE'
            if(itoken.eq.52) w='SFCONV'
            if(itoken.eq.53) w='SELF'
            if(itoken.eq.54) w='SFSE'
            if(itoken.eq.55) w='RCONV'
            if(itoken.eq.56) w='ELNES'
            if(itoken.eq.57) w='EXELFS'
            if(itoken.eq.58) w='MAGIC'
            if(itoken.eq.59) w='ABSOLUTE'
            if(itoken.eq.60) w='SYMMETRY'
            if(itoken.eq.61) w='REAL'
            if(itoken.eq.62) w='RECIPROCAL'
            if(itoken.eq.63) w='SGROUP'
            if(itoken.eq.64) w='LATTICE'
            if(itoken.eq.65) w='KMESH'
            if(itoken.eq.66) w='STRFAC'
            if(itoken.eq.67) w='BAND'
            if(itoken.eq.68) w='COREHOLE'
            if(itoken.eq.71) w='TARGET'
            if(itoken.eq.72) w='EGRID'
            if(itoken.eq.73) w='COORDINATES'
            if(itoken.eq.74) w='EXTPOT'
            if(itoken.eq.75) w='CHBROADENING'
            if(itoken.eq.76) w='CHSHIFT'
            if(itoken.eq.77) w='DIMS'
            if(itoken.eq.78) w='NRIXS'
            if(itoken.eq.79) w='LJMAX'
            if(itoken.eq.80) w='LDECMX'
            if(itoken.eq.81) w='SETE'
            if(itoken.eq.82) w='EPS0'
            if(itoken.eq.83) w='OPCONS'
            if(itoken.eq.84) w='NUMD'
            if(itoken.eq.85) w='PREP'
            if(itoken.eq.86) w='EGAP'
            if(itoken.eq.87) w='CHWIDTH'
            if(itoken.eq.88) w='MDFF'
            if(itoken.eq.89) w='RESTART'
            if(itoken.eq.90) w='CONFIGURATION'
            if(itoken.eq.91) w='SCREEN'
            if(itoken.eq.92) w='CIF'
            if(itoken.eq.93) w='EQUIVALENCE'
            if(itoken.eq.94) w='COMPTON'
            if(itoken.eq.95) w='RHOZZP'
            if(itoken.eq.96) w='CGRID'
            if(itoken.eq.97) w='CORVAL'
            if(itoken.eq.98) w='SIGGK'
            if(itoken.eq.99) w='TEMP'
            if(itoken.eq.100) w='DENS'
            if(itoken.eq.101) w='RIXS'
            if(itoken.eq.102) w='RLPR'
            if(itoken.eq.103) w='ICOR'
            if(itoken.eq.104) w='HUBBARD'
            if(itoken.eq.105) w='CRPA'
            if(itoken.eq.106) w='FULLSPECTRUM'
            if(itoken.eq.107) w='SCXC'
            if(itoken.eq.108) w='HIGHZ'
            if(itoken.eq.109) w='SCFTH'

      elseif (flname(1:10).eq.'spring.inp') then
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     These tokens are for spring.inp (input for eq of motion method)
            if(itoken.eq.1) w='STRETCH'
            if(itoken.eq.2) w='ANGLE'
            if(itoken.eq.3) w='VDOS'
            if(itoken.eq.4) w='PRDOS'
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      endif

      itoken_reverse=w
      
      return
      end

