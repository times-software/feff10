!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: pertab.f90,v $:
! $Revision: 1.2 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Periodic table of the elements
!     Written by Steven Zabinsky, Feb 1992.  Deo Soli Gloria

!     atwts(iz)  single precision fn, returns atomic weight
!     atwtd(iz)  double precision fn, returns atomic weight
!     atsym(iz)  character*2 fn, returns atomic symbol

      double precision function atwtd (iz)
      double precision weight
      common /atwtco/ weight(103)
      atwtd = weight(iz)
      return
      end

      real function atwts (iz)
      double precision weight
      common /atwtco/ weight(103)
      atwts = weight(iz)
      return
      end

      character*2 function atsym (iz)
      character*2 sym
      common /atsyco/ sym(103)
      atsym = sym(iz)
      return
      end

      block data prtbbd
!     PeRiodic TaBle Block Data

!     Atomic weights from inside front cover of Ashcroft and Mermin.

      double precision weight
      common /atwtco/ weight(103)

      character*2 sym
      common /atsyco/ sym(103)

      data weight /                                                     &
     &   1.0079, 4.0026, 6.941,  9.0122, 10.81,   12.01,                &
     &   14.007, 15.999, 18.998, 20.18,  22.9898, 24.305,               &
     &   26.982, 28.086, 30.974, 32.064, 35.453,  39.948,               &
     &   39.09,  40.08,  44.956, 47.90,  50.942,  52.00,                &
     &   54.938, 55.85,  58.93,  58.71,  63.55,   65.38,                &
     &   69.72,  72.59,  74.922, 78.96,  79.91,   83.80,                &
     &   85.47,  87.62,  88.91,  91.22,  92.91,   95.94,                &
     &   98.91,  101.07, 102.90, 106.40, 107.87,  112.40,               &
     &   114.82, 118.69, 121.75, 127.60, 126.90,  131.30,               &
     &   132.91, 137.34, 138.91, 140.12, 140.91,  144.24,               &
     &   145,    150.35, 151.96, 157.25, 158.92,  162.50,               &
     &   164.93, 167.26, 168.93, 173.04, 174.97,  178.49,               &
     &   180.95, 183.85, 186.2,  190.20, 192.22,  195.09,               &
     &   196.97, 200.59, 204.37, 207.19, 208.98,  210,                  &
     &   210,    222,    223,    226,    227,     232.04,               &
     &   231,    238.03, 237.05, 244,    243,     247,                  &
     &   247,    251,    254,    257,    256,     254,                  &
     &   257/

      data sym /  'H', 'He','Li','Be','B', 'C', 'N', 'O', 'F', 'Ne',    &
     &            'Na','Mg','Al','Si','P', 'S', 'Cl','Ar','K', 'Ca',    &
     &            'Sc','Ti','V', 'Cr','Mn','Fe','Co','Ni','Cu','Zn',    &
     &            'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y', 'Zr',    &
     &            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',    &
     &            'Sb','Te','I', 'Xe','Cs','Ba','La','Ce','Pr','Nd',    &
     &            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',    &
     &            'Lu','Hf','Ta','W', 'Te','Os','Ir','Pt','Au','Hg',    &
     &            'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',    &
     &            'Pa','U', 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm',    &
     &            'Md','No','Lw'/

      end
