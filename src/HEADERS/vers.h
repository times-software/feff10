      character*12 vfeff
      CHARACTER(300) revision
!                       123456789012  
!      PARAMETER(vfeff='FEFF 9.7.1 ')
      PARAMETER (vfeff='FEFF 10.0.0 ')
      PARAMETER (revision='Revision 1')

      ! 9.00 prepared by Yoshi and Josh from feff86
      ! 9.01 JP starts dynamic allocation, KJ merges and cleans up
      ! 9.02 KJ integrates feff8q, revision of rdinp and modules, etc.
      ! 9.03 JK a few bug fixes.
      ! 9.05 JK: first official release of feff9
      ! 9.1  KJ decides it's time to move on to the next digit.
      ! 9.5  KJ: big clean-up, integration with JFEFF
      ! 9.6  KJ: First real public release since 9.1
      ! 9.7  KJ: Big clean-up of output; MKL integration; new modules RHORRP/density
      !10.0  KJ: New major version.  NRIXS, HUBBARD-U (d/f), loads of smaller updates and bugfixes.
      !9.9   KJ: Whoops back to "9" for political reasons
      !10.0  JK: Added finite temperature and RIXS that will work with corvus.
