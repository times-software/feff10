!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: m_atomicpotio.f90,v $:
! $Revision: 1.16 $
! $Author: jorissen $
! $Date: 2010/12/14 00:22:37 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE AtomicPotIO
  USE ErrorMod
  USE IOMod
  IMPLICIT NONE
  
  CHARACTER(3),PRIVATE :: FileType
  !PARAMETER(FileType = 'PAD')
  PARAMETER(FileType = 'TXT')
  
  
  CONTAINS
    SUBROUTINE WriteAtomicPots(nph, iz, ihole, rho, dmag, rhoval, vcoul, dgc0, dpc0, dgc, dpc, &
           & adgc, adpc, erelax, emu, xnvmu, xnval, norb, eorb, drho, dvcoul, iphat,    &
           & rat, iatph, novr, iphovr, nnovr, rovr, nat, edens, edenvl, vclap, rnrm,    &
           & kappa, iorb, s02)
      
      ! Scalar data
      INTEGER,INTENT(IN) :: nph, nat, ihole
      DOUBLE PRECISION, INTENT(IN) :: emu, erelax, s02

      ! 1D data
      INTEGER,INTENT(IN) :: iz(:), norb(:), iphat(:), iatph(:), novr(:)
      DOUBLE PRECISION,INTENT(IN) :: dgc0(:), dpc0(:), drho(:), dvcoul(:), rnrm(:)

      ! 2D data
      INTEGER,INTENT(IN) :: iphovr(:,:), nnovr(:,:), kappa(:,:), iorb(:,:)
      DOUBLE PRECISION,INTENT(IN) :: rho(:,:), dmag(:,:), rhoval(:,:), vcoul(:,:), &
           & xnvmu(:,:), xnval(:,:), eorb(:,:), rat(:,:), rovr(:,:), edens(:,:),   &
           & edenvl(:,:), vclap(:,:)

      ! 3D data
      DOUBLE PRECISION,INTENT(IN) :: dgc(:,:,:), dpc(:,:,:), adgc(:,:,:), adpc(:,:,:)


         

      CHARACTER*80 Headers(20)
      CHARACTER*20 ColumnLabels(20)

      INTEGER iph

      ! Initialize stuff.
      Headers(:) = ' '
      ColumnLabels(:) = ' '

      ! Define some headers.
      Headers(1) = 'This file contains information about the free atom potentials.'
      Headers(2) = 'nat    - number of atoms.'
      Headers(3) = 'nph    - number of unique potentials'      
      Headers(3) = 'ihole  - hole index'      
      Headers(4) = 'erelax - relaxation energy for each occupied orbital.'
      Headers(5) = 'emu - edge energy for each occupied orbital.'
      Headers(6) = 's02 - many body amplitude reduction.'
      ColumnLabels(1) = 'nph'
      ColumnLabels(2) = 'nat'
      ColumnLabels(3) = 'ihole'
      ColumnLabels(4) = 'erelax'
      ColumnLabels(5) = 'emu'
      ColumnLabels(6) = 's02'
      ! Write headers along with the first line and columnlabels.
      ! Write scalar data.
      CALL WriteData('apot.bin', Int1 = nph, Int2 = nat, Int3 = ihole, Double4 = erelax, Double5 = emu, &
           & Double6 = s02, Headers = Headers, ColumnLabels = ColumnLabels, FileType = FileType, &
		   ForceNewSection = .TRUE. ) !KJ 7-09 added ForceNewSection bc otherwise doesn't work on my Windows pc
      Headers(:) = ' '
      ColumnLabels(:) = ' '

      ! Write iz, iatph, novr, rnrm.
      Headers(1) = 'iz(0:nphx)    - atomic number for each unique potential'
      Headers(2) = 'iatph(0:nphx) - given unique pot, which atom is model?'
      Headers(3) = 'novr(0:nphx)  - number of overlap shells for each unique pot'
      Headers(4) = 'rnrm(0:nphx)  - norman radius for each unique potential.'
      ColumnLabels(1) = 'iz'
      ColumnLabels(2) = 'iatph'
      ColumnLabels(3) = 'novr'
      ColumnLabels(4) = 'rnrm'
      CALL WriteArrayData('apot.bin', Int1 = iz, Int2 = iatph, Int3 = novr, Double4 = rnrm, &
           & Headers = Headers, ColumnLabels = ColumnLabels, FileType = FileType,           &
           & ForceNewSection = .TRUE.)
      Headers(:) = ' '
      ColumnLabels(:) = ' '

      ! Write norb
      Headers(1) = 'norb(0:nphx+1) - number of occupied orbitals for each unique potential'
      CALL WriteArrayData('apot.bin', Int1 = norb, Headers = Headers, &
           & FileType = FileType, ForceNewSection = .TRUE.)

      ! Write iphat
      Headers(1) = 'iphat(natx) - given specific atom, which unique pot?'
      CALL WriteArrayData('apot.bin', Int1 = iphat, Headers = Headers, &
           & FileType = FileType, ForceNewSection = .TRUE.)

      ! Write dgc0, dpc0, drho, and dvcoul
      Headers(1) = 'dgc0   - upper component of core hole orbital'
      Headers(2) = 'dpc0   - lower component of core hole orbital'
      Headers(3) = 'drho   - core hole density.'
      Headers(4) = 'dvcoul - core hole coulomb potential.'
      ColumnLabels(1) = 'dgc0'
      ColumnLabels(2) = 'dpc0'
      ColumnLabels(3) = 'drho'
      ColumnLabels(4) = 'dvcoul'
      CALL WriteArrayData('apot.bin', Double1 = dgc0, Double2 = dpc0, Double3 = drho, &
           & Double4 = dvcoul, Headers = Headers, ColumnLabels = ColumnLabels,          &
           & FileType = FileType, ForceNewSection = .TRUE.)
      Headers(:) = ' '
      ColumnLabels(:) = ' '

      ! Write iphovr.
      Headers(1) = 'iphovr(novrx,0:nphx) - unique pot for each overlap shell'
      CALL Write2D('apot.bin', iphovr, Headers = Headers, FileType = FileType)

      ! Write nnovr.
      Headers(1) = 'nnovr(novrx,0:nphx) - number of atoms in overlap shell'
      CALL Write2D('apot.bin', nnovr, Headers = Headers, FileType = FileType)

      ! Write rho.
      Headers(1) = 'rho(r,0:nphx+1) - atomic density for each unique potential'
      Headers(2) = '                  nph+1 holds final state potential for absorber'
      CALL Write2D('apot.bin', rho, Headers = Headers, FileType = FileType)
      Headers(:) = ' '

      ! Write dmag.
      Headers(1) = 'dmag(r,nph+1) - ?'
      CALL Write2D('apot.bin', dmag, Headers = Headers, FileType = FileType)

      ! Write rhoval.
      Headers(1) = 'rhoval(r,nph+1) - atomic valence density for each unique potential.'
      CALL Write2D('apot.bin', rhoval, Headers = Headers, FileType = FileType)

      ! Write vcoul.
      Headers(1) = 'vcoul(r,nph) - coulomb potential for each unique potential.'
      CALL Write2D('apot.bin', vcoul, Headers = Headers, FileType = FileType)

      ! Write xnvmu.
      Headers(1) = 'xnvmu(0:lx,0:nphx+1) - number of valence electron within norman sphere for each channel (?)'
      CALL Write2D('apot.bin', xnvmu, Headers = Headers, FileType = FileType)

      ! Write xnval.
      Headers(1) = 'xnval(30,nph) - occupation of each orbital.'
      CALL Write2D('apot.bin', xnval, Headers = Headers, FileType = FileType)

      ! Write eorb.
      Headers(1) = 'eorb(norb,iph)'
      CALL Write2D('apot.bin', eorb, Headers = Headers, FileType = FileType)

      ! Write rat.
      Headers(1) = 'rat(3,nat) - cartesian coordinates for each atom.'
      CALL Write2D('apot.bin', rat, Headers = Headers, FileType = FileType)

      ! Write rovr
      Headers(1) = 'rovr(novrx,0:nphx) - r for overlap shell'
      CALL Write2D('apot.bin', rovr, Headers = Headers, FileType = FileType)

      ! Write edens.
      Headers(1) = 'edens(r,0:nphx) - overlapped density for each unique potential.'
      CALL Write2D('apot.bin', edens, Headers = Headers, FileType = FileType)

      ! Write edenvl
      Headers(1) = 'edenvl(r,0:nphx) - overlapped density of valence electrons?'
      CALL Write2D('apot.bin', edenvl, Headers = Headers, FileType = FileType)


      ! Write vclap.
      Headers(1) = 'vclap(r,0:nphx) - overlapped coulomb potential.'
      CALL Write2D('apot.bin', vclap, Headers = Headers, FileType = FileType)

      ! Write kappa.
      Headers(1) = 'kappa(norb,0:nph) - quntum number kappa for each orbital and potential.'
      CALL Write2D('apot.bin', kappa, Headers = Headers, FileType = FileType)

      ! Write iorb.
      Headers(1) = 'iorb(-4:3,0:nphx+1) - last occupied orbital of a particular kappa or zero if none.'
      CALL Write2D('apot.bin', iorb, Headers = Headers, FileType = FileType)
      
      ! Write dgc.
      Headers(1) = 'dgc(r,30,nph) - upper component of each obital for each unique potential.'
      CALL WriteData('apot.bin',Headers = Headers)
      DO iph = 1, SIZE(dgc,3)
         IF(iph.eq.1) THEN
             Headers(1) = 'dgc(r,30,nph) - upper component of each obital for each unique potential.'
             Headers(2) = 'dgc(r,norb,' // ACHAR(iph+48) // ')'
          ELSEIF(iph.lt.10) THEN
             Headers(1) = 'dgc(r,norb,' // ACHAR(iph+48) // ')'
		  ELSEIF(iph.lt.100) THEN
             Headers(1) = 'dgc(r,norb,' // ACHAR((iph/10)+48) // ACHAR(mod(iph,10)+48) // ')'
		  ELSE
		     stop 'ERROR iph too large in WriteAtomicPots'
          END IF
         CALL Write2D('apot.bin', dgc(:,:,iph), Headers = Headers, FileType = FileType)
      END DO

      ! Write dpc.
      Headers(1) = 'dpc(r,30,nph) - lower component of each obital and unique potential.'
      CALL WriteData('apot.bin',Headers = Headers)
      DO iph = 1, SIZE(dpc,3)
         IF(iph.eq.1) THEN
            Headers(1) = 'dpc(r,30,nph) - lower component of each obital and unique potential.'
             Headers(2) = 'dpc(r,norb,' // ACHAR(iph+48) // ')'
          ELSEIF(iph.lt.10) THEN
             Headers(1) = 'dpc(r,norb,' // ACHAR(iph+48) // ')'
		  ELSEIF(iph.lt.100) THEN
             Headers(1) = 'dpc(r,norb,' // ACHAR((iph/10)+48) // ACHAR(mod(iph,10)+48) // ')'
		  ELSE
		     stop 'ERROR iph too large in WriteAtomicPots'
          END IF
         CALL Write2D('apot.bin', dpc(:,:,iph), Headers = Headers, FileType = FileType)
      END DO

      ! Write adgc.
      DO iph = 1, SIZE(adgc,3)
         IF(iph.eq.1) THEN
            Headers(1) = 'adgc(r,30,nph) - upper development coeficients for each obital and unique potential.'            
            WRITE(Headers(2),'(A,I10,A)') 'adgc(r,norb,', iph, ')'
         ELSE
            WRITE(Headers(1),'(A,I10,A)') 'adgc(r,norb,', iph, ')'
         END IF
         CALL Write2D('apot.bin', adgc(:,:,iph), Headers = Headers, FileType = FileType)
      END DO

      ! Write adpc.
      Headers(1) = 'adpc(r,30,nph) - lower development coeficients for each obital and unique potential.'
      CALL WriteData('apot.bin',Headers = Headers)
      DO iph = 1, SIZE(adpc,3)
         IF(iph.eq.1) THEN
             Headers(1) = 'adpc(r,30,nph) - lower development coeficients for each obital and unique potential.'
             WRITE(Headers(2),'(A,I10,A)') 'adpc(r,norb,', iph, ')'
          ELSE
             WRITE(Headers(1),'(A,I10,A)') 'adpc(r,norb,', iph, ')'
          END IF
         CALL Write2D('apot.bin', adpc(:,:,iph), Headers = Headers, FileType = FileType)
      END DO

      ! Close the file so that we can read it later.
      CALL CloseFl('apot.bin')
    END SUBROUTINE WriteAtomicPots







    SUBROUTINE ReadAtomicPots(nph, iz, ihole, rho, dmag, rhoval, vcoul, dgc0, dpc0, dgc, dpc, &
         & adgc, adpc, erelax, emu, xnvmu, xnval, norb, eorb, drho, dvcoul, iphat,     &
         & rat, iatph, novr, iphovr, nnovr, rovr, nat, edens, edenvl, vclap, rnrm,     &
         & kappa, iorb, s02)

!  CALL ReadAtomicPots(nph = nph, iz = iz(0:nph), ihole = ihole, dmag = dmag(:,0:nph), &
!            dgc0 = dgc0, dpc0 = dpc0, dgc = dgc(:,:,0:nph+1), dpc = dpc(:,:,0:nph+1),  &
!            adgc = adgc(:,:,0:nph+1), adpc = adpc(:,:,0:nph+1), erelax = erelax,       &
!            emu = emu_apot, xnval = xnval(:,0:nph+1), eorb = eorbTmp, rnrm = rnrm(0:nph))

      ! Scalar data
      INTEGER,INTENT(OUT),OPTIONAL :: nph, nat, ihole
      DOUBLE PRECISION,INTENT(OUT),OPTIONAL :: emu, erelax, s02

      ! 1D data
      INTEGER,INTENT(OUT),OPTIONAL :: iz(:), norb(:), iphat(:), iatph(:), novr(:)
      DOUBLE PRECISION,INTENT(OUT),OPTIONAL :: dgc0(:), dpc0(:), drho(:), dvcoul(:), rnrm(:)

      ! 2D data
      INTEGER,INTENT(OUT),OPTIONAL :: iphovr(:,:), nnovr(:,:), kappa(:,:), iorb(:,:)
      DOUBLE PRECISION,INTENT(OUT),OPTIONAL :: rho(:,:), dmag(:,:), rhoval(:,:), vcoul(:,:), &
           & xnvmu(:,:), xnval(:,:), eorb(:,:), rat(:,:), rovr(:,:), edens(:,:),   &
           & edenvl(:,:), vclap(:,:)

      ! 3D data
      DOUBLE PRECISION,INTENT(OUT),OPTIONAL :: dgc(:,:,:), dpc(:,:,:), adgc(:,:,:), adpc(:,:,:)

      INTEGER n1, n2, nSect, iError, iph

      ! Temp variables.
      INTEGER iTmp1, iTmp2, iTmp3
      INTEGER,ALLOCATABLE :: iTmpArray1(:), iTmpArray2(:), iTmpArray3(:)
      DOUBLE PRECISION DTmp1, DTmp2, DTmp3
      DOUBLE PRECISION,ALLOCATABLE :: DTmpArray1(:), DTmpArray2(:), DTmpArray3(:), DTmpArray4(:)

      NSect = 1
      ! Read nph, nat
      ! PRINT*, 'Read nph, nat'
      IF(PRESENT(nph).or.PRESENT(nat).or.PRESENT(ihole).or.PRESENT(emu).or.PRESENT(erelax)) THEN
         CALL ReadData('apot.bin', Int1 = iTmp1, Int2 = iTmp2, Int3 = iTmp3, Double4 = DTmp1, & 
              & Double5 = DTmp2 , Double6 = DTmp3, FileType = FileType, SectionNumber = NSect)
         IF(PRESENT(nph)) nph = iTmp1
         IF(PRESENT(nat)) nat = iTmp2
         IF(PRESENT(ihole)) ihole = iTmp3
         IF(PRESENT(erelax)) erelax = DTmp1
         IF(PRESENT(emu)) emu = DTmp2
         IF(PRESENT(s02)) s02 = DTmp3
      END IF
      !write(*,*) NSect,' done'

      NSect = NSect + 1
      IF(PRESENT(iz).or.PRESENT(iatph).or.PRESENT(novr).or.PRESENT(rnrm)) THEN
         ! Read iz, iatph, novr, rnrm.
         IF(PRESENT(iz)) THEN
            n1 = SIZE(iz)
         ELSEIF(PRESENT(iatph)) THEN
            n1 = SIZE(iatph)
         ELSEIF(PRESENT(novr)) THEN
            n1 = SIZE(novr)
         ELSEIF(PRESENT(rnrm)) THEN
            n1 = SIZE(rnrm)
         END IF
         !write(*,*) 'size is ',n1,size(rnrm)

         ALLOCATE(iTmpArray1(n1),STAT = iError)
         CALL CheckAllocation(iError,'ERROR: Cannot allocate space for array iTmpArray1 in subroutine ReadAtomicPots.')
         ALLOCATE(iTmpArray2(n1),STAT = iError)
         CALL CheckAllocation(iError,'ERROR: Cannot allocate space for array iTmpArray2 in subroutine ReadAtomicPots.')
         ALLOCATE(iTmpArray3(n1),STAT = iError)
         CALL CheckAllocation(iError,'ERROR: Cannot allocate space for array iTmpArray3 in subroutine ReadAtomicPots.')
         ALLOCATE(DTmpArray1(n1),STAT = iError)
         CALL CheckAllocation(iError,'ERROR: Cannot allocate space for array DTmpArray1 in subroutine ReadAtomicPots.')
         CALL ReadArrayData('apot.bin', Int1 = iTmpArray1, Int2 = iTmpArray2, Int3 = iTmpArray3, Double4 = DTmpArray1, &
              & FileType = FileType, SectionNumber = NSect)
         IF(PRESENT(iz)) iz(:)    = iTmpArray1(:)
         IF(PRESENT(iatph)) iatph(:) = iTmpArray2(:)
         IF(PRESENT(novr)) novr(:)  = iTmpArray3(:)
         IF(PRESENT(rnrm)) rnrm(:)  = DTmpArray1(:)
         !write(*,*) 'sizes',size(DTmpArray1),size(rnrm)
         !write(*,*) 'content',DTmpArray1, 'and ',rnrm

         DEALLOCATE(iTmpArray1, iTmpArray2, iTmpArray3, DTmpArray1)
      END IF

      NSect = NSect + 1
      IF(PRESENT(norb)) THEN
         ! Read norb
         CALL ReadArrayData('apot.bin', Int1 = norb,  &
              & FileType = FileType, SectionNumber = NSect)
      END IF

      NSect = NSect + 1
      IF(PRESENT(iphat)) THEN
         ! Read iphat
         CALL ReadArrayData('apot.bin', Int1 = iphat,  &
              & FileType = FileType, SectionNumber = NSect)
      END IF
      !write(*,*) NSect,' done'

      NSect = NSect + 1
      IF(PRESENT(dgc0).or.PRESENT(dgc0).or.PRESENT(drho).or.PRESENT(dvcoul)) THEN
         ! Read dgc0, dpc0, drho, and dvcoul
         IF(PRESENT(dgc0)) THEN
            n1 = SIZE(dgc0)
         ELSEIF(PRESENT(dpc0)) THEN
            n1 = SIZE(dpc0)
         ELSEIF(PRESENT(drho)) THEN
            n1 = SIZE(drho)
         ELSEIF(PRESENT(dvcoul)) THEN
            n1 = SIZE(dvcoul)
         END IF

         ALLOCATE(DTmpArray1(n1),STAT = iError)
         CALL CheckAllocation(iError,'ERROR: Cannot allocate space for array iTmpArray1 in subroutine ReadAtomicPots.')
         ALLOCATE(DTmpArray2(n1),STAT = iError)
         CALL CheckAllocation(iError,'ERROR: Cannot allocate space for array iTmpArray2 in subroutine ReadAtomicPots.')
         ALLOCATE(DTmpArray3(n1),STAT = iError)
         CALL CheckAllocation(iError,'ERROR: Cannot allocate space for array iTmpArray3 in subroutine ReadAtomicPots.')
         ALLOCATE(DTmpArray4(n1),STAT = iError)
         CALL CheckAllocation(iError,'ERROR: Cannot allocate space for array DTmpArray1 in subroutine ReadAtomicPots.')

         CALL ReadArrayData('apot.bin', Double1 = DTmpArray1, Double2 = DTmpArray2, &
              & Double3 = DTmpArray3, Double4 = DTmpArray4, FileType = FileType, SectionNumber = NSect)
         IF(PRESENT(dgc0))   dgc0    = DTmpArray1
         IF(PRESENT(dpc0))   dpc0    = DTmpArray2
         IF(PRESENT(drho))   drho    = DTmpArray3
         IF(PRESENT(dvcoul)) dvcoul  = DTmpArray4

         DEALLOCATE(DTmpArray1, DTmpArray2, DTmpArray3, DTmpArray4)
      END IF

      NSect = NSect + 1
      IF(PRESENT(iphovr)) THEN
         ! Read iphovr.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', iphovr,  n1, n2, FileType = FileType)
      END IF

      NSect = NSect + 1
      IF(PRESENT(nnovr)) THEN
         ! Read nnovr.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', nnovr,  n1, n2, FileType = FileType)
      END IF

      NSect = NSect + 1
      IF(PRESENT(rho)) THEN
         ! Read rho.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', rho,  n1, n2, FileType = FileType)
      END IF
      !write(*,*) NSect,' done'

      NSect = NSect + 1
      IF(PRESENT(dmag)) THEN
         ! Read dmag.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', dmag,  n1, n2, FileType = FileType)
      END IF

      NSect = NSect + 1
      IF(PRESENT(rhoval)) THEN
         ! Read rhoval.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', rhoval, n1, n2, FileType = FileType)
      END IF

      NSect = NSect + 1
      IF(PRESENT(vcoul)) THEN
         ! Read vcoul.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', vcoul, n1, n2, FileType = FileType)
      END IF

      NSect = NSect + 1
      IF(PRESENT(xnvmu)) THEN
         ! Read xnvmu.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', xnvmu, n1, n2, FileType = FileType)
      END IF

      NSect = NSect + 1
      IF(PRESENT(xnval)) THEN
         ! Read xnval.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', xnval, n1, n2, FileType = FileType)
      END IF

      NSect = NSect + 1
      IF(PRESENT(eorb)) THEN
         ! Read eorb.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', eorb, n1, n2, FileType = FileType)
      END IF
      !write(*,*) NSect,' done'

      NSect = NSect + 1
      ! PRINT*, 'Read rat'
      IF(PRESENT(rat)) THEN
         ! Read rat.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', rat, n1, n2, FileType = FileType)
      END IF

      NSect = NSect + 1
      IF(PRESENT(rovr)) THEN
         ! Read rovr
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', rovr, n1, n2, FileType = FileType)
      END IF

      NSect = NSect + 1
      IF(PRESENT(edens)) THEN
         ! Read edens.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', edens, n1, n2, FileType = FileType)
      END IF

      NSect = NSect + 1
      IF(PRESENT(edenvl)) THEN
         ! Read edenvl
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', edenvl, n1, n2, FileType = FileType)
      END IF

      NSect = NSect + 1
      IF(PRESENT(vclap)) THEN
         ! Read vclap.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', vclap, n1, n2, FileType = FileType)
      END IF
      !write(*,*) NSect,' done'

      NSect = NSect + 1
      IF(PRESENT(kappa)) THEN
         ! Read kappa.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', kappa, n1, n2, FileType = FileType)
      END IF

      NSect = NSect + 1
      IF(PRESENT(iorb)) THEN
         ! Read iorb.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         CALL Read2D('apot.bin', iorb, n1, n2, FileType = FileType)
      END IF
      
      NSect = NSect + 1
      IF(PRESENT(dgc)) THEN
         ! Read dgc.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         DO iph = 1, SIZE(dgc,3)
            CALL Read2D('apot.bin', dgc(:,:,iph), n1, n2, FileType = FileType)
         END DO
      END IF

      NSect = NSect + nph + 2
      IF(PRESENT(dpc)) THEN
         ! Read dpc.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         DO iph = 1, SIZE(dpc,3)
            CALL Read2D('apot.bin', dpc(:,:,iph), n1, n2, FileType = FileType)
         END DO
      END IF

      NSect = NSect + nph + 2
      IF(PRESENT(adgc)) THEN
         ! Read adgc.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         DO iph = 1, SIZE(adgc,3)
            CALL Read2D('apot.bin', adgc(:,:,iph), n1, n2, FileType = FileType)
         END DO
      END IF
      !write(*,*) NSect,' done'

      NSect = NSect + nph + 2
      IF(PRESENT(adpc)) THEN
         ! Read adpc.
         CALL ReadToSection('apot.bin',SectionNumber = NSect)
         DO iph = 1, SIZE(adpc,3)         
            CALL Read2D('apot.bin', adpc(:,:,iph), n1, n2, FileType = FileType)
         END DO
      END IF

      CALL CloseFl('apot.bin')
      ! PRINT*, 'Done reading information from atomic calculation.' 
    END SUBROUTINE ReadAtomicPots
    
        ! Hack for Fer's mt potentials
    SUBROUTINE WriteExternalPot(vtot, vint, edens, rhoint, rat, xmu, imt, rmt)
      DOUBLE PRECISION vtot(:,:), vint, edens(:,:), rhoint, rat(:,:), xmu, rmt(:), RadialGrid(251), xNElOut
      DOUBLE PRECISION,ALLOCATABLE :: ratTranspose(:,:)
      INTEGER nat, NRPts, i1, i2, imt(:), ncoord
      INTEGER,ALLOCATABLE :: iat(:)
      CHARACTER(30) FileName

      ncoord = 3
      FileName = 'extpot.aip'
      nat = 2
      ! Read number of atoms
      CALL WriteData(FileName, Int1 = nat, ForceNewSection = .TRUE.)

      ! Read atomic numbers
      ALLOCATE(iat(nat))
      iat(1) = 6
      iat(2) = 8
      CALL WriteArrayData(FileName, Int1 = iat, ForceNewSection = .TRUE.)
      
      ! Read coordinates
      ALLOCATE(ratTranspose(nat,3))
      DO i1 = 1, nat
         DO i2 = 1, 3
            ratTranspose(i1,i2) = rat(i2, i1)
         END DO
      END DO
      CALL Write2D(FileName, ratTranspose)

      ! Read number of radial points
      NRPts = 210
      CALL WriteData(FileName, Int1 = NRPts, ForceNewSection = .TRUE.)
      
      ! Read Interstitial potential
      CALL WriteData(FileName, Double1 = vint, ForceNewSection = .TRUE.)

      ! Read Fermi energy
      CALL WriteData(FileName, Double1 = xmu, ForceNewSection = .TRUE.)

      ! Read number of electrons outside muffin tins.
      xNElOut = 1
      CALL WriteData(FileName, Double1 = xNElOut, ForceNewSection = .TRUE.)

      ! Read muffin tin radii and jri.
      CALL WriteArrayData(FileName, Int1 = imt, Double2 = rmt, ForceNewSection = .TRUE.)
      
      ! Read radial grid, vtot, and edens
      DO i1 = 1, nat
         DO i2 = 1, NRPts
            CALL WriteData(FileName, Double1 = RadialGrid(i2), Double2 = vtot(i2,i1), Double3 = edens(i2,i1))
         END DO
      END DO

      CALL CloseFl(FileName)
    END SUBROUTINE WriteExternalPot


    ! Hack for Fer's mt potentials
    SUBROUTINE ReadExternalPot(vtot, vint, edens, rhoint, rat, xmu, imt, rmt)
      USE Mtdp
      USE IOFiles
      
      DOUBLE PRECISION vtot(:,:), vint, edens(:,:), rhoint, rat(:,:), xmu, rmt(:), RadialGrid(251), vTmp, xNElOut
      DOUBLE PRECISION :: rmt2(SIZE(rmt)), vtot2(SIZE(vtot,1),SIZE(vtot,2)), edens2(SIZE(edens,1),SIZE(edens,2)) 
      INTEGER nat, NRPts, i1, i2, imt(:), ncoord, iError, iSort(200), iunit
      INTEGER :: iz(SIZE(imt)), imt2(SIZE(imt))
      CHARACTER(30) PotFile, SortFile, mtdpFile
      TYPE(Mtdp_Data_Type) :: Mtdp_Data
      ncoord = 3
      PotFile  = 'extpot.aip'
      SortFile = 'sort.aip'
      mtdpFile = 'GeCl4.04.dft.mtdp'

      ! Read mtdp file.
      CALL OpenFl(mtdpFile,FileStatus='OLD',FileAction='READ')
      CALL GetIOFileInfo(mtdpFile, UnitNumber = iunit)
      CALL Read_Mtdp(iunit,Mtdp_Data)

      ! Number of points in the radial grid
      NRPts = Mtdp_Data%nR
      
      ! Number of atoms
      nat = Mtdp_Data%nAt
      
      ! For now atomic numbers are not needed.

      ! Atomic coordinates - ignore for now
!      rat(:,:) = Mtdp_Data%At_XYZ

      ! Muffin-Tin Radii
      rmt2(:nat) = Mtdp_Data%At_R(:)

      ! Index of the radii on the radial grid
      imt2(:nat) = Mtdp_Data%At_iR(:)

      ! Electron density inside each Muffin-Tin
      edens2(:,:nat) = Mtdp_Data%At_Den(:,:)

      ! Potential (H+XC) inside each Muffin-Tin
      vtot2(:,:nat) = Mtdp_Data%At_Pot(:,:)

      ! Empty spheres

      ! muffin tin radii
      rmt2(nat+1:Mtdp_Data%nESph+nat) = Mtdp_Data%ESph_R(:)
      
      ! muffin tin index
      imt2(nat+1:Mtdp_Data%nESph+nat) = Mtdp_Data%ESph_iR(:)

      ! Electron density inside each empty sphere
      edens2(:,nat+1:Mtdp_Data%nESph+nat) = Mtdp_Data%ESph_Den(:,:)

      ! Potential (H+XC) inside each empty sphere
      vtot2(:,nat+1:Mtdp_Data%nESph+nat) = Mtdp_Data%ESph_Pot(:,:)

      !PRINT*, 'Old vint = ', vint
      ! Interstitial potential
      vint = Mtdp_Data%V_Int
      !PRINT*, 'Enter vint: '
      !READ*, vint    
    
      ! HOMO energy
      xmu = (Mtdp_Data%V_HOMO + Mtdp_Data%V_LUMO)/2.d0
!      xmu = Mtdp_Data%V_HOMO
      nat = nat + Mtdp_Data%nESph
      ! Read isort out of sort.aip. This tells which potentials go with which atoms in feff.
      iSort(:) = -1
      ! PRINT*, 'Reading iSort'
      CALL ReadArrayData(SortFile, Int1 = iSort)
      ! iSort(i1) is the unique potential given the i1th potential defined in the extpot.aip file.
      ! iSort can be zero, but arrays here are defined from 1, so we need to shift.
      DO i1 = 1, SIZE(iSort)
         iSort(i1) = iSort(i1) + 1
      END DO

      ! Fill vtot, edens, rmt, and imt with the proper values based on iSort
      DO i1 = 1, SIZE(iSort)
         IF(iSort(i1).gt.0) THEN
            ! Fill rmt and imt
            IF(iSort(i1).gt.nat) CALL Error('ERROR: Number of potentials defined '// &
              & 'in sort.aip is greater than number defined in extpot.aip.')
            rmt(iSort(i1)) = rmt2(i1)
            imt(iSort(i1)) = imt2(i1)
            
            ! Fill vtot and edens
            DO i2 = 1, NRPts
               vtot(i2,iSort(i1)) = vtot2(i2,i1)
               edens(i2,iSort(i1)) = edens2(i2,i1)
            END DO
         END IF         
      END DO
      
      ! If vint is too close to zero, reset it.
!      IF(vint.gt.-0.1d0) vint = -0.1d0
      DO i1 = 1, nat
         DO i2 = NRPts + 1, 251
            vtot(i2,i1) = vint
            edens(i2,i1) = rhoint
         END DO
      END DO

!      DO i1 = 1, 251
!         vtot(i2,nat+1) = vtot(NRPts,1)
!         edens(i2,nat+1) = edens(NRPts,1)
!      END DO

      CALL CloseFl(SortFile)
      CALL CloseFl(mtdpFile)
      ! PRINT*, 'Done reading external potential.'
    END SUBROUTINE ReadExternalPot

!     SUBROUTINE SplitMtDP(PotFiles,mtdpFile)
!       USE Mtdp
!       USE IOFiles
!       ! This subroutine reads an mtdp file and splits it into nph pot files.
!       ! Input is the array of potential file names.
!       CHARACTER PotFiles(:)

!       INTEGER iunit, iat
!       CHARACTER AtomLabel

!       CALL OpenFl(mtdpFile)
!       CALL GetIOFileInfo(UnitNumber = iunit)
!       CALL Read_Mtdp(iunit,Mtdp_Data)

!       DO iat = 1, Mtdp_Data%nAt
!          ! Number of points in the radial grid
!          Mtdp_Data%nR

!          CALL WriteData



    SUBROUTINE ReadExternalPotWien2k   ! (vtot, vint, edens, rhoint, rat, xmu, imt, rmt)

    END SUBROUTINE ReadExternalPotWien2k


      
END MODULE AtomicPotIO
