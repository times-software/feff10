      subroutine rdop (cntrl, emin, emax, iepts, ncomps, cmpnm, nz,     &
     &                 numden, numedg, alledg, conv, ldos, finest,      &
     &                 drude,tau, ndrude,valence,eels)

      !reads core.inp
      use constants
      use dimsmod, only: nheadx,nex,nphx=>nphu
      implicit none
      include 'HEADERS/params.h'
      character*2 alledg(2*mxedgs,nphx),ename,names(2*mxedgs),          &
     &         flnm*20, l1*75, l2*75, l3*75,flnm2*20,cmd*80
      character*3 cmnt,cmpnm(nphx),key
      integer ie, j, iepts, numedg(nphx),iv,ic,i,                       &
     &        cntrl(1:6),ios,                                           &
     &        ncomps,iunit,ierr, nedges, ii, ihole
      real*8 emin, emax, emaxi, myemin, myemax, myestp, numden(nphx), tau
      real*8 val, emu,lowedg,hiedg,ndrude
      double precision emudp
      logical, external :: isedge
      logical drude,  rdflag, conv(2*mxedgs,nphx),convt(2*mxedgs)
      logical finest (2*mxedgs,nphx),ldos(nphx), valence
!     input and output strings
      character*512 slog, str
      character*512 line, keywrd
      integer mwords, nwords,ilen,k,l,m, nz(nphx),num
      parameter (mwords = 20)
      character*512 words(mwords),wrdsor(mwords)
      integer iz
      logical eels !KJ 6/10

!*************************************     INITIALIZE TO DEFAULT VALUES
      ncomps=0
      !do everything if no control card
      cntrl(:)=1
      !do not convolute or compute fine structure unless user says to.
      conv(:,:)=.false.
      finest(:,:)=.false.
      !set energy grid params negative as a flag to set them based on
      !which edges will be in the spectrum if the user doesn't set the
      !grid.
      iepts=-1
      emax=-1.0
      emin=-2.0
      drude=.false. !no drude term 
      valence=.false. !do not include xmu.val by default
      tau=0.0 
      ndrude=0.0
      eels=.false. !KJ 7/2011 bugfix

!*************************************     READ OPTIONS FROM CORE.INP
      ii=-1 !init flag for getcom. -1 means open file.
      line=infile !put filename to open in line; first call to getcom opens.
      rdflag=.true. 
      !rdflag is set to false if str already holds a line from the input file that needs to be processed.
! Modivied by FDV:
! Avoid warning in gfortran.
! This is pretty crappy Fortran programming, by the way, whoever is
! responsible, shame on you.
! One little caveat, now we have an infinite loop, but the exit condition
! should be fullfilled at end of file.
!     do l=1,1e6
      do
        if (rdflag) then !get new line if str does not already hold one.
          call getcom(ii,line) !read one line of input file, returns it in line
          call fixstr (line,str,ilen,words,wrdsor,mwords,nwords)
          !takes line and returns a fixed-up version in str; the individual
          !words are returned in words (lower case) and wrdsor (original case);
          !nwords is the number of words, mwords is max # of words.
        endif
        if (str.eq.'getcom_nofile') then !test for input file non-existence
          call wlog ('Input file does not exist.')
          call wlog ('Filename: '//infile)
          call wlog ('Stopping.')
          stop
        endif 
        !pick out first 3 letters of first word of line to determine which card has been read. 
        keywrd=words(1)
        key=keywrd(1:3)

        !read information appropriate to each card
        
        if(key.eq.'con') then !CONTROL card
          do i=1,max(nwords,6)
            keywrd=words(i+1)
            call str2in(keywrd,j,ierr)
            if (ierr.eq.0) then
              cntrl(i)=j
            else
              call wlog('Error reading CONTROL card. Check syntax')
              call wlog ('Stopping.')
              stop
            endif
          enddo
          rdflag=.true.

!       EELS CARD BY KEVIN JORISSEN 6/2010
        elseif (key.eq.'eel' .or. key.eq.'EEL') then
           eels=.true.

        elseif (key.eq.'val') then 
           !VALENCE card: set flag to include xmu.val in edge list
           valence=.true.
        elseif (key.eq.'dru') then !DRUDE card: look for tau and ndrude
          drude=.true.
          keywrd=words(2)
          call str2re(keywrd,val,ierr) 
          if (ierr.eq.0) then
            tau=val
          else
            call wlog('Error reading DRUDE card. Check syntax (tau).')
            call wlog ('Stopping.')
            stop
          endif
          keywrd=words(3)
          ierr=0
          if (keywrd.ne.' ') then 
            call str2re(keywrd,val,ierr) 
            if (ierr.eq.0) then
              ndrude=val
            else
              call wlog                                                 &
     &              ('Error reading DRUDE card. Check syntax (ndrude).')
              call wlog ('Stopping.')
              stop
            endif
          endif
          rdflag=.true.

        elseif (key.eq.'det') then !DETAIL card: set finest to true
          do i=1,2*mxedgs
            do j=1,nphx
              finest(i,j)=.true.
            enddo
          enddo
          rdflag=.true.
          
          
        elseif (key.eq.'egr') then !EGRID card
          keywrd=words(2)
          call str2re(keywrd,val,ierr) 
          if (val.ge.0.and.ierr.eq.0) then
            emin=val/hart
          else
            slog=                                                       &
     &      'Error reading EGRID card: emin is negative or nonsensical!'
            call wlog (slog)
            write  (slog,fmt="('    emax: ',e10.3)") val
            call wlog ('Stopping.')
            stop
          endif
          keywrd=words(3)
          call str2re(keywrd,val,ierr) 
          if (val.gt.0.and.ierr.eq.0) then
            emax=val/hart
          else
            slog=                                                       &
     &      'Error reading EGRID card: emax is negative or nonsensical!'
            call wlog (slog)
            write  (slog,fmt="('    emax: ',e10.3)") val
            call wlog ('Stopping.')
            stop
          endif
          keywrd=words(4)
          call str2in(keywrd,j,ierr) 
          if (j.gt.0.and.ierr.eq.0) then
            iepts=j
          else
          endif
          rdflag=.true.

        elseif (key.eq.'com') then !COMPONENT card: get numden and read edge names
          rdflag=.true.
          if (nwords.lt.3) then 
            call wlog                                                   &
     &        ('Problem with COMPONENT card; bad line:')
            call wlog('Stopping.')
            stop
          endif
          ncomps=ncomps+1
          cmpnm(ncomps)=wrdsor(2)

          keywrd=words(3) !atomic number
          call str2in(keywrd,j,ierr)
          if (ierr.eq.0) then
            nz(ncomps)=j
          else
            call wlog                                                   &
     &        ('Problem with COMPONENT card: what is Z?')
            call wlog(line)
            call wlog('Stopping.')
          endif
          keywrd=words(4) !numden
          call str2re(keywrd,val,ierr) 
          if (ierr.eq.0.and.val.gt.0.0) then
            numden(ncomps)=val*bohr**3 !convert from angs**3 to bohr**3
          else
            call wlog                                                   &
     &      ('No number density listed; ' //                            &
     &       'number density will be estimated from feff.inp.')
             numden(ncomps)=-1.0
             rdflag=.true.
          endif
          if (ierr.eq.0) then !numden was read, look for EDGES card in words(5)
            keywrd=words(5)
          else !look for EDGES card in words(4).
            keywrd=words(4)
          endif
          key=keywrd(1:3)
          if (key.eq.'edg') then !read the edge names
            do i=1,2*mxedgs
              finest(i,ncomps)=.true.
            enddo
            numedg(ncomps)=0
            k=1
! Modivied by FDV:
! Avoid warning in gfortran.
! This is pretty crappy Fortran programming, by the way, whoever is
! responsible, shame on you.
! One little caveat, now we have an infinite loop, but the exit condition
! should be fullfilled at end of file.
!           do  m=1,1e6
            do
              call getcom(k,line) !read one line of input file
              call fixstr (line,str,ilen,words,wrdsor,mwords,nwords)
              if (str.eq.'getcom_end') exit 
              !EOF, we are done. This exit escapes inner do loop
              !(over edges in this component);
              !there is an identical test at the bottom of the  
              !outer do loop (over cards in the input file). If
              !the last card is COMPONENT, both will be used.
              if (str.eq.'getcom_error') then
                call wlog ('Error in input file '//infile)
                call wlog ('Stopping.')
                stop
              endif
              ename=words(1)
              call upper(ename)
              if(isedge(ename)) then 
                call stdnm(ename) !convert edge number to label if needed.
                !If we read a valid name, save it.
                numedg(ncomps)=numedg(ncomps)+1
                alledg(numedg(ncomps),ncomps)=ename
                !Determine if edge is core or valence, and note this in 
                ! the logical array conv.  
                keywrd=words(2)
                key=keywrd(1:3)
                if (key.eq.'con') then
                  !current edge is valence edge
                  conv(numedg(ncomps),ncomps)=.true.
                  write(slog,fmt="('Edge ',a2,' will be convolved.')")  &
     &              ename
                call wlog (slog)
                elseif (key.eq.'det') then
                  !current edge is core edge, check to see if we should
                  !compute fine structure (detail)
                  finest(numedg(ncomps),ncomps)=.true.
                elseif (key.eq.'bac') then
                  !user wants background only calculation
                  finest(numedg(ncomps),ncomps)=.false.
                endif

              else
                !Current line does not contain a valid edge name, so we 
                !assume list of edge names has ended for this EDGES card.
                rdflag=.false.
                exit
              endif
            enddo
          else !no EDGES card, so get edges for nz(ncomps) from getorb
            iz=nz(ncomps)
            call gtedgs(iz,num,names,convt)
            numedg(ncomps)=num
            do i=1,num
              alledg(i,ncomps)=names(i)
              !don't do convolutions for automatic run since straight
              !feff seems better
              conv(i,ncomps)=.false.
            enddo 
            if (key.eq.'det') then !read edges for fine structure calc.

              !set all edges for background calculation only; we will
              !turn on fine structure for each edge listed.
              do i=1,2*mxedgs
                finest(i,ncomps)=.false.
              enddo

              do  m=1,mxedgs
                call getcom(k,line) !read one line of input file
                call fixstr (line,str,ilen,words,wrdsor,mwords,nwords)
                if (str.eq.'getcom_end') exit 
                !EOF, we are done. This exit escapes inner do loop
                !(over edges in this DETAIL card);
                !there is an identical test at the bottom of the  
                !outer do loop (over cards in the input file). If
                !the last card is COMPONENT with DETAIL, both will be used.
                if (str.eq.'getcom_error') then
                  call wlog ('Error in input file '//infile)
                  call wlog ('Stopping.')
                  stop
                endif
                ename=words(1)
                call upper(ename)
                if(isedge(ename)) then 
                  call stdnm(ename) !convert edge number to label if needed.
                  !If we read a valid name, find it in the list of names
                  !and mark it for fine structure calculation.
                  do i=1,numedg(ncomps)
                    if (alledg(i,ncomps).eq.ename) then
                      finest(i,ncomps)=.true.
                    endif
                  enddo
                else
                  !Current line does not contain a valid edge name, so we 
                  !assume list of edge names has ended for this DETAIL card.
                  rdflag=.false.
                  exit
                endif
              enddo


            endif
          endif
        
        elseif (key.eq.'get') then !getcom_end or getcom_err
          if (str.eq.'getcom_end') then !EOF, we are done.
            exit 
          end if
          if (str.eq.'getcom_error') then
            call wlog ('Error in input file '//infile)
            call wlog ('Stopping.')
            stop
          endif
  
        else !Line not recognized as comment or card.
          call wlog('Ignoring unrecognized line: ')
          call wlog(line)
          rdflag=.true. !read next line

        endif

!       Need this test in case COMPONENT is last card in file
        if (str.eq.'getcom_end') exit !EOF, we are done.
        if (str.eq.'getcom_error') then
          call wlog ('Error in input file '//infile)
          call wlog ('Stopping.')
          stop
        endif

      enddo 
!*************************************     DONE READING CORE.INP

      
!     make sure at least one component is specified.
      if (ncomps.eq.0) then
        call wlog ('No Components specified.')
        call wlog ('Stopping.')
        stop
      endif

      do i=1,ncomps
        ldos(i)=.false.
        do j=1,numedg(i)
          if (conv(j,i)) then
            ldos(i)=.true.
          end if
        enddo
      enddo

!*************************************     INITIALIZE ENERGY GRID
      !Set egrid based on which edges we are calculating with fine
      !structure.
      lowedg=1.0e6
      hiedg=0.0
      if (emin.lt.0.0.or.emax.lt.0.0) then
        !This loop finds the highest and lowest energies of edges to be computed with fine structure.
        do i=1,ncomps
          iz=nz(i)
          do j=1,numedg(i)
            if (finest(j,i)) then
              ename=alledg(j,i)
              !get index of hole from the name of the edge
              call setedg (ename,ihole)
              !get energy of edge from hole index
              call getedg (ihole,iz,emudp)
              !convert to single precision
              emu=real(emudp)
              !revise highest and lowest edge energies
              if (emu.lt.lowedg) lowedg=emu
              if (emu.gt.hiedg) hiedg=emu
            endif
          enddo
        enddo
        !if we found some edges with fine structure, set the grid
        !boundaries 50 eV before the first edge to 1000 eV beyond the
        !last edge.
        if (lowedg.le.hiedg) then
          emin=max(lowedg-50/hart,0.1/hart)
          emax=hiedg+1000/hart
          !0.5 eV spacing
          iepts=(emax-emin)*hart*2
        end if
      endif

      !If the above loop did not set the grid, set egrid based on all edges we are calculating.
      lowedg=1.0e6
      hiedg=0.0
      if (emin.lt.0.0.or.emax.lt.0.0) then
        !This loop finds the highest and lowest energies of any edge.
        do i=1,ncomps
          iz=nz(i)
          do j=1,numedg(i)
            ename=alledg(j,i)
            !get index of hole from the name of the edge
            call setedg (ename,ihole)
            !get energy of edge from hole index
            call getedg (ihole,iz,emudp)
            !convert to single precision
            emu=real(emudp)
            !revise highest and lowest edge energies
            if (emu.lt.lowedg) lowedg=emu
            if (emu.gt.hiedg) hiedg=emu
          enddo
        enddo
        !if we found some edges, set the grid
        !boundaries 50 eV before the first edge to 1000 eV beyond the
        !last edge.
        if (lowedg.le.hiedg) then
          emin=max(lowedg-50/hart,0.1/hart)
          emax=hiedg+1000/hart
          !0.5 eV spacing
          iepts=(emax-emin)*hart*2
        end if
      endif

      if (emin.lt.0.0.or.emax.lt.0.0) then
        call wlog ('Trouble setting energy grid.')
        call wlog ('Stopping.')
        stop
      end if


      return
      end
