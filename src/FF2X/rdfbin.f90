!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: rdfbin.f90,v $:
! $Revision: 1.4 $
! $Author: jorissen $
! $Date: 2011/12/11 00:17:33 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine rdfbin(fbfile, nphx, nex, npathx, nlegx,              &
     &      npaths, ne, npot, ihole, iorder, ilinit,                    &
     &      rnrmav, xmu, edge,  potlbl, iz, phc, ck, xk, index,         &
     &      nleg, deg, reff, crit, ipot,                                &
     &      rat, beta, eta, ri, achi, phchi)

!
! read path information from PAD-format feff.bin
!  arguments:
!   fbile   name of feff.bin file                             [in]
!   nphx    dimension of  potlbl, iz (both (0:nphx)           [in]
!             max # of potentials
!   nex     dimension of many energy arrays                   [in]
!             max # of energy points
!   npathx     dimension of index,nleg,ipot,deg,reff,crit,    [in]
!           rat,beta, eta, ri, achi, phchi
!             max # of paths
!   nlegx  dimension of  ipot, rat, beta,eta,ri              [in]
!            max # of legs in a path
!   npaths  number of paths read                             [out]
!   ne      maximum number of energy points read             [out]
!   npot    number of potentials read                        [out]
!   rnrmav  average norman radius                            [out]
!   edge    shift in edge energy (?)                         [out]
!   iorder  order of genfmt matrix used                      [out]
!   potlbl  array of potential labels                        [out]
!   iz      array of atomic numbers for potentials           [out]
!   phc     array of central atom phase shift (complex)      [out]
!   ck      array of wavenumbers/momentum (complex)          [out]
!   xk      array of wavenumbers/momentum (real)             [out]
!   index   array of path indices                            [out]
!   nleg    array:  number of legs in path                   [out]
!   deg     array:  path degeneracy                          [out]
!   reff    array:  half path length of path                 [out]
!   crit    array:  importance factor for path               [out]
!   ipot    array:  pots, in order, that make up the path    [out]
!   rat     array:  atomic positions of atoms in path        [out]
!   beta    array:  euler angle for path                     [out]
!   eta     array:  second euler angle for path              [out]
!   ri      array:  path leg distances for path              [out]
!   achi    array:  amplitude of chi for path                [out]
!   phchi   array:  phase of chi for path                    [out]
!
! notes:
!   the data in feff.bin is written completely in printable
!   ascii characters.  The file is however, highly formatted
!   and kept fairly small.  all text is stored as is, but most
!   numerical data in arrays (both real and complex) is stored
!   in a special Packed Ascii Data (PAD) format which uses 6
!   printable characters to represent a real number.
!
!   special markers in the first 1 or 2 characters of each line
!   give hints about the contents of the line:
!      #_    top 2 lines.  The first line must begin "#_feff.bin"
!      #"    title lines / plain text
!      #&    misc info about potentials, calc method
!      #@    potential labels and iz
!      ##    path index,  deg, reff, crit, ipots involved
!      !     PAD characters to be read as a real array
!      $     PAD characters to be read as a complex array
!
! copyright (c) 1999  matt newville:  jan 1999
! modified by alex ankudinov: feb 2000; few fixes for feff8.2
!------------------------------------------------------------------
       use dimsmod,only:nwordx,nexdims=>nex
       integer nphx, nex, npathx, nlegx, npaths
       integer i, j, ivers, nexmax
       character*(*) fbfile
       character*(*) potlbl(0:nphx), filnam*128, str*512, msg*256
       integer iorder,  index(npathx), nleg(npathx)
       integer ne, npot, ipot(nlegx,npathx), iz(0:nphx)
       integer istrln, ier1, ier2, ier3, nwords, npadx !KJ , nwordx
       real    deg(npathx), reff(npathx), crit(npathx)
       real    rnrmav, edge, xk(nex)
       double precision bohr, tmpdp
       parameter (bohr = 0.529177249d0)  !KJ , nwordx = 100)
!       parameter(nexmax = 256)
       parameter(nexmax = nexdims)
       character*20 words(nwordx)
       real     rat(3,nlegx,npathx), beta(nlegx,npathx)
       real     eta(nlegx,npathx),  ri(nlegx,npathx)
       real     achi(nex,npathx), phchi(nex,npathx), tmpr(nexmax)
       complex  phc(nex), ck(nex), tmpc(nexmax)
       external  istrln

! open feff.bin
       filnam = ' '
       filnam = fbfile
       call triml(filnam)
       il     = istrln(filnam)
       if (filnam.eq.' ') filnam = 'feff.bin'
!c       print*, ' RDFBIN!  ', filnam(1:il),':',il
       open (unit=3, file=filnam, status='old', err=450)
 10    format(a)
! first line, must match  "#_feff.bin"
       read(3,10,err=920) str
       call triml(str)
       if ((str(1:10).ne.'#_feff.bin')) go to 900
! check version of feff.bin : only support version 3 here!!
       ivers = 0
       if ((str(1:14).eq.'#_feff.bin fil')) ivers = 1
       if ((str(1:14).eq.'#_feff.bin v02')) ivers = 2
       if ((str(1:14).eq.'#_feff.bin v03')) ivers = 3
       if (ivers.ne.3) go to 930
! second line:   npot, ne
       read(3,10,err=920) str
       call triml(str)
       if ((str(1:2).ne.'#_')) go to 900
       nwords = 3
       str    = str(3:)
       call bwords(str,nwords,words)
       if (nwords.ne.3) go to 905
       call str2in(words(1), npot,  ier1)
       call str2in(words(2), ne,    ier2)
       call str2in(words(3), npadx, ier3)
       if ((ier1.ne.0).or.(ier2.ne.0).or.(ier3.ne.0)) go to 910

!  read in misc stuff:  (rnrmav, edge, iorder )
       read(3,10,err=920) str
       call triml(str)
       if (str(1:2).ne.'#=') go to 900
       nwords = 6
       str    = str(3:)
       call bwords(str,nwords,words)
       if (nwords.ne.6) go to 905
       call str2in(words(1), ihole,  ier1)
       call str2in(words(2), iorder, ier2)
       call str2in(words(3), ilinit, ier3)
       if ((ier1.ne.0).or.(ier2.ne.0).or.(ier3.ne.0)) go to 910
       call str2re(words(4), rnrmav, ier1)
       call str2re(words(5), xmu   , ier2)
       call str2re(words(6), edge  , ier3)
       if ((ier1.ne.0).or.(ier2.ne.0).or.(ier3.ne.0)) go to 910
!  read pot label and iz line
       read(3,'(a512)',err=920) str
       call triml(str)
       if (str(1:2).ne.'#@') go to 900
       nwords = 2 * npot + 2
!  note: potlbl cannot be blank!!
       str    = str(3:)
       call bwords(str, nwords, words)
       if (nwords.ne.(2 + 2*npot)) go to 905
       do 200 i = 0, npot
          potlbl(i) = words(i+1)
          iz(i) = -1
          call str2in(words(2+npot+i),iz(i),ier1)
          if (ier1.ne.0)  go to 910
 200   continue

! read  numerical data that are the same for all paths
       call rdpadc(3,npadx, phc, ne)
       call rdpadc(3,npadx, ck,ne)
       call rdpadr(3,npadx, xk,ne)
       npaths = 0
! now, for each path:
       do 300  i = 1, npathx
          index(i) = 0
!  read path  info "##" line  and retrieve all the stuff from it
          read(3,10,end=450,err=920) str
          call triml(str)
          if (str(1:2).ne.'##')   go to 900
          nwords = nwordx
          str    = str(3:)
          call bwords(str,nwords,words)
          call str2in(words(1),  index(i), ier1)
          call str2in(words(2),  nleg(i), ier2)
          call str2re(words(3),  deg(i), ier3)
          if ((ier1.ne.0).or.(ier2.ne.0).or.(ier3.ne.0)) go to 910
          call str2dp(words(4),  tmpdp, ier2)
          reff(i) = tmpdp / bohr
          call str2dp(words(5),  tmpdp, ier3)
          crit(i) = tmpdp
          if ((ier1.ne.0).or.(ier2.ne.0).or.(ier3.ne.0)) go to 910
          npaths = npaths + 1
          do 230 j = 1, nleg(i)
             call str2in(words(5+j),ipot(j,i),ier1)
             if (ier1.ne.0) go to 910
 230      continue
!
!  next, read padded arrays for rat,beta, ..., achi, phchi
          call rdpadr(3,npadx, rat(1,1,i),3*nleg(i))
          call rdpadr(3,npadx, beta(1,i),   nleg(i))
          call rdpadr(3,npadx, eta(1,i),    nleg(i))
          call rdpadr(3,npadx, ri(1,i),     nleg(i))
          call rdpadr(3,npadx, achi(1, i),  ne)
          call rdpadr(3,npadx, phchi(1, i), ne)
!  fill in rest of achi and phchi with zeros
          do 270 j = ne+1, nex
             achi(j,i)  = 0
             phchi(j,i) = 0
 270      continue
 300    continue
 450    continue
       close(3)
!c       print*, ' RDFBIN done!'
       return
 900   call wlog (' -- rdfbin error: wrong format : at line')
       go to 990
 905   call wlog (' -- rdfbin error: missing data : at line')
       go to 990
 910   call wlog (' -- rdfbin error:   bad data   : at line')
       go to 990
 920   call wlog (' -- rdfbin error: unknown error: at line')
       go to 990
 930   call wlog (' -- rdfbin error: unknown version of feff.bin')
       go to 990

 990   call wlog (str)
       call par_stop(' -- fatal error reading feff.bin -- ')
       end


