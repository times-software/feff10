       subroutine fixstr(string,str,ilen,words,wrdsor,mwords,nwords)
!  simple preparation of string for reading of keywords
       integer       ilen, mwords, nwords, i, lenp1
       integer       iexcla, iperct, ihash, ieolc, istrln
       character*(*) string, str, words(mwords), wrdsor(mwords)
!
!  fix-up string: untab, left-justify, make a lower-case version
       nwords = 0
       call untab(string)
       str   = string
       call triml(str)
       call smcase( str, 'case')
!  remove comments from str:
!   '!', '#', and '%' are end of line comments
!   '*' is a complete comment line if in col 1
       lenp1  = len(str) + 1
       iexcla = index(str,'!')
       if (iexcla.eq.0)  iexcla = lenp1
       iperct = index(str,'%')
       if (iperct.eq.0)  iperct = lenp1
       ihash  = index(str,'#')
       if (ihash.eq.0)  ihash = lenp1
       ieolc  = min(iperct,iexcla,ihash) - 1
       if ((ieolc.lt.1).or.(str(1:1).eq.'*')) ieolc = 1
       str    = str(1:ieolc)
       ilen   = max(1, istrln(str))
!      if (ilen.le.2)  return
!  break string into words (up to mwords)
!  words is in lower case,   wrdsor is in original case
       do 120 i = 1, mwords
          words(i)   =  ' '
          wrdsor(i) =  ' '
 120   continue
       nwords = mwords
       call bwords(str   , nwords, words)
       call bwords(string, nwords, wrdsor)
! end  subroutine fixstr
       return
       end
