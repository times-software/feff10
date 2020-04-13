!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Information about last revision of $RCSfile: str.f90,v $:
! $Revision: 1.4 $
! $Author: hebhop $
! $Date: 2010/02/23 23:52:06 $
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION ISTRLN (STRING)  Returns index of last non-blank
!                           character.  Returns zero if string is
!                           null or all blank.

      FUNCTION ISTRLN (STRING)
      CHARACTER*(*)  STRING
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
!     there is a tab character here  ^

!  -- If null string or blank string, return length zero.
      ISTRLN = 0
      IF (STRING (1:1) .EQ. CHAR(0))  RETURN
      IF (STRING .EQ. ' ')  RETURN

!  -- Find rightmost non-blank character.
      ILEN = LEN (STRING)
      DO 20  I = ILEN, 1, -1
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 30
   20 CONTINUE
   30 ISTRLN = I

      RETURN
      END
! SUBROUTINE TRIML (STRING)  Removes leading blanks.

      SUBROUTINE TRIML (STRING)
      CHARACTER*(*)  STRING
      CHARACTER*512  TMP !KJ 7-09 increased from 200 to 512.  Would prefer a more dynamic solution FIX LATER.  Important for header feff.bin
      CHARACTER BLANK, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
!     there is a tab character here  ^

      JLEN = ISTRLN (STRING)

!  -- All blank and null strings are special cases.
      IF (JLEN .EQ. 0)  RETURN

!  -- FInd first non-blank char
      DO 10  I = 1, JLEN
         IF (STRING(I:I).NE.BLANK .AND. STRING(I:I).NE.TAB)  GOTO 20
   10 CONTINUE
   20 CONTINUE

!  -- If I is greater than JLEN, no non-blanks were found.
      IF (I .GT. JLEN)  RETURN

!  -- Remove the leading blanks.
      TMP = STRING (I:)
      STRING = TMP
      RETURN
      END
! SUBROUTINE UPPER (STRING)  Changes a-z to upper case.

      SUBROUTINE UPPER (STRING)
      CHARACTER*(*)  STRING

      JLEN = ISTRLN (STRING)

      DO 10  I = 1, JLEN
         IC = ICHAR (STRING (I:I))
         IF ((IC .LT. 97)  .OR.  (IC .GT. 122))  GOTO 10
         STRING (I:I) = CHAR (IC - 32)
   10 CONTINUE

      RETURN
      END
! SUBROUTINE LOWER (STRING)  Changes A-Z to lower case.

      SUBROUTINE LOWER (STRING)
      CHARACTER*(*)  STRING

      JLEN = ISTRLN (STRING)

      DO 10  I = 1, JLEN
         IC = ICHAR (STRING (I:I))
         IF ((IC .LT. 65) .OR.  (IC .GT. 90))  GOTO 10
         STRING (I:I) = CHAR (IC + 32)
   10 CONTINUE

      RETURN
      END
!***********************************************************************
!
      SUBROUTINE BWORDS (S, NWORDS, WORDS)
!
!     Breaks string into words.  Words are seperated by one or more
!     blanks or tabs, or a comma and zero or more blanks.
!
!     ARGS        I/O      DESCRIPTION
!     ----        ---      -----------
!     S            I       CHAR*(*)  String to be broken up
!     NWORDS      I/O      Input:  Maximum number of words to get
!                          Output: Number of words found
!     WORDS(NWORDS) O      CHAR*(*) WORDS(NWORDS)
!                          Contains words found.  WORDS(J), where J is
!                          greater then NWORDS found, are undefined on
!                          output.
!
!      Written by:  Steven Zabinsky, September 1984
!      Tab char added July 1994.
!
!**************************  Deo Soli Gloria  **************************

!  -- No floating point numbers in this routine.
      IMPLICIT INTEGER (A-Z)

      CHARACTER*(*) S, WORDS(NWORDS)

      CHARACTER BLANK, COMMA, TAB
      PARAMETER (BLANK = ' ', COMMA = ',', TAB = '	')
!     there is a tab character here               ^.

!  -- BETW    .TRUE. if between words
!     COMFND  .TRUE. if between words and a comma has already been found
      LOGICAL BETW, COMFND

!  -- Maximum number of words allowed
      WORDSX = NWORDS

!  -- SLEN is last non-blank character in string
      SLEN = ISTRLN (S)
!  -- All blank string is special case
      IF (SLEN .EQ. 0)  THEN
         NWORDS = 0
         RETURN
      ENDIF

!  -- BEGC is beginning character of a word
      BEGC = 1
      NWORDS = 0

      BETW   = .TRUE.
      COMFND = .TRUE.

      DO 10  I = 1, SLEN
         IF (S(I:I) .EQ. BLANK .OR. S(I:I) .EQ. TAB)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S (BEGC : I-1)
               BETW = .TRUE.
               COMFND = .FALSE.
            ENDIF
         ELSEIF (S(I:I) .EQ. COMMA)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S(BEGC : I-1)
               BETW = .TRUE.
            ELSEIF (COMFND)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = BLANK
            ENDIF
            COMFND = .TRUE.
         ELSE
            IF (BETW)  THEN
               BETW = .FALSE.
               BEGC = I
            ENDIF
         ENDIF

         IF (NWORDS .GE. WORDSX)  RETURN

   10 CONTINUE

      IF (.NOT. BETW  .AND.  NWORDS .LT. WORDSX)  THEN
         NWORDS = NWORDS + 1
         WORDS (NWORDS) = S (BEGC :SLEN)
      ENDIF

      RETURN
      END

!***********************************************************************
!
      SUBROUTINE BWORDS2 (S, NWORDS, WORDS)
!
!     Breaks string into words.  Words are seperated by one or more
!     blanks or tabs.
!
!     ARGS        I/O      DESCRIPTION
!     ----        ---      -----------
!     S            I       CHAR*(*)  String to be broken up
!     NWORDS      I/O      Input:  Maximum number of words to get
!                          Output: Number of words found
!     WORDS(NWORDS) O      CHAR*(*) WORDS(NWORDS)
!                          Contains words found.  WORDS(J), where J is
!                          greater than NWORDS found, are undefined on
!                          output.
!
!      Written by:  Steven Zabinsky, September 1984
!      Tab char added July 1994.
!
!**************************  Deo Soli Gloria  **************************

!  -- No floating point numbers in this routine.
      IMPLICIT INTEGER (A-Z)

      CHARACTER*(*) S, WORDS(NWORDS)

      CHARACTER BLANK, COMMA, TAB
      PARAMETER (BLANK = ' ', TAB = '	')
!     there is a tab character here               ^.

!  -- BETW    .TRUE. if between words
!     COMFND  .TRUE. if between words and a comma has already been found
      LOGICAL BETW, COMFND

!  -- Maximum number of words allowed
      WORDSX = NWORDS

!  -- SLEN is last non-blank character in string
      SLEN = ISTRLN (S)

!  -- All blank string is special case
      IF (SLEN .EQ. 0)  THEN
         NWORDS = 0
         RETURN
      ENDIF

!  -- BEGC is beginning character of a word
      BEGC = 1
      NWORDS = 0

      BETW   = .TRUE.
      COMFND = .TRUE.

      DO 10  I = 1, SLEN
         IF (S(I:I) .EQ. BLANK .OR. S(I:I) .EQ. TAB)  THEN
            IF (.NOT. BETW)  THEN
               NWORDS = NWORDS + 1
               WORDS (NWORDS) = S (BEGC : I-1)
               BETW = .TRUE.
               COMFND = .FALSE.
            ENDIF
         ELSE
            IF (BETW)  THEN
               BETW = .FALSE.
               BEGC = I
            ENDIF
         ENDIF

         IF (NWORDS .GE. WORDSX)  RETURN

   10 CONTINUE

      IF (.NOT. BETW  .AND.  NWORDS .LT. WORDSX)  THEN
         NWORDS = NWORDS + 1
         WORDS (NWORDS) = S (BEGC :SLEN)
      ENDIF

      RETURN
      END


      SUBROUTINE UNTAB (STRING)
! REPLACE TABS WITH BLANKS :    TAB IS ASCII DEPENDENT
      INTEGER        ITAB , I, ILEN, ISTRLN
      PARAMETER      (ITAB = 9)
      CHARACTER*(*)  STRING, TAB*1
      EXTERNAL ISTRLN
      TAB  = CHAR(ITAB)
      ILEN = MAX(1, ISTRLN(STRING))
 10   CONTINUE
        I = INDEX(STRING(:ILEN), TAB )
        IF (I .NE. 0) THEN
            STRING(I:I) = ' '
            GO TO 10
        END IF
      RETURN
! END SUBROUTINE UNTAB
      END

      logical function iscomm (line)
!     returns true if line is a comment or blank line, false otherwise
!#mn{ rewritten to allow ";*%#" as comment characters
       character*(*) line
       iscomm = ((line.eq.' ').or.(index(';*%#',line(1:1)).ne.0))
!#mn}
      return
      end
      subroutine smcase (str, contrl)
!  convert case of string *str*to be the same case
!  as the first letter of string *contrl*
!  if contrl(1:1) is not a letter, *str* will be made lower case.
      character*(*) str, contrl, s1*1, t1*1
      s1 = contrl(1:1)
      t1 = s1
      call lower(t1)
      if (t1.eq.s1)  call lower(str)
      if (t1.ne.s1)  call upper(str)
      return
! end subroutine smcase
      end
       subroutine unblnk (string)
!
! remove blanks from a string
       integer        i, ilen, j
       character*(*)  string, str*256, blank*1
       parameter (blank = ' ')
       ilen = min(256, max(1, istrln(string)))
       j   = 0
       str = blank
       do 10 i = 1, ilen
         if (string(i:i).ne.blank) then
            j = j+1
            str(j:j) = string(i:i)
         end if
 10   continue
      string = blank
      string = str(1:j)
      return
! end subroutine unblnk
      end
      subroutine uncomm(str)
!
! purpose: remove comments from a string
!
! arguments:
!      str  string to modify        [in/out]
! notes:
!   1. '*' is a comment iff it occurs in col 1
!   2. char(10) and char(12) are end-of-line comments
!   3. '!', '#', and '%'  are end-of-line comments that
!       can be protected by matching " ", ' ', ( ), [], or {}
!
! requires:  istrln, triml, echo
!
! copyright 1997  matt newville
       integer i, istrln, ilen, iprot
       character*(*) str, copen*5, cclose*5, eol*3, spec*2, s*1
       character*1 blank, star
       parameter(blank = ' ',star = '*')
       external  istrln
       data copen, cclose, eol  / '[{"''(',  ']}"'')', '!#%' /
!
       spec(1:2) = char(10)//char(12)
       call triml(str)
       ilen = istrln(str)
       if ((ilen.le.0).or.(str(1:1).eq.star)) then
          str = blank
          i   = 1
       else
          iprot = 0
          do 50 i = 1, ilen
             s  = str(i:i)
             if (iprot.le.0) then
                iprot = index(copen,s)
             elseif (iprot.le.5) then
                if (s.eq.cclose(iprot:iprot)) iprot = 0
             else
!c                call echo('** uncomm confusion: iprot out of range')
                return
             end if
! if the string is unprotected, look for end-of-line comment characters
             if (((iprot.eq.0).and.(index(eol,s).ne.0)).or.             &
     &             index(spec,s).ne.0)  go to 60
 50       continue
          i  = ilen + 1
 60       continue
       end if
       str  = str(1:i-1)
! end subroutine uncomm
       return
       end
      subroutine strclp(str,str1,str2,strout)
!
!  a rather complex way of clipping a string:
!      strout = the part of str that begins with str2.
!  str1 and str2 are subsrtings of str, (str1 coming before str2),
!  and even if they are similar, strout begins with str2
!  for example:
!   1.  str =  "title title my title" with  str1 = str2 = "title"
!       gives strout = "title my title"
!   2.  str =  "id  1  1st path label" with str1 = "1", str2 = "1st"
!       gives strout = "1st path label"
!
      character*(*)  str, str1, str2, strout
      integer  i1, i2, ibeg, iend, istrln, ilen
      external istrln
      ilen   = len(strout)
      i1     = max(1, istrln(str1))
      i2     = max(1, istrln(str2))
      i1e    = index(str,str1(1:i1)) + i1
      ibeg   = index(str(i1e:),str2(1:i2) ) + i1e - 1
      iend   = min(ilen+ibeg, istrln(str) )
      strout = str(ibeg:iend)
      return
! end subroutine strclp
      end
       subroutine rmdels(s,s1,s2)
!
!  remove general enclosing delimeters from a string
       character*(*) s, s1, s2, t*512
       call triml(s)
       i  = istrln(s)
       t  = s
       if ((s(1:1).eq.s1) .and. (s(i:i).eq.s2)) s = t(2:i-1)
       return
       end
!
!        subroutine rmpars(str)
! c  remove enclosing parentheses for a string
!        character*(*) str
!        call rmdels(str,'(',')')
!        return
!        end

       subroutine rmquot(str)
!  remove enclosing single or double quotes from a string
       character*(*) str
       call rmdels(str,'''','''')
       call rmdels(str,'"','"')
       return
       end
       subroutine undels(s)
!  remove an enclosing delimiter from a string
       character*(*) s, op*5, cl*5
       integer j
       data op, cl / '[{"''(',  ']}"'')'/
       j = index(op,s(1:1))
       if (j.ne.0) then
          call rmdels(s, op(j:j), cl(j:j) )
       end if
       return
       end
      subroutine str2lg(str,flag,ierr)
!  return logical "flag" from character string "str".
!  flag is true unless the first character is
!     '0', 'f' or 'n' (not case-sensitive)
      character*(*) str, test*5
      parameter (test = 'fnFN0')
      logical    flag
      integer    ierr
      ierr  = 0
      flag  = index(test,str(1:1)).eq.0
      return
! end subroutine str2lg
      end
!        subroutine str2il(str,miar,niar,iar,ierr)
! !  convert a string into an integer _list_,
! !  supporting syntax like '1-2,12,4,6-8' returns
! !  iar =   1,2,4,6,7,8,12    niar = 7
! !
! !  returns ierr = -1 if string clearly non-integer
!        character*(*) str , s*128, sint*32
!        integer  miar, niar, iar(miar), ierr, istrln
!        integer  i, ibeg
!        logical  dash
!        external  istrln
!        s    = str
!        call triml(s)
!        ilen = istrln(s)+1
!        s    = s(1:ilen-1)//'^'
!        do 20 i = 1, miar
!           iar(i) = 0
!  20    continue
!        niar =  0
!        ierr = -1
!        ix1  =  0
!        dash = .false.
!        if (ilen.gt.1) then
!           i    = 1
!           ibeg = 1
!  100      continue
!           i = i + 1
!           if ((s(i:i).eq.',') .or. (s(i:i).eq.'^')) then
!              sint = s(ibeg:i-1)
!              ibeg = i+1
!              if (dash) then
!                 call str2in(sint,ix,ierr)
!                 do 130 j = ix1, ix
!                    niar = niar + 1
!                    iar(niar) = j
!  130            continue
!              else
!                 call str2in(sint,ix,ierr)
!                 niar = niar + 1
!                 iar(niar) = ix
!              end if
!              dash = .false.
!           elseif (s(i:i).eq.'-') then
!              sint = s(ibeg:i-1)
!              dash = .true.
!              call str2in(sint,ix1,ierr)
!              ibeg = i+1
!           end if
!           if (s(i:i).ne.'^') go to 100
!        end if
! ! now remove the zeroth one!
!        niar = niar - 1
! !
!        return
! ! end subroutine str2il
!        end

       logical function isnum (string)
!  tests whether a string can be a number. not foolproof!
!  to return true, string must contain:
!    - only characters in  'deDE.+-, 1234567890' (case is checked)
!    - no more than one 'd' or 'e'
!    - no more than one '.'
!    - if '+' or '-' is seen after a digit, 'deDE' must be seen.
!  matt newville
       character*(*)  string,  number*20
! note:  layout and case of *number* is important: do not change!
       parameter (number = 'deDE.,+- 1234567890')
       integer   iexp, idec, i, j, istrln, isign
       integer   jexp, jsign
       logical   ldig, l_op
       external  istrln
!       str   = string
!       call triml(str)
       iexp  = 0
       jexp  = 0
       idec  = 0
       isign = 0
       ldig  = .false.
       l_op  = .false.
       isnum = .false.
       do 100  i = 1, max(1, istrln(string))
          j = index(number,string(i:i))
!c          print*, 'X  ' , i, j, ' : ' , str(i:i)
          if (j.le.0)               go to 200
          if (j.ge.10)              ldig = .true.
          if((j.ge.1).and.(j.le.4)) then
             iexp = iexp + 1
             jexp = i
          endif
          if (j.eq.5)               idec = idec + 1
          if ((j.eq.7).or.(j.eq.8)) then
             isign= isign +1
             if ((i .gt. 1) .and. (i .ne. (jexp+1))) then
                l_op = .true.
             endif
          endif
 100   continue
!  every character in "string" is also in "number".  so, if there are
!  not more than one exponential and decimal markers, it's a number
       if ((iexp.le.1).and.(idec .le.1)) isnum = .true.
       if ((iexp.eq.0).and.(isign.gt.1)) isnum = .false.
       if (jexp.eq.1)  isnum = .false.
       isnum = isnum .and. (.not.l_op)
!c       print*, 'ISNUM: ', string(1:istrln(string))
!c       print*, '       ', isnum, l_op, iexp, idec, isign
 200   continue
       return
!  end logical function isnum
       end
       logical function isdat(string)
!  tests if string contains numerical data
!    returns true if the first (up to eight) words in string can
!    all be numbers. requires at least two words, and tests only
!    the first eight columns
       integer nwords, mwords, i
       parameter (mwords = 8)
       character*(30)  string*(*), words(mwords), line*(256)
       logical isnum
       external isnum
!
       isdat = .false.
       do 10 i = 1, mwords
          words(i) = 'no'
 10    continue
!
       nwords = mwords
       line   = string
       call triml(line)
       call untab(line)
       call bwords(line, nwords, words)
       if (nwords.ge.1) then
          isdat = .true.
          do 50 i = 1, nwords
             isdat = isdat .and. isnum(words(i))
 50       continue
       end if
       return
       end
       subroutine bkeys(str, mkeys, keys, values, nkeys)
!
! purpose:  break a string into {key,value} pairs.
! arguments:
!      str     string to break into pairs           [in]
!      mkeys   dimension of arrays keys and values  [in]
!      keys    character array of keys              [out]
!      values  character array of values            [out]
!      nkeys   number of keys found                 [out]
!
! parsing rules:
!  1. a key is a word terminated by whitespace, an equal sign,
!     a comma, or the final close paren.  keys are converted to
!     lower case before returning.
!
!  2. a value is a more general string, terminated by either
!     an "unprotected" comma or the final "unprotected" close paren.
!     Any part of the string can be "protected" by either matching
!     single quotes, double quotes, parens, braces, or brackets.
!     In fact, *all* of these pairs must be matched for the
!     value to terminate.  the values are left in their original case.
!
!  3. If a key does not have a value (because a comma or the last close
!     paren gets in the way) the value will be set to '%undef%'.
!     note that str2lg will interpret this as "true"!, and that it
!     will never make sense as any other value.
!
! example:  x =13.214, File = B.dat, Verbose, sig = sqrt(A + min(b,c))
!   will return these pairs:
!        key        value
!        x          13.214
!        file       B.dat
!        verbose    %undef%
!        sig        sqrt(A + min(b,c))
!
!  routines needed: istrln, triml, lower, rmdels, echo
!
!  copyright (c) 1998  matt newville
!
       integer   istrln, i, j, ilen, ibeg
       integer   nkeys, mkeys, nk, jprot
       character*(*) str, keys(mkeys), values(mkeys), tmp*64
       character s, t, u, blank, comma, equal, semicl
       character copen*3, cclose*3, undef*8
       logical   lcomma, seek_key, have_key
       parameter (blank = ' ',comma = ',',equal = '=',semicl = ';')
       parameter (undef = '%undef%')
       external istrln
       data copen, cclose / '[{(',  ']})'/
!
! initialize
       nkeys = 0
       do 10 i = 1, mkeys
          keys(i)   = blank
          values(i) = undef
 10    continue
       have_key = .false.
       seek_key = .true.
       lcomma   = .false.
       ibeg     = 1
       iprot    = 0
       jprot    = 0
!
! check for valid string to parse
       ilen = istrln(str)
!c       print*,'BKEYS:',str(1:ilen),':'
       if (ilen .eq. 0)  return
!
! loop through string
       i = 0
 100   continue
       i = i + 1
       s  = str(i:i)
! test for opening/closing delimiters
! and march over protected strings
       if ((s.eq.'''').or.(s.eq.'"')) then
          t = s
!c          print*, ' quote: ', t
 120      continue
          i  = i + 1
          if ((str(i:i).ne.t).and.(i.lt.ilen)) goto 120
       else
          iprot = index(copen,s)
          if ((iprot.ge.1).and.(iprot.le.3)) then
!c             print*, ' iprot = ',iprot , s, i
             jprot= jprot + 1
             t = copen(iprot:iprot)
             u = cclose(iprot:iprot)
 130         continue
             i  = i + 1
             if (str(i:i).eq.t)  jprot = jprot + 1
             if (str(i:i).eq.u)  jprot = jprot - 1
             if ((i.lt.ilen).and.(jprot.ne.0)) goto 130
          end if
       endif
       lcomma = s.eq.comma
! looking for keyword:
!   we've seen the beginning of a keyword, and now we see the end:
!   keyword  ends at "=",","," ", or the final positon
       if (seek_key) then
          if (((s.eq.equal).or.lcomma.or.(i.eq.ilen))) then
             nkeys  = nkeys + 1
             if (nkeys .ge. mkeys) go to 150
             keys(nkeys) = str(ibeg:i-1)
             if ((i.eq.ilen).and.(.not.lcomma).and.(s.ne.equal))        &
     &            keys(nkeys) = str(ibeg:i)
!c             print*, 'found key : ', nkeys, ' ', keys(nkeys)(1:32)
             ibeg   = min(i + 1, ilen)
             seek_key = .false.
             have_key = .false.
!      a bare word counts as a key with value= undefined (as above)
             if (lcomma .or.(i.eq.ilen) ) then
                seek_key = .true.
                call triml(keys(nkeys))
                ij = istrln(keys(nkeys))
                if  (index(keys(nkeys)(1:ij),blank).ne.0) then
                   tmp = keys(nkeys)(1:ij)
!      c                        call echo(' syntax error: '//tmp)
                   keys(nkeys)  = blank
                end if
             end if
          elseif (.not.have_key) then
             have_key = s.ne.blank
          end if
!      looking for a value:  ends at a comma or the final postion
       else
          if (lcomma.or.(i.eq.ilen)) then
             values(nkeys) = str(ibeg:i-1)
             if ((i.eq.ilen).and.(.not.lcomma))                         &
     &            values(nkeys) = str(ibeg:)
             ibeg   = min( i + 1, ilen)
             seek_key = .true.
          end if
       end if
       if (i.le.ilen) goto 100
 150   continue
!
!  finally, we may have ended with a one-letter keyword, in which case
!   have_key is true
       if (have_key) then
          nkeys       = nkeys + 1
          keys(nkeys) = str(ibeg:)
          call triml(keys(nkeys))
       end if
!
! now clean up keys and values, eliminate blank and invalid keys
       nk = nkeys
       nkeys = 0
       do 500 i = 1, nk
          if (keys(i).ne.blank .and. keys(i).ne.comma .and.             &
     &         keys(i).ne.equal .and. keys(i).ne.semicl) then
             nkeys = nkeys + 1
             keys(nkeys) = keys(i)
             call triml( values(i))
             if (values(i)(1:1).eq.equal) then
                values(i) = values(i)(2:)
                call triml(values(i) )
             end if
             call rmquot(values(i))
             do 470 j = 1, 3
                call rmdels(values(i),copen(j:j),cclose(j:j))
 470         continue
             call triml( values(i))
             values(nkeys) = values(i)
             if (values(nkeys).ne.undef) call lower(keys(nkeys))
             call triml(keys(nkeys))
          end if
          lk = istrln(keys(i))
          lv = istrln(values(i))
!c          print*, i,' |', keys(i)(1:lk),' | ', values(i)(1:lv), '|'
 500   continue
       return
! end subroutine bkeys
       end
