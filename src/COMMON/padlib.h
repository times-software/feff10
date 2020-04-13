! padlib.h -*-fortran-*-
!  header of parameters for packed-ascii-data (pad) routines
       implicit none
       character cpadr, cpadi, cpadc
       integer   maxlen, ibase, ioff, ihuge, ibas2
       double precision  ten, tenlog, huge, tiny, one, zero, base
       parameter(cpadr = '!', cpadc = '$', cpadi = '%')
       parameter(ibase = 90, ioff = 37, ihuge = 38, maxlen = 82)
       parameter(ibas2 = ibase/2, zero=0d0, one=1.d0, ten = 10.d0)
       parameter(tenlog= 2.302585092994045684d0)
       parameter(huge = ten**ihuge, tiny = one/huge)
       parameter(base = ibase*one)
!
