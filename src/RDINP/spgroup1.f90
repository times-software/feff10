 subroutine getlattype(sgname, lattyp)
! This routine gives the Bravais lattice type for a given spacegroup.
     implicit none
	 character*3,intent(out) :: lattyp
	 character*8,intent(inout) :: sgname
	 integer i

!    Correct for spellings like 'C2/C' which should be 'C2/c'
     do i=2,8
	    if (sgname(i:i).eq.'M') sgname(i:i)='m'
	    if (sgname(i:i).eq.'A') sgname(i:i)='a'
	    if (sgname(i:i).eq.'B') sgname(i:i)='b'
	    if (sgname(i:i).eq.'C') sgname(i:i)='c'
	 enddo
	 
     lattyp='   '
     if (sgname.eq.'P1      ') lattyp='P  '
     if (sgname.eq.'P-1     ') lattyp='P  '
     if (sgname.eq.'P2      ') lattyp='P  '
     if (sgname.eq.'P2      ') lattyp='P  '
     if (sgname.eq.'P2      ') lattyp='P  '
     if (sgname.eq.'P21     ') lattyp='P  '
     if (sgname.eq.'P21     ') lattyp='P  '
     if (sgname.eq.'P21     ') lattyp='P  '
     if (sgname.eq.'C2      ') lattyp='CXY'
     if (sgname.eq.'A2      ') lattyp='CYZ'
     if (sgname.eq.'B2      ') lattyp='CXZ'
     if (sgname.eq.'B2      ') lattyp='CXZ'
     if (sgname.eq.'C2      ') lattyp='CXY'
     if (sgname.eq.'A2      ') lattyp='CYZ'
     if (sgname.eq.'Pm      ') lattyp='P  '
     if (sgname.eq.'Pm      ') lattyp='P  '
     if (sgname.eq.'Pm      ') lattyp='P  '
     if (sgname.eq.'Pc      ') lattyp='P  '
     if (sgname.eq.'Pa      ') lattyp='P  '
     if (sgname.eq.'Pb      ') lattyp='P  '
     if (sgname.eq.'Pb      ') lattyp='P  '
     if (sgname.eq.'Pc      ') lattyp='P  '
     if (sgname.eq.'Pa      ') lattyp='P  '
     if (sgname.eq.'Pn      ') lattyp='P  '
     if (sgname.eq.'Pn      ') lattyp='P  '
     if (sgname.eq.'Pn      ') lattyp='P  '
     if (sgname.eq.'Cm      ') lattyp='CXY'
     if (sgname.eq.'Am      ') lattyp='CYZ'
     if (sgname.eq.'Bm      ') lattyp='CXZ'
     if (sgname.eq.'Bm      ') lattyp='CXZ'
     if (sgname.eq.'Cm      ') lattyp='CXY'
     if (sgname.eq.'Am      ') lattyp='CYZ'
     if (sgname.eq.'Cc      ') lattyp='CXY'
     if (sgname.eq.'Aa      ') lattyp='CYZ'
     if (sgname.eq.'Bb      ') lattyp='CXZ'
     if (sgname.eq.'Bb      ') lattyp='CXZ'
     if (sgname.eq.'Cc      ') lattyp='CXY'
     if (sgname.eq.'Aa      ') lattyp='CYZ'
     if (sgname.eq.'P2/m    ') lattyp='P  '
     if (sgname.eq.'P2/m    ') lattyp='P  '
     if (sgname.eq.'P2/m    ') lattyp='P  '
     if (sgname.eq.'P21/m   ') lattyp='P  '
     if (sgname.eq.'P21/m   ') lattyp='P  '
     if (sgname.eq.'P21/m   ') lattyp='P  '
     if (sgname.eq.'C2/m    ') lattyp='CXY'
     if (sgname.eq.'A2/m    ') lattyp='CYZ'
     if (sgname.eq.'B2/m    ') lattyp='CXZ'
     if (sgname.eq.'B2/m    ') lattyp='CXZ'
     if (sgname.eq.'C2/m    ') lattyp='CXY'
     if (sgname.eq.'A2/m    ') lattyp='CYZ'
     if (sgname.eq.'P2/c    ') lattyp='P  '
     if (sgname.eq.'P2/a    ') lattyp='P  '
     if (sgname.eq.'P2/b    ') lattyp='P  '
     if (sgname.eq.'P2/b    ') lattyp='P  '
     if (sgname.eq.'P2/c    ') lattyp='P  '
     if (sgname.eq.'P2/a    ') lattyp='P  '
     if (sgname.eq.'P2/n    ') lattyp='P  '
     if (sgname.eq.'P2/n    ') lattyp='P  '
     if (sgname.eq.'P2/n    ') lattyp='P  '
     if (sgname.eq.'P21/c   ') lattyp='P  '
     if (sgname.eq.'P21/a   ') lattyp='P  '
     if (sgname.eq.'P21/b   ') lattyp='P  '
     if (sgname.eq.'P21/b   ') lattyp='P  '
     if (sgname.eq.'P21/c   ') lattyp='P  '
     if (sgname.eq.'P21/a   ') lattyp='P  '
     if (sgname.eq.'P21/n   ') lattyp='P  '
     if (sgname.eq.'P21/n   ') lattyp='P  '
     if (sgname.eq.'P21/n   ') lattyp='P  '
     if (sgname.eq.'C2/c    ') lattyp='CXY'
     if (sgname.eq.'A2/a    ') lattyp='CYZ'
     if (sgname.eq.'B2/b    ') lattyp='CXZ'
     if (sgname.eq.'B2/b    ') lattyp='CXZ'
     if (sgname.eq.'C2/c    ') lattyp='CXY'
     if (sgname.eq.'A2/a    ') lattyp='CYZ'
     if (sgname.eq.'P222    ') lattyp='P  '
     if (sgname.eq.'P2221   ') lattyp='P  '
     if (sgname.eq.'P2122   ') lattyp='P  '
     if (sgname.eq.'P2212   ') lattyp='P  '
     if (sgname.eq.'P21212  ') lattyp='P  '
     if (sgname.eq.'P22121  ') lattyp='P  '
     if (sgname.eq.'P21221  ') lattyp='P  '
     if (sgname.eq.'P212121 ') lattyp='P  '
     if (sgname.eq.'C2221   ') lattyp='CXY'
     if (sgname.eq.'A2122   ') lattyp='CYZ'
     if (sgname.eq.'B2212   ') lattyp='CXZ'
     if (sgname.eq.'C222    ') lattyp='CXY'
     if (sgname.eq.'A222    ') lattyp='CYZ'
     if (sgname.eq.'B222    ') lattyp='CXZ'
     if (sgname.eq.'F222    ') lattyp='F  '
     if (sgname.eq.'I222    ') lattyp='B  '
     if (sgname.eq.'I212121 ') lattyp='B  '
     if (sgname.eq.'Pmm2    ') lattyp='P  '
     if (sgname.eq.'P2mm    ') lattyp='P  '
     if (sgname.eq.'Pm2m    ') lattyp='P  '
     if (sgname.eq.'Pmc21   ') lattyp='P  '
     if (sgname.eq.'P21ma   ') lattyp='P  '
     if (sgname.eq.'Pb21m   ') lattyp='P  '
     if (sgname.eq.'Pcm21   ') lattyp='P  '
     if (sgname.eq.'P21am   ') lattyp='P  '
     if (sgname.eq.'Pm21b   ') lattyp='P  '
     if (sgname.eq.'Pcc2    ') lattyp='P  '
     if (sgname.eq.'P2aa    ') lattyp='P  '
     if (sgname.eq.'Pb2b    ') lattyp='P  '
     if (sgname.eq.'Pma2    ') lattyp='P  '
     if (sgname.eq.'P2mb    ') lattyp='P  '
     if (sgname.eq.'Pc2m    ') lattyp='P  '
     if (sgname.eq.'Pbm2    ') lattyp='P  '
     if (sgname.eq.'P2cm    ') lattyp='P  '
     if (sgname.eq.'Pm2a    ') lattyp='P  '
     if (sgname.eq.'Pca21   ') lattyp='P  '
     if (sgname.eq.'P21ab   ') lattyp='P  '
     if (sgname.eq.'Pc21b   ') lattyp='P  '
     if (sgname.eq.'Pbc21   ') lattyp='P  '
     if (sgname.eq.'P21ca   ') lattyp='P  '
     if (sgname.eq.'Pb21a   ') lattyp='P  '
     if (sgname.eq.'Pnc2    ') lattyp='P  '
     if (sgname.eq.'P2na    ') lattyp='P  '
     if (sgname.eq.'Pb2n    ') lattyp='P  '
     if (sgname.eq.'Pcn2    ') lattyp='P  '
     if (sgname.eq.'P2an    ') lattyp='P  '
     if (sgname.eq.'Pn2b    ') lattyp='P  '
     if (sgname.eq.'Pmn21   ') lattyp='P  '
     if (sgname.eq.'P21mn   ') lattyp='P  '
     if (sgname.eq.'Pn21m   ') lattyp='P  '
     if (sgname.eq.'Pnm21   ') lattyp='P  '
     if (sgname.eq.'P21nm   ') lattyp='P  '
     if (sgname.eq.'Pm21n   ') lattyp='P  '
     if (sgname.eq.'Pba2    ') lattyp='P  '
     if (sgname.eq.'P2cb    ') lattyp='P  '
     if (sgname.eq.'Pc2a    ') lattyp='P  '
     if (sgname.eq.'Pna21   ') lattyp='P  '
     if (sgname.eq.'P21nb   ') lattyp='P  '
     if (sgname.eq.'Pc21n   ') lattyp='P  '
     if (sgname.eq.'Pbn21   ') lattyp='P  '
     if (sgname.eq.'P21cn   ') lattyp='P  '
     if (sgname.eq.'Pn21a   ') lattyp='P  '
     if (sgname.eq.'Pnn2    ') lattyp='P  '
     if (sgname.eq.'P2nn    ') lattyp='P  '
     if (sgname.eq.'Pn2n    ') lattyp='P  '
     if (sgname.eq.'Cmm2    ') lattyp='CXY'
     if (sgname.eq.'A2mm    ') lattyp='CYZ'
     if (sgname.eq.'Bm2m    ') lattyp='CXZ'
     if (sgname.eq.'Cmc21   ') lattyp='CXY'
     if (sgname.eq.'A21ma   ') lattyp='CYZ'
     if (sgname.eq.'Bb21m   ') lattyp='CXZ'
     if (sgname.eq.'Ccm21   ') lattyp='CXY'
     if (sgname.eq.'A21am   ') lattyp='CYZ'
     if (sgname.eq.'Bm21b   ') lattyp='CXZ'
     if (sgname.eq.'Ccc2    ') lattyp='CXY'
     if (sgname.eq.'A2aa    ') lattyp='CYZ'
     if (sgname.eq.'Bb2b    ') lattyp='CXZ'
     if (sgname.eq.'Amm2    ') lattyp='CYZ'
     if (sgname.eq.'B2mm    ') lattyp='CXZ'
     if (sgname.eq.'Cm2m    ') lattyp='CXY'
     if (sgname.eq.'Abm2    ') lattyp='CYZ'
     if (sgname.eq.'B2cm    ') lattyp='CXZ'
     if (sgname.eq.'Cm2a    ') lattyp='CXY'
     if (sgname.eq.'Bma2    ') lattyp='CXZ'
     if (sgname.eq.'C2mb    ') lattyp='CXY'
     if (sgname.eq.'Ac2m    ') lattyp='CYZ'
     if (sgname.eq.'Ama2    ') lattyp='CYZ'
     if (sgname.eq.'B2mb    ') lattyp='CXZ'
     if (sgname.eq.'Cc2m    ') lattyp='CXY'
     if (sgname.eq.'Bbm2    ') lattyp='CXZ'
     if (sgname.eq.'C2cm    ') lattyp='CXY'
     if (sgname.eq.'Am2a    ') lattyp='CYZ'
     if (sgname.eq.'Aba2    ') lattyp='CYZ'
     if (sgname.eq.'B2cb    ') lattyp='CXZ'
     if (sgname.eq.'Cc2a    ') lattyp='CXY'
     if (sgname.eq.'Bba2    ') lattyp='CXZ'
     if (sgname.eq.'C2cb    ') lattyp='CXY'
     if (sgname.eq.'Ac2a    ') lattyp='CYZ'
     if (sgname.eq.'Fmm2    ') lattyp='F  '
     if (sgname.eq.'F2mm    ') lattyp='F  '
     if (sgname.eq.'Fm2m    ') lattyp='F  '
     if (sgname.eq.'Fdd2    ') lattyp='F  '
     if (sgname.eq.'F2dd    ') lattyp='F  '
     if (sgname.eq.'Fd2d    ') lattyp='F  '
     if (sgname.eq.'Imm2    ') lattyp='B  '
     if (sgname.eq.'I2mm    ') lattyp='B  '
     if (sgname.eq.'Im2m    ') lattyp='B  '
     if (sgname.eq.'Iba2    ') lattyp='B  '
     if (sgname.eq.'I2cb    ') lattyp='B  '
     if (sgname.eq.'Ic2a    ') lattyp='B  '
     if (sgname.eq.'Ima2    ') lattyp='B  '
     if (sgname.eq.'I2mb    ') lattyp='B  '
     if (sgname.eq.'Ic2m    ') lattyp='B  '
     if (sgname.eq.'Ibm2    ') lattyp='B  '
     if (sgname.eq.'I2cm    ') lattyp='B  '
     if (sgname.eq.'Im2a    ') lattyp='B  '
     if (sgname.eq.'Pmmm    ') lattyp='P  '
     if (sgname.eq.'Pnnn    ') lattyp='P  '
     if (sgname.eq.'Pccm    ') lattyp='P  '
     if (sgname.eq.'Pmaa    ') lattyp='P  '
     if (sgname.eq.'Pbmb    ') lattyp='P  '
     if (sgname.eq.'Pban    ') lattyp='P  '
     if (sgname.eq.'Pncb    ') lattyp='P  '
     if (sgname.eq.'Pcna    ') lattyp='P  '
     if (sgname.eq.'Pmma    ') lattyp='P  '
     if (sgname.eq.'Pbmm    ') lattyp='P  '
     if (sgname.eq.'Pmcm    ') lattyp='P  '
     if (sgname.eq.'Pmam    ') lattyp='P  '
     if (sgname.eq.'Pmmb    ') lattyp='P  '
     if (sgname.eq.'Pcmm    ') lattyp='P  '
     if (sgname.eq.'Pnna    ') lattyp='P  '
     if (sgname.eq.'Pbnn    ') lattyp='P  '
     if (sgname.eq.'Pncn    ') lattyp='P  '
     if (sgname.eq.'Pnan    ') lattyp='P  '
     if (sgname.eq.'Pnnb    ') lattyp='P  '
     if (sgname.eq.'Pcnn    ') lattyp='P  '
     if (sgname.eq.'Pmna    ') lattyp='P  '
     if (sgname.eq.'Pbmn    ') lattyp='P  '
     if (sgname.eq.'Pncm    ') lattyp='P  '
     if (sgname.eq.'Pman    ') lattyp='P  '
     if (sgname.eq.'Pnmb    ') lattyp='P  '
     if (sgname.eq.'Pcnm    ') lattyp='P  '
     if (sgname.eq.'Pcca    ') lattyp='P  '
     if (sgname.eq.'Pbaa    ') lattyp='P  '
     if (sgname.eq.'Pbcb    ') lattyp='P  '
     if (sgname.eq.'Pbab    ') lattyp='P  '
     if (sgname.eq.'Pccb    ') lattyp='P  '
     if (sgname.eq.'Pcaa    ') lattyp='P  '
     if (sgname.eq.'Pbam    ') lattyp='P  '
     if (sgname.eq.'Pmcb    ') lattyp='P  '
     if (sgname.eq.'Pcma    ') lattyp='P  '
     if (sgname.eq.'Pccn    ') lattyp='P  '
     if (sgname.eq.'Pnaa    ') lattyp='P  '
     if (sgname.eq.'Pbnb    ') lattyp='P  '
     if (sgname.eq.'Pbcm    ') lattyp='P  '
     if (sgname.eq.'Pmca    ') lattyp='P  '
     if (sgname.eq.'Pbma    ') lattyp='P  '
     if (sgname.eq.'Pcmb    ') lattyp='P  '
     if (sgname.eq.'Pcam    ') lattyp='P  '
     if (sgname.eq.'Pmab    ') lattyp='P  '
     if (sgname.eq.'Pnnm    ') lattyp='P  '
     if (sgname.eq.'Pmnn    ') lattyp='P  '
     if (sgname.eq.'Pnmn    ') lattyp='P  '
     if (sgname.eq.'Pmmn    ') lattyp='P  '
     if (sgname.eq.'Pnmm    ') lattyp='P  '
     if (sgname.eq.'Pmnm    ') lattyp='P  '
     if (sgname.eq.'Pbcn    ') lattyp='P  '
     if (sgname.eq.'Pnca    ') lattyp='P  '
     if (sgname.eq.'Pbna    ') lattyp='P  '
     if (sgname.eq.'Pcnb    ') lattyp='P  '
     if (sgname.eq.'Pcan    ') lattyp='P  '
     if (sgname.eq.'Pnab    ') lattyp='P  '
     if (sgname.eq.'Pbca    ') lattyp='P  '
     if (sgname.eq.'Pcab    ') lattyp='P  '
     if (sgname.eq.'Pnma    ') lattyp='P  '
     if (sgname.eq.'Pbnm    ') lattyp='P  '
     if (sgname.eq.'Pmcn    ') lattyp='P  '
     if (sgname.eq.'Pnam    ') lattyp='P  '
     if (sgname.eq.'Pmnb    ') lattyp='P  '
     if (sgname.eq.'Pcmn    ') lattyp='P  '
     if (sgname.eq.'Cmcm    ') lattyp='CXY'
     if (sgname.eq.'Amma    ') lattyp='CYZ'
     if (sgname.eq.'Bbmm    ') lattyp='CXZ'
     if (sgname.eq.'Bmmb    ') lattyp='CXZ'
     if (sgname.eq.'Ccmm    ') lattyp='CXY'
     if (sgname.eq.'Amam    ') lattyp='CYZ'
     if (sgname.eq.'Cmca    ') lattyp='CXY'
     if (sgname.eq.'Cmce    ') lattyp='CXY'
     if (sgname.eq.'Abma    ') lattyp='CYZ'
     if (sgname.eq.'Bbcm    ') lattyp='CXZ'
     if (sgname.eq.'Bmab    ') lattyp='CXZ'
     if (sgname.eq.'Ccmb    ') lattyp='CXY'
     if (sgname.eq.'Acam    ') lattyp='CYZ'
     if (sgname.eq.'Cmmm    ') lattyp='CXY'
     if (sgname.eq.'Ammm    ') lattyp='CYZ'
     if (sgname.eq.'Bmmm    ') lattyp='CXZ'
     if (sgname.eq.'Cccm    ') lattyp='CXY'
     if (sgname.eq.'Amaa    ') lattyp='CYZ'
     if (sgname.eq.'Bbmb    ') lattyp='CXZ'
     if (sgname.eq.'Cmma    ') lattyp='CXY'
     if (sgname.eq.'Abmm    ') lattyp='CYZ'
     if (sgname.eq.'Bmcm    ') lattyp='CXZ'
     if (sgname.eq.'Bmam    ') lattyp='CXZ'
     if (sgname.eq.'Cmmb    ') lattyp='CXY'
     if (sgname.eq.'Acmm    ') lattyp='CYZ'
     if (sgname.eq.'Ccca    ') lattyp='CXY'
     if (sgname.eq.'Abaa    ') lattyp='CYZ'
     if (sgname.eq.'Bbcb    ') lattyp='CXZ'
     if (sgname.eq.'Bbab    ') lattyp='CXZ'
     if (sgname.eq.'Cccb    ') lattyp='CXY'
     if (sgname.eq.'Acaa    ') lattyp='CYZ'
     if (sgname.eq.'Fmmm    ') lattyp='F  '
     if (sgname.eq.'Fddd    ') lattyp='F  '
     if (sgname.eq.'Immm    ') lattyp='B  '
     if (sgname.eq.'Ibam    ') lattyp='B  '
     if (sgname.eq.'Imcb    ') lattyp='B  '
     if (sgname.eq.'Icma    ') lattyp='B  '
     if (sgname.eq.'Ibca    ') lattyp='B  '
     if (sgname.eq.'Icab    ') lattyp='B  '
     if (sgname.eq.'Imma    ') lattyp='B  '
     if (sgname.eq.'Ibmm    ') lattyp='B  '
     if (sgname.eq.'Imcm    ') lattyp='B  '
     if (sgname.eq.'Imam    ') lattyp='B  '
     if (sgname.eq.'Immb    ') lattyp='B  '
     if (sgname.eq.'Icmm    ') lattyp='B  '
     if (sgname.eq.'P4      ') lattyp='P  '
     if (sgname.eq.'P41     ') lattyp='P  '
     if (sgname.eq.'P42     ') lattyp='P  '
     if (sgname.eq.'P43     ') lattyp='P  '
     if (sgname.eq.'I4      ') lattyp='B  '
     if (sgname.eq.'I41     ') lattyp='B  '
     if (sgname.eq.'P-4     ') lattyp='P  '
     if (sgname.eq.'I-4     ') lattyp='B  '
     if (sgname.eq.'P4/m    ') lattyp='P  '
     if (sgname.eq.'P42/m   ') lattyp='P  '
     if (sgname.eq.'P4/n    ') lattyp='P  '
     if (sgname.eq.'P42/n   ') lattyp='P  '
     if (sgname.eq.'I4/m    ') lattyp='B  '
     if (sgname.eq.'I41/a   ') lattyp='B  '
     if (sgname.eq.'P422    ') lattyp='P  '
     if (sgname.eq.'P4212   ') lattyp='P  '
     if (sgname.eq.'P4122   ') lattyp='P  '
     if (sgname.eq.'P41212  ') lattyp='P  '
     if (sgname.eq.'P4222   ') lattyp='P  '
     if (sgname.eq.'P42212  ') lattyp='P  '
     if (sgname.eq.'P4322   ') lattyp='P  '
     if (sgname.eq.'P43212  ') lattyp='P  '
     if (sgname.eq.'I422    ') lattyp='B  '
     if (sgname.eq.'I4122   ') lattyp='B  '
     if (sgname.eq.'P4mm    ') lattyp='P  '
     if (sgname.eq.'P4bm    ') lattyp='P  '
     if (sgname.eq.'P42cm   ') lattyp='P  '
     if (sgname.eq.'P42nm   ') lattyp='P  '
     if (sgname.eq.'P4cc    ') lattyp='P  '
     if (sgname.eq.'P4nc    ') lattyp='P  '
     if (sgname.eq.'P42mc   ') lattyp='P  '
     if (sgname.eq.'P42bc   ') lattyp='P  '
     if (sgname.eq.'I4mm    ') lattyp='B  '
     if (sgname.eq.'I4cm    ') lattyp='B  '
     if (sgname.eq.'I41md   ') lattyp='B  '
     if (sgname.eq.'I41cd   ') lattyp='B  '
     if (sgname.eq.'P-42m   ') lattyp='P  '
     if (sgname.eq.'P-42c   ') lattyp='P  '
     if (sgname.eq.'P-421m  ') lattyp='P  '
     if (sgname.eq.'P-421c  ') lattyp='P  '
     if (sgname.eq.'P-4m2   ') lattyp='P  '
     if (sgname.eq.'P-4c2   ') lattyp='P  '
     if (sgname.eq.'P-4b2   ') lattyp='P  '
     if (sgname.eq.'P-4n2   ') lattyp='P  '
     if (sgname.eq.'I-4m2   ') lattyp='B  '
     if (sgname.eq.'I-4c2   ') lattyp='B  '
     if (sgname.eq.'I-42m   ') lattyp='B  '
     if (sgname.eq.'I-42d   ') lattyp='B  '
     if (sgname.eq.'P4/mmm  ') lattyp='P  '
     if (sgname.eq.'P4/mcc  ') lattyp='P  '
     if (sgname.eq.'P4/nbm  ') lattyp='P  '
     if (sgname.eq.'P4/nnc  ') lattyp='P  '
     if (sgname.eq.'P4/mbm  ') lattyp='P  '
     if (sgname.eq.'P4/mnc  ') lattyp='P  '
     if (sgname.eq.'P4/nmm  ') lattyp='P  '
     if (sgname.eq.'P4/ncc  ') lattyp='P  '
     if (sgname.eq.'P42/mmc ') lattyp='P  '
     if (sgname.eq.'P42/mcm ') lattyp='P  '
     if (sgname.eq.'P42/nbc ') lattyp='P  '
     if (sgname.eq.'P42/nnm ') lattyp='P  '
     if (sgname.eq.'P42/mbc ') lattyp='P  '
     if (sgname.eq.'P42/mnm ') lattyp='P  '
     if (sgname.eq.'P42/nmc ') lattyp='P  '
     if (sgname.eq.'P42/ncm ') lattyp='P  '
     if (sgname.eq.'I4/mmm  ') lattyp='B  '
     if (sgname.eq.'I4/mcm  ') lattyp='B  '
     if (sgname.eq.'I41/amd ') lattyp='B  '
     if (sgname.eq.'I41/acd ') lattyp='B  '
     if (sgname.eq.'P3      ') lattyp='H  '
     if (sgname.eq.'P31     ') lattyp='H  '
     if (sgname.eq.'P32     ') lattyp='H  '
     if (sgname.eq.'R3      ') lattyp='R  '
     if (sgname.eq.'P-3     ') lattyp='H  '
     if (sgname.eq.'R-3     ') lattyp='R  '
     if (sgname.eq.'P312    ') lattyp='H  '
     if (sgname.eq.'P321    ') lattyp='H  '
     if (sgname.eq.'P3112   ') lattyp='H  '
     if (sgname.eq.'P3121   ') lattyp='H  '
     if (sgname.eq.'P3212   ') lattyp='H  '
     if (sgname.eq.'P3221   ') lattyp='H  '
     if (sgname.eq.'R32     ') lattyp='R  '
     if (sgname.eq.'P3m1    ') lattyp='H  '
     if (sgname.eq.'P31m    ') lattyp='H  '
     if (sgname.eq.'P3c1    ') lattyp='H  '
     if (sgname.eq.'P31c    ') lattyp='H  '
     if (sgname.eq.'R3m     ') lattyp='R  '
     if (sgname.eq.'R3c     ') lattyp='R  '
     if (sgname.eq.'P-31m   ') lattyp='H  '
     if (sgname.eq.'P-31c   ') lattyp='H  '
     if (sgname.eq.'P-3m1   ') lattyp='H  '
     if (sgname.eq.'P-3c1   ') lattyp='H  '
     if (sgname.eq.'R-3m    ') lattyp='R  '
     if (sgname.eq.'R-3c    ') lattyp='R  '
     if (sgname.eq.'P6      ') lattyp='H  '
     if (sgname.eq.'P61     ') lattyp='H  '
     if (sgname.eq.'P65     ') lattyp='H  '
     if (sgname.eq.'P62     ') lattyp='H  '
     if (sgname.eq.'P64     ') lattyp='H  '
     if (sgname.eq.'P63     ') lattyp='H  '
     if (sgname.eq.'P-6     ') lattyp='H  '
     if (sgname.eq.'P6/m    ') lattyp='H  '
     if (sgname.eq.'P63/m   ') lattyp='H  '
     if (sgname.eq.'P622    ') lattyp='H  '
     if (sgname.eq.'P6122   ') lattyp='H  '
     if (sgname.eq.'P6522   ') lattyp='H  '
     if (sgname.eq.'P6222   ') lattyp='H  '
     if (sgname.eq.'P6422   ') lattyp='H  '
     if (sgname.eq.'P6322   ') lattyp='H  '
     if (sgname.eq.'P6mm    ') lattyp='H  '
     if (sgname.eq.'P6cc    ') lattyp='H  '
     if (sgname.eq.'P63cm   ') lattyp='H  '
     if (sgname.eq.'P63mc   ') lattyp='H  '
     if (sgname.eq.'P-6m2   ') lattyp='H  '
     if (sgname.eq.'P-6c2   ') lattyp='H  '
     if (sgname.eq.'P-62m   ') lattyp='H  '
     if (sgname.eq.'P-62c   ') lattyp='H  '
     if (sgname.eq.'P6/mmm  ') lattyp='H  '
     if (sgname.eq.'P6/mcc  ') lattyp='H  '
     if (sgname.eq.'P63/mcm ') lattyp='H  '
     if (sgname.eq.'P63/mmc ') lattyp='H  '
     if (sgname.eq.'P23     ') lattyp='P  '
     if (sgname.eq.'F23     ') lattyp='F  '
     if (sgname.eq.'I23     ') lattyp='B  '
     if (sgname.eq.'P213    ') lattyp='P  '
     if (sgname.eq.'I213    ') lattyp='B  '
     if (sgname.eq.'Pm-3    ') lattyp='P  '
     if (sgname.eq.'Pn-3    ') lattyp='P  '
     if (sgname.eq.'Fm-3    ') lattyp='F  '
     if (sgname.eq.'Fd-3    ') lattyp='F  '
     if (sgname.eq.'Im-3    ') lattyp='B  '
     if (sgname.eq.'Pa-3    ') lattyp='P  '
     if (sgname.eq.'Ia-3    ') lattyp='B  '
     if (sgname.eq.'P432    ') lattyp='P  '
     if (sgname.eq.'P4232   ') lattyp='P  '
     if (sgname.eq.'F432    ') lattyp='F  '
     if (sgname.eq.'F4132   ') lattyp='F  '
     if (sgname.eq.'I432    ') lattyp='B  '
     if (sgname.eq.'P4332   ') lattyp='P  '
     if (sgname.eq.'P4132   ') lattyp='P  '
     if (sgname.eq.'I4132   ') lattyp='B  '
     if (sgname.eq.'P-43m   ') lattyp='P  '
     if (sgname.eq.'F-43m   ') lattyp='F  '
     if (sgname.eq.'I-43m   ') lattyp='B  '
     if (sgname.eq.'P-43n   ') lattyp='P  '
     if (sgname.eq.'F-43c   ') lattyp='F  '
     if (sgname.eq.'I-43d   ') lattyp='B  '
     if (sgname.eq.'Pm-3m   ') lattyp='P  '
     if (sgname.eq.'Pn-3n   ') lattyp='P  '
     if (sgname.eq.'Pm-3n   ') lattyp='P  '
     if (sgname.eq.'Pn-3m   ') lattyp='P  '
     if (sgname.eq.'Fm-3m   ') lattyp='F  '
     if (sgname.eq.'Fm-3c   ') lattyp='F  '
     if (sgname.eq.'Fd-3m   ') lattyp='F  '
     if (sgname.eq.'Fd-3c   ') lattyp='F  '
     if (sgname.eq.'Im-3m   ') lattyp='B  '
     if (sgname.eq.'Ia-3d   ') lattyp='B  '


   end subroutine getlattype
