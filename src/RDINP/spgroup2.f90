 subroutine getsgnum(sgname,sgnum)
! This routine translates between two different naming conventions for spacegroups.
implicit none
character*8,intent(in) :: sgname
integer,intent(out) :: sgnum

     if (sgname.eq.'P1      ')  sgnum = 1   
     if (sgname.eq.'P-1     ')  sgnum = 2   
     if (sgname.eq.'P2      ')  sgnum = 3   
     if (sgname.eq.'P2      ')  sgnum = 3   
     if (sgname.eq.'P2      ')  sgnum = 3   
     if (sgname.eq.'P21     ')  sgnum = 4   
     if (sgname.eq.'P21     ')  sgnum = 4   
     if (sgname.eq.'P21     ')  sgnum = 4   
     if (sgname.eq.'C2      ')  sgnum = 5   
     if (sgname.eq.'A2      ')  sgnum = 5   
     if (sgname.eq.'B2      ')  sgnum = 5   
     if (sgname.eq.'B2      ')  sgnum = 5   
     if (sgname.eq.'C2      ')  sgnum = 5   
     if (sgname.eq.'A2      ')  sgnum = 5   
     if (sgname.eq.'Pm      ')  sgnum = 6   
     if (sgname.eq.'Pm      ')  sgnum = 6   
     if (sgname.eq.'Pm      ')  sgnum = 6   
     if (sgname.eq.'Pc      ')  sgnum = 7   
     if (sgname.eq.'Pa      ')  sgnum = 7   
     if (sgname.eq.'Pb      ')  sgnum = 7   
     if (sgname.eq.'Pb      ')  sgnum = 7   
     if (sgname.eq.'Pc      ')  sgnum = 7   
     if (sgname.eq.'Pa      ')  sgnum = 7   
     if (sgname.eq.'Pn      ')  sgnum = 7   
     if (sgname.eq.'Pn      ')  sgnum = 7   
     if (sgname.eq.'Pn      ')  sgnum = 7   
     if (sgname.eq.'Cm      ')  sgnum = 8   
     if (sgname.eq.'Am      ')  sgnum = 8   
     if (sgname.eq.'Bm      ')  sgnum = 8   
     if (sgname.eq.'Bm      ')  sgnum = 8   
     if (sgname.eq.'Cm      ')  sgnum = 8   
     if (sgname.eq.'Am      ')  sgnum = 8   
     if (sgname.eq.'Cc      ')  sgnum = 9   
     if (sgname.eq.'Aa      ')  sgnum = 9   
     if (sgname.eq.'Bb      ')  sgnum = 9   
     if (sgname.eq.'Bb      ')  sgnum = 9   
     if (sgname.eq.'Cc      ')  sgnum = 9   
     if (sgname.eq.'Aa      ')  sgnum = 9   
     if (sgname.eq.'P2/m    ')  sgnum = 10  
     if (sgname.eq.'P2/m    ')  sgnum = 10  
     if (sgname.eq.'P2/m    ')  sgnum = 10  
     if (sgname.eq.'P21/m   ')  sgnum = 11  
     if (sgname.eq.'P21/m   ')  sgnum = 11  
     if (sgname.eq.'P21/m   ')  sgnum = 11  
     if (sgname.eq.'C2/m    ')  sgnum = 12  
     if (sgname.eq.'A2/m    ')  sgnum = 12  
     if (sgname.eq.'B2/m    ')  sgnum = 12  
     if (sgname.eq.'B2/m    ')  sgnum = 12  
     if (sgname.eq.'C2/m    ')  sgnum = 12  
     if (sgname.eq.'A2/m    ')  sgnum = 12  
     if (sgname.eq.'P2/c    ')  sgnum = 13  
     if (sgname.eq.'P2/a    ')  sgnum = 13  
     if (sgname.eq.'P2/b    ')  sgnum = 13  
     if (sgname.eq.'P2/b    ')  sgnum = 13  
     if (sgname.eq.'P2/c    ')  sgnum = 13  
     if (sgname.eq.'P2/a    ')  sgnum = 13  
     if (sgname.eq.'P2/n    ')  sgnum = 13  
     if (sgname.eq.'P2/n    ')  sgnum = 13  
     if (sgname.eq.'P2/n    ')  sgnum = 13  
     if (sgname.eq.'P21/c   ')  sgnum = 14  
     if (sgname.eq.'P21/a   ')  sgnum = 14  
     if (sgname.eq.'P21/b   ')  sgnum = 14  
     if (sgname.eq.'P21/b   ')  sgnum = 14  
     if (sgname.eq.'P21/c   ')  sgnum = 14  
     if (sgname.eq.'P21/a   ')  sgnum = 14  
     if (sgname.eq.'P21/n   ')  sgnum = 14  
     if (sgname.eq.'P21/n   ')  sgnum = 14  
     if (sgname.eq.'P21/n   ')  sgnum = 14  
     if (sgname.eq.'C2/c    ')  sgnum = 15  
     if (sgname.eq.'A2/a    ')  sgnum = 15  
     if (sgname.eq.'B2/b    ')  sgnum = 15  
     if (sgname.eq.'B2/b    ')  sgnum = 15  
     if (sgname.eq.'C2/c    ')  sgnum = 15  
     if (sgname.eq.'A2/a    ')  sgnum = 15  
     if (sgname.eq.'P222    ')  sgnum = 16  
     if (sgname.eq.'P2221   ')  sgnum = 17  
     if (sgname.eq.'P2122   ')  sgnum = 17  
     if (sgname.eq.'P2212   ')  sgnum = 17  
     if (sgname.eq.'P21212  ')  sgnum = 18  
     if (sgname.eq.'P22121  ')  sgnum = 18  
     if (sgname.eq.'P21221  ')  sgnum = 18  
     if (sgname.eq.'P212121 ')  sgnum = 19  
     if (sgname.eq.'C2221   ')  sgnum = 20  
     if (sgname.eq.'A2122   ')  sgnum = 20  
     if (sgname.eq.'B2212   ')  sgnum = 20  
     if (sgname.eq.'C222    ')  sgnum = 21  
     if (sgname.eq.'A222    ')  sgnum = 21  
     if (sgname.eq.'B222    ')  sgnum = 21  
     if (sgname.eq.'F222    ')  sgnum = 22  
     if (sgname.eq.'I222    ')  sgnum = 23  
     if (sgname.eq.'I212121 ')  sgnum = 24  
     if (sgname.eq.'Pmm2    ')  sgnum = 25  
     if (sgname.eq.'P2mm    ')  sgnum = 25  
     if (sgname.eq.'Pm2m    ')  sgnum = 25  
     if (sgname.eq.'Pmc21   ')  sgnum = 26  
     if (sgname.eq.'P21ma   ')  sgnum = 26  
     if (sgname.eq.'Pb21m   ')  sgnum = 26  
     if (sgname.eq.'Pcm21   ')  sgnum = 26  
     if (sgname.eq.'P21am   ')  sgnum = 26  
     if (sgname.eq.'Pm21b   ')  sgnum = 26  
     if (sgname.eq.'Pcc2    ')  sgnum = 27  
     if (sgname.eq.'P2aa    ')  sgnum = 27  
     if (sgname.eq.'Pb2b    ')  sgnum = 27  
     if (sgname.eq.'Pma2    ')  sgnum = 28  
     if (sgname.eq.'P2mb    ')  sgnum = 28  
     if (sgname.eq.'Pc2m    ')  sgnum = 28  
     if (sgname.eq.'Pbm2    ')  sgnum = 28  
     if (sgname.eq.'P2cm    ')  sgnum = 28  
     if (sgname.eq.'Pm2a    ')  sgnum = 28  
     if (sgname.eq.'Pca21   ')  sgnum = 29  
     if (sgname.eq.'P21ab   ')  sgnum = 29  
     if (sgname.eq.'Pc21b   ')  sgnum = 29  
     if (sgname.eq.'Pbc21   ')  sgnum = 29  
     if (sgname.eq.'P21ca   ')  sgnum = 29  
     if (sgname.eq.'Pb21a   ')  sgnum = 29  
     if (sgname.eq.'Pnc2    ')  sgnum = 30  
     if (sgname.eq.'P2na    ')  sgnum = 30  
     if (sgname.eq.'Pb2n    ')  sgnum = 30  
     if (sgname.eq.'Pcn2    ')  sgnum = 30  
     if (sgname.eq.'P2an    ')  sgnum = 30  
     if (sgname.eq.'Pn2b    ')  sgnum = 30  
     if (sgname.eq.'Pmn21   ')  sgnum = 31  
     if (sgname.eq.'P21mn   ')  sgnum = 31  
     if (sgname.eq.'Pn21m   ')  sgnum = 31  
     if (sgname.eq.'Pnm21   ')  sgnum = 31  
     if (sgname.eq.'P21nm   ')  sgnum = 31  
     if (sgname.eq.'Pm21n   ')  sgnum = 31  
     if (sgname.eq.'Pba2    ')  sgnum = 32  
     if (sgname.eq.'P2cb    ')  sgnum = 32  
     if (sgname.eq.'Pc2a    ')  sgnum = 32  
     if (sgname.eq.'Pna21   ')  sgnum = 33  
     if (sgname.eq.'P21nb   ')  sgnum = 33  
     if (sgname.eq.'Pc21n   ')  sgnum = 33  
     if (sgname.eq.'Pbn21   ')  sgnum = 33  
     if (sgname.eq.'P21cn   ')  sgnum = 33  
     if (sgname.eq.'Pn21a   ')  sgnum = 33  
     if (sgname.eq.'Pnn2    ')  sgnum = 34  
     if (sgname.eq.'P2nn    ')  sgnum = 34  
     if (sgname.eq.'Pn2n    ')  sgnum = 34  
     if (sgname.eq.'Cmm2    ')  sgnum = 35  
     if (sgname.eq.'A2mm    ')  sgnum = 35  
     if (sgname.eq.'Bm2m    ')  sgnum = 35  
     if (sgname.eq.'Cmc21   ')  sgnum = 36  
     if (sgname.eq.'A21ma   ')  sgnum = 36  
     if (sgname.eq.'Bb21m   ')  sgnum = 36  
     if (sgname.eq.'Ccm21   ')  sgnum = 36  
     if (sgname.eq.'A21am   ')  sgnum = 36  
     if (sgname.eq.'Bm21b   ')  sgnum = 36  
     if (sgname.eq.'Ccc2    ')  sgnum = 37  
     if (sgname.eq.'A2aa    ')  sgnum = 37  
     if (sgname.eq.'Bb2b    ')  sgnum = 37  
     if (sgname.eq.'Amm2    ')  sgnum = 38  
     if (sgname.eq.'B2mm    ')  sgnum = 38  
     if (sgname.eq.'Cm2m    ')  sgnum = 38  
     if (sgname.eq.'Abm2    ')  sgnum = 39  
     if (sgname.eq.'B2cm    ')  sgnum = 39  
     if (sgname.eq.'Cm2a    ')  sgnum = 39  
     if (sgname.eq.'Bma2    ')  sgnum = 39  
     if (sgname.eq.'C2mb    ')  sgnum = 39  
     if (sgname.eq.'Ac2m    ')  sgnum = 39  
     if (sgname.eq.'Ama2    ')  sgnum = 40  
     if (sgname.eq.'B2mb    ')  sgnum = 40  
     if (sgname.eq.'Cc2m    ')  sgnum = 40  
     if (sgname.eq.'Bbm2    ')  sgnum = 40  
     if (sgname.eq.'C2cm    ')  sgnum = 40  
     if (sgname.eq.'Am2a    ')  sgnum = 40  
     if (sgname.eq.'Aba2    ')  sgnum = 41  
     if (sgname.eq.'B2cb    ')  sgnum = 41  
     if (sgname.eq.'Cc2a    ')  sgnum = 41  
     if (sgname.eq.'Bba2    ')  sgnum = 41  
     if (sgname.eq.'C2cb    ')  sgnum = 41  
     if (sgname.eq.'Ac2a    ')  sgnum = 41  
     if (sgname.eq.'Fmm2    ')  sgnum = 42  
     if (sgname.eq.'F2mm    ')  sgnum = 42  
     if (sgname.eq.'Fm2m    ')  sgnum = 42  
     if (sgname.eq.'Fdd2    ')  sgnum = 43  
     if (sgname.eq.'F2dd    ')  sgnum = 43  
     if (sgname.eq.'Fd2d    ')  sgnum = 43  
     if (sgname.eq.'Imm2    ')  sgnum = 44  
     if (sgname.eq.'I2mm    ')  sgnum = 44  
     if (sgname.eq.'Im2m    ')  sgnum = 44  
     if (sgname.eq.'Iba2    ')  sgnum = 45  
     if (sgname.eq.'I2cb    ')  sgnum = 45  
     if (sgname.eq.'Ic2a    ')  sgnum = 45  
     if (sgname.eq.'Ima2    ')  sgnum = 46  
     if (sgname.eq.'I2mb    ')  sgnum = 46  
     if (sgname.eq.'Ic2m    ')  sgnum = 46  
     if (sgname.eq.'Ibm2    ')  sgnum = 46  
     if (sgname.eq.'I2cm    ')  sgnum = 46  
     if (sgname.eq.'Im2a    ')  sgnum = 46  
     if (sgname.eq.'Pmmm    ')  sgnum = 47  
     if (sgname.eq.'Pnnn    ')  sgnum = 48  
     if (sgname.eq.'Pccm    ')  sgnum = 49  
     if (sgname.eq.'Pmaa    ')  sgnum = 49  
     if (sgname.eq.'Pbmb    ')  sgnum = 49  
     if (sgname.eq.'Pban    ')  sgnum = 50  
     if (sgname.eq.'Pncb    ')  sgnum = 50  
     if (sgname.eq.'Pcna    ')  sgnum = 50  
     if (sgname.eq.'Pmma    ')  sgnum = 51  
     if (sgname.eq.'Pbmm    ')  sgnum = 51  
     if (sgname.eq.'Pmcm    ')  sgnum = 51  
     if (sgname.eq.'Pmam    ')  sgnum = 51  
     if (sgname.eq.'Pmmb    ')  sgnum = 51  
     if (sgname.eq.'Pcmm    ')  sgnum = 51  
     if (sgname.eq.'Pnna    ')  sgnum = 52  
     if (sgname.eq.'Pbnn    ')  sgnum = 52  
     if (sgname.eq.'Pncn    ')  sgnum = 52  
     if (sgname.eq.'Pnan    ')  sgnum = 52  
     if (sgname.eq.'Pnnb    ')  sgnum = 52  
     if (sgname.eq.'Pcnn    ')  sgnum = 52  
     if (sgname.eq.'Pmna    ')  sgnum = 53  
     if (sgname.eq.'Pbmn    ')  sgnum = 53  
     if (sgname.eq.'Pncm    ')  sgnum = 53  
     if (sgname.eq.'Pman    ')  sgnum = 53  
     if (sgname.eq.'Pnmb    ')  sgnum = 53  
     if (sgname.eq.'Pcnm    ')  sgnum = 53  
     if (sgname.eq.'Pcca    ')  sgnum = 54  
     if (sgname.eq.'Pbaa    ')  sgnum = 54  
     if (sgname.eq.'Pbcb    ')  sgnum = 54  
     if (sgname.eq.'Pbab    ')  sgnum = 54  
     if (sgname.eq.'Pccb    ')  sgnum = 54  
     if (sgname.eq.'Pcaa    ')  sgnum = 54  
     if (sgname.eq.'Pbam    ')  sgnum = 55  
     if (sgname.eq.'Pmcb    ')  sgnum = 55  
     if (sgname.eq.'Pcma    ')  sgnum = 55  
     if (sgname.eq.'Pccn    ')  sgnum = 56  
     if (sgname.eq.'Pnaa    ')  sgnum = 56  
     if (sgname.eq.'Pbnb    ')  sgnum = 56  
     if (sgname.eq.'Pbcm    ')  sgnum = 57  
     if (sgname.eq.'Pmca    ')  sgnum = 57  
     if (sgname.eq.'Pbma    ')  sgnum = 57  
     if (sgname.eq.'Pcmb    ')  sgnum = 57  
     if (sgname.eq.'Pcam    ')  sgnum = 57  
     if (sgname.eq.'Pmab    ')  sgnum = 57  
     if (sgname.eq.'Pnnm    ')  sgnum = 58  
     if (sgname.eq.'Pmnn    ')  sgnum = 58  
     if (sgname.eq.'Pnmn    ')  sgnum = 58  
     if (sgname.eq.'Pmmn    ')  sgnum = 59  
     if (sgname.eq.'Pnmm    ')  sgnum = 59  
     if (sgname.eq.'Pmnm    ')  sgnum = 59  
     if (sgname.eq.'Pbcn    ')  sgnum = 60  
     if (sgname.eq.'Pnca    ')  sgnum = 60  
     if (sgname.eq.'Pbna    ')  sgnum = 60  
     if (sgname.eq.'Pcnb    ')  sgnum = 60  
     if (sgname.eq.'Pcan    ')  sgnum = 60  
     if (sgname.eq.'Pnab    ')  sgnum = 60  
     if (sgname.eq.'Pbca    ')  sgnum = 61  
     if (sgname.eq.'Pcab    ')  sgnum = 61  
     if (sgname.eq.'Pnma    ')  sgnum = 62  
     if (sgname.eq.'Pbnm    ')  sgnum = 62  
     if (sgname.eq.'Pmcn    ')  sgnum = 62  
     if (sgname.eq.'Pnam    ')  sgnum = 62  
     if (sgname.eq.'Pmnb    ')  sgnum = 62  
     if (sgname.eq.'Pcmn    ')  sgnum = 62  
     if (sgname.eq.'Cmcm    ')  sgnum = 63  
     if (sgname.eq.'Amma    ')  sgnum = 63  
     if (sgname.eq.'Bbmm    ')  sgnum = 63  
     if (sgname.eq.'Bmmb    ')  sgnum = 63  
     if (sgname.eq.'Ccmm    ')  sgnum = 63  
     if (sgname.eq.'Amam    ')  sgnum = 63  
     if (sgname.eq.'Cmca    ')  sgnum = 64  
     if (sgname.eq.'Cmce    ')  sgnum = 64  
     if (sgname.eq.'Abma    ')  sgnum = 64  
     if (sgname.eq.'Bbcm    ')  sgnum = 64  
     if (sgname.eq.'Bmab    ')  sgnum = 64  
     if (sgname.eq.'Ccmb    ')  sgnum = 64  
     if (sgname.eq.'Acam    ')  sgnum = 64  
     if (sgname.eq.'Cmmm    ')  sgnum = 65  
     if (sgname.eq.'Ammm    ')  sgnum = 65  
     if (sgname.eq.'Bmmm    ')  sgnum = 65  
     if (sgname.eq.'Cccm    ')  sgnum = 66  
     if (sgname.eq.'Amaa    ')  sgnum = 66  
     if (sgname.eq.'Bbmb    ')  sgnum = 66  
     if (sgname.eq.'Cmma    ')  sgnum = 67  
     if (sgname.eq.'Abmm    ')  sgnum = 67  
     if (sgname.eq.'Bmcm    ')  sgnum = 67  
     if (sgname.eq.'Bmam    ')  sgnum = 67  
     if (sgname.eq.'Cmmb    ')  sgnum = 67  
     if (sgname.eq.'Acmm    ')  sgnum = 67  
     if (sgname.eq.'Ccca    ')  sgnum = 68  
     if (sgname.eq.'Abaa    ')  sgnum = 68  
     if (sgname.eq.'Bbcb    ')  sgnum = 68  
     if (sgname.eq.'Bbab    ')  sgnum = 68  
     if (sgname.eq.'Cccb    ')  sgnum = 68  
     if (sgname.eq.'Acaa    ')  sgnum = 68  
     if (sgname.eq.'Fmmm    ')  sgnum = 69  
     if (sgname.eq.'Fddd    ')  sgnum = 70  
     if (sgname.eq.'Immm    ')  sgnum = 71  
     if (sgname.eq.'Ibam    ')  sgnum = 72  
     if (sgname.eq.'Imcb    ')  sgnum = 72  
     if (sgname.eq.'Icma    ')  sgnum = 72  
     if (sgname.eq.'Ibca    ')  sgnum = 73  
     if (sgname.eq.'Icab    ')  sgnum = 73  
     if (sgname.eq.'Imma    ')  sgnum = 74  
     if (sgname.eq.'Ibmm    ')  sgnum = 74  
     if (sgname.eq.'Imcm    ')  sgnum = 74  
     if (sgname.eq.'Imam    ')  sgnum = 74  
     if (sgname.eq.'Immb    ')  sgnum = 74  
     if (sgname.eq.'Icmm    ')  sgnum = 74  
     if (sgname.eq.'P4      ')  sgnum = 75  
     if (sgname.eq.'P41     ')  sgnum = 76  
     if (sgname.eq.'P42     ')  sgnum = 77  
     if (sgname.eq.'P43     ')  sgnum = 78  
     if (sgname.eq.'I4      ')  sgnum = 79  
     if (sgname.eq.'I41     ')  sgnum = 80  
     if (sgname.eq.'P-4     ')  sgnum = 81  
     if (sgname.eq.'I-4     ')  sgnum = 82  
     if (sgname.eq.'P4/m    ')  sgnum = 83  
     if (sgname.eq.'P42/m   ')  sgnum = 84  
     if (sgname.eq.'P4/n    ')  sgnum = 85  
     if (sgname.eq.'P42/n   ')  sgnum = 86  
     if (sgname.eq.'I4/m    ')  sgnum = 87  
     if (sgname.eq.'I41/a   ')  sgnum = 88  
     if (sgname.eq.'P422    ')  sgnum = 89  
     if (sgname.eq.'P4212   ')  sgnum = 90  
     if (sgname.eq.'P4122   ')  sgnum = 91  
     if (sgname.eq.'P41212  ')  sgnum = 92  
     if (sgname.eq.'P4222   ')  sgnum = 93  
     if (sgname.eq.'P42212  ')  sgnum = 94  
     if (sgname.eq.'P4322   ')  sgnum = 95  
     if (sgname.eq.'P43212  ')  sgnum = 96  
     if (sgname.eq.'I422    ')  sgnum = 97  
     if (sgname.eq.'I4122   ')  sgnum = 98  
     if (sgname.eq.'P4mm    ')  sgnum = 99  
     if (sgname.eq.'P4bm    ')  sgnum = 100 
     if (sgname.eq.'P42cm   ')  sgnum = 101 
     if (sgname.eq.'P42nm   ')  sgnum = 102 
     if (sgname.eq.'P4cc    ')  sgnum = 103 
     if (sgname.eq.'P4nc    ')  sgnum = 104 
     if (sgname.eq.'P42mc   ')  sgnum = 105 
     if (sgname.eq.'P42bc   ')  sgnum = 106 
     if (sgname.eq.'I4mm    ')  sgnum = 107 
     if (sgname.eq.'I4cm    ')  sgnum = 108 
     if (sgname.eq.'I41md   ')  sgnum = 109 
     if (sgname.eq.'I41cd   ')  sgnum = 110 
     if (sgname.eq.'P-42m   ')  sgnum = 111 
     if (sgname.eq.'P-42c   ')  sgnum = 112 
     if (sgname.eq.'P-421m  ')  sgnum = 113 
     if (sgname.eq.'P-421c  ')  sgnum = 114 
     if (sgname.eq.'P-4m2   ')  sgnum = 115 
     if (sgname.eq.'P-4c2   ')  sgnum = 116 
     if (sgname.eq.'P-4b2   ')  sgnum = 117 
     if (sgname.eq.'P-4n2   ')  sgnum = 118 
     if (sgname.eq.'I-4m2   ')  sgnum = 119 
     if (sgname.eq.'I-4c2   ')  sgnum = 120 
     if (sgname.eq.'I-42m   ')  sgnum = 121 
     if (sgname.eq.'I-42d   ')  sgnum = 122 
     if (sgname.eq.'P4/mmm  ')  sgnum = 123 
     if (sgname.eq.'P4/mcc  ')  sgnum = 124 
     if (sgname.eq.'P4/nbm  ')  sgnum = 125 
     if (sgname.eq.'P4/nnc  ')  sgnum = 126 
     if (sgname.eq.'P4/mbm  ')  sgnum = 127 
     if (sgname.eq.'P4/mnc  ')  sgnum = 128 
     if (sgname.eq.'P4/nmm  ')  sgnum = 129 
     if (sgname.eq.'P4/ncc  ')  sgnum = 130 
     if (sgname.eq.'P42/mmc ')  sgnum = 131 
     if (sgname.eq.'P42/mcm ')  sgnum = 132 
     if (sgname.eq.'P42/nbc ')  sgnum = 133 
     if (sgname.eq.'P42/nnm ')  sgnum = 134 
     if (sgname.eq.'P42/mbc ')  sgnum = 135 
     if (sgname.eq.'P42/mnm ')  sgnum = 136 
     if (sgname.eq.'P42/nmc ')  sgnum = 137 
     if (sgname.eq.'P42/ncm ')  sgnum = 138 
     if (sgname.eq.'I4/mmm  ')  sgnum = 139 
     if (sgname.eq.'I4/mcm  ')  sgnum = 140 
     if (sgname.eq.'I41/amd ')  sgnum = 141 
     if (sgname.eq.'I41/acd ')  sgnum = 142 
     if (sgname.eq.'P3      ')  sgnum = 143 
     if (sgname.eq.'P31     ')  sgnum = 144 
     if (sgname.eq.'P32     ')  sgnum = 145 
     if (sgname.eq.'R3      ')  sgnum = 146 
     if (sgname.eq.'P-3     ')  sgnum = 147 
     if (sgname.eq.'R-3     ')  sgnum = 148 
     if (sgname.eq.'P312    ')  sgnum = 149 
     if (sgname.eq.'P321    ')  sgnum = 150 
     if (sgname.eq.'P3112   ')  sgnum = 151 
     if (sgname.eq.'P3121   ')  sgnum = 152 
     if (sgname.eq.'P3212   ')  sgnum = 153 
     if (sgname.eq.'P3221   ')  sgnum = 154 
     if (sgname.eq.'R32     ')  sgnum = 155 
     if (sgname.eq.'P3m1    ')  sgnum = 156 
     if (sgname.eq.'P31m    ')  sgnum = 157 
     if (sgname.eq.'P3c1    ')  sgnum = 158 
     if (sgname.eq.'P31c    ')  sgnum = 159 
     if (sgname.eq.'R3m     ')  sgnum = 160 
     if (sgname.eq.'R3c     ')  sgnum = 161 
     if (sgname.eq.'P-31m   ')  sgnum = 162 
     if (sgname.eq.'P-31c   ')  sgnum = 163 
     if (sgname.eq.'P-3m1   ')  sgnum = 164 
     if (sgname.eq.'P-3c1   ')  sgnum = 165 
     if (sgname.eq.'R-3m    ')  sgnum = 166 
     if (sgname.eq.'R-3c    ')  sgnum = 167 
     if (sgname.eq.'P6      ')  sgnum = 168 
     if (sgname.eq.'P61     ')  sgnum = 169 
     if (sgname.eq.'P65     ')  sgnum = 170 
     if (sgname.eq.'P62     ')  sgnum = 171 
     if (sgname.eq.'P64     ')  sgnum = 172 
     if (sgname.eq.'P63     ')  sgnum = 173 
     if (sgname.eq.'P-6     ')  sgnum = 174 
     if (sgname.eq.'P6/m    ')  sgnum = 175 
     if (sgname.eq.'P63/m   ')  sgnum = 176 
     if (sgname.eq.'P622    ')  sgnum = 177 
     if (sgname.eq.'P6122   ')  sgnum = 178 
     if (sgname.eq.'P6522   ')  sgnum = 179 
     if (sgname.eq.'P6222   ')  sgnum = 180 
     if (sgname.eq.'P6422   ')  sgnum = 181 
     if (sgname.eq.'P6322   ')  sgnum = 182 
     if (sgname.eq.'P6mm    ')  sgnum = 183 
     if (sgname.eq.'P6cc    ')  sgnum = 184 
     if (sgname.eq.'P63cm   ')  sgnum = 185 
     if (sgname.eq.'P63mc   ')  sgnum = 186 
     if (sgname.eq.'P-6m2   ')  sgnum = 187 
     if (sgname.eq.'P-6c2   ')  sgnum = 188 
     if (sgname.eq.'P-62m   ')  sgnum = 189 
     if (sgname.eq.'P-62c   ')  sgnum = 190 
     if (sgname.eq.'P6/mmm  ')  sgnum = 191 
     if (sgname.eq.'P6/mcc  ')  sgnum = 192 
     if (sgname.eq.'P63/mcm ')  sgnum = 193 
     if (sgname.eq.'P63/mmc ')  sgnum = 194 
     if (sgname.eq.'P23     ')  sgnum = 195 
     if (sgname.eq.'F23     ')  sgnum = 196 
     if (sgname.eq.'I23     ')  sgnum = 197 
     if (sgname.eq.'P213    ')  sgnum = 198 
     if (sgname.eq.'I213    ')  sgnum = 199 
     if (sgname.eq.'Pm-3    ')  sgnum = 200 
     if (sgname.eq.'Pn-3    ')  sgnum = 201 
     if (sgname.eq.'Fm-3    ')  sgnum = 202 
     if (sgname.eq.'Fd-3    ')  sgnum = 203 
     if (sgname.eq.'Im-3    ')  sgnum = 204 
     if (sgname.eq.'Pa-3    ')  sgnum = 205 
     if (sgname.eq.'Ia-3    ')  sgnum = 206 
     if (sgname.eq.'P432    ')  sgnum = 207 
     if (sgname.eq.'P4232   ')  sgnum = 208 
     if (sgname.eq.'F432    ')  sgnum = 209 
     if (sgname.eq.'F4132   ')  sgnum = 210 
     if (sgname.eq.'I432    ')  sgnum = 211 
     if (sgname.eq.'P4332   ')  sgnum = 212 
     if (sgname.eq.'P4132   ')  sgnum = 213 
     if (sgname.eq.'I4132   ')  sgnum = 214 
     if (sgname.eq.'P-43m   ')  sgnum = 215 
     if (sgname.eq.'F-43m   ')  sgnum = 216 
     if (sgname.eq.'I-43m   ')  sgnum = 217 
     if (sgname.eq.'P-43n   ')  sgnum = 218 
     if (sgname.eq.'F-43c   ')  sgnum = 219 
     if (sgname.eq.'I-43d   ')  sgnum = 220 
     if (sgname.eq.'Pm-3m   ')  sgnum = 221 
     if (sgname.eq.'Pn-3n   ')  sgnum = 222 
     if (sgname.eq.'Pm-3n   ')  sgnum = 223 
     if (sgname.eq.'Pn-3m   ')  sgnum = 224 
     if (sgname.eq.'Fm-3m   ')  sgnum = 225 
     if (sgname.eq.'Fm-3c   ')  sgnum = 226 
     if (sgname.eq.'Fd-3m   ')  sgnum = 227 
     if (sgname.eq.'Fd-3c   ')  sgnum = 228 
     if (sgname.eq.'Im-3m   ')  sgnum = 229 
     if (sgname.eq.'Ia-3d   ')  sgnum = 230 

   end subroutine getsgnum
