xsphSRC = \
./XSPH/xsph.f90 ./PAR/parallel.f90 ./COMMON/wlog.f90 \
./COMMON/str.f90 ./COMMON/chopen.f90 ./XSPH/rexsph.f90 \
./KSPACE/crystalstructure.f90 ./KSPACE/readcrystaldata.f90 ./KSPACE/pointgroup.f90 \
./MATH/invertmatrix.f90 ./KSPACE/spacegroup.f90 ./KSPACE/subtract_a.f90 \
./KSPACE/change_car.f90 ./KSPACE/symmetrycheck.f90 ./KSPACE/kmesh.f90 \
./COMMON/rdhead.f90 ./COMMON/setkap.f90 ./XSPH/xsphsub.f90 \
./COMMON/rdpot.f90 ./COMMON/padlib.f90 ./TDLDA/correorb.f90 \
./TDLDA/cdos.f90 ./MATH/besjn.f90 ./MATH/bjnser.f90 \
./FOVRG/inmuac.f90 ./COMMON/getorb.f90 ./RDINP/rdline.f90 \
./COMMON/nxtunt.f90 ./FOVRG/diff.f90 ./MATH/besjh.f90 \
./FOVRG/wfirdc.f90 ./FOVRG/nucdec.f90 ./ATOM/nucmass.f90 \
./FOVRG/potdvp.f90 ./FOVRG/aprdep.f90 ./FOVRG/solout.f90 \
./FOVRG/intout.f90 ./FOVRG/dfovrg.f90 ./FOVRG/muatcc.f90 \
./MATH/cwig3j.f90 ./FOVRG/potex.f90 ./FOVRG/aprdec.f90 \
./FOVRG/yzkrdc.f90 ./FOVRG/yzktec.f90 ./FOVRG/solin.f90 \
./MATH/phamp.f90 ./MATH/csomm2.f90 ./XSPH/getedg.f90 \
./COMMON/head.f90 ./XSPH/phmesh2T.f90 ./COMMON/getxk.f90 \
./XSPH/rdgrid.f90 ./COMMON/rdcmt.f90 ./XSPH/phmesh2.f90 \
./MATH/qsortd.f90 ./XSPH/phmeshjas.f90 ./XSPH/phmesh.f90 \
./TDLDA/meshlda.f90 ./FMS/kprep.f90 ./KSPACE/calccgc.f90 \
./KSPACE/strinit.f90 ./KSPACE/strvecgen.f90 ./KSPACE/strgaunt.f90 \
./KSPACE/cgcrac.f90 ./KSPACE/strsmat.f90 ./MATH/lu.f90 \
./MATH/seigen.f90 ./KSPACE/straa.f90 ./KSPACE/strconfra.f90 \
./KSPACE/strfunqjl.f90 ./KSPACE/strharpol.f90 ./KSPACE/strcc.f90 \
./KSPACE/change_eta.f90 ./KSPACE/bastrans.f90 ./KSPACE/makerotations.f90 \
./POT/istprm.f90 ./MATH/dist.f90 ./POT/sidx.f90 \
./COMMON/xx.f90 ./EXCH/vbh.f90 ./EXCH/edp.f90 \
./POT/movrlp.f90 ./POT/ovp2mt.f90 ./MATH/terp.f90 \
./MATH/polint.f90 ./MATH/somm2.f90 ./POT/fermi.f90 \
./XSPH/szlz.f90 ./XSPH/acoef.f90 ./POT/grids.f90 \
./COMMON/fixvar.f90 ./COMMON/fixdsx.f90 ./XSPH/rholat.f90 \
./MATH/exjlnl.f90 ./XSPH/rholsz.f90 ./XSPH/fmssz.f90 \
./FMS/yprep.f90 ./FMS/xstaff.f90 ./FMS/fmskspace.f90 \
./FMS/kkrintegral.f90 ./COMMON/writematrix.f90 ./KSPACE/structurefactor.f90 \
./KSPACE/strset.f90 ./KSPACE/strbbdd2.f90 ./KSPACE/strbbdd.f90 \
./FMS/fmspack.f90 ./MATH/ylm.f90 ./FMS/gglu.f90 \
./FMS/ggbi.f90 ./FMS/ggrm.f90 ./FMS/gggm.f90 \
./FMS/ggtf.f90 ./COMMON/fixdsp.f90 ./XSPH/getholeorb0.f90 \
./XSPH/xsectjas.f90 ./MATH/somm.f90 ./XSPH/bcoefjas.f90 \
./XSPH/mincalc.f90 ./XSPH/qbesselget.f90 ./XSPH/besjnjas.f90 \
./XSPH/radjas.f90 ./XSPH/xmultjas.f90 ./XSPH/csommjas.f90 \
./EXCH/xcpot.f90 ./SELF/csigz.f90 ./SELF/csigma.f90 \
./SELF/bpr1_2.f90 ./SELF/logi.f90 ./SELF/omegaq.f90 \
./MATH/xlogx.f90 ./SELF/bpr2_2.f90 ./SELF/bpr3_2.f90 \
./SELF/fndsng.f90 ./MATH/czeros.f90 ./EXCH/qsorti.f90 \
./DEBYE/sigms.f90 ./COMMON/pertab.f90 ./SELF/quinn.f90 \
./EXCH/rhl.f90 ./EXCH/imhl.f90 ./EXCH/ffq.f90 \
./EXCH/cubic.f90 ./EXCH/rhlbp.f90 ./XSPH/ljneeded0.f90 \
./XSPH/specupd.f90 ./MATH/cpl0.f90 ./XSPH/specupdatom.f90 \
./XSPH/specupdlg.f90 ./XSPH/xsect.f90 ./MATH/bcoef.f90 \
./MATH/rotwig.f90 ./TDLDA/phiscf.f90 ./TDLDA/lipman.f90 \
./TDLDA/chiklu.f90 ./XSPH/radint.f90 ./XSPH/xmult.f90 \
./MATH/csomm.f90 ./FF2X/xscorr.f90 ./MATH/terpc.f90 \
./TDLDA/rdpotp.f90 ./TDLDA/xsectd.f90 ./COMMON/setgam.f90 \
./TDLDA/getmat.f90 ./TDLDA/ridxmu.f90 ./TDLDA/getwf.f90 \
./TDLDA/ellfun.f90 ./TDLDA/yzktd.f90 ./MATH/conv.f90 \
./TDLDA/kkchi.f90 ./TDLDA/dmscf.f90 ./XSPH/phase.f90 \
./XSPH/phase_h.f90 ./XSPH/wphase.f90 ./XSPH/wrxsph.f90 \
./XSPH/axafs.f90 ./MATH/determ.f90 \
./EXCH/rgw.f90 ./EXCH/imgw.f90

xsph_MODULESRC = \
./PAR/m_par.f90 ./COMMON/m_constants.f90 ./KSPACE/m_struct.f90 \
./KSPACE/m_controls.f90 ./KSPACE/m_boundaries.f90 ./KSPACE/m_kklist.f90 \
./KSPACE/m_kgenwork.f90 ./KSPACE/m_tetrahedra.f90 ./KSPACE/m_controlkgen.f90 \
./COMMON/m_dimsmod.f90 ./KSPACE/m_strfacs.f90 ./COMMON/m_inpmodules.f90 \
./COMMON/m_nrixs.f90 ./COMMON/m_config.f90 ./ERRORMODS/m_errorfile.f90 \
./XSPH/m_elam.f90 ./KSPACE/m_energygrid.f90 ./KSPACE/m_workstrfacs2.f90 \
./KSPACE/m_workstrfacs.f90 ./KSPACE/m_workstrfacssimple.f90 ./KSPACE/m_trafo.f90 \
./KSPACE/m_wigner3j.f90 ./EXCH/m_pdw.f90 ./EXCH/m_pz.f90 \
./EXCH/m_ksdt.f90 ./COMMON/m_ifuns.f90 ./COMMON/m_rotx.f90 \
./COMMON/m_lnlm.f90 ./COMMON/m_xstruc.f90 ./COMMON/m_afctr.f90 \
./ERRORMODS/m_errormod.f90 ./IOMODS/m_iofiles.f90 ./IOMODS/m_padio.f90 \
./IOMODS/m_iomod.f90 ./COMMON/m_stkets.f90 ./FMS/m_fms.f90 \
./COMMON/m_t3j.f90 ./SELF/m_SelfEnergy.f90 ./TDLDA/getchi0.f90 \
./COMMON/m_kinds.f90 ./POT/m_mtdp.f90 ./POT/m_atomicpotio.f90 

