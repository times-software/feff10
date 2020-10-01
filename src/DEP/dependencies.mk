./ATOM/akeato.o: Compiler.mk
./ATOM/apot.o: ./POT/m_atomicpotio.o ./ERRORMODS/m_errormod.o ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./PAR/m_par.o ./COMMON/m_config.o ./COMMON/m_inpmodules.o Compiler.mk
./ATOM/aprdev.o: Compiler.mk
./ATOM/atomic.o: ./COMMON/m_dimsmod.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./ERRORMODS/m_errorfile.o Compiler.mk
./ATOM/bkmrdf.o: Compiler.mk
./ATOM/cofcon.o: Compiler.mk
./ATOM/dentfa.o: Compiler.mk
./ATOM/dsordf.o: ./ERRORMODS/m_errormod.o Compiler.mk
./ATOM/etotal.o: Compiler.mk
./ATOM/fdmocc.o: Compiler.mk
./ATOM/fdrirk.o: Compiler.mk
./ATOM/fpf0.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./PAR/m_par.o Compiler.mk
./ATOM/inmuat.o: Compiler.mk
./ATOM/intdir.o: Compiler.mk
./ATOM/lagdat.o: Compiler.mk
./ATOM/messer.o: Compiler.mk
./ATOM/muatco.o: Compiler.mk
./ATOM/nucdev.o: Compiler.mk
./ATOM/nucmass.o: Compiler.mk
./ATOM/ortdat.o: Compiler.mk
./ATOM/potrdf.o: Compiler.mk
./ATOM/potslw.o: Compiler.mk
./ATOM/s02at.o: Compiler.mk
./ATOM/scfdat.o: ./COMMON/m_dimsmod.o ./PAR/m_par.o Compiler.mk
./ATOM/soldir.o: Compiler.mk
./ATOM/tabrat.o: Compiler.mk
./ATOM/vlda.o: ./COMMON/m_constants.o Compiler.mk
./ATOM/wfirdf.o: Compiler.mk
./ATOM/yzkrdf.o: Compiler.mk
./ATOM/yzkteg.o: Compiler.mk
./BAND/band.o: ./COMMON/m_dimsmod.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./ERRORMODS/m_errorfile.o Compiler.mk
./BAND/bandblaslapack.o: Compiler.mk
./BAND/bandtot.o: ./COMMON/m_dimsmod.o ./KSPACE/m_boundaries.o ./COMMON/m_constants.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./BAND/m_fitting.o ./KSPACE/m_struct.o Compiler.mk
./BAND/brent.o: ./BAND/m_fitting.o Compiler.mk
./BAND/dbrent.o: Compiler.mk
./BAND/fmsband.o: ./COMMON/m_dimsmod.o ./KSPACE/m_struct.o ./KSPACE/m_strfacs.o ./KSPACE/m_boundaries.o ./KSPACE/m_wigner3j.o ./KSPACE/m_trafo.o ./KSPACE/m_controls.o ./IOMODS/m_iomod.o ./COMMON/m_constants.o Compiler.mk
./BAND/gauleg.o: Compiler.mk
./BAND/ibravais.o: Compiler.mk
./BAND/ikapmue.o: Compiler.mk
./BAND/invert.o: Compiler.mk
./BAND/kkrband.o: ./KSPACE/m_struct.o ./KSPACE/m_boundaries.o Compiler.mk
./BAND/kpath.o: Compiler.mk
./BAND/m_fitting.o: Compiler.mk
./BAND/mnbrak.o: ./BAND/m_fitting.o Compiler.mk
./BAND/polcoe.o: Compiler.mk
./COMMON/chopen.o: Compiler.mk
./COMMON/fixdsp.o: ./COMMON/m_dimsmod.o ./COMMON/m_ifuns.o ./COMMON/m_constants.o Compiler.mk
./COMMON/fixdsx.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_ifuns.o Compiler.mk
./COMMON/fixstr.o: Compiler.mk
./COMMON/fixvar.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_ifuns.o Compiler.mk
./COMMON/getcom.o: Compiler.mk
./COMMON/getorb.o: ./COMMON/m_config.o Compiler.mk
./COMMON/getxk.o: Compiler.mk
./COMMON/head.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./HEADERS/vers.h Compiler.mk
./COMMON/iniptz.o: Compiler.mk
./COMMON/isedge.o: Compiler.mk
./COMMON/itoken.o: Compiler.mk
./COMMON/m_afctr.o: Compiler.mk
./COMMON/m_config.o: ./COMMON/m_inpmodules.o Compiler.mk
./COMMON/m_constants.o: Compiler.mk
./COMMON/m_dimsmod.o: Compiler.mk
./COMMON/m_ifuns.o: Compiler.mk
./COMMON/m_inpmodules.o: ./COMMON/m_dimsmod.o ./KSPACE/m_controls.o ./KSPACE/m_struct.o ./KSPACE/m_kklist.o ./KSPACE/m_strfacs.o ./COMMON/m_constants.o Compiler.mk
./COMMON/m_kinds.o: Compiler.mk
./COMMON/m_lnlm.o: Compiler.mk
./COMMON/m_nrixs.o: ./COMMON/m_inpmodules.o ./COMMON/m_dimsmod.o Compiler.mk
./COMMON/m_readxmu.o: ./IOMODS/m_iomod.o Compiler.mk
./COMMON/m_rotx.o: Compiler.mk
./COMMON/m_stkets.o: ./COMMON/m_dimsmod.o Compiler.mk
./COMMON/m_t3j.o: Compiler.mk
./COMMON/m_xstruc.o: Compiler.mk
./COMMON/nxtunt.o: Compiler.mk
./COMMON/padlib.o: ./COMMON/padlib.h Compiler.mk
./COMMON/pertab.o: Compiler.mk
./COMMON/pijump.o: ./COMMON/m_constants.o Compiler.mk
./COMMON/rdcmt.o: Compiler.mk
./COMMON/rdhead.o: Compiler.mk
./COMMON/rdpot.o: ./COMMON/m_dimsmod.o Compiler.mk
./COMMON/rdxsph.o: ./COMMON/m_dimsmod.o Compiler.mk
./COMMON/rdxsph_h.o: ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o Compiler.mk
./COMMON/rdxsphjas.o: ./COMMON/m_dimsmod.o ./COMMON/m_nrixs.o ./COMMON/m_inpmodules.o Compiler.mk
./COMMON/setgam.o: Compiler.mk
./COMMON/setkap.o: Compiler.mk
./COMMON/stdnm.o: Compiler.mk
./COMMON/str.o: Compiler.mk
./COMMON/str2dp.o: Compiler.mk
./COMMON/wlog.o: ./PAR/m_par.o Compiler.mk
./COMMON/writematrix.o: Compiler.mk
./COMMON/xx.o: Compiler.mk
./COMPTON/compton.o: ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o ./COMPTON/m_compton.o ./ERRORMODS/m_errorfile.o ./PAR/m_par.o ./COMMON/m_constants.o ./RHORRP/m_rhorrp.o Compiler.mk
./COMPTON/m_compton.o: ./COMMON/m_constants.o ./COMPTON/m_rotation.o ./RHORRP/m_rhorrp.o ./COMMON/m_inpmodules.o ./PAR/m_par.o Compiler.mk
./COMPTON/m_rotation.o: Compiler.mk
./CRPA/chi_crpa.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./CRPA/crpa.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./PAR/m_par.o ./COMMON/m_inpmodules.o Compiler.mk
./DEBYE/sigcl.o: ./COMMON/m_constants.o Compiler.mk
./DEBYE/sigm3.o: Compiler.mk
./DEBYE/sigms.o: ./COMMON/m_constants.o Compiler.mk
./DEBYE/sigrem.o: ./PAR/m_par.o ./COMMON/m_dimsmod.o Compiler.mk
./DEBYE/sigte3.o: Compiler.mk
./DMDW/dmdw.o: ./COMMON/m_kinds.o ./DMDW/m_const_and_conv.o ./DMDW/m_dmdw.o ./COMMON/m_inpmodules.o ./ERRORMODS/m_errorfile.o ./PAR/m_par.o Compiler.mk
./DMDW/dym2feffinp.o: ./COMMON/m_kinds.o ./DMDW/m_const_and_conv.o ./DMDW/m_cmdline.o ./DMDW/m_dmdw.o Compiler.mk
./DMDW/m_cmdline.o: Compiler.mk
./DMDW/m_const_and_conv.o: ./COMMON/m_kinds.o Compiler.mk
./DMDW/m_dmdw.o: ./COMMON/m_kinds.o ./DMDW/m_const_and_conv.o ./INPGEN/m_ptable.o ./DMDW/m_math.o ./INPGEN/m_strings.o Compiler.mk
./DMDW/m_math.o: ./COMMON/m_kinds.o Compiler.mk
./EELS/angularmesh.o: ./COMMON/m_inpmodules.o ./EELS/m_qvectors.o ./EELS/m_work.o ./COMMON/m_constants.o Compiler.mk
./EELS/calculateweights.o: ./EELS/m_work.o ./EELS/m_qvectors.o ./COMMON/m_inpmodules.o ./COMMON/m_constants.o Compiler.mk
./EELS/concat.o: Compiler.mk
./EELS/eels.o: ./EELS/m_program_control.o ./EELS/m_qvectors.o ./COMMON/m_inpmodules.o ./EELS/m_work.o ./EELS/m_spectrum.o ./PAR/m_par.o ./COMMON/m_constants.o ./ERRORMODS/m_errorfile.o Compiler.mk
./EELS/euler.o: Compiler.mk
./EELS/m_program_control.o: Compiler.mk
./EELS/m_qvectors.o: Compiler.mk
./EELS/m_spectrum.o: Compiler.mk
./EELS/m_work.o: Compiler.mk
./EELS/productmatvect.o: Compiler.mk
./EELS/qmesh.o: ./EELS/m_program_control.o ./EELS/m_qvectors.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o Compiler.mk
./EELS/readsp.o: ./COMMON/m_inpmodules.o ./EELS/m_spectrum.o ./COMMON/m_constants.o Compiler.mk
./EELS/wavelength.o: ./COMMON/m_constants.o Compiler.mk
./EELS/writeangulardependence1.o: ./COMMON/m_inpmodules.o ./EELS/m_qvectors.o ./COMMON/m_constants.o ./KSPACE/m_energygrid.o ./EELS/m_program_control.o Compiler.mk
./EELS/writeangulardependence2.o: ./COMMON/m_inpmodules.o ./EELS/m_work.o ./EELS/m_qvectors.o ./EELS/m_program_control.o ./COMMON/m_constants.o ./EELS/m_spectrum.o Compiler.mk
./EELS/writeangulardependence3.o: ./COMMON/m_inpmodules.o ./EELS/m_work.o ./EELS/m_qvectors.o ./EELS/m_program_control.o ./COMMON/m_constants.o ./EELS/m_spectrum.o Compiler.mk
./EELSMDFF/mdff_angularmesh.o: ./COMMON/m_inpmodules.o ./EELSMDFF/mdff_m_qvectors.o ./EELSMDFF/mdff_m_work.o ./COMMON/m_constants.o Compiler.mk
./EELSMDFF/mdff_concat.o: Compiler.mk
./EELSMDFF/mdff_eels.o: ./EELSMDFF/mdff_m_program_control.o ./EELSMDFF/mdff_m_qvectors.o ./COMMON/m_inpmodules.o ./EELSMDFF/mdff_m_work.o ./EELSMDFF/mdff_m_spectrum.o ./PAR/m_par.o ./COMMON/m_constants.o Compiler.mk
./EELSMDFF/mdff_euler.o: Compiler.mk
./EELSMDFF/mdff_m_program_control.o: Compiler.mk
./EELSMDFF/mdff_m_qvectors.o: Compiler.mk
./EELSMDFF/mdff_m_spectrum.o: Compiler.mk
./EELSMDFF/mdff_m_work.o: Compiler.mk
./EELSMDFF/mdff_productmatvect.o: Compiler.mk
./EELSMDFF/mdff_qmesh.o: ./EELSMDFF/mdff_m_program_control.o ./EELSMDFF/mdff_m_qvectors.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o Compiler.mk
./EELSMDFF/mdff_readsp.o: ./COMMON/m_inpmodules.o ./EELSMDFF/mdff_m_spectrum.o ./COMMON/m_constants.o Compiler.mk
./EELSMDFF/mdff_wavelength.o: ./COMMON/m_constants.o Compiler.mk
./ERRORMODS/m_errorfile.o: Compiler.mk
./ERRORMODS/m_errormod.o: Compiler.mk
./ERRORMODS/m_errorstackmod.o: ./ERRORMODS/m_errormod.o ./IOMODS/m_iomod.o Compiler.mk
./EXCH/cubic.o: Compiler.mk
./EXCH/edp.o: ./PAR/m_par.o ./COMMON/m_constants.o Compiler.mk
./EXCH/ffq.o: Compiler.mk
./EXCH/imhl.o: ./COMMON/m_constants.o Compiler.mk
./EXCH/m_ksdt.o: Compiler.mk
./EXCH/m_pdw.o: ./COMMON/m_constants.o Compiler.mk
./EXCH/m_pz.o: Compiler.mk
./EXCH/qsorti.o: Compiler.mk
./EXCH/rhl.o: ./PAR/m_par.o ./COMMON/m_constants.o Compiler.mk
./EXCH/rhlbp.o: ./COMMON/m_constants.o Compiler.mk
./EXCH/vbh.o: Compiler.mk
./EXCH/xcpot.o: ./IOMODS/m_iomod.o ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./FF2X/dwaddl.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./FF2X/exconv.o: Compiler.mk
./FF2X/feffdt.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./HEADERS/vers.h Compiler.mk
./FF2X/ff2afs.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./FF2X/ff2afsjas.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./FF2X/ff2chi.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./FF2X/ff2chijas.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./FF2X/ff2gen.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./DMDW/m_dmdw.o ./HEADERS/vers.h Compiler.mk
./FF2X/ff2x.o: ./PAR/m_par.o ./COMMON/m_inpmodules.o ./COMMON/m_nrixs.o ./ERRORMODS/m_errorfile.o ./COMMON/m_dimsmod.o Compiler.mk
./FF2X/ff2xmu.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./FF2X/m_thermal_xscorr.o Compiler.mk
./FF2X/ff2xmujas.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./FF2X/fprime.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./PAR/m_par.o Compiler.mk
./FF2X/m_thermal_xscorr.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./FF2X/rdfbin.o: ./COMMON/m_dimsmod.o Compiler.mk
./FF2X/rdfbinl.o: ./COMMON/m_dimsmod.o Compiler.mk
./FF2X/reff2x.o: ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./COMMON/m_nrixs.o Compiler.mk
./FF2X/xscorr.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./FF2X/xscorratan.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./FF2X/xscorrjas.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./FMS/fms.o: ./COMMON/m_dimsmod.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./COMMON/m_stkets.o ./COMMON/m_rotx.o ./COMMON/m_lnlm.o ./COMMON/m_xstruc.o ./ERRORMODS/m_errorfile.o ./COMMON/m_t3j.o Compiler.mk
./FMS/fmskspace.o: ./COMMON/m_dimsmod.o ./KSPACE/m_struct.o ./KSPACE/m_kklist.o ./KSPACE/m_strfacs.o ./KSPACE/m_boundaries.o ./KSPACE/m_wigner3j.o ./KSPACE/m_trafo.o ./KSPACE/m_kgenwork.o ./KSPACE/m_controls.o ./IOMODS/m_iomod.o ./COMMON/m_constants.o Compiler.mk
./FMS/fmspack.o: ./COMMON/m_dimsmod.o ./KSPACE/m_controls.o ./COMMON/m_stkets.o ./COMMON/m_rotx.o ./COMMON/m_lnlm.o ./COMMON/m_xstruc.o ./COMMON/m_t3j.o ./PAR/m_par.o ./FMS/m_fms.o ./COMMON/m_constants.o Compiler.mk
./FMS/fmspack_h.o: ./COMMON/m_dimsmod.o ./COMMON/m_stkets.o ./COMMON/m_rotx.o ./COMMON/m_lnlm.o ./COMMON/m_xstruc.o ./COMMON/m_t3j.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./COMMON/m_constants.o Compiler.mk
./FMS/fmstot.o: ./COMMON/m_dimsmod.o ./KSPACE/m_controls.o ./KSPACE/m_kklist.o ./KSPACE/m_kgenwork.o ./KSPACE/m_boundaries.o ./IOMODS/m_iomod.o ./ERRORMODS/m_errormod.o ./COMMON/m_constants.o ./PAR/m_par.o ./COMMON/m_stkets.o ./COMMON/m_rotx.o ./COMMON/m_inpmodules.o ./COMMON/m_nrixs.o ./ERRORMODS/m_errorfile.o ./FMS/m_fms.o Compiler.mk
./FMS/ggbi.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_stkets.o Compiler.mk
./FMS/gggm.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_stkets.o Compiler.mk
./FMS/gglu.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./IOMODS/m_iomod.o ./COMMON/m_stkets.o ./FMS/m_fms.o Compiler.mk
./FMS/gglu_h.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./IOMODS/m_iomod.o ./COMMON/m_stkets.o ./KSPACE/m_controls.o Compiler.mk
./FMS/gglufullpot.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_stkets.o Compiler.mk
./FMS/ggrm.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_stkets.o Compiler.mk
./FMS/ggtf.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_stkets.o Compiler.mk
./FMS/kkrintegral.o: ./KSPACE/m_struct.o ./KSPACE/m_kklist.o ./KSPACE/m_boundaries.o ./KSPACE/m_strfacs.o ./KSPACE/m_kgenwork.o ./KSPACE/m_trafo.o ./KSPACE/m_controls.o Compiler.mk
./FMS/kprep.o: ./KSPACE/m_struct.o ./KSPACE/m_kklist.o ./KSPACE/m_strfacs.o ./KSPACE/m_workstrfacs.o ./KSPACE/m_workstrfacs2.o ./KSPACE/m_workstrfacssimple.o ./KSPACE/m_boundaries.o ./KSPACE/m_wigner3j.o ./KSPACE/m_trafo.o ./KSPACE/m_controls.o ./KSPACE/m_kgenwork.o ./KSPACE/m_energygrid.o Compiler.mk
./FMS/m_fms.o: Compiler.mk
./FMS/reafms.o: ./COMMON/m_dimsmod.o ./KSPACE/m_controls.o ./KSPACE/m_struct.o ./KSPACE/m_kklist.o ./KSPACE/m_strfacs.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./COMMON/m_nrixs.o Compiler.mk
./FMS/xprep.o: ./IOMODS/m_iomod.o ./COMMON/m_dimsmod.o ./PAR/m_par.o ./COMMON/m_constants.o ./COMMON/m_rotx.o ./COMMON/m_lnlm.o ./COMMON/m_xstruc.o ./COMMON/m_t3j.o ./DMDW/m_dmdw.o Compiler.mk
./FMS/xstaff.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./COMMON/m_rotx.o ./COMMON/m_lnlm.o ./COMMON/m_xstruc.o ./COMMON/m_afctr.o Compiler.mk
./FMS/yprep.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_rotx.o ./COMMON/m_lnlm.o ./COMMON/m_xstruc.o Compiler.mk
./FOVRG/aprdec.o: Compiler.mk
./FOVRG/aprdep.o: Compiler.mk
./FOVRG/dfovrg.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./FOVRG/diff.o: ./COMMON/m_dimsmod.o Compiler.mk
./FOVRG/dsordc.o: ./COMMON/m_dimsmod.o Compiler.mk
./FOVRG/inmuac.o: ./COMMON/m_dimsmod.o Compiler.mk
./FOVRG/intout.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./FOVRG/muatcc.o: ./COMMON/m_dimsmod.o Compiler.mk
./FOVRG/nucdec.o: ./COMMON/m_dimsmod.o Compiler.mk
./FOVRG/ortdac.o: ./COMMON/m_dimsmod.o Compiler.mk
./FOVRG/potdvp.o: ./COMMON/m_dimsmod.o Compiler.mk
./FOVRG/potex.o: ./COMMON/m_dimsmod.o Compiler.mk
./FOVRG/solin.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./FOVRG/solout.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./FOVRG/wfirdc.o: ./COMMON/m_dimsmod.o Compiler.mk
./FOVRG/yzkrdc.o: ./COMMON/m_dimsmod.o Compiler.mk
./FOVRG/yzktec.o: ./COMMON/m_dimsmod.o Compiler.mk
./FULLSPECTRUM/addedg.o: ./COMMON/m_constants.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/drdtrm.o: ./COMMON/m_constants.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/egrid.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/egrid_lin.o: ./COMMON/m_constants.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/fullspectrum.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o ./PAR/m_par.o ./ERRORMODS/m_errorfile.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/gtedgs.o: ./COMMON/m_constants.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/hamaker.o: ./COMMON/m_constants.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/kk.o: ./COMMON/m_constants.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/opcons.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/qsum.o: ./COMMON/m_constants.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/rdbkg.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/rddens.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/rdldos.o: ./COMMON/m_constants.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/rdop.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/rdpotp_fs.o: ./COMMON/m_dimsmod.o Compiler.mk
./FULLSPECTRUM/rdst.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/rdval.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/rdxmu.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/rdxmunorm.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./FULLSPECTRUM/sumrules.o: ./COMMON/m_constants.o ./FULLSPECTRUM/HEADERS/params.h Compiler.mk
./GENFMT/fmtrxi.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./GENFMT/m_genfmt.o Compiler.mk
./GENFMT/genfmt.o: ./COMMON/m_dimsmod.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./COMMON/m_nrixs.o ./ERRORMODS/m_errorfile.o Compiler.mk
./GENFMT/genfmtjas.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_nrixs.o ./COMMON/m_inpmodules.o ./GENFMT/m_genfmt.o ./HEADERS/vers.h Compiler.mk
./GENFMT/genfmtsub.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./GENFMT/m_genfmt.o ./HEADERS/vers.h Compiler.mk
./GENFMT/m_genfmt.o: ./COMMON/m_dimsmod.o Compiler.mk
./GENFMT/mmtr.o: ./COMMON/m_dimsmod.o ./GENFMT/m_genfmt.o ./COMMON/m_constants.o Compiler.mk
./GENFMT/mmtrjas.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_nrixs.o ./COMMON/m_inpmodules.o ./GENFMT/m_genfmt.o Compiler.mk
./GENFMT/mmtrjas0.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_nrixs.o ./GENFMT/m_genfmt.o Compiler.mk
./GENFMT/mmtrxi.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./GENFMT/m_genfmt.o Compiler.mk
./GENFMT/mmtrxijas.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./GENFMT/m_genfmt.o Compiler.mk
./GENFMT/mmtrxijas0.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./GENFMT/m_genfmt.o Compiler.mk
./GENFMT/rdpath.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./GENFMT/m_genfmt.o Compiler.mk
./GENFMT/regenf.o: ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./COMMON/m_nrixs.o ./GENFMT/m_genfmt.o Compiler.mk
./GENFMT/rot3i.o: ./COMMON/m_dimsmod.o ./GENFMT/m_genfmt.o Compiler.mk
./GENFMT/sclmz.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./GENFMT/m_genfmt.o Compiler.mk
./GENFMT/setlam.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./GENFMT/m_genfmt.o Compiler.mk
./GENFMT/snlm.o: ./COMMON/m_dimsmod.o ./GENFMT/m_genfmt.o Compiler.mk
./GENFMT/xstar.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./GENFMT/m_genfmt.o Compiler.mk
./HEADERS/feff.o: Compiler.mk
./INPGEN/m_pot_generator.o: ./COMMON/m_kinds.o ./INPGEN/m_strings.o ./INPGEN/m_ptable.o Compiler.mk
./INPGEN/m_ptable.o: ./COMMON/m_kinds.o ./INPGEN/m_strings.o Compiler.mk
./INPGEN/m_strings.o: Compiler.mk
./INPGEN/pot_generator_test.o: ./COMMON/m_kinds.o ./INPGEN/m_pot_generator.o ./INPGEN/m_ptable.o Compiler.mk
./IOMODS/m_iofiles.o: ./ERRORMODS/m_errormod.o Compiler.mk
./IOMODS/m_iomod.o: ./ERRORMODS/m_errormod.o ./IOMODS/m_iofiles.o ./IOMODS/m_padio.o Compiler.mk
./IOMODS/m_padio.o: ./ERRORMODS/m_errormod.o Compiler.mk
./KSPACE/bastrans.o: Compiler.mk
./KSPACE/calccgc.o: Compiler.mk
./KSPACE/cgcrac.o: Compiler.mk
./KSPACE/change_car.o: Compiler.mk
./KSPACE/change_eta.o: ./KSPACE/m_struct.o ./KSPACE/m_workstrfacs2.o ./KSPACE/m_boundaries.o Compiler.mk
./KSPACE/crystalstructure.o: ./KSPACE/m_struct.o ./COMMON/m_constants.o ./KSPACE/m_controls.o Compiler.mk
./KSPACE/kmesh.o: ./KSPACE/m_struct.o ./KSPACE/m_kklist.o ./KSPACE/m_kgenwork.o ./KSPACE/m_tetrahedra.o ./KSPACE/m_controlkgen.o Compiler.mk
./KSPACE/m_boundaries.o: ./KSPACE/m_controls.o Compiler.mk
./KSPACE/m_controlkgen.o: Compiler.mk
./KSPACE/m_controls.o: Compiler.mk
./KSPACE/m_energygrid.o: Compiler.mk
./KSPACE/m_kgenwork.o: Compiler.mk
./KSPACE/m_kklist.o: ./KSPACE/m_boundaries.o Compiler.mk
./KSPACE/m_strfacs.o: Compiler.mk
./KSPACE/m_struct.o: Compiler.mk
./KSPACE/m_tetrahedra.o: Compiler.mk
./KSPACE/m_trafo.o: Compiler.mk
./KSPACE/m_wigner3j.o: Compiler.mk
./KSPACE/m_workstrfacs.o: ./KSPACE/m_boundaries.o Compiler.mk
./KSPACE/m_workstrfacs2.o: ./KSPACE/m_boundaries.o Compiler.mk
./KSPACE/m_workstrfacssimple.o: ./KSPACE/m_boundaries.o ./KSPACE/m_workstrfacs.o ./KSPACE/m_workstrfacs2.o Compiler.mk
./KSPACE/makerotations.o: ./KSPACE/m_boundaries.o ./KSPACE/m_trafo.o ./KSPACE/m_kklist.o Compiler.mk
./KSPACE/pointgroup.o: Compiler.mk
./KSPACE/readcrystaldata.o: Compiler.mk
./KSPACE/spacegroup.o: ./COMMON/m_constants.o Compiler.mk
./KSPACE/straa.o: ./KSPACE/m_boundaries.o ./KSPACE/m_controls.o ./KSPACE/m_workstrfacs.o ./KSPACE/m_workstrfacs2.o Compiler.mk
./KSPACE/strbbdd.o: ./KSPACE/m_boundaries.o ./KSPACE/m_workstrfacs.o ./KSPACE/m_workstrfacs2.o Compiler.mk
./KSPACE/strbbdd2.o: ./KSPACE/m_boundaries.o ./KSPACE/m_workstrfacssimple.o ./KSPACE/m_workstrfacs2.o Compiler.mk
./KSPACE/strcc.o: ./KSPACE/m_boundaries.o ./KSPACE/m_workstrfacs.o ./KSPACE/m_workstrfacs2.o ./KSPACE/m_workstrfacssimple.o ./KSPACE/m_controls.o Compiler.mk
./KSPACE/strconfra.o: Compiler.mk
./KSPACE/strfunqjl.o: Compiler.mk
./KSPACE/strgaunt.o: ./KSPACE/m_controls.o Compiler.mk
./KSPACE/strharpol.o: ./KSPACE/m_boundaries.o ./KSPACE/m_workstrfacs.o ./KSPACE/m_workstrfacs2.o ./KSPACE/m_controls.o Compiler.mk
./KSPACE/strinit.o: ./KSPACE/m_workstrfacs.o ./KSPACE/m_workstrfacs2.o ./KSPACE/m_boundaries.o ./KSPACE/m_controls.o Compiler.mk
./KSPACE/strset.o: ./KSPACE/m_workstrfacs.o ./KSPACE/m_boundaries.o ./KSPACE/m_controls.o Compiler.mk
./KSPACE/strsmat.o: Compiler.mk
./KSPACE/structurefactor.o: ./KSPACE/m_struct.o ./KSPACE/m_trafo.o ./KSPACE/m_boundaries.o ./COMMON/m_constants.o Compiler.mk
./KSPACE/strvecgen.o: ./KSPACE/m_energygrid.o ./KSPACE/m_controls.o ./KSPACE/m_boundaries.o ./KSPACE/m_workstrfacs2.o ./KSPACE/m_workstrfacs.o Compiler.mk
./KSPACE/subtract_a.o: ./COMMON/m_constants.o Compiler.mk
./KSPACE/symmetrycheck.o: ./COMMON/m_constants.o Compiler.mk
./LDOS/ff2rho.o: ./COMMON/m_dimsmod.o ./KSPACE/m_controls.o ./COMMON/m_constants.o ./PAR/m_par.o Compiler.mk
./LDOS/ff2rho_h.o: ./PAR/m_par.o ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o Compiler.mk
./LDOS/ff2rho_h_step1.o: ./PAR/m_par.o ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o Compiler.mk
./LDOS/ff2rho_h_step2.o: ./PAR/m_par.o ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./LDOS/fmsdos.o: ./COMMON/m_dimsmod.o ./KSPACE/m_controls.o ./KSPACE/m_kklist.o ./COMMON/m_constants.o ./PAR/m_par.o ./ERRORMODS/m_errorfile.o Compiler.mk
./LDOS/fmsdos_h.o: ./COMMON/m_dimsmod.o ./KSPACE/m_controls.o ./COMMON/m_constants.o ./PAR/m_par.o ./ERRORMODS/m_errorfile.o ./COMMON/m_inpmodules.o Compiler.mk
./LDOS/fmsdos_h_step1.o: ./COMMON/m_dimsmod.o ./KSPACE/m_controls.o ./COMMON/m_constants.o ./PAR/m_par.o ./ERRORMODS/m_errorfile.o ./COMMON/m_inpmodules.o Compiler.mk
./LDOS/fmsdos_h_step2.o: ./COMMON/m_dimsmod.o ./KSPACE/m_controls.o ./COMMON/m_constants.o ./PAR/m_par.o ./ERRORMODS/m_errorfile.o ./COMMON/m_inpmodules.o Compiler.mk
./LDOS/ldos.o: ./COMMON/m_dimsmod.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./ERRORMODS/m_errorfile.o Compiler.mk
./LDOS/ldos_driver.o: ./COMMON/m_dimsmod.o ./COMMON/m_rotx.o ./COMMON/m_stkets.o ./COMMON/m_lnlm.o ./COMMON/m_xstruc.o ./COMMON/m_t3j.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./ERRORMODS/m_errorfile.o Compiler.mk
./LDOS/ldossub.o: ./KSPACE/m_controls.o ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./PAR/m_par.o ./COMMON/m_inpmodules.o Compiler.mk
./LDOS/ldossub_h.o: ./KSPACE/m_controls.o ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./PAR/m_par.o ./COMMON/m_inpmodules.o Compiler.mk
./LDOS/ldossub_h_unrolled.o: ./KSPACE/m_controls.o ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./PAR/m_par.o ./COMMON/m_inpmodules.o Compiler.mk
./LDOS/reldos.o: ./COMMON/m_dimsmod.o ./KSPACE/m_controls.o ./KSPACE/m_struct.o ./KSPACE/m_kklist.o ./KSPACE/m_strfacs.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o Compiler.mk
./LDOS/rhol.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./LDOS/rhol_h.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o Compiler.mk
./LDOS/rhol_h_step1.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./LDOS/rhol_h_step2.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o Compiler.mk
./MATH/bcoef.o: ./COMMON/m_dimsmod.o Compiler.mk
./MATH/besjh.o: ./COMMON/m_dimsmod.o Compiler.mk
./MATH/besjn.o: ./COMMON/m_dimsmod.o Compiler.mk
./MATH/bjnser.o: Compiler.mk
./MATH/conv.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./MATH/cpl0.o: Compiler.mk
./MATH/csomm.o: Compiler.mk
./MATH/csomm2.o: Compiler.mk
./MATH/cwig3j.o: Compiler.mk
./MATH/czeros.o: Compiler.mk
./MATH/determ.o: Compiler.mk
./MATH/dist.o: Compiler.mk
./MATH/exjlnl.o: Compiler.mk
./MATH/invertmatrix.o: Compiler.mk
./MATH/lint.o: Compiler.mk
./MATH/lu.o: Compiler.mk
./MATH/m_polyfit.o: Compiler.mk
./MATH/phamp.o: ./COMMON/m_constants.o Compiler.mk
./MATH/polint.o: Compiler.mk
./MATH/qsortd.o: Compiler.mk
./MATH/quartc.o: Compiler.mk
./MATH/rotwig.o: Compiler.mk
./MATH/sdist.o: Compiler.mk
./MATH/seigen.o: Compiler.mk
./MATH/somm.o: Compiler.mk
./MATH/somm2.o: Compiler.mk
./MATH/strap.o: Compiler.mk
./MATH/terp.o: Compiler.mk
./MATH/terpc.o: Compiler.mk
./MATH/trap.o: Compiler.mk
./MATH/xlogx.o: Compiler.mk
./MATH/ylm.o: Compiler.mk
./MKGTR/calclbcoef.o: Compiler.mk
./MKGTR/getgtr.o: ./COMMON/m_dimsmod.o ./IOMODS/m_iomod.o ./ERRORMODS/m_errormod.o ./COMMON/m_constants.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./COMMON/m_nrixs.o Compiler.mk
./MKGTR/getgtrjas.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./IOMODS/m_iomod.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./COMMON/m_nrixs.o Compiler.mk
./MKGTR/mkgtr.o: ./COMMON/m_dimsmod.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./ERRORMODS/m_errorfile.o Compiler.mk
./MKGTR/rotgmatrix.o: Compiler.mk
./OPCONSAT/addeps.o: ./IOMODS/m_iomod.o ./IOMODS/m_iofiles.o Compiler.mk
./OPCONSAT/epsdb.o: ./IOMODS/m_iomod.o Compiler.mk
./OPCONSAT/getelement.o: Compiler.mk
./OPCONSAT/opconsat.o: ./COMMON/m_inpmodules.o ./COMMON/m_constants.o ./POT/m_atomicpotio.o ./PAR/m_par.o ./ERRORMODS/m_errorfile.o ./COMMON/m_dimsmod.o Compiler.mk
./PAR/m_par.o: Compiler.mk
./PAR/parallel.o: ./PAR/m_par.o Compiler.mk
./PATH/ccrit.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./PATH/heap.o: Compiler.mk
./PATH/ipack.o: Compiler.mk
./PATH/mcrith.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./PATH/mcritk.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./PATH/mpprmd.o: ./COMMON/m_dimsmod.o Compiler.mk
./PATH/mpprmp.o: ./COMMON/m_dimsmod.o Compiler.mk
./PATH/mrb.o: ./COMMON/m_dimsmod.o Compiler.mk
./PATH/outcrt.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./PATH/path.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./ERRORMODS/m_errorfile.o Compiler.mk
./PATH/paths.o: ./COMMON/m_dimsmod.o Compiler.mk
./PATH/pathsd.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o Compiler.mk
./PATH/phash.o: ./COMMON/m_dimsmod.o Compiler.mk
./PATH/prcrit.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./COMMON/m_nrixs.o Compiler.mk
./PATH/repath.o: ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./COMMON/m_nrixs.o Compiler.mk
./PATH/sortix.o: Compiler.mk
./PATH/timrep.o: ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o Compiler.mk
./POT/afolp.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./POT/broydn.o: ./POT/m_broydn_workspace.o ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./POT/corval.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o Compiler.mk
./POT/coulom.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./POT/fermi.o: ./COMMON/m_constants.o Compiler.mk
./POT/ff2g.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./POT/fmsie.o: ./KSPACE/m_controls.o ./COMMON/m_dimsmod.o Compiler.mk
./POT/frnrm.o: ./COMMON/m_dimsmod.o Compiler.mk
./POT/grids.o: ./COMMON/m_constants.o Compiler.mk
./POT/importpot.o: ./COMMON/m_dimsmod.o Compiler.mk
./POT/inipot.o: ./COMMON/m_dimsmod.o Compiler.mk
./POT/istprm.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./EXCH/m_pdw.o ./COMMON/m_inpmodules.o ./EXCH/m_pz.o ./EXCH/m_ksdt.o Compiler.mk
./POT/istval.o: ./COMMON/m_dimsmod.o Compiler.mk
./POT/kpreppot.o: ./KSPACE/m_struct.o ./KSPACE/m_kklist.o ./KSPACE/m_strfacs.o ./KSPACE/m_workstrfacs.o ./KSPACE/m_workstrfacs2.o ./KSPACE/m_boundaries.o ./KSPACE/m_wigner3j.o ./KSPACE/m_trafo.o ./KSPACE/m_controls.o ./KSPACE/m_kgenwork.o ./KSPACE/m_energygrid.o Compiler.mk
./POT/m_atomicpotio.o: ./ERRORMODS/m_errormod.o ./IOMODS/m_iomod.o ./POT/m_mtdp.o ./IOMODS/m_iofiles.o Compiler.mk
./POT/m_broydn_workspace.o: ./COMMON/m_inpmodules.o ./COMMON/m_dimsmod.o Compiler.mk
./POT/m_mtdp.o: ./COMMON/m_kinds.o Compiler.mk
./POT/m_thermal_scf.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./PAR/m_par.o Compiler.mk
./POT/moveh.o: ./COMMON/m_dimsmod.o Compiler.mk
./POT/movrlp.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./POT/ovp2mt.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./POT/ovrlp.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./POT/pot.o: ./COMMON/m_dimsmod.o ./COMMON/m_stkets.o ./COMMON/m_rotx.o ./COMMON/m_lnlm.o ./COMMON/m_xstruc.o ./COMMON/m_t3j.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./ERRORMODS/m_errorfile.o ./POT/m_broydn_workspace.o Compiler.mk
./POT/potsub.o: ./POT/m_atomicpotio.o ./COMMON/m_dimsmod.o ./PAR/m_par.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./KSPACE/m_workstrfacs2.o ./KSPACE/m_controls.o ./POT/m_thermal_scf.o Compiler.mk
./POT/reapot.o: ./KSPACE/m_controls.o ./KSPACE/m_struct.o ./KSPACE/m_kklist.o ./KSPACE/m_strfacs.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o Compiler.mk
./POT/rhofmslie.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o Compiler.mk
./POT/rholie.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./POT/scmt.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./POT/scmtmp.o: ./COMMON/m_dimsmod.o ./PAR/m_par.o ./COMMON/m_constants.o Compiler.mk
./POT/sidx.o: Compiler.mk
./POT/sumax.o: Compiler.mk
./POT/wpot.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./POT/wrpot.o: ./COMMON/m_dimsmod.o Compiler.mk
./RDINP/cif2feffatoms.o: Compiler.mk
./RDINP/ciftbx.o: ./RDINP/ciftbx.sys Compiler.mk
./RDINP/clearfp.o: Compiler.mk
./RDINP/consistency.o: Compiler.mk
./RDINP/ffsort.o: ./PAR/m_par.o ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o Compiler.mk
./RDINP/findz.o: Compiler.mk
./RDINP/fixlinenow.o: ./COMMON/m_dimsmod.o Compiler.mk
./RDINP/hash_funcs.o: Compiler.mk
./RDINP/importcif.o: ./COMMON/m_dimsmod.o ./KSPACE/m_struct.o ./COMMON/m_inpmodules.o ./PAR/m_par.o ./COMMON/m_constants.o ./RDINP/ciftbx.cmn Compiler.mk
./RDINP/iniall.o: ./COMMON/m_inpmodules.o Compiler.mk
./RDINP/mkptz.o: ./COMMON/m_constants.o ./COMMON/m_inpmodules.o Compiler.mk
./RDINP/rdinp.o: ./PAR/m_par.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./COMMON/m_dimsmod.o ./ERRORMODS/m_errorfile.o ./HEADERS/vers.h Compiler.mk
./RDINP/rdline.o: Compiler.mk
./RDINP/setedg.o: Compiler.mk
./RDINP/spgroup1.o: Compiler.mk
./RDINP/spgroup2.o: Compiler.mk
./RDINP/spgroup3.o: Compiler.mk
./RDINP/wrtall.o: ./PAR/m_par.o ./COMMON/m_inpmodules.o Compiler.mk
./RHORRP/m_density_inp.o: ./COMMON/m_constants.o Compiler.mk
./RHORRP/m_rhorrp.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./PAR/m_par.o ./COMMON/m_ifuns.o ./MATH/m_polyfit.o Compiler.mk
./RHORRP/rhorrp.o: ./RHORRP/m_rhorrp.o ./RHORRP/m_density_inp.o ./COMMON/m_constants.o ./PAR/m_par.o ./ERRORMODS/m_errorfile.o Compiler.mk
./RIXS/blinterp2d.o: Compiler.mk
./RIXS/doublelorentz.o: Compiler.mk
./RIXS/kkint.o: Compiler.mk
./RIXS/rdxsphrxs.o: ./COMMON/m_dimsmod.o Compiler.mk
./RIXS/rixs.o: ./COMMON/m_dimsmod.o ./ERRORMODS/m_errorfile.o ./IOMODS/m_iomod.o ./PAR/m_par.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./COMMON/m_nrixs.o ./POT/m_atomicpotio.o ./COMMON/m_readxmu.o Compiler.mk
./RIXS/test.o: Compiler.mk
./SCREEN/fegrid.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./SCREEN/frgrid.o: ./COMMON/m_dimsmod.o Compiler.mk
./SCREEN/fxc.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./SCREEN/getph.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./SCREEN/prep.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o Compiler.mk
./SCREEN/rdgeom.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./COMMON/m_stkets.o ./COMMON/m_rotx.o ./COMMON/m_lnlm.o ./COMMON/m_xstruc.o ./COMMON/m_t3j.o Compiler.mk
./SCREEN/screen.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./ERRORMODS/m_errorfile.o Compiler.mk
./SCREEN/screensub.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./SELF/bpr1_2.o: Compiler.mk
./SELF/bpr2_2.o: Compiler.mk
./SELF/bpr3_2.o: Compiler.mk
./SELF/calcse.o: ./SELF/m_SelfEnergy.o Compiler.mk
./SELF/csigma.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./SELF/csigz.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./SELF/eps2exc.o: ./IOMODS/m_iomod.o ./COMMON/m_dimsmod.o ./SELF/m_SelfEnergy.o Compiler.mk
./SELF/fndsng.o: Compiler.mk
./SELF/logi.o: Compiler.mk
./SELF/m_SelfEnergy.o: ./IOMODS/m_iomod.o ./ERRORMODS/m_errormod.o ./COMMON/m_constants.o Compiler.mk
./SELF/omegaq.o: Compiler.mk
./SELF/quinn.o: ./COMMON/m_constants.o Compiler.mk
./SELF/testbp.o: Compiler.mk
./SELF/testrs.o: Compiler.mk
./SELF/testwq.o: Compiler.mk
./SFCONV/croots.o: Compiler.mk
./SFCONV/grater.o: Compiler.mk
./SFCONV/interpsf.o: Compiler.mk
./SFCONV/mkrmu.o: Compiler.mk
./SFCONV/mksat.o: Compiler.mk
./SFCONV/mkspectf.o: Compiler.mk
./SFCONV/plset.o: Compiler.mk
./SFCONV/ppole.o: Compiler.mk
./SFCONV/qlimits.o: Compiler.mk
./SFCONV/rdeps.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./SFCONV/senergies.o: Compiler.mk
./SFCONV/sfconv.o: ./PAR/m_par.o ./COMMON/m_inpmodules.o ./ERRORMODS/m_errorfile.o Compiler.mk
./SFCONV/sfconvsub.o: Compiler.mk
./SFCONV/so2conv.o: ./COMMON/m_inpmodules.o ./ERRORMODS/m_errorfile.o Compiler.mk
./TDLDA/cdos.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./TDLDA/chiklu.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./TDLDA/correorb.o: ./ERRORMODS/m_errorfile.o ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./TDLDA/dmscf.o: Compiler.mk
./TDLDA/ellfun.o: Compiler.mk
./TDLDA/getchi0.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./TDLDA/getmat.o: Compiler.mk
./TDLDA/getwf.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./TDLDA/kkchi.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./TDLDA/lipman.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./TDLDA/meshlda.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./TDLDA/phiscf.o: ./COMMON/m_dimsmod.o ./ERRORMODS/m_errorfile.o ./COMMON/m_constants.o Compiler.mk
./TDLDA/rdpotp.o: ./COMMON/m_dimsmod.o Compiler.mk
./TDLDA/ridxmu.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./TDLDA/xsectd.o: ./IOMODS/m_iomod.o ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./SELF/m_SelfEnergy.o Compiler.mk
./TDLDA/yzktd.o: ./COMMON/m_dimsmod.o Compiler.mk
./XSPH/acoef.o: ./COMMON/m_dimsmod.o Compiler.mk
./XSPH/axafs.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./XSPH/bcoefjas.o: ./COMMON/m_nrixs.o Compiler.mk
./XSPH/besjnjas.o: Compiler.mk
./XSPH/csommjas.o: Compiler.mk
./XSPH/fmssz.o: ./COMMON/m_dimsmod.o ./KSPACE/m_controls.o Compiler.mk
./XSPH/getedg.o: ./COMMON/m_constants.o ./XSPH/m_elam.o ./COMMON/m_dimsmod.o Compiler.mk
./XSPH/getholeorb0.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./XSPH/getoccnorm.o: Compiler.mk
./XSPH/ljneeded0.o: Compiler.mk
./XSPH/m_elam.o: Compiler.mk
./XSPH/mincalc.o: Compiler.mk
./XSPH/phase.o: ./IOMODS/m_iomod.o ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./SELF/m_SelfEnergy.o ./COMMON/m_inpmodules.o Compiler.mk
./XSPH/phase_h.o: ./IOMODS/m_iomod.o ./COMMON/m_dimsmod.o ./COMMON/m_constants.o ./SELF/m_SelfEnergy.o Compiler.mk
./XSPH/phmesh.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./XSPH/phmesh2.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o Compiler.mk
./XSPH/phmesh2T.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o Compiler.mk
./XSPH/phmeshjas.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./XSPH/qbesselget.o: Compiler.mk
./XSPH/radint.o: ./COMMON/m_dimsmod.o Compiler.mk
./XSPH/radjas.o: ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o Compiler.mk
./XSPH/rdgrid.o: ./COMMON/m_constants.o Compiler.mk
./XSPH/rexsph.o: ./KSPACE/m_controls.o ./KSPACE/m_struct.o ./KSPACE/m_kklist.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./COMMON/m_nrixs.o ./COMMON/m_dimsmod.o Compiler.mk
./XSPH/rholat.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./XSPH/rholsz.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./XSPH/specupd.o: ./COMMON/m_inpmodules.o Compiler.mk
./XSPH/specupdatom.o: ./COMMON/m_inpmodules.o Compiler.mk
./XSPH/specupdlg.o: ./COMMON/m_inpmodules.o Compiler.mk
./XSPH/szlz.o: ./COMMON/m_constants.o ./COMMON/m_dimsmod.o Compiler.mk
./XSPH/wphase.o: ./COMMON/m_dimsmod.o ./COMMON/m_constants.o Compiler.mk
./XSPH/wrxsph.o: ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o Compiler.mk
./XSPH/xmult.o: ./COMMON/m_constants.o Compiler.mk
./XSPH/xmultjas.o: Compiler.mk
./XSPH/xsect.o: ./IOMODS/m_iomod.o ./COMMON/m_constants.o ./COMMON/m_dimsmod.o ./SELF/m_SelfEnergy.o Compiler.mk
./XSPH/xsectjas.o: ./COMMON/m_dimsmod.o ./IOMODS/m_iomod.o ./SELF/m_SelfEnergy.o ./COMMON/m_constants.o ./COMMON/m_inpmodules.o ./COMMON/m_nrixs.o Compiler.mk
./XSPH/xsph.o: ./COMMON/m_dimsmod.o ./COMMON/m_stkets.o ./COMMON/m_rotx.o ./COMMON/m_lnlm.o ./COMMON/m_xstruc.o ./COMMON/m_t3j.o ./PAR/m_par.o ./COMMON/m_inpmodules.o ./ERRORMODS/m_errorfile.o Compiler.mk
./XSPH/xsphsub.o: ./KSPACE/m_controls.o ./IOMODS/m_iomod.o ./POT/m_atomicpotio.o ./COMMON/m_dimsmod.o ./COMMON/m_inpmodules.o ./COMMON/m_constants.o ./COMMON/m_nrixs.o ./ERRORMODS/m_errorfile.o Compiler.mk
