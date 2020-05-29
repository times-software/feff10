comptonSRC = \
./COMPTON/compton.f90 ./PAR/parallel.f90 ./COMMON/wlog.f90 \
./COMMON/str.f90 ./COMMON/chopen.f90 ./COMMON/rdhead.f90 \
./KSPACE/readcrystaldata.f90 ./RDINP/mkptz.f90 ./COMMON/rdpot.f90 \
./COMMON/padlib.f90 ./COMMON/rdxsph.f90 ./COMMON/fixvar.f90 \
./MATH/terp.f90 ./MATH/polint.f90 ./COMMON/fixdsx.f90 \
./FOVRG/dfovrg.f90 ./FOVRG/inmuac.f90 ./COMMON/getorb.f90 \
./RDINP/rdline.f90 ./COMMON/nxtunt.f90 ./FOVRG/diff.f90 \
./FOVRG/wfirdc.f90 ./FOVRG/nucdec.f90 ./FOVRG/potdvp.f90 \
./FOVRG/aprdep.f90 ./FOVRG/solout.f90 ./FOVRG/intout.f90 \
./FOVRG/solin.f90 ./MATH/besjh.f90 ./MATH/bjnser.f90 \
./FOVRG/muatcc.f90 ./MATH/cwig3j.f90 ./FOVRG/potex.f90 \
./FOVRG/aprdec.f90 ./FOVRG/yzkrdc.f90 ./FOVRG/yzktec.f90 \
./MATH/besjn.f90 ./MATH/phamp.f90 ./MATH/exjlnl.f90 \
./MATH/cpl0.f90 ./MATH/ylm.f90 ./MATH/terpc.f90 \
./MATH/lu.f90 ./MATH/seigen.f90 
compton_MODULESRC = \
./PAR/m_par.f90 ./COMMON/m_dimsmod.f90 ./KSPACE/m_controls.f90 \
./KSPACE/m_struct.f90 ./KSPACE/m_boundaries.f90 ./KSPACE/m_kklist.f90 \
./KSPACE/m_strfacs.f90 ./COMMON/m_constants.f90 ./COMMON/m_inpmodules.f90 \
./COMPTON/m_rotation.f90 ./COMMON/m_ifuns.f90 ./COMMON/m_config.f90 \
./MATH/m_polyfit.f90 ./RHORRP/m_rhorrp.f90 ./COMPTON/m_compton.f90 \
./ERRORMODS/m_errorfile.f90 
