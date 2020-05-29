atomicSRC = \
./ATOM/atomic.f90 ./PAR/parallel.f90 ./COMMON/wlog.f90 \
./COMMON/str.f90 ./COMMON/chopen.f90 ./POT/reapot.f90 \
./KSPACE/crystalstructure.f90 ./KSPACE/readcrystaldata.f90 ./KSPACE/pointgroup.f90 \
./MATH/invertmatrix.f90 ./KSPACE/spacegroup.f90 ./KSPACE/subtract_a.f90 \
./KSPACE/change_car.f90 ./KSPACE/symmetrycheck.f90 ./KSPACE/kmesh.f90 \
./POT/kpreppot.f90 ./MATH/cwig3j.f90 ./KSPACE/calccgc.f90 \
./KSPACE/strinit.f90 ./KSPACE/strvecgen.f90 ./KSPACE/strgaunt.f90 \
./KSPACE/cgcrac.f90 ./KSPACE/strsmat.f90 ./MATH/seigen.f90 \
./KSPACE/straa.f90 ./KSPACE/strconfra.f90 ./KSPACE/strfunqjl.f90 \
./KSPACE/strharpol.f90 ./KSPACE/strcc.f90 ./KSPACE/change_eta.f90 \
./KSPACE/bastrans.f90 ./KSPACE/makerotations.f90 ./COMMON/rdhead.f90 \
./ATOM/apot.f90 ./RDINP/setedg.f90 ./POT/moveh.f90 \
./ATOM/scfdat.f90 ./ATOM/dsordf.f90 ./ATOM/aprdev.f90 \
./ATOM/inmuat.f90 ./COMMON/getorb.f90 ./RDINP/rdline.f90 \
./COMMON/nxtunt.f90 ./ATOM/wfirdf.f90 ./ATOM/dentfa.f90 \
./ATOM/nucdev.f90 ./ATOM/soldir.f90 ./ATOM/intdir.f90 \
./ATOM/messer.f90 ./ATOM/muatco.f90 ./ATOM/ortdat.f90 \
./ATOM/lagdat.f90 ./ATOM/akeato.f90 ./ATOM/fdrirk.f90 \
./ATOM/yzkrdf.f90 ./ATOM/yzkteg.f90 ./ATOM/potrdf.f90 \
./ATOM/vlda.f90 ./EXCH/vbh.f90 ./EXCH/edp.f90 \
./ATOM/cofcon.f90 ./ATOM/tabrat.f90 ./ATOM/etotal.f90 \
./ATOM/fdmocc.f90 ./ATOM/bkmrdf.f90 ./MATH/somm.f90 \
./ATOM/potslw.f90 ./ATOM/fpf0.f90 ./ATOM/s02at.f90 \
./MATH/determ.f90 ./POT/ovrlp.f90 ./MATH/dist.f90 \
./POT/sumax.f90 ./COMMON/xx.f90 ./POT/frnrm.f90 \
./MATH/somm2.f90 ./COMMON/setgam.f90 ./MATH/terp.f90 \
./MATH/polint.f90 ./COMMON/padlib.f90 
atomic_MODULESRC = \
./PAR/m_par.f90 ./COMMON/m_constants.f90 ./KSPACE/m_struct.f90 \
./KSPACE/m_controls.f90 ./KSPACE/m_boundaries.f90 ./KSPACE/m_kklist.f90 \
./KSPACE/m_kgenwork.f90 ./KSPACE/m_tetrahedra.f90 ./KSPACE/m_controlkgen.f90 \
./KSPACE/m_energygrid.f90 ./KSPACE/m_workstrfacs2.f90 ./KSPACE/m_workstrfacs.f90 \
./KSPACE/m_workstrfacssimple.f90 ./KSPACE/m_trafo.f90 ./KSPACE/m_strfacs.f90 \
./KSPACE/m_wigner3j.f90 ./COMMON/m_dimsmod.f90 ./COMMON/m_inpmodules.f90 \
./ERRORMODS/m_errormod.f90 ./COMMON/m_config.f90 ./IOMODS/m_iofiles.f90 \
./IOMODS/m_padio.f90 ./IOMODS/m_iomod.f90 ./COMMON/m_kinds.f90 \
./POT/m_mtdp.f90 ./POT/m_atomicpotio.f90 ./ERRORMODS/m_errorfile.f90 
