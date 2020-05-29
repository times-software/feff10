bandSRC = \
./BAND/band.f90 ./PAR/parallel.f90 ./COMMON/wlog.f90 \
./COMMON/str.f90 ./COMMON/chopen.f90 ./FMS/reafms.f90 \
./KSPACE/crystalstructure.f90 ./KSPACE/readcrystaldata.f90 ./KSPACE/pointgroup.f90 \
./MATH/invertmatrix.f90 ./KSPACE/spacegroup.f90 ./KSPACE/subtract_a.f90 \
./KSPACE/change_car.f90 ./KSPACE/symmetrycheck.f90 ./KSPACE/kmesh.f90 \
./COMMON/rdhead.f90 ./COMMON/setkap.f90 ./BAND/bandtot.f90 \
./BAND/ibravais.f90 ./COMMON/rdxsph.f90 ./COMMON/padlib.f90 \
./FMS/kprep.f90 ./MATH/cwig3j.f90 ./KSPACE/calccgc.f90 \
./KSPACE/strinit.f90 ./KSPACE/strvecgen.f90 ./KSPACE/strgaunt.f90 \
./KSPACE/cgcrac.f90 ./KSPACE/strsmat.f90 ./MATH/seigen.f90 \
./KSPACE/straa.f90 ./KSPACE/strconfra.f90 ./KSPACE/strfunqjl.f90 \
./KSPACE/strharpol.f90 ./KSPACE/strcc.f90 ./KSPACE/change_eta.f90 \
./KSPACE/bastrans.f90 ./KSPACE/makerotations.f90 ./BAND/kpath.f90 \
./MATH/terp.f90 ./MATH/polint.f90 ./BAND/fmsband.f90 \
./BAND/kkrband.f90 ./KSPACE/structurefactor.f90 ./KSPACE/strset.f90 \
./KSPACE/strbbdd2.f90 ./KSPACE/strbbdd.f90 ./FMS/kkrintegral.f90 \
./COMMON/writematrix.f90 ./BAND/bandblaslapack.f90 ./RDINP/rdline.f90 \
./COMMON/nxtunt.f90 
band_MODULESRC = \
./PAR/m_par.f90 ./COMMON/m_constants.f90 ./KSPACE/m_struct.f90 \
./KSPACE/m_controls.f90 ./KSPACE/m_boundaries.f90 ./KSPACE/m_kklist.f90 \
./KSPACE/m_kgenwork.f90 ./KSPACE/m_tetrahedra.f90 ./KSPACE/m_controlkgen.f90 \
./COMMON/m_dimsmod.f90 ./KSPACE/m_strfacs.f90 ./COMMON/m_inpmodules.f90 \
./COMMON/m_nrixs.f90 ./KSPACE/m_energygrid.f90 ./KSPACE/m_workstrfacs2.f90 \
./KSPACE/m_workstrfacs.f90 ./KSPACE/m_workstrfacssimple.f90 ./KSPACE/m_trafo.f90 \
./KSPACE/m_wigner3j.f90 ./ERRORMODS/m_errormod.f90 ./IOMODS/m_iofiles.f90 \
./IOMODS/m_padio.f90 ./IOMODS/m_iomod.f90 ./BAND/m_fitting.f90 \
./ERRORMODS/m_errorfile.f90 
