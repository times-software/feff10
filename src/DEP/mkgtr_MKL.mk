mkgtrSRC = \
./MKGTR/mkgtr.f90 ./PAR/par.f90 ./COMMON/wlog.f90 \
./COMMON/str.f90 ./COMMON/chopen.f90 ./FMS/reafms.f90 \
./KSPACE/crystalstructure.f90 ./KSPACE/readcrystaldata.f90 ./KSPACE/pointgroup.f90 \
./MATH/invertmatrix.f90 ./KSPACE/spacegroup.f90 ./KSPACE/subtract_a.f90 \
./KSPACE/change_car.f90 ./KSPACE/symmetrycheck.f90 ./KSPACE/kmesh.f90 \
./COMMON/rdhead.f90 ./COMMON/setkap.f90 ./MKGTR/getgtrjas.f90 \
./COMMON/rdxsphjas.f90 ./COMMON/padlib.f90 ./MATH/cpl0.f90 \
./XSPH/bcoefjas.f90 ./MATH/cwig3j.f90 ./MKGTR/calclbcoef.f90 \
./MKGTR/rotgmatrix.f90 ./MATH/rotwig.f90 ./RDINP/rdline.f90 \
./COMMON/nxtunt.f90 ./MKGTR/getgtr.f90 ./COMMON/rdxsph.f90 \
./COMMON/iniptz.f90 ./MATH/bcoef.f90 
mkgtr_MODULESRC = \
./PAR/m_par.f90 ./COMMON/m_constants.f90 ./KSPACE/m_struct.f90 \
./KSPACE/m_controls.f90 ./KSPACE/m_boundaries.f90 ./KSPACE/m_kklist.f90 \
./KSPACE/m_kgenwork.f90 ./KSPACE/m_tetrahedra.f90 ./KSPACE/m_controlkgen.f90 \
./COMMON/m_dimsmod.f90 ./KSPACE/m_strfacs.f90 ./COMMON/m_inpmodules.f90 \
./COMMON/m_nrixs.f90 ./ERRORMODS/m_errormod.f90 ./IOMODS/m_iofiles.f90 \
./IOMODS/m_padio.f90 ./IOMODS/m_iomod.f90 ./ERRORMODS/m_errorfile.f90 
