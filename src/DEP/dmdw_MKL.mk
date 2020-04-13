dmdwSRC = \
./DMDW/dmdw.f90 ./PAR/par.f90 ./COMMON/wlog.f90 \
./COMMON/str.f90 ./COMMON/str2dp.f90 ./COMMON/rdhead.f90 \
./COMMON/chopen.f90 
dmdw_MODULESRC = \
./PAR/m_par.f90 ./COMMON/m_kinds.f90 ./DMDW/m_const_and_conv.f90 \
./INPGEN/m_strings.f90 ./INPGEN/m_ptable.f90 ./DMDW/m_math.f90 \
./DMDW/m_dmdw.f90 ./COMMON/m_dimsmod.f90 ./KSPACE/m_controls.f90 \
./KSPACE/m_struct.f90 ./KSPACE/m_kklist.f90 ./KSPACE/m_strfacs.f90 \
./COMMON/m_constants.f90 ./COMMON/m_inpmodules.f90 ./ERRORMODS/m_errorfile.f90 
