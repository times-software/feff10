eelsSRC = \
./EELS/eels.f90 ./EELS/wavelength.f90 ./PAR/parallel.f90 \
./COMMON/wlog.f90 ./COMMON/str.f90 ./COMMON/chopen.f90 \
./EELS/readsp.f90 ./EELS/concat.f90 ./COMMON/rdhead.f90 \
./EELS/calculateweights.f90 ./EELS/angularmesh.f90 ./EELS/qmesh.f90 \
./EELS/euler.f90 ./EELS/productmatvect.f90 ./EELS/writeangulardependence3.f90 \
./EELS/writeangulardependence2.f90 
eels_MODULESRC = \
./COMMON/m_constants.f90 ./PAR/m_par.f90 ./COMMON/m_dimsmod.f90 \
./KSPACE/m_controls.f90 ./KSPACE/m_struct.f90 ./KSPACE/m_kklist.f90 \
./KSPACE/m_strfacs.f90 ./COMMON/m_inpmodules.f90 ./EELS/m_spectrum.f90 \
./EELS/m_qvectors.f90 ./EELS/m_work.f90 ./EELS/m_program_control.f90 \
./ERRORMODS/m_errorfile.f90 
