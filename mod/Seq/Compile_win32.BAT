set comp="ifort.exe"
%comp% %flags% rdinp_tot.f90 -o ..\win32\rdinp
%comp% %flags% dmdw_tot.f90 -o ..\win32\dmdw
%comp% %flags% dym2feffinp_tot.f90 -o ..\win32\dym2feffinp
%comp% %flags% atomic_tot.f90 -o ..\win32\atomic
%comp% %flags% opconsat_tot.f90 -o ..\win32\opconsat
%comp% %flags% eels_tot.f90 -o ..\win32\eels
%comp% %flags% compton_tot.f90 -o ..\win32\compton
%comp% %flags% ff2x_tot.f90 -o ..\win32\ff2x
%comp% %flags% fms_tot.f90 -o ..\win32\fms
%comp% %flags% genfmt_tot.f90 -o ..\win32\genfmt
%comp% %flags% ldos_tot.f90 -o ..\win32\ldos
%comp% %flags% mkgtr_tot.f90 -o ..\win32\mkgtr
%comp% %flags% path_tot.f90 -o ..\win32\path
%comp% %flags% pot_tot.f90 -o ..\win32\pot
%comp% %flags% screen_tot.f90 -o ..\win32\screen
%comp% %flags% sfconv_tot.f90 -o ..\win32\sfconv
%comp% %flags% xsph_tot.f90 -o ..\win32\xsph
del *.obj
del *.mod
