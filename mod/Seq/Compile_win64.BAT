set comp="ifort.exe"
%comp% %flags% rdinp_tot.f90 -o ..\win64\rdinp
%comp% %flags% dmdw_tot.f90 -o ..\win64\dmdw
%comp% %flags% dym2feffinp_tot.f90 -o ..\win64\dym2feffinp
%comp% %flags% atomic_tot.f90 -o ..\win64\atomic
%comp% %flags% opconsat_tot.f90 -o ..\win64\opconsat
%comp% %flags% eels_tot.f90 -o ..\win64\eels
%comp% %flags% compton_tot.f90 -o ..\win64\compton
%comp% %flags% ff2x_tot.f90 -o ..\win64\ff2x
%comp% %flags% fms_tot.f90 -o ..\win64\fms
%comp% %flags% genfmt_tot.f90 -o ..\win64\genfmt
%comp% %flags% ldos_tot.f90 -o ..\win64\ldos
%comp% %flags% mkgtr_tot.f90 -o ..\win64\mkgtr
%comp% %flags% path_tot.f90 -o ..\win64\path
%comp% %flags% pot_tot.f90 -o ..\win64\pot
%comp% %flags% screen_tot.f90 -o ..\win64\screen
%comp% %flags% sfconv_tot.f90 -o ..\win64\sfconv
%comp% %flags% xsph_tot.f90 -o ..\win64\xsph
del *.obj
del *.mod
