      integer mxedgs, maxpts, fullpts, lmax, nexts, nlog, numcnv,vfeff
      real econv, ecnv, bigemax, bigemin
      character*8 infile
      real khi, klo, delta, minres, trsize, xkstep

      parameter(infile = 'fullspectrum.inp')
      parameter ( lmax = 3)
!     parameter ( fullpts = 100001)
      parameter ( fullpts = 200001)
      parameter(mxedgs = 26)
!     parameter(maxpts = 602)
      parameter(maxpts = 6002)
      parameter(nexts = 2)
      parameter ( nlog = 1000)
      parameter ( econv = 1.83746545) !50/hart, max. edge onset for conv.
      parameter(numcnv = 2400) !no. of points on linear grid for convolutions
      parameter(ecnv = 3.6749309) !100/hart energy at which convolution stops 
      parameter(khi = 4.0)
      parameter(klo = 3.0)
!     parameter(minres=0.0367493) !1.0/hart smallest energy step for fprime grids.
      parameter(minres=0.015)
      !factor by which coarsness of fprime grids differ in consecutive files
      parameter(delta=10.0)  
      !length of interval to use in the smooth transition at end of exafs
      !region in hartrees
      parameter(trsize=0.05)
!     parameter(trsize=0.367)
!     parameter(trsize=1.0)
      !energy bounds for sumrule integration
!     parameter (bigemin=0.000367493090027428) !bigemin=0.01/hart
      parameter (bigemin=0.0)
!     parameter (bigemax=3674.93090027428) !bigemax=100000.0/hart
!     parameter (bigemax=5512.39635041142) !bigemax=150000.0/hart
      parameter (bigemax=18383.0) !about 0.5 GeV
!     parameter (xkstep=0.01)
      parameter (xkstep=0.005)
