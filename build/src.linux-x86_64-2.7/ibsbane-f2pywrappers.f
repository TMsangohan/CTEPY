C     -*- fortran -*-
C     This file is autogenerated with f2py (version:2)
C     It contains Fortran 77 wrappers to fortran functions.

      subroutine f2pywrapgbane (gbanef2pywrap, alpha, gtab, nx, ny
     &)
      external gbane
      double precision alpha
      integer nx
      integer ny
      double precision gtab(nx,ny)
      double precision gbanef2pywrap, gbane
      gbanef2pywrap = gbane(alpha, gtab, nx, ny)
      end

