C     -*- fortran -*-
C     This file is autogenerated with f2py (version:2)
C     It contains Fortran 77 wrappers to fortran functions.

      subroutine f2pywrapfmohl (fmohlf2pywrap, a, b, q, np)
      external fmohl
      double precision a
      double precision b
      double precision q
      integer np
      double precision fmohlf2pywrap, fmohl
      fmohlf2pywrap = fmohl(a, b, q, np)
      end


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


      subroutine f2pywraprds (rdsf2pywrap, x, y, z)
      external rds
      double precision x
      double precision y
      double precision z
      double precision rdsf2pywrap, rds
      rdsf2pywrap = rds(x, y, z)
      end

