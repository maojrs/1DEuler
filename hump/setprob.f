      subroutine setprob
      implicit double precision (a-h,o-z)
      character*12 fname
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ omegas,omeplas,omewat
      common /param/ rhog,rhop,rhow
c
c     # Set gamma and gamma1 = gamma-1 for Euler equations
c     # Passed to the Riemann solver rp1.f in a common block
c
c
      iunit = 7
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)

      read(7,*) gammagas
      read(7,*) gammaplas
      read(7,*) gammawat
      read(7,*) pinfgas
      read(7,*) pinfplas
      read(7,*) pinfwat
      read(7,*) omegas
      read(7,*) omeplas
      read(7,*) omewat
      read(7,*) rhog
      read(7,*) rhop
      read(7,*) rhow
      

      return
      end

