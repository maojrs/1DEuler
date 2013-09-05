c     ============================================
      subroutine b4step1(mbc,mx,meqn,q,
     &            xlower,dx,t,dt,maux,aux)
c     ============================================
c
c     # called from claw1 before each call to step1.
c     # use to set time-dependent aux arrays or perform other tasks
c     # which must be done every time step.
c
c
c     
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
      dimension P(1-mbc:mx+mbc)
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ omegas,omeplas,omewat
      common /param/ rhog,rhop,rhow
c
      ! GAUGE for pressure (choose point to see pressure as a function of time)
      ! Choose x value to obtain pressure as a function of time
      xcell = -2.0
      i = floor((xcell - xlower)/dx + 0.5)

      ! Calculate pressure at point xcell
      gamma = aux(1,i)
      gamma1 = aux(1,i) - 1.0
      pinf = aux(2,i)
      omega = aux(3,i)
      rho = q(1,i)           ! density
      mom = q(2,i)           ! momentum
      ene = q(3,i)           ! energy
      P = gamma1*(ene - 0.5*mom*mom/rho)/(1.0 - omega*rho)
      P = P - gamma*pinf 
     
      ! Write Gauge data to file
      write (22,*) t, P
      
      ! Save solution on aux array
!       j = t/dt
!       if (j < mx+mbc) then
! 	aux(3,j+1) = P
!       end if

      return
      end

