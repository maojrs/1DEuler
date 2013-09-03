c
c
c =========================================================
       subroutine qinit(meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c     # Smooth entropy wave hitting a shock
c
c
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc)
      dimension p(1-mbc:mx+mbc)
      dimension aux(maux,1-mbc:mx+mbc)
      dimension dist(500), pressure(500)
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ omegas,omeplas,omewat
      common /param/ rhog,rhop,rhow
c
c     
      ! Open file for pressure Gauge plot
      open (22,file="a-pgauge.dat",action="write",
     & status="replace")
      open (23,file="a-pIC.dat",action="read",status="old")
!       
      ! Read pressure data from file
      do i=1,500
        read(23,*) dist(i), pressure(i)
        !print*, dist(i),pressure(i)
      enddo


      ! Write initial conditions
      rhog = 1.0 !kg/m^3
      rhop = 1050.0 !kg/m^3
      rhow = 1000.0 !kg/m^3
      p0 = 101325.0 !
      p = p0
      c0 = sqrt(gammagas*p0/rhog) 
      Egas0 = p0/(gammagas - 1.d0)  
      do 150 i=1,mx
         xcell = xlower + (i-0.5d0)*dx

         ! Look for correct value for pressure in data file
         ddx = dist(101) - dist(100)
         do j=1,500
          dist2 = (dist(j) - 16.0) 
          if (abs(dist2 - xcell) <  ddx/2) then
            p(i) = 6894.75729*pressure(j) + p0
!           exit
          end if
         end do

         if (aux(1,i) == gammagas) then
          q(1,i) = rhog*(p(i)/p0)**(1/gammagas) !+ 5.0d0*dexp(-200.d0*(xcell+1.0)**2)
          q(2,i) = (2/(gammagas - 1.0))*(-c0 + 
     & sqrt(gammagas*p(i)/q(1,i)))
!          q(2,i) = 0.d0
          q(3,i) = p(i)/(aux(1,i) - 1.0) + q(2,i)**2/(2.0*q(1,i))
         else if (aux(1,i) == gammaplas) then
          q(1,i) = rhop
          q(2,i) = 0.0
          ! make sure pressure jump is zero across interface using SGEOS
          q(3,i) = gammaplas*pinfplas - gammagas*pinfgas
          q(3,i) = q(3,i) + Egas0*(gammagas - 1.d0)/(1.d0 - omegas*rhog)
          q(3,i) = q(3,i)*(1.d0 - omeplas*rhop)/(gammaplas - 1.d0)
          Eplas0 = 1.0*q(3,i) ! Needed to compute energy in water state 
         else
          q(1,i) = rhow
          q(2,i) = 0.0
          ! Make sure jump in pressure is zero again
          q(3,i) = -gammaplas*pinfplas + gammawat*pinfwat
          q(3,i) = q(3,i)+Eplas0*(gammaplas - 1.d0)/(1.d0 - omewat*rhow)
          q(3,i) = q(3,i)*(1.d0 - omeplas*rhop)/(gammawat - 1.d0)
         end if
         
  150    continue
c
      return
      end

