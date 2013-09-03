c
c
c
c =========================================================
      subroutine rp1(maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &           wave,s,amdq,apdq,num_aux)
c =========================================================
c
c     # solve Riemann problems for the 1D Euler equations using Roe's 
c     # approximate Riemann solver.  
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c     # On output, wave contains the waves, 
c     #            s the speeds, 
c     #            amdq the  left-going flux difference  A^- \Delta q
c     #            apdq the right-going flux difference  A^+ \Delta q
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routine step1, rp is called with ql = qr = q.
c
c
      implicit double precision (a-h,o-z)
      dimension   ql(meqn,1-mbc:mx+mbc)
      dimension   qr(meqn,1-mbc:mx+mbc)
      dimension   qml(meqn,1-mbc:mx+mbc)
      dimension   qmr(meqn,1-mbc:mx+mbc)
      dimension    s(mwaves,1-mbc:mx+mbc)
      dimension wave(meqn, mwaves,1-mbc:mx+mbc)
      dimension amdq(meqn,1-mbc:mx+mbc)
      dimension apdq(meqn,1-mbc:mx+mbc)
      dimension auxl(num_aux,1-mbc:mx+mbc)
      dimension auxr(num_aux,1-mbc:mx+mbc)
      dimension   flux_l(meqn)
      dimension   flux_r(meqn)
c
c     # local storage
c     ---------------
      parameter (max2 = 2002)  !# assumes at most 2000 grid points with mbc=2
      dimension delta(3)
      dimension u(-1:max2),enth(-1:max2),a(-1:max2)
      logical efix
!       common /param/  gamma,gamma1
c
      data efix /.false./    !# use entropy fix for transonic rarefactions
c
c     # Compute Roe-averaged quantities:
c
      do 20 i=2-mbc,mx+mbc
         gammal = auxr(1,i-1)
         gammar = auxl(1,i)
         gamma1l = gammal - 1.0
         gamma1r = gammar - 1.0
         pinfl = auxr(2,i-1)
         pinfr = auxl(2,i)
         omel = auxr(3,i-1)
         omer = auxl(3,i)
         ! Densities
         rho_l = qr(1,i-1)
         rho_r = ql(1,i)
         ! Velocities
         ul = qr(2,i-1)/rho_l
         ur = ql(2,i)/rho_r
         ! Kinetic Energy
         ek_l = 0.5*rho_l*ul**2
         ek_r = 0.5*rho_r*ur**2
         ! Pressures
         pl = gamma1l*(qr(3,i-1) - ek_l) 
         pl = pl/(1.0 - omel*rho_l) - pinfl*gammal
         pr = gamma1r*(ql(3,i) - ek_r) 
         pr = pr/(1.0 - omer*rho_r) - pinfr*gammar

C        ! Compute left and right fluxes to simplify HLLC calculations
         flux_l(1) = rho_l*ul
         flux_l(2) = rho_l*ul**2 + pl 
         flux_l(3) = ul*(qr(3,i-1) + pl)
         ! ...
         flux_r(1) = rho_r*ur
         flux_r(2) = rho_r*ur**2 + pr 
         flux_r(3) = ur*(ql(3,i) + pr)
         ! Compute left and right speeds
!          enth_l = (qr(i-1,3) + pl)/rho_l
!          enth_r = (ql(i,3) + pr)/rho_r
!          cl = dsqrt(gamma1l*(enth_l - .5d0*ul**2))
!          cr = dsqrt(gamma1r*(enth_r - .5d0*ur**2))
         cl = dsqrt(gammal*(pl + pinfl)/rho_l)
         cr = dsqrt(gammar*(pr + pinfr)/rho_r)

         ! Compute the speed of left and right HLLC wave
         Sl = min(ul - cl,ur - cr) ! u(i) - a(i)
         Sr = max(ul + cl,ur + cr) ! u(i) + a(i),
         s(1,i) = 1.0*Sl
         s(3,i) = 1.0*Sr

         ! Compute HLLC middle speed state (see research notebook)
         Sm = pr - pl + rho_r*ur*(ur - Sr) - rho_l*ul*(ul - Sl) 
         Sm = Sm/(rho_r*(ur - Sr) - rho_l*(ul - Sl))
         s(2,i) = 1.0*Sm

         ! Force zero speed at contact discontinuity
         if (gammal .ne. gammar ) then
!           Sm = 0.0
          s(2,i) = 0.0
         end if
         
         ! Compute middle state pressure (both formulas yield the same value except maybe on interface)
         pstar_l = pl + rho_l*(ul - Sm)*(ul - Sl)
         pstar_r = pr + rho_r*(ur - Sm)*(ur - Sr)

!          ! Calculate ql* and qr* HLLC middle states (without pressure term)
!          do j=1,3
! 	  qml(i,j) = (flux_l(j) - Sl*qr(i-1,j))/(Sm - Sl) 
!           qmr(i,j) = (flux_r(j) - Sr*ql(i,j))/(Sm - Sr)
! 	 end do
!          ! Add pressure to momentum ones
! 	 qml(i,2) = qml(i,2) - pstar_l/(Sm - Sl)
!          qmr(i,2) = qmr(i,2) - pstar_r/(Sm - Sr)
!          ! Add pressure to energy ones
! 	 qml(i,3) = qml(i,3) - Sm*pstar_l/(Sm - Sl)
!          qmr(i,3) = qmr(i,3) - Sm*pstar_r/(Sm - Sr)

         ! Calculate ql* and qr* HLLC middle states (without pressure term)
         do j=1,3
           qml(j,i) = rho_l*(Sl - ul)/(Sl - Sm)
           qmr(j,i) = rho_r*(Sr - ur)/(Sr - Sm)
         end do
         ! Add second terms to momentum one (see Toro pg. 325)
         qml(2,i) = Sm*qml(2,i)
         qmr(2,i) = Sm*qmr(2,i)
         ! Add second terms to energy one (see Toro pg. 325)
         qml(3,i) = qml(3,i)*(qr(3,i-1)/rho_l + 
     & (Sm - ul)*(Sm + pl/(rho_l*(Sl - ul))))
         qmr(3,i) = qmr(3,i)*(ql(3,i)/rho_r + 
     & (Sm - ur)*(Sm + pr/(rho_r*(Sr - ur))))

c        # Compute the 3 waves.
c        j index over q variables
         do j=1,3
           q_l = qr(j,i-1)
           q_r = ql(j,i)
           wave(j,1,i) = qml(j,i) - q_l
           wave(j,2,i) = qmr(j,i) - qml(j,i)
           wave(j,3,i) = q_r - qmr(j,i) 
         end do

   20    continue

c!          if (gammal .ne. gammar ) then
! 	  Sm = 0.0
!           s(i,2) = 0.0
! 	 end if
c     # compute Godunov flux f0:
c     --------------------------
c
c
      if (efix) go to 110
c
c     # no entropy fix
c     ----------------
c
c     # amdq = SUM s*wave   over left-going waves
c     # apdq = SUM s*wave   over right-going waves
c
      do 100 m=1,3
         do 100 i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do 90 mw=1,mwaves
               if (s(mw,i) .lt. 0.d0) then
                   amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                 else
                   apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                 endif
   90          continue
  100       continue
      go to 900
c
c-----------------------------------------------------
c
  110 continue


c
  900 continue
      return
      end
