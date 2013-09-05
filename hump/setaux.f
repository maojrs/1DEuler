c     ============================================
      subroutine setaux(mbc,mx,xlower,dx,maux,aux)
c     ============================================
c
c     # set auxiliary arrays 
c     # dummy routine when no auxiliary arrays
c
c     
      implicit double precision (a-h,o-z)
      dimension aux(maux,1-mbc:mx+mbc)
      common /param/ gammagas, gammaplas, gammawat
      common /param/ pinfgas,pinfplas,pinfwat
      common /param/ omegas,omeplas,omewat
      common /param/ rhog,rhop,rhow
      
! !       aux(1,i) plays the role of gamma, aux(2,i) could be another parameter
      ! Aux arrays establish values for the three free parameters in polystyrene EOS
      ! Polystyrene EOS based on Van der Waals approximation in Spender and Gilmore paper
      
      pwidth = 1.3 !0.0 !0.1 ! 0.05 !0.1 !0.3 ! 0.7 !1.3
      do i=1-mbc,mx + mbc
!         indexi = i + mbc
        xcell = xlower + (i-0.5d0)*dx
        if (abs(xcell) < pwidth) then
          aux(1,i) = gammaplas
          aux(2,i) = pinfplas
          aux(3,i) = omeplas
        else if (xcell < -pwidth) then
          aux(1,i) = gammagas
          aux(2,i) = pinfgas
          aux(3,i) = omegas
        else
          aux(1,i) = gammawat
          aux(2,i) = pinfwat
          aux(3,i) = omewat
        end if
      end do

!       do i=1-mbc,mx + mbc - 1
! 	if (aux(1,i) .ne. aux(1,i+1)) then
! 	    aux(1,i) = 0.5d0*(gamma + gamma2)
! 	end if
!       end do

!       ! Smooth varying gamma
!       gamma = 1.4
!       gamma2 = 3.0
!       do i=1-mbc,mx + mbc
!         indexi = i + mbc
! 	xcell = xlower + (indexi-0.5d0)*dx
!         aux(2,i) = 0.0
!         aux(1,i) = gamma + (gamma2 - gamma)*ATAN(10*xcell)/3.1416
!       end do


c
      return
      end


