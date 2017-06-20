!     =========================================================================
!     subroutine :: pipiscat
!     -------------------------------------------------------------------------
!     Summary :: Code produces pion-pion scattering phase shifts and
!                inelasticities up to 1.42 GeV based on KPY III [1]
!
!     Input  :: s     - Center-of-Momentum energy squared [GeV^2]
!               wav   - Specifies LI-wave (L = Angular momentum, I = Isospin)
!                       wav = { 'S0', 'S2', 'P1', 'D0', 'D2', 'F1' }
!     Output :: delta - Phase Shift for LI-wave 
!               eta   - Inelasticity for LI-wave
!
!     Type declarations - real (8)      :: s, delta, eta
!                       - character (2) :: wav
!
!     To call :: call pipiscat (s,wav,delta,eta)
!     
!     Parameterization taken from Pion-Pion Scattering III
!     [1] R. Kamiński, J. R. Peláez, and F. J. Ynduráin, 
!         Phys. Rev. D 77, 054015 (2008).
!     Parameters represent Appendix A: (UFD) Set
!     
!     Author      :: Andrew W. Jackura
!     Affiliation :: Joint Physics Analysis Center (JPAC)
!     Contact     :: ajackura@indiana.edu
!     Last Update :: Last Update :: 20/06/2016   
!     =========================================================================
      subroutine pipiscat (s,wav,delta,eta)
!     -------------------------------------------------------------------------
      implicit none
      real (8)               :: s, delta, eta
      character (2)          :: wav
      real (8), parameter    :: xmpi = 0.1396d0, xmK = 0.4936d0
      real (8), parameter    :: pi = 4.d0*datan(1.d0)
      real (8), parameter    :: sth0 = 4.d0 * xmpi * xmpi
      real (8), parameter    :: sth1 = 0.932d0**2
      real (8), parameter    :: sth2 = 1.420d0**2
      real (8), parameter    :: sthK = 4.d0 * xmK * xmK
      complex (8), parameter :: xr = (1.d0,0.d0)
      complex (8), parameter :: xi = (0.d0,1.d0)
!     Parameter declarations common to all waves
      real (8) :: s0, B0, B1, B2, cot_delta, xkmom, xkmom1, xkmom2, w
      real (8) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
      real (8), external :: mom_el, w_conf
!     Parameter declaration for S0
      real (8) :: z0, alpha1, alpha2, beta1, beta2
      real (8) :: gamma11, gamma12, gamma22, xM1, xM2, mu, tan_delta
      real (8) :: detK, K11, K12, K21, K22
!     Parameter declaration for S2
      real (8) :: z2, sl, sh, sm, Bh0, Bh1, Bh2, wl, wh, wlsm, whsm, eps
      real (8) :: derv
!     Parameter declaration for P1
      real (8) :: xmrho, lambda0, lambda1, lambda2, eps1, eps2
!     Parameter declaration for D0
      real (8) :: xmf2, sqxmf2, r, xkmomK
!     Parameter declaration for D2
      real (8) :: del
!     Parameter declaration for F1
      real (8) :: lambda
!     -------------------------------------------------------------------------
      delta = 0.d0
      eta   = 1.d0
      if ( s .lt. sth0 ) then
         delta = 0.d0
         return
      endif
      xkmom = mom_el (s,sth0)
!     ----------------------
!     Wave Selection
!     ----------------------
      select case(wav)
!     ----------------------
!     case ('S0')
!     ----------------------
      case ('S0')
!
         if ( ( s .le. sth1 ) .and. ( s .gt. sth0 ) ) then
            s0 = 4.d0 * xmK**2
            z0 = xmpi
            w  = w_conf (s,s0)
!            z0 = xmpi
!            B0 = 4.3d0
!            B1 = -26.7d0
!            B2 = -14.1d0
            B0 = 4.41d0
            B1 = -26.25d0
            B2 = -15.8d0
            z0 = .1661d0
            tmp1 = ( dsqrt(s) / 2.d0 / xkmom ) 
            tmp2 = ( xmpi * xmpi / ( s - 0.5d0 * z0 * z0 ) )
            tmp3 = z0 * z0 / xmpi / dsqrt (s) + B0 + B1 * w + B2 * w * w
            cot_delta = tmp1 * tmp2 * tmp3
            delta = atan2 (1.d0,cot_delta)
         elseif ( ( s .le. sth2 ) .and. ( s .gt. sth1 ) ) then
!     
            alpha1  = 0.843d0
            alpha2  = 0.2d0
            beta1   = 1.02d0
            beta2   = 1.33d0
            gamma11 = 3.10d0
            gamma12 = 1.82d0
            gamma22 = -7.00d0
            xM1     = 0.888d0
            xM2     = 1.327d0
            mu      = 1.d0
            K11     = mu * alpha1 * alpha1 / ( xM1 * xM1 - s ) 
     c                + mu * beta1 * beta1 / ( xM2 * xM2 - s )
     c                + gamma11 / mu
            K12     = mu * alpha1 * alpha2 / ( xM1 * xM1 - s )
     c                + mu * beta1 * beta2 / ( xM2 * xM2 - s )
     c                + gamma12 / mu
            K22     = mu * alpha2 * alpha2 / ( xM1 * xM1 - s )
     c                + mu * beta2 * beta2 / ( xM2 * xM2 - s )
     c                + gamma22 / mu
            K21     = K12
            detK    = K11 * K22 - K12 * K21
            xkmom1  = mom_el (s,sth0)
!            xkmom2  = mom_el (s,sthK)
            if ( s .lt. sthK ) then
               xkmom2 = 0.5d0 * dsqrt ( sthK - s )
               tmp1 = xkmom1 * abs (xkmom2) * detK + xkmom1 * K11
               tmp2 = 1.d0 + abs (xkmom2) * K22
               tan_delta = tmp1 / tmp2
!               delta = datan2(tmp1,tmp2)
            elseif ( s .ge. sthK ) then
               xkmom2 = mom_el (s,sthK)
               tmp1 = 2.d0 * xkmom1 * (K11 + xkmom2**2 * K22 * detK)
               tmp2 = xkmom1**2 * K11**2 - xkmom2**2 * K22**2 
     c                + xkmom1**2 * xkmom2**2 * detK**2 - 1.d0
               tmp3 = xkmom1**2* K11**2 + xkmom2**2 * K22**2
     c                + xkmom1**2 * xkmom2**2 *detK**2 + 1.d0
               tmp4 = 4.d0 * xkmom1**2 * xkmom2**2 * K12**4
               tmp5 = dsqrt ( tmp3 * tmp3 - tmp4 )
               tan_delta = ( tmp2 + tmp5 ) / tmp1
!               delta = datan2(tmp2 + tmp5 , tmp1)
!     Inelasticity
               tmp6 = ( 1.d0 + xkmom1 * xkmom2 * detK )**2
               tmp7 = ( xkmom1 * K11 - xkmom2 * K22 )**2
               tmp8 = ( 1.d0 - xkmom1 * xkmom2 * detK )**2
               tmp9 = ( xkmom1 * K11 + xkmom2 * K22 )**2
               eta  = dsqrt( (tmp6 + tmp7) / (tmp8 + tmp9) )
            endif
!            cot_delta = 1./tan_delta
!            delta = atan2(1.d0,cot_delta)
            delta = atan(tan_delta)
            if  (s.lt.1.2405**2) then
               delta = delta + pi;
            else
               delta = delta + 2*pi;
            endif
         else 
            delta = 0.d0
            eta   = 0.d0
         endif
!     
!     ----------------------
!     case ('S2')
!     ----------------------
      case ('S2')
         if ( ( s .le. sth1 ) .and. ( s .gt. sth0 ) ) then
            z2   = xmpi
            sl   = 1.05d0**2
            B0   = -80.4d0
            B1   = -73.6d0
            w    = w_conf (s,sl)
            tmp1 = dsqrt (s) / 2.d0 / xkmom
            tmp2 = xmpi * xmpi / ( s - 2.d0 * z2 * z2 )
            tmp3 = B0 + B1 * w
            cot_delta = tmp1 * tmp2 * tmp3
         elseif ( ( s .le. sth2 ) .and. ( s .gt. sth1 ) ) then
            sl   = (1.050d0)**2
            sh   = (1.450d0)**2
            sm   = 0.932d0
            wl   = w_conf (s,sl)
            wh   = w_conf (s,sh)
            wlsm = w_conf (sm,sl)
            whsm = w_conf (sm,sh)
            B0   = -80.4d0
            B1   = -73.6d0
            Bh0  = B0 + B1 * wlsm
            derv = 2.50814d0
            Bh1  = B1 * derv
            Bh2  = 112.d0
            tmp1 = sqrt (s) / 2.d0 / xkmom
            tmp2 = xmpi * xmpi / ( s - 2.d0 * xmpi * xmpi )
            tmp3 = Bh0 + Bh1 * (wh - whsm)
            tmp4 = Bh2 * ( wh - whsm )**2
            tmp5 = tmp3 + tmp4
            cot_delta = tmp1 * tmp2 * tmp5
            if ( s.gt. sl ) then
               eps  = 0.17d0
               eta  =  1.d0 - eps * ( 1.d0 - sl / s )**(1.5d0)
            endif
         endif
         delta = atan2(1.d0,cot_delta)
!     ----------------------
!     case ('P1')
!     ----------------------
      case ('P1')
         if ( ( s .gt. sth0 ) .and. ( s .le. sthK ) ) then
            s0    = 1.05d0**2
            B0    = 1.055d0
            B1    = 0.15d0
            xmrho = 0.7736d0
            w     = w_conf (s,s0)
            tmp1  = dsqrt (s) / 2.d0 / xkmom / xkmom / xkmom
            tmp2  = xmrho * xmrho - s
            tmp3  = 2.d0 * xmpi**3 / xmrho / xmrho / dsqrt (s)
            tmp4  = tmp3 + B0 + B1 * w         
            cot_delta = tmp1 * tmp2 * tmp4
            delta = atan2 (1.d0,cot_delta)
         elseif ( ( s .gt. sthK ) .and. ( s .le. sth2 ) ) then
            lambda0 = 2.681d0
            lambda1 = 1.57d0
            lambda2 = -1.96d0
            eps1    = 0.10d0
            eps2    = 0.11d0
            tmp1    = dsqrt ( s / sthK ) - 1.d0
            xkmom2  = dsqrt ( 1.d0 - sthK / s )
            delta   = lambda0 + lambda1 * tmp1 + lambda2 * tmp1 * tmp1
            eta     = 1.d0 - eps1 * xkmom2 - eps2 * xkmom2 * xkmom2
         endif
!     ----------------------
!     case ('D0')
!     ----------------------
      case ('D0')
         if ( ( s .le. sthK ) .and. ( s .gt. sth0 ) ) then
            s0   = 1.05d0**2
            B0   = 12.47d0
            B1   = 10.12d0
            xmf2 = 1.2754d0
            w    = w_conf (s,s0)
            tmp1 = dsqrt (s) / 2.d0 / xkmom**5
            tmp2 = xmf2 * xmf2 - s
            tmp3 = B0 + B1 * w            
            cot_delta = tmp1 * tmp2 * xmpi * xmpi * tmp3
         elseif ( ( s .gt. sthK ) .and. ( s .le. sth2) ) then
            xmf2 = 1.2754d0
            sh   = 1.45d0**2
            w    = w_conf (s,sh)
            Bh0  = 18.77d0
            Bh1  = 43.7d0
            tmp1 = dsqrt (s) / 2.d0 / xkmom**5
            tmp2 = xmf2 * xmf2 - s
            tmp3 = Bh0 + Bh1 * w
            cot_delta = tmp1 * tmp2 * xmpi * xmpi * tmp3
         
            eps  = 0.284d0
            r    = 2.54d0
            tmp4 = ( 1.d0 - sthK / s )**(2.5d0)
            tmp5 = ( 1.d0 - sthK / xmf2 / xmf2 )**(-2.5d0)
            xkmom2 = mom_el (s,sthK)
            sqxmf2 = xmf2 * xmf2
            xkmomK = mom_el (sqxmf2,sthK)
            tmp6   = 1.d0 + r * ( 1.d0 - xkmom2 / xkmomK )
            eta = 1.d0 - eps * tmp4 * tmp5 * tmp6
         endif
         delta = atan2 (1.d0,cot_delta)
!     ----------------------
!     case ('D2')
!     ----------------------
      case ('D2')
         if ( ( s .gt. sth0 ) .and. ( s .le. sth2 ) ) then
            B0   = 2.4d3
            B1   = 7.8d3
            B2   = 23.7d3
            Del  = 0.196d0
            s0   = 1.45d0**2
            w    = w_conf (s,s0)
            tmp1 = dsqrt (s) / 2.d0 / xkmom**5
            tmp2 = B0 + B1 * w + B2 * w * w
            tmp3 = xmpi**4 * s / ( 4.d0 * ( xmpi**2 + del**2 ) - s )
            cot_delta = tmp1 * tmp2 * tmp3
            delta = atan (1.d0/cot_delta)
            if ( s .gt. 1.05d0 ) then
               sl  = 1.05d0**2
               eps = 0.2d0
               eta = 1.d0 - eps * ( 1.d0 - sl / s )**3
            endif
         endif


!     ----------------------
!     case ('F1')
!     ----------------------
      case ('F1')
         if ( ( s .gt. sth0 ) .and. ( s .le. sth2 ) ) then
            s0     = 1.45d0**2
            B0     = 1.09d5
            B1     = 1.41d5
            lambda = 0.051d5
            w      = w_conf (s,s0)
            tmp1   = dsqrt (s) / 2.d0 / xkmom**7
            tmp2   = 2.d0 * lambda * xmpi / dsqrt (s)
            tmp3   = B0 + B1 * w
            tmp4   = tmp2 + tmp3
            cot_delta = tmp1 * xmpi**6 * tmp4
         endif
         delta = atan (1.d0 / cot_delta)
!     ----------------------
!     case default
!     ----------------------
      case default
         delta = 0.d0
      end select
      return
      end subroutine pipiscat

!     =========================================================================
!     function :: w_conf
!     -------------------------------------------------------------------------
!     Summary  :: Returns conformal map given 's' and branching point 's0'
!
!     =========================================================================
      real (8) function w_conf (s,s0)
      implicit none
      real (8) :: s, s0, tmp1, tmp2
      tmp1   = dsqrt ( s ) - dsqrt ( s0 - s )
      tmp2   = dsqrt ( s ) + dsqrt ( s0 - s )
      w_conf = tmp1 / tmp2
      return
      end function w_conf

!     =========================================================================
!     function :: mom_el
!     -------------------------------------------------------------------------
!     Summary  :: Returns momentum for elastic scattering given 's' and 'sth'
!                 where sth = 4 m ^2 is the threshold 
!
!     =========================================================================
      real (8) function mom_el (s,sth)
      implicit none
      real (8) :: s, sth
      mom_el = dsqrt (s - sth) / 2.d0
      return
      end function mom_el
