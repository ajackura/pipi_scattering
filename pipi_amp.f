!     =========================================================================
!     subroutine :: pipi_pw_amp
!     -------------------------------------------------------------------------
!     Summary :: Code produces pion-pion scattering partial wave amplitudes
!                given phase shifts and inelasticities.
!                Normalization based on GKPRY IV [1]
!
!     Input  :: s   - Center-of-Momentum energy squared [GeV^2]
!               wav - Specifies LI-wave (L = Angular momentum, I = Isospin)
!                     wav = { 'S0', 'S2', 'P1', 'D0', 'D2', 'F1' }
!     Output :: amp - complex LI partial wave amplitude
!
!     Type declarations - real (8)      :: s
!                       - character (2) :: wav
!                       - complex (8)   :: amp
!
!     To call :: call pipi_pw_amp (s,wav,amp)
!     
!     t_{L,I}(s) = ( sqrt(s) / 2 k ) f_{L,I}(s)
!     
!     f_{L,I}(s) = ( eta_{L,I}(s) exp ( 2 i delta_{L,I}(s) ) - 1 ) / 2i
!
!     Amplitude normalization taken from Eq. (1) from Pion-Pion Scattering IV
!     [1] R. García-Martín, R. Kamiński, J. R. Peláez,
!         J. Ruiz de Elvira, and F. J. Ynduráin, 
!         Phys. Rev. D 83, 074004 (2011).
!     
!     Author      :: Andrew W. Jackura
!     Affiliation :: Joint Physics Analysis Center (JPAC)
!     Contact     :: ajackura@indiana.edu
!     Last Update :: 20/06/2016
!     =========================================================================
      subroutine pipi_pw_amp (s,wav,amp)
      implicit none
      real (8)               :: s
      complex (8)            :: amp
      character(len=2)       :: wav
      complex (8), parameter :: xi = (0.d0,1.d0)
      real (8)               :: delta, eta, xkmom
      real (8), parameter    :: pi = 4.d0*datan(1.d0)
      real (8), parameter    :: xmpi = 0.1396d0
      call pipiscat (s,wav,delta,eta)
      xkmom = 0.5d0 * sqrt(s - 4.d0*xmpi**2)
      amp   = ( eta * cdexp (2*xi*delta) - 1.d0 ) / ( 2*xi )
      amp   = sqrt(s) * amp / 2.d0 / xkmom      
      return
      end subroutine pipi_pw_amp
!     =========================================================================
!     subroutine :: pipi_amp
!     -------------------------------------------------------------------------
!     Summary :: Code produces pion-pion scattering amplitudes given
!                phase shifts and inelasticities based on GKPRY IV [1]
!
!     Input  :: s     - Center-of-Momentum energy squared [GeV^2]
!               theta - Scattering angle in Center-of-Momentum frame
!               Iso   - Isospin quantum number
!                     
!     Output :: amp   - complex isospin  amplitude
!
!     Type declarations - real (8)      :: s, theta
!                       - integer       :: Iso
!                       - complex (8)   :: amp
!
!     To call :: call pipi_amp (s,theta,Iso,amp)
!     
!     A^{I} (s,z) = 32 pi Sum_{l = 0}^{\infty} ( 2 l + 1 ) t_{L,I} P_{L} (z)
!
!     Amplitude normalization taken from Eq. (1) from Pion-Pion Scattering IV
!     [1] R. García-Martín, R. Kamiński, J. R. Peláez,
!         J. Ruiz de Elvira, and F. J. Ynduráin, 
!         Phys. Rev. D 83, 074004 (2011).
!     
!     Author      :: Andrew W. Jackura
!     Affiliation :: Joint Physics Analysis Center (JPAC)
!     Contact     :: ajackura@indiana.edu
!     Last Update :: 20/06/2016
!     =========================================================================
      subroutine pipi_amp (s,theta,Iso,amp)
      implicit none
      integer                         :: Iso, i
      real (8)                        :: s, theta, z
      complex (8)                     :: amp, amp1, amp2
      integer, parameter              :: nset = 6
      real (8), dimension (0:9)       :: pol, pol1p
      complex (8), dimension (nset)   :: pw_amp
      character (2), dimension (nset) :: waveset
      real (8), parameter             :: xmpi = 0.1396d0
      real (8), parameter             :: pi = 4.d0 * datan (1.d0)
      waveset (1) = 'S0'
      waveset (2) = 'S2'
      waveset (3) = 'P1'
      waveset (4) = 'D0'
      waveset (5) = 'D2'
      waveset (6) = 'F1'
      do i = 1, nset
         call pipi_pw_amp (s,waveset(i),pw_amp(i))
      enddo
      z = dcos (theta)
      call legendreP (z,pol,pol1p)
      if     ( Iso .eq. 0 ) then
         amp1 = ( 2.d0 * 0.d0 + 1.d0 ) * pw_amp (1) * pol (0)
         amp2 = ( 2.d0 * 2.d0 + 1.d0 ) * pw_amp (4) * pol (2)
         amp  = amp1 + amp2
      elseif ( Iso .eq. 1 ) then
         amp1 = ( 2.d0 * 1.d0 + 1.d0 ) * pw_amp (3) * pol (1)
         amp2 = ( 2.d0 * 3.d0 + 1.d0 ) * pw_amp (6) * pol (3)
         amp  = amp1 + amp2
      elseif ( Iso .eq. 2 ) then
         amp1 = ( 2.d0 * 0.d0 + 1.d0 ) * pw_amp (2) * pol (0)
         amp2 = ( 2.d0 * 2.d0 + 1.d0 ) * pw_amp (5) * pol (2)
         amp  = amp1 + amp2
      endif
      amp = 32.d0 * pi * amp
      return
      end subroutine pipi_amp

!     =========================================================================
!     subroutine :: pipi_amp_charged
!     -------------------------------------------------------------------------
!     Summary :: Code produces pion-pion scattering amplitudes for 
!                pi^+ pi^- --> pi^+ pi^-
!
!     Input  :: s     - Center-of-Momentum energy squared [GeV^2]
!               theta - Scattering angle in Center-of-Momentum frame 
!               
!     Output :: amp   - complex amplitude for charged scattering
!
!     Type declarations - real (8)      :: s, theta
!                       - complex (8)   :: amp
!
!     To call :: call pipi_amp_charged (s,theta,amp)
!     
!     Amplitude for :: pi^+ pi^- --> pi^+ pi^-
!
!     A (s,z) = sum_{I=0}^{2} < I 0 | 1 + 1 - > A^{I} (s,z)
!
!     Amplitude normalization taken from Eq. (1) from Pion-Pion Scattering IV
!     [1] R. García-Martín, R. Kamiński, J. R. Peláez,
!         J. Ruiz de Elvira, and F. J. Ynduráin, 
!         Phys. Rev. D 83, 074004 (2011).
!     
!     Author      :: Andrew W. Jackura
!     Affiliation :: Joint Physics Analysis Center (JPAC)
!     Contact     :: ajackura@indiana.edu
!     Last Update :: 20/06/2016
!     =========================================================================
      subroutine pipi_amp_charged (s,theta,amp)
      implicit none
      integer                      :: I
      real (8)                     :: s, theta
      complex (8)                  :: amp
      real (8), dimension (0:2)    :: clebsch
      complex (8), dimension (0:2) :: amp_pipi
      clebsch (0) = 1.d0 / 3.d0
      clebsch (1) = 1.d0 / 2.d0
      clebsch (2) = 1.d0 / 6.d0
      amp = (0.d0,0.d0)
      do I = 0,2
         call pipi_amp (s,theta,I,amp_pipi (I))
         amp = amp + clebsch (I) * amp_pipi (I)
      enddo
!      write (0,*) 'amp', s, amp
      return
      end subroutine pipi_amp_charged


!     =========================================================================
!     subroutine :: pipi_cross_section
!     -------------------------------------------------------------------------
!     Summary :: Code produces pion-pion scattering cross section
!                for the charged reaction pi^+ pi^- --> pi^+ pi^-
!
!     Input  :: s     - Center-of-Momentum energy squared [GeV^2]

!     Output :: sigma - Cross Section for charged pi^+ pi^- --> pi^+ pi^-
!
!     Type declarations - real (8)      :: s, sigma
!
!     To call :: call pipi_cross_section (s,sigma)
!     
!     sigma = Im ( A (s,z=0) ) / ( 2 sqrt(s) k )
!
!     Amplitude normalization taken from Eq. (1) from Pion-Pion Scattering IV
!     [1] R. García-Martín, R. Kamiński, J. R. Peláez,
!         J. Ruiz de Elvira, and F. J. Ynduráin, 
!         Phys. Rev. D 83, 074004 (2011).
!     
!     Author      :: Andrew W. Jackura
!     Affiliation :: Joint Physics Analysis Center (JPAC)
!     Contact     :: ajackura@indiana.edu
!     Last Update :: 20/06/2016
!     =========================================================================
      subroutine pipi_cross_section (s,sigma)
      implicit none
      real (8)            :: s, sigma, theta, xkmom
      complex (8)         :: amp
      real (8), parameter :: xmpi = 0.1396d0
      xkmom = 0.5d0 * sqrt(s - 4.d0 * xmpi**2)
      theta = 0.d0
      call pipi_amp_charged (s,theta,amp)
      sigma = dimag (amp) / ( 2.d0 * xkmom * dsqrt(s) )
      return
      end subroutine pipi_cross_section


      subroutine pipi_intensities (s,sigma)
      implicit none
      real (8)                        :: s, xkmom
      integer                         :: i
      integer, parameter              :: nset = 6
      real (8), parameter             :: xmpi = 0.1396d0
      real (8), parameter             :: pi = 4.d0 * datan (1.d0)
      character (2), dimension (nset) :: waveset
      complex (8), dimension (nset)   :: amp
      real (8), dimension (nset)      :: clebsch, sigma
      waveset (1) = 'S0'
      waveset (2) = 'S2'
      waveset (3) = 'P1'
      waveset (4) = 'D0'
      waveset (5) = 'D2'
      waveset (6) = 'F1'
                                !0 1.d0 / 3.d0
      clebsch (1) = 1.d0 / 3.d0 !1 1.d0 / 2.d0
      clebsch (2) = 1.d0 / 6.d0 !2 1.d0 / 6.d0
      clebsch (3) = 1.d0 / 2.d0
      clebsch (4) = 1.d0 / 3.d0
      clebsch (5) = 1.d0 / 6.d0
      clebsch (6) = 1.d0 / 2.d0
      xkmom = 0.5d0 * sqrt(s - 4.d0 * xmpi**2)
      do i = 1, nset
         call pipi_pw_amp (s,waveset (i),amp (i))
         sigma (i) = 32.d0 * pi * clebsch (i) * dimag (amp (i)) /
     c               ( 2.d0 * xkmom * dsqrt(s) )
      enddo
      return
      end subroutine pipi_intensities

!     =========================================================================
!     subroutine :: legendreP
!     -------------------------------------------------------------------------
!     Summary :: Code produces Legendre polynomials P_{l} (z) and their first
!                derivative w.r.t 'z' up to l = 9
!
!     Input  :: z     - Argument of Legendre polynomial (z = cos (theta))
!
!     Output :: pol   - Legendre polynomial
!               pol1p - Derivative of Legendre polynomial w.r.t z
!
!     Type declarations - real (8)      :: z
!                       - real (8), dimension (0:9) :: pol, pol1p
!
!     To call :: call legendreP (z,pol,pl1p)
!
!     Author      :: Andrew W. Jackura
!     Affiliation :: Joint Physics Analysis Center (JPAC)
!     Contact     :: ajackura@indiana.edu
!     Last Update :: 20/06/2016
!     =========================================================================
      subroutine legendreP(z,pol,pol1p)
      implicit none
      real (8) :: z
      real (8), dimension (0:9) :: pol, pol1p
      pol(0) = 1.d0
      pol(1) = z
      pol(2) = (3*z*z - 1.d0)/2.d0
      pol(3) = (-3*z + 5*z*z*z)/2.d0
      pol(4) = (3.d0 - 30*z*z + 35*z**4)/8.d0
      pol(5) = (15*z - 70*z*z*z + 63*z**5)/8.d0
      pol(6) = (-5.d0 +105*z*z - 315*z**4 + 231*z**6)/16.d0
      pol(7) = (-35*z + 315*z*z*z - 693*z**5 + 429*z**7 )/16.d0
      pol(8) = (35.d0 - 1260*z*z + 6934*z**4 - 12012*z**6 +
     c     6435*z**8)/128.d0
      pol(9) = (315*z - 4620*z*z*z + 18018*z**5 - 25740*z**7 +
     c     12155*z**9 )/129.d0

      pol1p(0) = 0.d0
      pol1p(1) = 1.d0
      pol1p(2) = 3.d0*z
      pol1p(3) = (-3.d0 + 15*z*z)/2.d0
      pol1p(4) = (-60*z + 140*z*z*z)/8.d0
      pol1p(5) = (15.d0 - 210*z*z + 315*z**4)/8.d0
      pol1p(6) = (210*z - 1260*z*z*z + 1386*z**5)/16.d0
      pol1p(7) = (-35.d0 + 945*z*z - 3465*z**4 + 3003*z**6)/16.d0
      pol1p(8) = (-2520*z + 27720*z*z*z - 72072*z**5 +
     c     51480*z**7)/128.d0
      pol1p(9) = (315.d0 - 13860*z*z + 90090*z**4 - 180180*z**6 +
     c     109395*z**8)/128.d0
      return
      end subroutine legendreP
