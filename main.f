      program run_pipi
      implicit none
      integer                :: i, j, nstep
      real (8)               :: s, delta, eta, phi
      character (len=2)      :: wav
      complex (8), parameter :: xi = (0.d0,1.d0)
      real (8)               :: sf, step, sigma
      real (8), dimension(6) :: intensity
      real (8), parameter    :: pi = 4.d0*datan(1.d0)
      real (8), parameter    :: xmpi = 0.1396d0, xmK = 0.4936d0
      complex (8)            :: amp

      write (0,*) 'Wave?'
      read *, wav
      s = (2.d0 * xmpi)**2
      sf = 1.42d0**2
      step = 1.d-4
      nstep = nint ((sf-s)/step)

      do i = 1, nstep
!         call pipi_intensities (s,intensity)
         call pipiscat (s,wav,delta,eta)
!         write (*,*) dsqrt (s), delta, eta
         call pipi_pw_amp (s,wav,amp)
!         write (*,*) dsqrt (s), dble (amp), dimag (amp)
         call pipi_cross_section (s,sigma)
!         write (*,*) dsqrt(s), sigma
         write (*,*) dsqrt (s), delta, eta, dble(amp), dimag(amp), sigma
         write (111,*) dsqrt (s), (intensity (j), j=1,6)
         s = s + step
      enddo
      stop
      end program run_pipi
