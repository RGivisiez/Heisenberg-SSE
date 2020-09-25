Module mean_and_error

  Implicit None

  Contains

  Subroutine basic(measurements, mean, mean_err)

    ! Given n measurements, returns the mean and the error.
    !
    ! Equations used:
      ! mean = \sum_{i=0}^N measurements(i)
      ! mean_err = \sqrt{\frac{1}{n(n-1)} \sum_{i=0}^N (measurements(i) - mean)^2}
      !
      ! Parameters:
      !             measurements:   Real*8, Dimension(:)
      !                             Measurements.
      ! Returns:
      !             mean:           Real*8
      !                             Mean of the measurements.
      !
      !             mean_err:       Real*8
      !                             Mean error.

    Implicit None

    Real*8, Dimension (:), intent(in) :: measurements
    Real*8, intent(out) :: mean, mean_err
    Integer*4 N, i

    N = size(measurements)

    mean = sum(measurements) / N

    mean_err = 0.0d0

    Do i = 1, N

      mean_err = mean_err + (measurements(i) - mean)**2

    end do

    mean_err = mean_err / (N * (N - 1))

  End Subroutine basic

  Subroutine bootstrap(measurements, mean, mean_err, n_resample)

    ! Given n measurements, returns the mean and the error using the
    ! bootstrap method. This method is a resampling method.
    !
      !
      ! Parameters:
      !             measurements:   Real*8, Dimension(:)
      !                             Measurements.
      !             
      !             n_resample:     Integer*4
      !                             Number of times a resample is made.
      !
      ! Returns:
      !             mean:           Real*8
      !                             Mean of the measurements.
      !
      !             mean_err:       Real*8
      !                             Mean error.

    Implicit None

    Real*8, Dimension (:), intent(in) :: measurements
    Integer*4, intent(in) :: n_resample
    Real*8, intent(out) :: mean, mean_err
    Real*8, Allocatable, Dimension (:) :: res_mean, ran
    Integer*4 N, i
    Real*8 mean2

    N = size(measurements)

    Allocate(res_mean(n_resample), ran(N))

    res_mean = 0.0d0
    ran = 0.0d0



    Do i = 1, n_resample

      Call random_number(ran)

      ran = floor(ran * N + 1)

      res_mean(i) = sum(measurements(int(ran))) / N

    end do

    mean = sum(res_mean) / n_resample
    
    mean2 = sum(res_mean**2) / n_resample

    mean_err = sqrt(mean2 - mean**2)

  End Subroutine bootstrap

  Subroutine jackknife(measurements, mean, mean_err)

    ! Given n measurements, returns the mean and the error.
    !
    ! Equations used:
      ! mean = \sum_{i=0}^N measurements(i)
      ! mean_err = \sqrt{\frac{1}{n(n-1)} \sum_{i=0}^N (measurements(i) - mean)^2}
      !
      ! Parameters:
      !             measurements:   Real*8, Dimension(:)
      !                             Measurements.
      ! Returns:
      !             mean:           Real*8
      !                             Mean of the measurements.
      !
      !             mean_err:       Real*8
      !                             Mean error.

    Implicit None

    Real*8, Dimension (:), intent(in) :: measurements
    Real*8, intent(out) :: mean, mean_err
    Real*8, Allocatable, Dimension (:) :: res_mean
    Real*8 summ
    Integer*4 N, i

    N = size(measurements)

    Allocate(res_mean(N))
    
    summ = sum(measurements)

    Do i = 1, N

      res_mean(i) = (summ - measurements(i)) / (N - 1)

    end do

    mean = sum(res_mean) / N

    mean_err = 0.0d0

    Do i = 1, N

      mean_err = mean_err + (res_mean(i) - mean)**2

    end do

    mean_err = mean_err * (dble(N - 1) / dble(N))

  End Subroutine jackknife

End Module mean_and_error

Program Main

  Use mean_and_error
  Implicit None

  Real*8 beta, mean, mean_err
  Character*60 arq
  Integer*4 lx, ly, lz, mcsteps, Nbins, i
  Integer*4, Allocatable, Dimension (:)   :: raw_NH

  Open(20, file='parameters.dat')
  
  read(20, '(A)') arq
  print*, arq

  read(20, *) lx, ly, lz, beta, mcsteps, Nbins
  print*, lx, ly, lz, beta, mcsteps

  write(arq, '("raw_",I0,"x",I0,"x",I0,"_T=",F6.4)') lx, ly, lz, 1.0d0 / beta
  Open(20, file=arq, form='UNFORMATTED')

  Allocate(raw_NH(Nbins * mcsteps))

  Do i = 1, Nbins * mcsteps

    read(20) raw_NH(i)

  end do

  Call basic(dble(raw_NH), mean, mean_err)

  print*, 'Basic:'
  print*, mean, mean_err

  Call bootstrap(dble(raw_NH), mean, mean_err, 50)

  print*, 'Bootstrap:'
  print*, mean, mean_err

  Call jackknife(dble(raw_NH), mean, mean_err)

  print*, 'Jackknife:'
  print*, mean, mean_err

End Program Main