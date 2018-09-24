MODULE fft_call
  USE nrtype

  INTERFACE
     SUBROUTINE rdft(n,isgn,a,ip,w)
       USE nrtype
       INTEGER :: n,isgn
       INTEGER, DIMENSION(0:*) :: ip
       REAL(DP), DIMENSION(0:n-1) :: a
       REAL(DP), DIMENSION(0:*) :: w
     END SUBROUTINE rdft
  END INTERFACE
  INTERFACE
     SUBROUTINE cdft(n,isgn,a,ip,w)
       USE nrtype
       INTEGER :: n,isgn
       INTEGER, DIMENSION(0:*) :: ip
       REAL(DP), DIMENSION(0:n-1) :: a
       REAL(DP), DIMENSION(0:*) :: w
     END SUBROUTINE cdft
  END INTERFACE

END MODULE fft_call
