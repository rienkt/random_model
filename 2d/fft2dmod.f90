MODULE fft_call_2d
  USE nrtype

  INTERFACE
     SUBROUTINE rdft2d(n1max,n1,n2,isgn,a,t,ip,w)
       USE nrtype
       INTEGER, INTENT(IN) :: n1max,n1,n2,isgn
       INTEGER, DIMENSION(0:*) :: ip
       REAL(DP), DIMENSION(0:n1max-1,0:n2-1) :: a
       REAL(DP), DIMENSION(0:*) :: w,t
     END SUBROUTINE rdft2d
  END INTERFACE
  INTERFACE
     SUBROUTINE rdft2dsort(n1max,n1,n2,isgn,a)
       USE nrtype
       INTEGER,INTENT(IN) :: n1max,n1,n2,isgn
       REAL(DP), DIMENSION(0:n1max-1,0:n2-1) :: a
     END SUBROUTINE rdft2dsort
  END INTERFACE


END MODULE fft_call_2d

MODULE fft_params_2d
  USE nrtype
  USE header2d
  !
  !----------------------------------------------------------
  ! for fftsg2d.f90
  !  w  : cos, sin table for n-data length
  !  ip : work area for bit reversal for n-data length
  !  w2 : for n*2-data length
  !  ip2: 
  !----------------------------------------------------------
  INTEGER,DIMENSION(:),ALLOCATABLE :: ip
  REAL(DP), DIMENSION(:),ALLOCATABLE :: w,t
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: aa
!  INTEGER,DIMENSION(:),ALLOCATABLE :: ip2
!  REAL(DP), DIMENSION(:),ALLOCATABLE :: w2,t2
  !
END MODULE fft_params_2d
