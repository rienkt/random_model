! This program is tests for 1-D autocorrelation function 
!
USE nrtype
USE acorr
USE nr, ONLY : realft
IMPLICIT NONE
INTEGER :: n,lagmax
PARAMETER(n=1024, lagmax=200)
INTEGER :: i,j,no2
REAL,DIMENSION(n) :: org_data
! DIRECT METHOD
REAL, DIMENSION(-n:n*2) :: a
REAL, DIMENSION(-n:n) :: b
REAL, DIMENSION(n*2-1) :: bs_full
REAL, DIMENSION(lagmax*2+1) :: bs_lag
! FFT METHOD
REAL, DIMENSION(n*2) :: a2, b2, bsf
COMPLEX, DIMENSION(n) :: cdat

! Data Input
open(100,FILE="temp.dat")
do i=1,n
   read(100,*) org_data(i)
enddo
close(100)
print *, size(org_data)

!Direct Method by acorr1d (subroutine including in acorr.f90)
bs_full=0.d0
print *, "===== Full acorr1d ===="
call acorr1d(org_data, bs_full)
call acorr1d(org_data, bs_lag, lagmax)

!open(100,FILE="temp_corr_s_full.dat")
!write(100,'(F9.5)') (bs_full(i),i=n,2*n-1)
!close(100)

open(100,FILE="calc_acorr.dat")
write(100,'(F9.5)') (bs_lag(i),i=lagmax+1,2*lagmax+1)
close(100)


end
