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
!   print *, org_data(i)
enddo
close(100)
print *, size(org_data)

!Direct Method
a=0.d0
a(1:n)=org_data
!do i=-0,n-1
!   do j=1,n
!      b(i+1)=b(i+1)+a(j)*a(j+i)
!   enddo
!enddo
do i=-n+1,n-1
   do j=1,n
      b(i+1)=b(i+1)+a(j)*a(j+i)
   enddo
enddo

open(100,FILE="temp_corr.dat")
!write(100,'(F9.5)') (b(i),i=1,n)
write(100,'(F9.5)') (b(i),i=-n,n)
close(100)

!Direct Method by acorr1d (subroutine including in acorr.f90)
print *, size(org_data), size(bs_full)
bs_full=0.d0
print *, size(org_data), size(bs_full),org_data(1)
call test_scaler(org_data(1))
call test(org_data)
print *, "===== Full acorr1d ===="
call acorr1d(org_data, bs_full)
call acorr1d(org_data, bs_lag, lagmax)

open(100,FILE="temp_corr_s_full.dat")
write(100,'(F9.5)') (bs_full(i),i=n,2*n-1)
close(100)

open(100,FILE="temp_s_lag.dat")
write(100,'(F9.5)') (bs_lag(i),i=lagmax+1,2*lagmax+1)
close(100)


!FFT Method
a2=0.d0
a2(1:n)=org_data
no2=n
call realft(a2,1,cdat)
cdat(1)=cmplx(real(cdat(1))*real(cdat(1))/no2,aimag(cdat(1))**2/no2,kind=spc)
cdat(2:)=cdat(2:)*conjg(cdat(2:))/no2
call realft(b2,-1,cdat)

open(100,FILE="temp_corr_fft.dat")
write(100,'(F9.5)') (b2(i),i=1,n)
close(100)

! FFT METHOD byy acorr1df
call acorr1df(org_data,bsf)
open(100,FILE="temp_corr_fft_s.dat")
write(100,'(F9.5)') (bsf(i),i=1,n)
close(100)


end
