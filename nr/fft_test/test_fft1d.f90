PROGRAM test_fft1d
USE nrtype 
USE nr, ONLY : realft
IMPLICIT NONE
! 5T
INTEGER :: n
REAL :: dt,f0
PARAMETER(n=1024,dt=1./20,f0=1.0)
INTEGER :: i
REAL :: data_t(n),data_org_t(n)
COMPLEX(SPC),DIMENSION(1:n/2) :: data_f
do i=1,n
   data_t(i)=sin(2.*PI*f0*(i-1)*dt)
!   write(*,*) i,PI,f0,dt,data_t(i)
enddo
data_org_t=data_t
call realft(data_t,1,data_f)
!call realft_sp(data_t,1)
open(100,file="time.dat")
write(100,'(f9.4)') (data_t(i),i=1,n)
close(100)
open(100,file="freq_cabs.dat")
write(100,'(f9.4)') (cabs(data_f(i)),i=1,n/2)
close(100)
call realft(data_t,-1,data_f)
data_t=data_t*2./n
open(100,file="inverse_time.dat")
write(100,'(f9.4)') (data_t(i),i=1,n)
close(100)
open(100,file="inverse_time_dif.dat")
write(100,'(f9.4)') (data_t(i)-data_org_t(i),i=1,n)
close(100)

end PROGRAM test_fft1d
