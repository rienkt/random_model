PROGRAM test_fft2d
USE nrtype 
USE nr, ONLY : rlft2
IMPLICIT NONE
! 5T
INTEGER :: n1,n2
REAL :: dx,dz,f1,f2
PARAMETER(n1=128,n2=128,dx=1./20,dz=1./20, f1=1.0, f2=2.0)
INTEGER :: i,j
REAL, DIMENSION(n1,n2):: data_xz,data_org_xz
COMPLEX(SPC),DIMENSION(1:n1/2,1:n2) :: data_kk
COMPLEX(SPC),DIMENSION(1:n1) :: data_nyq
do i=1,n1
   do j=1,n2
      data_xz(i,j)=sin(2.*PI*f1*(i-1)*dx)+sin(2.*PI*f2*(j-1)*dz)
!   write(*,*) i,PI,f0,dt,data_t(i)
   enddo
enddo
data_org_xz=data_xz
call rlft2(data_xz,data_kk,data_nyq,1)
!call realft_sp(data_t,1)
open(100,file="time.dat")
do j=1,n2
   write(100,'(f9.4)') (data_xz(i,j),i=1,n1)
enddo
close(100)
open(100,file="freq_cabs.dat")
do j=1,n2
   write(100,'(f9.4)') (cabs(data_kk(i,j)),i=1,n1/2)
enddo
close(100)
call rlft2(data_xz,data_kk,data_nyq,-1)
data_xz=data_xz*2./n1/n2
open(100,file="inverse_time.dat")
do j=1,n2
   write(100,'(f9.4)') (data_xz(i,j),i=1,n1)
enddo
close(100)
open(100,file="inverse_time_dif.dat")
do j=1,n2
   write(100,'(f9.4)') (data_xz(i,j)-data_org_xz(i,j),i=1,n1)
enddo
close(100)

end PROGRAM test_fft2d
