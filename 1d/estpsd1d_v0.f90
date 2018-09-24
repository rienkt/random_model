PROGRAM estpsd1d
!**************************************************************
! This program will generate PDF for bi-modal random field.
!**************************************************************
!
USE acflib, ONLY : gauss1df, exp1df,vk1df
USE nrtype
!USE nr, ONLY : ran3,moment
!USE ran_state
USE header_stat
USE header1d
USE modfft4gsub
USE subrnd1d
USE mt19937

!=============================================================
! Define variables
!=============================================================
IMPLICIT NONE
REAL(DP), PARAMETER :: eps=0.01_dp
REAL(DP), PARAMETER :: large=5.0_dp ! large*std = max(gaussian)
REAL(DP), PARAMETER :: lim=4.0_dp
REAL(DP), PARAMETER :: lw=0.15
INTEGER,  PARAMETER :: iter_max=20
! define general temporary variables 
INTEGER :: i,j,m, ireal,iter
! Tolerability
REAL(DP) :: tol
! for error analysis
REAL(DP) :: err,err2,err3
! define 
REAL(DP), DIMENSION(:), ALLOCATABLE :: zb,zg,zge,zbe
! define
CHARACTER(LEN=20) :: fname


!================================================================
! PREPARATION
!===============================================================

!----------------------------------------------
! * Read parameters and allocate allays
!-----------------------------------------------
call read_params
! Allocate arrays
allocate(ip(0:2+int(sqrt(real(n/2)))),w(0:n/2-1))
allocate(ip2(0:2+int(sqrt(real(n)))),w2(0:n-1))
allocate(sg0(n/2+1))
ip(0)=0; ip2(0)=0 ! initialization for fft4g.f

!--------------------------------------------
! * Check length of PSD is large enough
!      for exponential type :-P
!--------------------------------------------
tol=1.0_dp/a*tan((1.0_dp-eps)*PI_D/2.0_dp)
write(*,1001) atan(PI_D/dz)/PIO2_D, PI_D/dz
1001 format(f6.4, " of total spectral is contained at Nyquist ", f8.4)
write(*,1002) tol
1002 format( "recommended upper cut frequency ", f8.4)

!------------------------------------
! 0. set target PSD for Bi-modal
!------------------------------------
allocate(s0(n/2+1))
if (type ==1) then
   call gauss1df(s0,a,dk)
elseif (type==2) then
   call exp1df(s0,a,dk)
elseif (type==3) then
   call vk1df(s0,a,dk,nu)
else 
   stop
end if

! Write target PSD
open(100,file="tmp/s0.dat")
write(100,'(F14.8)') (s0(i)/s0(2),i=1,n/2+1)
close(100)

!--------------------------------------
! 1. Read input PDF for Gaussian
!    - set initial PSD for Gaussian
!---------------------------------------


! Initialize random generator
!call ran_seed(iseed)      ! N.R. 
call init_genrand(iseed)   ! MT

sg0=s0

do iter=1,iter_max
! For error analysis
   allocate(zge(n),zbe(n),sge(n/2+1),sbe(n/2+1))
   ! >>
   do ireal=1,nreal
      print *, ireal
      allocate(theta(n/2+1),zg(n),zb(n),sb(n/2+1),sg(n/2+1))  
      ! * Phase angle
      theta=0.0
      do i=1,n/2+1
         theta(i)=grnd()
      end do
      theta=theta*TWOPI_D
      print*, theta(1), theta(30)
      
      ! * Pseudo-Standard Gaussian field
      call std_gauss_field(zg,sg0,theta(1:n/2),dk,lim,lw) 
      call calc_psd(zg,sg,dk,dz)
      
      ! * Bi-modal field
      call cnv_bimodal(zb,zg,lim)
      call calc_psd(zb-sum(zb)/n,sb,dk,dz)
      
      ! * Calculate emsemble 
      zge=zge+zg
      zbe=zbe+zb
      sge=sge+sg
      sbe=sbe+sb 

      if (ireal .eq. 1) then
         write(fname,'("tmp/zg",i2.2,".dat")') iter
         open(100,file=fname)
         write(100,'(D14.8)') (zg(i), i=1,n)
         close(100)
         write(fname,'("tmp/zb",i2.2,".dat")') iter
         open(100,file=fname)
         write(100,'(D14.8)') (zb(i), i=1,n)
         close(100)
      end if

      
      deallocate(theta,zg,zb,sb,sg)
   enddo
   
   write(fname,'("tmp/sg",i2.2,".dat")') iter
   open(100,file=fname)
   write(100,'(D14.8)') (sge(i)/sge(2), i=1,n/2+1)
   close(100)
   write(fname,'("tmp/zge",i2.2,".dat")') iter
   open(100,file=fname)
   write(100,'(D14.8)') (zge(i)/real(nreal),i=1,n)
   close(100)
   write(fname,'("tmp/zbe",i2.2,".dat")') iter
   open(100,file=fname)
   write(100,'(D14.8)') (zbe(i)/real(nreal),i=1,n)
   close(100)
   write(fname,'("tmp/sge",i2.2,".dat")') iter
   open(100,file=fname)
   write(100,'(D14.8)') (sge(i)/sge(2), i=1,n/2+1)
   close(100)
   write(fname,'("tmp/sbe",i2.2,".dat")') iter
   open(100,file=fname)
   write(100,'(D14.8)') (sbe(i)/sbe(2),i=1,n/2+1)
   close(100)

   ! * Calculate error with three different ways
   !   - err   = sum(sbe-s0)
   !   - err2  = sum(sbe-s0)**2
   !   - err3  = sum(log(sbe)-log(s0))**2
   ! 
   sbe=sbe/(sbe(2)+sbe(3))*2.0_dp ! Normalise PSD
   sge=sge/sge(2)
   s0=s0/s0(2)     
   err=0.0_dp     ! Initialize errors
   err2=0.0_dp
   err3=0.0_dp
   do i=2,n/2+1-lw*n
      err=err+abs(sbe(i)-s0(i))
      err2=err2+(sbe(i)-s0(i))**2
     if ((sbe(i).gt.1.0e-5).and. (s0(i) .gt. 1.0e-5)) then 
     err3=err3+(log(sbe(i))-log(s0(i)))**2
     end if
   end do
   write(*,'(i2.2," : error ", d12.5, d12.5,d12.5)') iter, err, err2,err3
   
   
   do i=2,n/2
      sg0(i)=sg0(i)/sbe(i)*s0(i)
   enddo
   

   deallocate(zge,zbe,sge,sbe)

enddo

END PROGRAM estpsd1d
