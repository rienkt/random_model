PROGRAM rnd2d_v0
!OCL INDEPENDENT(out)
!OCL INDEPENDENT(cpd_std_gauss,cpd_bi_gauss,erf)
!**************************************************************
! This program will generate bi-modal random field.
!**************************************************************
!
USE acflib, ONLY : gauss2df, exp2df,vk2df
USE nrtype
!USE nr, ONLY : ran3,moment
!USE ran_state
USE header_stat
USE header2d
USE fft_params_2d
USE subrnd2d
USE mt19937

!=============================================================
! Define variables
!=============================================================
IMPLICIT NONE
REAL(DP), PARAMETER :: eps=0.01_dp
REAL(DP), PARAMETER :: large=5.0_dp ! large*std = max(gaussian)
REAL(DP), PARAMETER :: lim=4.0_dp
REAL(DP), PARAMETER :: lw=0.1
INTEGER,  PARAMETER :: iter_max=10
! define general temporary variables 
INTEGER :: i,j, ireal,iter,nw
! Tolerability
REAL(DP) :: tol
! for error analysis
REAL(DP) :: err,err2,err3
! define 
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zb,zg,zge,zbe
! define
CHARACTER(LEN=20) :: fname
character(LEN=10) ::  time

!================================================================
! PREPARATION
!===============================================================

!----------------------------------------------
! * Read parameters and allocate allays
!-----------------------------------------------
call read_params
! Allocate arrays
allocate(ip(0:2+nint(sqrt(float(max(nx/2,nz))))),w(0:max(nx/4,nz/2)+nx/4))
allocate(t(0:8*nz-1),aa(0:nx+1,0:nz-1))
allocate(sg0(nx/2+1,nz))
ip(0)=0;! initialization for fft4g.f

!--------------------------------------------
! * Check length of PSD is large enough
!      for 1D exponential type :-P
!--------------------------------------------
tol=1.0_dp/min(ax,az)*tan((1.0_dp-eps)*PI_D/2.0_dp)
write(*,1001) atan(PI_D/dz)/PIO2_D, PI_D/dz
1001 format(f7.4, " of total spectral is contained at Nyquist ", f8.4)
write(*,1002) tol
1002 format( "recommended upper cut frequency ", f8.4)

!------------------------------------
! 0. set target PSD for Bi-modal
!------------------------------------
allocate(s0(nx,nz))
if (type ==1) then
   call gauss2df(s0,ax,az,dkx,dkz,0.0_dp)
elseif (type==2) then
   call exp2df(s0,ax,az,dkx,dkz,0.0_dp)
elseif (type==3) then
   call vk2df(s0,ax,az,nu,dkx,dkz,0.0_dp)
else 
   stop
end if

! Write target PSD
s0(1,1)=0
s0=s0/s0(1,2)
print *, s0(1,2)
s0(1,1)=0.0_dp
open(unit=100,file="tmp/s0.bin",form='unformatted')
write(100) ((s0(i,j)/s0(1,2),i=1,nx/2+1),j=1,nz)
close(100)

!--------------------------------------
! 1. Read input PDF for Gaussian
!    - set initial PSD for Gaussian
!---------------------------------------

sg0=s0(1:nx/2+1,1:nz);
sg0(1,1)=0.0_dp
!deallocate(s0)

! Initialize random generator
!call ran_seed(iseed)      ! N.R. 
call init_genrand(iseed)   ! MT

print *, iter_max
 
do iter=1,iter_max
   print *, "----------------------------------"
   call DATE_AND_TIME (TIME=time)
   print *, "* Iteration : ", iter,time 
!   allocate(zge(nx,nz),zbe(nx,nz),sge(nx/2+1,nz),sbe(nx/2+1,nz))
   allocate(sbe(nx/2+1,nz))
   sbe=0.0d0
   do ireal=1,nreal
      allocate(theta(nx/2+1,nz))
      allocate(zg(nx,nz))
      allocate(zb(nx,nz))
!      if (mod(ireal,10) == 0) then
      call DATE_AND_TIME (TIME=time)
         print *, "  realization ", ireal,time
!      end if
         zb=0.0d0
         zg=0.0d0
      ! * Phase angle
      print *, 'calculate phase'
      theta=0.0d0
      do i=1,nx/2+1
         do j=1,nz
            !zg(i,j)=grnd()
            theta(i,j)=grnd()
         end do
      end do
      theta=theta*TWOPI_D
      
      ! * Pseudo-Standard Gaussian field
      print *, "std_gauss_field aho"
      call std_gauss_field(zg,sg0,theta,lim,lw) 
!      call std_gauss_field(zg,sg0,theta,lim,lw) 
!      deallocate(theta)
!      print *, "calc_psd"
!      call calc_psd(zg,sg,dz)
      ! * Bi-modal field

      print *, "* cnv_bimodal"
      call cnv_bimodal(zb,zg,lim)


!      deallocate(zg)
!      allocate(sb(nx/2+1,nz))
!      sb=sg0
!      print *, "* calc_psd"
      tol=sum(zb)/float(n)
      zb=zb-tol
      call calc_psd(zb,zg(1:nx/2+1,1:nz),dz)
!      zg=zg/zg(1,2)
      zg=zg/1000.0_dp;
      ! * Write data
      !   write(fname,'("tmp/zg",i5.5,".dat")') ireal
      !   open(100,file=fname)
      !   write(100,'(D14.8)') ((zg(i,j),j=1,nz),i=1,nx)
      !   close(100)
      !   write(fname,'("tmp/zb",i5.5,".dat")') ireal
      !   open(100,file=fname)
      !   write(100,'(D14.8)') ((zb(i,j),j=1,nz),i=1,nx)
      !   close(100)
      !   write(fname,'("tmp/sb",i5.5,".dat")') ireal
      !   open(100,file=fname)
      !   write(100,'(D14.8)') ((sb(i,j)/sb(2,2),j=1,nz),i=1,nx/2+1)
      !   close(100)
      !   write(fname,'("tmp/sg",i5.5,".dat")') ireal
      !   open(100,file=fname)
      !   write(100,'(D14.8)') ((sg(i,j)/sg(2,2),j=1,nz),i=1,nx/2+1)
      !   close(100)
      
      ! * Calculate emsemble 
!      zge=zge+zg
!      zbe=zbe+zb
!      sge=sge+sg
!      sbe=sbe+sb
    !  print *, zg(20,20), sbe(20,20)
      sbe=sbe+zg(1:nx/2+1,1:nz)
!      if (mod(ireal,10)==0) then
!         write(fname,'("tmp/sbe",i3.3,".bin")') ireal
!         open(unit=100,file=fname,form='unformatted')
!         write(100) ((sbe(i,j)/sbe(1,2),j=1,nz),i=1,nx/2+1)
!         close(100)
!      end 
 if
  
     ! >>
 
      deallocate(zb,theta,zg)
!      deallocate(zb,sb)
   enddo

   sbe(1,1)=0.25_dp*(sbe(2,2)+sbe(1,2)+sbe(2,1)+sbe(2,nz))
   
 
!  writing data is time-consuming  
!   print *, 'write_data'
!   write(fname,'("tmp/zge",i2.2,".dat")') iter
!   open(100,file=fname)
!   write(100,'(D14.8)') ((zge(i,j)/real(nreal),j=1,nz),i=1,nx)
!   close(100)
!   write(fname,'("tmp/zbe",i2.2,".bin")') iter
!   open(unit=100,file=fname,form='unformatted')
!   write(100) ((zbe(i,j)/real(nreal),j=1,nz),i=1,nx)
!   close(100)

   write(fname,'("tmp/sg",i2.2,".bin")') iter
   open(unit=100,file=fname,form='unformatted')
   write(100) ((sg0(i,j), j=1,nz),i=1,nx/2+1)
   close(100)
   write(fname,'("tmp/sbe",i2.2,".bin")') iter
   open(unit=100,file=fname,form='unformatted')
   write(100) ((sbe(i,j)/sbe(1,1),j=1,nz),i=1,nx/2+1)
   close(100)
   
   ! * Calculate error with three different ways
   !   - err   = sum(sbe-s0)
   !   - err2  = sum(sbe-s0)**2
   !   - err3  = sum(log(sbe)-log(s0))**2
   ! 
   !print *, 'prepare for error analysis'
   sbe=sbe/sbe(1,1) ! Normalise PSD
 !   sge=sge/sge(2,2)
 !  sge(1,1)=0.0_dp
   sbe(1,1)=0.0_dp

   
   !print *, 'calculate error'
   err=0.0_dp     ! Initialize errors
   err2=0.0_dp
   err3=0.0_dp    
   nw=int(min(nx,nz)/2*lw)
   !print *, nw
   do i=1,nx/2+1-nw
      do j=1,nz-nw
         err=err+abs(sbe(i,j)-s0(i,j))
         err2=err2+(sbe(i,j)-s0(i,j))**2
         if ( i>1 .or. j>1 ) then
            err3=err3+(log(sbe(i,j))-log(s0(i,j)))**2
         end if
      end do
   enddo
   write(*,'(i2.2," : error ", d12.5, d12.5,d12.5)') iter, err/sum(s0), err2,err3
!   close(100)


   do i=1,nx/2+1-nw
      do j=1,nz-nw
         if ( i>1 .or. j>1) then
            sg0(i,j)=sg0(i,j)/sbe(i,j)*s0(i,j)
         endif
      end do
   end do
   !print *, 'update data'
   deallocate(sbe)

!   deallocate(zge,zbe,sge,sbe)
enddo
!-----------------------------------------------------
! << OMAKE
!-----------------------------------------------------  

EN D PROGRAM rnd2d_v0
