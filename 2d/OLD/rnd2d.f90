PROGRAM rnd2d
!**************************************************************
! This program will generate bi-modal random field.
! 
! 2007/10/03
!   Revision
!   To converge faster, we will change the way to find the bimodal
!   point.
!   1. prepare the array
!   2. 
!
!**************************************************************
!
USE acflib, ONLY : gauss2df, exp2df,vk2df
USE nrtype
!!USE nr, ONLY : ran3,moment
!!USE ran_state
USE header_stat
USE header2d
USE fft_params_2d
USE subrnd2d
USE mt19937

!=============================================================
! Define variables
!=============================================================
IMPLICIT NONE
REAL(DP), PARAMETER :: eps=0.0001_dp
REAL(DP), PARAMETER :: large=5.0_dp ! large*std = max(gaussian)
REAL(DP), PARAMETER :: lim=4.0_dp
REAL(DP), PARAMETER :: lw=0.15
INTEGER,  PARAMETER :: iter_max=10
! define general temporary variables 
INTEGER :: i,j, ireal,iter,icount,iolength
! Tolerability
REAL(DP) :: tol
! for error analysis
REAL(DP) :: err,err2,err3
! define 
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zb,zg,zge,zbe
! define
REAL(SP) :: tmp
CHARACTER(LEN=20) :: fname
character(LEN=10) ::  time
REAL(SP) :: ihi,tmpin

!================================================================
! PREPARATION
!===============================================================
print *, 'in'
!----------------------------------------------
! * Read parameters and allocate allays
!-----------------------------------------------
call read_params
! Allocate arrays
allocate(ip(0:2+nint(sqrt(float(max(nx/2,nz))))),w(0:max(nx/4,nz/2)+nx/4))
allocate(t(0:8*nz-1))
print *, 'set'
allocate(sg0(nx/2+1,nz))
ip(0)=0; !ip2(0)=0 ! initialization for fft4g.f

print *, 'set2'
!--------------------------------------------
! * Check length of PSD is large enough
!      for 1D exponential type :-P
!--------------------------------------------
tol=1.0_dp/min(ax,az)*tan((1.0_dp-eps)*PI_D/2.0_dp)
write(*,1001) atan(PI_D/dz)/PIO2_D, PI_D/dz
1001 format(f6.4, " of total spectral is contained at Nyquist ", f8.4)
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
!open(unit=100,file="tmp/s0.bin",form="unformatted")
!write(100) ((s0(i,j)/s0(2,2),i=1,nx/2+1),j=1,nz)
!do i=1,nx/2+1
!   do j=1,nz
!      tmp=s0(i,j)/s0(2,2)
!      write(100) tmp
!   end do
!enddo
!close(100)

!--------------------------------------
! 1. Read input PDF for Gaussian
!    - set initial PSD for Gaussian
!---------------------------------------
if ( ratio1 < 1.0_dp )then
   icount=1
   inquire(iolength=iolength) float(1)
!   print *, iolength
   open(100,file="sg.bin",form='unformatted',access='direct',recl=iolength)
      do i=1,nx/2+1
         do j=1,nz
         ! read(100) ((sg0(i,j),j=1,nz),i=1,nx/2+1)
         read(100,rec=icount) tmpin
         sg0(i,j)=dble(tmpin)
         icount=icount+1
      enddo
   enddo
   close(100)
!   sg0=s0(1:n/2+1,1:nz)
else
   sg0=s0(1:n/2+1,1:nz);
end if
sg0(1,1)=0.0_dp
deallocate(s0);

! Initialize random generator
!call ran_seed(iseed)      ! N.R. 
call init_genrand(iseed)   ! MT

! << for OMAKE (error analysis)
!allocate(zge(nx,nz),zbe(nx,nz),sge(nx/2+1,nz),sbe(nx/2+1,nz))
! >>
do ireal=1,nreal
   print *, ""
   print *, "[REALIZATION : ", ireal, " ]"
   allocate(theta(nx/2+1,nz),zg(nx,nz),zb(nx,nz))
   !allocate(sb(nx/2+1,nz),sg(nx/2+1,nz))  
   ! * Phase angle
   theta=0.0
   open(1000,file='theta.dat')
   do i=1,nx/2+1
      do j=1,nz
         theta(i,j)=grnd()
         write(1000,*) i,j,theta(i,j)
      end do
   end do
   close(1000)
   theta=theta*TWOPI_D
   
   ! * Pseudo-Standard Gaussian field
!   print *, "std_gauss_field",size(zg,1),size(zg,2),size(sg0,1),size(sg0,2)
!   print *, size(theta,1), size(theta,2),dkx,dkz,lim,lw
   print *, "- generating Gaussian random field"
   call std_gauss_field(zg,sg0,theta,lim,lw) 
   !print *, "-cal_psd"
   !call calc_psd(zg,sg,dz)
   
   ! * Bi-modal field
   if (ratio1 < 1) then
      print *, "- converting to bimodal field"
      call cnv_bimodal(zb,zg,lim)
      !print *, "calc_psd"
      !call calc_psd(zb-sum(zb)/n,sb,dz)
   else
      zb=zg*std1+mu1
!      sb=sg
   end if
   ! * Write data
   call DATE_AND_TIME (TIME=time)
   print *, 'time = ', time

   print *, "- writing to file"
   write(fname,'("zz",i5.5,".bin")') ireal
   call write_direct(fname,zb(1:nx0,1:nz0),nx0,nz0)
   
!         write(100) ((zb(i,j),j=1,nz0),i=1,nx0)
!   open(100,file=fname,form='unformatted')
!    do i=1,nx0
!      do j=1,nz0
!         print *, i,j
!         tmp=zb(i,j)
!         write(100) tmp
!      end do
!   end do
!   close(100)


!   call DATE_AND_TIME (TIME=time)
!   print *, 'time = ', time

   !write(fname,'("tmp/zg",i5.5,".dat")') ireal
   !open(100,file=fname)
   !write(100,'(D14.8)') ((zg(i,j),j=1,nz0),i=1,nx0)
   !close(100)

   !call DATE_AND_TIME (TIME=time)
   !print *, 'time = ', time

   !write(fname,'("tmp/sg",i5.5,".dat")') ireal
   !open(100,file=fname)
   !write(100,'(D14.8)') ((sg(i,j)/sg(2,2),j=1,nz),i=1,nx/2+1)
   !close(100)

!   if (ratio1 < 1) then
!      write(fname,'("tmp/zb",i5.5,".bin")') ireal
!      open(unit=100,file=fname,form='unformatted')
!      write(100) ((zb(i,j),j=1,nz0),i=1,nx0)
!      close(100)

!      write(fname,'("tmp/zb",i5.5,".dat")') ireal
!      open(100,file=fname)
!      write(100,'(D14.8)') ((zb(i,j),j=1,nz),i=1,nx)
!      close(100)
      !write(fname,'("tmp/sb",i5.5,".dat")') ireal
      !open(100,file=fname)
      !write(100,'(D14.8)') ((sb(i,j)/sb(2,2),j=1,nz),i=1,nx/2+1)
      !close(100)
!   else
!      write(fname,'("tmp/zg",i5.5,".bin")') ireal
!      open(unit=100,file=fname,form='unformatted')
!      write(100) ((zg(i,j),j=1,nz0),i=1,nx0)
!      close(100)

!   end if


   ! << for OMAKE (error analysis) 
   ! * Calculate emsemble 
!   zge=zge+zg
!   sge=sge+sg
!   if (ratio1 < 1) then
!      zbe=zbe+zb
!      sbe=sbe+sb
!   end if
   ! >>
   
   deallocate(theta,zg,zb)
enddo

!------------------------------------------------
! >> OMAKE
!-------------------------------------------------
!write(fname,'("tmp/zge.dat")') 
!open(100,file=fname)
!write(100,'(D14.8)') ((zge(i,j)/real(nreal),j=1,nz),i=1,nx)
!close(100)
!write(fname,'("tmp/sge.dat")') 
!open(100,file=fname)
!write(100,'(D14.8)') ((sge(i,j)/sge(2,2), j=1,nz),i=1,nx/2+1)
!close(100)
!if ( ratio1 < 1) then
!   write(fname,'("tmp/zbe.dat")') 
!   open(100,file=fname)
!   write(100,'(D14.8)') ((zbe(i,j)/real(nreal),j=1,nz),i=1,nx)
!   close(100)
!   write(fname,'("tmp/sbe.dat")') 
!   open(100,file=fname)
!   write(100,'(D14.8)') ((sbe(i,j)/sbe(2,2),j=1,nz),i=1,nx/2+1)
!   close(100)
!end if

! * Calculate error with three different ways
!   - err   = sum(sbe-s0)
!   - err2  = sum(sbe-s0)**2
!   - err3  = sum(log(sbe)-log(s0))**2
! 
!if ( ratio1 < 1 ) then
!   sbe=sbe/sbe(2,2) ! Normalise PSD
!   sge=sge/sge(2,2)
!   s0=s0/s0(2,2)
!   s0(1,1)=0.0_dp
!   sge(1,1)=0.0_dp
!   sbe(1,1)=0.0_dp
!   
!   err=0.0_dp     ! Initialize errors
!   err2=0.0_dp
!   err3=0.0_dp
!   do i=1,nx/2+1
!      do j=1,nz
!         err=err+abs(sbe(i,j)-s0(i,j))
!         err2=err2+(sbe(i,j)-s0(i,j))**2
!         if ( i==1 .and. j==1 ) then
!            err3=err3+(log(sbe(i,j))-log(s0(i,j)))**2
!         end if
!      end do
!   enddo
!   write(*,'(i2.2," : error ", d12.5, d12.5,d12.5)') iter, err/sum(sge), err2,err3
!end if
!deallocate(zge,zbe,sge,sbe)
 
!-----------------------------------------------------
! << OMAKE
!-----------------------------------------------------  

END PROGRAM rnd2d

SUBROUTINE write_direct(fname,data,nx,nz)
USE nrtype

character(LEN=20) :: fname
REAL(DP) :: data(nx,nz)
REAL(SP) :: tmp
INTEGER :: nx,nz,count,irecl
tmp=data(1,1)
inquire(iolength=irecl) tmp
!print *, irecl

count=1
!print *, "writing",nx,nz
open(100,file=fname,access="direct",form="unformatted",recl=irecl)
do i=1,nx
   do j=1,nz
      tmp=data(i,j)
      write(100,rec=count) tmp
      count=count+1
   enddo
end do
close(100)
END SUBROUTINE write_direct
