PROGRAM extractModel
!********************************************************************
! This program will generate single presicion P-wave and S-wave data 
! and density data from P-wave data
! [OPTION for S-wave velocity]
!  1) fixing Vp/Vs ratio
!  2) fixing Poisson's ratio
! [OPTION for density]
!  1) homogeneous
!  2) Gardner's formula (1974)  
!*********************************************************************
USE nrtype
!USE header_stat
!USE header2d

!=====================================================================
! Define variables
!=====================================================================
IMPLICIT NONE
REAL(SP),dimension(:,:), allocatable  :: vporg
REAL(SP) :: vpout
INTEGER  :: type,option,index,nx,nz
INTEGER  :: i, j, icount, iolen_in, iolen_out, option_pfile
INTEGER, dimension(2) :: ix,iz
CHARACTER(LEN=20) :: fname_vpin, fname_vpout, fname_vsout, fname_rhoout

!=====================================================================
! Program start
!=====================================================================
inquire(iolength=iolen_in) 1.0_sp
inquire(iolength=iolen_out) 1.0_sp
!----------------------------------------------
! * Set nx and nz
!-----------------------------------------------
write(*,*) "Input nx and nz"
read(*,*) nx,nz

!-------------------------------------------------
! * Decide Input file and read
!-------------------------------------------------
write(*,*) "Choose type of original P-wave velocity data file name"
write(*,*) "1. generated from random program (zz00000??.bin)"
write(*,*) "2. arbitrary"
read(*,*) option_pfile
if (option_pfile==1) then
   write(*,*) "Which file do you want to convert ?(zz0000??.bin)?"
   read(*,*) index
   write(fname_vpin,'("zz",i5.5,".bin")') index
else
   write(*,*) "Input input P-wave file name"
   read(*,*) fname_vpin
endif
allocate(vporg(nx,nz))
open(100,file=fname_vpin,form='unformatted',access='direct',recl=iolen_in)
icount=1
do i=1,nx
   do j=1,nz
      read(100,rec=icount) vporg(i,j)
      icount=icount+1
   end do
enddo
close(100)

!-------------------------------------------------------------
! Action
!------------------------------------------------------------
write(*,*) "Input coordinates to extract (x1,x2,z1,z2)"
read(*,*) ix(1),ix(2),iz(1),iz(2)
if (option_pfile==1) then 
   write(fname_vpout,'("vp",i5.5,".bin")') index
else
   write(*,*) "Input output P-wave file name"
   read(*,*) fname_vpout
end if

write(*,*) "Now extract", ix(1), ix(2), iz(1), iz(2)
open(101,file=fname_vpout,form='unformatted',access='direct',recl=iolen_out)
icount=1
do i=ix(1),ix(2)
   do j=iz(1),iz(2)
      vpout=vporg(i,j)
      write(101,rec=icount) vpout
      icount=icount+1
   enddo
end do
close(101)


END PROGRAM extractModel

