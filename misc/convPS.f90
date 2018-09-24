PROGRAM convPS
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
REAL(SP),dimension(:), allocatable  :: vporg
REAL(SP) :: vpout, vsout, q,qinv
REAL(SP) :: poisson, vsvp, vpvs,density,vsfixed
INTEGER  :: type,option,index,nx,nz
INTEGER  :: option_pwave,option_density, option_swave, option_pfile
INTEGER  :: i, icount, iolen_in, iolen_out
CHARACTER(LEN=80) :: fname_vpin, fname_vpout, fname_vsout, fname_rhoout
CHARACTER(LEN=80) :: fname_qout, fname_qinvout
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
allocate(vporg(nx*nz))
open(100,file=fname_vpin,form='unformatted',access='direct',recl=iolen_in)
icount=1
do i=1,nx*nz
   read(100,rec=icount) vporg(i)
   icount=icount+1
enddo
close(100)

!-------------------------------------------------------------
! Action
!------------------------------------------------------------

type=1
do while ( type < 10 )
   write(*,*) "Choose action"
   write(*,*) "1. Create P-wave file"
   write(*,*) "2. Create S-wave file"
   write(*,*) "3. Create density file"
   write(*,*) "4. Create Q and 1/Q files"
   write(*,*) "99. Exit"
   read(*,*) type
   select case(type)
   case(1)
      if (option_pfile==1) then 
         write(fname_vpout,'("vp",i5.5,".bin")') index
      else
         write(*,*) "Input output P-wave file name"
         read(*,*) fname_vpout
      end if
      open(101,file=fname_vpout,form='unformatted',access='direct',recl=iolen_out)
      icount=1
      do i=1,nx*nz
         vpout=vporg(i)
         write(101,rec=icount) vpout
      end do
      close(101)
   case(2)
      if (option_pfile==1) then 
         write(fname_vsout,'("vs",i5.5,".bin")') index
      else
         write(*,*) "Input output S-wave file name"
         read(*,*) fname_vsout
      end if
      write(*,*) "Choose conversion method from P-wave velocity to S-wave velocity"
      write(*,*) "1. fixing Vp/Vs ratio"
      write(*,*) "2. fixing Poisson's ratio"
      write(*,*) "3. fixing values"
      read(*,*) option
      if (option == 1) then
         write(*,*) "Input Vp/Vs ratio"
         read(*,*) vpvs
         vsvp=1./vpvs
      else if (option == 2) then
         write(*,*) "Input Poisson's ratio"
         read(*,*) poisson
         vsvp=sqrt((1.-2.*poisson)*0.5/(1.-poisson))
      else if (option == 3) then
         write(*,*) "Input fixed value"
         read(*,*) vsfixed
      end if
      open(101,file=fname_vsout,form='unformatted',access='direct',recl=iolen_out)
      icount=1
      if ( option == 3 ) then
         do i=1,nx*nz
            vpout=vsfixed
            write(101,rec=icount) vsout
            icount=icount+1
         end do
      else
         do i=1,nx*nz
            vpout=vporg(i)
            if (abs(vpout-1500.0)<1e-4) then
               vsout=0
            else
               vsout=vpout*vsvp
            endif
            write(101,rec=icount) vsout
            icount=icount+1
         end do
      endif
      close(101)
   case(3)
      if (option_pfile==1) then 
         write(fname_rhoout,'("rho",i5.5,".bin")') index
      else
         write(*,*) "Input output Density file name"
         read(*,*) fname_rhoout
      end if

      write(*,*) "Choose method to decide density "
      write(*,*) "1. homogeneous values"
      write(*,*) "2. Gardner's formula"
      read(*,*) option
      if ( option == 1) then
         write(*,*) "Input density values"
         read(*,*) density
      endif
      open(101,file=fname_rhoout,form='unformatted',access='direct',recl=iolen_out)
      icount=1
      do i=1,nx*nz
         vpout=vporg(i)
         if (option ==1) then
            write(101,rec=icount) density
         else
            if (abs(vpout-1500.0) < 1e-4) then
               density=1050
            else
               density=1000*0.31*vpout**0.25
            endif
            write(101,rec=icount) density
         end if
         icount=icount+1
      end do
      close(101)
   case(4)
      if (option_pfile==1) then 
         write(fname_qout,'("q",i5.5,".bin")') index
         write(fname_qinvout,'("qinv",i5.5,".bin")') index
      else
         write(*,*) "Input output Q file name"
         read(*,*) fname_qout
         write(*,*) "Input output 1/Q file name"
         read(*,*) fname_qinvout
      end if

!      write(*,*) "Choose method to decide density "
!      write(*,*) "1. homogeneous values"
!      write(*,*) "2. Gardner's formula"
!      read(*,*) option
!      if ( option == 1) then
         write(*,*) "Input Q value"
         read(*,*) q
         qinv=1./q
!      endif
      open(101,file=fname_qout,form='unformatted',access='direct',recl=iolen_out)
      open(102,file=fname_qinvout,form='unformatted',access='direct',recl=iolen_out)
      icount=1
      do i=1,nx*nz
         
         if (abs(vporg(i)-1500.0) < 1e-4) then
!            print *, "here"
            write(101,rec=icount) 1e5
            write(102,rec=icount) 1e-4
         else
!            print *, 
            write(101,rec=icount) q
            write(102,rec=icount) qinv
         endif
         icount=icount+1
      end do
      close(101)
      close(102)



   end select
end do

!------------------------------------------------------------------
! possible improvement
!    - total number of record can be determined while reading data
!       -> no need to input nx,nz
!    - adopting input file for single precision
!-------------------------------------------------------------------

END PROGRAM convPS

