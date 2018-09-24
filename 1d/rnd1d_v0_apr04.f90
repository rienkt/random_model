PROGRAM rnd1d
! Creating 1-D random field from auto-correlation function.
! Not minimize memory costs.
!
USE acorr
USE nrtype
USE nr, ONLY : ran3,moment
USE ran_state
USE header_stat
USE header1d
USE modfft4gsub
USE subrnd1d
USE mt19937

IMPLICIT NONE
INTEGER :: nbin,nlen
real(DP) :: eps,tiny,large,lim,lw
PARAMETER(eps=0.001_dp,tiny=1.0_dp,large=5.0_dp,lim=4.0_dp,&
     nbin=20, nlen=5, lw=0.15) 
! define general temporary variables 
INTEGER :: i,j,m, index, itol,ireal
REAL(DP) :: tmp,tol
INTEGER, DIMENSION(:), ALLOCATABLE :: itmp_array
REAL(DP),DIMENSION(:), ALLOCATABLE :: tmp_array
! for error analysis
REAL(DP) :: total_error
! define 
REAL(DP), DIMENSION(:), ALLOCATABLE :: outdata, tmpdata,sumdata
REAL(DP), DIMENSION(:), ALLOCATABLE :: error
! define
CHARACTER(LEN=20) :: fname
CHARACTER(LEN=nlen) :: flen,flen2
!-------
! PARAMETERS  --> Subroutine
!------
call read_params
! ** Allocate arrays **
nlagmax=n-1
allocate(acf0(nlagmax+1),error(nlagmax+1))
allocate(ip(0:2+int(sqrt(real(n/2)))),w(0:n/2-1))
allocate(ip2(0:2+int(sqrt(real(n)))),w2(0:n-1))
ip(0)=0; ip2(0)=0
!----
! Define autocorrelation function & power spectrum
!----

tol=1.0_dp/a*tan((1.0_dp-eps)*PI_D/2.0_dp)
write(*,1001) atan(PI_D/dz)/PIO2_D, PI_D/dz
1001 format(f6.4, " of total spectral is contained at Nyquist ", f8.4)
write(*,1002) tol
1002 format( "recommended upper cut frequency ", f8.4)

!itol=int(tol/dk)
!dk=tol/n*2
!do while (tol > PI_D/dz)
!   n=n*2
!   dz=dz/2.0_dp
!end do
!print *, "dk", dk, "n", n, "dz",dz
!if (itol>n/2+1) itol=n/2+1



allocate(s0(n/2+1))
if (type ==1) then
   call gauss1df(s0,a,dk)
   call gauss1d(acf0,a,dz)
elseif (type==2) then
   call exp1df(s0,a,dk)
   call exp1d(acf0,a,dz)
elseif (type==3) then
   call vk1df(s0,a,dk,nu)
   call vk1d(acf0,a,dz,nu)
else 
   stop
end if


open(100,file="tmp/s0.dat")
write(100,'(F15.8)') (s0(i)/s0(2),i=1,n/2+1)
close(100)
!s0(1)=0.0_dp
!print *, s0(itol)

!s0(itol:n/2+1)=0.0_dp

!----------
! Realization
!----------
ireal=1
!call ran_seed(iseed)
call init_genrand(iseed)
allocate(sumdata(n))
sumdata=0.0_dp
200 continue
write(*,*)
write(*,*) "================================="
write(*,*) "   Realization No. ", ireal
write(*,*) "================================="
write(flen,'(i5.5)') ireal
allocate(acf(nlagmax*2+1),error(nlagmax+1)) 
allocate(tmp_array(nbin+1), itmp_array(nbin+1), theta(n/2+1))
allocate(tmpdata(n),outdata(n),sb(n/2+1))

!---
! Define Phase angle
!---
write(*,*)
write(*,*) "==== Phase ====="
theta=0.0

do i=1,n/2+1
   theta(i)=grnd()
end do

! replaced by MT
!call ran3(theta)
!do while ( dabs(sum(theta)/size(theta)-0.5_dp) > 0.00001_dp)
!   theta=theta-sum(theta)/size(theta,1)+0.50_dp
!   do i=1,n/2-1
!      if ( theta(i)<0.0_dp ) call ran3(theta(i))
!      if ( theta(i)>1.0_dp ) call ran3(theta(i))
!   end do
!end do
call calc_stat_values(theta,"angle")
call hist(theta,itmp_array,nbin)

write(fname,'("tmp/",a5,"ht.dat")') flen
open(100,file=fname)
write(100, *) (itmp_array(i), i=1,nbin)
close(100)

write(fname,'("tmp/",a5,"zt.dat")') flen
open(100,file=fname)
write(100,'(f12.8)') (theta(i)-0.50_dp, i=1,n/2-1)
close(100)
theta=theta*TWOPI_D


! **** ITERATION WILL START FROM HERE *****
allocate(sg0(n/2+1))
sg0=s0

!mu1=100.0; mu2=300.0; std1=10.0; std2=30.0; ratio1=0.8_dp


!------------------------------------------------
! GENERATE RANDOM MEDIA
!-------------------------------------------------
! * Standard Gaussian distribution
write(*,*)
write(*,*) "==== Gaussian field ===="
call std_gauss_field(tmpdata,sg0,theta(1:n/2),dk,lim,lw) 
call calc_stat_values(tmpdata,'normal distribution')
call check_gauss(tmpdata)

call hist(tmpdata,itmp_array,nbin,-lim,lim)
write(fname,'("tmp/",a5,"hg.dat")') flen
open(100,file=fname)
write(100, *) (itmp_array(i), i=1,nbin)
close(100)

sumdata=sumdata+tmpdata

! * Bi-modal distribution
write(*,*)
write(*,*) "==== Bi-modal field ===="
call cnv_bimodal(outdata,tmpdata,lim)
call calc_stat_values(outdata,'bi-modal distribution')
call hist(outdata,itmp_array,nbin,min(mu1-lim*std1,mu2-lim*std2),&
     max(mu1+lim*std1, mu2+lim*std2))
write(fname,'("tmp/",a5,"hb.dat")') flen
open(100,file=fname)
write(100, *) (itmp_array(i), i=1,nbin)
close(100)
call get_norm_power(tmpdata,sb,dk,dz)
write(fname,'("tmp/",a5,"sg.dat")') flen
open(100,file=fname)
write(100,'(F12.5)') (sb(i)/sb(2),i=1,n/2+1)
close(100)

! * calculate power spectrul normalized by its summation
call get_norm_power(outdata-sum(outdata)/n,sb,dk,dz)
write(fname,'("tmp/",a5,"zg.dat")') flen
open(100,file=fname)
write(100,'(F12.5)') (tmpdata(i),i=1,n)
close(100)

write(fname,'("tmp/",a5,"zb.dat")') flen
open(100,file=fname)
write(100,'(F12.5)') (outdata(i),i=1,n)
close(100)

write(fname,'("tmp/",a5,"sb.dat")') flen
open(100,file=fname)
write(100,'(F12.5)') (sb(i)/sb(2),i=1,n/2+1)
close(100)


!---
! CHECK ACF
!----
write(fname,'("tmp/",a5,"acf.dat")') flen
open(100,file=fname)
call acorr1df(tmpdata,acf,nlagmax)
!write(100,'(F12.5)') (acf(i),i=nlagmax+1,nlagmax*2+1)
write(100,'(F12.5)') (acf(i),i=1,nlagmax+1)

call acorr1df(outdata-sum(outdata)/n,acf,nlagmax)

!write(100,'(F12.5)') (acf(i),i=nlagmax+1,nlagmax*2+1)
write(100,'(F12.5)') (acf(i),i=1,nlagmax+1)
write(100,'(F12.5)') (acf0(i),i=1,nlagmax+1,2)
close(100)

!--
! Calculate Error 
!   (to be OPTIONAL)
!--
total_error=0.0
do i=1,nlagmax+1
   if (acf0(i)>tiny) then
      error(i)=abs((acf(i+nlagmax)-acf0(i)))/acf0(i)
   else
      error(i)=abs((acf(i+nlagmax)-acf0(i)))
   end if
   total_error=total_error+(acf(i+nlagmax)-acf0(i))**2
end do
print *, "Total Square Error", total_error
write(fname,'("tmp/",a5,"err.dat")') flen
open(100,file=fname)
write(100,'(F12.5)') (error(i),i=1,nlagmax+1)
close(100)

deallocate(acf,error) 
deallocate(tmp_array, itmp_array, theta)
deallocate(tmpdata,outdata,sb)

ireal=ireal+1
if (ireal <= nreal) goto 200

!write(fname,'("tmp/sum.dat")') flen
open(100,file="tmp/sum.dat")
write(100,'(F12.5)') (sumdata(i)/real(nreal),i=1,n)
close(100)


END PROGRAM rnd1d
