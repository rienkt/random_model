PROGRAM rnd1d
! Creating 1-D random field from auto-correlation function.
! Not minimize memory costs.
!
USE acorr
USE nrtype
USE nr, ONLY : realft, ran3,moment,rtbis
USE ran_state
USE header_stat
USE subrnd1d
IMPLICIT NONE
real(DP) :: tiny,large,limit,lambda
PARAMETER(tiny=0.01_dp,large=5.0_dp,limit=4.0_dp,lambda=0.1_dp) 
! define general temporary variables 
INTEGER :: i,j,m, index, ipos,ipos2
REAL(DP) :: tmp,tmp2
REAL(DP),DIMENSION(:), ALLOCATABLE :: tmp_array
! fpr error analysis
REAL(DP) :: ave,adev,sdev,var,skew,curt, total_error
! define 
INTEGER type,clagmax,n,iseed,iter
REAL(DP), DIMENSION(:), ALLOCATABLE :: outdata, tmpdata
REAL(DP), DIMENSION(:), ALLOCATABLE :: auto_corr
REAL(DP), DIMENSION(:), ALLOCATABLE ::  tr_auto_corr, error
COMPLEX(DPC), DIMENSION(:), ALLOCATABLE  :: k_data,k2_data
REAL(DP), DIMENSION(:), ALLOCATABLE :: p1,p1_new,p1t
REAL(DP), DIMENSION(:), ALLOCATABLE   :: theta
REAL(DP) :: nu, a,dk, dz
! distribution
!REAL(DP) :: mu1, mu2, std1, std2, ratio1

! define for bi-peak conversion
REAL(DP), DIMENSION(:), ALLOCATABLE :: sb,sg0,sg1,sg2,st, serr1, serr2
REAL(DP), DIMENSION(:), ALLOCATABLE :: dev
REAL(DP) :: lambda1(3),tserr1,tserr2,tserr(3)

!-------
! PARAMETERS  --> Subroutine
!------
call read_params(n,dz,dk,type,nu,a,clagmax,iseed)
! ** Allocate arrays **
allocate(p1(n/2+1), st(n/2+1),theta(n/2-1))
allocate(auto_corr(clagmax*2+1), tr_auto_corr(clagmax+1), &
     error(clagmax+1))

!----
! Define autocorrelation function & power spectrum
!----
m=n*2
allocate(p1t(m))

if (type ==1) then
   call gauss1df(p1,a,dk)
   call gauss1df(p1t,a,dk)
   call gauss1d(tr_auto_corr,a,dz)
elseif (type==2) then
   call exp1df(p1,a,dk)
   call exp1df(p1t,a,dk)
   call exp1d(tr_auto_corr,a,dz)
elseif (type==3) then
   call vk1df(p1,a,dk,nu)
   call vk1df(p1t,a,dk,nu)
   call vk1d(tr_auto_corr,a,dz,nu)
else 
   stop
end if

open(100,file="p1t.dat")
write(100,'(F12.5)') (p1t(i)/p1t(1),i=1,m)
close(100)

open(100,file="p1.dat")
write(100,'(F12.5)') (p1(i)/p1(1),i=1,n/2+1)
!write(100,'(F12.5)') ((1.0_dp-p1t(i)),i=1,n/2+1)
p1t=p1t/sum(p1t)
tmp2=1.0_dp-p1t(1)
write(100,'(F12.5)') tmp2
ipos=n/2+1
do i=2,n/2
   tmp2=tmp2-p1t(i)
   write(100,'(F12.5)') tmp2
   if (tmp2 >2d-1) ipos2=i
   if ( tmp2 < 1d-1 ) then
      ipos=i
      exit
   end if
enddo
print *, "IPOS",ipos,tmp2,ipos2


close(100)


!---
! Define k-domain components.
!---
theta=0.0
call ran_seed(iseed)
call ran3(theta)
print *, "Statistical values of Phase/2/PI"
call moment(theta,ave,adev,sdev,var,skew,curt)
print *, "ave", ave, "adev", adev,"sdev",sdev,"var",var
print *, "skew",skew,"curt",curt
print *, ""
theta=theta*TWOPI


! **** ITERATION WILL START FROM HERE *****
allocate(sg0(n/2+1),sg1(n/2+1),sg2(n/2+1))
p1(1)=0.0_dp; st=p1/sum(p1(2:n/2+1)); sg1=st
iter=1
open(201,file="sg.dat")
open(202,file="sb.dat")
write(201,'(F15.8)') (st(i)/sum(st(2:n/2+1)),i=2,n/2+1)
write(202,'(F15.8)') (st(i)/sum(st(2:n/2+1)),i=2,n/2+1)
open(301,file="ag.dat")
open(302,file="ab.dat")
write(301,'(F12.5)') (tr_auto_corr(i),i=1,clagmax+1)
write(302,'(F12.5)') (tr_auto_corr(i),i=1,clagmax+1)


allocate(serr1(n/2+1),serr2(n/2+1))!,lmbd(n/2+1))
serr1=10000.0_dp ! initialize
!lmbd=lambda
!lmbd(ipos2+1:n/2+1)=0.0_dp
   mu1=100.0; mu2=300.0; std1=10.0; std2=30.0; ratio1=0.8_dp

lambda1(1)=lambda/10.0_dp;lambda1(2)=lambda;lambda1(3)=lambda*10.0_dp

!*****
! Calculate Derivatives (Oh No.......)
!****
allocate(dev(n/2+1))
allocate(tmpdata(n),outdata(n),sb(n/2+1))
dev=0.0_dp
! * Standard Gaussian distribution
call std_gauss_field(tmpdata,sg1,theta,limit) 
! * Bi-Peak distribution
call cnv_bipeak(outdata,tmpdata,limit)
! * calculate power spectrul normalized by its summation
call moment(outdata,ave,adev,sdev,var,skew,curt)
call get_norm_power(outdata-ave,sb)
tserr1=sum(abs(sb(2:ipos)-st(2:ipos)))/ipos
sg2=sg1
do i=2,ipos
   sg2(i)=sg1(i)*1.05
   call std_gauss_field(tmpdata,sg2,theta,limit) 
   ! * Bi-Peak distribution
   call cnv_bipeak(outdata,tmpdata,limit)
   ! * calculate power spectrul normalized by its summation
   call moment(outdata,ave,adev,sdev,var,skew,curt)
   call get_norm_power(outdata-ave,sb)
   !* calculate error function using spectral information
   tserr2=sum(abs(sb(2:ipos)-st(2:ipos)))/ipos
   dev(i)=(tserr2-tserr1)/(sg2(i)-sg1(i))
   write (102,*) tserr2
!   print *, sg1(i)
   print *, sg2(i)-sg1(i),tserr2-tserr1,dev(i)
   sg2(i)=sg1(i)
enddo
deallocate(sb,tmpdata,outdata)
tserr1=0.0_dp
tserr2=0.0_dp

open(102,file="error_sp.dat")
do iter=1,3
   !***********
   !  Forward Modeling 1
   !***********
   allocate(tmpdata(n),outdata(n),sb(n/2+1))
   ! * Standard Gaussian distribution
   call std_gauss_field(tmpdata,sg1,theta,limit) 
   ! * Bi-Peak distribution
   call cnv_bipeak(outdata,tmpdata,limit)
   ! * calculate power spectrul normalized by its summation
   call moment(outdata,ave,adev,sdev,var,skew,curt)
   call get_norm_power(outdata-ave,sb)
   write(201,'(F15.8)') (sg1(i),i=2,n/2+1)
   write(202,'(F15.8)') (sb(i),i=2,n/2+1)
   !* calculate error function using spectral information
   serr2=0.0
   serr2=abs(sb-st)
   tserr2=sum(serr2(2:ipos))/ipos
   write (102,*) tserr2

   if ( tserr2 < 1e-4 ) exit
   !**********
   ! Decide Lambda
   !**********
   do i=1,3
   ! Decide Lambda
      do j=2,n/2+1
         !         if ( serr2(j) < serr1(j)) then
         if (st(j) > 1d-5) then
            tmp=(sg1(j)-sb(j))/sb(j)
            if (tmp > 1.0) tmp=1.0
            sg2(j)=st(j)-dev(j)*lambda1(i)*st(j)
         else
            sg2(j)=st(j)
         endif
         !         else
         !            sg2(j)=sg0(j)
         !         endif
      enddo
      ! * Standard Gaussian distribution
      call std_gauss_field(tmpdata,sg2,theta,limit) 
      ! * Bi-Peak distribution
      call cnv_bipeak(outdata,tmpdata,limit)
      ! * calculate power spectrul normalized by its summation
      call moment(outdata,ave,adev,sdev,var,skew,curt)
      call get_norm_power(outdata-ave,sb)
      !* calculate error function using spectral information
      tserr(i)=sum(abs(sb(2:ipos)-st(2:ipos)))/ipos
   enddo
   
   write (103,*) tserr(1),tserr(2),tserr(3)
   lambda1(2)=log(lambda1(2))-0.5_dp*&
        ( &
        (log(lambda1(2))-log(lambda1(1)))**2*(tserr(2)-tserr(1))&
        -(log(lambda1(2))-log(lambda1(3)))**2*(tserr(2)-tserr(3))&
        )/&
        ( &
        (log(lambda1(2))-log(lambda1(1)))*(tserr(2)-tserr(1))&
        -(log(lambda1(2))-log(lambda1(3)))*(tserr(2)-tserr(3))&
        ) 

!   lambda1(2)=lambda1(2)-0.5_dp*&
!        ( &
!        (lambda1(2)-lambda1(1))**2*(tserr(2)-tserr(1))&
!        -(lambda1(2)-lambda1(3))**2*(tserr(2)-tserr(3))&
!        )/&
!        ( &
!        (lambda1(2)-lambda1(1))*(tserr(2)-tserr(1))&
!        -(lambda1(2)-lambda1(3))*(tserr(2)-tserr(3))&
!        ) 

   lambda1(2)=exp(lambda1(2))
   lambda1(1)=lambda1(2)/10.0_dp
   lambda1(3)=lambda1(2)*10.0_dp
!   if ( lambda1(3)>1.0 ) lambda1(3)=0.9_dp
   print *, lambda1


   !* update spectral information p1 and goes back to pos. 100
   sg2=sg1
!   if ((iter == 6).and.(ipos2.ne. ipos)) lmbd(ipos2+1:ipos)=lambda1(2)
   do i=2,ipos
      if ( serr2(i) < serr1(i)) then
         if (st(i) > 1d-4) then
            tmp=(sg1(i)-sb(i))/sb(i)
            if (tmp > 1.0) tmp=1.0
            sg2(i)=st(i)*(1.0_dp+tmp*lambda1(2))
         else
            sg2(i)=st(i)
         endif
      else
         sg2(i)=sg0(i)
      endif
   enddo
   print *,"P1", sg2(ipos), sg1(ipos),sg0(ipos)
   !      print *, p1_new(i),sg(i),sb(i)
   sg0=sg1
   sg1=sg2

   call acorr1d(tmpdata,auto_corr,clagmax)
   write(301,'(F12.5)') (auto_corr(i),i=clagmax+1,clagmax*2+1)
   call moment(outdata,ave,adev,sdev,var,skew,curt)
   call acorr1d(outdata-ave,auto_corr,clagmax)
   write(302,'(F12.5)') (auto_corr(i),i=clagmax+1,clagmax*2+1)


   total_error=0.0
   do i=1,clagmax+1
      if (tr_auto_corr(i)>tiny) then
         error(i)=abs((auto_corr(i+clagmax)-tr_auto_corr(i)))/tr_auto_corr(i)
      else
!         error(i)=abs((auto_corr(i+clagmax)-tr_auto_corr(i)))
      end if
      total_error=total_error+(auto_corr(i+clagmax)-tr_auto_corr(i))**2
   end do
   print *, "Total Square Error", total_error


   if (iter < 3) deallocate(sb,tmpdata,outdata)
end do

!   open(100,file="tmpdata.dat")
!   write(100,'(F12.5)') (tmpdata(i),i=1,n)
!   close(100)


open(100,file="outdata.dat")
write(100,'(F12.5)') (outdata(i),i=1,n)
close(100)


!---
! Finaly,
! Calculate Autocorrelation Function 
!  (to be OPTIONAL)
!
!----
if (type ==1) then
   call gauss1d(tr_auto_corr,a,dz)
elseif (type==2) then
   call exp1d(tr_auto_corr,a,dz)
!   call exp1df(p1,a,dk)
elseif (type==3) then
   print *, "TYPE=3"
   call vk1d(tr_auto_corr,a,dz,nu)
!   call vk1df(p1,a,dk)
else 
   stop
end if
open(100,file="acorr.dat")
call acorr1d(tmpdata,auto_corr,clagmax)
write(100,'(F12.5)') (auto_corr(i),i=clagmax+1,clagmax*2+1)
print *, ""
print *, "Stastical values of normalized data after converting"
call moment(outdata,ave,adev,sdev,var,skew,curt)
print *, "ave", ave, "adev", adev,"sdev",sdev,"var",var
print *, "skew",skew,"curt",curt
call acorr1d(outdata-ave,auto_corr,clagmax)

write(100,'(F12.5)') (auto_corr(i),i=clagmax+1,clagmax*2+1)
write(100,'(F12.5)') (tr_auto_corr(i),i=1,clagmax+1)
close(100)

!--
! Calculate Error 
!   (to be OPTIONAL)
!--
total_error=0.0
do i=1,clagmax+1
!   print *, tr_auto_corr(i), auto_corr(i+clagmax),tr_auto_corr(i)
   if (tr_auto_corr(i)>tiny) then
      error(i)=abs((auto_corr(i+clagmax)-tr_auto_corr(i)))/tr_auto_corr(i)
   else
      error(i)=abs((auto_corr(i+clagmax)-tr_auto_corr(i)))
   end if
   total_error=total_error+(auto_corr(i+clagmax)-tr_auto_corr(i))**2
!   print *, tr_auto_corr(i), auto_corr(i+clagmax),tr_auto_corr(i),error(i)
end do
print *, "Total Square Error", total_error
open(100,file="error.dat")
write(100,'(F12.5)') (error(i),i=1,clagmax+1)
close(100)
CONTAINS
  SUBROUTINE read_params(n,dz,dk,type,nu,a,clagmax,iseed)
    INTEGER  :: type, clagmax, iseed
    REAL(DP) :: dz,dk,nu,a
    REAL(DP) :: clagmax_value
    INTEGER :: i,j,n
    REAL(DP) :: tmp
    
    print *, "n, dz"
!tmp    read (*,*) n, dz
    read (*,*) n, dz, type,nu, a, clagmax_value, iseed
    dk=1.0_dp/dz/float(n)*TWOPI ! Nyquist : 1/2dz >> kan
    print *, "n=", n, "dz=", dz, "dk=",dk
    print *, "Nyqist"

    ! ** Parameters related to auto correlation function
10  print *, "Type [1] Gaussian [2] Exponential [3] Von Karman ";
!tmp    read (*,*) type
           print *, "Nu "
!tmp       read (*,*) nu
    if (type == 3) then
       print *, "Nu "
    elseif (type == 2) then
       nu=0.5_dp
    elseif (type == 1) then
       nu=100.0_dp
    else
       print *, "Wrong type"
       goto 10
    end if
11  print *, "Auto correlation length";
!tmpread(*,*) a

    print *, "Lagmax for calculation maximum"
!tmpread(*,*) clagmax_value
    clagmax=int(clagmax_value/dz)
    
    print *, "Random seed"
!tmpread(*,*) iseed
  END SUBROUTINE read_params

END PROGRAM rnd1d
