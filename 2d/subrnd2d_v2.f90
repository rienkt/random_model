MODULE subrnd2d
! This module is now as same as subrnd2d_hpc.f90
  USE nrtype
  REAL(DP) :: tmp_for_min
  PRIVATE :: tmp_for_min

CONTAINS
  SUBROUTINE std_gauss_field(data,p1in,ttheta,limit,lwin)
    !==================================================================
    ! STD_GAUSS_FIELD(DATA,P1IN,THETA,DK,LIMIT,lWIN)
    !    Generating standard gaussian field 
    !
    ! * Input data
    !    P1IN(NX/2+1,NZ)  : spectral
    !    THETA(NX/2+1,NZ) : phase angle
    !    DKX,DKX      : dk
    !    LIMIT        : max value of output data
    !   <optional>
    !    LWIN         : length of window (%)
    !    
    ! * Output data
    !    DATA(NX,NZ)      : 
    !
    ! * Required subroutines
    !    NRTYPE.F90
    !    NR.F90
    !    FFT.F90
    !    STATIC.F90
    !    RAN3.F90
    !=================================================================

    USE nrtype
    USE nr, ONLY : moment
    USE fft_params_2d
    !USE fft_call_2d
    USE mt19937
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: limit
    REAL(DP), DIMENSION(:,:),INTENT(IN) :: p1in,ttheta
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: data
    REAL(DP), OPTIONAL, TARGET :: lwin
    REAL(DP), POINTER :: lw

    INTEGER :: new_index
    REAL(DP) :: ave,sdev!tmp,ave,adev,sdev,var,skew,curt
    INTEGER :: i,j,mx,mz,nw
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a,p1
!    REAL(DP), DIMENSION(:), ALLOCATABLE :: tmparray
    mx=size(data,1)
    mz=size(data,2)
!    print *,'ahou'
    allocate(p1(mx/2+1,mz))
    allocate(a(0:mx+1,0:mz-1))
!    print *, 'windowing'
    p1=p1in
    if (present(lwin)) then
       lw=>lwin
    else
       allocate(lw)
       lw=0.1
    end if

    !--
    ! Window
    !--
    if (lw > 0) then
       nw=int(min(mx,mz)/2*lw)
       p1(mx/2+1,mz/2+1)=0.0_dp
       do i=1,nw
          p1(mx/2+1-i,1:mz)=sin(dble(i)/nw*PIO2_D)*p1(mx/2+1-i,1:mz)
       enddo
    end if


!    print *, 'preparing data'
    !--
    ! Prepare data
    !--
    data=0.d0
    p1=sqrt(p1)
    a(0,0)=0.0_dp
    do i=1,mx/2+1
       do j=1,mz
          a(2*(i-1),j-1)=p1(i,j)*cos(ttheta(i,j))
          a(2*(i-1)+1,j-1)=p1(i,j)*sin(ttheta(i,j))
!          a(2*(i-1),j-1)=p1(i,j)*cos(data(i,j))
!          a(2*(i-1)+1,j-1)=p1(i,j)*sin(data(i,j))
       enddo
    enddo
!    a(1,0)=p1(nx/2+1,1)
!    a(0,nz/2)=p1(1,nz/2+1)
!    a(1,nz/2)=p1(nx/2+1,nz/2+1)
!    do i=2,nx/2
!       do j=2,nz
!          a(2*(i-1),j-1)=p1(i,j)*cos(theta(i,j))
!          a(2*(i-1)+1,j-1)=p1(i,j)*sin(theta(i,j))
!       enddo
!       a(2*(i-1),0)=p1(i,1)*cos(i,1)
!       a(2*(i-1)+1,0)=p1(i,1)*sin(i,1)
!    enddo
!    do j=2,nz/2
!       a(0,j-1)=p1(1,j)*cos(1,j)
!       a(1,j-1)=p1(1,j)*sin(1,j)
!       a(1,nz-j+1)=p1(nx/2+1,j)*cos(nx/2+1,j)
!       a(0,nz-j+1)=-p1(nx/2+1,j)*sin(nx/2+1,j)
!    enddo

!    a=a*sqrt(dkx)*sqrt(dkz)
    call rdft2dsort(mx+2,mx,mz,-1,a(0:mx+1,0:mz-1))
    call rdft2d(mx+2,mx,mz,-1,a(0:mx+1,0:mz-1),t,ip,w)
    data(1:mx,1:mz)=a(0:mx-1,0:mz-1)
!    print *, data(100,100),a(40,40)
!    data=2.0_dp*data/sqrt(TWOPI_D)
    !data=data/TWOPI_D
    !call calc_stat_values(data,'gaussian distribution')

    !--
    ! Normalization
    !--
    print *, 'normalization'
    !allocate(tmparray(nx*nz))
    !call cnv1d2d(tmparray,data,-1)
    !print *, 'normalization2',size(tmparray,1)
    !call moment(tmparray,ave,adev,sdev,var,skew,curt)
    !print *, 'normalization3'
    !print *, sdev,ave
    ave=0; sdev=0;
    !print *, maxval(data), minval(data)
    call calc_mean_var(data,ave,sdev)
    !print *, sdev,ave
    data=(data-ave)/sdev ! modify to normal distribution
    !---
    ! It's better to be cut larger than 4*std
    !---
    !print *, ""
    !print *, "Max and Min of temporary data"
    !print *, maxval(data), minval(data)

    !print *, 'abon'
    do i=1,mx
       do j=1,mz
          do while ( abs(data(i,j)) > limit ) 
             !          call ran3(tmp)
             new_index=int(grnd()*nz)
!             print *, new_index
             if (new_index == 0) new_index=1
             !print *, "Value is grater than ", limit
             !print *, i, j, data(i,j), "-->",new_index, data(i,new_index)
             data(i,j)=data(i,new_index)
          enddo
       enddo
    enddo
    deallocate(p1,a)
  END SUBROUTINE std_gauss_field

  SUBROUTINE calc_mean_var(data,ave,sdev)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: data
    REAL(DP), INTENT(OUT) :: ave,sdev
    INTEGER :: i,j,nx,nz
    nx=size(data,1)
    nz=size(data,2)
    ave=0.0d0
    do i=1,nx
       do j=1,nz
          ave=ave+data(i,j)
       enddo
    enddo
    ave=ave/float(nx)/float(nz)
    sdev=0.0d0
    do i=1,nx
       do j=1,nz
          sdev=sdev+(data(i,j)-ave)*(data(i,j)-ave)
       enddo
    enddo
    sdev=sdev/float(nx)/float(nz)
    sdev=dsqrt(sdev)
  end SUBROUTINE calc_mean_var
    

  SUBROUTINE cnv_bimodal(bidata,orgdata,limit)
    USE header_stat
    USE nr, ONLY : rtbis
!    USE acorr, ONLY : hist
    USE acflib, ONLY : cpd_std_gauss
    IMPLICIT NONE 
    REAL(DP),DIMENSION(:,:), INTENT(IN) :: orgdata
    REAL(DP),DIMENSION(:,:), INTENT(OUT) :: bidata
    REAL(DP), INTENT(IN) :: limit
!    REAL(DP), DIMENSION(:), ALLOCATABLE :: tmparray
    REAL(DP) :: bind1, bind2
    INTEGER :: i,j,nx,nz
    nx=size(bidata,1)
    nz=size(bidata,2)
   ! convert to Bi-Peak distribution
    bind1=min(mu1-std1*limit,mu2-std2*limit)
    bind2=max(mu1+std1*limit,mu2+std2*limit)
!    print *, bind1, bind2
!OCL INDEPENDENT(out)
    do i=1,nx
       do j=1,nz
!          tmp_for_min=cpd_std_gauss(orgdata(i,j))
!          print *, tmp_for_min
!          if ( tmp_for_min == 1 ) then
!             bidata(i,j)=bind2
!          else
!             bidata(i,j)=rtbis(func2min,bind1,bind2,1.0d-2)
!          end if
          bidata(i,j)=out(orgdata(i,j),bind1,bind2,1.0d-2)
       enddo
    enddo
  END SUBROUTINE cnv_bimodal
  
  RECURSIVE FUNCTION test(in)
    real(DP):: in,test
    test=in+1.0d0
  end FUNCTION test

  RECURSIVE FUNCTION out(in,bind1,bind2,lim)
    USE header_stat
    USE nr, ONLY : rtbis
    USE acflib, ONLY : cpd_std_gauss, cpd_bi_gauss
    IMPLICIT NONE 
    REAL(DP) :: in, bind1, bind2,lim,out,mine,ddx,xmid,fmid,f
    INTEGER :: j
    mine=cpd_std_gauss(in)
    if ( mine == 1 ) then
       out=bind2
    else
       !out=rtbis(func2min,bind1,bind2,lim)
       fmid=cpd_bi_gauss(bind2,mu1,std1,mu2,std2,ratio1)-mine
       f=cpd_bi_gauss(bind1,mu1,std1,mu2,std2,ratio1)-mine
       !       if (f*fmid >= 0.0) call nrerror('rtbis: root must be bracketed')
       if (f < 0.0) then ! Orient the search so that f>0 lies at x+dx.
          out=bind1
          ddx=bind2-bind1
       else
          out=bind2
          ddx=bind1-bind2
       end if
       do j=1,10000! Bisection loop.
          ddx=ddx*0.5_dp
          xmid=out+ddx
          fmid=cpd_bi_gauss(xmid,mu1,std1,mu2,std2,ratio1)-mine
          if (fmid <= 0.0) out=xmid
          if (abs(ddx) < lim .or. fmid == 0.0) goto 77
       end do
    end if
77  continue
  end FUNCTION out

  SUBROUTINE calc_psd(data,power,ddz)
    USE fft_params_2d
    !USE fft_call_2d
    IMPLICIT NONE
    REAL(DP) :: ddz!,dk
    REAL(DP), DIMENSION(:,:), INTENT(IN):: data
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: power
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a
    INTEGER :: mx,mz,i,j
    print *, 'check0'
    mx=size(data,1)
    mz=size(data,2)
    print *, mx,mz
    allocate(a(0:mx+1,0:mz-1))
    a(0:mx-1,0:mz-1)=data(1:mx,1:mz)
    !print *, 'check1'
    call rdft2d(mx+2,mx,mz,1,a,t,ip,w)
    !print *, 'check2'
    call rdft2dsort(mx+2,mx,mz,1,a)
    !print *, 'check3'
    do i=1,mx/2+1
       do j=1,mz
          power(i,j)=a(2*(i-1),j-1)**2+a(2*(i-1)+1,j-1)**2
       enddo
    enddo
    power(1,1)=0.25_dp*(power(2,1)+power(1,2)+power(2,2)+power(2,mz))
    deallocate(a)
    allocate(a(0:mx/2+2,0:mz+1))
    a(1:mx/2+1,1:mz)=power
    a(0,2:mz)=a(2,mz:2:-1)
    a(mx/2+2,2:mz)=a(mx/2,mz:2:-1)
    a(0,1)=a(2,1)
    a(mx/2+2,1)=a(mx/2,1)
    a(:,0)=a(:,mz)
    a(:,mz+1)=a(:,1)
    do i=1,mx/2+1
       do j=1,mz
          power(i,j)=1.0_dp/9.0_dp*(a(i-1,j)+a(i,j)+a(i+1,j)+a(i,j+1)+a(i,j-1)+a(i-1,j-1)+a(i+1,j+1)+a(i-1,j+1)+a(i+1,j-1))
       enddo
    enddo
    power=power*ddz**2/mx/mz/2.0_dp
    !print *, dz**2/mx/mz, maxval(power)    
    deallocate(a)
  END SUBROUTINE calc_psd

!  RECURSIVE FUNCTION func2min(x)
!    USE header_stat
!    USE acflib, ONLY : cpd_bi_gauss
!!   REAL(DP) :: x
!    REAL(DP) :: func2min
!    func2min=cpd_bi_gauss(x,mu1,std1,mu2,std2,ratio1)&
!         -tmp_for_min
! END FUNCTION func2min
  
  SUBROUTINE calc_stat_values(data,comment)
    USE nr, ONLY : moment
    REAL(DP), DIMENSION(:), INTENT(IN) :: data
    CHARACTER(LEN=*) :: comment
    REAL(DP) :: ave,adev,sdev,var,skew,curt
    call moment(data,ave,adev,sdev,var,skew,curt)
    write(*,*)
    write(*,*) "Stastical values of ", comment
    write(*,1001) ave,var
    write(*,1002) adev,sdev
    write(*,1003) skew,curt
1001 format("ave =", f9.4, " var =",f9.4)
1002 format("adev=", f9.4, " sdev=", f9.4)
1003 format("skew=", f9.4, " curt=",f9.4)
  END SUBROUTINE calc_stat_values

  SUBROUTINE check_gauss(indata)
    USE nr, ONLY : sort
    IMPLICIT NONE
    REAL(DP),  DIMENSION(:), INTENT(IN) :: indata ! mean=0, std=1
    REAL(DP), DIMENSION(:), ALLOCATABLE :: data
    REAL(DP) :: exdprb(5)
    INTEGER :: n,nstd,i
    n=size(indata,1)
    allocate(data(n))
    data=indata
    call sort(data)
    nstd=5
    do while (data(n) < nstd)
       exdprb(nstd)=0.0
       nstd=nstd-1
    enddo
    do i=1,n
       if (data(n-i+1) <= nstd) then
          exdprb(nstd)=float(i)/float(n)
          nstd=nstd-1
          if (nstd < 1) goto 100
       endif
    end do
100 continue
    
    print *, "exceedence probability of pseudo-gaussian field"
    write(*,'(f12.5)') (exdprb(i),i=1,5)
  END SUBROUTINE check_gauss


  SUBROUTINE cnv1d2d(a,b,flag)
    USE nrutil, ONLY : nrerror
    ! 1d -> 2d when flag==1
    ! 2d -> 1d when flag==-1
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: a
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: b
    INTEGER, INTENT(IN) :: flag
    INTEGER :: n1,n2,i,j
    n1=size(b,1)
    n2=size(b,2)
    if ( flag==1) then
       do j=1,n2
          do i=1,n1
             b(i,j)=a((j-1)*n1+i)
          enddo
       enddo
    elseif (flag==-1) then
       !print *, "matrix to array"
       do j=1,n2
          do i=1,n1
             !print *, i,j,n1,n2,(j-1)*n1+j,b(i,j)
             a((j-1)*n1+i)=b(i,j)
          enddo
       enddo
    else 
       call nrerror('cnv1d2d: invalid value for flag')
    end if

  END SUBROUTINE cnv1d2d


END MODULE subrnd2d
