MODULE subrnd1d
  USE nrtype
  REAL(DP) :: tmp_for_min
  PRIVATE :: tmp_for_min

CONTAINS
  SUBROUTINE std_gauss_field(data,p1in,theta,dk,limit,lwin)
    !==================================================================
    ! STD_GAUSS_FIELD(DATA,P1IN,THETA,DK,LIMIT,lWIN)
    !    Generating standard gaussian field 
    !
    ! * Input data
    !    P1IN(N/2+1)  : spectral
    !    THETA(N/2+1) : phase angle
    !    DK           : dk
    !    LIMIT        : max value of output data
    !   <optional>
    !    LWIN         : length of window (%)
    !    
    ! * Output data
    !    DATA(N)      : 
    !
    ! * Required subroutines
    !    NRTYPE.F90
    !    NR.F90
    !    FFT.F90
    !    STATIC.F90
    !    RAN3.F90
    !=================================================================

    USE nrtype
    USE nr, ONLY : moment,ran3 !, realft
    USE modfft4g
    USE modfft4gsub
    USE mt19937
!    USE header_stat
    IMPLICIT NONE 
    REAL(DP),INTENT(IN) :: limit,dk
    REAL(DP),DIMENSION(:),INTENT(IN) :: p1in, theta
    REAL(DP), DIMENSION(:), INTENT(OUT) :: data
    REAL(DP), OPTIONAL, TARGET :: lwin
    REAL(DP), POINTER :: lw

    ! << Numerical Recipe
    !COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: k_data
    ! >> N.R.
    INTEGER :: new_index
    REAL(DP) :: tmp,ave,adev,sdev,var,skew,curt
    INTEGER :: i,j,n,nw
    REAL(DP), DIMENSION(:), ALLOCATABLE :: a!,w
    !INTEGER,DIMENSION(:),ALLOCATABLE :: ip
    REAL(DP), DIMENSION(:), ALLOCATABLE :: p1

    n=size(data,1)
    !allocate(k_data(n/2))
    allocate(p1(n/2+1))
    ! << Real FFT
    allocate(a(0:n-1))!,ip(0:2+int(sqrt(real(n/2)))),w(0:n/2-1))
    ! >> R.F.
    ! << Complex FFT
    ! allocate(a(0:2*n-1),ip(0:2+int(sqrt(real(n)))),w(0:n/2-1))
    ! >> C.F.
    p1=p1in
!    p1(1)=0.0_dp
!    print *, " applying window"
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
       nw=int(n/2*lw)
     !  print *, "nw", nw
       p1(n/2+1)=0.0_dp
       do i=1,nw
          p1(n/2+1-i)=sin(i/nw*PIO2_D)*p1(n/2+1-i)
       enddo
    end if
    ! << Numerical Recipe
    !k_data(2:n/2)=sqrt(p1(2:n/2))*exp(cmplx(0.0_dp,1.0_dp)*theta(1:n/2-1))
    !k_data(1)=cmplx(0.0_dp,sqrt(p1(n/2+1)))
    !k_data=k_data*2.0_dp*sqrt(dk)
    ! >> N.R.
    !---
    ! FFT and return to x-domain.
    !---
    data=0.d0
    ! Uncomment following lines when using fft4g.f90
    ! << Real FFT
    a(0)=0.0_dp
    a(1)=p1(n/2+1)
    do i=2,n/2
       a(2*(i-1))=sqrt(p1(i))*cos(theta(i))
       a(2*(i-1)+1)=sqrt(p1(i))*sin(theta(i))
    end do
    !ip(0)=0
    a=a*sqrt(dk)
    !a=sqrt(dk)
    call rdft(n,-1,a,ip,w)
    data(1:n)=a(0:n-1)

    !<< Compex FFT (fft4g.f90)
    !a(0:1)=0.0_dp
    !do i=2,n/2
    !   a(2*(i-1))=sqrt(p1(i))*cos(theta(i))
    !   a(2*(i-1)+1)=sqrt(p1(i))*sin(theta(i))
    !   a(2*(n-i+1))=sqrt(p1(i))*cos(theta(i))
    !   a(2*(n-i+1)+1)=-sqrt(p1(i))*sin(theta(i))
    !end do
    !a(n)=sqrt(p1(n/2+1))
    !a(n+1)=0.0_dp
    !ip(0)=0
    !call cdft(2*n,-1,a,ip,w)
    !data(1:n)=a(1:2*n-1:2)
    ! >>
    !<< Numerical Recipe
    ! call realft(data,-1,k_data)
    ! deallocate(k_data)!,a,ip,w)
    ! << 
    data=2.0_dp*data/sqrt(TWOPI_D)
    !data=data/TWOPI_D
    !call calc_stat_values(data,'gaussian distribution')

    !--
    ! Normalization
    !--
    call moment(data,ave,adev,sdev,var,skew,curt)
    data=(data-ave)/sdev ! modify to normal distribution
    !---
    ! It's better to be cut larger than 4*std
    !---
    !print *, ""
    !print *, "Max and Min of temporary data"
    !print *, maxval(data), minval(data)
    do i=1,n
       !call init_genrand(100)
       do while ( abs(data(i)) > limit ) 
!          call ran3(tmp)
          new_index=int(grnd()*n)
          if (new_index == 0) new_index=1
          print *, "Value is grater than ", limit
          print *, i, data(i), "-->",new_index, data(new_index)
          data(i)=data(new_index)
       enddo
    enddo
  END SUBROUTINE std_gauss_field

  SUBROUTINE cnv_bimodal(bidata,orgdata,limit)
    USE header_stat
    USE nr, ONLY : rtbis
!    USE acorr, ONLY : hist
    USE acflib, ONLY : cpd_std_gauss
    IMPLICIT NONE 
    REAL(DP),DIMENSION(:), INTENT(IN) :: orgdata
    REAL(DP),DIMENSION(:), INTENT(OUT) :: bidata
    REAL(DP), INTENT(IN) :: limit
    REAL(DP) :: bndlow,bndhigh
    INTEGER :: i,n
    n=size(bidata,1)
    bndlow=min(mu1-std1*limit, mu2-std2*limit)
    bndhigh=max(mu1+std1*limit,mu2+std2*limit)
   ! convert to Bi-Peak distribution
    do i=1,n
       tmp_for_min=cpd_std_gauss(orgdata(i))
       if ( tmp_for_min == 1 ) then
          bidata(i)=bndlow
       else
          bidata(i)=rtbis(func2min,bndlow,bndhigh,1.0d-3)
       end if
    enddo
  END SUBROUTINE cnv_bimodal
  
  SUBROUTINE calc_psd(data,power,dk,dz)
    !USE nr, ONLY : realft
    USE modfft4g
    USE modfft4gsub
    IMPLICIT NONE
    REAL(DP) :: dz,dk
    REAL(DP), DIMENSION(:), INTENT(IN):: data
    REAL(DP), DIMENSION(:), INTENT(OUT) :: power
    REAL(DP), DIMENSION(:), ALLOCATABLE :: a
    ! << N.R
    !COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: k_data
    ! >> N.R.
    INTEGER :: n
    n=size(data,1)
    ! << N.R. 
    ! allocate(k_data(n/2))
    ! >> N.R.
    allocate(a(0:n-1))
    a(0:n-1)=data(1:n)
    !ip(0)=0
    call rdft(n,1,a,ip,w)
    power(1)=a(0)**2
    power(2:n/2)=a(2:n-2:2)**2+a(3:n-1:2)**2
    power(n/2+1)=a(1)**2
    power=power*dz/n/2
!    print *, dz/n, maxval(power)
    
    deallocate(a)
    ! << N.R.
    !call realft(data,1,k_data)
    !power(1)=real(k_data(1))
    !power(2:n/2)=cdabs(k_data(2:n/2))
    !power(n/2+1)=imag(k_data(1))
    !power=power**2
    !power=power/power(2)
    !deallocate(k_data)
    ! >> N.R.
  END SUBROUTINE calc_psd

  FUNCTION func2min(x)
    USE header_stat
    USE acflib, ONLY : cpd_bi_gauss
    REAL(DP) :: x
    REAL(DP) :: func2min
    func2min=cpd_bi_gauss(x,mu1,std1,mu2,std2,ratio1)&
         -tmp_for_min
  END FUNCTION func2min
  
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
    REAL(DP) :: value,exdprb(5)
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

END MODULE subrnd1d
