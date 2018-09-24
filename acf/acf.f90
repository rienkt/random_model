MODULE acflib
! Subroutines related to autocorrelation function.
! 
!
!
  USE nrtype
CONTAINS
!==================
! 1-D autocorrelation functions.
!    Gaussian, Exponential, and von Karman
!==================
  SUBROUTINE gauss1d(xnr,a,dr) ! Function would be much better
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: xnr
    REAL(DP) :: dr,a,r
    INTEGER :: i,n
    n=size(xnr)
    do i=1,n
       r=float(i-1)*dr
       xnr(i)=dexp(-r**2/a**2)
    enddo
  END SUBROUTINE gauss1d

  SUBROUTINE gauss1df(p1,a,dk)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: p1
    REAL(DP) :: a,dk,k
    INTEGER :: i,n
    n=size(p1)
    do i=1,n
       k=float(i-1)*dk
       p1(i)=dsqrt(PI_D)*a*dexp(-(k*a)**2/4.0_dp)
    enddo
  END SUBROUTINE gauss1df

  SUBROUTINE exp1d(xnr,a,dr)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: xnr
    REAL(DP) :: dr,a,r
    INTEGER :: i,n
    n=size(xnr)
    do i=1,n
       r=float(i-1)*dr
       xnr(i)=dexp(-(float(i-1)*dr)/a)
    enddo
  END SUBROUTINE exp1d

  SUBROUTINE exp1df(p1,a,dk)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: p1
    REAL(DP) :: a,dk,k
    INTEGER :: i,n
    n=size(p1)
    do i=1,n
        k=float(i-1)*dk
       p1(i)=2.0_dp*a/(1.0_dp+k**2*a**2)
 !      print *, i,p1(i),k,a
    enddo
  END SUBROUTINE exp1df

   SUBROUTINE vk1d(xnr,a,dr,nu)
    USE nr, ONLY : bessik, gammln
    IMPLICIT NONE 
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: xnr
    REAL(DP), INTENT(IN) :: dr,a,nu
    REAL(DP) :: r,k_nu,dummy
    INTEGER :: i,n
    n=size(xnr)
    xnr(1)=1.0
    do i=2,n
       r=float(i-1)*dr
       call bessik(r/a,nu,dummy,k_nu,dummy,dummy)
       xnr(i)=2.0**(1.0-nu)/exp(gammln(nu))*(r/a)**nu*k_nu
    enddo
  END SUBROUTINE vk1d

  SUBROUTINE vk1df(p1,a,dk,nu)
    USE nr, ONLY : gammln
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: p1
    REAL(DP), INTENT(IN) :: a,dk,nu
    REAL(DP) :: k
    INTEGER :: i,n
    n=size(p1)
!    print *, a,dk,nu,n
!    p1(1)=0.0_dp
    do i=1,n
       k=float(i-1)*dk
!       print *, exp(gammln(nu+0.5_dp))/exp(gammln(nu))*sqrt(PI)
       p1(i)=exp(gammln(nu+0.5_dp))/exp(gammln(nu))*2.0_dp*&
            sqrt(PI_D)*a/(1.0_dp+k**2*a**2)**(nu+0.5_dp)
!       print *, i, p1(i)
    enddo
  END SUBROUTINE vk1df

!=====
! 2-D autocorrelation function
!    Gaussian, Exponential, von Karman
!    considers anisotropic case 
!              by substituting r/a --> sqrt((x/a)^2+(z/b)^2) 
!=======
  SUBROUTINE gauss2d(xnr,a,b,dx,dz,theta)
    IMPLICIT NONE 
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: xnr
    REAL(DP), INTENT(IN) :: dx,dz,a,b,theta
    REAL(DP) :: twora,x,z
    INTEGER :: i,j,nx,nz,nnx,nnz
    ! number of elements of xnr must be odd number
    nx=size(xnr,1)
    nz=size(xnr,2)
    nnx=(nx+1)/2
    nnz=(nz+1)/2
    print *, nx,nnx
    do i=nnx,nx
       do j=1,nz
          x=float(i-nnx)*dx
          z=float(j-nnz)*dz
          twora=((x*cos(theta)-z*sin(theta))/a)**2&
               +((x*sin(theta)+z*cos(theta))/b)**2
          xnr(i,j)=exp(-twora)
       enddo
    enddo
    xnr(1:nnx-1,1:nnz)=xnr(nx:nnx+1:-1,nz:nnz:-1)
    xnr(1:nnx-1,nnz+1:nz)=xnr(nx:nnx+1:-1,nnz-1:1:-1)
  END SUBROUTINE gauss2d

  SUBROUTINE gauss2df(p2,a,b,dkx,dkz,theta)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: p2
    REAL(DP), INTENT(IN) :: a,b,dkx,dkz,theta
    INTEGER :: i,j,nkx,nkz
    REAL(DP) :: ka
    nkx=size(p2,1)
    nkz=size(p2,2)
    do i=1,nkx
       do j=1,nkz/2+1
          ka=((float(i-1)*dkx*cos(theta)-float(j-1)*dkz*sin(theta))*a)**2&
               +((float(i-1)*dkx*sin(theta)+float(j-1)*dkz*cos(theta))*b)**2
          p2(i,j)=a**2/2.0_dp*exp(-ka/4.0_dp)
       enddo
       do j=nkz/2+2,nkz
          ka=((float(i-1)*dkx*cos(theta)-float(nkz-j+1)*dkz*sin(theta))*a)**2&
             +((float(i-1)*dkx*sin(theta)+float(nkz-j+1)*dkz*cos(theta))*b)**2
          p2(i,j)=a**2/2.0_dp*exp(-ka/4.0_dp)
       enddo
    enddo
  END SUBROUTINE gauss2df



  SUBROUTINE exp2d(xnr,a,b,dx,dz,theta)
    IMPLICIT NONE 
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: xnr
    REAL(DP), INTENT(IN) :: dx,dz,a,b,theta
    REAL(DP) :: twora,x,z
    INTEGER :: i,j,nx,nz,nnx,nnz
    nx=size(xnr,1)
    nz=size(xnr,2)
    nnx=(nx+1)/2
    nnz=(nz+1)/2
    print *, nx,nnx
    do i=nnx,nx
       do j=1,nz
          x=float(i-nnx)*dx
          z=float(j-nnz)*dz
          twora=((x*cos(theta)-z*sin(theta))/a)**2&
               +((x*sin(theta)+z*cos(theta))/b)**2
          xnr(i,j)=exp(-sqrt(twora))
       enddo
    enddo
    xnr(1:nnx-1,1:nnz)=xnr(nx:nnx+1:-1,nz:nnz:-1)
    xnr(1:nnx-1,nnz+1:nz)=xnr(nx:nnx+1:-1,nnz-1:1:-1)
  END SUBROUTINE exp2d


  SUBROUTINE exp2df(p2,a,b,dkx,dkz,theta)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: p2
    REAL(DP), INTENT(IN) :: a,b,dkx,dkz,theta
    INTEGER :: i,j,nkx,nkz
    REAL(DP) :: ka
    nkx=size(p2,1)
    nkz=size(p2,2)
    do i=1,nkx
       do j=1,nkz/2+1
          ka=((float(i-1)*dkx*cos(theta)-float(j-1)*dkz*sin(theta))*a)**2&
               +((float(i-1)*dkx*sin(theta)+float(j-1)*dkz*cos(theta))*b)**2
          p2(i,j)=a**2/(1.0_dp+ka)**1.5_dp
       enddo
       do j=nkz/2+2,nkz
          ka=((float(i-1)*dkx*cos(theta)-float(nkz-j+1)*dkz*sin(theta))*a)**2&
             +((float(i-1)*dkx*sin(theta)+float(nkz-j+1)*dkz*cos(theta))*b)**2
          p2(i,j)=a**2/(1.0_dp+ka)**1.5_dp
       enddo
    enddo
  END SUBROUTINE exp2df


  SUBROUTINE vk2d(xnr,a,b,nu,dx,dz,theta)
    USE nr, ONLY : bessik, gammln
    IMPLICIT NONE 
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: xnr
    REAL(DP), INTENT(IN) :: dx,dz,a,b,nu,theta
    REAL(DP) :: ra,k_nu,dummy,x,z
    INTEGER :: i,j,nx,nz,nnx,nnz
    nx=size(xnr,1)
    nz=size(xnr,2)
    nnx=(nx+1)/2
    nnz=(nz+1)/2
    print *, nx,nnx
    
    ! at x=0,z=0
    xnr(nnx,nnz)=1.0_dp
    ! this increase computing time, but I believe it may reduce bugs.
    do i=nnx,nx
       do j=1,nz
          if (.not.(j.eq.nnx .and. i.eq.nnz)) then
             x=float(i-nnx)*dx
             z=float(j-nnz)*dz
             ra=sqrt(((x*cos(theta)-z*sin(theta))/a)**2&
                  +((x*sin(theta)+z*cos(theta))/b)**2)
             call bessik(ra,nu,dummy,k_nu,dummy,dummy)
             xnr(i,j)=2.0_dp**(1.0-nu)/exp(gammln(nu))*(ra)**nu*k_nu
          endif
       enddo
    enddo
    xnr(1:nnx-1,1:nnz)=xnr(nx:nnx+1:-1,nz:nnz:-1)
    xnr(1:nnx-1,nnz+1:nz)=xnr(nx:nnx+1:-1,nnz-1:1:-1)


  END SUBROUTINE vk2d

  SUBROUTINE vk2df(p2,a,b,nu,dkx,dkz,theta)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: p2
    REAL(DP), INTENT(IN) :: a,b,dkx,dkz,theta,nu
    INTEGER :: i,j,nkx,nkz
    REAL(DP) :: ka
    nkx=size(p2,1)
    nkz=size(p2,2)
    do i=1,nkx
       do j=1,nkz/2+1
          ka=((float(i-1)*dkx*cos(theta)-float(j-1)*dkz*sin(theta))*a)**2&
               +((float(i-1)*dkx*sin(theta)+float(j-1)*dkz*cos(theta))*b)**2
          p2(i,j)=2.0_dp*nu*a**2/(1.0_dp+ka)**(1.0_dp+nu)
       enddo
       do j=nkz/2+2,nkz
          ka=((float(i-1)*dkx*cos(theta)-float(nkz-j+1)*dkz*sin(theta))*a)**2&
             +((float(i-1)*dkx*sin(theta)+float(nkz-j+1)*dkz*cos(theta))*b)**2
          p2(i,j)=2.0_dp*nu*a**2/(1.0_dp+ka)**(1.0_dp+nu)
       enddo
    enddo
  END SUBROUTINE vk2df

!====
! Function for cumulative probability density function
!   Bi-Peak Gaussian
!===
  FUNCTION cpd_std_gauss(x)
    USE nr, ONLY : erf
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: cpd_std_gauss
    cpd_std_gauss=0.5_dp*(1.0_dp+sign(erf(abs(x)/sqrt(2.0_dp)),x))
  END FUNCTION cpd_std_gauss

  FUNCTION cpd_bi_gauss(x,mu1,std1,mu2,std2,ratio1)
    REAL(DP), INTENT(IN) :: x, mu1, std1,mu2, std2,ratio1
    REAL(DP) :: cpd_bi_gauss
    cpd_bi_gauss=ratio1*cpd_std_gauss((x-mu1)/std1)+&
         (1.0_dp-ratio1)*cpd_std_gauss((x-mu2)/std2)
  END FUNCTION cpd_bi_gauss

 
  SUBROUTINE hist(data,histd,nbin,inlim1,inlim2,outindex_data,outcum_data)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN) :: data
    INTEGER, INTENT(IN) :: nbin
    INTEGER,DIMENSION(:), INTENT(OUT):: histd
    REAL(DP), INTENT(IN), OPTIONAL, TARGET :: inlim1,inlim2
    REAL(DP), DIMENSION(:), INTENT(OUT),OPTIONAL,TARGET :: outcum_data
    INTEGER, DIMENSION(:), INTENT(OUT),OPTIONAL, TARGET :: outindex_data
    ! define
    REAL(DP), POINTER :: lim1, lim2
    REAL(DP), POINTER, DIMENSION(:) :: cum_data
    INTEGER, POINTER, DIMENSION(:) :: index_data
    REAL(DP) :: dbin
    INTEGER :: i,n,index

    n=size(data,1)
    ! set lim
    if (present(inlim1)) then
       lim1=>inlim1
       lim2=>inlim2
    else
       allocate(lim1,lim2)
       lim1=minval(data)
       lim2=maxval(data)
    end if
    ! set index_data
    if (present(outindex_data)) then
       index_data=> outindex_data
    else
       allocate(index_data(n))
    end if
    dbin=(lim2-lim1)/nbin
    histd=0.0
    do i=1,n
       index=nint((data(i)-lim1)/dbin+0.5)
       index_data(i)=index
       histd(index)=histd(index)+1
    enddo
    if (present(outcum_data)) then
       cum_data=>outcum_data
       cum_data(1)=histd(1)
       do i=2,nbin
          cum_data(i)=cum_data(i-1)+real(histd(i))
       enddo
       cum_data=cum_data/float(n)
    end if
    !    print *, 'check6'

  END SUBROUTINE hist

!=========
! Evaluating Autocorrelation function of given data
!    1-D
!    2-D
!=========

  SUBROUTINE acorr1d(data,acorr,lagmaxin)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(IN):: data
    REAL(DP), DIMENSION(:),INTENT(OUT) :: acorr
    INTEGER, OPTIONAL, TARGET :: lagmaxin
    !
    ! This subroutine calculates the 1-D autocorrelation function
    ! by a direct method. 
    ! data                       : N real-valued data points
    ! acorr(1:lagmax-1)          : -lagmax =< lag < 0
    ! acorr(lagmax)              :  lag = 0
    ! acorr(lagmax+1:2*lagmax-1) :  0 < lag <= lagmax
    ! lagmax                     : max of lags
    !     NOTE: When lagmax is not present, lagmax is set to N-1.
    !
    INTEGER :: i,j,n
    INTEGER, POINTER :: lagmax
    REAL(DP), DIMENSION(:), ALLOCATABLE :: cdata
    n=size(data)
    if (present(lagmaxin)) then
       lagmax=>lagmaxin
    else
       allocate(lagmax)
       lagmax=n-1
    endif
    print *, "Maximum of lags is set to", lagmax
    allocate(cdata(n+lagmax))
    cdata(:)=0.d0
    cdata(1:n)=data(1:n)
    !Direct Method
    do i=0,lagmax,1
       do j=1,n
          acorr(i+lagmax+1)=acorr(i+lagmax+1)+cdata(j)*cdata(j+i)
       enddo
    enddo
    acorr=acorr/acorr(lagmax+1)
    acorr(1:lagmax)=acorr(2*lagmax+1:lagmax:-1)
    if (present(lagmaxin)) then
    else
       deallocate(lagmax)
    endif
    !  deallocate(cdata)
  END SUBROUTINE acorr1d
  
!  SUBROUTINE acorr1df(data,acorr,lagmaxin)
!    USE nr, ONLY : realft
!    USE nrtype
!
!    IMPLICIT NONE
!    REAL(DP), DIMENSION(:), INTENT(IN)  :: data
!    REAL(DP), DIMENSION(:), INTENT(INOUT) :: acorr
!    INTEGER, OPTIONAL, TARGET :: lagmaxin
!    !
!    ! This subroutine calculates the 1-D autocorrelation function
!    ! by a FFT method. 
!    ! N must be power of 2
!    ! data                       : N real-valued data points
!    ! acorr(1:lagmax-1)          : -lagmax =< lag < 0
!    ! acorr(lagmax)              :  lag = 0
!    ! acorr(lagmax+1:2*lagmax-1) :  0 < lag <= lagmax
!    ! lagmax                     : N-1
!    REAL(DP), DIMENSION(size(data)*2) :: cdata ! zero patch too large?
!    COMPLEX(DPC), DIMENSION(size(data)) :: cfdata
!    INTEGER, POINTER :: lagmax
!    INTEGER :: n, no2,nw,i
!    n=size(data)
!    print *, n
!   if (present(lagmaxin)) then
!       lagmax=>lagmaxin
!    else
!       allocate(lagmax)
!       lagmax=n-1
!    endif!
!
!    nw=int(n*0.05)
!    print *, "nw", nw
!    cdata=0.d0
!    cdata(1:n)=data
!!    cdata(1)=0.0_dp
!!    cdata(n)=0.0_dp
!!    do i=1,nw
!!       cdata(1+i)=sin(i/nw*PIO2_D)*cdata(1+i)
!!       cdata(n-i)=sin(i/nw*PIO2_D)*cdata(n-i)
!!    enddo
!
!    no2=n
!
!    call realft(cdata,1,cfdata)
!    print *, real(cfdata(1))
!    cfdata(1)=cmplx(real(cfdata(1))**2/no2,&
!         aimag(cfdata(1))**2/no2,kind=dpc)
!    cfdata(2:)=cfdata(2:)*conjg(cfdata(2:))/no2
!    call realft(cdata,-1,cfdata)
!     acorr=cdata(1:lagmax*2+1)/cdata(1)
!  END SUBROUTINE acorr1df

!================
! ACORR2D
!   calculate 2-D auto correlation function
!=============
  SUBROUTINE acorr2d(data,acorr,lag1maxin,lag2maxin)
    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(IN):: data
    REAL, DIMENSION(:,:),INTENT(INOUT) :: acorr
    INTEGER, OPTIONAL, TARGET :: lag1maxin, lag2maxin
    !
    ! This subroutine calculates the 1-D autocorrelation function
    ! by a direct method. 
    ! data                       : N real-valued data points
    ! acorr(1:lagmax-1)          : -lagmax =< lag < 0
    ! acorr(lagmax)              :  lag = 0
    ! acorr(lagmax+1:2*lagmax-1) :  0 < lag <= lagmax
    ! lagmax                     : max of lags
    !     NOTE: When lagmax is not present, lagmax is set to N-1.
    !
    INTEGER :: i,j,n1,n2,ilag,jlag
    INTEGER, POINTER :: lag1max,lag2max
    REAL, DIMENSION(:,:), ALLOCATABLE :: cdata
    REAL :: acorr0
    !
    n1=size(data,1)
    n2=size(data,2)
    print *, "data size ", n1, n2
    print *, "acorr size", size(acorr,1), size(acorr,2)
    if (present(lag1maxin)) then
       lag1max=>lag1maxin
       if (present(lag2maxin)) then
          lag2max=>lag2maxin
       else
          allocate(lag2max); lag2max=n2-1
       endif
    else
       allocate(lag1max); lag1max=n1-1
       allocate(lag2max); lag2max=n2-1
    endif
    print *, "Maximum of lags is set to", lag1max,lag2max
    allocate(cdata(n1+lag1max,-lag2max:n2+lag2max))

    !==
    ! START calculating autocorrelation function
    cdata=0.d0
    cdata(1:n1,1:n2)=data(1:n1,1:n2)
    acorr=0.d0
    !Direct Method
    do i=1,n1
       do j=1,n2
          do jlag=-lag2max,lag2max
             do ilag=0,lag1max
                acorr(lag1max+ilag+1,lag2max+jlag+1)=&
                     acorr(lag1max+1+ilag,lag2max+1+jlag)&
                     +cdata(i,j)*cdata(i+ilag,j+jlag)
             enddo
          enddo
       enddo
    enddo
    acorr0=acorr(lag1max+1,lag2max+1)
    acorr=acorr/acorr0
    ! Using symmetry at Lag0
    print *, "ehehe "
    acorr(1:lag1max,1:lag2max+1)=&
         acorr(lag1max*2+1:lag1max+2:-1,lag2max*2+1:lag2max+1:-1)
    acorr(1:lag1max,lag2max+2:lag2max*2+1)=&
         acorr(lag1max*2+1:lag1max+2:-1, lag2max:1:-1)
    print *, "degege"
    print *, size(acorr(1:lag1max,1:lag2max+1),1)
    print *, size(acorr(1:lag1max,1:lag2max+1),2)
    print *, size(acorr(lag1max*2+1:lag1max+2:-1,lag2max*2+1:lag2max+1:-1),1)
    print *, size(acorr(lag1max*2+1:lag1max+2:-1,lag2max*2+1:lag2max+1:-1),2)
    print *, size(acorr(1:lag1max,lag2max+1:lag2max*2+1),1)
    print *, size(acorr(1:lag1max,lag2max+1:lag2max*2+1),2)
    print *, size(acorr(lag1max*2+1:lag1max+2:-1, lag2max+1:1:-1),1)
    print *, size(acorr(lag1max*2+1:lag1max+2:-1, lag2max+1:1:-1),2)

    ! END calculating autocorrelation function
    !==
    print *, "ehehehe"

    if (.not.present(lag1maxin)) deallocate(lag1max)
    if (.not.present(lag2maxin)) deallocate(lag2max)
    deallocate(cdata)

  END SUBROUTINE acorr2d

END MODULE acflib
