! real FFT subroutines
! from Numerical Recipe
!MODULE fft
!CONTAINS
SUBROUTINE fourrow_sp(data,isign)
  USE nrtype; USE nrutil, ONLY : assert,swap
  IMPLICIT NONE
  COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
  INTEGER(I4B), INTENT(IN) :: isign
  !Replaces each row (constant first index) of data(1:M,1:N) 
  !by its discrete Fourier transform (transform on second index), 
  !if isign is input as 1; or replaces each row of data by N times 
  !its inverse discrete Fourier transform, if isign is input 
  !as .1. N must be an integer power of 2. Parallelism is M-fold 
  !on the first index of data.
  INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
  REAL(DP) :: theta
  COMPLEX(SPC), DIMENSION(size(data,1)) :: temp
  COMPLEX(DPC) :: w,wp !Double precision for the trigonometric recurrences.
  COMPLEX(SPC) :: ws
  n=size(data,2)
  call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_sp')
  n2=n/2
  j=n2
  !This is the bit-reversal section of the routine.
  do i=1,n-2
     if (j > i) call swap(data(:,j+1),data(:,i+1))
     m=n2
     do
        if (m < 2 .or. j < m) exit
        j=j-m
        m=m/2
     end do
     j=j+m
  end do
  mmax=1
  !Here begins the Danielson-Lanczos section of the routine.
  do !Outer loop executed log2 N times.
     if (n <= mmax) exit
     istep=2*mmax
     theta=PI_D/(isign*mmax) !Initialize for the trigonometric recurrence.
     wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
     w=cmplx(1.0_dp,0.0_dp,kind=dpc)
     do m=1,mmax !Here are the two nested inner loops.
        ws=w
        do i=m,n,istep
           j=i+mmax
           temp=ws*data(:,j) !This is the Danielson-Lanczos formula.
           data(:,j)=data(:,i)-temp
           data(:,i)=data(:,i)+temp
        end do
        w=w*wp+w !Trigonometric recurrence.
     end do
     mmax=istep
  end do
END SUBROUTINE fourrow_sp

SUBROUTINE four1_sp(data,isign)
  USE nrtype; USE nrutil, ONLY : arth,assert
  USE nr, ONLY : fourrow
  IMPLICIT NONE
  COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
  INTEGER(I4B), INTENT(IN) :: isign
  !Replaces a complex array data by its discrete Fourier transform, 
  !if isign is input as 1; or replaces data by its inverse discrete 
  !Fourier transform times the size of data, if isign is input as -ss1. 
  !The size of data must be an integer power of 2. Parallelism is achieved
  !by internally reshaping the input array to two dimensions. (Use 
  !this version if fourrow is faster than fourcol on your machine.)
  COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
  COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
  REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
  INTEGER(I4B) :: n,m1,m2,j
  n=size(data)
  call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_sp')
  !Find dimensions as close to square as possible, allocate space, 
  !and reshape the input array.
  m1=2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
  m2=n/m1
  allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
  dat=reshape(data,shape(dat))
  call fourrow(dat,isign) !Transform on second index.
  theta=arth(0,isign,m1)*TWOPI_D/n !Set up recurrence.
  wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
  w=cmplx(1.0_dp,0.0_dp,kind=dpc)
  do j=2,m2 !Multiply by the extra phase factor.
     w=w*wp+w
     dat(:,j)=dat(:,j)*w
  end do
  temp=transpose(dat) !Transpose, and transform on (original) first index.
  call fourrow(temp,isign)
  data=reshape(temp,shape(data)) !Reshape the result back to one dimension.
  deallocate(dat,w,wp,theta,temp)
END SUBROUTINE four1_sp

SUBROUTINE four2(data,isign)
  USE nrtype
  USE nr, ONLY : fourrow
  IMPLICIT NONE
  COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
  INTEGER(I4B), INTENT(IN) :: isign
  !Replaces a 2-d complex array data by its discrete 2-d Fourier 
  !transform, if isign is input as 1; or replaces data by its inverse 
  !2-d discrete Fourier transform times the product of its two sizes, 
  !if isign is input as .1. Both of data¡Çs sizes must be integer powers 
  !of 2(this is checked for in fourrow). Parallelism is by use of fourrow.
  COMPLEX(SPC), DIMENSION(size(data,2),size(data,1)) :: temp
  call fourrow(data,isign) !Transform in second dimension.
  temp=transpose(data) !Tranpose.
  call fourrow(temp,isign) !Transform in (original) first dimension.
  data=transpose(temp) !Transpose into data.
END SUBROUTINE four2


!========================================================================
! REALFT_SP
!     Real Fourier Transform, single prescion.
!=======================================================================
SUBROUTINE realft_sp(data,isign,zdata)
  USE nrtype; USE nrutil, ONLY : assert,assert_eq,zroots_unity
  USE nr, ONLY : four1
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
  INTEGER(I4B), INTENT(IN) :: isign
  COMPLEX(SPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
  !When isign = 1, calculates the Fourier transform of a set of 
  !N real-valued data points, input in the array data. 
  !If the optional argument zdata is not present, the data are replaced
  !by the positive frequency half of its complex Fourier transform. 
  !The real-valued first and last components of the complex transform 
  !are returned as elements data(1) and data(2), respectively. If the 
  !complex array zdata of length N/2 is present, data is unchanged and
  !the transform is returned in zdata. N must be a power of 2. 
  !If isign = -1, this routine calculates the inverse transform of a complex 
  !data array if it is the transform of real data. (Result in this case 
  !must be multiplied by 2/N.) The data can be supplied either in data,
  !with zdata absent, or in zdata.
  INTEGER(I4B) :: n,ndum,nh,nq
  COMPLEX(SPC), DIMENSION(size(data)/4) :: w
  COMPLEX(SPC), DIMENSION(size(data)/4-1) :: h1,h2
  COMPLEX(SPC), DIMENSION(:), POINTER :: cdata !Used for internal complex computations.
  COMPLEX(SPC) :: z
  REAL(SP) :: c1=0.5_sp,c2
  n=size(data)
  call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_sp')
  nh=n/2
  nq=n/4
  if (present(zdata)) then
     ndum=assert_eq(n/2,size(zdata),'realft_sp')
     cdata=>zdata ! Use zdata as cdata.
     if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
  else
     allocate(cdata(n/2)) !Have to allocate storage ourselves.
     cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
  end if
  if (isign == 1) then
     c2=-0.5_sp
     !print *, size(cdata)
     call four1(cdata,+1) !The forward transform is here.
  else !Otherwise set up for an inverse transform.
     c2=0.5_sp
  end if

  w=zroots_unity(sign(n,isign),n/4)
  w=cmplx(-aimag(w),real(w),kind=spc)
  h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1))) 
  !The two separate transforms are separated out of cdata. 
  h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
  !Next they are recombined to form the true transform of the original real data:
  cdata(2:nq)=h1+w(2:nq)*h2
  cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
  z=cdata(1) 
  !Squeeze the first and last data together to get them all within the original array.
  if (isign == 1) then
     cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=spc)
  else
     cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=spc)
     call four1(cdata,-1) !This is the inverse transform for the case isign=-1.
  end if
  if (present(zdata)) then !Ship out answer in data if required.
     if (isign /= 1) then
        data(1:n-1:2)=real(cdata)
        data(2:n:2)=aimag(cdata)
     end if
  else
     data(1:n-1:2)=real(cdata)
     data(2:n:2)=aimag(cdata)
     deallocate(cdata)
  end if
END SUBROUTINE realft_sp

!=================================================================
! RLFT2
!     Real 2D Fourier Transform, single prescion.
!================================================================ 
SUBROUTINE rlft2(data,spec,speq,isign)
  USE nrtype; USE nrutil, ONLY : assert,assert_eq
  USE nr, ONLY : four2
  REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: data
  COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: spec
  COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: speq
  INTEGER(I4B), INTENT(IN) :: isign
  !Given a two-dimensional real array data(1:M,1:N), this routine 
  !returns (for isign=1) the complex fast Fourier transform as 
  !two complex arrays: On output, spec(1:M/2,1:N) contains the 
  !zero and positive frequency values of the first frequency component, 
  !while speq(1:N) contains the Nyquist critical frequency values of 
  !the first frequency component. The second frequency components are 
  !stored for zero, positive, and negative frequencies, in standard 
  !wrap-around order. For isign=-1, the inverse transform (times M * N/2 as
  !a constant multiplicative factor) is performed, with output data 
  !deriving from input spec and speq. For inverse transforms on data 
  !not generated first by a forward transform, make
  !sure the complex input data array satisfies property (12.5.2). The 
  !size of all arrays must always be integer powers of 2.
  INTEGER :: i1,j1,nn1,nn2
  REAL(DP) :: theta
  COMPLEX(SPC) :: c1=(0.5_sp,0.0_sp),c2,h1,h2,w
  COMPLEX(SPC), DIMENSION(size(data,2)-1) :: h1a,h2a
  COMPLEX(DPC) :: ww,wp
  nn1=assert_eq(size(data,1),2*size(spec,1),'rlft2: nn1')
  nn2=assert_eq(size(data,2),size(spec,2),size(speq),'rlft2: nn2')
  call assert(iand((/nn1,nn2/),(/nn1,nn2/)-1)==0, &
       'dimensions must be powers of 2 in rlft2')
  c2=cmplx(0.0_sp,-0.5_sp*isign,kind=spc)
  theta=TWOPI_D/(isign*nn1)
  wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=spc)
  if (isign == 1) then !Case of forward transform.
     spec(:,:)=cmplx(data(1:nn1:2,:),data(2:nn1:2,:),kind=spc)
     call four2(spec,isign) !Here is where most all of the compute time is spent. 
     speq=spec(1,:)
  end if
  h1=c1*(spec(1,1)+conjg(speq(1)))
  h1a=c1*(spec(1,2:nn2)+conjg(speq(nn2:2:-1)))
  h2=c2*(spec(1,1)-conjg(speq(1)))
  h2a=c2*(spec(1,2:nn2)-conjg(speq(nn2:2:-1)))
  spec(1,1)=h1+h2
  spec(1,2:nn2)=h1a+h2a
  speq(1)=conjg(h1-h2)
  speq(nn2:2:-1)=conjg(h1a-h2a)
  ww=cmplx(1.0_dp,0.0_dp,kind=dpc) !Initialize trigonometric recurrence.
  do i1=2,nn1/4+1
     j1=nn1/2-i1+2 !Corresponding negative frequency.
     ww=ww*wp+ww !Do the trig recurrence.
     w=ww
     h1=c1*(spec(i1,1)+conjg(spec(j1,1))) !Equation (12.3.5).
     h1a=c1*(spec(i1,2:nn2)+conjg(spec(j1,nn2:2:-1)))
     h2=c2*(spec(i1,1)-conjg(spec(j1,1)))
     h2a=c2*(spec(i1,2:nn2)-conjg(spec(j1,nn2:2:-1)))
     spec(i1,1)=h1+w*h2
     spec(i1,2:nn2)=h1a+w*h2a
     spec(j1,1)=conjg(h1-w*h2)
     spec(j1,nn2:2:-1)=conjg(h1a-w*h2a)
  end do
  if (isign == -1) then !Case of reverse transform.
     call four2(spec,isign)
     data(1:nn1:2,:)=real(spec)
     data(2:nn1:2,:)=aimag(spec)
  end if
END SUBROUTINE rlft2

