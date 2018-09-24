FUNCTION chebev_s(a,b,c,x)
  USE nrtype; USE nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,b,x
  REAL(SP), DIMENSION(:), INTENT(IN) :: c
  REAL(SP) :: chebev_s
  ! Chebyshev evaluation: All arguments are input. c is an array of length M 
  ! of Chebyshev coefficients, the first M elements of c output from chebft 
  ! (which must have been called with the same a and b). The Chebyshev 
  ! polynomial Sum_k=1_M ckT_(k-1)(y)-c1/2 is evaluated at a point 
  ! y = [x-(b+a)/2]/[(b-a)/2], and the result is returned as the function value.
  INTEGER(I4B) :: j,m
  REAL(SP) :: d,dd,sv,y,y2
  if ((x-a)*(x-b) > 0.0) call nrerror('x not in range in chebev_s')
  m=size(c)
  d=0.0
  dd=0.0
  y=(2.0_sp*x-a-b)/(b-a)                    !Change of variable.
  y2=2.0_sp*y
  do j=m,2,-1                               !Clenshaw’s recurrence.
     sv=d
     d=y2*d-dd+c(j)
     dd=sv
  end do
  chebev_s=y*d-dd+0.5_sp*c(1)               ! Last step is different.
END FUNCTION chebev_s
FUNCTION chebev_v(a,b,c,x)
  USE nrtype; USE nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,b
  REAL(SP), DIMENSION(:), INTENT(IN) :: c,x
  REAL(SP), DIMENSION(size(x)) :: chebev_v
  INTEGER(I4B) :: j,m
  REAL(SP), DIMENSION(size(x)) :: d,dd,sv,y,y2
  if (any((x-a)*(x-b) > 0.0)) call nrerror('x not in range in chebev_v')
  m=size(c)
  d=0.0
  dd=0.0
  y=(2.0_sp*x-a-b)/(b-a)
  y2=2.0_sp*y
  do j=m,2,-1
     sv=d
     d=y2*d-dd+c(j)
     dd=sv
  end do
  chebev_v=y*d-dd+0.5_sp*c(1)
END FUNCTION chebev_v
