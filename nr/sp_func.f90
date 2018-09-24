FUNCTION gammln_s_sp(xx)
  USE nrtype; USE nrutil, ONLY : arth,assert
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: xx
  REAL(SP) :: gammln_s_sp
  ! Returns the value ln[Gamma(xx)] for xx > 0.
  REAL(DP) :: tmp,x
  ! Internal arithmetic will be done in double precision, a nicety 
  ! that you can omit if five-figure accuracy is good enough.
  REAL(DP) :: stp = 2.5066282746310005_dp
  REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
       -86.50532032941677_dp,24.01409824083091_dp,&
       -1.231739572450155_dp,0.1208650973866179e-2_dp,&
       -0.5395239384953e-5_dp/)
  call assert(xx > 0.0, 'gammln_s arg')
  x=xx
  tmp=x+5.5_dp
  tmp=(x+0.5_dp)*log(tmp)-tmp
  gammln_s_sp=tmp+dlog(stp*(1.000000000190015_dp+&
       sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
END FUNCTION gammln_s_sp

FUNCTION gammln_v_sp(xx)
  USE nrtype; USE nrutil, ONLY: assert
  IMPLICIT NONE
  INTEGER(I4B) :: i
  REAL(SP), DIMENSION(:), INTENT(IN) :: xx
  REAL(SP), DIMENSION(size(xx)) :: gammln_v_sp
  REAL(DP), DIMENSION(size(xx)) :: ser,tmp,x,y
  REAL(DP) :: stp = 2.5066282746310005_dp
  REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
       -86.50532032941677_dp,24.01409824083091_dp,&
       -1.231739572450155_dp,0.1208650973866179e-2_dp,&
       -0.5395239384953e-5_dp/)
  if (size(xx) == 0) RETURN
  call assert(all(xx > 0.0), 'gammln_v arg')
  x=xx
  tmp=x+5.5_dp
  tmp=(x+0.5_dp)*log(tmp)-tmp
  ser=1.000000000190015_dp
  y=x
  do i=1,size(coef)
     y=y+1.0_dp
     ser=ser+coef(i)/y
  end do
  gammln_v_sp=tmp+log(stp*ser/x)
END FUNCTION gammln_v_sp

FUNCTION gammln_s_dp(xx)
  USE nrtype; USE nrutil, ONLY : arth,assert
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: xx
  REAL(DP) :: gammln_s_dp
  ! Returns the value ln[Gamma(xx)] for xx > 0.
  REAL(DP) :: tmp,x
  ! Internal arithmetic will be done in double precision, a nicety 
  ! that you can omit if five-figure accuracy is good enough.
  REAL(DP) :: stp = 2.5066282746310005_dp
  REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
       -86.50532032941677_dp,24.01409824083091_dp,&
       -1.231739572450155_dp,0.1208650973866179e-2_dp,&
       -0.5395239384953e-5_dp/)
  call assert(xx > 0.0, 'gammln_s arg')
  x=xx
  tmp=x+5.5_dp
  tmp=(x+0.5_dp)*dlog(tmp)-tmp
  gammln_s_dp=tmp+dlog(stp*(1.000000000190015_dp+&
       sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
END FUNCTION gammln_s_dp

FUNCTION gammln_v_dp(xx)
  USE nrtype; USE nrutil, ONLY: assert
  IMPLICIT NONE
  INTEGER(I4B) :: i
  REAL(DP), DIMENSION(:), INTENT(IN) :: xx
  REAL(DP), DIMENSION(size(xx)) :: gammln_v_dp
  REAL(DP), DIMENSION(size(xx)) :: ser,tmp,x,y
  REAL(DP) :: stp = 2.5066282746310005_dp
  REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
       -86.50532032941677_dp,24.01409824083091_dp,&
       -1.231739572450155_dp,0.1208650973866179e-2_dp,&
       -0.5395239384953e-5_dp/)
  if (size(xx) == 0) RETURN
  call assert(all(xx > 0.0), 'gammln_v arg')
  x=xx
  tmp=x+5.5_dp
  tmp=(x+0.5_dp)*log(tmp)-tmp
  ser=1.000000000190015_dp
  y=x
  do i=1,size(coef)
     y=y+1.0_dp
     ser=ser+coef(i)/y
  end do
  gammln_v_dp=tmp+log(stp*ser/x)
END FUNCTION gammln_v_dp


!========

FUNCTION gcf_s_sp(a,x,gln)
  USE nrtype; USE nrutil, ONLY : nrerror
  USE nr, ONLY : gammln
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,x
  REAL(SP), OPTIONAL, INTENT(OUT) :: gln
  REAL(SP) :: gcf_s_sp
  INTEGER(I4B), PARAMETER :: ITMAX=100
  REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
  ! Returns the incomplete gamma function Q(a, x) evaluated by its 
  ! continued fraction representation as gammcf. Also optionally returns 
  ! ln Gamma(a) as gln.
  ! Parameters: ITMAX is the maximum allowed number of iterations; 
  ! EPS is the relative accuracy;
  ! FPMIN is a number near the smallest representable floating-point number.
  INTEGER(I4B) :: i
  REAL(SP) :: an,b,c,d,del,h
  if (x == 0.0) then
     gcf_s_sp=1.0
     RETURN
  end if
  b=x+1.0_sp-a ! Set up for evaluating continued fraction by modified
               ! Lentz¡Çs method (¡ø5.2) with b0 = 0. 
  c=1.0_sp/FPMIN
  d=1.0_sp/b
  h=d
  do i=1,ITMAX ! Iterate to convergence.
     an=-i*(i-a)
     b=b+2.0_sp
     d=an*d+b
     if (abs(d) < FPMIN) d=FPMIN
     c=b+an/c
     if (abs(c) < FPMIN) c=FPMIN
     d=1.0_sp/d
     del=d*c
     h=h*del
     if (abs(del-1.0_sp) <= EPS) exit
  end do
  if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_s')
  if (present(gln)) then
     gln=gammln(a)
     gcf_s_sp=exp(-x+a*log(x)-gln)*h !Put factors in front.
  else
     gcf_s_sp=exp(-x+a*log(x)-gammln(a))*h
  end if
END FUNCTION gcf_s_sp

FUNCTION gcf_v_sp(a,x,gln)
  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  USE nr, ONLY : gammln
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
  REAL(SP), DIMENSION(size(a)) :: gcf_v_sp
  INTEGER(I4B), PARAMETER :: ITMAX=100
  REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
  INTEGER(I4B) :: i
  REAL(SP), DIMENSION(size(a)) :: an,b,c,d,del,h
  LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
  i=assert_eq(size(a),size(x),'gcf_v')
  zero=(x == 0.0)
  where (zero)
     gcf_v_sp=1.0
  elsewhere
     b=x+1.0_sp-a
     c=1.0_sp/FPMIN
     d=1.0_sp/b
     h=d
  end where
  converged=zero
  do i=1,ITMAX
     where (.not. converged)
        an=-i*(i-a)
        b=b+2.0_sp
        d=an*d+b
        d=merge(FPMIN,d, abs(d)<FPMIN )
        c=b+an/c
        c=merge(FPMIN,c, abs(c)<FPMIN )
        d=1.0_sp/d
        del=d*c
        h=h*del
        converged = (abs(del-1.0_sp)<=EPS)
     end where
     if (all(converged)) exit
  end do
  if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_v')
  if (present(gln)) then
     if (size(gln) < size(a)) call &
          nrerror('gser: Not enough space for gln')
     gln=gammln(a)
     where (.not. zero) gcf_v_sp=exp(-x+a*log(x)-gln)*h
  else
     where (.not. zero) gcf_v_sp=exp(-x+a*log(x)-gammln(a))*h
  end if
END FUNCTION gcf_v_sp


FUNCTION gcf_s_dp(a,x,gln)
  USE nrtype; USE nrutil, ONLY : nrerror
  USE nr, ONLY : gammln
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a,x
  REAL(DP), OPTIONAL, INTENT(OUT) :: gln
  REAL(DP) :: gcf_s_dp
  INTEGER(I4B), PARAMETER :: ITMAX=100
  REAL(DP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
  ! Returns the incomplete gamma function Q(a, x) evaluated by its 
  ! continued fraction representation as gammcf. Also optionally returns 
  ! ln Gamma(a) as gln.
  ! Parameters: ITMAX is the maximum allowed number of iterations; 
  ! EPS is the relative accuracy;
  ! FPMIN is a number near the smallest representable floating-point number.
  INTEGER(I4B) :: i
  REAL(DP) :: an,b,c,d,del,h
  if (x == 0.0) then
     gcf_s_dp=1.0
     RETURN
  end if
  b=x+1.0_dp-a ! Set up for evaluating continued fraction by modified
               ! Lentz¡Çs method (¡ø5.2) with b0 = 0. 
  c=1.0_dp/FPMIN
  d=1.0_dp/b
  h=d
  do i=1,ITMAX ! Iterate to convergence.
     an=-i*(i-a)
     b=b+2.0_dp
     d=an*d+b
     if (abs(d) < FPMIN) d=FPMIN
     c=b+an/c
     if (abs(c) < FPMIN) c=FPMIN
     d=1.0_dp/d
     del=d*c
     h=h*del
     if (abs(del-1.0_dp) <= EPS) exit
  end do
  if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_s')
  if (present(gln)) then
     gln=gammln(a)
     gcf_s_dp=exp(-x+a*log(x)-gln)*h !Put factors in front.
  else
     gcf_s_dp=exp(-x+a*log(x)-gammln(a))*h
  end if
END FUNCTION gcf_s_dp

FUNCTION gcf_v_dp(a,x,gln)
  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  USE nr, ONLY : gammln
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
  REAL(DP), DIMENSION(size(a)) :: gcf_v_dp
  INTEGER(I4B), PARAMETER :: ITMAX=100
  REAL(DP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
  INTEGER(I4B) :: i
  REAL(DP), DIMENSION(size(a)) :: an,b,c,d,del,h
  LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
  i=assert_eq(size(a),size(x),'gcf_v')
  zero=(x == 0.0)
  where (zero)
     gcf_v_dp=1.0
  elsewhere
     b=x+1.0_sp-a
     c=1.0_dp/FPMIN
     d=1.0_dp/b
     h=d
  end where
  converged=zero
  do i=1,ITMAX
     where (.not. converged)
        an=-i*(i-a)
        b=b+2.0_dp
        d=an*d+b
        d=merge(FPMIN,d, abs(d)<FPMIN )
        c=b+an/c
        c=merge(FPMIN,c, abs(c)<FPMIN )
        d=1.0_dp/d
        del=d*c
        h=h*del
        converged = (abs(del-1.0_dp)<=EPS)
     end where
     if (all(converged)) exit
  end do
  if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_v')
  if (present(gln)) then
     if (size(gln) < size(a)) call &
          nrerror('gser: Not enough space for gln')
     gln=gammln(a)
     where (.not. zero) gcf_v_dp=exp(-x+a*log(x)-gln)*h
  else
     where (.not. zero) gcf_v_dp=exp(-x+a*log(x)-gammln(a))*h
  end if
END FUNCTION gcf_v_dp


!====

FUNCTION gser_s_sp(a,x,gln)
  USE nrtype; USE nrutil, ONLY : nrerror
  USE nr, ONLY : gammln
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,x
  REAL(SP), OPTIONAL, INTENT(OUT) :: gln
  REAL(SP) :: gser_s_sp
  INTEGER(I4B), PARAMETER :: ITMAX=100
  REAL(SP), PARAMETER :: EPS=epsilon(x)
  ! Returns the incomplete gamma function P(a, x) evaluated by 
  ! its series representation as gamser. Also optionally returns 
  ! ln Gamma (a) as gln.
  INTEGER(I4B) :: n
  REAL(SP) :: ap,del,summ
  if (x == 0.0) then
     gser_s_sp=0.0
     RETURN
  end if
  ap=a
  summ=1.0_sp/a
  del=summ
  do n=1,ITMAX
     ap=ap+1.0_sp
     del=del*x/ap
     summ=summ+del
     if (abs(del) < abs(summ)*EPS) exit
  end do
  if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_s')
  if (present(gln)) then
     gln=gammln(a)
     gser_s_sp=summ*exp(-x+a*log(x)-gln)
  else
     gser_s_sp=summ*exp(-x+a*log(x)-gammln(a))
  end if
END FUNCTION gser_s_sp

FUNCTION gser_v_sp(a,x,gln)
  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  USE nr, ONLY : gammln
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
  REAL(SP), DIMENSION(size(a)) :: gser_v_sp
  INTEGER(I4B), PARAMETER :: ITMAX=100
  REAL(SP), PARAMETER :: EPS=epsilon(x)
  INTEGER(I4B) :: n
  REAL(SP), DIMENSION(size(a)) :: ap,del,summ
  LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
  n=assert_eq(size(a),size(x),'gser_v')
  zero=(x == 0.0)
  where (zero) gser_v_sp=0.0
  ap=a
  summ=1.0_sp/a
  del=summ
  converged=zero
  do n=1,ITMAX
     where (.not. converged)
        ap=ap+1.0_sp
        del=del*x/ap
        summ=summ+del
        converged = (abs(del) < abs(summ)*EPS)
     end where
     if (all(converged)) exit
  end do
  if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_v')
  if (present(gln)) then
     if (size(gln) < size(a)) call &
          nrerror('gser: Not enough space for gln')
     gln=gammln(a)
     where (.not. zero) gser_v_sp=summ*exp(-x+a*log(x)-gln)
  else
     where (.not. zero) gser_v_sp=summ*exp(-x+a*log(x)-gammln(a))
  end if
END FUNCTION gser_v_sp

FUNCTION gser_s_dp(a,x,gln)
  USE nrtype; USE nrutil, ONLY : nrerror
  USE nr, ONLY : gammln
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,x
  REAL(SP), OPTIONAL, INTENT(OUT) :: gln
  REAL(SP) :: gser_s_dp
  INTEGER(I4B), PARAMETER :: ITMAX=100
  REAL(SP), PARAMETER :: EPS=epsilon(x)
  ! Returns the incomplete gamma function P(a, x) evaluated by 
  ! its series representation as gamser. Also optionally returns 
  ! ln Gamma (a) as gln.
  INTEGER(I4B) :: n
  REAL(SP) :: ap,del,summ
  if (x == 0.0) then
     gser_s_dp=0.0
     RETURN
  end if
  ap=a
  summ=1.0_dp/a
  del=summ
  do n=1,ITMAX
     ap=ap+1.0_dp
     del=del*x/ap
     summ=summ+del
     if (abs(del) < abs(summ)*EPS) exit
  end do
  if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_s')
  if (present(gln)) then
     gln=gammln(a)
     gser_s_dp=summ*exp(-x+a*log(x)-gln)
  else
     gser_s_dp=summ*exp(-x+a*log(x)-gammln(a))
  end if
END FUNCTION gser_s_dp

FUNCTION gser_v_dp(a,x,gln)
  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  USE nr, ONLY : gammln
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
  REAL(DP), DIMENSION(size(a)) :: gser_v_dp
  INTEGER(I4B), PARAMETER :: ITMAX=100
  REAL(DP), PARAMETER :: EPS=epsilon(x)
  INTEGER(I4B) :: n
  REAL(DP), DIMENSION(size(a)) :: ap,del,summ
  LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
  n=assert_eq(size(a),size(x),'gser_v')
  zero=(x == 0.0)
  where (zero) gser_v_dp=0.0
  ap=a
  summ=1.0_dp/a
  del=summ
  converged=zero
  do n=1,ITMAX
     where (.not. converged)
        ap=ap+1.0_dp
        del=del*x/ap
        summ=summ+del
        converged = (abs(del) < abs(summ)*EPS)
     end where
     if (all(converged)) exit
  end do
  if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_v')
  if (present(gln)) then
     if (size(gln) < size(a)) call &
          nrerror('gser: Not enough space for gln')
     gln=gammln(a)
     where (.not. zero) gser_v_dp=summ*exp(-x+a*log(x)-gln)
  else
     where (.not. zero) gser_v_dp=summ*exp(-x+a*log(x)-gammln(a))
  end if
END FUNCTION gser_v_dp


FUNCTION gammp_s_sp(a,x)
  USE nrtype; USE nrutil, ONLY : assert
  USE nr, ONLY : gcf,gser
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,x
  REAL(SP) :: gammp_s_sp
  !Returns the incomplete gamma function P(a, x).
  call assert( x >= 0.0, a > 0.0, 'gammp_s args')
  if (x<a+1.0_sp) then !Use the series representation.
     gammp_s_sp=gser(a,x)
  else                       !Use the continued fraction representation
     gammp_s_sp=1.0_sp-gcf(a,x) !and take its complement.
  end if
END FUNCTION gammp_s_sp

FUNCTION gammp_v_sp(a,x)
  USE nrtype; USE nrutil, ONLY : assert,assert_eq
  USE nr, ONLY : gcf,gser
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(SP), DIMENSION(size(x)) :: gammp_v_sp
  LOGICAL(LGT), DIMENSION(size(x)) :: mask
  INTEGER(I4B) :: ndum
  ndum=assert_eq(size(a),size(x),'gammp_v')
  call assert( all(x >= 0.0), all(a > 0.0), 'gammp_v args')
  mask = (x<a+1.0_sp)
  gammp_v_sp=merge(gser(a,merge(x,0.0_sp,mask)), &
       1.0_sp-gcf(a,merge(x,0.0_sp,.not. mask)),mask)
END FUNCTION gammp_v_sp

FUNCTION gammp_s_dp(a,x)
  USE nrtype; USE nrutil, ONLY : assert
  USE nr, ONLY : gcf,gser
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a,x
  REAL(DP) :: gammp_s_dp
  !Returns the incomplete gamma function P(a, x).
  call assert( x >= 0.0, a > 0.0, 'gammp_s args')
  if (x<a+1.0_dp) then !Use the series representation.
     gammp_s_dp=gser(a,x)
  else                       !Use the continued fraction representation
     gammp_s_dp=1.0_dp-gcf(a,x) !and take its complement.
  end if
END FUNCTION gammp_s_dp

FUNCTION gammp_v_dp(a ,x)
  USE nrtype; USE nrutil, ONLY : assert,assert_eq
  USE nr, ONLY : gcf,gser
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
  REAL(DP), DIMENSION(size(x)) :: gammp_v_dp
  LOGICAL(LGT), DIMENSION(size(x)) :: mask
  INTEGER(I4B) :: ndum
  ndum=assert_eq(size(a),size(x),'gammp_v')
  call assert( all(x >= 0.0), all(a > 0.0), 'gammp_v args')
  mask = (x<a+1.0_dp)
  gammp_v_dp=merge(gser(a,merge(x,0.0_dp,mask)), &
       1.0_dp-gcf(a,merge(x,0.0_dp,.not. mask)),mask)
END FUNCTION gammp_v_dp



FUNCTION erf_s_sp(x)
  USE nrtype
  USE nr, ONLY : gammp
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: x
  REAL(SP) :: erf_s_sp
  !Returns the error function erf(x).
  erf_s_sp=gammp(0.5_sp,x**2)
  if (x < 0.0) erf_s_sp=-erf_s_sp
END FUNCTION erf_s_sp

FUNCTION erf_v_sp(x)
  USE nrtype
  USE nr, ONLY : gammp
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: x
  REAL(SP), DIMENSION(size(x)) :: erf_v_sp
  erf_v_sp=gammp(spread(0.5_sp,1,size(x)),x**2)
  where (x < 0.0) erf_v_sp=-erf_v_sp
END FUNCTION erf_v_sp


FUNCTION erf_s_dp(x)
  USE nrtype
  USE nr, ONLY : gammp
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x
  REAL(DP) :: erf_s_dp
  !Returns the error function erf(x).
  erf_s_dp=gammp(0.5_dp,x**2)
  if (x < 0.0) erf_s_dp=-erf_s_dp
END FUNCTION erf_s_dp

FUNCTION erf_v_dp(x)
  USE nrtype
  USE nr, ONLY : gammp
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: x
  REAL(DP), DIMENSION(size(x)) :: erf_v_dp
  erf_v_dp=gammp(spread(0.5_dp,1,size(x)),x**2)
  where (x < 0.0) erf_v_dp=-erf_v_dp
END FUNCTION erf_v_dp


SUBROUTINE bessik_sp(x,xnu,ri,rk,rip,rkp)
  USE nrtype; USE nrutil, ONLY : assert,nrerror
  USE nr, ONLY : beschb
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: x,xnu
  REAL(SP), INTENT(OUT) :: ri,rk,rip,rkp
  INTEGER(I4B), PARAMETER :: MAXIT=10000
  REAL(SP), PARAMETER :: XMIN=2.0
  REAL(DP), PARAMETER :: EPS=1.0e-10_dp,FPMIN=1.0e-30_dp
  ! Returns the modified Bessel functions ri = I_nu, rk = K_nu and their 
  ! derivatives rip = I_nu',rkp = K_nu' , for positive x and for xnu = ge 0. 
  ! The relative accuracy is within one ortwo significant digits of EPS. 
  ! FPMIN is a number close to the machine¡Çs smallest floatingpoint
  ! number. All internal arithmetic is in double precision. To convert the 
  ! entire routineto double precision, change the REAL declaration above and 
  ! decrease EPS to 10.16. Also convert the subroutine beschb.
  INTEGER(I4B) :: i,l,nl
  REAL(DP) :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,&
       gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,&
       ril,ril1,rimu,rip1,ripl,ritemp,rk1,rkmu,rkmup,rktemp,&
       s,sum,sum1,x2,xi,xi2,xmu,xmu2
  call assert(x > 0.0, xnu >= 0.0, 'bessik args')
  nl=int(xnu+0.5_dp)              ! nl is the number of downward recurrences
                                  ! of the I's and upward recurrences
                                  ! of K's. xmu lies between .1/2 and 1/2.
  xmu=xnu-nl
  xmu2=xmu*xmu
  xi=1.0_dp/x
  xi2=2.0_dp*xi
  h=xnu*xi                        ! Evaluate CF1 by modified Lentz's method
                                  !(sec.5.2). 
  if (h < FPMIN) h=FPMIN
  b=xi2*xnu
  d=0.0
  c=h
  do i=1,MAXIT
     b=b+xi2
     d=1.0_dp/(b+d)               ! Denominators cannot be zero here, so no
                                  ! need for special precautions. 
     c=b+1.0_dp/c
     del=c*d
     h=del*h
     if (abs(del-1.0_dp) < EPS) exit
  end do
  if (i > MAXIT) call nrerror('x too large in bessik; try asymptotic expansion')
  ril=FPMIN                ! Initialize I_nu and I_nu for downward recurrence.
  ripl=h*ril
  ril1=ril                                 !Store values for later rescaling.
  rip1=ripl
  fact=xnu*xi
  do l=nl,1,-1
     ritemp=fact*ril+ripl
     fact=fact-xi
     ripl=fact*ritemp+ril
     ril=ritemp
  end do
  f=ripl/ril                     ! Now have unnormalized I_nu and I_nu.
  if (x < XMIN) then               ! Use series.
     x2=0.5_dp*x
     pimu=PI_D*xmu
     if (abs(pimu) < EPS) then
        fact=1.0
     else
        fact=pimu/sin(pimu)
     end if
     d=-log(x2)
     e=xmu*d
     if (abs(e) < EPS) then
        fact2=1.0
     else
        fact2=sinh(e)/e
     end if
     call beschb(xmu,gam1,gam2,gampl,gammi) ! Chebyshev evaluation of Gamma1 
                                            ! and Gamma2.
     ff=fact*(gam1*cosh(e)+gam2*fact2*d)    ! f0.
     sum=ff
     e=exp(e)
     p=0.5_dp*e/gampl                       ! p0.
     q=0.5_dp/(e*gammi)                     ! q0.
     c=1.0
     d=x2*x2
     sum1=p
     do i=1,MAXIT
        ff=(i*ff+p+q)/(i*i-xmu2)
        c=c*d/i
        p=p/(i-xmu)
        q=q/(i+xmu)
        del=c*ff
        sum=sum+del
        del1=c*(p-i*ff)
        sum1=sum1+del1
        if (abs(del) < abs(sum)*EPS) exit
     end do
     if (i > MAXIT) call nrerror('bessk series failed to converge')
     rkmu=sum
     rk1=sum1*xi2
  else                            ! Evaluate CF2 by Steed's algorithm (sec5.2),
                                  !which is OK because there can be no
                                  !zero denominators.
     b=2.0_dp*(1.0_dp+x)
     d=1.0_dp/b
     delh=d
     h=delh
     q1=0.0                       !Initializations for recurrence (6.7.35).
     q2=1.0
     a1=0.25_dp-xmu2
     c=a1
     q=c                          !First term in equation (6.7.34).
     a=-a1
     s=1.0_dp+q*delh
     do i=2,MAXIT
        a=a-2*(i-1)
        c=-a*c/i
        qnew=(q1-b*q2)/a
        q1=q2
        q2=qnew
        q=q+c*qnew
        b=b+2.0_dp
        d=1.0_dp/(b+a*d)
        delh=(b*d-1.0_dp)*delh
        h=h+delh
        dels=q*delh
        s=s+dels
        if (abs(dels/s) < EPS) exit !Need only test convergence of sum, since
                                    !CF2 itself converges more quickly. 
     end do
     if (i > MAXIT) call nrerror('bessik: failure to converge in cf2')
     h=a1*h
     rkmu=sqrt(PI_D/(2.0_dp*x))*exp(-x)/s !Omit the factor exp(-x) to scale 
                                          !all the returned functions by exp(x) 
                                          !for x le XMIN.
     rk1=rkmu*(xmu+x+0.5_dp-h)*xi
  end if
  rkmup=xmu*xi*rkmu-rk1
  rimu=xi/(f*rkmu-rkmup)                  ! Get I_nu from Wronskian.
  ri=(rimu*ril1)/ril                      ! Scale original I_nu and I_nu' .
  rip=(rimu*rip1)/ril
  do i=1,nl                               !Upward recurrence of K_nu.
     rktemp=(xmu+i)*xi2*rk1+rkmu
     rkmu=rk1
     rk1=rktemp
  end do
  rk=rkmu
  rkp=xnu*xi*rkmu-rk1
END SUBROUTINE bessik_sp

SUBROUTINE beschb_s(x,gam1,gam2,gampl,gammi)
  USE nrtype
  USE nr, ONLY : chebev
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x
  REAL(DP), INTENT(OUT) :: gam1,gam2,gampl,gammi
  INTEGER(I4B), PARAMETER :: NUSE1=5,NUSE2=5
  !Evaluates Gamma1 and Gamma2 by Chebyshev expansion for |x| < 1/2. 
  !Also returns 1/Gamma(1 + x) and 1/Gamma(1 - x). If converting to double 
  !precision, set NUSE1 = 7, NUSE2 = 8.
  REAL(SP) :: xx
  REAL(SP), DIMENSION(7) :: c1=(/-1.142022680371168_sp,&
       6.5165112670737e-3_sp,3.087090173086e-4_sp,-3.4706269649e-6_sp,&
       6.9437664e-9_sp,3.67795e-11_sp,-1.356e-13_sp/)
  REAL(SP), DIMENSION(8) :: c2=(/1.843740587300905_sp,&
       -7.68528408447867e-2_sp,1.2719271366546e-3_sp,&
       -4.9717367042e-6_sp, -3.31261198e-8_sp,2.423096e-10_sp,&
       -1.702e-13_sp,-1.49e-15_sp/)
  xx=8.0_dp*x*x-1.0_dp !Multiply x by 2 to make range be .1 to 1, and then apply
                       !transformation for evaluating even Chebyshev series.
  gam1=chebev(-1.0_sp,1.0_sp,c1(1:NUSE1),xx)
  gam2=chebev(-1.0_sp,1.0_sp,c2(1:NUSE2),xx)
  gampl=gam2-x*gam1
  gammi=gam2+x*gam1
END SUBROUTINE beschb_s

SUBROUTINE beschb_v(x,gam1,gam2,gampl,gammi)
  USE nrtype
  USE nr, ONLY : chebev
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: x
  REAL(DP), DIMENSION(:), INTENT(OUT) :: gam1,gam2,gampl,gammi
  INTEGER(I4B), PARAMETER :: NUSE1=5,NUSE2=5
  REAL(SP), DIMENSION(size(x)) :: xx
  REAL(SP), DIMENSION(7) :: c1=(/-1.142022680371168_sp,&
       6.5165112670737e-3_sp,3.087090173086e-4_sp,-3.4706269649e-6_sp,&
       6.9437664e-9_sp,3.67795e-11_sp,-1.356e-13_sp/)
  REAL(SP), DIMENSION(8) :: c2=(/1.843740587300905_sp,&
       -7.68528408447867e-2_sp,1.2719271366546e-3_sp,&
       -4.9717367042e-6_sp, -3.31261198e-8_sp,2.423096e-10_sp,&
       -1.702e-13_sp,-1.49e-15_sp/)
  xx=8.0_dp*x*x-1.0_dp
  gam1=chebev(-1.0_sp,1.0_sp,c1(1:NUSE1),xx)
  gam2=chebev(-1.0_sp,1.0_sp,c2(1:NUSE2),xx)
  gampl=gam2-x*gam1
  gammi=gam2+x*gam1
END SUBROUTINE beschb_v

SUBROUTINE bessik_dp(x,xnu,ri,rk,rip,rkp)
  USE nrtype; USE nrutil, ONLY : assert,nrerror
  USE nr, ONLY : beschb
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x,xnu
  REAL(DP), INTENT(OUT) :: ri,rk,rip,rkp
  INTEGER(I4B), PARAMETER :: MAXIT=10000
  REAL(DP), PARAMETER :: XMIN=2.0
  REAL(DP), PARAMETER :: EPS=1.0e-10_dp,FPMIN=1.0e-30_dp
  ! Returns the modified Bessel functions ri = I_nu, rk = K_nu and their 
  ! derivatives rip = I_nu',rkp = K_nu' , for positive x and for xnu = ge 0. 
  ! The relative accuracy is within one ortwo significant digits of EPS. 
  ! FPMIN is a number close to the machine¡Çs smallest floatingpoint
  ! number. All internal arithmetic is in double precision. To convert the 
  ! entire routineto double precision, change the REAL declaration above and 
  ! decrease EPS to 10.16. Also convert the subroutine beschb.
  INTEGER(I4B) :: i,l,nl
  REAL(DP) :: a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,&
       gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,&
       ril,ril1,rimu,rip1,ripl,ritemp,rk1,rkmu,rkmup,rktemp,&
       s,sum,sum1,x2,xi,xi2,xmu,xmu2
  call assert(x > 0.0, xnu >= 0.0, 'bessik args')
  nl=int(xnu+0.5_dp)              ! nl is the number of downward recurrences
                                  ! of the I's and upward recurrences
                                  ! of K's. xmu lies between .1/2 and 1/2.
  xmu=xnu-nl
  xmu2=xmu*xmu
  xi=1.0_dp/x
  xi2=2.0_dp*xi
  h=xnu*xi                        ! Evaluate CF1 by modified Lentz's method
                                  !(sec.5.2). 
  if (h < FPMIN) h=FPMIN
  b=xi2*xnu
  d=0.0
  c=h
  do i=1,MAXIT
     b=b+xi2
     d=1.0_dp/(b+d)               ! Denominators cannot be zero here, so no
                                  ! need for special precautions. 
     c=b+1.0_dp/c
     del=c*d
     h=del*h
     if (abs(del-1.0_dp) < EPS) exit
  end do
  if (i > MAXIT) call nrerror('x too large in bessik; try asymptotic expansion')
  ril=FPMIN                ! Initialize I_nu and I_nu for downward recurrence.
  ripl=h*ril
  ril1=ril                                 !Store values for later rescaling.
  rip1=ripl
  fact=xnu*xi
  do l=nl,1,-1
     ritemp=fact*ril+ripl
     fact=fact-xi
     ripl=fact*ritemp+ril
     ril=ritemp
  end do
  f=ripl/ril                     ! Now have unnormalized I_nu and I_nu.
  if (x < XMIN) then               ! Use series.
     x2=0.5_dp*x
     pimu=PI_D*xmu
     if (abs(pimu) < EPS) then
        fact=1.0
     else
        fact=pimu/sin(pimu)
     end if
     d=-log(x2)
     e=xmu*d
     if (abs(e) < EPS) then
        fact2=1.0
     else
        fact2=sinh(e)/e
     end if
     call beschb(xmu,gam1,gam2,gampl,gammi) ! Chebyshev evaluation of Gamma1 
                                            ! and Gamma2.
     ff=fact*(gam1*cosh(e)+gam2*fact2*d)    ! f0.
     sum=ff
     e=exp(e)
     p=0.5_dp*e/gampl                       ! p0.
     q=0.5_dp/(e*gammi)                     ! q0.
     c=1.0
     d=x2*x2
     sum1=p
     do i=1,MAXIT
        ff=(i*ff+p+q)/(i*i-xmu2)
        c=c*d/i
        p=p/(i-xmu)
        q=q/(i+xmu)
        del=c*ff
        sum=sum+del
        del1=c*(p-i*ff)
        sum1=sum1+del1
        if (abs(del) < abs(sum)*EPS) exit
     end do
     if (i > MAXIT) call nrerror('bessk series failed to converge')
     rkmu=sum
     rk1=sum1*xi2
  else                            ! Evaluate CF2 by Steed's algorithm (sec5.2),
                                  !which is OK because there can be no
                                  !zero denominators.
     b=2.0_dp*(1.0_dp+x)
     d=1.0_dp/b
     delh=d
     h=delh
     q1=0.0                       !Initializations for recurrence (6.7.35).
     q2=1.0
     a1=0.25_dp-xmu2
     c=a1
     q=c                          !First term in equation (6.7.34).
     a=-a1
     s=1.0_dp+q*delh
     do i=2,MAXIT
        a=a-2*(i-1)
        c=-a*c/i
        qnew=(q1-b*q2)/a
        q1=q2
        q2=qnew
        q=q+c*qnew
        b=b+2.0_dp
        d=1.0_dp/(b+a*d)
        delh=(b*d-1.0_dp)*delh
        h=h+delh
        dels=q*delh
        s=s+dels
        if (abs(dels/s) < EPS) exit !Need only test convergence of sum, since
                                    !CF2 itself converges more quickly. 
     end do
     if (i > MAXIT) call nrerror('bessik: failure to converge in cf2')
     h=a1*h
     rkmu=sqrt(PI_D/(2.0_dp*x))*exp(-x)/s !Omit the factor exp(-x) to scale 
                                          !all the returned functions by exp(x) 
                                          !for x le XMIN.
     rk1=rkmu*(xmu+x+0.5_dp-h)*xi
  end if
  rkmup=xmu*xi*rkmu-rk1
  rimu=xi/(f*rkmu-rkmup)                  ! Get I_nu from Wronskian.
  ri=(rimu*ril1)/ril                      ! Scale original I_nu and I_nu' .
  rip=(rimu*rip1)/ril
  do i=1,nl                               !Upward recurrence of K_nu.
     rktemp=(xmu+i)*xi2*rk1+rkmu
     rkmu=rk1
     rk1=rktemp
  end do
  rk=rkmu
  rkp=xnu*xi*rkmu-rk1
END SUBROUTINE bessik_dp

!=================================================
! Subroutine in the Chapter "Function Evaluation"
!==============================

FUNCTION chebev_s_sp(a,b,c,x)
  USE nrtype; USE nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,b,x
  REAL(SP), DIMENSION(:), INTENT(IN) :: c
  REAL(SP) :: chebev_s_sp
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
  do j=m,2,-1                               !Clenshaw¡Çs recurrence.
     sv=d
     d=y2*d-dd+c(j)
     dd=sv
  end do
  chebev_s_sp=y*d-dd+0.5_sp*c(1)               ! Last step is different.
END FUNCTION chebev_s_sp

FUNCTION chebev_v_sp(a,b,c,x)
  USE nrtype; USE nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: a,b
  REAL(SP), DIMENSION(:), INTENT(IN) :: c,x
  REAL(SP), DIMENSION(size(x)) :: chebev_v_sp
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
  chebev_v_sp=y*d-dd+0.5_sp*c(1)
END FUNCTION chebev_v_sp

FUNCTION chebev_s_dp(a,b,c,x)
  USE nrtype; USE nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a,b,x
  REAL(DP), DIMENSION(:), INTENT(IN) :: c
  REAL(DP) :: chebev_s_dp
  ! Chebyshev evaluation: All arguments are input. c is an array of length M 
  ! of Chebyshev coefficients, the first M elements of c output from chebft 
  ! (which must have been called with the same a and b). The Chebyshev 
  ! polynomial Sum_k=1_M ckT_(k-1)(y)-c1/2 is evaluated at a point 
  ! y = [x-(b+a)/2]/[(b-a)/2], and the result is returned as the function value.
  INTEGER(I4B) :: j,m
  REAL(DP) :: d,dd,sv,y,y2
  if ((x-a)*(x-b) > 0.0) call nrerror('x not in range in chebev_s')
  m=size(c)
  d=0.0
  dd=0.0
  y=(2.0_dp*x-a-b)/(b-a)                    !Change of variable.
  y2=2.0_dp*y
  do j=m,2,-1                               !Clenshaw¡Çs recurrence.
     sv=d
     d=y2*d-dd+c(j)
     dd=sv
  end do
  chebev_s_dp=y*d-dd+0.5_dp*c(1)               ! Last step is different.
END FUNCTION chebev_s_dp
FUNCTION chebev_v_dp(a,b,c,x)
  USE nrtype; USE nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: a,b
  REAL(DP), DIMENSION(:), INTENT(IN) :: c,x
  REAL(DP), DIMENSION(size(x)) :: chebev_v_dp
  INTEGER(I4B) :: j,m
  REAL(DP), DIMENSION(size(x)) :: d,dd,sv,y,y2
  if (any((x-a)*(x-b) > 0.0)) call nrerror('x not in range in chebev_v')
  m=size(c)
  d=0.0_dp
  dd=0.0_dp
  y=(2.0_dp*x-a-b)/(b-a)
  y2=2.0_dp*y
  do j=m,2,-1
     sv=d
     d=y2*d-dd+c(j)
     dd=sv
  end do
  chebev_v_dp=y*d-dd+0.5_dp*c(1)
END FUNCTION chebev_v_dp
