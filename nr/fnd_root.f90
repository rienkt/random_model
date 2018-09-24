FUNCTION rtbis_sp(func,x1,x2,xacc)
  USE nrtype; USE nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: x1,x2,xacc
  REAL(SP) :: rtbis_sp
  INTERFACE
     FUNCTION func(x)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: func
     END FUNCTION func
  END INTERFACE
  INTEGER(I4B), PARAMETER :: MAXIT=40
  ! Using bisection, find the root of a function func known to lie 
  ! between x1 and x2. The root, returned as rtbis, will be refined 
  ! until its accuracy is }xacc.
  ! Parameter: MAXIT is the maximum allowed number of bisections.
  INTEGER(I4B) :: j
  REAL(SP) :: dx,f,fmid,xmid
  fmid=func(x2)
  f=func(x1)
  if (f*fmid >= 0.0) call nrerror('rtbis: root must be bracketed')
  if (f < 0.0) then ! Orient the search so that f>0 lies at x+dx.
     rtbis_sp=x1
     dx=x2-x1
  else
     rtbis_sp=x2
     dx=x1-x2
  end if
  do j=1,MAXIT ! Bisection loop.
     dx=dx*0.5_sp
     xmid=rtbis_sp+dx
     fmid=func(xmid)
     if (fmid <= 0.0) rtbis_sp=xmid
     if (abs(dx) < xacc .or. fmid == 0.0) RETURN
  end do
  call nrerror('rtbis: too many bisections')
END FUNCTION rtbis_sp


FUNCTION rtbis_dp(func,x1,x2,xacc)
  USE nrtype; USE nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL(DP), INTENT(IN) :: x1,x2,xacc
  REAL(DP) :: rtbis_dp
  INTERFACE
     FUNCTION func(x)
       USE nrtype
       IMPLICIT NONE
       REAL(DP), INTENT(IN) :: x
       REAL(DP) :: func
     END FUNCTION func
  END INTERFACE
  INTEGER(I4B), PARAMETER :: MAXIT=40
  ! Using bisection, find the root of a function func known to lie 
  ! between x1 and x2. The root, returned as rtbis, will be refined 
  ! until its accuracy is }xacc.
  ! Parameter: MAXIT is the maximum allowed number of bisections.
  INTEGER(I4B) :: j
  REAL(DP) :: dx,f,fmid,xmid
  fmid=func(x2)
  f=func(x1)
  if (f*fmid >= 0.0) call nrerror('rtbis: root must be bracketed')
  if (f < 0.0) then ! Orient the search so that f>0 lies at x+dx.
     rtbis_dp=x1
     dx=x2-x1
  else
     rtbis_dp=x2
     dx=x1-x2
  end if
  do j=1,MAXIT ! Bisection loop.
     dx=dx*0.5_dp
     xmid=rtbis_dp+dx
     fmid=func(xmid)
     if (fmid <= 0.0) rtbis_dp=xmid
     if (abs(dx) < xacc .or. fmid == 0.0) RETURN
  end do
  call nrerror('rtbis: too many bisections')
END FUNCTION rtbis_dp
