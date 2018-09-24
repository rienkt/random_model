SUBROUTINE moment_sp(data,ave,adev,sdev,var,skew,curt)
  USE nrtype; USE nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL(SP), INTENT(OUT) :: ave,adev,sdev,var,skew,curt
  REAL(SP), DIMENSION(:), INTENT(IN) :: data
  !Given an array of data, this routine returns its mean ave, 
  !average deviation adev, standard deviation sdev, variance var, 
  !skewness skew, and kurtosis curt.
  INTEGER(I4B) :: n
  REAL(SP) :: ep
  REAL(SP), DIMENSION(size(data)) :: p,s
  n=size(data)
  if (n <= 1) call nrerror('moment: n must be at least 2')
  ave=sum(data(:))/n !First pass to get the mean.
  s(:)=data(:)-ave   !Second pass to get the first (absolute), 
                     !second, third, and fourth moments of the deviation 
                     !from the mean. 
  ep=sum(s(:))
  adev=sum(abs(s(:)))/n
  p(:)=s(:)*s(:)
  var=sum(p(:))
  p(:)=p(:)*s(:)
  skew=sum(p(:))
  p(:)=p(:)*s(:)
  curt=sum(p(:))
  var=(var-ep**2/n)/(n-1) !Corrected two-pass formula.
  sdev=sqrt(var)
  if (var /= 0.0) then
     skew=skew/(n*sdev**3)
     curt=curt/(n*var**2)-3.0_sp
  else
     call nrerror('moment: no skew or kurtosis when zero variance')
  end if
END SUBROUTINE moment_sp

SUBROUTINE moment_dp(data,ave,adev,sdev,var,skew,curt)
  USE nrtype; USE nrutil, ONLY : nrerror
  IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: ave,adev,sdev,var,skew,curt
  REAL(DP), DIMENSION(:), INTENT(IN) :: data
  !Given an array of data, this routine returns its mean ave, 
  !average deviation adev, standard deviation sdev, variance var, 
  !skewness skew, and kurtosis curt.
  INTEGER(I4B) :: n
  REAL(DP) :: ep
  REAL(DP), DIMENSION(size(data)) :: p,s
  n=size(data)
  if (n <= 1) call nrerror('moment: n must be at least 2')
  ave=sum(data(:))/n !First pass to get the mean.
  s(:)=data(:)-ave   !Second pass to get the first (absolute), 
                     !second, third, and fourth moments of the deviation 
                     !from the mean. 
  ep=sum(s(:))
  adev=sum(abs(s(:)))/n
  p(:)=s(:)*s(:)
  var=sum(p(:))
  p(:)=p(:)*s(:)
  skew=sum(p(:))
  p(:)=p(:)*s(:)
  curt=sum(p(:))
  var=(var-ep**2/n)/(n-1) !Corrected two-pass formula.
  sdev=sqrt(var)
  if (var /= 0.0) then
     skew=skew/(n*sdev**3)
     curt=curt/(n*var**2)-3.0_sp
  else
     call nrerror('moment: no skew or kurtosis when zero variance')
  end if
END SUBROUTINE moment_dp

