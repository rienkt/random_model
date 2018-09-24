SUBROUTINE ran3_s_sp(harvest)
  USE nrtype
  USE ran_state, ONLY: K4B,amm,lenran,ran_init,ran_hash,mran0,nran0,rans
  IMPLICIT NONE
  REAL(SP), INTENT(OUT) :: harvest
  !Random number generation by DES-like hashing of two 32-bit words, 
  !using the algorithm ran hash. Returns as harvest a uniform random 
  !deviate between 0.0 and 1.0 (exclusive of the endpoint values).
  INTEGER(K4B) :: temp
!  print *, "ran3_s"
  if (lenran < 1) call ran_init(1)  !Initialize.
  nran0=ieor(nran0,ishft(nran0,13)) !Two Marsaglia shift sequences are
                                    !maintained as input to the hashing.
                                    !The period of the combined
                                    !generator is about 1.8 x 1019.
  nran0=ieor(nran0,ishft(nran0,-17))
  nran0=ieor(nran0,ishft(nran0,5))
  if (nran0 == 1) nran0=270369_k4b
  rans=nran0
  mran0=ieor(mran0,ishft(mran0,5))
  mran0=ieor(mran0,ishft(mran0,-13))
  mran0=ieor(mran0,ishft(mran0,6))
  temp=mran0
  call ran_hash(temp,rans) !Hash.
  harvest=amm*merge(rans,not(rans), rans<0 ) 
  !Make the result positive definite (note that amm is negative). 
END SUBROUTINE ran3_s_sp

SUBROUTINE ran3_v_sp(harvest)
  USE nrtype
  USE ran_state, ONLY: K4B,amm,lenran,ran_init,ran_hash,mran,nran,ranv
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
  INTEGER(K4B), DIMENSION(size(harvest)) :: temp
  INTEGER(K4B) :: n
  n=size(harvest)
!  print *, "ran3_v", lenran, n+1
  if (lenran < n+1) call ran_init(n+1)
  nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
  nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
  nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
  where (nran(1:n) == 1) nran(1:n)=270369_k4b
  ranv(1:n)=nran(1:n)
  mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
  mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
  mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
  temp=mran(1:n)
  call ran_hash(temp,ranv(1:n))
  harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
END SUBROUTINE ran3_v_sp


SUBROUTINE ran3_s_dp(harvest)
  USE nrtype
  USE ran_state, ONLY: K4B,amm,lenran,ran_init,ran_hash,mran0,nran0,rans
  IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: harvest
  !Random number generation by DES-like hashing of two 32-bit words, 
  !using the algorithm ran hash. Returns as harvest a uniform random 
  !deviate between 0.0 and 1.0 (exclusive of the endpoint values).
  INTEGER(K4B) :: temp
!  print *, "ran3_s"
  if (lenran < 1) call ran_init(1)  !Initialize.
  nran0=ieor(nran0,ishft(nran0,13)) !Two Marsaglia shift sequences are
                                    !maintained as input to the hashing.
                                    !The period of the combined
                                    !generator is about 1.8 x 1019.
  nran0=ieor(nran0,ishft(nran0,-17))
  nran0=ieor(nran0,ishft(nran0,5))
  if (nran0 == 1) nran0=270369_k4b
  rans=nran0
  mran0=ieor(mran0,ishft(mran0,5))
  mran0=ieor(mran0,ishft(mran0,-13))
  mran0=ieor(mran0,ishft(mran0,6))
  temp=mran0
  call ran_hash(temp,rans) !Hash.
  harvest=amm*merge(rans,not(rans), rans<0 ) 
  !Make the result positive definite (note that amm is negative). 
END SUBROUTINE ran3_s_dp

SUBROUTINE ran3_v_dp(harvest)
  USE nrtype
  USE ran_state, ONLY: K4B,amm,lenran,ran_init,ran_hash,mran,nran,ranv
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(OUT) :: harvest
  INTEGER(K4B), DIMENSION(size(harvest)) :: temp
  INTEGER(K4B) :: n
  n=size(harvest)
!  print *, "ran3_v", lenran, n+1
  if (lenran < n+1) call ran_init(n+1)
  nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
  nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
  nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
  where (nran(1:n) == 1) nran(1:n)=270369_k4b
  ranv(1:n)=nran(1:n)
  mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
  mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
  mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
  temp=mran(1:n)
  call ran_hash(temp,ranv(1:n))
  harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
END SUBROUTINE ran3_v_dp
