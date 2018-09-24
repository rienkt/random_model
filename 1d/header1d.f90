MODULE header1d
  USE nrtype
  USE header_stat
  !-----------------------------------------------------------
  ! TYPE     : ACF 1. Gauss, 2. Exp, 3. vK
  ! NU       : Hurst number
  ! N        : Size of random media
  ! A        : Autocorrelation length
  ! DZ       : dz
  ! DK       : dk
  !
  ! ISEED    : random seed for generating THETA
  ! NLAGMAX  : Max of lag for calculating ACF
  !-----------------------------------------------------------
  INTEGER  :: type,nlagmax,n,iseed,nreal
  REAL(DP) :: nu, a,dk,dz
  !
  !----------------------------------------------------------
  ! S0       : Desired power spectral density (PSD)
  ! SG0      : Input PSD of Gaussian
  ! SG       : Output PSD of Gaussian
  ! SB       : Output PSD of Bi-modal
  ! SGE      : ensamble of SG (useless in future)
  ! SBE      : ensamble of SB
  ! THETA    : Phase angle
  !----------------------------------------------------------
  REAL(DP), DIMENSION(:), ALLOCATABLE :: s0,sb,sg,sg0,theta
  REAL(DP), DIMENSION(:), ALLOCATABLE :: sge, sbe
  !
  !----------------------------------------------------------
  ! ACF0     : Desired ACF
  ! ACF      : ACF of generated media
  !----------------------------------------------------------
  REAL(DP), DIMENSION(:), ALLOCATABLE :: acf, acf0
  !



CONTAINS
  SUBROUTINE read_params
    IMPLICIT NONE
    CHARACTER(64) :: buf
    CHARACTER(64) :: val
    REAL(DP) :: lagmax

    write(*,*) "Reading Parameters ......"
   
    do while( (buf .ne. "END").or.(buf .ne. "end"))
       read(5,*,END=100) buf,val
       write(*, '(2A20)') buf,val
       
       if ((buf .eq. "n")     .or.(buf .eq. "N"))     read(val,*) n
       if ((buf .eq. "dz")    .or.(buf .eq. "DZ"))    read(val,*) dz
       if ((buf .eq. "type")  .or.(buf .eq. "TYPE"))  read(val,*) type
       if ((buf .eq. "nu")    .or.(buf .eq. "NU"))    read(val,*) nu
       if ((buf .eq. "a")     .or.(buf .eq. "A"))     read(val,*) a
       if ((buf .eq. "iseed") .or.(buf .eq. "ISEED")) read(val,*) iseed
       if ((buf .eq. "lagmax").or.(buf .eq. "LAGMAX"))read(val,*) lagmax
       if ((buf .eq. "mu1")   .or.(buf .eq. "MU1"))   read(val,*) mu1
       if ((buf .eq. "mu2")   .or.(buf .eq. "MU2"))   read(val,*) mu2
       if ((buf .eq. "std1")  .or.(buf .eq. "STD1"))  read(val,*) std1
       if ((buf .eq. "std2")  .or.(buf .eq. "STD2"))  read(val,*) std2
       if ((buf .eq. "ratio1").or.(buf .eq. "RATIO1"))read(val,*) ratio1
       if ((buf .eq. "nreal").or.(buf .eq. "NREAL"))read(val,*) nreal
    enddo
100 continue
    
    dk=1.0_dp/dz/float(n)*TWOPI ! Nyquist : 1/2dz >> kan

    if (type == 3) then
    elseif (type==2) then
       nu=0.5_dp
    elseif (type == 1) then
       nu=100.0_dp
    else
       print *, "Wrong type"
       stop
    end if
    nlagmax=int(lagmax/dz)

    write(*,*) "===PARAMETERS===="
    write(*,*) "* random field"
    write(*,1000) nreal
    write(*,1001) n, dz, real(n)*dz
    write(*,1002) dk, PI_D/dz
    write(*,1003) type, nu, a
    write(*,1004) iseed
    write(*,*) "* stat value of (bi-modal) field"
    write(*,1005) mu1, std1
    write(*,1006) mu2, std2
    write(*,1007) ratio1
    write(*,*) "* check"
    write(*,1008) lagmax,nlagmax
    write(*,*)
1000 format('  relization ', i10, 'times')
1001 format('  n =', i10,   ' dz  =', f8.4, ' size=',f8.4)
1002 format('  dk=', f10.4, ' nqst=', f8.4)
1003 format('  type=',i8,   ' nu  =', f8.4, ' a   =',f8.4)
1004 format('  seed=',i10)
1005 format('  unit 1 : mu=',f9.4,' std=',f9.4)
1006 format('  unit 2 : mu=',f9.4,' std=',f9.4)
1007 format('  ratio of unit 1 : 'f4.2)
1008 format('  lag max :',f8.4 , ' (', i10,')')

  END SUBROUTINE read_params

END MODULE header1d



MODULE header_stat
  USE nrtype
  REAL(DP) :: mu1,mu2,std1,std2,ratio1
END MODULE header_stat


MODULE modfft4gsub
  USE nrtype
  !
  !----------------------------------------------------------
  ! for fft4g.f90
  !  w  : cos, sin table for n-data length
  !  ip : work area for bit reversal for n-data length
  !  w2 : for n*2-data length
  !  ip2: 
  !----------------------------------------------------------
  INTEGER,DIMENSION(:),ALLOCATABLE :: ip
  REAL(DP), DIMENSION(:),ALLOCATABLE :: w
  INTEGER,DIMENSION(:),ALLOCATABLE :: ip2
  REAL(DP), DIMENSION(:),ALLOCATABLE :: w2
  !

END MODULE modfft4gsub

 
