MODULE ran_state
  ! This module supports the random number routines ran0, 
  ! ran1, ran2, and ran3. It provides each generator with five integers 
  ! (for vector versions, five vectors of integers), for use as internal 
  ! state space. The first three integers (iran, jran, kran) are maintained
  ! as nonnegative values, while the last two (mran, nran) have 32-bit 
  ! nonzero values. Also provided by this module is support for initializing 
  ! or reinitializing the state space to a desired standard sequence number, 
  ! hashing the initial values to random values, and allocating and
  ! deallocating the internal workspace.
  USE nrtype
  IMPLICIT NONE
  INTEGER, PARAMETER :: K4B=selected_int_kind(9)
  !Independent of the usual integer kind I4B, we need a kind value 
  !for (ideally) 32-bit integers.
  INTEGER(K4B), PARAMETER :: hg=huge(1_K4B), hgm=-hg, hgng=hgm-1
  INTEGER(K4B), SAVE :: lenran=0, seq=0
  INTEGER(K4B), SAVE :: iran0,jran0,kran0,nran0,mran0,rans
  INTEGER(K4B), DIMENSION(:,:), POINTER, SAVE :: ranseeds
  INTEGER(K4B), DIMENSION(:), POINTER, SAVE :: iran,jran,kran, &
       nran,mran,ranv
  REAL(SP), SAVE :: amm
  INTERFACE ran_hash !Scalar and vector versions of the hashing procedure.
     MODULE PROCEDURE ran_hash_s, ran_hash_v
  END INTERFACE
CONTAINS
!
!==
! RAN_INIT(LENGTH)
!==
  SUBROUTINE ran_init(length)
    USE nrtype; USE nrutil, ONLY : arth,nrerror,reallocate
    IMPLICIT NONE
    INTEGER(K4B), INTENT(IN) :: length
    ! Initialize or reinitialize the random generator state 
    ! space to vectors of size length. The saved variable seq is 
    ! hashed (via calls to the module routine ran hash) to create unique
    ! starting seeds, different for each vector component.
    INTEGER(K4B) :: new,j,hgt
    if (length < lenran) RETURN  ! Simply return if enough space is already 
                                 !allocated.
    hgt=hg
    ! The following lines check that kind value K4B is in fact 
    ! a 32-bit integer with the usual properties that we expect it to have 
    ! (under negation and wrap-around addition). If all of these tests are 
    ! satisfied, then the routines that use this module are portable, even 
    ! though they go beyond Fortran 90Åfs integer model.
    if (hg /= 2147483647) call nrerror('ran_init: arith assump 1 fails')
    if (hgng >= 0) call nrerror('ran_init: arith assump 2 fails')
    if (hgt+1 /= hgng) call nrerror('ran_init: arith assump 3 fails')
    if (not(hg) >= 0) call nrerror('ran_init: arith assump 4 fails')
    if (not(hgng) < 0) call nrerror('ran_init: arith assump 5 fails')
    if (hg+hgng >= 0) call nrerror('ran_init: arith assump 6 fails')
    if (not(-1_k4b) < 0) call nrerror('ran_init: arith assump 7 fails')
    if (not(0_k4b) >= 0) call nrerror('ran_init: arith assump 8 fails')
    if (not(1_k4b) >= 0) call nrerror('ran_init: arith assump 9 fails')
    if (lenran > 0) then ! Reallocate space, or ...
       ranseeds=>reallocate(ranseeds,length,5)
       ranv=>reallocate(ranv,length-1)
       new=lenran+1
    else                 ! allocate space.
       allocate(ranseeds(length,5))
       allocate(ranv(length-1))
       new=1                 !Index of first location not yet initialized.
       amm=nearest(1.0_sp,-1.0_sp)/hgng ! Use of nearest is to ensure that
                                        ! returned random deviates are 
                                        ! strictly less than 1.0.
       if (amm*hgng >= 1.0 .or. amm*hgng <= 0.0) &
            call nrerror('ran_init: arth assump 10 fails')
    end if
    !Set starting values, unique by seq and vector component.
    ranseeds(new:,1)=seq
    ranseeds(new:,2:5)=spread(arth(new,1,size(ranseeds(new:,1))),2,4)
    do j=1,4 ! Hash them.
       call ran_hash(ranseeds(new:,j),ranseeds(new:,j+1))
    end do
    where (ranseeds(new:,1:3) < 0) &              ! Enforce nonnegativity.
    ranseeds(new:,1:3)=not(ranseeds(new:,1:3))
    where (ranseeds(new:,4:5) == 0) ranseeds(new:,4:5)=1 ! Enforce nonzero.
    if (new == 1) then !Set scalar seeds.
       iran0=ranseeds(1,1)
       jran0=ranseeds(1,2)
       kran0=ranseeds(1,3)
       mran0=ranseeds(1,4)
       nran0=ranseeds(1,5)
       rans=nran0
    end if
    if (length > 1) then !Point to vector seeds.
       iran => ranseeds(2:,1)
       jran => ranseeds(2:,2)
       kran => ranseeds(2:,3)
       mran => ranseeds(2:,4)
       nran => ranseeds(2:,5)
       ranv = nran
    end if
    lenran=length
  END SUBROUTINE ran_init
!
!===
! RAN_DEALLOCATE
!==
  SUBROUTINE ran_deallocate
    !User interface to release the workspace used by the random number 
    !routines.
    if (lenran > 0) then
       deallocate(ranseeds,ranv)
       nullify(ranseeds,ranv,iran,jran,kran,mran,nran)
       lenran = 0
    end if
  END SUBROUTINE ran_deallocate
!===
! RAN_SEED
!===
  SUBROUTINE ran_seed(sequence,size,put,get)
    IMPLICIT NONE
    INTEGER, OPTIONAL, INTENT(IN) :: sequence
    INTEGER, OPTIONAL, INTENT(OUT) :: size
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN) :: put
    INTEGER, DIMENSION(:), OPTIONAL, INTENT(OUT) :: get
    !User interface for seeding the random number routines. Syntax is 
    !exactly like Fortran 90Åfs random seed routine, with one additional 
    !argument keyword: sequence, set to any integer value, causes an 
    !immediate new initialization, seeded by that integer.
    if (present(size)) then
       size=5*lenran
    else if (present(put)) then
       if (lenran == 0) RETURN
       ranseeds=reshape(put,shape(ranseeds))
       where (ranseeds(:,1:3) < 0) ranseeds(:,1:3)=not(ranseeds(:,1:3))
          !Enforce nonnegativity and nonzero conditions on any 
          !user-supplied seeds.
       where (ranseeds(:,4:5) == 0) ranseeds(:,4:5)=1
       iran0=ranseeds(1,1)
       jran0=ranseeds(1,2)
       kran0=ranseeds(1,3)
       mran0=ranseeds(1,4)
       nran0=ranseeds(1,5)
    else if (present(get)) then
       if (lenran == 0) RETURN
       ranseeds(1,1:5)=(/ iran0,jran0,kran0,mran0,nran0 /)
       get=reshape(ranseeds,shape(get))
    else if (present(sequence)) then
       call ran_deallocate
       seq=sequence
    end if
  END SUBROUTINE ran_seed
!===
! RAN_HASH_S(IL,IR)
!===
  SUBROUTINE ran_hash_s(il,ir)
    IMPLICIT NONE
    INTEGER(K4B), INTENT(INOUT) :: il,ir
    !DES-like hashing of two 32-bit integers, using shifts, 
    !xorÅfs, and adds to make the internal nonlinear function.
    INTEGER(K4B) :: is,j
    do j=1,4
       is=ir
       ir=ieor(ir,ishft(ir,5))+1422217823 !
       !The various constants are chosen to give good bit mixing and 
       !should not be changed.
       ir=ieor(ir,ishft(ir,-16))+1842055030
       ir=ieor(ir,ishft(ir,9))+80567781
       ir=ieor(il,ir)
       il=is
    end do
  END SUBROUTINE ran_hash_s
!===
! RAN_HASH_V(IL,IR)
!===
  SUBROUTINE ran_hash_v(il,ir)
    IMPLICIT NONE
    INTEGER(K4B), DIMENSION(:), INTENT(INOUT) :: il,ir
    !Vector version of ran hash s.
    INTEGER(K4B), DIMENSION(size(il)) :: is
    INTEGER(K4B) :: j
    do j=1,4
       is=ir
       ir=ieor(ir,ishft(ir,5))+1422217823
       ir=ieor(ir,ishft(ir,-16))+1842055030
       ir=ieor(ir,ishft(ir,9))+80567781
       ir=ieor(il,ir)
       il=is
    end do
  END SUBROUTINE ran_hash_v
END MODULE ran_state
