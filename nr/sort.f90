SUBROUTINE sort_sp(arr)
  USE nrtype; USE nrutil, ONLY : swap,nrerror
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
  INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
  ! Sorts an array arr into ascending numerical order using the Quicksort 
  ! algorithm. arr is replaced on output by its sorted rearrangement.
  ! Parameters: NN is the size of subarrays sorted by straight insertion and 
  ! NSTACK is the required auxiliary storage.
  REAL(SP) :: a
  INTEGER(I4B) :: n,k,i,j,jstack,l,r
  INTEGER(I4B), DIMENSION(NSTACK) :: istack
  n=size(arr)
  jstack=0
  l=1
  r=n
  do
     if (r-l < NN) then !Insertion sort when subarray small enough.
        do j=l+1,r
           a=arr(j)
           do i=j-1,l,-1
              if (arr(i) <= a) exit
              arr(i+1)=arr(i)
           end do
           arr(i+1)=a
        end do
        if (jstack == 0) RETURN
        r=istack(jstack) !Pop stack and begin a new round of partitioning.
        l=istack(jstack-1)
        jstack=jstack-2
     else ! Choose median of left, center, and right elements
          ! as partitioning element a. Also rearrange so
          ! that a(l) =< a(l+1) =< a(r).
        k=(l+r)/2
        call swap(arr(k),arr(l+1))
        call swap(arr(l),arr(r),arr(l)>arr(r))
        call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
        call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
        i=l+1 ! Initialize pointers for partitioning.
        j=r
        a=arr(l+1) !Partitioning element.
        do !Here is the meat.
           do !Scan up to find element >= a.
              i=i+1
              if (arr(i) >= a) exit
           end do
           do !Scan down to find element <= a.
              j=j-1
              if (arr(j) <= a) exit
           end do
           if (j < i) exit !Pointers crossed. Exit with partitioning complete.
           call swap(arr(i),arr(j)) !Exchange elements.
        end do
        arr(l+1)=arr(j) !Insert partitioning element.
        arr(j)=a
        jstack=jstack+2
        ! Push pointers to larger subarray on stack; 
        ! process smaller subarray immediately.
        if (jstack > NSTACK) call nrerror('sort: NSTACK too small')
        if (r-i+1 >= j-l) then
           istack(jstack)=r
           istack(jstack-1)=i
           r=j-1
        else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
        end if
     end if
  end do
END SUBROUTINE sort_sp
SUBROUTINE sort_dp(arr)
  USE nrtype; USE nrutil, ONLY : swap,nrerror
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
  INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
  ! Sorts an array arr into ascending numerical order using the Quicksort 
  ! algorithm. arr is replaced on output by its sorted rearrangement.
  ! Parameters: NN is the size of subarrays sorted by straight insertion and 
  ! NSTACK is the required auxiliary storage.
  REAL(DP) :: a
  INTEGER(I4B) :: n,k,i,j,jstack,l,r
  INTEGER(I4B), DIMENSION(NSTACK) :: istack
  n=size(arr)
  jstack=0
  l=1
  r=n
  do
     if (r-l < NN) then !Insertion sort when subarray small enough.
        do j=l+1,r
           a=arr(j)
           do i=j-1,l,-1
              if (arr(i) <= a) exit
              arr(i+1)=arr(i)
           end do
           arr(i+1)=a
        end do
        if (jstack == 0) RETURN
        r=istack(jstack) !Pop stack and begin a new round of partitioning.
        l=istack(jstack-1)
        jstack=jstack-2
     else ! Choose median of left, center, and right elements
          ! as partitioning element a. Also rearrange so
          ! that a(l) =< a(l+1) =< a(r).
        k=(l+r)/2
        call swap(arr(k),arr(l+1))
        call swap(arr(l),arr(r),arr(l)>arr(r))
        call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
        call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
        i=l+1 ! Initialize pointers for partitioning.
        j=r
        a=arr(l+1) !Partitioning element.
        do !Here is the meat.
           do !Scan up to find element >= a.
              i=i+1
              if (arr(i) >= a) exit
           end do
           do !Scan down to find element <= a.
              j=j-1
              if (arr(j) <= a) exit
           end do
           if (j < i) exit !Pointers crossed. Exit with partitioning complete.
           call swap(arr(i),arr(j)) !Exchange elements.
        end do
        arr(l+1)=arr(j) !Insert partitioning element.
        arr(j)=a
        jstack=jstack+2
        ! Push pointers to larger subarray on stack; 
        ! process smaller subarray immediately.
        if (jstack > NSTACK) call nrerror('sort: NSTACK too small')
        if (r-i+1 >= j-l) then
           istack(jstack)=r
           istack(jstack-1)=i
           r=j-1
        else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
        end if
     end if
  end do
END SUBROUTINE sort_dp



SUBROUTINE indexx_sp(arr,index)
  USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror,swap
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: arr
  INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
  INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
  ! Indexes an array arr, i.e., outputs the array index of length N such 
  ! that arr(index(j)) is in ascending order for j = 1, 2, . . . ,N. The 
  ! input quantity arr is not changed.
  REAL(SP) :: a
  INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
  INTEGER(I4B), DIMENSION(NSTACK) :: istack
  n=assert_eq(size(index),size(arr),'indexx_sp')
  index=arth(1,1,n)
  jstack=0
  l=1
  r=n
  do
     if (r-l < NN) then
        do j=l+1,r
           indext=index(j)
           a=arr(indext)
           do i=j-1,l,-1
              if (arr(index(i)) <= a) exit
              index(i+1)=index(i)
           end do
           index(i+1)=indext
        end do
        if (jstack == 0) RETURN
        r=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
     else
        k=(l+r)/2
        call swap(index(k),index(l+1))
        call icomp_xchg(index(l),index(r))
        call icomp_xchg(index(l+1),index(r))
        call icomp_xchg(index(l),index(l+1))
        i=l+1
        j=r
        indext=index(l+1)
        a=arr(indext)
        do
           do
              i=i+1
              if (arr(index(i)) >= a) exit
           end do
           do
              j=j-1
              if (arr(index(j)) <= a) exit
           end do
           if (j < i) exit
           call swap(index(i),index(j))
        end do
        index(l+1)=index(j)
        index(j)=indext
        jstack=jstack+2
        if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
        if (r-i+1 >= j-l) then
           istack(jstack)=r
           istack(jstack-1)=i
           r=j-1
        else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
        end if
     end if
  end do
CONTAINS
  SUBROUTINE icomp_xchg(i,j)
    INTEGER(I4B), INTENT(INOUT) :: i,j
    INTEGER(I4B) :: swp
    if (arr(j) < arr(i)) then
       swp=i
       i=j
       j=swp
    end if
  END SUBROUTINE icomp_xchg
END SUBROUTINE indexx_sp

SUBROUTINE indexx_i4b(iarr,index)
  USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror,swap
  IMPLICIT NONE
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
  INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: index
  INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
  INTEGER(I4B) :: a
  INTEGER(I4B) :: n,k,i,j,indext,jstack,l,r
  INTEGER(I4B), DIMENSION(NSTACK) :: istack
  n=assert_eq(size(index),size(iarr),'indexx_sp')
  index=arth(1,1,n)
  jstack=0
  l=1
  r=n
  do
     if (r-l < NN) then
        do j=l+1,r
           indext=index(j)
           a=iarr(indext)
           do i=j-1,l,-1
              if (iarr(index(i)) <= a) exit
              index(i+1)=index(i)
           end do
           index(i+1)=indext
        end do
        if (jstack == 0) RETURN
        r=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
     else
        k=(l+r)/2
        call swap(index(k),index(l+1))
        call icomp_xchg(index(l),index(r))
        call icomp_xchg(index(l+1),index(r))
        call icomp_xchg(index(l),index(l+1))
        i=l+1
        j=r
        indext=index(l+1)
        a=iarr(indext)
        do
           do
              i=i+1
              if (iarr(index(i)) >= a) exit
           end do
           do
              j=j-1
              if (iarr(index(j)) <= a) exit
           end do
           if (j < i) exit
           call swap(index(i),index(j))
        end do
        index(l+1)=index(j)
        index(j)=indext
        jstack=jstack+2
        if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
        if (r-i+1 >= j-l) then
           istack(jstack)=r
           istack(jstack-1)=i
           r=j-1
        else
           istack(jstack)=j-1
           istack(jstack-1)=l
           l=i
        end if
     end if
  end do
CONTAINS
  SUBROUTINE icomp_xchg(i,j)
    INTEGER(I4B), INTENT(INOUT) :: i,j
    INTEGER(I4B) :: swp
    if (iarr(j) < iarr(i)) then
       swp=i
       i=j
       j=swp
    end if
  END SUBROUTINE icomp_xchg
END SUBROUTINE indexx_i4b


FUNCTION select_s(k,arr)
  USE nrtype; USE nrutil, ONLY : assert,swap
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: k
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr
  REAL(SP) :: select_s
  ! Returns the kth smallest value in the array arr. The input array will 
  ! be rearranged to have this value in location arr(k), with all smaller 
  ! elements moved to arr(1:k-1) (in arbitrary order) and all larger elements 
  ! in arr(k+1:) (also in arbitrary order).
  INTEGER(I4B) :: i,r,j,l,n
  REAL(SP) :: a
  n=size(arr)
  call assert(k >= 1, k <= n, 'select args')
  l=1
  r=n
  do
     if (r-l <= 1) then !Active partition contains 1 or 2 elements.
        if (r-l == 1) call swap(arr(l),arr(r),arr(l)>arr(r)) 
                        !Active partition contains 2 elements. 
        select_s=arr(k)
        RETURN
     else !Choose median of left, center, and right elements
          ! as partitioning element a. Also rearrange so that 
          ! arr(l) <= arr(l+1) <= arr(r).
        i=(l+r)/2
        call swap(arr(i),arr(l+1))
        call swap(arr(l),arr(r),arr(l)>arr(r))
        call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
        call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
        i=l+1 ! Initialize pointers for partitioning.
        j=r
        a=arr(l+1) ! Partitioning element.
        do !Here is the meat.
           do !Scan up to find element > a.
              i=i+1
              if (arr(i) >= a) exit
           end do
           do !Scan down to find element < a.
              j=j-1
              if (arr(j) <= a) exit
           end do
           if (j < i) exit ! Pointers crossed. Exit with partitioning complete.
           call swap(arr(i),arr(j)) !Exchange elements.
        end do
        arr(l+1)=arr(j) ! Insert partitioning element.
        arr(j)=a
        if (j >= k)  r=j-1 ! Keep active the partition that contains 
                           ! the kth element.
        if (j <= k) l=i
     end if
  end do
END FUNCTION select_s

FUNCTION select_inplace_s(k,arr)
  USE nrtype
  USE nr, ONLY : select
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: k
  REAL(SP), DIMENSION(:), INTENT(IN) :: arr
  REAL(SP) :: select_inplace_s
  ! Returns the kth smallest value in the array arr, without altering 
  ! the input array. In Fortran 90fs assumed memory-rich environment, 
  ! we just call select in scratch space.
  REAL(SP), DIMENSION(size(arr)) :: tarr
  tarr=arr
  select_inplace_s=select(k,tarr)
END FUNCTION select_inplace_s



FUNCTION select_d(k,arr)
  USE nrtype; USE nrutil, ONLY : assert,swap
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: k
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: arr
  REAL(DP) :: select_d
  ! Returns the kth smallest value in the array arr. The input array will 
  ! be rearranged to have this value in location arr(k), with all smaller 
  ! elements moved to arr(1:k-1) (in arbitrary order) and all larger elements 
  ! in arr(k+1:) (also in arbitrary order).
  INTEGER(I4B) :: i,r,j,l,n
  REAL(DP) :: a
  n=size(arr)
  call assert(k >= 1, k <= n, 'select args')
  l=1
  r=n
  do
     if (r-l <= 1) then !Active partition contains 1 or 2 elements.
        if (r-l == 1) call swap(arr(l),arr(r),arr(l)>arr(r)) 
                        !Active partition contains 2 elements. 
        select_d=arr(k)
        RETURN
     else !Choose median of left, center, and right elements
          ! as partitioning element a. Also rearrange so that 
          ! arr(l) <= arr(l+1) <= arr(r).
        i=(l+r)/2
        call swap(arr(i),arr(l+1))
        call swap(arr(l),arr(r),arr(l)>arr(r))
        call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
        call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
        i=l+1 ! Initialize pointers for partitioning.
        j=r
        a=arr(l+1) ! Partitioning element.
        do !Here is the meat.
           do !Scan up to find element > a.
              i=i+1
              if (arr(i) >= a) exit
           end do
           do !Scan down to find element < a.
              j=j-1
              if (arr(j) <= a) exit
           end do
           if (j < i) exit ! Pointers crossed. Exit with partitioning complete.
           call swap(arr(i),arr(j)) !Exchange elements.
        end do
        arr(l+1)=arr(j) ! Insert partitioning element.
        arr(j)=a
        if (j >= k)  r=j-1 ! Keep active the partition that contains 
                           ! the kth element.
        if (j <= k) l=i
     end if
  end do
END FUNCTION select_d

FUNCTION select_inplace_d(k,arr)
  USE nrtype
  USE nr, ONLY : select
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: k
  REAL(DP), DIMENSION(:), INTENT(IN) :: arr
  REAL(DP) :: select_inplace_d
  ! Returns the kth smallest value in the array arr, without altering 
  ! the input array. In Fortran 90fs assumed memory-rich environment, 
  ! we just call select in scratch space.
  REAL(DP), DIMENSION(size(arr)) :: tarr
  tarr=arr
  select_inplace_d=select(k,tarr)
END FUNCTION select_inplace_d
