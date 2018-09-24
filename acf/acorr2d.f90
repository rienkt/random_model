! Testing program  for 2-D autocorrelation function 
!
!USE nrtype
USE acorr
!USE nr, ONLY : realft
IMPLICIT NONE
INTEGER :: n1,n2,lag1max,lag2max
PARAMETER(n1=128,n2=64,lag1max=127, lag2max=63)
!PARAMETER(n1=1,n2=128,lag1max=1, lag2max=127)
INTEGER :: ilag,jlag,i,j,no2
REAL,DIMENSION(n1,n2) :: org_data
! DIRECT METHOD
REAL, DIMENSION(n1*2,-n2:n2*2) :: a
REAL, DIMENSION(-lag1max:lag1max,-lag2max:lag2max) :: b!,b2
REAL, DIMENSION(lag1max*2+1,lag2max*2+1) :: b2
REAL :: b0
! Data Input
open(100,FILE="temp2d.dat")
do j=1,n2
   do i=1,n1
      read(100,*) org_data(i,j)
   enddo
enddo
close(100)
print *, size(org_data)

!Direct Method
a=0.d0
a(1:n1,1:n2)=org_data
!a(1:n1,-1:-n2:-1)=org_data
!print *, (a(1,i),i=-4,5)
!print *, (a(2,i),i=-4,5)
b=0.d0

do ilag=0,lag1max
!   jlag=0
!   do i=1,n1
!      do j=1,n2
!         b(ilag,jlag)=b(ilag,jlag)+a(i,j)*a(i+ilag,j+jlag)
!      enddo
!   end do
!    if (lag2max .gt. 0) then
       do jlag=-lag2max,lag2max
          do i=1,n1
             do j=1,n2
                b(ilag,jlag)=b(ilag,jlag)+a(i,j)*a(i+ilag,j+jlag)
!                b(ilag,-jlag)=b(ilag,-jlag)+a(i,j)*a(i+ilag, j-jlag)
             enddo
          enddo
       enddo
!    endif
enddo
b0=b(0,0)
b=b/b0
!write(*,'(F12.5)') (b(i,0),i=0,lag1max)

print *, size(b,1),size(b,2),-1+lag1max, -1+lag2max
print *, size(b,1),size(b,2)

b(-lag1max:-1,0:lag2max)=b(lag1max:1:-1,0:-lag2max:-1)
b(-lag1max:-1,-lag2max:-1)=b(lag1max:1:-1,lag2max:1:-1)
print *, size(b,1),size(b,2)
open(100,FILE="test2dcorr1.dat")
!write(100,'(F12.5)') (b(i,0),i=0,lag1max)
write(100,'(F12.5)') (b(0,j),j=-lag2max,lag2max)
close(100)
open(100,FILE="test2dcorr.dat")
write(100,'(F12.5)') ((b(i,j),i=-lag1max,lag1max),j=-lag2max,lag2max)
!write(100,'(F9.5)') (b(i,0),i=-lag1max,lag1max)
close(100)

! TEST AGAIN
b2=0.d0
call acorr2d(org_data,b2,lag1max,lag2max)
open(100,FILE="test2dcorrsub.dat")
write(100,'(F12.5)') ((b2(i,j),i=1,lag1max*2+1),j=1,lag2max*2+1)
!write(100,'(F9.5)') (b(i,0),i=-lag1max,lag1max)
close(100)

end
