program main
  implicit none
  integer,parameter :: n=20,itmax=n
  integer :: i,itnum
  real(8) :: adim(n,n),cdim(n),xsol(n),xgauss(n),work(n)
  real(8),parameter::zero=0.0d0,one=1.0d0,two=2.0d0,epsln=1.0d-12
  !
  ! setting matrix and vectors.
  !
  adim(:,:) = zero
  adim(1,1) =  one
  adim(1,2) = -one
  do i=2, n-1
     adim(i,i  ) =  two
     adim(i,i-1) = -one
     adim(i,i+1) = -one
  end do
  adim(n,n-1) = -one
  adim(n,n  ) =  two
  ! ... setting inhomogenious vector cdim.
  do i=1, n
     work(i) = dble(mod(i,6)+1)
  end do
  cdim = matmul(adim,work)
  !
  ! solution of linear system
  !
  xsol(:) = one
  call gauss(n,adim,cdim,xgauss)
  call cgm(n,adim,cdim,itmax,epsln,.true.,xsol,itnum)
  !
  ! output solotions into system standard output.
  !
  write(*,'(/a)') '**** solution of simultaneous equation ****'
  write(*,'(3x,a7,4x,a7)') 'cgm','gauss'
  do i=1,n
     write(*,'(2(1pe10.3,1x))') xsol(i),xgauss(i)
  end do
  !
  stop
end program main
