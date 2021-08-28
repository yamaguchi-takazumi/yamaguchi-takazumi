program main
  implicit none
  integer,parameter :: n=20
  integer :: i,k,kmax,n1,n2,ip(n)
  real(8) :: adim(n,n),cdim(n),xsol(n),xgauss(n)
  real(8) :: harv1,harv2,work(n)
  real(8),parameter::zero=0.0d0,one=1.0d0,two=2.0d0
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
  ! ... exchanging two different rows.
  kmax = n/2
  do k=1,kmax
     call random_number(harv1)
     call random_number(harv2)
     n1 = int(dble(n-1)*harv1) + 1
     n2 = int(dble(n-1)*harv2) + 1
     if(n1 /= n2) then
        work(:)      = adim(n1,:)
        adim(n1,:) = adim(n2,:)
        adim(n2,:) = work
     end if
  end do
  ! ... setting inhomogenious vector cdim.
  do i=1, n
     work(i) = dble(mod(i,6)+1)
  end do
  cdim = matmul(adim,work)
  !
  ! solution of linear system
  !
  call gauss(n,adim,cdim,xgauss)
  call lu_decomp(n,adim,ip)
  call lu_solve(n,adim,cdim,ip,xsol)
  !
  ! output solotions into system standard output.
  !
  write(*,'(/a)') '**** solution of simultaneous equation ****'
  write(*,'(3x,a7,4x,a7)') 'lu_decomp','gauss'
  do i=1,n
     write(*,'(2(1pe10.3,1x))') xsol(i),xgauss(i)
  end do
  write(*,*) 'ip:'
  write(*,*) ip(1:n)
  !
  stop
end program main
