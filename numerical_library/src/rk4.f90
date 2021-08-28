!======================================================================
!
! SUBROUTINE: Solve ordinary differential equations by using 
!             4nd-order Runge-Kutta method
!
!      [in ] neq : number of equations.
!            t0  : time at initial condition.
!            tn  : time at netx time step.
!            x0  : inidial condition.
!            fsub: subroutine for rigit-hand-side vector.
!
!      [out] xn  : solution of ordinary differential equations.
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine rk4(neq,t0,tn,x0,xn,fsub)
  real(8),intent(in ) :: t0,tn,x0(neq)
  real(8),intent(out) :: xn(neq)
  integer :: i,j
  real(8) :: dt,c(4),tx,tmp(neq),K(neq,0:4),f(neq)
  interface
     subroutine fsub(neq,t,x,f)
       integer,intent(in ) :: neq
       real(8),intent(in ) :: t,x(neq)
       real(8),intent(out) :: f(neq)
     end subroutine fsub
  end interface
  !
  dt = tn - t0
  c(1:4) = (/0d0,0.5d0,0.5d0,1d0/)
  tmp = 0.0d0
  K = 0d0
  f = 0d0  
  !
  do j=1, 4
     do i=1, neq
        tmp(i) = x0(i) + K(i,j-1)*c(j)
     end do
     !
     tx = t0 + c(j)*dt
     call fsub(neq,tx,tmp,f)
     !
     do i=1, neq
        K(i,j) = dt*f(i)
     end do
  end do
  !
  do i=1,neq
     xn(i) = x0(i) + (K(i,1)+K(i,4))/6d0+(K(i,2)+K(i,3))/3d0
  enddo
  !
  return
end subroutine rk4
