!======================================================================
!
! SUBROUTINE: Solve ordinary differential equations by using 
!             2nd-order Runge-Kutta method
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
subroutine rk2(neq,t0,tn,x0,xn,fsub)
  real(8),intent(in ) :: t0,tn,x0(neq)
  real(8),intent(out) :: xn(neq)
  real(8) :: dt,k1(neq),k2(neq)
  real(8),parameter   :: two=2.0D0
  interface
     subroutine fsub(neq,t,x,f)
       integer,intent(in ) :: neq
       real(8),intent(in ) :: t,x(neq)
       real(8),intent(out) :: f(neq)
     end subroutine fsub
  end interface
  !
  dt = tn - t0
  ! ... Calc. coefficient for runge-kutta method
  call fsub(neq,t0,x0,k1)
  k1 = dt*k1(:)
  call fsub(neq,tn,x0,k2)
  k2 = dt*k2(:)
  !
  xn(:) = x0(:) + (k1(:) + k2(:))/two
  !
  return
end subroutine rk2
