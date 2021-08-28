!======================================================================
!
! SUBROUTINE: Notify end of a subroutine
!
!      [in ] strnm: subroutine name
!
! Programmed by Takazumi Yamaguchi
!======================================================================
subroutine sucmp(strnm)
  implicit none
  character(len=*),intent(in) :: strnm
  !
  write(*,'(/a,1x,a,1x,a/)') &
       '.... Subroutine', trim(strnm), 'has been completed ! ....'
  !
  return
end subroutine sucmp
