include 'blck0'
include 'input'
include 'model'
include 'poisson_eq'
include 'solve'
include 'matrix_A'
include 'matrix_B'
include 'vector_C'
include 'output'
include 'result_error'
include 'pydata'
!
program main
  call input
  call model
  call solve
  call output
  !
  stop
end program main
