module vmat_module
!  include 'par.h'
integer, parameter :: nch_array_max=2000 
real*8, dimension(:,:), allocatable:: vmat,vmat01, vmat0,vmat1
complex*16, dimension(:,:), allocatable:: CVmat
integer, dimension(nch_array_max) :: nchistart, nchistop
logical:: scalapack
end module vmat_module
