module vmat_exch_module
  use sturmian_class
  use state_class

  implicit none

  type(basis_sturmian_nr):: bst_orth
  type(basis_state), target, public :: TargState_Orthonorm ! Used to build projection operator for H2 
  integer, allocatable, dimension(:) :: is_core_MO ! core type one-electron target states
  integer, allocatable, dimension(:) :: is_core_POMO ! Projection operator core type one-electron target states
  real*8, allocatable, dimension(:,:) :: ovlp_OrtnrmMO_k  ! Projection operator overlaps <varphi|k>

end module vmat_exch_module


