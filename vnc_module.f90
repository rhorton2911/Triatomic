module vnc_module
  use sturmian_class
  use state_class
  integer:: lamtop_vc ! limit for vnc
  integer:: ltop_vc  ! limit for vdcore
  real*8, dimension(:,:), allocatable:: vnc
  integer, dimension(:), allocatable:: minvnc, maxvnc
  real*8, dimension(:,:), allocatable:: vdcore
  integer:: minvdc, maxvdc
  real*8, dimension(:), allocatable:: vcentre_jn  ! joint nuclear potential
  type(basis_sturmian_nr):: wfcore_nr
  real*8 :: core_energy
  type(basis_state) :: CoreOrbitals !Liam added: to store the orbitals describing the inert core

end module vnc_module
