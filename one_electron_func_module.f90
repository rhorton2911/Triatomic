module one_electron_func

  use spheroidal_class
  use sturmian_class

  integer:: basis_type ! =0 if  nonrel. basis is in use, or =1 if rel.
  type(basis_sturmian_nr) :: bst_nr   ! this is Sturmian basis  (depend on k and l)
  type(spheroidal_basis) :: bst_oid   ! More general basis for spheroidal calculations: can include sturmians as well as separated one-electron states.
  integer:: nspm

end module one_electron_func

