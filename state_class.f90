module state_class
	 use numbers  !Reese: contains format specifier 'dpf' for double precision reals

!   private
  type, public :: state
     logical:: hlike
     character(len=5) :: label ! Label representing the quantum numbers (e.g. ground state is 1sSg)
     character(len=4) :: domconfig
     real(dpf):: M ! total angular momentum projection of the state, mmajor config for non-linear targets
     integer:: parity  ! parity (+1, -1) of the state
     real(dpf):: spin  ! spin = 0,1 for two-electron targets, 1/2 for one-electron targets
     real(dpf):: energy
     real(dpf), dimension(:), allocatable:: CI
     integer, dimension(:), allocatable:: na, nb
     integer, dimension(:), allocatable:: ma, mb  ! magnetic sublevels
     integer:: nam   !  size of array na(:), 
     integer:: nusemax ! number of orbitals used to describe this state (useful for two-electron states only)
     integer, dimension(:), allocatable:: nuse   ! idex to orbitals used to describe the state
     real(dpf), dimension(:), allocatable :: SF_core

     integer:: l   ! major config l value
     integer:: n   ! major config n value
		 ! major config m value stored in M above

     integer:: inum ! state index for given target symmetry

  end type state
!
!
  type, public:: basis_state
     integer:: basis_type ! =0 if  nonrel. basis is in use, or =1 if rel.
     logical:: hlike
     real(dpf):: Mmax
     type(state), dimension(:), pointer :: b  ! array of states
     integer:: Nmax    ! number of states (size of the basis)
     integer:: Nstates    ! last initialized state
     integer:: nicm  ! number of core orbitals, set in H12.f90
     integer, dimension(:), allocatable:: ncore, mncore    ! index to core orbitals, allocated in H12.f90
     real(dpf):: en_ion   ! energy of one-el. ion ground state
     integer:: Nmax_bound
     integer:: Nmax_open
     real(dpf), dimension(:,:), allocatable:: br_ratio

  end type basis_state

  public :: get_ang_mom_proj, get_par, get_spin, get_energy, new_basis_st, destruct_basis_st, sort_by_energy_basis_st, basis_size,print_energy_basis_st, calc_spectro_factors,  construct_st, destruct_st, copy_st, get_nusemax, get_nam, get_nuse, get_na, get_ma,get_CI, ovlp_st, modify_CI_state, get_n_majconf, set_n_majconf, get_l_majconf, set_l_majconf, set_label, get_label, write_1elState
  
  
  interface get_n_majconf
     module procedure get_n_majconf
  end interface

  interface get_l_majconf
     module procedure get_l_majconf
  end interface

  !interface get_m_majconf
  !   module procedure get_m_majconf
  !end interface

  interface get_ang_mom_proj
     module procedure get_ang_mom_proj_st
  end interface

  interface get_par
     module procedure get_par_st
  end interface

  interface get_spin
     module procedure get_spin_st
  end interface

  interface get_nusemax
     module procedure get_nusemax_st
  end interface

  interface get_nam
     module procedure get_nam_st
  end interface

  interface get_nuse
     module procedure get_nuse_st
  end interface     

  interface get_na
     module procedure get_na_st
  end interface

  interface get_ma
     module procedure get_ma_st
  end interface

  interface set_label
     module procedure set_label_1e, set_label_2e
  end interface

  interface get_label
     module procedure get_label_st
  end interface

  interface get_energy
     module procedure get_energy_st
  end interface

  interface get_elec_energy
     module procedure get_elec_energy_st
  end interface
  
  interface print_energy
     module procedure print_energy_basis_st
  end interface

  interface new_basis
     module procedure new_basis_st
  end interface

  interface basis_size
     module procedure basis_size_st
  end interface

  interface copy
     module procedure copy_st, copy_basis_st
  end interface


contains
!
!
  subroutine  construct_st(self,hlike,m,parity,spin,energy,inum,ncm,CI,no1,mo1,no2,mo2,phase)
	  use input_data  !data_in
    implicit none

    type(state), intent(inout):: self
    logical, intent(in):: hlike
    real(dpf), intent(in):: m  
    integer, intent(in):: parity 
    real(dpf), intent(in):: spin
    real(dpf), intent(in):: energy
    integer, intent(in):: inum
    integer, intent(in):: ncm
    real(dpf), dimension(ncm), intent(in) :: CI
    integer, dimension(ncm), intent(in) :: no1, no2, phase, mo1, mo2
    optional:: no2, phase, mo2
    integer:: i
		logical:: goodm

    self%hlike = hlike
    self%m = m
    self%parity = parity
    self%spin = spin
    self%energy = energy

    self%l = -1
    self%n = -1

    self%inum = inum

    if(hlike) then
       allocate(self%CI(ncm),self%na(ncm),self%ma(ncm))
       self%nam=ncm
       self%na(1:ncm) = no1(1:ncm)
       self%ma(1:ncm) = mo1(1:ncm)
       do i = 1, ncm
          if ( data_in%good_m .and. (self%m /= self%ma(i)) ) then
             print*,"Error construct_st: angular projection for target state m and mo1 different"
             print*, 'self%m, self%ma(i) =', self%m, self%ma(i)
             stop
          end if
       end do
       self%CI(1:ncm) = CI(1:ncm)
    else
       self%label = "  -  "
       call setupCI(self,ncm,CI,no1,mo1,no2,mo2,phase)
    endif

!    print*,'created a state with energy:', energy
   
  end   subroutine  construct_st
  !
  !
  subroutine  setupCI(self,ncm,CI,no1,mo1,no2,mo2,phase)    !  to be modified for new definition of CI
    implicit none

    type(state), intent(inout):: self
    integer, intent(in):: ncm
    real(dpf), dimension(ncm), intent(in) :: CI
    integer, dimension(ncm), intent(in) :: no1, no2, mo1,mo2, phase

    integer:: i, nc, imax, ji, n1, itest
    integer, dimension(:), allocatable:: ntmp

     
    ! How many unsymetrized configurations are used to describe this state
    i = 0    
    do nc=1,ncm
       i = i + 1
       if( no1(nc) .eq. no2(nc) .and. mo1(nc) .eq. mo2(nc)) then   ! |aa>
          continue
       else   ! |ab>
          i = i + 1   ! also include |ba>
       endif
    enddo
    imax = i

    self%nam = i
    if(allocated(self%na)) then
       deallocate(self%CI,self%na,self%ma,self%nb,self%mb)
    endif
    allocate(self%na(i),self%ma(i),self%nb(i),self%mb(i))
    allocate(self%CI(i))

    self%CI(:) = 0d0

    i = 0
    do nc=1,ncm
       i = i + 1
       self%na(i) = no1(nc)
       self%nb(i) = no2(nc)
       self%ma(i) = mo1(nc)
       self%mb(i) = mo2(nc)
       if( no1(nc) .eq. no2(nc) .and. mo1(nc) .eq. mo2(nc)) then   ! |aa>
          self%CI(i) = CI(nc)
       else   ! |ab> (unsymmetrised). Therefore we want [ |ab> + (-1)^S |ba> ] / sqrt(2)
          self%CI(i) = CI(nc)/sqrt(2d0)

          i = i + 1
          self%CI(i) = phase(nc)*CI(nc)/sqrt(2d0)
          self%na(i) = no2(nc)
          self%nb(i) = no1(nc)
          self%ma(i) = mo2(nc)
          self%mb(i) = mo1(nc)

       endif       
    enddo

!!$  define array that has all orbitals that are used to describe this state
    ! first define auxiliary array
    allocate(ntmp(1:imax))
    ntmp(1:imax) = 0
    ji = 0
    do nc=1,imax
       n1 = self%na(nc)
       itest = 1
       do i=1,ji
          if(n1 .eq. ntmp(i)) then
             itest = 0
             exit
          endif
       enddo
       if(itest .eq. 1) then
          ji = ji + 1
          ntmp(ji) = n1
       endif
    enddo
!    print*, 'ji=', ji

    self%nusemax = ji
    if(allocated(self%nuse)) deallocate(self%nuse)
    allocate(self%nuse(1:ji))
    self%nuse(1:ji) = ntmp(1:ji)
    deallocate(ntmp)
    
  end   subroutine setupCI
  !
  ! 
  subroutine modify_CI_state(self,nam_in,na_in,C_in, ma_in)
    implicit none  
    
    type(state), intent(inout):: self
    integer, intent(in):: nam_in
    integer, dimension(nam_in), intent(in):: na_in
    real(dpf), dimension(nam_in), intent(in):: C_in
    integer, dimension(nam_in), intent(in):: ma_in
    integer:: nam

    if(nam_in .ne. self%nam) then
       if(allocated(self%CI)) then
          deallocate(self%CI,self%na, self%ma)
       endif
       nam = nam_in
       allocate(self%CI(nam),self%na(nam),self%ma(nam))       
    endif
    self%nam = nam_in
    self%na(:) = na_in(:)
    self%CI(:) = C_in(:)
    self%ma(:) = ma_in(:)

  end subroutine modify_CI_state
!
! Shellsort algorithm for integer array na(:)
  subroutine sort_by_value(N,na)
    implicit none

    integer, intent(in):: N
    integer, dimension(N), intent(inout):: na
    integer:: gap, i, j
    integer:: Tmp

    gap = N/2
    do 
       if(gap .le. 0) exit
       do i=gap+1,N
          Tmp = na(i)
          do j=i,gap+1,-gap
             if(Tmp.lt. na(j-gap)) then
                na(j) = na(j-gap)
             else
                exit
             endif
             na(j-gap) = Tmp
          enddo
       enddo
       if ( gap .eq. 2) then
          gap = 1
       else
          gap = nint( real(gap)/2.2 )
       endif       
    enddo

  end subroutine sort_by_value
  !
  !
  subroutine  destruct_st(self)
    implicit none
    type(state), intent(inout):: self
    
    if(allocated(self%CI)) deallocate(self%CI)
    if(allocated(self%na)) deallocate(self%na)
    if(allocated(self%nb)) deallocate(self%nb)
    if(allocated(self%ma)) deallocate(self%ma)
    if(allocated(self%mb)) deallocate(self%mb)
    self%label = 'xxxxx'
    
  end subroutine destruct_st
  !
  subroutine copy_st(state_l,state_r)
    implicit none

    type(state), intent(out):: state_l
    type(state), intent(in):: state_r
    integer:: i1,i2, nam
    

    state_l%hlike = state_r%hlike

    nam = state_r%nam
    
    state_l%energy = state_r%energy
    state_l%M = state_r%M
    state_l%parity = state_r%parity
    state_l%spin = state_r%spin
    
    state_l%l = state_r%l
    state_l%n = state_r%n
    state_l%label = state_r%label
    state_l%domconfig = state_r%domconfig
    state_l%inum = state_r%inum

    if(allocated(state_l%CI)) then
       deallocate(state_l%CI,state_l%na,state_l%ma)
    endif
    state_l%nam = state_r%nam

    allocate(state_l%na(nam))
    allocate(state_l%CI(nam))
    allocate(state_l%ma(nam))

    state_l%CI = state_r%CI
    state_l%na = state_r%na
    state_l%ma = state_r%ma

    if(allocated(state_l%nb)) then
       deallocate(state_l%nb)
       deallocate(state_l%mb)
    endif

    if(allocated(state_r%SF_core)) then
      allocate(state_l%SF_core(size(state_r%SF_core)))
      state_l%SF_core = state_r%SF_core
    endif

    if(allocated(state_l%nuse)) deallocate(state_l%nuse)
    if(.not. state_l%hlike) then
       allocate(state_l%nb(nam))
       state_l%nb = state_r%nb    
       allocate(state_l%mb(nam))
       state_l%mb = state_r%mb    
       if(allocated(state_r%nuse)) then !nuse not always going to be used because it is made during rearrange
         state_l%nusemax = state_r%nusemax
         allocate(state_l%nuse(state_l%nusemax))
         state_l%nuse = state_r%nuse
       endif
    endif
    
  end subroutine copy_st
!
  subroutine copy_H2plus_st(state_l,state_r,Nmax)
    ! Mark: This subroutine copies the H2+ ground state, 
    ! which is built from Laguerre functions in the data.in file.
    ! We then create a basis from the F5 file to describe the excited A.O.
    ! The data.in Lagueere functions are added ontop of the F5 file basis.
    implicit none

    type(state), intent(out):: state_l
    type(state), intent(in):: state_r
    integer, intent(in):: Nmax
    integer:: i1,i2, nam
    

    state_l%hlike = state_r%hlike

    nam = state_r%nam
    
    state_l%energy = state_r%energy
    state_l%M = state_r%M
    state_l%parity = state_r%parity
    state_l%spin = state_r%spin
    
    state_l%l = state_r%l
    state_l%n = state_r%n

    state_l%inum = state_r%inum

    if(allocated(state_l%CI)) then
       deallocate(state_l%CI,state_l%na,state_l%ma)
    endif
    state_l%nam = state_r%nam
    state_l%label = state_r%label

    allocate(state_l%na(nam))
    allocate(state_l%CI(nam))
    allocate(state_l%ma(nam))

    state_l%CI = state_r%CI
    state_l%na = state_r%na + Nmax ! F5 file functions copied
    state_l%ma = state_r%ma
    
  end subroutine copy_H2plus_st
!
  function get_energy_st(self)
    implicit none
    real(dpf):: get_energy_st
    type(state), intent(in):: self
    
    get_energy_st = self%energy
  end function get_energy_st
!
  function get_elec_energy_st(self,n)
    implicit none
    real(dpf) :: get_elec_energy_st
    type(basis_state), intent(in) :: self
    integer, intent(in) :: n
    
    get_elec_energy_st = self%b(n)%energy - self%en_ion
  end function get_elec_energy_st
!
!
  function get_ang_mom_proj_st(self)
    implicit none
    integer :: get_ang_mom_proj_st
    type(state), intent(in) :: self

    get_ang_mom_proj_st = nint(self%M)

  end function get_ang_mom_proj_st
!
!
  function get_par_st(self)
    implicit none
    integer:: get_par_st
    type(state), intent(in) :: self
    
    get_par_st = self%parity

  end function get_par_st
!
!
  function get_spin_st(self)
    implicit none
    real(dpf):: get_spin_st
    type(state), intent(in) :: self
    
    get_spin_st = self%spin

  end function get_spin_st
!
!
  subroutine set_label_1e(self,n,l,m)
    implicit none
    type(state), intent(inout) :: self
    integer, intent(in) :: n, l, m

    self%label = make_label_1e(n,l,m)

  end subroutine set_label_1e
!
  subroutine set_label_2e(self)
    implicit none
    type(state), intent(inout) :: self

    self%label = make_label_2e(self%inum, self%parity, nint(self%M), nint(self%spin))

  end subroutine set_label_2e
!
  function get_label_st(self)
    implicit none
    character(len=5) :: get_label_st
    type(state), intent(in) :: self
    
    get_label_st = self%label

  end function get_label_st
!
!
  function get_inum_st(self)
    implicit none
    integer:: get_inum_st
    type(state), intent(in):: self
    
    get_inum_st = self%inum

  end function get_inum_st
!
!
  function get_n_majconf(self)
    implicit none
    integer:: get_n_majconf
    type(state), intent(in):: self
    
    get_n_majconf = self%n

  end function get_n_majconf!
!
  subroutine set_n_majconf(self,n)
    implicit none    
    type(state):: self
    integer, intent(in):: n

    self%n = n

  end subroutine set_n_majconf
!
!
  function get_l_majconf(self)
    implicit none
    integer:: get_l_majconf
    type(state), intent(in):: self
    
    get_l_majconf = self%l

  end function get_l_majconf!
!
  subroutine set_l_majconf(self,l)
    implicit none    
    type(state):: self
    integer, intent(in):: l

    self%l = l

  end subroutine set_l_majconf
!
!
!  function get_m_majconf(self)
!    implicit none
!    integer:: get_m_majconf
!    type(state), intent(in):: self
!    
!    get_m_majconf = self%mmajconf
!
!  end function get_m_majconf!
!
!
!  subroutine set_m_majconf(self,m)
!    implicit none    
!    type(state):: self
!    integer, intent(in):: m
!
!    self%mmajconf = m
!
!  end subroutine set_m_majconf
!!
!
  function get_nam_st(self)
    implicit none
    integer:: get_nam_st
    type(state), intent(in):: self
    
    get_nam_st = self%nam

  end function get_nam_st
!
!
  function get_nusemax_st(self)
    implicit none
    integer:: get_nusemax_st
    type(state), intent(in):: self
    
    get_nusemax_st = self%nusemax

  end function get_nusemax_st
!
!
  function get_nuse_st(self,n)
    implicit none
    type(state), intent(in) :: self
    integer, intent(in) :: n
    integer :: get_nuse_st

    get_nuse_st = self%nuse(n)

  end function get_nuse_st
!
!
  function get_na_st(self,n,i)
    implicit none
    integer:: get_na_st
    type(state), intent(in):: self
    integer, intent(in):: n
    integer, intent(in):: i   ! = 1 or 2 - which electron 'a' or 'b'
    optional:: i


    if(n .le. self%nam .and. n.ge. 1) then
       if(present(i)) then
          if(i.eq. 1) then
             get_na_st = self%na(n)
          elseif(i .eq. 2) then
             get_na_st = self%nb(n)
          else
             print*,'Error: state_class.f90, get_na_st(self,n,i): value of i is not 1 or 2: i=', i
          endif
       else
          get_na_st = self%na(n)
       endif
    else
       print*,'Error: state_class.f90, get_na_st(self,n): value of n is out of bounds: n=', n
       error stop
    endif

  end function get_na_st
!
!  
  function get_ma_st(self,n,i)
    implicit none
    integer:: get_ma_st
    type(state), intent(in):: self
    integer, intent(in):: n
    integer, intent(in):: i   ! = 1 or 2 - which electron 'a' or 'b'
    optional:: i


    if(n .le. self%nam .and. n.ge. 1) then
       if(present(i)) then
          if(i.eq. 1) then
             get_ma_st = self%ma(n)
          elseif(i .eq. 2) then
             get_ma_st = self%mb(n)
          else
             print*,'Error: state_class.f90, get_ma_st(self,n,i): value of i is not 1 or 2: i=', i
          endif
       else
          get_ma_st = self%ma(n)
       endif
    else
       print*,'Error: state_class.f90, get_ma_st(self,n): value of n is out of bounds: n=', n
       stop
    endif

  end function get_ma_st
!
  function get_CI(self,i)
    implicit none
    real(dpf):: get_CI
    type(state), intent(in):: self
    integer, intent(in):: i

    get_CI = self%CI(i)

  end function get_CI
!
!
  subroutine new_basis_st(self,n,hlike,basis_type)
    implicit none

    type(basis_state), intent(inout):: self
    integer, intent(in):: n  ! number of states
    logical, intent(in):: hlike 
    integer, intent(in) :: basis_type ! =0 if  nonrel. basis is in use, or =1 if rel.
    !

    self%basis_type = basis_type
    self%hlike = hlike

    self%Nstates = 0
    self%Nmax = n
    ! create array of n states
    if(n.ne. 0) allocate( self%b(n) )!, self%l(n) )


  end subroutine new_basis_st

  subroutine copy_basis_st(basis_l, basis_r)
    implicit none
    type(basis_state), intent(in) :: basis_r
    type(basis_state), intent(inout) :: basis_l
    integer :: i

    if(associated(basis_l%b)) call destruct_basis_st(basis_l)
    call new_basis_st(basis_l,basis_r%Nmax,basis_r%hlike,basis_r%basis_type)
    do i=1, basis_r%Nmax
      call copy_st(basis_l%b(i),basis_r%b(i))
    enddo

    basis_l%Nmax = basis_r%Nmax
    basis_l%basis_type = basis_r%basis_type
    basis_l%Mmax = basis_r%Mmax
    basis_l%Nstates = basis_r%Nstates
    basis_l%nicm = basis_r%nicm
    if(allocated(basis_r%ncore)) then
      allocate(basis_l%ncore(basis_r%nicm))
      basis_l%ncore = basis_r%ncore
    endif
    if(allocated(basis_r%mncore)) then
      allocate(basis_l%mncore(basis_r%nicm))
      basis_l%mncore = basis_r%mncore
    endif
    basis_l%en_ion = basis_r%en_ion
    basis_l%Nmax_bound = basis_r%Nmax_bound
    basis_l%Nmax_open = basis_r%Nmax_open

  end subroutine copy_basis_st
    
 subroutine destruct_basis_st(self)
    implicit none
    type(basis_state), intent(inout):: self
    integer:: n
    integer:: i
    integer:: stat

    n= self%Nmax
    do i=1,n
       !       print*,'dealocate: i=', i
       call  destruct_st(self%b(i))
    enddo
    deallocate(self%b, STAT=stat)

    if(allocated(self%ncore)) deallocate(self%ncore)
    if(allocated(self%mncore)) deallocate(self%mncore)
    if(allocated(self%br_ratio)) deallocate(self%br_ratio)


  end subroutine destruct_basis_st
!
!
  function basis_size_st(self)
    implicit none
    integer:: basis_size_st
    type(basis_state), intent(in):: self
    
    basis_size_st = self%Nmax

  end function basis_size_st
!
!
  subroutine  print_energy_basis_st(self)
    use input_data
    use target_data
    implicit none

    type(basis_state), intent(in):: self
    type(state), pointer :: stateObj
    real(dpf):: en_ion, Energy, Hartree
    integer:: nc, n, l, m, Ref_E
    real(dpf):: ioniz_en, ioniz_en_au, exit_en, two_electron_en


    Hartree = data_in%eV   ! 27.211385
    
    write(*,'("#********* Energy levels: **********")')
    print*, 'Nmax = ', self%Nmax
    if(self%hlike) then
       
       en_ion = self%en_ion
       print*, 'en_ion =', en_ion, '  NOTE: Z^2/R term might or might not be included and the internuclear distance R could be redefined from the value in the data.in'
       write(*,'(4X,"Exact ionization energy: Ref[0]: None, Ref[1]: T.E. Sharp 1970. Ref[2]: D.R. Bates et al 1953. " )')
       write(*,'(4X,"N    n    l    m   par state label    Excitation energy(eV) Ionization energy(au) Ionization energy(eV) Exact energy(eV)  Ref " )')
       
       !       write(*,'(4X,"                Excitation   Ionization  Ionization   Ionization )   " )')
       !       write(*,'(4X,"N    M    par   Energy(eV)   Energy(au)  Energy(eV) Exact Energy(eV)   Ref" )')
       
       do nc=1,self%Nmax
          stateObj => self%b(nc)
          n = get_n_majconf(stateObj)
          l = get_l_majconf(stateObj)
          m = get_ang_mom_proj(stateObj)
          exit_en =  real(Hartree*( get_energy(stateObj) -  get_energy(self%b(1))))
!          ioniz_en_au = real(( get_energy(stateObj) - en_ion)) 
          ioniz_en_au = get_elec_energy(self,nc)
          ioniz_en = Hartree*ioniz_en_au
          call get_target_energy( data_target, n,l,m, Energy, Ref_E )
          
          Energy = Hartree * Energy
          if (Ref_E == 0 ) Energy = 0d0
          write(*,'(5I5,A10,4F22.5,I5)') nc, n,l,m,get_par(stateObj), get_label(stateObj), exit_en, ioniz_en_au, ioniz_en, Energy, Ref_E
!          write(*,'(I5,F6.1,I5,3X,4(F10.5,3X),1X,I5)') nc, stateObj%M, stateObj%parity, exit_en, ioniz_en_au, ioniz_en, real(Hartree * H2Ieny(nc)), Ref_H2Ieny(nc)
       end do

    else   ! not hlike       
       en_ion = self%en_ion
       
       write(*,'(4X,"N    M   par   S   label?     two-electron energy(au) excitation energy(eV) ionization energy(au) ionization energy(eV)" )')

       do nc=1,self%Nmax
          stateObj => self%b(nc)
          exit_en =  real(Hartree*( get_energy(stateObj) -  get_energy(self%b(1))))
          ioniz_en_au = real(( get_energy(stateObj) - en_ion)) 
          ioniz_en = real(Hartree*( get_energy(stateObj) - en_ion)) 
          two_electron_en = get_energy(stateObj)
!if (data_in%Rd > 0) two_electron_en = two_electron_en - data_in%Z*data_in%Z/data_in%Rd
          write(*,'(4I5,3X,A5,4F22.5)') nc, get_ang_mom_proj(stateObj), get_par(stateObj), nint(get_spin(stateObj)), get_label(stateObj), two_electron_en, exit_en, ioniz_en_au, ioniz_en
          
!          write(101,'(F20.12)',advance='no') two_electron_en*Hartree
       end do

!!$ Mark: Neeed to remove the below. Just there for testing
!       open(302,file='Energy_R')
!       write(302,'(F12.6,3X,F12.6,3X,F12.6)') data_in%Rd,real(Hartree*get_energy(self%b(1)) ), real(en_ion*Hartree)
!      write(302,'(4F22.12)') data_in%Rd, get_energy(self%b(1)), en_ion
!      close(302)
!       stop 


    endif
    print*

  end subroutine print_energy_basis_st
  !
  ! Shellsort algorithm
  subroutine sort_by_energy_basis_st(self)
    use MPI_module
    use input_data
    implicit none

    type(basis_state), intent(inout):: self
    integer:: gap, i, j, N, jjj, iii
    type(state):: Tmp
    logical:: oldyesno, yesno
    logical :: input_order

    N = self%Nmax
    gap = N/2
    do 
       if(gap .le. 0) exit
       do i=gap+1,N
          call copy_st(Tmp,self%b(i))
          do j=i,gap+1,-gap
             if(Tmp%energy .lt. self%b(j-gap)%energy ) then
                call copy_st(self%b(j),self%b(j-gap))
             else
                exit
             endif
             call copy_st(self%b(j-gap),Tmp)
          enddo
       enddo
       if ( gap .eq. 2) then
          gap = 1
       else
          gap = gap/2!.2
       endif       
    enddo

!!$  Additional sorting to ensure that postive M is before negative M for the degenerate states
!!$ do it as default
!!$ for old runs (potl is already exists) we must call  check_potl_ordering(self) to use the same state ordering as in the potl file  - called from H12.f90
    if (data_in%good_m) then
       inquire(file='old_state_ordering', exist=oldyesno)
       if(oldyesno) then
          ! do nothing
          if (myid==0) print*,'use old state ordering'
       else
          if (myid==0) print*,'use stable state ordering'
          do i=1,N
             if(self%b(i)%M .gt. 0) then
                if(self%b(i)%M .eq. abs(self%b(i-1)%M) .and. self%b(i)%parity .eq. self%b(i-1)%parity .and. self%b(i)%spin .eq. self%b(i-1)%spin) then
                   if( abs(self%b(i)%energy - self%b(i-1)%energy) .lt. 1e-6) then
                      call copy_st(Tmp,self%b(i))
                      call copy_st(self%b(i),self%b(i-1))
                      call copy_st(self%b(i-1),Tmp)
                   endif
                endif
                
             endif
          enddo
       endif

       inquire(file='state_order.in',exist=input_order)
       if(input_order .and. .not. self%hlike) call reorder_states_2el(self)
    end if

  end subroutine sort_by_energy_basis_st

  subroutine reorder_states_2el(self)
    !Liam
    implicit none
    type(basis_state), intent(inout):: self
    type(state):: temp, temp_j
    integer :: i, j, k, nunit, N, par_in, M_in, iS_in
    real(dpf) :: S_in, E_in
    character(len=5) :: label_in
    integer :: io
    
    N = self%Nmax

    open(newunit=nunit,file='state_order.in',action='read')
    !read first 5 lines in state_order.in file - they are not important
    do i=1,5        !TODO: this is only valid if the file has only ONE H2+ state in it
      read(nunit,*)
    enddo
    do i=1, N
      
      read(nunit,'(4I4,F12.6,2X,A5)',iostat=io) j, par_in, M_in, iS_in, E_in, label_in
      if(io /= 0) then
        label_in = ''
        backspace(nunit)
        read(nunit,*) j, par_in, M_in, iS_in, E_in
      endif
      
      if(j/=i) error stop 'error: j/=i in reorder_states_2el'

      S_in = iS_in

      !print*, 'i=', i
      !print*, 'sym_in:', M_in, par_in, S_in
      !print*, 'sym:   ', self%b(i)%M, self%b(i)%parity, self%b(i)%spin
      !print*, 'label_in:', label_in

      !First try to sort using state labels

      if(trim(adjustl(label_in)) /= '?' .and. trim(adjustl(label_in)) /= '') then
      
        if(trim(adjustl(self%b(i)%label)) == trim(adjustl(label_in))) cycle
        do j=i+1,N
          if(trim(adjustl(self%b(j)%label)) == trim(adjustl(label_in))) exit 
        enddo
        if(j>N) error stop 'error: incorrect structure in reorder_states_2el (1)'

      else

        !If the i-th state has the same symmetry as the i-th state in the state_order.in file then goto next state
        if(self%b(i)%M==M_in .and. self%b(i)%parity==par_in .and. self%b(i)%spin==S_in) cycle 

        !Search through states until we find one matching the symmetry of the i-th state in the state_order.in file
        ! - we start looking at j=i+1 because the lower states have already been set correctly at this point
        do j=i+1,N
          if(self%b(j)%M==M_in .and. self%b(j)%parity==par_in .and. self%b(j)%spin==S_in) exit
        enddo

        !If the loop got to the end without finding a state of the right symmetry then the present structure
        !  calculation does not match the one used to generate the state_order.in file
        if(j>N) error stop 'error: incorrect structure in reorder_states_2el (2)'

      endif

      !Now j is the index of the state which should be in index i
      !Make room for state j in index i by shifting all states from i to j-1 up one index
      !This is necessary rather than just swapping i and j to maintain the order of the states above i

      call copy_st(temp_j,self%b(j)) !now temp_j contains the j-th state
      do k=j,i+1,-1
        call copy_st(self%b(k),self%b(k-1)) !shift state up one index
      enddo
      call copy_st(self%b(i),temp_j) !copy the original j-th state into index i

    enddo !i


    print*, ''; print*, 'STATES REORDERED TO MATCH STATE_ORDER.IN FILE'; print*, ''

  end subroutine reorder_states_2el

!
!
! this code will check the state ordering by comparing it with the potl file (if it exists)
! Aim is to correct an error in the state order related with degenerate M=+-|M| states
! no other errors are allowed
  subroutine check_potl_ordering(self,icheck_Nmax_open)
    
    use input_data    
    use MPI_module
    
    implicit none
 
    type(basis_state), intent(inout):: self
    integer, intent(inout):: icheck_Nmax_open

    real(dpf):: eV
    character(len=80):: potlfile
    logical:: ex
    character(len=12):: ench
    integer:: npotlfile
    type(input):: old_data_in

    real(dpf), dimension(:), allocatable:: energy
    integer, dimension(:), allocatable:: parity
    integer, dimension(:), allocatable:: M
    integer:: n, Nmax, Nmax_open, iprint, Mself, Mselfp1, iostat_data
    character(len=10):: text
    type(state):: Tmp

    eV = data_in%eV  ! 27.2116
    if(myid .eq. 0) print*, 'Entering  check_potl_ordering(self)'




!!$  Check if potl file is there
!!$ if no do nothing, if yes adopt the state ordering that is already in the potl file
!!$ open potl file

    if(myid .eq. 0) write(ench,'(1p,"_",e11.5)') eV*data_in%energy
    potlfile = 'potl'//ench
    potlfile = trim(adjustl(potlfile))
    npotlfile = 110
    inquire(file=potlfile,exist=ex)
    if(ex) then
       if(myid .eq. 0) print*, 'check_potl_ordering(self): potl file exists: check state ordering in potl with the current state ordering '
       open(npotlfile,file=potlfile)
       if(myid .eq. 0) print*,'Start reading potl file:',  potlfile
       iprint = 1
       !       call  readin(old_data_in,basis_1e,basis_2e,sub_basis_2e,npotlfile,iprint)
       !call  readin_skip(npotlfile,iprint)
    else
       if(myid .eq. 0) print*,'check_potl_ordering(): no potl file, can continue with present state ordering'
       return
    endif

    if(myid .eq. 0)  print*,'check_potl_ordering(): Reading energies,parity, ang.mom.proj.'
    read(npotlfile,'(5X,I5)',ADVANCE='NO',iostat=iostat_data)  Nmax
    if(myid .eq. 0)  print*, 'Nmax=', Nmax
    read(npotlfile,'(12X,I5)',END=467,iostat=iostat_data) Nmax_open
    !read(npotlfile,'(5X,I5,12X,I5)')  Nmax, Nmax_open
    if(Nmax_open .eq. 0) Nmax_open = -1    ! DF 5-10-2018 this is due to a wrong behaviour of compiler, EOR is not working
    if(myid .eq. 0) print*, 'Nmax_open=', Nmax_open
    goto 468
467 Nmax_open = -1
468 if(Nmax_open .gt. 0) then
       if(Nmax_open .ne. self%Nmax_open) then
          if(myid .eq. 0) print*,'!!!!!>>>>>   state_class.f90: check_potl_ordering(self):  (Nmax_open .ne. self%Nmax_open', Nmax_open, self%Nmax_open
          if(myid .eq. 0) print*, ' trying to fix the problem by redefining Nmax_open in Jloop.f90'
          icheck_Nmax_open = Nmax_open
          self%Nmax_open = Nmax_open
       endif
    endif

    allocate(energy(Nmax),parity(Nmax),M(Nmax))
    read(npotlfile,*) (energy(n), n=1,Nmax)
    read(npotlfile,*) (parity(n), n=1,Nmax)
    read(npotlfile,*) (M(n), n=1,Nmax)

    close(npotlfile)

!    print*
!    do n=1,10
!       print*, n, parity(n), M(n), energy(n)
!    enddo
!    print*


    if(Nmax .ne. self%Nmax) then
       if(myid .eq. 0) print*,'check_potl_ordering(): different Nmax values: Nmax, self%Nmax=', Nmax, self%Nmax
       stop
    endif

    do n=1,Nmax-1       
       if(parity(n) .ne. self%b(n)%parity) then ! errors like this are not allowed
          if(myid .eq. 0) print*,'check_potl_ordering(): parity missmatch when trying to correct the state ordering:n,parity(n),self%parity(n):', n,parity(n),self%b(n)%parity
          stop
       endif

       Mself = self%b(n)%M
       if(M(n) .ne. Mself) then  
          if(M(n) .ne. -Mself) then  ! errors like this are not allowed
             if(myid .eq. 0) print*,'check_potl_ordering(): wrong M value: n, M(n), -self%M(n)=', n, M(n), -Mself
             stop
          endif

          if(myid .eq. 0) print*,'check_potl_ordering(): correcting wrong M:  n,M(n), self%b(n)%M=', n,M(n), Mself

          ! check that the next state has -M and the same parity for both new and old ordering
          if(M(n) .ne. -M(n+1)) then
             stop 'M(n) .ne. -M(n+1)'
          endif
          Mselfp1 = self%b(n+1)%M
          if(self%b(n)%M .ne. -Mselfp1) then
             stop 'self%b(n)%M .ne. self%b(n+1)%M'
          endif
          if(parity(n) .ne. parity(n+1)) then
             stop 'parity(n) .ne. parity(n+1)'
          endif
          if(self%b(n)%parity .ne. self%b(n+1)%parity) then
             stop 'self%b(n)%parity .ne. self%b(n+1)%parity'
          endif
          !
          
          ! here we have: same parity for n and n+1 for new and old ordering; values of |M| are the same
          ! need just to swap states in new ordering
          call copy_st(Tmp,self%b(n))
          call copy_st(self%b(n),self%b(n+1))
          call copy_st(self%b(n+1),Tmp)
          
       endif
    enddo
    ! check the last state
    n = Nmax
    if(parity(n) .ne. self%b(n)%parity .or. M(n) .ne. self%b(n)%M) then ! errors like this are not allowed
       if(myid .eq. 0)  print*,'check_potl_ordering(): parity missmatch when trying to correct the state ordering:n,parity(n),self%parity(n):', n,parity(n),self%b(n)%parity
       if(myid .eq. 0)  print*,'check_potl_ordering(): M missmatch when trying to correct the state ordering:n,M(n), self%b(n)%M:', n, M(n), self%b(n)%M
       stop
    endif

    if(myid .eq. 0) print*, 'Finish check_potl_ordering()'
    return
  end subroutine check_potl_ordering
  !
  !
  subroutine calc_spectro_factors(self, bst, labot,latop)
    !
    ! Calculates the spectroscopic factors of each state for each subset of basis functions
    ! and guesses the l (in TargetState%l array) and thus n (in each state%n) of the state. 
    !
    ! JS 6/12/11
    !
    use input_data
    use MPI_module
    use sturmian_class

    implicit none

    type(basis_state), intent(inout) :: self
    type(basis_sturmian_nr), intent(in) :: bst
    integer, intent(in):: labot, latop

    type(state), pointer :: stateObj, stateObj2
    type(sturmian_nr), pointer :: sturmObj
    logical :: isDebug
    character(len=13) :: fileName
    integer :: stateNum, n,l,m,par, nFunc, i,j
    integer, dimension(:), allocatable :: naVec
    real(dpf) :: energy, CI
    real(dpf), dimension(0:latop) :: spectroVec
    real(dpf), dimension(:), allocatable :: sumVec
    real(dpf), dimension(:,:), allocatable :: spectroMat



    ! Set up the file.
    open(10, file='states.core_parts')
    write(10,'(A)',advance='no') '#state n   l   m  par  label  Energy    S(l)   sum(S)'
    do i = labot, min(latop,9)
       write(10,'(3X,A,I1,A)',advance='no') 'S(', i, ') '
    enddo
    do i = 10, latop
       write(10,'(3X,A,I2,A)',advance='no') 'S(', i, ')'
    enddo
    write(10,*)


    do stateNum = 1, basis_size_st(self)   ! Loop through all the states.

       stateObj => self%b(stateNum)   ! Pointer to the current state object.
       energy = get_energy(stateObj)   ! Energy of the state.
       m = get_ang_mom_proj(stateObj)   ! Angular momentum projection.
       par = get_par(stateObj)   ! Parity (+1 or -1) of the state.
       nFunc = get_nam_st(stateObj)   

       allocate( naVec(nFunc), spectroMat(nFunc,nFunc), sumVec(nFunc) )
       naVec = stateObj%na
       spectroMat(:,:) = bst%ortint(naVec(:),naVec(:))

       do i = 1, nFunc
          CI = get_CI(stateObj,i)
          spectroMat(i,:) = CI * spectroMat(i,:)
          spectroMat(:,i) = CI * spectroMat(:,i)
       end do

       spectroVec(:) = 0d0
       sumVec = sum( spectroMat(:,:), 1 )
       do i = 1, nFunc
          l = get_ang_mom( bst%b(naVec(i)) )
          spectroVec(l) = spectroVec(l) + sumVec(i)
       end do

       l = maxloc(spectroVec, 1) - 1   ! Largest spectroscopic factor.
       n = l + 1   ! Ensures n is greater than l.
       stateObj%l = l   ! Assigns a value of l to the state.

       do i = 1, stateNum-1   ! Search through previous states to generate n.
          stateObj2 => self%b(i)
          if ( get_l_majconf(stateObj2)==l .and. get_ang_mom_proj(stateObj2)==m ) n=n+1
       enddo
       stateObj%n = n
       !Liam commented below so inum retains definition of index within target symmetry
       !stateObj%inum = n - l ! Mark addition inum = k of Laguerre function 
       stateObj%label = make_label_1e(n,l,m)

       write(10,'(5I4,A7,F10.5,99F8.4)') stateNum, n,l,m,par, stateObj%label, energy, spectroVec(l), sum(spectroVec), spectroVec
       deallocate(naVec, spectroMat, sumVec)

    enddo ! stateNum

    close(10)


  end subroutine calc_spectro_factors
  !
  !

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Subroutine: calc_spectro_factors_group
	!Purpose: calculates the spectroscopic factors for a set of states,
	!         which are used to assign a principal klm to a given electronic
	!         state.
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
  subroutine calc_spectro_factors_group(self, bst, labot,latop, sturm_ind_list, m_list)
    !
    ! Calculates the spectroscopic factors of each state for each subset of basis functions
    ! and guesses the l (in TargetState%l array) and thus n (in each state%n) of the state. 
    !
    ! JS 6/12/11
    !
    use input_data
    use MPI_module
    use sturmian_class

    implicit none

    type(basis_state), intent(inout) :: self
    type(basis_sturmian_nr), intent(in) :: bst
    integer, intent(in):: labot, latop

    type(state), pointer :: stateObj, stateObj2
    type(sturmian_nr), pointer :: sturmObj
    logical :: isDebug
    character(len=13) :: fileName
    integer :: stateNum, n,l,m,par, nFunc, i,j
    integer, dimension(:), allocatable :: naVec
    real(dpf) :: energy, CI
    !real(dpf), dimension(0:latop) :: spectroVec
    real(dpf), dimension(:), allocatable :: spectroVec	!Will have dimension (latop +1)^2
    real(dpf), dimension(:), allocatable :: sumVec
    real(dpf), dimension(:,:), allocatable :: spectroMat

		integer, dimension(:):: sturm_ind_list   !length of nam=num_func
		integer, dimension(:):: m_list   !length of nam=num_func
		integer:: lPrev, mPrev, specInd
		integer:: lind, mind, lmind, lSearch, mSearch, maxInd
		integer:: varl, varm


    ! Set up the file.
    open(10, file='states.core_parts')
    write(10,'(A)',advance='no') '#state n   l   m  par  label  Energy    S(lm)   sum(S)'
    do i = labot, min(latop,9)
			 do j = -i, -1
          write(10,'(3X,A,I1,A,I2,A)',advance='no') 'S(', i,',', j, ') '
			 end do
			 do j = 0, i
          write(10,'(3X,A,I1,A,I1,A)',advance='no') 'S(', i,',', j, ') '
			 end do
    enddo
    do i = 10, latop
			 do j = -i, -1
          write(10,'(3X,A,I2,A,I3,A)',advance='no') 'S(', i,',', j, ')'
		   end do
			 do j = 0, i
          write(10,'(3X,A,I2,A,I2,A)',advance='no') 'S(', i,',', j, ')'
		   end do
    enddo
    write(10,*)


		!$OMP PARALLEL DO DEFAULT(SHARED) &
		!$OMP PRIVATE(stateNum,energy,m,par,nFunc,i,j,CI,spectroMat,spectroVec,naVec,sumVec,l,specInd,lPrev,mPrev,maxInd,lSearch,mSearch,lmind,lind,mind,n,stateObj,stateObj2)
    do stateNum = 1, basis_size_st(self)   ! Loop through all the states.

       stateObj => self%b(stateNum)   ! Pointer to the current state object.
       energy = get_energy(stateObj)   ! Energy of the state.
       m = get_ang_mom_proj(stateObj)   ! Angular momentum projection.
       par = get_par(stateObj)   ! Parity (+1 or -1) of the state.
       nFunc = get_nam_st(stateObj)   

       allocate( naVec(nFunc), spectroMat(nFunc,nFunc), sumVec(nFunc) )
       naVec = stateObj%na
       !spectroMat(:,:) = bst%ortint(naVec(:),naVec(:))
			 !SpectroMat should be initialised to overlap matrix of sturmian basis
			 spectroMat(:,:) = 0d0
			 do i = 1, nFunc
					do j = 1, nFunc
						 if (m_list(i) .eq. m_list(j)) then
						    spectroMat(i,j) = bst%ortint(sturm_ind_list(i),sturm_ind_list(j))
						 end if
				  end do
			 end do

			 !Spectrofactor_lm = sum_{i,j st li=lj=l mi=mj=m} c_i c_j S_ij
       do i = 1, nFunc
          CI = get_CI(stateObj,i)
          spectroMat(i,:) = CI * spectroMat(i,:)
          spectroMat(:,i) = CI * spectroMat(:,i)
       end do

			 !sum_{l=0^L} (2l+1) = (L+1)^2   = number of lm pairs
			 !Spectrovec now loops over lm, so: 1=(0,0), 2=(1,-1), 3=(1,0), 4=(1,1), 5=(2,-2), etc.
			 allocate(spectroVec((latop+1)**2))
       spectroVec(:) = 0d0
       sumVec = sum( spectroMat(:,:), 1 )
			 specInd = 0
			 lPrev = -1
			 mPrev = -1000
       !do i = 1, nFunc
       !   l = get_ang_mom( bst%b(sturm_ind_list(naVec(i))) )
			 ! 	m = m_list(naVec(i))
			 ! 	if (l .ne. lPrev) then
			 !      specInd = specInd + 1
			 ! 		 lPrev = l
			 ! 		 mPrev = m
			 !   else
			 ! 		 if (m .ne. mPrev) then
			 ! 		    specInd = specInd + 1
			 ! 				mPrev = m
			 ! 		 end if
			 !   end if

       !   spectroVec(specInd) = spectroVec(specInd) + sumVec(i)
       !end do

			 specInd = 0
			 do l = labot, latop
					do m = -l, l
						 specInd = specInd + 1
						 do i = 1, nFunc
						    do j = 1, nFunc
								   if ((get_ang_mom(bst%b(sturm_ind_list(naVec(i)))) .eq. l) &
									     .and. (m_list(i) .eq. m)) then
								      if ((get_ang_mom(bst%b(sturm_ind_list(naVec(j)))) .eq. l) &
									       .and. (m_list(j) .eq. m)) then
							           spectroVec(specInd) = spectroVec(specInd) + spectroMat(i,j)
										  end if
							     end if
						    end do
						 end do
				  end do
			 end do
 
       !l = maxloc(spectroVec, 1) - 1   ! Largest spectroscopic factor.
       maxInd = maxloc(spectroVec, 1)    ! Largest spectroscopic factor.
			 !Get lm value corresponding to this index 
       lmind = 1
			 lSearch = 0
			 mSearch = 0
			 do lind = labot, latop
					do mind = -lind, lind
						 if (lmind .eq. maxInd) then
								lSearch = lind
								mSearch = mind
						 end if
						 lmind = lmind + 1
				  end do
			 end do
			 l = lSearch
			 m = mSearch 

       n = l + 1   ! Ensures n is greater than l.
       stateObj%l = l   ! Assigns a value of l to the state.
			 stateObj%M = m   ! Assigns a value of m to the state.

       do i = 1, stateNum-1   ! Search through previous states to generate n.
          stateObj2 => self%b(i)
          if ( get_l_majconf(stateObj2)==l .and. get_ang_mom_proj(stateObj2)==m ) n=n+1
       enddo
       stateObj%n = n
       !Liam commented below so inum retains definition of index within target symmetry
       !stateObj%inum = n - l ! Mark addition inum = k of Laguerre function 
			 if (data_in%good_m) then
          stateObj%label = make_label_1e(n,l,m)
			 else
					stateObj%label = "   -   "
			 end if

       !write(10,'(5I4,A7,F10.5,99F8.4)') stateNum, n,l,m,par, stateObj%label, energy, spectroVec(l), sum(spectroVec), spectroVec
       write(10,'(5I4,A7,F10.5,99F9.5)') stateNum, n,l,m,par, stateObj%label, energy, spectroVec(maxInd), sum(spectroVec), spectroVec
       deallocate(naVec, spectroMat, sumVec, spectroVec)

    enddo ! stateNum
		!$OMP END PARALLEL DO

    close(10)


  end subroutine calc_spectro_factors_group




  subroutine spectroscopic_2e(TargetStates2el)
    !
    implicit none

    type(basis_state), intent(in) :: TargetStates2el

    type(state), pointer :: state2e
    integer :: j2e


    open(11, file='states.core_parts', action='write', position='append') !Cray compiler doesn't like access='append'
    write(11,'(/A)') '# TWO-ELECTRON STATES'
    write(11,'(4A4,A12)') '#  N', 'par','M','S', 'Energy (au)'
    
    do j2e = 1, TargetStates2el%Nmax
       state2e => TargetStates2el%b(j2e)
       write(11,'(4I4,F12.6,2X,A)') j2e, state2e%parity,nint(state2e%M),nint(state2e%spin), state2e%energy, state2e%label
    end do

    close(11)

  end subroutine spectroscopic_2e
  !
  !
  subroutine spectroscopic_2e_old(self,TargetStates, bst, lmax, nConfigs,ind1Vec,ind2Vec,CIVec,bMat)
    use sturmian_class
    !
    implicit none

    type(state), intent(in) :: self
    type(basis_state), intent(in) :: TargetStates
    type(basis_sturmian_nr), intent(in) :: bst
    integer, intent(in) :: lmax, nConfigs
    integer, dimension(nConfigs), intent(in) :: ind1Vec, ind2Vec
    real(dpf), dimension(nConfigs), intent(in) :: CIVec
    real(dpf), dimension(nConfigs,nConfigs), intent(inout) :: bMat

    type(state), pointer :: orb
    type(sturmian_nr), pointer :: spf
    integer :: indf,indi, indorb, indspf,nspfs, indff,indii, l
    real(dpf) :: sumf,sumff, CI
    real(dpf), dimension(0:lmax) :: SF1Vec,SF2Vec
    real(dpf), dimension(:,:), allocatable :: bbMat


    SF1Vec(:) = 0d0
    SF2Vec(:) = 0d0
    open(11,file='states.core_parts',position='append')
    write(11,'(/4A4,3A10)',advance='no') '#', 'Par','M','S', 'E (au)','S(l)','sum(S)'
    do l = 0, lmax
       write(11,'(A8,I1,A1)',advance='no') 'S(',l,')'
    end do
    write(11,*)

    do indf = 1, nConfigs
       sumf = 0d0
       do indi = 1, nConfigs
          sumf = sumf + CIVec(indf)*CIVec(indi)*bMat(indf,indi)
       end do

       ! Electron 1.
       indorb = ind1Vec(indf); orb => TargetStates%b(indorb)
       nspfs = orb%nam
       allocate( bbMat(nspfs,nspfs) )
       bbMat(:,:) = bst%ortint( orb%na(:), orb%na(:) )
       
       do indff = 1, nspfs
          sumff = 0d0
          do indii = 1, nspfs
             sumff = sumff + orb%CI(indff)*orb%CI(indii)*bbMat(indff,indii)
          end do

          indspf = orb%na(indff); spf => bst%b(indspf)
          l = get_ang_mom(spf)
          SF1Vec(l) = SF1Vec(l) + sumf*sumff
       end do
       deallocate(bbMat)

       ! Electron 2.
       indorb = ind2Vec(indf); orb => TargetStates%b(indorb)
       nspfs = orb%nam
       allocate( bbMat(nspfs,nspfs) )
       bbMat(:,:) = bst%ortint( orb%na(:), orb%na(:) )

       do indff = 1, nspfs
          sumff = 0d0
          do indii = 1, nspfs
             sumff = sumff + orb%CI(indff)*orb%CI(indii)*bbMat(indff,indii)
          end do

          indspf = orb%na(indff); spf => bst%b(indspf)
          l = get_ang_mom(spf)
          SF2Vec(l) = SF2Vec(l) + sumf*sumff
       end do
       deallocate(bbMat)

    end do

    l = maxloc(SF1Vec(:),1) - 1
    write(11,'(4I4,99F10.5)') self%inum, self%parity,nint(self%M),nint(self%spin), self%energy, SF1Vec(l),sum(SF1Vec(:)),SF1Vec(:)
    l = maxloc(SF2Vec(:),1) - 1
    write(11,'(26X,99F10.5)') SF2Vec(l), sum(SF2Vec(:)), SF2Vec(:)
    close(11)


  end subroutine spectroscopic_2e_old
  !
  !

  subroutine spectro_core(TargetStates2el, TargetStates)
    use ovlpste1me
    use input_data
    use MPI_module
    implicit none

    type(basis_state), intent(inout) :: TargetStates2el
    type(basis_state), intent(in) :: TargetStates
    type(state), pointer :: st
    integer :: i, j, jcore, c, cc, ca, cb, cca, ccb, nicm, Nmax
    integer, allocatable :: ncore(:)
    real(dpf) :: B, tol
    logical :: reduce_pseudo
    
    reduce_pseudo = (data_in%SF_tol < 1.0d0)

    if(.not. reduce_pseudo) return !temp because there is an error with SF_core bounds

    nicm = TargetStates2el%nicm
    allocate(ncore(nicm))
    ncore = TargetStates2el%ncore
 
    do i=1, nicm
      ncore(i) = i
    enddo
    
    !print*, 'ncore:', ncore
    !print*, 'label:', TargetStates%b(ncore)%label

    if(myid==0) open(unit=123,file='spectro_core.out',action='write',status='replace')

    if(myid==0) write(123,'("State    sum     ", 1000(A4,4X))') (TargetStates%b(ncore(i))%label, i=1, nicm)
    do i=1, TargetStates2el%Nmax
      st => TargetStates2el%b(i)

      allocate(st%SF_core(nicm))

      st%SF_core = 0.0d0  
      do c=1, st%nam
            ca = st%na(c)
            cb = st%nb(c)
            j = min(ca,cb)
            do cc=1, st%nam
              cca = st%na(cc)
              ccb = st%nb(cc)
              B = ovlpst(ca,cca) * ovlpst(cb,ccb)
              st%SF_core(j) = st%SF_core(j) + st%CI(c) * st%CI(cc) * B
            enddo !cc
        enddo !c

        if(myid==0) write(123,'(A5,3X,1000(F6.3,2X))') st%label, sum(st%SF_core(:)), st%SF_core(:)

    enddo !i (state)

    if(myid==0) close(123)

    if(reduce_pseudo) then

      tol = data_in%SF_tol

      jcore = 0
      do i=1, nicm
        st => TargetStates%b(i)
        if(st%n == 1 .and. st%l == 0 .and. st%m == 0.0d0) jcore = i
      enddo
      if(jcore == 0) error stop '*** ERROR in spectro_core: couldn''t find 1sS core orbital'

      Nmax = TargetStates2el%Nmax
      i=1
    
      do while (i <= Nmax)
  
        st => TargetStates2el%b(i)
 
        !print*, 'N, LABEL, M, SF:', i, st%label, st%m, abs(st%SF_core(jcore))

        if(st%energy > TargetStates2el%en_ion .and. abs(st%SF_core(jcore)) < tol) then
          !print*, '   ..remove'
          do j=i+1, Nmax
            call copy_st(TargetStates2el%b(j-1),TargetStates2el%b(j))
          enddo
          Nmax = Nmax - 1
        else
          i = i + 1
        endif
  
      enddo
  
     TargetStates2el%Nmax = Nmax
  
     if(myid==0) then
       call print_energy_basis_st(TargetStates2el)
       call oscstr_2e()
     endif

   endif

  end subroutine spectro_core
  !
  function make_label_1e(n,l,m)
    use input_data
    !
    ! Converts a set of quantum numbers (n,l,m,parity) into a 5-character label.
    !
    ! First is n; second is l (in lowercase s,p,d,f notation);
    ! third is m (in uppercase S,P,D,F to symbolise sigma,pi,delta,phi); and
    ! fourth is parity (g = gerade = even = +1, u = ungerade = odd = -1).
    !
    ! For example, 1sSg is the ground state. The next energetic is 2pSu.
    !
    ! Note that m-degenerate states are labelled the same.
    !
    ! JS 6/12/11
    !
    implicit none

    character(len=5) :: make_label_1e
    integer, intent(in) :: n, l, m

    make_label_1e = ''
    make_label_1e(1:2) = char(n+48)   ! 48->0, 49->1, etc. Google "ASCII table".
    
    select case(l)
       case(0)
          make_label_1e(2:3) = 's'
       case(1)
          make_label_1e(2:3) = 'p'
       case(2)
          make_label_1e(2:3) = 'd'
       case default
          make_label_1e(2:3) = char(l+99)   ! 102->f, 103->g, etc. Good till l>12 .
    end select

    select case( abs(m) )   ! Degenerate +/- m states are treated as the same.
       case(0)
          make_label_1e(3:4) = 'S'
       case(1)
          make_label_1e(3:4) = 'P'
       case(2)
          make_label_1e(3:4) = 'D'
       case default
          make_label_1e(3:4) = char(abs(m)+67)   ! 70->F, 71->G, etc. Till m>12.
    end select

    if(data_in%good_parity) then  !Liam added if statement: don't include u/g in labels for heteronuclear diatomics
      if ( (-1)**l == 1 ) then
         make_label_1e(4:) = 'g'
      elseif ( (-1)**l == -1 ) then
         make_label_1e(4:) = 'u'
      endif
    endif !good parity


  end function make_label_1e

  function make_label_2e(j, Par, M, Spin)
    use input_data
    !
    ! Convert set of quantum numbers Par, M, Spin, and state number within given symmetry into a 5-character label
    !
    ! Liam 15/12/2020
    !
    implicit none
    
    character(len=5) :: make_label_2e
    integer, intent(in) :: j, Par,M,Spin
    character(len=30), parameter :: symchars = 'u gSPDFGHIKLMNOQRTUV??????????'
       
    write(make_label_2e,'(I1,"_",I1,A1,A1)') j, 2*Spin+1, symchars(abs(M)+4:abs(M)+4),symchars(Par+2:Par+2)

    select case(trim(adjustl(data_in%target)))
      case ('H2')
        call add_state_letter_H2(j,make_label_2e)
      case ('HeH+')
        call add_state_letter_HeHplus(j,make_label_2e)
    end select


  end function make_label_2e


  subroutine add_state_letter_H2(j,label)
    !
    ! Determines Latin character to add to symmetry label for the H2 molecule
    !
    ! Liam 15/12/2020
    !
    implicit none
    integer, intent(in) :: j
    integer :: jmax, col
    character(len=5), intent(inout) :: label
    character(len=3) :: sym
    character(len=2), dimension(6,12) :: letters
    character(len=3), dimension(12) :: syms

    sym = label(3:5)
    
    jmax = size(letters, 1)

    if(j>jmax) then
      if(j<10) then
        write(label,'(I1,"_",A3)') j, sym
      else
        write(label,'(I2,A3)') j, sym
      endif
      return
    endif

    syms    =          (/   '1Sg', '1Su',  '1Pg',  '1Pu',  '1Dg',  '1Du',  '3Sg',  '3Su',  '3Pg',  '3Pu',  '3Dg',  '3Du' /)
    letters = transpose( &
            & reshape( (/   ' X',  ' B' ,  ' I',   ' C',   ' J',   '4f',   ' a',   ' b',   ' i',   ' c',   ' j',   '4f', &  !j=1
                          & 'EF',  "B'" ,  ' R',   ' D',   ' S',   '2_',   ' h',   ' e',   ' r',   ' d',   ' s',   '2_', &  !j=2
                          & 'GK',  'B"',   '3_',   ' V',   '3_',   '3_',   ' g',   ' f',   '3_',   ' k',   '3_',   '3_', &  !j=3
                          & ' H',  '4f' ,  '4_',   "D'",   '4_',   '4_',   ' p',   '4f',   '4_',   '4f',   '4_',   '4_', &  !j=4
                          & ' P',  '5_' ,  '5_',   '5_',   '5_',   '5_',   '4s',   '5_',   '5_',   '5_',   '5_',   '5_', &  !j=5
                          & ' O',  '6_' ,  '6_',   '6_',   '6_',   '6_',   '6_',   '6_',   '6_',   '6_',   '6_',   '6_'  &  !j=6
                          &/), shape(transpose(letters)) ) )


    !col = findloc(syms, sym, dim=1)
    !Above annoying doesn't work - I think the FINDLOC function is not implemented in IFORT?
    !Have to do this instead:
    select case(sym)
      case ('1Sg')
        col = 1
      case ('1Su')
        col = 2
      case ('1Pg')
        col = 3
      case ('1Pu')
        col = 4
      case ('1Dg')
        col = 5
      case ('1Du')
        col = 6
      case ('3Sg')
        col = 7
      case ('3Su')
        col = 8
      case ('3Pg')
        col = 9
      case ('3Pu')
        col = 10
      case ('3Dg')
        col = 11
      case ('3Du')
        col = 12
      case default
        col = 0
    end select

    if(col == 0) then ! Haven't coded letters for this symmetry
      !label = '__'//sym
      if(j<10) then
        write(label,'(I1,"_",A3)') j, sym
      else
        write(label,'(I2,A3)') j, sym
      endif
    else
      label = letters(j,col)//sym
    endif

  end subroutine add_state_letter_H2
  
  subroutine add_state_letter_HeHplus(j,label)
    !
    ! Determines Latin character to add to symmetry label for the HeH+ molecule
    !
    ! Liam 15/12/2020
    !
    implicit none
    integer, intent(in) :: j
    integer :: jmax, col
    character(len=5), intent(inout) :: label
    character(len=3) :: sym
    character(len=2), dimension(6,6) :: letters
    character(len=3), dimension(6) :: syms

    sym = label(3:5)
    
    jmax = size(letters, 1)

    if(j>jmax) then
      if(j<10) then
        write(label,'(I1,"_",A3)') j, sym
      else
        write(label,'(I2,A3)') j, sym
      endif
      return
    endif

    syms    =          (/   '1S ',  '1P ',  '1D ',  '3S ',  '3P ',   '3D ' /)
    !letters = transpose( &
    !        & reshape( (/   ' X',  '2p' ,  '3d',   '2s',   '2p',   '3d', &  !j=1
    !                      & '2s',  '3d' ,  '2_',   '2p',   '3p',   '2_', &  !j=2
    !                      & '2p',  '3p',   '3_',   '3s',   '3d',   '3_', &  !j=3
    !                      & '3s',  '4_' ,  '4_',   '3p',   '4_',   '4_', &  !j=4
    !                      & '3d',  '5_' ,  '5_',   '3d',   '5_',   '5_', &  !j=5
    !                      & '3p',  '6_' ,  '6_',   '6_',   '6_',   '6_'  &  !j=6
    !                      &/), shape(transpose(letters)) ) )
    
    letters = transpose( &
            & reshape( (/   ' X',  ' C' ,  ' I',   ' a',   ' c',   ' i', &  !j=1
                          & ' A',  ' F' ,  '2_',   ' b',   ' e',   '2_', &  !j=2
                          & ' B',  ' H',   '3_',   ' d',   ' h',   '3_', &  !j=3
                          & ' D',  '4_' ,  '4_',   ' f',   '4_',   '4_', &  !j=4
                          & ' E',  '5_' ,  '5_',   ' g',   '5_',   '5_', &  !j=5
                          & ' G',  '6_' ,  '6_',   '6_',   '6_',   '6_'  &  !j=6
                          &/), shape(transpose(letters)) ) )


    !col = findloc(syms, sym, dim=1)
    !Above annoying doesn't work - I think the FINDLOC function is not implemented in IFORT?
    !Have to do this instead:
    select case(sym)
      case ('1S ')
        col = 1
      case ('1P ')
        col = 2
      case ('1D ')
        col = 3
      case ('3S ')
        col = 4
      case ('3P ')
        col = 5
      case ('3D ')
        col = 6
      case default
        col = 0
    end select

    if(col == 0) then ! Haven't coded letters for this symmetry
      !label = '__'//sym
      if(j<10) then
        write(label,'(I1,"_",A3)') j, sym
      else
        write(label,'(I2,A3)') j, sym
      endif
    else
      label = letters(j,col)//sym
    endif

  end subroutine add_state_letter_HeHplus


  subroutine convert_from_oid_to_st(oid, bst, states)
    !
    use grid_radial
    use ovlpste1me
    use spheroidal_class
    use sturmian_class
    implicit none

    type(spheroidal_basis), intent(in) :: oid
    type(basis_sturmian_nr), intent(out) :: bst
    type(basis_state), intent(out) :: states

    type(spheroidal_fn), pointer :: oidFn
    type(spheroidal_fn) :: oid1,oid2
    integer :: nBst, jFn,nFn, jTerm,nTerm,nPrev, l,m,i1,i2
    integer, dimension(:), allocatable :: indVec, mVec
    real(dpf) :: D
    real(dpf), dimension(:), allocatable :: CIVec
    real(dpf), dimension(:), pointer :: oldVec
    real(dpf), dimension(1:grid%nr) :: newVec


    ! Create the state orbitals, which will point to various s.p. sturmians.
    nFn = get_num_fns(oid)
    call new_basis(states, nFn, .true., 2)

    ! Create the basis of single particle "sturmian" functions.
    nBst = 0
    do jFn = 1, nFn
       oidFn => get_spheroidal_fn(oid,jFn)
       nBst = nBst + get_num_terms(oidFn)
    end do
    call new_basis(bst,nBst)

    nBst = 0; nTerm = 0
    do jFn = 1, nFn   ! Loop through each spheroidal function.
       oidFn => get_spheroidal_fn(oid,jFn)
       nPrev = nTerm; nTerm = get_num_terms(oidFn)
       m = get_m(oidFn)
       i1 = get_rad_min(oidFn); i2 = get_rad_max(oidFn)
       oldVec => get_rad_vec(oidFn)

       if (nTerm /= nPrev) then
          if ( allocated(indVec) ) deallocate(indVec, mVec, CIVec)
          allocate( indVec(1:nTerm), mVec(1:nTerm), CIVec(1:nTerm) )
          CIVec(:) = 1d0
       end if
       mVec(:) = m

       ! Each term in the angular expansion will become an individual sturmian.
       do jTerm = 1, nTerm
          nBst = nBst + 1   ! Running tally of the number of sturmian functions.
          indVec(jTerm) = nBst

          l = get_term_l(oidFn,jTerm); D = get_term_D(oidFn,jTerm)
          newVec(:) = 0d0; newVec(i1:i2) = D*oldVec(i1:i2)
          call init_function(bst%b(nBst), l,m,jFn,0d0, i1,i2,newVec,grid%nr)
       end do

       ! The spheroidal function is now a state with links to several sturmians.
       call construct_st(states%b(jFn),.true., dble(m),(-1)**get_lambda(oidFn),0.5d0,0d0,jFn, nTerm,CIVec,indVec,mVec)
    end do

    ! Calculate the sturmian overlap matrix from scratch.
    do jFn = 1, nBst
       call convert_from_sturm_to_oid( bst%b(jFn), oid1 )
       bst%ortint(jFn,jFn) = oid_overlap(oid1,oid1)

       do jTerm = jFn+1, nBst
          call convert_from_sturm_to_oid( bst%b(jTerm), oid2 )
          bst%ortint(jFn,jTerm) = oid_overlap(oid1,oid2)
          bst%ortint(jTerm,jFn) = bst%ortint(jFn,jTerm)
       end do
    end do

    write(*,'(/3(A,I4)/)') 'convert_from_oid_to_st: converted ', nFn, ' orbitals from spheroidal functions to target states, with ', nBst, ' underlying single particle sturmians.'
    

  end subroutine convert_from_oid_to_st


!!$
!!$
!!$  state_l,state_r  are one-electron target states !>> to do make it work with two-electron states too
  function ovlp_st(state_l,state_r,bst)

    use sturmian_class
    use ovlpste1me
    use input_data

    implicit none

    real(dpf):: ovlp_st

    type(state), intent(in):: state_l, state_r
    type(basis_sturmian_nr), intent(in) :: bst

    real(dpf), dimension(:,:), pointer:: ortint
    real(dpf):: tmp_l, tmp_r, tmpsum, ovlpsp1, ovlpsp2
    integer:: i, j, m, n, ni, nm,  n1l, n2l, n1r, n2r, li, lm

    
    ovlp_st = 0d0


    if(nint(2*state_l%m) .ne. nint(2*state_r%m)) then
!       print*, 'm: l:', state_l%m
!       print*, '   r:', state_r%m
       return
    endif
    if(data_in%good_parity .and. (state_l%parity .ne. state_r%parity)) then
!       print*, 'par: l:',  state_l%parity
!       print*, '     r:',  state_r%parity
       return
    endif
    if(nint(2*state_l%spin) .ne. nint(2*state_r%spin)) then
!       print*, 'spin: l:', state_l%spin
!       print*, '      r:', state_r%spin
       return
    endif

    if( (state_r%hlike.and..not.state_l%hlike) .or. (.not.state_r%hlike.and.state_l%hlike) ) then
       print*,'state_class.f90:  ovlp_st(): state_r%hlike .ne. statel%hlike'
       print*,'stop'
       stop
    endif

!    if(state_r%hlike) then
       ortint => bst%ortint
!    endif


    tmpsum = 0d0

    if(state_r%hlike) then ! one-electron states
        do i=1,state_l%nam
          ni = state_l%na(i)   ! index to a Sturmian type function with fixe dvalue of angular moemntum
!          print*, 'ni=', ni
          li = get_ang_mom(bst%b(ni))  
          tmp_l = get_CI(state_l,i)
          do m=1,state_r%nam
             nm = state_r%na(m)  ! index to a Sturmian type function with fixe dvalue of angular moemntum
!             print*, 'nm=', nm
             lm = get_ang_mom(bst%b(nm))
             if ((data_in%calculation_type == 0 .or. data_in%calculation_type == 1).and. li.ne.lm) cycle   ! Spherical l'=l.
             ovlpsp1 = ortint(ni,nm)
             if(ovlpsp1 .eq. 0d0) cycle             
             tmp_r = get_CI(state_r,m)
             tmpsum = tmpsum + tmp_l * tmp_r * ovlpsp1
          enddo
       enddo
    else  ! two-electron states
       do i=1,state_l%nam
          n1l = state_l%na(i)
          n2l = state_l%nb(i)
          tmp_l = get_CI(state_l,i)
          do j=1,state_r%nam
             n1r = state_r%na(j)
             n2r = state_r%nb(j)
             tmp_r = get_CI(state_r,j)
             
             tmpsum = tmpsum + ortint(n1l,n1r)*ortint(n2l,n2r)*tmp_l*tmp_r
             !    tmpsum = tmpsum + ovlpst(n1l,n1r)*ovlpst(n2l,n2r)*tmp_l*tmp_r
             !  print'(">>>",4E15.6)', ovlpst(n1l,n1r),ovlpst(n2l,n2r),tmp_l,tmp_r
!             print'(2i3,2(2i5,F15.5))', i,j, n1l, n1r,ortint(n1l,n1r), n2l, n2r,ortint(n2l,n2r)
             
          enddo
       enddo
       
    endif
    ovlp_st = tmpsum
    
  end function ovlp_st


!Function: ovlp_st_group
!Purpose: calculates the overlap of two states expanded in a basis with sturmian parts in bst
!         this version allows states to have different m values in their expansion. Stored in
!         m_list_nr, index of radial functions in basis stored in sturm_ind_list
!!$  state_l,state_r  are one-electron target states !>> to do make it work with two-electron states too
  function ovlp_st_group(state_l,state_r,bst, m_list_nr, sturm_ind_list)

    use sturmian_class
    use ovlpste1me
    use input_data

    implicit none

    real(dpf):: ovlp_st_group

    type(state), intent(in):: state_l, state_r
    type(basis_sturmian_nr), intent(in) :: bst

    real(dpf), dimension(:,:), pointer:: ortint
    real(dpf):: tmp_l, tmp_r, tmpsum, ovlpsp1, ovlpsp2
    integer:: i, j, m, n, ni, nm,  n1l, n2l, n1r, n2r, li, lm

		integer, dimension(:):: m_list_nr, sturm_ind_list
		integer:: m1, m2  !m values in expansion
		integer:: indi, indj
    
    ovlp_st_group = 0d0


    if(nint(2*state_l%m) .ne. nint(2*state_r%m)) then
!       print*, 'm: l:', state_l%m
!       print*, '   r:', state_r%m
       return
    endif
    if(data_in%good_parity .and. (state_l%parity .ne. state_r%parity)) then
!       print*, 'par: l:',  state_l%parity
!       print*, '     r:',  state_r%parity
       return
    endif
    if(nint(2*state_l%spin) .ne. nint(2*state_r%spin)) then
!       print*, 'spin: l:', state_l%spin
!       print*, '      r:', state_r%spin
       return
    endif

    if( (state_r%hlike.and..not.state_l%hlike) .or. (.not.state_r%hlike.and.state_l%hlike) ) then
       print*,'state_class.f90:  ovlp_st(): state_r%hlike .ne. statel%hlike'
       print*,'stop'
       stop
    endif

!    if(state_r%hlike) then
       ortint => bst%ortint
!    endif


    tmpsum = 0d0

    if(state_r%hlike) then ! one-electron states
        do i=1,state_l%nam
          indi = state_l%na(i)  !Index to sturmian functions resolved in klm 
          ni = sturm_ind_list(indi)  !phi_{kl} for given index
					m1 = m_list_nr(indi)

!         print*, 'ni=', ni
          li = get_ang_mom(bst%b(ni))  
          tmp_l = get_CI(state_l,i)
          do m=1,state_r%nam
             indj = state_r%na(m)  ! index to a Sturmian type function resolved in klm 
						 nm = sturm_ind_list(indj)
						 m2 = m_list_nr(indj)

!            print*, 'nm=', nm
             lm = get_ang_mom(bst%b(nm))
             if ((data_in%calculation_type == 0 .or. data_in%calculation_type == 1).and. li.ne.lm) cycle   ! Spherical l'=l.
						 if ((.not. data_in%good_m) .and. (m1 .ne. m2)) cycle   !different m orthogonal
             ovlpsp1 = ortint(ni,nm)
             if(ovlpsp1 .eq. 0d0) cycle             
             tmp_r = get_CI(state_r,m)
             tmpsum = tmpsum + tmp_l * tmp_r * ovlpsp1
          enddo
       enddo
    else  ! two-electron states
			 if (.not. data_in%good_m) then
					print*, "ERROR: ovlp_st_group not yet implemented for 2e states in H3 mode"
					stop
			 end if
  
       do i=1,state_l%nam
          n1l = state_l%na(i)
          n2l = state_l%nb(i)
          tmp_l = get_CI(state_l,i)
          do j=1,state_r%nam
             n1r = state_r%na(j)
             n2r = state_r%nb(j)
             tmp_r = get_CI(state_r,j)
             
             tmpsum = tmpsum + ortint(n1l,n1r)*ortint(n2l,n2r)*tmp_l*tmp_r
             !    tmpsum = tmpsum + ovlpst(n1l,n1r)*ovlpst(n2l,n2r)*tmp_l*tmp_r
             !  print'(">>>",4E15.6)', ovlpst(n1l,n1r),ovlpst(n2l,n2r),tmp_l,tmp_r
!             print'(2i3,2(2i5,F15.5))', i,j, n1l, n1r,ortint(n1l,n1r), n2l, n2r,ortint(n2l,n2r)
             
          enddo
       enddo
       
    endif
    ovlp_st_group = tmpsum
    
  end function ovlp_st_group







!
!!$  this is a function that calculate matrix element for one-electron Hamiltonian (H2+)
!!$  state_l,state_r  are one-electron target states
!!$  it relies on special form of underlying Laguerre basis (2l+1) 
  function H1el_st(state_l,state_r,bst)

    use sturmian_class

    implicit none

    real(dpf):: H1el_st
    type(state), intent(in):: state_l, state_r  ! it is one-electron states
    type(basis_sturmian_nr), intent(in) :: bst

    real(dpf), dimension(:,:), pointer:: ortint
    real(dpf):: tmp_l, tmp_r, tmpsum, result
    integer:: ma, i, m, n, ni, nm


    ortint => bst%ortint
    
    H1el_st = 0d0

    if(state_l%m .ne. state_r%m) return
    if(state_l%parity .ne. state_r%parity) return
    
    ma = state_l%m

    if(state_l%hlike) then
    else
       print*,'>>> state_class.f90: H1el_st(...)'
       print*,'Two-electron states have not been coded yet.'
       print*, 'stop'
       stop
    endif

    tmpsum = 0d0
    do i=1,state_l%nam
       ni = state_l%na(i)
       tmp_l = get_CI(state_l,i)
       do m=1,state_r%nam
          nm = state_r%na(m)
          tmp_r = get_CI(state_r,m)

          call Hlagorb(bst,ni,nm,ma,result)

          tmpsum = tmpsum +  result * tmp_l * tmp_r
!          print'(4i5,3E15.5)', i,ni,m,nm, tmp_l, tmp_r,result
       enddo
    enddo
    
     H1el_st = tmpsum
    
   end function H1el_st

 !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
!!$ write a one-electron state in a file = State_1el
   subroutine write_1elState(bst_nr, state1el)
     !
     use grid_radial
     use sturmian_class

     implicit none

     type(basis_sturmian_nr), intent(in) :: bst_nr
     type(state), intent(in) :: state1el
     
     integer:: M, ipar, nam
     integer:: n, in, l, i1, i2, i

     type(sturmian_nr), pointer :: spFunction
     character(len=16) :: filename
     real(dpf), dimension(:), pointer :: rVec,spVec
     real(dpf):: CI, spectr_factor, sum_spectr_factor, tmp, tmpsum

     rVec => grid%gridr

     nam = get_nam(state1el)

     filename = 'State_1el'

     open(11,file=filename)

     write(11,*) '# EXPANSION OF ONE-ELECTRON TARGET STATE'
     write(11,*) '# IN TERMS OF ONE-ELECTRON ORBITALS WITH DEFINITE ANGULAR MOMENTUM'     
     M = get_ang_mom_proj(state1el)
     ipar = get_par(state1el)
     write(11,'("# M, ipar",2i5)') M, ipar

     write(11,'("# nam=",i5,"  CI:")',ADVANCE='NO') nam
     do n=1,nam
        CI = get_CI(state1el,n)
        write(11,'(1P,E15.5)',ADVANCE='NO') 
     enddo
     write(11,*)

     write(11,'("# r            (n l spectr.-factor): ")')
        write(11,'("# ")',ADVANCE='NO')
     sum_spectr_factor = 0d0
     do n=1,nam
        CI = get_CI(state1el,n)
        in = get_na(state1el,n)
        spFunction => bst_nr%b(in)
        l = get_ang_mom(spFunction)
        spVec => fpointer(spFunction) 
        i1 = get_minf(spFunction)
        i2 = get_maxf(spFunction)
        spectr_factor = CI*CI*SUM(spVec(i1:i2)*spVec(i1:i2)*grid%weight(i1:i2))
        write(11,'(2i5,1P,E15.5)',ADVANCE='NO') in,l,spectr_factor
        sum_spectr_factor = sum_spectr_factor + spectr_factor
     enddo
     tmp = ovlp_st(state1el,state1el,bst_nr)
     write(11,'("   sum:",1P,2E15.5)') sum_spectr_factor, tmp

!     write(*,*) ( bst_nr%b(get_na(state1el,n)) )
     
     do i=1,grid%nr
        write(11,'(1P,E15.5,2X)',ADVANCE='NO') grid%gridr(i)
        tmpsum = 0d0
        do n=1,nam
           in = get_na(state1el,n)
           spFunction => bst_nr%b(in)
           l = get_ang_mom(spFunction)
           CI = get_CI(state1el,n)
           tmp = CI*value(spFunction,i) 
           tmpsum = tmpsum + tmp
           write(11,'(1P,E15.5,12X)',ADVANCE='NO') tmp
        enddo
        write(11,*)
        if(tmpsum .eq. 0 .and. nam .gt. 1) exit        

     enddo

     write(*,*) 'Made ', filename
     close(11)

   end subroutine write_1elState

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

   subroutine write_configs(bst_nr, TargetStates, TargetStates2el, rearranged)
     !
     use grid_radial
     use sturmian_class
     implicit none

     type(basis_sturmian_nr), intent(in) :: bst_nr
     type(basis_state), intent(in) :: TargetStates, TargetStates2el
     logical, intent(in) :: rearranged

     type(sturmian_nr), pointer :: spFunction
     type(state), pointer :: orbital, state2e
     character(len=16) :: filename
     integer :: ind1,ind2,ind, j,jOrb,jCore,jState
     real(dpf), dimension(:), pointer :: rVec,spVec


     rVec => grid%gridr

     if (rearranged) then
        filename = 'configs_rearr'
     else
        filename = 'configs_original'
     end if
     open(11,file=filename)

     write(11,*) '# CI EXPANSIONS OF TWO-ELECTRON TARGET STATES'
     write(11,*) '# IN TERMS OF ONE-ELECTRON ORBITALS (STATES)'
     write(11,*) '# AND UNDERLYING SINGLE PARTICLE FUNCTIONS (STURMIANS).'
     write(11,*) '#'
     if (rearranged) then
        write(11,*) '# TWO-ELECTRON STATES HAVE GONE THROUGH REARRANGE12.'
     else
        write(11,*) '# TWO-ELECTRON STATES HAVE NOT BEEN REARRANGED.'
     end if
     write(11,*) '#'

     write(11,*) '# CORE ORBITALS:'
     write(11,*) '# no.  par  m  :  ind  l   1st value   magnitude'
     do jOrb = 1, TargetStates2el%nicm

        jCore = TargetStates2el%ncore(jOrb)
!        jCore = jOrb

        orbital => TargetStates%b(jCore)
        write(11,'(I6,SP,2I4)',advance='no') jCore, get_par(orbital), get_ang_mom_proj(orbital)

        do j = 1, get_nam(orbital)
           ind = get_na(orbital,j)
           spFunction => bst_nr%b(ind)
           spVec => fpointer(spFunction)
           write(11,'(I8,I3,2Es12.4/14X)',advance='no') ind, get_ang_mom(spFunction), spVec(1), sum(spVec(:)*spVec(:))
        end do
        write(11,*)

     end do

     write(11,'(/A/)') ' # TWO-ELECTRON TARGET STATES:'
     do jState = 1, TargetStates2el%Nmax
        state2e => TargetStates2el%b(jState)
        write(11,'(A,I4,SP,3(A,I2),SS,A,F12.8,A)') ' # State no.', jState, ' with (Par,M,S) = (', get_par(state2e),',',get_ang_mom_proj(state2e),',',nint(get_spin(state2e)),') and energy =', get_energy(state2e), ' a.u.'

        write(11,*) '# CORE  FCO  par  m  :  ind  l   1st value   magnitude'
        do jOrb = 1, get_nam(state2e)
           ind1 = get_na(state2e,jOrb,1); ind2 = get_na(state2e,jOrb,2)
           if (ind1 > ind2) then
              write(11,'(2I6,A/)') ind1,ind2, ' antisymmetrisation'
              cycle
           end if

           orbital => TargetStates%b(ind2)
           write(11,'(2I6,SP,2I4)',advance='no') ind1,ind2, get_par(orbital), get_ang_mom_proj(orbital)
           
           do j = 1, get_nam(orbital)
              ind = get_na(orbital,j)
              spFunction => bst_nr%b(ind)
              spVec => fpointer(spFunction)
              write(11,'(I8,I3,2Es12.4/20X)',advance='no') ind, get_ang_mom(spFunction), spVec(1), sum(spVec(:)*spVec(:))
           end do
           write(11,*)
        end do
        write(11,*)

     end do

     write(*,*) 'Made ', filename
     close(11)

   end subroutine write_configs

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

   subroutine orthonormal12(TargetStates, TargetStates2el, bst)
     !
     use sturmian_class
     implicit none

     type(basis_state), intent(in) :: TargetStates, TargetStates2el
     type(basis_sturmian_nr), intent(in) :: bst

     type(basis_state) :: tmpStates2el
     type(state), pointer :: statef,statef1,statef2, statei,statei1,statei2
     integer :: n2e,j2ef,j2ei, jf,ji,stf1,sti1, ncff,ncfi,jcff,jcfi, spf1,spf2,spi1,spi2, timeStart,timeStop,timeRate
     real(dpf) :: overlap,overlap2, CIf,CIi, acc


     open(11, file='orthonormality')
     write(11,*) '# CHECKING ORTHONORMALITY OF TWO-ELECTRON TARGET STATES'
     write(11,*) '# IN FULL, NUSE_ST, AND NUSE_SP BASIS REPRESENTATIONS.'

     n2e = TargetStates2el%Nmax
     call new_basis(tmpStates2el, n2e, .true., TargetStates2el%basis_type)
     do j2ef = 1, n2e
        call copy( tmpStates2el%b(j2ef), TargetStates2el%b(j2ef) )
     end do
     call convert_from_st_to_sp(bst,TargetStates,tmpStates2el, basis_size(bst),TargetStates%Nmax,n2e)
     
     do j2ef = 1, n2e
        write(11,'(/A,I4)') '# final state index:', j2ef
        write(11,'(4A28)') '# initial state index', 'full    time', 'nuse_st    time', 'nuse_sp    time'

        acc = 0d0
        do j2ei = 1, n2e
           statef => TargetStates2el%b(j2ef); ncff = statef%nam
           statei => TargetStates2el%b(j2ei); ncfi = statei%nam
           write(11,'(I28)',advance='no') j2ei

           !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
           ! FULL REPRESENTATION.
           ! List of configurations, i.e. two one-electron orbital states.
           ! Each orbital state links to one or more underlying sp functions.

           call system_clock(timeStart,timeRate)
           overlap = 0d0
           do jcff = 1, ncff   ! Configurations of the final 2e state.
              statef1 => TargetStates%b( statef%na(jcff) )   ! One-electron
              statef2 => TargetStates%b( statef%nb(jcff) )   ! orbital states.
              CIf = statef%CI(jcff)   ! CI coefficient for this configuration.
              do jcfi = 1, ncfi   ! Configurations of the initial 2e state.
                 statei1 => TargetStates%b( statei%na(jcfi) )
                 statei2 => TargetStates%b( statei%nb(jcfi) )
                 CIi = statei%CI(jcfi)

                 overlap = overlap + CIf*CIi * ovlp_st(statef1,statei1,bst) * ovlp_st(statef2,statei2,bst)
              end do
           end do
           call system_clock(timeStop)

           write(11,'(F20.12,F8.4)',advance='no') overlap, dble(timeStop-timeStart)/dble(timeRate)
           if (j2ef == j2ei) overlap = overlap - 1d0
           if (abs(overlap) > abs(acc)) acc = overlap

           !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
           ! NUSE_ST REPRESENTATION.
           ! List of one-electron orbital states used by the two-electron state.
           ! Avoids re-evaluating r1 orbitals from multiple configurations.

           call system_clock(timeStart)
           overlap = 0d0
           do jf = 1, statef%nusemax
              stf1 = statef%nuse(jf); statef1 => TargetStates%b(stf1)
              do ji = 1, statei%nusemax
                 sti1 = statei%nuse(ji); statei1 => TargetStates%b(sti1)

                 overlap2 = 0d0   ! Overlap of the r2 orbital states.
                 do jcff = 1, ncff
                    if (statef%na(jcff) /= stf1) cycle
                    statef2 => TargetStates%b( statef%nb(jcff) )
                    CIf = statef%CI(jcff)
                    do jcfi = 1, ncfi
                       if (statei%na(jcfi) /= sti1) cycle
                       statei2 => TargetStates%b( statei%nb(jcfi) )
                       CIi = statei%CI(jcfi)

                       overlap2 = overlap2 + CIf*CIi * ovlp_st(statef2,statei2,bst)
                    end do
                 end do
                 if (overlap2 == 0d0) cycle

                 overlap = overlap + overlap2 * ovlp_st(statef1,statei1,bst)
              end do
           end do
           call system_clock(timeStop)

           write(11,'(F20.12,F8.4)',advance='no') overlap, dble(timeStop-timeStart)/dble(timeRate)
           if (j2ef == j2ei) overlap = overlap - 1d0
           if (abs(overlap) > abs(acc)) acc = overlap

           !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
           ! NUSE_SP REPRESENTATION.
           ! List of single particle functions used by the two-electron state.
           ! Single-orbital-angular-momentum functions for simple integrals.

           statef => tmpStates2el%b(j2ef)
           statei => tmpStates2el%b(j2ei)
           
           call system_clock(timeStart)
           overlap = 0d0
           do jf = 1, statef%nusemax
              spf1 = statef%nuse(jf)
              do ji = 1, statei%nusemax
                 spi1 = statei%nuse(ji)

                 overlap2 = 0d0
                 do jcff = 1, statef%nam
                    if (statef%na(jcff) /= spf1) cycle
                    spf2 = statef%nb(jcff)
                    CIf = statef%CI(jcff)
                    do jcfi = 1, statei%nam
                       if (statei%na(jcfi) /= spi1) cycle
                       spi2 = statei%nb(jcfi)
                       CIi = statei%CI(jcfi)

                       overlap2 = overlap2 + CIf*CIi * bst%ortint(spf2,spi2)
                    end do
                 end do
                 if (overlap2 == 0d0) cycle
                 
                 overlap = overlap + overlap2 * bst%ortint(spf1,spi1)
              end do
           end do
           call system_clock(timeStop)

           write(11,'(F20.12,F8.4)') overlap, dble(timeStop-timeStart)/dble(timeRate)
           if (j2ef == j2ei) overlap = overlap - 1d0
           if (abs(overlap) > abs(acc)) acc = overlap

        end do ! j2ei
        write(11,'(A,Es20.12)') '# Largest error for this state is', acc

     end do ! j2ef

     print*, 'Made orthonormality'
     close(11)


   end subroutine orthonormal12

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

 end module state_class


