module input_data

  use numbers
  public:: readin, read_debug

  type input
     character(len=10):: target        ! target label
     real(dpf):: energy              ! incident energy in au
     integer:: calculation_type   ! Coordinate system/radial basis type. 0: spherical/nonorthogonal; 1: spherical/orthogonal; 2: spheroidal/orthogonal
     integer :: calculation_mode  ! 0:structure only; 1:scattering; 2:scattering+DCS
     real(dpf):: Rd,origin           ! distance between two nuclei
     real(dpf):: Z1,Z2            ! (one) nuclear charge
     real(dpf):: Zasym               ! molecule charge asymmtery (2*Z(nuclear) - N(electrons))
     real(dpf):: Zplus, Zminus       !spheroidal Zplus and Zminus parameters
     real(dpf):: Zproj               ! The charge of projectile. 1.0 = positron, -1.0 = electron
     real(dpf):: C                   ! the speed of light
     integer:: la_core            ! maximum value of l for atom core orbitals 
     integer, dimension(:), pointer:: npr_core  ! principal quantum number  for last core orbital for each l

!!$ Laguerre basis configurations: 
!!$ MSC_nconfig = 0 Single diagonalisation only using below Laguerre function configurations.
!!$ MSC_nconfig > 0 Used for second diagonalisation.
     integer:: labot, latop       ! minimum value of atom l, maximum value of atom l 
     integer, dimension(:), allocatable :: nps     ! number of functions for each l
     real(dpf),  dimension(:), allocatable ::  alpha  ! exp.fall-off for each l
!!$ Molecular State Configuration: 
!!$ Diagonalisation of the one-electron ion produced MS to replace one-electron configurations.
!!$ For second diagonalisation Molelcar states replaced corresponding Laguerre basis functions below.
     integer:: MSC_nconfig                  ! Number of one-electron configurations replaced with Molecular ion States
     integer :: use_MSC, use_sub
!!$ Different alpha
     real(dpf), allocatable, dimension(:,:):: alpha_nl  ! Storage of alpha n,l functions
     integer :: labot_diff, latop_diff
     integer, allocatable, dimension(:):: nps_diff
     real(dpf), allocatable, dimension(:):: alpha_diff
!!$
     integer::  Mtot_start, Mtot_stop  ! starting and finishing value for the total ang. mom
     integer:: ipar               ! 1 or -1   only this parity to be calculated, 0 if both are calc.
     integer:: nent             ! number if incident states
     integer:: iborn              ! 0 if no analytical Born subtraction is used, 1 if use it, 2 if want to output Born amplitudes
     real(dpf):: rmax                ! max value of radial grid
     real(dpf):: qmax                ! max value of momentum that can be integrated accurately
     integer:: ndouble            ! number of doubling in rgrid
     integer:: npdbl              ! the number of points with the same dx per interval
     integer:: npwave             ! the number of points per oscillation
     integer:: ltmax              ! maximum l in v(1,2) expansion: rpow module
     real(dpf):: formcut
     real(dpf):: regcut
     real(dpf):: expcut
     integer :: iweight
     real(dpf):: corep               ! one-electron pol.potential parameters
     real(dpf):: r0                  ! one-electron pol.potential parameters
     real(dpf):: gamma               ! di-electron pol.potential parameters
     real(dpf):: rho                 ! di-electron pol.potential parameters
     integer :: irec              ! =1 to reconstruct the wf and calculate Zeff
     integer:: iosc               ! =0 no osc.str., =1 absorption osc. str., =2 emission osc.str
     integer:: iSlater  ! =1 if Slater core exchange, =0 if full f.c. exchange  
     integer:: Mt_min,Mt_max      ! max and min values of Mt
     integer, dimension(:,:), pointer:: nst    ! nst(m) number of one-electron-states for each value of m
     integer:: idistpot           ! =0 for plane-waves, =N for distorted-waves, N is a state number
     integer:: Ldw           ! max L value for DW
     integer:: ndw           ! number of functions to be used in diogonalization for DW
     real(dpf):: aldw           ! exp. fall-off for DW
     integer:: ifirst              ! =0 for no exchange, =1 for exchange, =-1 for first order calc., =-10 this sets iAnalitBorn =1 and ifirst=-1
     logical :: exchange
     integer:: iAnalitBorb   ! =0, =1 to stop teh code after Analitical Born cross sections
     real(dpf):: theta                ! theta - nonuniquence
     integer:: inc                ! inc:  =1 for calling rearrange and =0 for not calling it
     integer:: Lpmax                ! max value of L for projectile electron
     integer:: Lmaxkg                ! kgrid: max projectile orbital ang. mom.
     integer, dimension(:), allocatable :: NBND               ! kgrid: number of bound states in dist. potential
     integer, dimension(:,:), allocatable :: nk ! nk(i,l) : kgrid: interval description, number of points
     real(dpf), dimension(:,:), allocatable ::  sk  !  sk(i,l) : kgrid: interval description, boundary value
     integer:: orient              ! 1 for fixed orientation, 0 for orientation averaging of molecule
     real(dpf):: thrad, phrad   ! incident electron angles (degrees) in the body frame

     logical :: combine_potl_files

     logical:: hlike, good_parity           ! true for quasi one-electron atoms, false for quasitwo-electron atoms
     integer:: iBorn_amp       ! 0 if no analytical Born amp files to be written to disk, 1 if amp files to be written to disk
     integer:: lbnd_max   ! max value of L for which bound states of asymptotic potential (Zasym/r + V_{dw}(r) should be calculated
     logical:: non_uniq   ! We solve for non-uniqueness. ifirst = 1, theta /= 0 and Zproj == -1d0
     
     real(dpf):: CI_min   ! minimum value for CI coef.: in 2-elelctron structure code: if |CI| < CI_min it will be set to zero 

     integer :: l_ion_core, l12max, M12max, Mmax2el
     integer, dimension(:), allocatable :: n_ion_core, nkin, nkout
     integer, dimension(:,:,:), allocatable :: nst_2e, nst_2e_bound
  
     integer:: N_core_el
     
! Rav: Below for input data for Ps-channels:
     integer:: l_Ps ! added by Rav, for Ps-channels
     integer:: lpbot, lptop                     ! min and max l of Ps
     integer, dimension(:), allocatable:: npbot, nptop  ! n for each l of Ps
     integer, dimension(:), allocatable:: npsp          ! basis size N for each l of Ps
     real(dpf), dimension(:), allocatable::  alphap        ! fall of parameter of Ps states
     integer::   igz, igp                           ! number of points for z-integration, ipg - for composite mesh xgp    
     logical:: analyticd, numericalv                !  to switch between analytic and numerical calculations
     integer:: lstoppos ! no calculations of Ps-formation for for ! J>lsstoppos
!!! end of Ps input data

     real(dpf) :: eV
     real(dpf) :: eV_old
     
     logical :: print_1el_basis, print_2el_config, print_CI, print_dipole, print_channel_timing, print_pol_contribution

     integer :: num_nat_orb, natorb_GS_M, natorb_GS_par, natorb_GS_spin
     logical :: only_nat_orbs
     logical :: new_input_format
     logical :: pwborn, UBA
     logical :: pseudo_pot
     real(dpf), dimension(0:2) :: pseudo_pot_B, pseudo_pot_Beta 
     real(dpf) :: core_energy

     logical :: skip_vmat
     logical :: optical
     logical :: no_second_spin_channel
     logical :: weak_coupling, weak_coupling_exchange
     logical :: reduce_bound
!     integer :: max_bound
     logical :: load_balance, read_channel_timing
     integer :: NP
     logical :: spectro_core
     real(dpf) :: SF_tol
     logical :: print_vmatrix
     logical :: print_BornME
     logical :: print_halfk
     logical :: print_ionisation

  end type input


  type basis_input
     integer :: labot,latop, mabot,matop
     integer, dimension(:), allocatable :: nps
     real(dpf), dimension(:), allocatable :: alpha
  end type basis_input


  type debug
     integer :: list_sturm_1e = 0
     integer :: n_sturm_1e = 0
     integer, dimension(:), allocatable :: v_sturm_1e

     integer :: list_state_1e = 0
     integer :: n_state_1e = 0
     integer, dimension(:), allocatable :: v_state_1e

     integer :: list_orb_2e = 0
     integer :: n_orb_2e = 0
     integer, dimension(:), allocatable :: v_orb_2e

     integer :: list_state_2e = 0
     integer :: n_state_2e = 0
     integer, dimension(:), allocatable :: v_state_2e

  end type debug


!****

  type(input):: data_in     ! contains all input data
  type(basis_input), target :: basis_1e, basis_2e, sub_basis_2e
  type(basis_input), target :: core_basis_1e, core_basis_2e, core_sub_basis_2e
  type(input):: dataMSC     ! Contains Laguerre basis used in second diagonalisation with Molecular State Configurations
  type(debug) :: data_debug

!****
    
  interface invalid_argument
    module procedure invalid_argument_i, invalid_argument_r, invalid_argument_c
  end interface

contains
  
  subroutine read_debug(self)
    !
    use MPI_module
    implicit none
    type(debug), intent(inout) :: self

    integer :: j
    logical :: ex

    inquire( file='debug.in', exist=ex )
    if ( .not. ex ) then
       write(*,'(/A/)') 'WARNING: debug.in does not exist'
       return
    end if
       
    open(11,file='debug.in')

    read(11,*) self%list_sturm_1e
    read(11,*) self%n_sturm_1e
    allocate( self%v_sturm_1e(self%n_sturm_1e) )
    read(11,*) ( self%v_sturm_1e(j), j=1,self%n_sturm_1e )

    read(11,*) self%list_state_1e
    read(11,*) self%n_state_1e
    allocate( self%v_state_1e(self%n_state_1e) )
    read(11,*) ( self%v_state_1e(j), j=1,self%n_state_1e )

    read(11,*) self%list_orb_2e
    read(11,*) self%n_orb_2e
    allocate( self%v_orb_2e(self%n_orb_2e) )
    read(11,*) ( self%v_orb_2e(j), j=1,self%n_orb_2e )

    read(11,*) self%n_state_2e
    allocate( self%v_state_2e(self%n_state_2e) )
    read(11,*) ( self%v_state_2e(j), j=1,self%n_state_2e )

    close(11)

  end subroutine read_debug

  subroutine readin( self, nfile,iwrite )
    use target_data
    use MPI_module
    use input_parser
    implicit none
    type(input), intent(out) :: self

    integer, intent(in):: nfile
    integer:: nfile_Ps
    integer, intent(in):: iwrite  ! if =1 write to screen the read data, =0 no writing to screen

    logical:: ex
    integer:: iostat_data, icheck
    integer:: l, latop,labot
    real(dpf):: en_eV, eV
    integer:: Ltt, LSW, LSW_prev, i, NBNDt, Lbot
    integer :: nkt
    real(dpf) :: skt
    integer:: n, n_elec 
    character(len=50) :: tempchar
    character(len=9) :: projectile
    character(len=:), allocatable :: form
    character(len=8) :: header

    logical :: origin_not_specified = .false.
    logical :: l_ps_not_specified = .false.
    
    character(len=255) :: block, blockname, label, w1, w2, w3
    logical :: eof, templ, templ1, templ2
    logical :: scattering, use_sub, born, exchange, inc, distpot
    real(dpf) :: tempr
    integer :: tempi, tempi1, M, par_start, par_stop
    integer, dimension(:), allocatable :: L_k_input

    logical :: scattering_system_read, output_read, first_basis_read, core_first_basis_read, one_el_states_read, &
              &core_one_el_states_read, second_basis_read, core_second_basis_read, two_el_configs_read, two_el_states_read, &
              &radial_grid_read, kgrid_read, calculation_options_read, bound_two_el_states_read

    logical :: read_B0, read_B1, read_B2, read_Beta0, read_Beta1, read_Beta2

    logical :: star, star_prev, use_second_basis, core_input

    type(basis_input), pointer :: basis_1e_pntr, basis_2e_pntr, sub_basis_2e_pntr

    self%eV=27.21138505d0
    self%eV_old=27.2116d0

    eV = self%eV  !  27.2116
    
    !default values for parameters which don't need to be specified:
      self%origin = 0.5d0
      self%calculation_mode = 0 !default to structure only
      self%combine_potl_files = .false.
      self%use_sub = 0
      self%iweight = 0 ! 0 = Simpson
      self%irec = 0
      self%iosc = 0
      self%idistpot = 0
      self%Ldw = 0
      self%ndw = 0
      self%aldw = 0.0d0
      self%ifirst = 1
      self%iAnalitBorb = 0
      self%theta = 0.0d0
      self%inc = 1
      self%orient = 0
      self%thrad = 0.0d0
      self%phrad = 0.0d0
      self%iBorn_amp = 0
      self%print_1el_basis = .false.
      self%print_2el_config = .false.
      self%print_CI = .false.
      self%print_dipole = .false.
      self%print_pol_contribution = .false.
      self%num_nat_orb = 0
      self%only_nat_orbs = .false.
      self%good_parity = .true.
      self%l_ps = -1
      self%optical = .false.
      self%no_second_spin_channel = .false.
      self%weak_coupling = .false.
      self%weak_coupling_exchange = .false.
      self%reduce_bound = .false.
!      self%max_bound = -1
      self%load_balance = .false.
      self%read_channel_timing = .false.
      self%NP = -1 !-1 = ALL
      self%CI_min = 1.0D-6
      self%spectro_core = .false.
      self%SF_tol = 1.0d0
      self%print_vmatrix = .false.
      self%print_BornME = .false.
      self%print_halfk = .false.
      self%skip_vmat = .false.
      self%pwborn = .false.
      self%N_core_el = 0
      self%corep=0.0d0
      self%r0   =0.0d0
      self%non_uniq = .false.
      self%pseudo_pot = .false.

    inquire(file='data.in',exist=ex)
    if(.not. ex) then
      error stop '***ERROR: Cannot find data.in file'
    endif

    if(myid==0) write(*,*) 'Found data.in file'

    open(unit=nfile,file='data.in',action='read')
    read(nfile,'(A8)') header
    if(header == '#MCCC-FN') then
      if(myid==0) write(*,*) ' -- New data.in format'
      self%new_input_format = .true.
    elseif(header == '#CONVERT') then

      if(myid==0) then
        write(*,*) ' -- Converting old data.in to new format'
        call readin_old(self, basis_1e,basis_2e,sub_basis_2e, nfile,iwrite)
        call convert_input_file(self, basis_1e, basis_2e, sub_basis_2e)
      endif

      !if(mpi) call MPI_Finalize(ierr)
      stop
    else
      if(myid==0) write(*,*) ' -- Old data.in format - calling readin_old'
      rewind(nfile)
      call readin_old(self, basis_1e,basis_2e,sub_basis_2e, nfile,iwrite)
      self%new_input_format = .false.
      return
    endif

    !Temporary values for variables not presently used
    self%la_core = 0
    allocate(self%npr_core(0:0))
    self%npr_core = 0
    self%gamma=0.0d0
    self%rho  =0.0d0
    
    call input_options(skip_blank_lines=.true., echo_lines=.false.)
    
    core_input = .false.

    !set default values for input parameters
    self%calculation_type = -1 ! -1 = not set
    scattering = .false.

    !First scan through file to determine coordinate system because that effects how some input is handled
    do while (self%calculation_type == -1)
      call read_line(eof,nfile)
      if(eof) then
        if(myid==0) write(*,*) '*** ERROR reading input: Could not find CALCULATION_MODES block'
        error stop
      endif
      call read_block(block,blockname,core_input) !read the block heading
      select case (block)
        case('CALCULATION_MODES') !found the CALCULATION_MODES block
          if(nitems>1) call heading_error(block) 
          do !read lines inside CALCULATION_MODES block
            call read_line(eof,nfile)
            if(eof) call eof_during_block(block)
            call read_label(label,block)
            select case (label)
              case('COORDINATE_SYSTEM')
                if(nitems /= 3) call nitem_error(label,block)
                call readu(w1)
                select case(w1)
                  case('SPHERICAL')
                    self%calculation_type = 0
                  case('SPHEROIDAL')
                    self%calculation_type = 2
                  case default
                    if(myid==0) call report('invalid coordinate system')
                    error stop
                end select
              case('SCATTERING')
                call switch(scattering,label,block)
                if(scattering) then
                  self%calculation_mode = 1 !default to no DCS or TV files
                else
                  self%calculation_mode = 0
                endif
              case('COMBINE_POTL_FILES')
                call switch(self%combine_potl_files,label,block)
              case('ANALYTICAL_BORN_ONLY')
                call switch(born,label,block)
                if(born) then
                  self%iAnalitBorb = 1
                else
                  self%iAnalitBorb = 0
                endif
              case('END') !Found end of CALCULATION_MODES block
                if(self%calculation_type == -1) then
                  if(myid==0) write(*,*) '*** ERROR reading input: coordinate system not specified in <CALCULATION_MODES> block'
                  error stop
                endif
                if(self%combine_potl_files .and. .not. scattering) then
                  if(myid==0) write(*,*) '*** ERROR reading input: turn SCATTERING ON if you want to combine potl files'
                  error stop
                endif
                exit
              case default
                call invalid_label_in_block(label,block)
            end select 
          enddo !end reading lines in CALCULATION_MODES block

        case default 
          !Skip any blocks which are not CALCULATION_MODES
      end select
    enddo

    rewind(nfile) !rewind data.in file in case there were blocks above the CALCULATION_MODES block
    call rewind_file
    blocks_read = '' !reset blocks_read array
       

    !logical variables to keep track of whether or not required input blocks are present in data.in file
    !Don't bother with CALCULATION_MODES because that was explicitly checked above
    scattering_system_read= .false.
    output_read= .false.
    first_basis_read= .false.
    core_first_basis_read= .false.
    one_el_states_read= .false.
    core_one_el_states_read= .false.
    second_basis_read= .false.
    core_second_basis_read= .false.
    two_el_configs_read= .false.
    two_el_states_read= .false.
    bound_two_el_states_read= .false.
    radial_grid_read= .false.
    kgrid_read= .false.
    calculation_options_read = .false.
    
      
      !default "bad" values for parameters to check if they are provided 
      self%target = 'UNDEFINED'
      projectile = 'UNDEFINED'
      self%energy = -1000d0
      self%Rd = -1000d0
      self%labot = -1000
      self%latop = -1000
      self%labot_diff = -1000
      self%latop_diff = -1000
      self%Mtot_start = -1000
      self%Mtot_stop = -1000
      self%ipar = -1000
      self%nent = -1000
      self%iborn = -1000
      self%rmax = -1000d0
      self%qmax = -1000d0
      self%ndouble = -1000
      self%npdbl = -1000
      self%npwave = -1000
      self%ltmax = -1000
      self%formcut = -1000d0
      self%regcut = -1000d0
      self%expcut = -1000d0
      self%Mt_min = -1000
      self%Mt_max = -1000
      self%Lpmax = -1000
      self%Lmaxkg = -1000
      self%l_ion_core = -1000
      self%l12max = -1000
      self%M12max = -1000
      self%Mmax2el = -1000
      self%natorb_GS_M = -1000
      self%natorb_GS_par = -1000
      self%natorb_GS_spin = -1000

      !Doesn't do anything but the Cray compiler complains about them potentially being referenced later on without being defined
      par_start = -1000
      par_stop = -1000
      en_eV = -1000.0d0

    !Main loop over data.in file
    do 
      call read_line(eof,nfile)
      if(eof) exit
      call read_block(block,blockname,core_input) !read the block heading

      select case (block)
        case('CALCULATION_MODES') !found the CALCULATION_MODES block
          !skip since we already read this in earlier
          do !read lines inside CALCULATION_MODES block
            call read_line(eof,nfile)
            call read_label(label,block)
            select case (label)
              case('END') !Found end of CALCULATION_MODES block
                exit
            end select 
          enddo !end reading lines in CALCULATION_MODES block
       
        case('CORE')
          if(core_input) then
            if(myid==0) call report('found CORE while already in a CORE block')
            error stop
          else
            core_input = .true.
          endif

        case('END_CORE')
          if(.not.core_input) then
            if(myid==0) call report('found END_CORE without corresponding CORE')
            error stop
          else
            core_input = .false.
          endif

        case('SCATTERING_SYSTEM') !found the SCATTERING_SYSTEM block
          scattering_system_read = .true.
          if(nitems>1) call heading_error(block) 
          do !read lines inside SCATTERING_SYSTEM block
            call read_line(eof,nfile)
            if(eof) call eof_during_block(block)
            call read_label(label,block)
            select case (label)
              case('TARGET')
                call reada(self%target)
                select case(self%target)
                  case('H2')
                    self%Z1 = 1
                    self%Z2 = 1
                    self%Zasym = 0
                    self%N_core_el = 0
                  case('H2+')
                    self%Z1 = 1
                    self%Z2 = 1
                    self%Zasym = 1
                    self%N_core_el = 0
                  case('HeH')
                    self%Z1 = 2
                    self%Z2 = 1
                    self%Zasym = 0
                    self%N_core_el = 2
                  case('HeH+')
                    self%Z1 = 2
                    self%Z2 = 1
                    self%Zasym = 1
                    self%N_core_el = 0
                  case('HeH++')
                    self%Z1 = 2
                    self%Z2 = 1
                    self%Zasym = 2
                    self%N_core_el = 0
                  case('LiH+')
                    self%Z1 = 3
                    self%Z2 = 1
                    self%Zasym = 1
                    self%N_core_el = 2
                  case('LiH++')
                    self%Z1 = 3
                    self%Z2 = 1
                    self%Zasym = 2
                    self%N_core_el = 0
                  case('LiH+++')
                    self%Z1 = 3
                    self%Z2 = 1
                    self%Zasym = 3
                    self%N_core_el = 0
                  case('He2+++')
                    self%Z1 = 2
                    self%Z2 = 2
                    self%Zasym = 3
                    self%N_core_el = 0
                  case('He2++')
                    self%Z1 = 2
                    self%Z2 = 2
                    self%Zasym = 2
                    self%N_core_el = 0
                  case default
                    call invalid_argument(label, block, self%target)
                    error stop
                end select
                if (self%Z1+self%Z2 - self%Zasym - self%N_core_el == 1) then   ! One-electron target.
                  self%hlike = .true.
                elseif (self%Z1+self%Z2 - self%Zasym - self%N_core_el == 2) then   ! Two-electron target.
                  self%hlike = .false.
                endif

              case('INTERNUCLEAR_SEPARATION')
                call readf(self%Rd)
                if(self%Rd < 0) call invalid_argument(label,block,self%Rd) 
              case('ORIGIN')
                call readf(self%origin)
                if(self%origin < 0 .or. self%origin > 1) call invalid_argument(label,block,self%origin) 
                if(self%origin /= 0.5d0) self%good_parity = .false.
              case('PROJECTILE')
                call readu(projectile)
                select case(projectile)
                  case('ELECTRON')
                    self%Zproj = -1
                  case('POSITRON')
                    self%Zproj = 1
                  case default
                    call invalid_argument(label,block,projectile)
                end select
              case('ENERGY')
                if(nitems == 3) then !energy given without unit - assume eV
                  call readf(tempr)
                  self%energy = tempr / eV
                  en_eV = tempr
                elseif(nitems == 4) then !energy given with unit
                  call readf(tempr)
                  if(tempr <= 0.0d0) call invalid_argument(label,block,tempr)
                  call readu(w1)
                  select case (w1)
                    case('EV')
                      self%energy = tempr / eV
                  en_eV = tempr
                    case('HARTREE','HA')
                      self%energy = tempr
                      en_eV = tempr * eV
                    case default
                      if(tempr <= 0.0d0) call invalid_argument(label,block,w1)
                      error stop
                  end select
                else
                  call nitem_error(label,block)
                endif

              case('M_TOT')
                if(nitems == 3) then
                  call readi(self%Mtot_start)
                  self%Mtot_stop = self%Mtot_start
                  if(self%Mtot_stop < 0) call invalid_argument(label,block,self%Mtot_stop)
                elseif(nitems == 4) then
                  call readi(self%Mtot_start)
                  call readi(self%Mtot_stop)
                  if(self%Mtot_start < 0 .or. self%Mtot_stop < 0 .or. self%Mtot_stop < self%Mtot_start) then
                    write(w1,'(I0,", ",I0)') self%Mtot_start, self%Mtot_stop
                    call invalid_argument(label, block, w1)
                  endif
                else 
                  call nitem_error(label,block)
                endif

              case('PARITY')
                if(nitems == 3) then
                  call readi(tempi)
                  tempi1 = tempi
                elseif(nitems == 4) then
                  call readi(tempi)
                  call readi(tempi1)
                else 
                  call nitem_error(label,block)
                endif
                if(tempi==-1.and.tempi1==1 .or. tempi==1.and.tempi1==-1) then !both parities
                  self%ipar = 0
                  par_start = -1; par_stop = 1
                elseif(tempi==1.and.tempi1==1) then
                  self%ipar = 1
                  par_start = 1; par_stop = 1
                elseif(tempi==-1.and.tempi1==-1) then
                  self%ipar = -1
                  par_start = -1; par_stop = -1
                elseif(tempi==0.and.tempi1==0) then
                  self%ipar = 0
                  self%good_parity = .false.
                  par_start = 0; par_stop = 0
                else
                  if(myid==0) call report('invalid argument(s) to the <PARITY> keyword')
                  error stop
                endif

              case('NUMBER_ENTRANCE_STATES')
                call readi(self%nent)
                if(self%nent <= 0) call invalid_argument(label, block, self%nent)

              case('PROJECTILE_LMAX')
                call readi(self%Lpmax)
                if(self%Lpmax < 0) call invalid_argument(label, block, self%Lpmax)

              case('END')
                !Check for missing input:
                if(trim(adjustl(self%target)) == 'UNDEFINED') then
                  if(myid==0) write(*,*) '*** ERROR reading input: missing TARGET in the <SCATTERING_SYSTEM> block'
                  error stop
                endif
                if(trim(adjustl(projectile)) == 'UNDEFINED') then
                  if(myid==0) write(*,*) '*** ERROR reading input: missing PROJECTILE in the <SCATTERING_SYSTEM> block'
                  error stop
                endif
                if(self%Rd == -1000) then
                  if(myid==0) write(*,*) '*** ERROR reading input: missing INTERNUCLEAR_SEPARATION in the <SCATTERING_SYSTEM> block'
                  error stop
                endif
                if(self%origin == -1000 .and. self%calculation_type == 0) then 
                  if(myid==0) write(*,*) '*** WARNING reading input: missing ORIGIN in the <SCATTERING_SYSTEM> block'
                  if(myid==0) write(*,*) '                                --> setting origin = 0.5'
                  self%origin = 0.5d0
                endif
                if(scattering.and.self%energy == -1000d0) then
                  if(myid==0) write(*,*) '*** ERROR reading input: missing ENERGY in the <SCATTERING_SYSTEM> block'
                  error stop
                endif
                if(scattering.and.self%Mtot_stop == -1000) then
                  if(myid==0) write(*,*) '*** ERROR reading input: missing M_TOT in the <SCATTERING_SYSTEM> block'
                  error stop
                endif
                if(scattering.and.self%ipar == -1000) then
                  if(myid==0) error stop '*** ERROR reading input: missing PARITY in the <SCATTERING_SYSTEM> block'
                  error stop
                endif
                if(self%nent == -1000) then 
                  write(*,*) '*** WARNING reading input: missing NUMBER_ENTRANCE_STATES in the <SCATTERING_SYSTEM> block'
                  write(*,*) '                                --> setting NUMBER_ENTRANCE_STATES = 1'
                  self%nent = 1
                endif
                if(scattering.and.self%Lpmax == -1000) then
                  if(myid==0) write(*,*) '*** ERROR reading input: missing PROJECTILE_LMAX in the <SCATTERING_SYSTEM> block'
                  error stop
                endif
                
                if(self%Z1 /= self%Z2 .and. self%Rd > 0.0d0) self%good_parity = .false.
                self%Zplus = self%Zproj * (self%Z1 + self%Z2)
                self%Zminus = self%Zproj * (self%Z1 - self%Z2)

                exit

              case default
                call invalid_label_in_block(label,block)
            end select 
          enddo !end reading lines in SCATTERING_SYSTEM block
        
        case('OUTPUT') !found the OUTPUT block
          output_read = .true.
          if(nitems>1) call heading_error(block) 
          templ1 = .false. !DCS
          templ2 = .false. !TV_FILES
          do !read lines inside OUTPUT block
            call read_line(eof,nfile)
            if(eof) call eof_during_block(block)
            call read_label(label,block)
            select case (label)
              case('OSCILLATOR_STRENGTHS')
                call switch(templ,label,block)
                if(templ) then
                  self%iosc = 1
                else
                  self%iosc = 0
                endif
              case('DIPOLE_MOMENTS')
                call switch(self%print_dipole,label,block)
              case('POL_CONTRIBUTIONS')
                call switch(self%print_pol_contribution,label,block)
              case('DCS')
                call switch(templ1,label,block)
              case('TV_FILES')
                call switch(templ2,label,block)
              case('IONISATION_RUNNING_SUM')
                call switch(self%print_ionisation,label,block)
              case('ONE_EL_BASIS_LIST')
                call switch(self%print_1el_basis,label,block)
              case('TWO_EL_CONFIG_LIST')
                call switch(self%print_2el_config,label,block)
              case('CI_COEFFICIENTS')
                call switch(self%print_CI,label,block)
              case('CHANNEL_TIMING')
                call switch(self%print_channel_timing,label,block)
              case('VMATRIX')
                call switch(self%print_vmatrix,label,block)
                if(self%print_vmatrix .and. ntasks > 1) then
                  if(myid==0) call report('can only print V matrix with single-task calculations')
                  error stop
                endif
              case('KMATRIX')
                call switch(self%print_halfk,label,block)
              case('BORNME')
                call switch(self%print_BornME,label,block)
              case('END')
                !Validate output options
                if(scattering) then
                  if(.not.templ1 .and. .not.templ2) then !no DCS or TV files
                    self%calculation_mode = 1
                  elseif(templ1 .and. templ2) then !DCS and TV files
                    self%calculation_mode = 2
                  elseif(.not.templ1 .and. templ2) then !TV files and no DCS
                    self%calculation_mode = 3
                  else
                    if(myid==0) write(*,*) '*** ERROR reading input: cannot output DCS but not TV files'
                    error stop
                  endif
                else
                  self%calculation_mode = 0
                endif

                if(self%print_dipole .and. self%iosc == 0) then
                  if(myid==0) error stop '*** ERROR reading input: cannot print dipole moments but not oscillator strengths'
                  error stop
                endif

                exit
              case default
                call invalid_label_in_block(label,block)
            end select 
          enddo !end reading lines in OUTPUT block

        case('FIRST_DIAGONALISATION_BASIS') !found the FIRST_DIAGONALISATION_BASIS block
          if(core_input) then
            core_first_basis_read = .true.
            basis_1e_pntr => core_basis_1e
          else
            first_basis_read = .true.
            basis_1e_pntr => basis_1e
          endif

          if(nitems>1) call heading_error(block) 
          basis_1e_pntr%latop = -1
          do !read lines inside FIRST_DIAGONALISATION_BASIS block
            call read_line(eof,nfile)
            if(eof) call eof_during_block(block)
            call read_label(label,block)
            select case (label)
              case('LMAX')
                basis_1e_pntr%labot = 0
                call readi(basis_1e_pntr%latop)
                if(basis_1e_pntr%latop < 0) call invalid_argument(label,block,basis_1e_pntr%latop)
                allocate(basis_1e_pntr%nps(0:basis_1e_pntr%latop), basis_1e_pntr%alpha(0:basis_1e_pntr%latop))
              case('BASIS_SIZE')
                call read_basis_size(basis_1e_pntr%latop, basis_1e_pntr%nps, block)
              case('EXPONENTIAL_FALLOFF')
                call read_alpha(basis_1e_pntr%latop, basis_1e_pntr%alpha, block)
              case('END')
                !Check for missing input
                if(basis_1e_pntr%latop == -1000) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing LMAX in the <'//trim(adjustl(blockname))//'> block'
                  error stop
                endif
                if(.not.allocated(basis_1e_pntr%nps)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing BASIS_SIZE in the <'//trim(adjustl(blockname))//'> block'
                  error stop
                endif
                if(.not.allocated(basis_1e_pntr%alpha)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing EXPONENTIAL_FALLOFF in the <'//trim(adjustl(blockname))//'> block'
                  error stop
                endif
                exit
              case default
                call invalid_label_in_block(label,block)
            end select 
          enddo !end reading lines in OUTPUT block
        
        case('ONE_ELECTRON_STATES') !found the ONE_ELECTRON_STATES block
          one_el_states_read = .true.
          if(nitems>1) call heading_error(block) 
          do !read lines inside ONE_ELECTRON_STATES block
            call read_line(eof,nfile)
            if(eof) call eof_during_block(block)
            call read_label(label,block)
            select case (label)
              case('M_MAX')
                self%Mt_min = 0
                call readi(self%Mt_max)
                if(self%Mt_max < 0) call invalid_argument(label,block,self%Mt_max)
                if ( self%calculation_type==2 .or. self%calculation_type==3 ) then
                   basis_1e%mabot = self%Mt_min
                   basis_1e%matop = self%Mt_max
                end if   
                allocate(self%nst(0:self%Mt_max,-1:1))
                self%nst = -1000

              case('NST_POS_PAR')
                call read_ALL_or_N(self%nst(0,1),label,block)
                if(nitems == 2 + self%Mt_max+1) then
                  do M=1, self%Mt_max
                    call read_ALL_or_N(self%nst(M,1),label,block)
                  enddo
                else
                  call nitem_error(label,block)
                endif

              case('NST_NEG_PAR')
                call read_ALL_or_N(self%nst(0,-1),label,block)
                if(nitems == 2 + self%Mt_max+1) then
                  do M=1, self%Mt_max
                    call read_ALL_or_N(self%nst(M,-1),label,block)
                  enddo
                else
                  call nitem_error(label,block)
                endif

              case('END')
                !check for missing input:
                if(self%Mt_max == -1000) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing M_MAX in the <ONE_ELECTRON_STATES> block'
                  error stop
                endif
                if(any(self%nst(:,1) == -1000)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NST_POS_PAR in the <ONE_ELECTRON_STATES> block'
                  error stop
                endif
                if(any(self%nst(:,-1) == -1000)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NST_NEG_PAR in the <ONE_ELECTRON_STATES> block'
                  error stop
                endif
                exit
              case default
                call invalid_label_in_block(label,block)
            end select 
          enddo !end reading lines in ONE_ELECTRON_STATES block
        
        case('SECOND_DIAGONALISATION_BASIS') !found the SECOND_DIAGONALISATION_BASIS block
          second_basis_read = .true.
          if(nitems>1) call heading_error(block) 
          use_second_basis = .true.
          do !read lines inside SECOND_DIAGONALISATION_BASIS block
            call read_line(eof,nfile)
            if(eof) call eof_during_block(block)
            call read_label(label,block)
            select case (label)
              case('USE_SECOND_BASIS')
                call switch(use_second_basis,label,block)
              case('LMAX')
                basis_2e%labot = 0
                call readi(basis_2e%latop)
                if(basis_2e%latop < 0) call invalid_argument(label,block,basis_1e%latop)
                basis_2e%mabot = basis_2e%labot
                basis_2e%matop = basis_2e%latop
                allocate(basis_2e%nps(0:basis_2e%latop), basis_2e%alpha(0:basis_2e%latop))
              case('BASIS_SIZE')
                call read_basis_size(basis_2e%latop, basis_2e%nps, block)
              case('EXPONENTIAL_FALLOFF')
                call read_alpha(basis_2e%latop, basis_2e%alpha, block)
              case('INSERT_ONE_EL_STATES')
                call read_ALL_or_N(self%use_MSC,label,block)
              case('INSERT_NATURAL_ORBITALS')
                call read_ALL_or_N(self%num_nat_orb,label,block)
              case('ONLY_NATURAL_ORBITALS')
                call switch(self%only_nat_orbs,label,block)
              case('NATORB_GS_M')
                call readi(self%natorb_GS_M)
                if(self%natorb_GS_M < 0) call invalid_argument(label,block,self%natorb_GS_M)
              case('NATORB_GS_PAR')
                call readi(self%natorb_GS_PAR)
                if(self%natorb_GS_par < 0) call invalid_argument(label,block,self%natorb_GS_par)
              case('NATORB_GS_SPIN')
                call readi(self%natorb_GS_SPIN)
                if(self%natorb_GS_spin < 0) call invalid_argument(label,block,self%natorb_GS_spin)
              case('USE_SUB_BASIS')
                call switch(use_sub,label,block)
                if(use_sub) then
                  self%use_sub = 1
                else
                  self%use_sub = 0
                endif
              case('SUB_LMAX')
                sub_basis_2e%labot = 0
                call readi(sub_basis_2e%latop)
                if(sub_basis_2e%latop < 0) call invalid_argument(label,block,sub_basis_2e%latop)
                if(sub_basis_2e%latop > basis_2e%latop) then
                  if(myid==0) write(*,*) "***ERROR reading input: sub-basis latop not within the confines of the second Laguerre basis"
                  error stop 
                end if
                sub_basis_2e%mabot = sub_basis_2e%labot
                sub_basis_2e%matop = sub_basis_2e%latop
                  
                allocate(sub_basis_2e%nps(0:sub_basis_2e%latop), sub_basis_2e%alpha(0:sub_basis_2e%latop))
              case('SUB_BASIS_SIZE')
                call read_basis_size(sub_basis_2e%latop, sub_basis_2e%nps, block)
              case('SUB_BASIS_EXP_FALLOFF')
                call read_alpha(sub_basis_2e%latop, sub_basis_2e%alpha, block)
              case('END')
                !Check for missing input
                if(basis_2e%latop == -1000) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing LMAX in the <SECOND_DIAGONALISATION_BASIS> block'
                  error stop
                endif
                if(.not.allocated(basis_2e%nps)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing BASIS_SIZE in the <SECOND_DIAGONALISATION_BASIS> block'
                  error stop
                endif
                if(.not.allocated(basis_2e%alpha)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing EXPONENTIAL_FALLOFF in the <SECOND_DIAGONALISATION_BASIS> block'
                  error stop
                endif
                if(self%num_nat_orb == -1 .or. self%num_nat_orb > 0) then
                  if(self%natorb_GS_M == -1000 .or. self%natorb_GS_par == -1000 .or. self%natorb_GS_spin == -1000) then
                    if(myid==0) write(*,*) '*** ERROR reading input: Missing symmetry for ground state used to generate natural orbitals'
                    error stop
                  endif
                endif
                if(.not.use_second_basis) self%use_MSC = -1
                exit
              case default
                call invalid_label_in_block(label,block)
            end select 
          enddo !end reading lines in SECOND_DIAGONALISATION_BASIS block
        
        case('TWO_ELECTRON_CONFIGURATIONS') !found the TWO_ELECTRON_CONFIGURATIONS block
          if(.not.second_basis_read) then          
            if(myid==0) write(*,*) '*** ERROR reading input: place <TWO_ELECTRON_CONFIGURATIONS> after <SECOND_DIAGONALISATION_BASIS>' 
            error stop
          endif
          two_el_configs_read = .true.
          if(nitems>1) call heading_error(block) 
          do !read lines inside TWO_ELECTRON_CONFIGURATIONS block
            if(use_second_basis) then
              self%l12max = basis_2e%latop
            else
              self%l12max = basis_1e%latop
            endif
            call read_line(eof,nfile)
            if(eof) call eof_during_block(block)
            call read_label(label,block)
            select case (label)
              case('L_ION_CORE')
                if(nitems /= 3) call nitem_error(label,block)
                call readi(self%l_ion_core)
                if(self%l_ion_core < 0) call invalid_argument(label,block,self%l_ion_core)
                allocate(self%n_ion_core(0:self%l_ion_core))
              case('N_ION_CORE')
                if(nitems /= 2 + self%l_ion_core+1) call nitem_error(label,block)
                do i=0, self%l_ion_core
                  call readi(self%n_ion_core(i))
                  if(self%n_ion_core(i) < 0) call invalid_argument(label,block,self%n_ion_core(i))
                enddo
              !case('L12MAX')
              !  if(nitems /= 3) call nitem_error(label,block)
              !  call readi(self%l12max)
              !  if(self%l12max < 0) call invalid_argument(label,block,self%l12max)
              case('M12MAX')
                if(nitems /= 3) call nitem_error(label,block)
                call readi(self%M12max)
                if(self%M12max < 0) call invalid_argument(label,block,self%M12max)
              case('NKIN')
                !if(nitems < 2 + self%l12max+1) call nitem_error(label,block)
                allocate(self%nkin(0:self%l12max))
                !do i=0, self%l12max
                do i=0, min(self%l12max, nitems-3)
                  call readi(self%nkin(i))
                  if(self%nkin(i) < 0) call invalid_argument(label,block,self%nkin(i))
                enddo
                if(nitems-3 < self%l12max) self%nkin(nitems-2:self%l12max) = 0
                allocate(self%nkout(0:self%l12max))
                if(use_second_basis) then
                  self%nkout = maxval(basis_2e%nps)
                else
                  self%nkout = maxval(basis_1e%nps)
                endif
              !case('NKOUT')
              !  if(nitems /= 2 + self%l12max+1) call nitem_error(label,block)
              !  allocate(self%nkout(0:self%l12max))
              !  do i=0, self%l12max
              !    call readi(self%nkout(i))
              !    if(self%nkout(i) < 0) call invalid_argument(label,block,self%nkout(i))
              !  enddo

              case('END')
                !Check for missing input
                if(self%l_ion_core == -1000) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing L_ION_CORE in the <TWO_ELECTRON_CONFIGURATIONS> block'
                  error stop
                endif
                !if(self%l12max == -1000) then
                !  if(myid==0) write(*,*) '*** ERROR reading input: Missing L12MAX in the <TWO_ELECTRON_CONFIGURATIONS> block'
                !  error stop
                !endif
                if(self%M12max == -1000) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing M12MAX in the <TWO_ELECTRON_CONFIGURATIONS> block'
                  error stop
                endif
                if(.not.allocated(self%n_ion_core)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing N_ION_CORE in the <TWO_ELECTRON_CONFIGURATIONS> block'
                  error stop
                endif
                if(.not.allocated(self%nkin)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NKIN in the <TWO_ELECTRON_CONFIGURATIONS> block'
                  error stop
                endif
                !if(.not.allocated(self%nkout)) then
                !  if(myid==0) write(*,*) '*** ERROR reading input: Missing NKOUT in the <TWO_ELECTRON_CONFIGURATIONS> block'
                !  error stop
                !endif

                exit
              case default
                call invalid_label_in_block(label,block)
            end select 
          enddo !end reading lines in TWO_ELECTRON_CONFIGURATIONS block
        
        case('TWO_ELECTRON_STATES') !found the TWO_ELECTRON_STATES block
          two_el_states_read = .true.
          if(nitems>1) call heading_error(block) 
          do !read lines inside TWO_ELECTRON_STATES block
            call read_line(eof,nfile)
            if(eof) call eof_during_block(block)
            call read_label(label,block)
            select case (label)
              case('M_MAX')
                call readi(self%Mmax2el)
                if(self%Mmax2el < 0) call invalid_argument(label,block,self%Mmax2el)
                allocate(self%nst_2e(0:self%Mmax2el,-1:1,0:1))
                self%nst_2e = -1000

              case('NST_POS_PAR_SINGLET')
                call read_ALL_or_N(self%nst_2e(0,1,0),label,block)
                if(nitems >= 2 + self%Mmax2el+1) then
                  do M=1, self%Mmax2el
                    call read_ALL_or_N(self%nst_2e(M,1,0),label,block)
                  enddo
                else
                  call nitem_error(label,block)
                endif

              case('NST_NEG_PAR_SINGLET')
                call read_ALL_or_N(self%nst_2e(0,-1,0),label,block)
                if(nitems >= 2 + self%Mmax2el+1) then
                  do M=1, self%Mmax2el
                    call read_ALL_or_N(self%nst_2e(M,-1,0),label,block)
                  enddo
                else
                  call nitem_error(label,block)
                endif
              
              case('NST_POS_PAR_TRIPLET')
                call read_ALL_or_N(self%nst_2e(0,1,1),label,block)
                if(nitems >= 2 + self%Mmax2el+1) then
                  do M=1, self%Mmax2el
                    call read_ALL_or_N(self%nst_2e(M,1,1),label,block)
                  enddo
                else
                  call nitem_error(label,block)
                endif

              case('NST_NEG_PAR_TRIPLET')
                call read_ALL_or_N(self%nst_2e(0,-1,1),label,block)
                if(nitems >= 2 + self%Mmax2el+1) then
                  do M=1, self%Mmax2el
                    call read_ALL_or_N(self%nst_2e(M,-1,1),label,block)
                  enddo
                else
                  call nitem_error(label,block)
                endif

              case('END')
                !Check for missing input
                if(self%Mmax2el == -1000) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing M_MAX in the <TWO_ELECTRON_STATES> block'
                  error stop
                endif
                if(any(self%nst_2e(:,1,0) == -1000)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NST_POS_PAR_SINGLET in the <TWO_ELECTRON_STATES> block'
                  error stop
                endif
                if(any(self%nst_2e(:,-1,0) == -1000)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NST_NEG_PAR_SINGLET in the <TWO_ELECTRON_STATES> block'
                  error stop
                endif
                if(any(self%nst_2e(:,1,1) == -1000)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NST_POS_PAR_TRIPLET in the <TWO_ELECTRON_STATES> block'
                  error stop
                endif
                if(any(self%nst_2e(:,-1,1) == -1000)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NST_NEG_PAR_TRIPLET in the <TWO_ELECTRON_STATES> block'
                  error stop
                endif
                exit
              case default
                call invalid_label_in_block(label,block)
            end select 
          enddo !end reading lines in TWO_ELECTRON_STATES block
        
        case('BOUND_TWO_ELECTRON_STATES') !found the BOUND_TWO_ELECTRON_STATES block
          bound_two_el_states_read = .true.
          self%reduce_bound = .true.
          if(.not.two_el_states_read) then
            if(myid==0) write(*,*) '*** ERROR reading input: place BOUND_TWO_ELECTRON_STATES after TWO_ELECTRON_STATES'
            error stop
          endif
          if(nitems>1) call heading_error(block) 
          allocate(self%nst_2e_bound(0:self%Mmax2el,-1:1,0:1))
          self%nst_2e_bound = -1000
          do !read lines inside TWO_ELECTRON_STATES block
            call read_line(eof,nfile)
            if(eof) call eof_during_block(block)
            call read_label(label,block)
            select case (label)

              case('NST_POS_PAR_SINGLET')
                call read_ALL_or_N(self%nst_2e_bound(0,1,0),label,block)
                if(nitems >= 2 + self%Mmax2el+1) then
                  do M=1, self%Mmax2el
                    call read_ALL_or_N(self%nst_2e_bound(M,1,0),label,block)
                  enddo
                else
                  call nitem_error(label,block)
                endif

              case('NST_NEG_PAR_SINGLET')
                call read_ALL_or_N(self%nst_2e_bound(0,-1,0),label,block)
                if(nitems >= 2 + self%Mmax2el+1) then
                  do M=1, self%Mmax2el
                    call read_ALL_or_N(self%nst_2e_bound(M,-1,0),label,block)
                  enddo
                else
                  call nitem_error(label,block)
                endif
              
              case('NST_POS_PAR_TRIPLET')
                call read_ALL_or_N(self%nst_2e_bound(0,1,1),label,block)
                if(nitems >= 2 + self%Mmax2el+1) then
                  do M=1, self%Mmax2el
                    call read_ALL_or_N(self%nst_2e_bound(M,1,1),label,block)
                  enddo
                else
                  call nitem_error(label,block)
                endif

              case('NST_NEG_PAR_TRIPLET')
                call read_ALL_or_N(self%nst_2e_bound(0,-1,1),label,block)
                if(nitems >= 2 + self%Mmax2el+1) then
                  do M=1, self%Mmax2el
                    call read_ALL_or_N(self%nst_2e_bound(M,-1,1),label,block)
                  enddo
                else
                  call nitem_error(label,block)
                endif

              case('END')
                !Check for missing input
                if(self%Mmax2el == -1000) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing M_MAX in the <BOUND_TWO_ELECTRON_STATES> block'
                  error stop
                endif
                if(any(self%nst_2e_bound(:,1,0) == -1000)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NST_POS_PAR_SINGLET in the <BOUND_TWO_ELECTRON_STATES> block'
                  error stop
                endif
                if(any(self%nst_2e_bound(:,-1,0) == -1000)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NST_NEG_PAR_SINGLET in the <BOUND_TWO_ELECTRON_STATES> block'
                  error stop
                endif
                if(any(self%nst_2e_bound(:,1,1) == -1000)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NST_POS_PAR_TRIPLET in the <BOUND_TWO_ELECTRON_STATES> block'
                  error stop
                endif
                if(any(self%nst_2e_bound(:,-1,1) == -1000)) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NST_NEG_PAR_TRIPLET in the <BOUND_TWO_ELECTRON_STATES> block'
                  error stop
                endif
                exit
              case default
                call invalid_label_in_block(label,block)
            end select 
          enddo !end reading lines in BOUND_TWO_ELECTRON_STATES block
        
        case('RADIAL_GRID') !found the RADIAL_GRID block
          radial_grid_read = .true.
          if(nitems>1) call heading_error(block) 
          do !read lines inside RADIAL GRID block
            call read_line(eof,nfile)
            if(eof) call eof_during_block(block)
            call read_label(label,block)
            if(nitems /= 3 .and. trim(adjustl(label)) /= 'END') call nitem_error(label,block)
            select case (label)
              case('RMAX')
                call readf(self%rmax)
                if(self%rmax < 0) call invalid_argument(label,block,self%rmax)
              case('QMAX')
                call readf(self%qmax)
                if(self%qmax < 0) call invalid_argument(label,block,self%qmax)
              case('NDOUBLE')
                call readi(self%ndouble)
                if(self%ndouble < 0) call invalid_argument(label,block,self%ndouble)
              case('NPDBL')
                call readi(self%npdbl)
                if(self%npdbl < 0) call invalid_argument(label,block,self%npdbl)
              case('NPWAVE')
                call readi(self%npwave)
                if(self%npwave < 0) call invalid_argument(label,block,self%npwave)
              case('FORMCUT')
                call readf(self%formcut)
                if(self%formcut < 0) call invalid_argument(label,block,self%formcut)
              case('REGCUT')
                call readf(self%regcut)
                if(self%regcut < 0) call invalid_argument(label,block,self%regcut)
              case('EXPCUT')
                call readf(self%expcut)
                if(self%expcut < 0) call invalid_argument(label,block,self%expcut)
              case('END')
                !Check for missing input
                if(self%rmax == -1000d0) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing RMAX in the <RADIAL_GRID> block'
                  error stop
                endif
                if(self%qmax == -1000d0) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing QMAX in the <RADIAL_GRID> block'
                  error stop
                endif
                if(self%ndouble == -1000) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NDOUBLE in the <RADIAL_GRID> block'
                  error stop
                endif
                if(self%npdbl == -1000) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NPDBLE in the <RADIAL_GRID> block'
                  error stop
                endif
                if(self%npwave == -1000) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing NPWAVE in the <RADIAL_GRID> block'
                  error stop
                endif
                if(self%formcut == -1000d0) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing FORMCUT in the <RADIAL_GRID> block'
                  error stop
                endif
                if(self%regcut == -1000d0) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing REGCUT in the <RADIAL_GRID> block'
                  error stop
                endif
                if(self%expcut == -1000d0) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Missing EXPCUT in the <RADIAL_GRID> block'
                  error stop
                endif
                exit
              case default
                call invalid_label_in_block(label,block)
            end select 
          enddo !end reading lines in RADIAL_GRID block
        
        case('KGRID') !found the KGRID block
          if(nitems>1) call heading_error(block) 
          if(self%Lpmax == -1000) then
            if(myid==0) write(*,*) '*** ERROR reading input: <KGRID> block needs to appear after the <SCATTERING_SYSTEM> block'
            error stop
          endif
          Lbot = -1
          star = .false.
          star_prev = .false.
          allocate(self%nk(4,0:self%Lpmax), self%sk(4,0:self%Lpmax), self%NBND(0:self%Lpmax))
          self%nk(:,:) = -1000
          self%sk(:,:) = -1000.0
          self%lmaxkg = self%Lpmax
          LSW_prev = -1
          do !read lines inside KGRID block
            call read_line(eof,nfile)
            if(eof) call eof_during_block(block)
            call readu(label)
            if(trim(adjustl(label)) == 'END') exit
            kgrid_read = .true.
            call reread(-1)
            if(nitems /= 10) then
              if(myid==0) call report('wrong number of items in kgrid')
              error stop
            endif
            call readu(tempchar)
            if(trim(adjustl(tempchar)) == '*') then
              if(star_prev) then
                if(myid==0) call report('cannot have two "*" lines in a row in <KGRID> block')
                error stop
              endif
              star = .true.
              LSW = self%Lpmax
            else
              star = .false.
              call reread(-1)
              call readi(LSW)
              write(tempchar,'(I0)') LSW
            endif
            !Ltt = min(LSW,self%Lpmax)
            if(LSW > self%Lpmax) then
              if(myid==0) then
                call report('L in <KGRID> block is greater than projectile LMAX')
                write(*,'("     L=",I0,", while projectile LMAX=",I0)') LSW, self%Lpmax
                write(*,'("     (if you''re trying to fill the kgrid out up to LMAX try using the ""*"" notation)")')
              endif
              error stop
            endif
            Ltt = LSW
            if(Ltt <= Lbot) then
              if(myid==0) call report('duplicate L values in <KGRID> block')
              error stop
            elseif(Ltt < Lbot) then
              if(myid==0) call report('L values out of order in <KGRID> block')
              error stop
            endif
            if(star) then
              Lbot = LSW_prev + 1
            else
              Lbot = LSW
            endif
            if(star_prev .and. Ltt < self%Lpmax) then
              self%nk(:,Ltt+1:) = -1000
              self%sk(:,Ltt+1:) = -1000.0
            endif

            call readi(NBNDt)
            
            self%NBND(Lbot:Ltt) = NBNDt
            do i=1, 4
              call readi(nkt)
              if(nkt <= 0 .and. i < 4) then
                if(myid==0) call report('invalid number of points in kgrid interval')
                error stop
              endif
              self%nk(i,Lbot:Ltt) = nkt
              call readf(skt)
              if(skt <= 0) then
                if(myid==0) call report('invalid boundary/width of kgrid interval')
                error stop
              elseif(i == 2) then
                if(skt <= self%sk(i-1,Lbot)) then
                  if(myid==0) call report('invalid boundary/width of kgrid interval')
                  error stop
                endif
              endif
              self%sk(i,Lbot:Ltt) = skt
            enddo !i
              
            if(star) then
              star_prev = .true.
              star = .false.
            else
              star_prev = .false.
            endif

            LSW_prev = LSW

          enddo !end reading lines in KGRID block
          if(any(self%nk(1,:) == -1000)) then
            if(myid==0) then
              write(*,*) '*** ERROR reading input: <KGRID> block is missing the following L values:'
              write(*,'(25X)',advance='no')
              do i=0, self%Lpmax
                if(self%nk(1,i) == -1000) write(*,'(I0,2X)',advance='no') i
              enddo
              write(*,*)
            endif
            error stop
          endif
        
        case('CALCULATION_OPTIONS') !found the CALCULATION_OPTIONS block
          calculation_options_read = .true.
          if(nitems>1) call heading_error(block) 
          do !read lines inside CALCULATIONS_OPTIONS block
            call read_line(eof,nfile)
            if(eof) call eof_during_block(block)
            call read_label(label,block)
            if(nitems /= 3 .and. trim(adjustl(label)) /= 'END') call nitem_error(label,block)
            select case (label)
              case('LTMAX')
                call readi(self%ltmax)
                if(self%ltmax < 0) call invalid_argument(label,block,self%ltmax)
              case('BORN_COMPLETION')
                call switch(born,label,block)
                if(born) then
                  self%iborn = 1
                else
                  self%iborn = 0
                endif
              case('UBA')
                call switch(self%UBA,label,block)
              case('PWBORN')
                call switch(self%pwborn,label,block)
              case('THETA')
                call readf(self%theta)
              case('EXCHANGE')
                call switch(exchange,label,block)
              case('REARRANGE')
                call switch(inc,label,block)
                if(inc) then
                  self%inc = 1
                else
                  self%inc = 0
                endif
              case('MPI_LOAD_BALANCING')
                call switch(self%load_balance,label,block)
                if(ntasks == 1) self%load_balance = .false.
              case('READ_CHANNEL_TIMING')
                call switch(self%read_channel_timing,label,block)
              case('DISTORTING_POTENTIAL')
                call switch(distpot,label,block)
                if(distpot) then
                  self%idistpot = 1
                else
                  self%idistpot = 0
                endif
                if(self%calculation_type == 2 .and. self%idistpot > 0) then
                  if(myid==0) call report('distorting potential not coded in Spheroidal implementation')
                  error stop
                endif
              case('DISTPOT_LMAX')
                call readi(self%ldw)
                if(self%ldw < 0) call invalid_argument(label,block,self%ldw)
              case('PROJ_BOUNDSTATE_NMAX')
                call readi(self%ndw)
                if(self%ndw < 0) call invalid_argument(label,block,self%ndw)
              case('PROJ_BOUNDSTATE_ALPHA')
                call readf(self%aldw)
                if(self%aldw <= 0.0d0) call invalid_argument(label,block,self%aldw)
              case('POSITRONIUM_FORMATION')
                call switch(templ,label,block)
                if(templ) then
                  self%l_Ps = 1
                else
                  self%l_Ps = -1
                endif
              case('CI_MIN')
                call readf(self%ci_min)
                if(self%ci_min <= 0.0d0) then
                  if(myid==0) call report('CI_MIN must be positive')
                  error stop
                endif
              case('OPTICAL_POTENTIAL')
                call switch(self%optical,label,block)
              case('WEAK_COUPLING')
                call switch(self%weak_coupling,label,block)
              case('WEAK_COUPLING_EXCHANGE')
                call switch(self%weak_coupling_exchange,label,block)
              case('P_SPACE')
                call read_ALL_or_N(self%NP,label,block)
              !case('MAX_BOUND_STATES')
              !  call readu(tempchar)
              !  if(trim(adjustl(tempchar)) /= 'OFF') then
              !    call reread(-1)
              !    call readi(self%max_bound)
              !    if(self%max_bound < 1) then
              !      if(myid==0) call report('MAX_BOUND_STATES must be greater than 0')
              !      error stop
              !    endif
              !    self%reduce_bound = .true.
              !  endif
              case('INTEGRATION_WEIGHTS')
                call readu(tempchar)
                select case (tempchar) 
                  case ('SIMPSON')
                    self%iweight = 0
                  case ('BOOL','BOOLE')
                    self%iweight = 1
                  case default
                    if(myid==0) call report ('Invalid integration weight')
                    error stop
                end select
              case('SPECTRO_CORE')
                call switch(self%spectro_core,label,block)
              case('SF_TOLERANCE')
                call readf(self%SF_tol)
                if(self%SF_tol < 0 .or. self%SF_tol > 1) then
                  if(myid==0) call report('SF_tol must be between 0 and 1')
                  error stop
                endif
              case('CORE_POLARISABILITY')
                call readf(self%corep)
                if(self%corep < 0.0d0) then
                  if(myid==0) call report('Core polarisation must be >= 0.0')
                  error stop
                endif
              case('CORE_POTENTIAL_CUTOFF')
                call readf(self%r0)
                if(self%r0 <= 0.0d0) then
                  if(myid==0) call report('Core potential cutoff must be > 0.0')
                  error stop
                endif
              case('NO_SECOND_SPIN_CHANNEL')
                call switch(self%no_second_spin_channel,label,block)
              case('SKIP_VMATRIX')
                call switch(self%skip_vmat,label,block)
              case('END')
                if(self%optical .and. self%NP == -1) then
                  if(myid==0) write(*,*) '*** ERROR reading input: you have OPTICAL ON but have not specified PSPACE < Nmax'
                  error stop
                endif
                if(self%weak_coupling .and. self%NP == -1) then
                  if(myid==0) write(*,*) '*** ERROR reading input: you have WEAK_COUPLING ON but have not specified PSPACE < Nmax'
                  error stop
                endif
                if(self%weak_coupling_exchange .and. self%NP == -1) then
                  if(myid==0) write(*,*) '*** ERROR reading input: you have WEAK_COUPLING ON but have not specified PSPACE < Nmax'
                  error stop
                endif
                if(self%weak_coupling .and. self%weak_coupling_exchange) then
                  if(myid==0) write(*,*) '*** ERROR reading input: Incompatible options WEAK_COUPLING and WEAK_COUPLING_EXCHANGE' 
                  error stop
                endif
                self%exchange = exchange
                
                !ifirst: used to control first order calcs AND determine if exchange is included
                !Now: we use data_in%UBA and data_in%pwborn to control two different types of first-order calcs
                !     and data_in%exchange to control exchange (so exchange can be turned on for FO calcs if desired
               
                
                if(self%UBA .and. self%pwborn) then
                  if(myid==0) write(*,*) '*** ERROR reading input: UBA and PWBORN are mutually exclusive'
                  error stop
                endif
                
                if ( self%exchange .AND. .not.(self%UBA.or.self%pwborn) .AND. self%theta /= 0d0 .AND. self%Zproj== -1d0 ) self%non_uniq = .TRUE.

                
                !if(exchange .and. .not.UBA) then
                !  self%ifirst = 1
                !elseif(UBA .and. .not.exchange) then
                !  self%ifirst = -1
                !elseif(.not.exchange .and. .not.UBA) then
                !  self%ifirst = 0
                !else
                !  if(myid==0) write(*,*) '*** ERROR reading input: cannot have UBA with exchange'
                !  error stop
                !endif
                !if(exchange .and. self%pwborn.and..false.) then
                !  if(myid==0) write(*,*) '*** ERROR reading input: cannot have PWBORN with exchange'
                !  error stop
                !endif
                !if(UBA .and. self%pwborn) then
                !  if(myid==0) write(*,*) '*** ERROR reading input: UBA and PWBORN are mutually exclusive'
                !  error stop
                !endif
                !if(self%pwborn) self%ifirst = -1 !ifirst = -1 sets up first-order calculation but pwborn=true stops UBA
                !if(self%iAnalitBorb == 1) self%iBorn = 1
    
                !if ( self%ifirst == 1 .AND. self%theta /= 0d0 .AND. self%Zproj== -1d0 ) self%non_uniq = .TRUE.
                
                exit
              case default
                call invalid_label_in_block(label,block)
            end select 
          enddo !end reading lines in CALCULATION_OPTIONS block

          case('PSEUDO_POTENTIAL')
            read_B0 = .false.; read_B1 = .false.; read_B2 = .false.
            read_Beta0 = .false.; read_Beta1 = .false.; read_Beta2 = .false.
            self%core_energy = 0.0d0
            do !read lines inside PSEUDO_POTENTIAL block
              call read_line(eof,nfile)
              if(eof) call eof_during_block(block)
              call read_label(label,block)
              if(nitems /= 3 .and. trim(adjustl(label)) /= 'END') call nitem_error(label,block)
              select case(label)
                case('USE_PSEUDO_POTENTIAL')
                  call switch(self%pseudo_pot,block,block)
                case('B0')
                  call readf(self%pseudo_pot_B(0))
                  read_B0 = .true.
                case('BETA0')
                  call readf(self%pseudo_pot_Beta(0))
                  read_Beta0 = .true.
                case('B1')
                  call readf(self%pseudo_pot_B(1))
                  read_B1 = .true.
                case('BETA1')
                  call readf(self%pseudo_pot_Beta(1))
                  read_Beta1 = .true.
                case('B2')
                  call readf(self%pseudo_pot_B(2))
                  read_B2 = .true.
                case('BETA2')
                  call readf(self%pseudo_pot_Beta(2))
                  read_Beta2 = .true.
                case('CORE_ENERGY')
                  call readf(self%core_energy)
                case('END')
                  if(.not. read_B0) then
                    if(myid==0) write(*,*) '*** ERROR reading input: B0 missing from PSEUDO_POTENTIAL_BLOCK'
                    error stop
                  elseif(.not. read_Beta0) then
                    if(myid==0) write(*,*) '*** ERROR reading input: BETA0 missing from PSEUDO_POTENTIAL_BLOCK'
                    error stop
                  elseif(.not. read_B1) then
                    if(myid==0) write(*,*) '*** ERROR reading input: B1 missing from PSEUDO_POTENTIAL_BLOCK'
                    error stop
                  elseif(.not. read_Beta1) then
                    if(myid==0) write(*,*) '*** ERROR reading input: BETA1 missing from PSEUDO_POTENTIAL_BLOCK'
                    error stop
                  elseif(.not. read_B2) then
                    if(myid==0) write(*,*) '*** ERROR reading input: B2 missing from PSEUDO_POTENTIAL_BLOCK'
                    error stop
                  elseif(.not. read_Beta2) then
                    if(myid==0) write(*,*) '*** ERROR reading input: BETA2 missing from PSEUDO_POTENTIAL_BLOCK'
                    error stop
                  endif

                  exit
                case default
                  call invalid_label_in_block(label,block)
              end select

            enddo


        case default 
          if(myid==0) write(*,*) '*** ERROR reading input: Invalid block "'//trim(adjustl(block))//'"'
          error stop
      end select

    enddo !end main loop over data.in file

    close(nfile)

    !Check that all nessesary input blocks have been provided. 
    if(.not.scattering_system_read) then
      if(myid==0) write(*,*) '*** ERROR reading input: SCATTERING_SYSTEM block missing'
      error stop
    endif
    if(.not.radial_grid_read) then
      if(myid==0) write(*,*) '*** ERROR reading input: RADIAL_GRID block missing'
      error stop
    endif
    if(.not.calculation_options_read) then
      if(myid==0) write(*,*) '*** ERROR reading input: CALCULATION_OPTIONS block missing'
      error stop
    endif
    
    if(.not.first_basis_read) then
      if(myid==0) write(*,*) '*** ERROR reading input: FIRST_DIAGONALISATION_BASIS block missing'
      error stop
    endif
    if(.not.core_first_basis_read) core_basis_1e = basis_1e

    if(.not.one_el_states_read) then
      if(myid==0) write(*,*) '*** ERROR reading input: ONE_ELECTRON_STATES block missing'
      error stop
    endif
    if(.not.second_basis_read) write(*,*) '*** WARNING: Second basis not specified: assuming you want to use only one-electron states for second diagonalisation'
    if(.not.core_second_basis_read) then
      core_basis_2e = basis_2e
      core_sub_basis_2e = sub_basis_2e
    endif
    
    if(.not.self%hlike) then
      if(.not.two_el_configs_read) then
        if(myid==0) write(*,*) '*** ERROR reading input: TWO_ELECTRON_CONFIGURATIONS block missing'
        error stop
      endif
      if(.not.two_el_states_read)then
        if(myid==0) write(*,*) '*** ERROR reading input: TWO_ELECTRON_STATES block missing'
        error stop
      endif
    endif

    if(scattering .and. .not.kgrid_read) then
      if(myid==0) write(*,*) '*** ERROR reading input: KGRID block missing'
      error stop
    endif
    
    if(.not.output_read) then
      if(myid==0) write(*,*) '*** ERROR reading input: OUTPUT block missing'
      error stop
    endif
    
    !-=-= Misc input validation which depends on values in multiple blocks so we do it at the very end -=-=
    
    if(self%use_MSC == 0 .and. self%calculation_type == 0) then 
      print*, 'Including no one-el states not possible in spherical - gotta trick it'
      dataMSC%MSC_nconfig = 1
      if(.not.first_basis_read) then
        if(myid==0) write(*,*) '*** ERROR reading input: place <FIRST_DIAGONALISATION_BASIS> before SECOND' 
        error stop
      endif
      basis_1e%latop = 0
      if(allocated(basis_1e%nps)) deallocate(basis_1e%nps)
      if(allocated(basis_1e%alpha)) deallocate(basis_1e%alpha)
      allocate(basis_1e%nps(0:0), basis_1e%alpha(0:0))
      basis_1e%nps(0) = 1
      basis_1e%alpha(0) = basis_2e%alpha(0)
      self%nst = 0
      self%nst(0,1) = 1
    endif

    if(.not.allocated(sub_basis_2e%nps)) then
      sub_basis_2e%labot=0
      sub_basis_2e%latop=0
      allocate(sub_basis_2e%nps(0:0))
      allocate(sub_basis_2e%alpha(0:0))
      sub_basis_2e%nps = 0
      sub_basis_2e%alpha = 0.0d0
    endif
          
    if(self%calculation_type == 0) then !spherical stores basis differently (gotta change this at some point)
      self%labot = 0
      self%latop = basis_1e%latop
      allocate(self%alpha(0:self%latop))
      allocate(self%nps(0:self%latop))
      self%nps = basis_1e%nps
      self%alpha = basis_1e%alpha
      deallocate(basis_1e%nps, basis_1e%alpha) !basis_1e won't be referred to in spherical code
      dataMSC%labot = 0
      dataMSC%latop = basis_2e%latop
      allocate(dataMSC%alpha(0:basis_2e%latop))
      allocate(dataMSC%nps(0:basis_2e%latop))
      dataMSC%nps = basis_2e%nps
      dataMSC%alpha = basis_2e%alpha
      if(self%use_MSC == -1) then
        dataMSC%MSC_nconfig = 0
      elseif(self%use_MSC > 0) then
        dataMSC%MSC_nconfig = self%use_MSC !case with use_MSC = 0 attended to above
      endif
      
      !Liam added .true. here so alpha_nl array is defined when using natural orbitals
      if (.true. .or. ABS(dataMSC%MSC_nconfig) /= 0)  allocate(dataMSC%alpha_nl(MAXVAL(dataMSC%nps(dataMSC%labot:dataMSC%latop)),dataMSC%labot:dataMSC%latop))
      
      dataMSC%labot_diff = sub_basis_2e%labot
      dataMSC%latop_diff = sub_basis_2e%latop
      allocate(dataMSC%nps_diff(dataMSC%labot_diff:dataMSC%latop_diff), dataMSC%alpha_diff(dataMSC%labot_diff:dataMSC%latop_diff))
      dataMSC%nps_diff = sub_basis_2e%nps
      dataMSC%alpha_diff = sub_basis_2e%alpha
      if (ABS(dataMSC%MSC_nconfig) /= 0) then
         do l = dataMSC%labot, dataMSC%latop
            dataMSC%alpha_nl(:,l) = dataMSC%alpha(l)
            if ( self%use_sub > 0 ) then ! overwrite alpha with different alpha's
               if ( l>=dataMSC%labot_diff .AND. l<=dataMSC%latop_diff ) dataMSC%alpha_nl(1:dataMSC%nps_diff(l),l) = dataMSC%alpha_diff(l)
            end if
         end do
      end if
    endif

    if(.not.self%good_parity .and. (par_start /= 0 .or. par_stop /= 0)) then
      if(myid==0) write(*,*) '*** ERROR reading input: input to the &
        &<PARITY> keyword in the <SCATTERING_SYSTEM> block not consistent with lack of parity for this system'
      error stop
    endif

    if(self%good_parity) then
      self%nst(:,0) = 0
      self%nst_2e(:,0,:) = 0
    else
      self%nst(:,0) = self%nst(:,1) + self%nst(:,-1)
      self%nst(:,1) = 0
      self%nst(:,-1) = 0
      self%nst_2e(:,0,:) = self%nst_2e(:,1,:) + self%nst_2e(:,-1,:)
      self%nst_2e(:,1,:) = 0
      self%nst_2e(:,-1,:) = 0
      where(self%nst < 0) self%nst = -1  !To keep functionality that "-1" = "all states"
      where(self%nst_2e < 0) self%nst_2e = -1
    endif
    

    if(self%use_sub > sum(self%nst) .and. all(self%nst >= 0)) then
      if(myid==0) write(*,*) '*** ERROR reading input: requested more 1-el states in the second basis than kept from the&
        & first diagonalisation'
      error stop
    endif

    if(self%iosc == 1) self%iosc = self%nent !reducing number of input parameters, calculate osc strengths for all entrance states

    self%lbnd_max = -1
    do i = 0, self%Lmaxkg
      if ( self%NBND(i) > 0 ) then
        self%lbnd_max = i
      end if
    end do

    if (self%idistpot <= 0) self%Ldw = -1
    
    if(self%Ldw .gt. self%Lmaxkg) then
       self%Ldw = self%Lmaxkg 
       if(myid.eq.0) print*,'*** WARNING: You have too large value of Ldw, changing it to maximum L of projectile:'
       if(myid.eq.0) print*,'             Ldw = ',self%Ldw
    endif

    if(self%l_Ps >= 0) call readin_Ps(self,myid+1)


    !Print calculation details
    if(myid==0) then
      
      write(*,'(A)') '##################### CALCULATION DETAILS ######################'
      if(ntasks > 1) then
        write(*,*) '[+] Parallelism mode: Hybrid MPI + OpenMP'
      elseif(nomp > 1) then
        write(*,*) '[+] Parallelism mode: OpenMP'
      else
        write(*,*) '[+] Parallelism mode: None'
      endif

      write(*,'("     > Number of tasks: ",I0)') ntasks
      if(ntasks == 1) then
        write(*,'("     > Number of OpenMP threads: ",I0)') nomp
      else
        write(*,'("     > Number of OpenMP threads per task: ",I0)') nomp
      endif
      write(*,'("     > Total number of processors: ",I0)') nomp*ntasks
      write(*,*)
    
      if(en_eV < 10.0d0) then
        form='F3.1'
      elseif(en_eV < 100.0d0) then
        form='F4.1'
      elseif(en_eV < 1000.0d0) then
        form='F5.1'
      else
        form='F6.1'
      endif
  
      if(self%Zproj == -1) then
        projectile = 'Electron'
      elseif(self%Zproj == 1) then
        projectile = 'Positron'
      else
        error stop '***ERROR: Invalid projectile charge'
      endif
        
      form='(" [+] Scattering system: ",A," (",'//form//'," eV) --> ",A)'
      write(*,form) projectile, en_eV, self%target
      write(*,'("     > Internuclear separation R = ",F9.5)') self%Rd
      write(*,'("     > Target charges: Z1=",F3.1,", Z2=",F3.1,", Zasym=",F3.1)') self%Z1, self%Z2, self%Zasym
      n_elec = self%Z1 + self%Z2 - self%Zasym
      if(self%N_core_el > 0) then
        write(*,'("       number of electrons: ",i0," (",i0," active and ",i0," in the inert core)")') n_elec, n_elec-self%n_core_el, self%n_core_el
        if(self%corep > 0.0d0) then
          write(*,'("     > Adding core polarisation potential with static dipole polarisability alpha=",F5.3)') self%corep
          write(*,'("       and cutoff radius r0=",F5.3)') self%r0
        else
          write(*,'("       WARNING: Core polarisation potential parameters not specified - not including it")')
        endif
      else
        write(*,'("       number of electrons: ",i0)') n_elec
      endif
      if(n_elec-self%N_core_el /= 1 .and. n_elec-self%N_core_el /= 2) error stop '***ERROR: Zasym invalid, or this target not coded yet'
      write(*,*) '    > Hlike: ', self%hlike
      write(*,*) '    > Good parity: ', self%good_parity
      write(*,*)

      if(self%calculation_type == 0) then
         write(*,'(" [+] Calculation type : 0 - Spherical with nonorthogonal (2l+1 type) radial basis functions")')
      elseif(self%calculation_type == 1) then
         write(*,'(" [+] Calculation type : 1 - Spherical with orthogonal (2l+2 type) radial basis functions")')
         error stop '***ERROR: This option is not coded'
      elseif (self%calculation_type == 2) then
         write(*,'(" [+] Calculation type : 2 - Spheroidal with orthogonal (m type) Hylleraas radial basis functions")') 
      elseif (self%calculation_type == 3) then
         write(*,'(" [+] Calculation type : 3 - Spheroidal with non-orthogonal (m type) Hylleraas radial basis functions")')
      endif

      if(self%calculation_type == 0 .or. self%calculation_type == 1) then !spherical - origin is specified
        write(*,'("     > Spherical coordinate-system origin: ",F5.3)') self%origin
      endif
      write(*,*)

      if (self%calculation_mode == 0) then
         write(*,*) '[+] Calculation mode : 0 - structure only'
      elseif (self%calculation_mode == 1) then
         write(*,*) '[+] Calculation mode : 1 - scattering, no DCS'
      elseif (self%calculation_mode == 2) then
         write(*,*) '[+] Calculation mode : 2 - scattering with DCS and TV files'
      elseif (self%calculation_mode == 3) then
         write(*,*) '[+] Calculation mode : 3 - scattering with TV files (no DCS)'
      end if
  
      if(self%calculation_mode > 0 .and. self%Zproj == 1) then !positron scattering

        if(self%calculation_mode > 0) then
          if(self%l_Ps < 0) then
            write(*,'("[+] Single-centre positron-scattering calculations, l_Ps =",4i0)') self%l_Ps
          else
            write(*,'("[+] Two-centre positron-scattering calculations, l_Ps =",4i0)') self%l_Ps
          endif
        endif
      endif !positron scattering
      write(*,*)
      
      if(self%calculation_type == 0) then !spherical
        if(self%N_core_el > 0) then
          write(*,'(" [+] Spherical Lageurre basis: (used for first diagonalisation for inert core structure)")')
          write(*,'("     > labot,latop: ",2I5)') core_basis_1e%labot, core_basis_1e%latop
          write(tempchar,'(I0)') core_basis_1e%latop-core_basis_1e%labot + 1
          write(*,'("     > nps(l),alpha(l): ",'//tempchar//'("(",I0,",",F7.4,"), "),"(",I0,",",F7.4,")")') (core_basis_1e%nps(l),core_basis_1e%alpha(l), l=core_basis_1e%labot,core_basis_1e%latop)
          write(*,*)
        endif
        
        if(self%hlike .and. dataMSC%MSC_nconfig == 0) then
          write(*,'(" [+] Spherical Lageurre basis: (used for final one-electron structure calculations)")')
        else
          write(*,'(" [+] Spherical Lageurre basis: (used to generate accurate MOs for second diagonalisation)")')
        endif
        write(*,'("     > labot,latop: ",2I5)') self%labot, self%latop
        write(tempchar,'(I0)') self%latop-self%labot + 1
        write(*,'("     > nps(l),alpha(l): ",'//tempchar//'("(",I0,",",F7.4,"), "),"(",I0,",",F7.4,")")') (self%nps(l),self%alpha(l), l=self%labot,self%latop)
        write(*,*)
         
        if(self%hlike .and. dataMSC%MSC_nconfig > 0) then
          write(*,*)   "[+] Second spherical basis: (Hybrid MO + Laguerre, used for final one-electron structure calculations)"
        elseif(dataMSC%MSC_nconfig > 0) then
          write(*,*)   "[+] Second spherical basis: (Hybrid MO + Laguerre, used for final one- and two-electron structure calculations)"
        endif

        if(dataMSC%MSC_nconfig == 0 .and. .not.self%hlike) then
          write(*,*) "    > MSC_nconfig = 0: basis will be formed only from the MOs saved from the first diagonalisation"
        elseif(dataMSC%MSC_nconfig > 0) then
          write(*,'("     > labot,latop: ",2I5)') dataMSC%labot, dataMSC%latop
          write(tempchar,'(I0)') dataMSC%latop-dataMSC%labot + 1
          write(*,'("     > nps(l),alpha(l): ",'//tempchar//'("(",I0,",",F7.4,"), "),"(",I0,",",F7.4,")")') (dataMSC%nps(l),dataMSC%alpha(l), l=dataMSC%labot,dataMSC%latop)
          write(*,'("     > MSC_nconfig = ",I0,": the first ",I0," MOs saved from first diagonalisation will be subtituted into the second basis")') dataMSC%MSC_nconfig, dataMSC%MSC_nconfig
        endif !nconfig
       
        if(self%use_sub > 0) then
          write(*,'("     > use_sub = ",I0,": functions from the following basis will replace Laguerres (not MOs) in the Hybrid basis")') self%use_sub
          write(tempchar,'(I0)') dataMSC%latop_diff-dataMSC%labot_diff + 1
          write(*,'("     > nps(l),alpha(l): ",'//tempchar//'("(",I0,",",F7.4,"), "),"(",I0,",",F7.4,")")') (dataMSC%nps_diff(l),dataMSC%alpha_diff(l), l=dataMSC%labot_diff,dataMSC%latop_diff)
        endif
      
      elseif(self%calculation_type == 2) then !spheroidal
        if(self%hlike) then
          write(*,'(" [+] Spheroidal Hylleraas basis:")')
        else !two-electron target
          write(*,'(" [+] First spheroidal Hylleraas basis: (used for one-electron structure")')
          write(*,'("                                        and to generate accurate MOs for second diagonalisation)")')
        endif
        write(*,'("     > labot,latop: ",2I5)') basis_1e%labot, basis_1e%latop
        write(tempchar,'(I0)') basis_1e%latop-basis_1e%labot+1
        write(*,'("     > nps(l),alpha(m): ",'//tempchar//'("(",I0,",",F7.4,"), "),"(",I0,",",F7.4,")")') (basis_1e%nps(l),basis_1e%alpha(l), l=basis_1e%labot,basis_1e%latop)
        write(*,*)
        if (.not.self%hlike) then !second basis only used for two-electron structure
          if(self%use_MSC == 0) then
            write(*,*)   "[+] Second spheroidal basis: (Pure Hylleraas, used for two-electron structure calculations)"
          else
            write(*,*)   "[+] Second spheroidal basis: (Hybrid MO + Hylleraas, used for two-electron structure calculations)"
          endif
          write(*,'("     > labot,latop: ",2I5)') basis_2e%labot, basis_2e%latop
          write(tempchar,'(I0)') basis_2e%latop-basis_2e%labot+1
          write(*,'("     > nps(l),alpha(m): ",'//tempchar//'("(",I0,",",F7.4,"), "),"(",I0,",",F7.4,")")') (basis_2e%nps(l),basis_2e%alpha(l), l=basis_2e%labot,basis_2e%latop)
          if(self%use_MSC > 0) then
            write(*,'("     > MSC_nconfig = ",I0,": the first ",I0," MOs saved from first diagonalisation will be subtituted into the second basis")') self%use_MSC, self%use_MSC
          endif
        end if
       
        if (.not.self%hlike) then !sub basis only used for two-electron structure
          write(*,'("     > use_sub = ",I0,": functions from the following basis will replace Laguerres (not MOs) in the Hybrid basis")') self%use_sub
          write(tempchar,'(I0)') sub_basis_2e%latop-sub_basis_2e%labot+1
          write(*,'("     > nps(l),alpha(m): ",'//tempchar//'("(",I0,",",F7.4,"), "),"(",I0,",",F7.4,")")') (sub_basis_2e%nps(l),sub_basis_2e%alpha(l), l=sub_basis_2e%labot,sub_basis_2e%latop)
        endif 
      
      endif !spherical / spheroidal
     
      if(self%calculation_mode > 0) then
        write(*,*)
        write(*,'("[+] Mtot_start, Mtot_stop, ipar, nent, iborn: ",6I5)') self%Mtot_start, self%Mtot_stop, self%ipar, self%nent, self%iborn
      endif
        
      write(*,*)
      write(*,'(" [+] Radial grid parameters:")')
      write(*,'("     > rmax,qmax,ndouble,npdbl,npwave,ltmax: ",2F15.5,4I5)') self%rmax, self%qmax, self%ndouble, self%npdbl, self%npwave, self%ltmax
      write(*,'("     > formcut,regcut,expcut: ",3E10.2)') self%formcut, self%regcut, self%expcut
      if(self%iweight == 0) then
        write(*,'("     > Integration weights: SIMPSON")')
      else
        write(*,'("     > Integration weights: BOOLE")')
      endif
      write(*,*)
      
      write(*,'(" [+] iSlater: ",I5,", iosc=",I5,", irec=",I5 )') self%iSlater,self%iosc, self%irec
      write(*,*)
      
      write(*,'(" [+] One-electron target states to be saved:")')
      write(*,'("     > Mt_min, Mt_max: ",I0,2X,I0 )') self%Mt_min, self%Mt_max
      if(self%good_parity) then
        write(*,'("     > nst(m,+1): ",20(I0,2X))') (self%nst(m,1), m=self%Mt_min,self%Mt_max)
        write(*,'("     > nst(m,-1): ",20(I0,2X))') (self%nst(m,-1), m=self%Mt_min,self%Mt_max)
      else
        write(*,'("     > nst(m): ",20(I0,2X))') (self%nst(m,0), m=self%Mt_min,self%Mt_max)
      endif
      write(*,*)
     
      if(self%calculation_mode > 0) then
        write(*,'(" [+] idistpot, LDW, nDW, alDW: ",3I5,F10.4)') self%idistpot, self%Ldw, self%ndw, self%aldw
        write(*,*)
        write(*,'(" [+] ifirst,theta: ",I5,F10.4)') self%ifirst, self%theta
        if(self%Zproj == -1) write(*,'("     >  Solving non-uniqueness: ", L)') self%non_uniq
        write(*,*)
        write(*,'(" [+] inc: ",I5)') self%inc
        write(*,*)
        write(*,'(" [+] Lpmax, orient, thrad, phrad: ",2I5,2F10.4)') self%Lpmax, self%orient, self%thrad, self%phrad
        write(*,*)
        
        write(*,*) '[+] Kgrid and bound projectile-state info:'
        do Ltt = 0,self%Lpmax
          write(*,'("     > LSW,NBND, nk(1:4),sk(1:4): ",2I5,4(I5,F7.2))') Ltt, self%NBND(Ltt), (self%nk(i,Ltt),self%sk(i,Ltt), i=1,4)
        enddo
      
      endif

      write(*,'(A)') '################### END CALCULATION DETAILS ####################'
      write(*,*)
    
    endif !myid==0



  end subroutine readin

  subroutine readin_old( self, basis_1e,basis_2e,sub_basis_2e, nfile,iwrite )
    use target_data
    use MPI_module
    implicit none
    type(input):: self
    type(basis_input) :: basis_1e,basis_2e,sub_basis_2e

    integer, intent(in):: nfile
    integer:: nfile_Ps
    integer, intent(in):: iwrite  ! if =1 write to screen the read data, =0 no writing to screen

    logical:: ex
    integer:: iostat_data, icheck, LSW_old
    integer:: l, m, mpiwrite, latop,labot
    real(dpf):: en_eV, eV
    integer:: Ltt, LSW, i, mint, NBNDt, Lmaxkgt, Lbot
    integer, dimension(:), allocatable:: nkt
    real, dimension(:), allocatable:: skt
!!$ Different alpha n,l 
    integer:: n, n_elec 
    character(len=50) :: tempchar
    character(len=8) :: projectile
    character(len=:), allocatable :: form

    logical :: origin_not_specified = .false.
    logical :: l_ps_not_specified = .false.

    self%eV=27.21138505d0
    self%eV_old=27.2116d0

    eV = self%eV  !  27.2116
    
    self%CI_min = 1.0D-6
    
    !New features disabled for old input files
    self%print_1el_basis = .false.
    self%print_2el_config = .false.
    self%print_CI = .false.
    self%print_dipole = .false.
    self%num_nat_orb = 0
    self%only_nat_orbs = .false.

    ! MPI print to screen
    if (iwrite == 1 .AND. myid==0 ) then
       mpiwrite = 1
    else
       mpiwrite = 0
    end if

    !
    if(nfile .eq. 10) then
       inquire(file='data.in',exist=ex)
       if(.not. ex) then
          if(myid.eq.0) print*,'File data.in does not exists. STOP'
          STOP
       end if
       open(nfile,file='data.in',iostat=iostat_data)
       if(iostat_data.ne.0) then
          if(myid.eq.0) print*, '********  File  canot open file data.in'
          stop
       end if
    endif
    
    if(mpiwrite==1) then
      write(*,'(A)') '##################### CALCULATION DETAILS ######################'
      if(ntasks > 1) then
        write(*,*) '[+] Parallelism mode: Hybrid MPI + OpenMP'
      elseif(nomp > 1) then
        write(*,*) '[+] Parallelism mode: OpenMP'
      else
        write(*,*) '[+] Parallelism mode: None'
      endif

      write(*,'("     > Number of tasks: ",I0)') ntasks
      if(ntasks == 1) then
        write(*,'("     > Number of OpenMP threads: ",I0)') nomp
      else
        write(*,'("     > Number of OpenMP threads per task: ",I0)') nomp
      endif
      write(*,'("     > Total number of processors: ",I0)') nomp*ntasks
      write(*,*)
    endif

    read(nfile,'(A10)') self%target

    en_eV = 0d0
    read(nfile,*) en_eV
    self%energy = en_eV / eV
    
    self%calculation_type = -1
    self%calculation_mode = 2   ! Default to structure + scattering + DCS.
    self%l_Ps = -1
    read(nfile,*,iostat=iostat_data) self%calculation_type, self%calculation_mode, self%l_Ps

    if(iostat_data.ne.0 ) l_ps_not_specified = .true.

    read(nfile,*,iostat=iostat_data) self%Rd, self%origin
    if (iostat_data /= 0) then
       if (self%calculation_type==0) origin_not_specified = .true.
       self%origin = 0.5d0
    end if
    dataMSC%Rd = self%Rd
    
    read(nfile,*,iostat=iostat_data) self%Z1,self%Z2, self%Zasym, self%Zproj 
    if (iostat_data /= 0) then   ! Compatibility with older input files.
       if (mpiwrite==1) write(*,'(A/A)') 'WARNING: Not enough parameters of charge. Should be Z1,Z2,Zasym,Zproj', 'Assuming Z2 was omitted and adjusting accordingly.'
       self%Zproj = self%Zasym
       self%Zasym = self%Z2
       self%Z2 = self%Z1
    end if
    self%Zplus = self%Zproj * (self%Z1 + self%Z2)
    self%Zminus = self%Zproj * (self%Z1 - self%Z2)
    
    self%good_parity = .true.
    if (self%Rd > 0d0) then
       if (self%Z1 /= self%Z2) self%good_parity = .false.
       if (self%calculation_type==0 .and. self%origin/=0.5d0) self%good_parity = .false.
    end if
    
!!$   Set hlike: correct for H, He, H2+, and H2. To be modified for complex molecules to account for core electrons
    if (self%Z1+self%Z2 - self%Zasym == 1) then   ! One-electron target.
       self%hlike = .true.
    elseif (self%Z1+self%Z2 - self%Zasym == 2) then   ! Two-electron target.
       self%hlike = .false.
    else
       stop 'Error in input data: check Z and Zasym'
    end if

    if(en_eV < 10.0d0) then
      form='F3.1'
    elseif(en_eV < 100.0d0) then
      form='F4.1'
    elseif(en_eV < 1000.0d0) then
      form='F5.1'
    else
      form='F6.1'
    endif

    if(self%Zproj == -1) then
      projectile = 'Electron'
    elseif(self%Zproj == 1) then
      projectile = 'Positron'
    else
      error stop '***ERROR: Invalid projectile charge'
    endif
      
    form='(" [+] Scattering system: ",A," (",'//form//'," eV) --> ",A)'
    if(mpiwrite==1) then
      write(*,form) projectile, en_eV, self%target
      write(*,'("     > Internuclear separation R = ",F9.5)') self%Rd
      write(*,'("     > Target charges: Z1=",F3.1,", Z2=",F3.1,", Zasym=",F3.1)') self%Z1, self%Z2, self%Zasym
      n_elec = self%Z1 + self%Z2 - self%Zasym
      write(*,'("       (Number of electrons: ",I0,")")') n_elec
      if(n_elec /= 1 .and. n_elec /= 2) error stop '***ERROR: Zasym invalid, or this target not coded yet'
      write(*,*) '    > Good parity: ', self%good_parity
      write(*,*)

      if(self%calculation_type == 0) then
         write(*,'(" [+] Calculation type : 0 - Spherical with nonorthogonal (2l+1 type) radial basis functions")')
      elseif(self%calculation_type == 1) then
         write(*,'(" [+] Calculation type : 1 - Spherical with orthogonal (2l+2 type) radial basis functions")')
         error stop '***ERROR: This option is not coded'
      elseif (self%calculation_type == 2) then
         write(*,'(" [+] Calculation type : 2 - Spheroidal with orthogonal (m type) Hylleraas radial basis functions")') 
      elseif (self%calculation_type == 3) then
         write(*,'(" [+] Calculation type : 3 - Spheroidal with non-orthogonal (m type) Hylleraas radial basis functions")')
      else
        error stop '***ERROR: Invalid calculation type in data.in'
      endif

      if(self%calculation_type == 0 .or. self%calculation_type == 1) then !spherical - origin is specified
        write(*,'("     > Spherical coordinate-system origin: ",F5.3)') self%origin
        if (self%origin<0d0 .or. 1d0<self%origin) then
          error stop '***ERROR: Origin value must be a scale factor between 0 (on Z1) and 1 (on Z2)'
        endif
      endif
      write(*,*)

      if (self%calculation_mode == 0) then
         write(*,*) '[+] Calculation mode : 0 - structure only'
      elseif (self%calculation_mode == 1) then
         write(*,*) '[+] Calculation mode : 1 - scattering, no DCS'
      elseif (self%calculation_mode == 2) then
         write(*,*) '[+] Calculation mode : 2 - scattering with DCS and TV files'
      elseif (self%calculation_mode == 3) then
         write(*,*) '[+] Calculation mode : 3 - scattering with TV files (no DCS)'
      else
         print*, '*** ERROR: Wrong calculation mode in data.in. It must be either:'
         print*, '        0: structure only'
         print*, '        1: structure + scattering, no DCS'
         print*, '        2: structure + scattering with DCS and TV files created'
         print*, '        3: structure + scattering with TV files created (skips DCS)'
         error stop
      end if
  
      if(self%calculation_mode > 0 .and. self%Zproj == 1) then !positron scattering

        if(l_ps_not_specified) write(*,*) '*** WARNING: l_Ps not specified in data.in, continuing without Ps-channels'

        if(self%calculation_mode > 0) then
          if(self%l_Ps < 0) then
            write(*,'("[+] Single-centre positron-scattering calculations, l_Ps =",4i0)') self%l_Ps
          else
            write(*,'("[+] Two-centre positron-scattering calculations, l_Ps =",4i0)') self%l_Ps
          endif
        endif
      endif !positron scattering
      write(*,*)
    endif !mpiwrite


    read(nfile,*) self%la_core
    !if(mpiwrite .eq. 1) write(*,'("la_core: ",I5)') self%la_core   !TODO print info when these are utilised
    allocate(self%npr_core(0:self%la_core))
    
    read(nfile,*) (self%npr_core(l), l=0,self%la_core)
    !if(mpiwrite .eq. 1) write(*,'("npr_core: ",20I5)') (self%npr_core(l), l=0,self%la_core)

    ! Read in the bases. This is handled differently in spherical to spheroidal.
    if ( self%calculation_type == 0 ) then   ! SPHERICAL.
!!$ Laguerre basis configurations     
      read(nfile,*) self%labot, self%latop
      allocate(self%alpha(0:self%latop))
      allocate(self%nps(0:self%latop))
      read(nfile,*) (self%nps(l),self%alpha(l), l=self%labot,self%latop)
       
      read(nfile,*) dataMSC%MSC_nconfig, dataMSC%labot, dataMSC%latop
      if(dataMSC%MSC_nconfig < 0) error stop "*** ERROR: MSC_nconfig < 0... Need to check what this really does in the code."
      allocate(dataMSC%alpha(0:dataMSC%latop))
      allocate(dataMSC%nps(0:dataMSC%latop))    
      read(nfile,*) (dataMSC%nps(l),dataMSC%alpha(l), l=dataMSC%labot,dataMSC%latop)
      
      if(mpiwrite==1)  then
        if(self%hlike .and. dataMSC%MSC_nconfig == 0) then
          write(*,'(" [+] Spherical Lageurre basis: (used for final one-electron structure calculations)")')
        else
          write(*,'(" [+] Spherical Lageurre basis: (used to generate accurate MOs for second diagonalisation)")')
        endif
        write(*,'("     > labot,latop: ",2I5)') self%labot, self%latop
        write(tempchar,'(I0)') self%latop-self%labot+1
        write(*,'("     > nps(l),alpha(l): ",'//tempchar//'("(",I0,",",F7.4,"), "),"(",I0,",",F7.4,")")') (self%nps(l),self%alpha(l), l=self%labot,self%latop)
        write(*,*)
      endif
!!$ Molecular State Configuration 
       
       if(mpiwrite==1) then
         if(self%hlike .and. dataMSC%MSC_nconfig > 0) then
           write(*,*)   "[+] Second spherical basis: (Hybrid MO + Laguerre, used for final one-electron structure calculations)"
         elseif(dataMSC%MSC_nconfig > 0) then
           write(*,*)   "[+] Second spherical basis: (Hybrid MO + Laguerre, used for final one- and two-electron structure calculations)"
         endif
         if(dataMSC%MSC_nconfig == 0 .and. .not.self%hlike) then
           write(*,*) "    > MSC_nconfig = 0: basis will be formed only from the MOs saved from the first diagonalisation"
         elseif(dataMSC%MSC_nconfig > 0) then
           write(*,'("     > labot,latop: ",2I5)') dataMSC%labot, dataMSC%latop
           write(tempchar,'(I0)') dataMSC%latop-dataMSC%labot+1
           write(*,'("     > nps(l),alpha(l): ",'//tempchar//'("(",I0,",",F7.4,"), "),"(",I0,",",F7.4,")")') (dataMSC%nps(l),dataMSC%alpha(l), l=dataMSC%labot,dataMSC%latop)
           write(*,'("     > MSC_nconfig = ",I0,": the first ",I0," MOs saved from first diagonalisation will be subtituted into the second basis")') dataMSC%MSC_nconfig, dataMSC%MSC_nconfig
         endif !nconfig
       endif

       if (allocated(dataMSC%alpha_nl)) deallocate(dataMSC%alpha_nl) 
       !Liam added .true. here so alpha_nl array is defined when using natural orbitals
       if (.true. .or. ABS(dataMSC%MSC_nconfig) /= 0)  allocate(dataMSC%alpha_nl(MAXVAL(dataMSC%nps(dataMSC%labot:dataMSC%latop)),dataMSC%labot:dataMSC%latop))
       read(nfile,*) self%use_sub, dataMSC%labot_diff, dataMSC%latop_diff

       if (self%use_sub > 0) then
          if ( dataMSC%labot_diff < dataMSC%labot .OR. dataMSC%latop_diff > dataMSC%latop) then 
             if (mpiwrite == 1 )  print*,"***ERROR: sub-basis labot or latop are not within the confines of the MSC Laguerre basis"
             error stop 
          end if
          if ( allocated(dataMSC%nps_diff) ) deallocate(dataMSC%nps_diff, dataMSC%alpha_diff)
          allocate(dataMSC%nps_diff(dataMSC%labot_diff:dataMSC%latop_diff), dataMSC%alpha_diff(dataMSC%labot_diff:dataMSC%latop_diff))
          read(nfile,*) (dataMSC%nps_diff(l), dataMSC%alpha_diff(l), l=dataMSC%labot_diff,dataMSC%latop_diff)
          if(mpiwrite==1) then
            write(*,'("     > use_sub = ",I0,": functions from the following basis will replace Laguerres (not MOs) in the Hybrid basis")') self%use_sub
            write(tempchar,'(I0)') dataMSC%latop_diff-dataMSC%labot_diff+1
            write(*,'("     > nps(l),alpha(l): ",'//tempchar//'("(",I0,",",F7.4,"), "),"(",I0,",",F7.4,")")') (dataMSC%nps_diff(l),dataMSC%alpha_diff(l), l=dataMSC%labot_diff,dataMSC%latop_diff)
          endif !mpiwrite
       else
          dataMSC%labot_diff = 0
          dataMSC%latop_diff = -1
          read(nfile,*)   ! Skip the different alpha line.
       end if
       if(mpiwrite==1) write(*,*)
       
       if ( ABS(dataMSC%MSC_nconfig) /= 0) then
          do l = dataMSC%labot, dataMSC%latop
             dataMSC%alpha_nl(:,l) = dataMSC%alpha(l)
             if ( self%use_sub > 0 ) then ! overwrite alpha with different alpha's
                if ( l>=dataMSC%labot_diff .AND. l<=dataMSC%latop_diff ) dataMSC%alpha_nl(1:dataMSC%nps_diff(l),l) = dataMSC%alpha_diff(l)
             end if
          end do
       end if

       !TODO: change spherical code to use the basis_1e, basis_2e, and sub_basis_2e structures
       !For now copying into these for use with the convert_input_file subroutine

       basis_1e%labot = self%labot ; basis_1e%latop = self%latop
       allocate( basis_1e%nps(0:self%latop), basis_1e%alpha(0:self%latop) )
       basis_1e%nps = self%nps; basis_1e%alpha = self%alpha
       
       basis_2e%labot = dataMSC%labot ; basis_2e%latop = dataMSC%latop
       allocate( basis_2e%nps(0:dataMSC%latop), basis_2e%alpha(0:dataMSC%latop) )
       basis_2e%nps = dataMSC%nps; basis_2e%alpha = dataMSC%alpha
       
       sub_basis_2e%labot = dataMSC%labot_diff ; sub_basis_2e%latop = dataMSC%latop_diff
       allocate( sub_basis_2e%nps(0:dataMSC%latop_diff), sub_basis_2e%alpha(0:dataMSC%latop_diff) )
       sub_basis_2e%nps = dataMSC%nps_diff; sub_basis_2e%alpha = dataMSC%alpha_diff


    ! SPHEROIDAL.
    elseif ( self%calculation_type==2 .or. self%calculation_type==3 ) then
       ! Read in basis for one-electron diagonalisation.
       labot = -1; latop = -1
       read(nfile,*) labot, latop
       basis_1e%labot = labot ; basis_1e%latop = latop
       allocate( basis_1e%nps(0:latop), basis_1e%alpha(0:latop) )
       read(nfile,*) ( basis_1e%nps(l),basis_1e%alpha(l), l=labot,latop )
      
       if(mpiwrite==1)  then
         if(self%hlike) then
           write(*,'(" [+] Spheroidal Hylleraas basis:")')
         else !two-electron target
           write(*,'(" [+] First spheroidal Hylleraas basis: (used for one-electron structure")')
           write(*,'("                                        and to generate accurate MOs for second diagonalisation)")')
         endif
         write(*,'("     > labot,latop: ",2I5)') basis_1e%labot, basis_1e%latop
         write(tempchar,'(I0)') basis_1e%latop-basis_1e%labot+1
         write(*,'("     > nps(l),alpha(m): ",'//tempchar//'("(",I0,",",F7.4,"), "),"(",I0,",",F7.4,")")') (basis_1e%nps(l),basis_1e%alpha(l), l=basis_1e%labot,basis_1e%latop)
         write(*,*)
       endif

       ! Read in basis for two-electron diagonalisation.
       read(nfile,*) self%use_MSC, labot,latop
       basis_2e%labot = labot ; basis_2e%latop = latop
       basis_2e%mabot = labot ; basis_2e%matop = latop
       allocate( basis_2e%nps(0:latop), basis_2e%alpha(0:latop) )
       read(nfile,*) ( basis_2e%nps(l),basis_2e%alpha(l), l=labot,latop )
       if(mpiwrite==1 .and. .not.self%hlike) then !second basis only used for two-electron structure
         if(self%use_MSC == 0) then
           write(*,*)   "[+] Second spheroidal basis: (Pure Hylleraas, used for two-electron structure calculations)"
         else
           write(*,*)   "[+] Second spheroidal basis: (Hybrid MO + Hylleraas, used for two-electron structure calculations)"
         endif
         write(*,'("     > labot,latop: ",2I5)') basis_2e%labot, basis_2e%latop
         write(tempchar,'(I0)') basis_2e%latop-basis_2e%labot+1
         write(*,'("     > nps(l),alpha(m): ",'//tempchar//'("(",I0,",",F7.4,"), "),"(",I0,",",F7.4,")")') (basis_2e%nps(l),basis_2e%alpha(l), l=basis_2e%labot,basis_2e%latop)
         if(self%use_MSC > 0) then
           write(*,'("     > MSC_nconfig = ",I0,": the first ",I0," MOs saved from first diagonalisation will be subtituted into the second basis")') self%use_MSC, self%use_MSC
         endif
       end if

       ! Read in sub-basis for two-electron diagonalisation
       read(nfile,*) self%use_sub, labot,latop
       sub_basis_2e%labot = labot ; sub_basis_2e%latop = latop
       sub_basis_2e%mabot = labot ; sub_basis_2e%matop = latop
       allocate( sub_basis_2e%nps(0:latop), sub_basis_2e%alpha(0:latop) )
       read(nfile,*) ( sub_basis_2e%nps(l),sub_basis_2e%alpha(l), l=labot,latop )
       
       if(mpiwrite==1 .and. .not.self%hlike) then !sub basis only used for two-electron structure
         write(*,'("     > use_sub = ",I0,": functions from the following basis will replace Laguerres (not MOs) in the Hybrid basis")') self%use_sub
         write(tempchar,'(I0)') sub_basis_2e%latop-sub_basis_2e%labot+1
         write(*,'("     > nps(l),alpha(m): ",'//tempchar//'("(",I0,",",F7.4,"), "),"(",I0,",",F7.4,")")') (sub_basis_2e%nps(l),sub_basis_2e%alpha(l), l=sub_basis_2e%labot,sub_basis_2e%latop)
       endif 
    end if ! Different ways spherical/spheroidal deals with MSC basis.
       

    read(nfile,*) self%Mtot_start, self%Mtot_stop, self%ipar, self%nent, self%iborn
    if(.not.self%good_parity) self%ipar = 0
    if(mpiwrite==1 .and. self%calculation_mode > 0) then
      write(*,'("[+] Mtot_start, Mtot_stop, ipar, nent, iborn: ",6I5)') self%Mtot_start, self%Mtot_stop, self%ipar, self%nent, self%iborn
      write(*,*)
    endif
    self%iBorn_amp = 0
    if(self%iBorn .eq. 2) then
       self%iBorn = 1
       self%iBorn_amp = 1
    endif
    
    read(nfile,*) self%rmax, self%qmax, self%ndouble, self%npdbl, self%npwave, self%ltmax
    read(nfile,*,iostat=iostat_data) self%formcut, self%regcut, self%expcut, self%iweight
    if(iostat_data/=0) then
      backspace(nfile)
      read(nfile,*) self%formcut, self%regcut, self%expcut
      self%iweight = 0
    endif
    if(mpiwrite==1) then 
      write(*,'(" [+] Radial grid parameters:")')
      write(*,'("     > rmax,qmax,ndouble,npdbl,npwave,ltmax: ",2F15.5,4I5)') self%rmax, self%qmax, self%ndouble, self%npdbl, self%npwave, self%ltmax
      write(*,'("     > formcut,regcut,expcut: ",3E10.2)') self%formcut, self%regcut, self%expcut
      if(self%iweight == 0) then
        write(*,'("     > Integration weights: SIMPSON")')
      else
        write(*,'("     > Integration weights: BOOLE")')
      endif
      write(*,*)
    endif
    
    read(nfile,*) self%corep, self%r0, self%gamma, self%rho
    !if(mpiwrite .eq. 1) write(*,'("corep,r0,gamma,rho: ",4F10.3)') self%corep, self%r0, self%gamma, self%rho !TODO: print when it is utilised
    
    read(nfile,*) self%iSlater, self%iosc, self%irec
    if(mpiwrite==1) then
      write(*,'(" [+] iSlater: ",I5,", iosc=",I5,", irec=",I5 )') self%iSlater,self%iosc, self%irec
      write(*,*)
    endif
    
    
    read(nfile,*) self%Mt_min, self%Mt_max
    if ( self%calculation_type==2 .or. self%calculation_type==3 ) then
       basis_1e%mabot = self%Mt_min
       basis_1e%matop = self%Mt_max
    end if    

    allocate(self%nst(0:self%Mt_max,-1:1))
    self%nst(:,:) = 0
    read(nfile,*) (self%nst(m,1), m=self%Mt_min,self%Mt_max)
    read(nfile,*) (self%nst(m,-1), m=self%Mt_min,self%Mt_max)


    if(.not.self%good_parity) then
      self%nst(:,0) = self%nst(:,1) + self%nst(:,-1)
      self%nst(:,1) = 0
      self%nst(:,-1) = 0
    endif
    
    if(mpiwrite==1) then
      write(*,'(" [+] One-electron target states to be saved:")')
      write(*,'("     > Mt_min, Mt_max: ",I0,2X,I0 )') self%Mt_min, self%Mt_max
      if(self%good_parity) then
        write(*,'("     > nst(m,+1): ",20(I0,2X))') (self%nst(m,1), m=self%Mt_min,self%Mt_max)
        write(*,'("     > nst(m,-1): ",20(I0,2X))') (self%nst(m,-1), m=self%Mt_min,self%Mt_max)
      else
        write(*,'("     > nst(m): ",20(I0,2X))') (self%nst(m,0), m=self%Mt_min,self%Mt_max)
      endif
      write(*,*)
    endif
    
    if(myid.eq.0) then
       if ( dataMSC%MSC_nconfig /=  0 ) then
!             if( self%target(1:3) == 'H2+' ) call set_data_H2plus( data_target, max(maxval(self%nst(:,:)),maxval(dataMSC%nst(:,:))), max(self%latop,dataMSC%latop), self%Mt_max )
          if (self%target(1:3) == 'H2+') then
             call set_data_H2plus( data_target,maxval(self%nst(:,:)),max(self%latop,dataMSC%latop), self%Mt_max )
          elseif (self%target(1:2) == 'H2') then
!             call set_data_H2( data_target,maxval(self%nst(:,:)),max(self%latop,dataMSC%latop), self%Mt_max )
          end if
       else   
          if( self%target(1:3) == 'H2+' ) then
             call set_data_H2plus( data_target, maxval(self%nst(:,:)), self%latop, self%Mt_max )
          elseif (self%target(1:2) == 'H2') then
!             call set_data_H2( data_target,maxval(self%nst(:,:)), self%latop, self%Mt_max )
          end if
       end if
    end if
    
    read(nfile,*) self%idistpot, self%Ldw, self%ndw, self%aldw

    if((self%calculation_type == 2 .or. self%calculation_type == 3) .and. self%idistpot > 0) then
      if(mpiwrite==1) print*, '***ERROR: Distorting potential not coded in Spheroidal implementation'
      error stop
    endif

    read(nfile,*) self%ifirst, self%theta
    
    ! Mark: If the projectile is a positron setting ifirst = 0, so we just include Direct Matrix elements
    if ( self%Zproj > 0d0 .AND. self%ifirst == 1 ) self%ifirst = 0


!!$ DF 8-10-2018
    self%iAnalitBorb = 0
    if( self%ifirst .eq. -10) then
       self%ifirst = 0
       self%iAnalitBorb = 1
    endif
    
    !We don't use the ifirst parameter anymore - here we set the replacement parameters according to the supplied ifirst 
    self%pwborn = .false. !this cannot be set in the old data.in file
    if(self%ifirst == -1) then
      self%exchange = .false.
      self%UBA = .true.
    elseif(self%ifirst == 0) then
      self%exchange = .false.
      self%UBA = .false.
    elseif(self%ifirst == 1) then
      self%exchange = .true.
      self%UBA = .false.
    endif
    
    read(nfile,*) self%inc
    self%non_uniq = .FALSE.
    if ( self%ifirst == 1 .AND. self%theta /= 0d0 .AND. self%Zproj== -1d0 ) self%non_uniq = .TRUE.
    
    ! Mark: Input of new variables for either fixed or average orientation. If fixed incoming angles read thrad and phrad.
    read(nfile,*) self%Lpmax, self%orient, self%thrad, self%phrad
    
    if(mpiwrite==1 .and. self%calculation_mode > 0) then
      write(*,'(" [+] idistpot, LDW, nDW, alDW: ",3I5,F10.4)') self%idistpot, self%Ldw, self%ndw, self%aldw
      write(*,*)
      write(*,'(" [+] ifirst,theta: ",I5,F10.4)') self%ifirst, self%theta
      if(self%Zproj == -1) write(*,'("     >  Solving non-uniqueness: ", L)') self%non_uniq
      write(*,*)
      write(*,'(" [+] inc: ",I5)') self%inc
      write(*,*)
      write(*,'(" [+] Lpmax, orient, thrad, phrad: ",2I5,2F10.4)') self%Lpmax, self%orient, self%thrad, self%phrad
      write(*,*)
    endif

    mint = 4
    Lmaxkgt =  self%Lpmax
    self%Lmaxkg = Lmaxkgt
    allocate(nkt(mint),skt(mint))
    allocate(self%nk(mint,0:Lmaxkgt),self%sk(mint,0:Lmaxkgt),self%NBND(0:Lmaxkgt))
    
    icheck = 0   ! check that there at least one line for kgrid
    Lbot = 0
    LSW_old = 0
    do 
       LSW = -1; NBNDt = -1
       read(nfile,*,err=5,end=5) LSW, NBNDt , (nkt(i),skt(i), i=1,mint)
       !if(mpiwrite.eq.1) print*,"****", LSW, NBNDt , (nkt(i),skt(i), i=1,mint)
       LSW_old = LSW
       icheck = 1
       if(LSW .le. Lmaxkgt) then
          Ltt = LSW
       else
          Ltt = Lmaxkgt
       endif
       if(Lbot .gt. Ltt) then
          if(myid.eq.0)   print*,"***ERROR in input_data.f90: readin(): Lbot .gt. Ltt", Lbot, Ltt
          stop
       endif
       do i=1,mint
          self%NBND(Lbot:Ltt) = NBNDt
          self%nk(i,Lbot:Ltt) = nkt(i)
          self%sk(i,Lbot:Ltt) = skt(i)
       enddo
       Lbot = Ltt +1
       if(Ltt .eq. Lmaxkgt) then
          ! read until the error or end of file just to mak esure that when reading frm potl file all required line are read.
          do   
             read(nfile,*,err=5,end=5) LSW, NBNDt , (nkt(i),skt(i), i=1,mint)
          enddo
          
          exit
       endif
    enddo
    
5   continue
    
    if(icheck .eq. 0) then
       if(myid.eq.0) print*,'***ERRROR in input_data.f90: icheck=0, no kgrid line!'
       stop
    endif

    Lmaxkgt = min(Lmaxkgt,LSW_old)   !!! check if this is OK, temporary coding
    self%Lmaxkg = Lmaxkgt
    
    !!! TODO


    ! Mark: Modified Ldw = Max L of projectile when bound states from Coulomb or Distorted waves are included. 
    self%lbnd_max = -1
    do i = 0, self%Lmaxkg
       if ( self%NBND(i) .gt. 0 ) then
          self%lbnd_max = i
       end if
    end do

    if (self%idistpot <= 0) self%Ldw = -1
    !if(mpiwrite.eq.1) print*,"Ldw=",self%Ldw
    !if(mpiwrite.eq.1) print*,"Lmaxkg=",Lmaxkgt
!!$ Make sure that Ldw is not larger that Lmaxkgt    
    if(self%Ldw .gt. self%Lmaxkg) then
       self%Ldw = self%Lmaxkg 
       if(myid.eq.0) print*,'***WARNING: You have too large value of Ldw, changing it to maximum L of projectile:'
       if(myid.eq.0) print*,'            Ldw = ',self%Ldw
    endif
    
    if(mpiwrite==1 .and. self%calculation_mode > 0) then
      write(*,*) '[+] Kgrid and bound projectile-state info:'
      do Ltt = 0,Lmaxkgt
        write(*,'("     > LSW,NBND, nk(1:4),sk(1:4): ",2I5,4(I5,F7.2))') Ltt, self%NBND(Ltt), (self%nk(i,Ltt),self%sk(i,Ltt), i=1,mint)
      enddo
    endif


    if(nfile .eq. 10)   close(nfile)
    
    deallocate(nkt,skt)

!!! Rav: following is to read Ps-channel input data
    if(self%l_Ps >= 0) call readin_Ps(self,mpiwrite)
!R     
    
  if(mpiwrite==1) then
    write(*,'(A)') '################### END CALCULATION DETAILS ####################'
    write(*,*)
  endif
    
  inquire(file='print_channel_timing',exist=ex)
  if(ex) self%print_channel_timing = .true.
  if(self%print_channel_timing) print*, '>>>>> PRINTING CHANNEL TIMING'

  end subroutine readin_old
  
  subroutine readin_Ps(self,mpiwrite)
    use MPI_module
    implicit none
    type(input), intent(inout) :: self
    integer, intent(in) :: mpiwrite
    integer :: nfile_Ps, iostat_data, l
    logical :: ex
  
    nfile_Ps=33
    inquire(file='Ps_data.in',exist=ex)
    if(ex) then
       print*,'----------------- input data for Ps-formation channels____________'
     else
       if(myid.eq.0) print*,' Ps_data.in does not exist. For 1-centre CCC set l_Ps=-1 in data.in. STOP'
       STOP
    end if
    open(nfile_Ps,file='Ps_data.in',iostat=iostat_data)
    if(iostat_data.ne.0) then
       if(myid.eq.0) print*, '********  File  canot open file Ps_data.in'
       stop
    end if
    read(nfile_Ps,*) self%lpbot, self%lptop
    allocate(self%npbot(0:self%lptop),self%nptop(0:self%lptop))
    read(nfile_Ps,*) (self%npbot(l),self%nptop(l),l=self%lpbot,self%lptop)
    if(mpiwrite.eq.1) write(*,'("lpbot,lptop,npbot(lp)nptop(lp):",10I10)') self%lpbot, self%lptop,&
      &                  (self%npbot(l),self%nptop(l),l=self%lpbot,self%lptop)
    allocate(self%npsp(0:self%lptop),self%alphap(0:self%lptop))
    read(nfile_Ps,*) (self%npsp(l),self%alphap(l), l = self%lpbot, self%lptop)
    if(mpiwrite.eq.1) write(*,'("npsp(l),alphap(l):",10(I5,F10.4))') (self%npsp(l),self%alphap(l), l = self%lpbot, self%lptop)
    read(nfile_Ps,*)  self%igz, self%igp, self%analyticd, self%numericalv,self%lstoppos
    if(mpiwrite.eq.1) write(*,'(" igz, igp, analyticd, numericalv,lstoppos:",(2i5,2L5,i5))') &
      &     self%igz, self%igp, self%analyticd, self%numericalv,self%lstoppos
    print*,'----------------- end of input data for Ps-formation channels____________'
 
    close(nfile_Ps)
  end subroutine readin_Ps

!!$-------------------------------------------------------------------------------------------

! Used in potl head. MPI directives when calling potl routines.  
  subroutine write_inputdata(self, basis_1e,basis_2e,sub_basis_2e, nfile)
    implicit none
    type(input), intent(in):: self
    type(basis_input), intent(in) :: basis_1e, basis_2e, sub_basis_2e
    integer, intent(in):: nfile

    logical:: ex, op
    integer:: i, l, m, mint, Lmaxkgt, Ltt, LL, lprint, Ltt_max, calc_type
    real(dpf):: eV

    eV = self%eV  !  27.2116
    calc_type = self%calculation_type

    inquire(nfile,exist=ex,opened=op)
    if(.not. ex) then
       print*,'File with unit number', nfile,'  does not exists. STOP'
       STOP
       elseif(.not. op) then
       print*,'File with unit number', nfile,'  is not opened. STOP'
       STOP
    end if

    write(nfile,'(A10)') self%target

    write(nfile,'(F15.5)') self%energy * self%eV

    write(nfile,'(3I5,2X,"!")') calc_type, self%calculation_mode, self%l_Ps

    write(nfile,*) self%Rd, self%origin

    write(nfile,'(4F5.1)') self%Z1,self%Z2, self%Zasym, self%Zproj

    write(nfile,'(I5)') self%la_core
    write(nfile,'(20I5)') (self%npr_core(l), l=0,self%la_core)


    if (calc_type == 0) then
       write(nfile,'(2I5)') self%labot, self%latop
       write(nfile,'(20(I5,F10.4))') (self%nps(l),self%alpha(l), l=self%labot,self%latop)

       write(nfile,'(3I5)') dataMSC%MSC_nconfig, dataMSC%labot, dataMSC%latop
       write(nfile,'(20(I5,F10.4))') (dataMSC%nps(l),dataMSC%alpha(l), l=dataMSC%labot,dataMSC%latop)

       write(nfile,'(3I5)') self%use_sub, dataMSC%labot_diff, dataMSC%latop_diff
       write(nfile,'(20(I5,F10.4))') (dataMSC%nps_diff(l),dataMSC%alpha_diff(l),l=dataMSC%labot_diff,dataMSC%latop_diff)

    elseif ( calc_type==2 .or. calc_type==3 ) then
       write(nfile,'(2I5)') basis_1e%labot, basis_1e%latop
       write(nfile,'(20(I5,F10.4))') (basis_1e%nps(l),basis_1e%alpha(l), l=basis_1e%labot,basis_1e%latop)

       write(nfile,'(3I5)') self%use_MSC, basis_2e%labot, basis_2e%latop
       write(nfile,'(20(I5,F10.4))') (basis_2e%nps(l),basis_2e%alpha(l), l=basis_2e%labot,basis_2e%latop)

       write(nfile,'(3I5)') self%use_sub, sub_basis_2e%labot, sub_basis_2e%latop
       write(nfile,'(20(I5,F10.4))') (sub_basis_2e%nps(l),sub_basis_2e%alpha(l), l=sub_basis_2e%labot,sub_basis_2e%latop)

    end if


    write(nfile,'(6I5)') self%Mtot_start, self%Mtot_stop, self%ipar, self%nent, self%iborn

    write(nfile,'(2F15.5,4I5)') self%rmax, self%qmax, self%ndouble, self%npdbl, self%npwave, self%ltmax

    write(nfile,'(3E10.2)') self%formcut, self%regcut, self%expcut

    write(nfile,'(4F10.3)') self%corep, self%r0, self%gamma, self%rho

    write(nfile,'(I5)') self%iSlater, self%iosc, self%irec

    write(nfile,'(20I5)') self%mt_min, self%mt_max
    write(nfile,'(20I5)') (self%nst(m,1), m=self%mt_min, self%mt_max)
    write(nfile,'(20I5)') self%mt_min, self%mt_max
    write(nfile,'(20I5)') (self%nst(m,-1), m=self%mt_min, self%mt_max)

    write(nfile,'(3I5,F15.5)') self%idistpot, self%Ldw, self%ndw, self%aldw

    write(nfile,'(I5,F10.4)') self%ifirst, self%theta

    write(nfile,'(I5)') self%inc

    write(nfile,'(I5,I5,2F10.5)') self%Lpmax, self%orient, self%thrad, self%phrad
     
    mint = 4
    Lmaxkgt = self%Lmaxkg 
    Ltt_max = max(0,Lmaxkgt-1)
    LL = 0
    do Ltt = 0,Ltt_max
       lprint = 0

       if(Ltt+1 .le. Lmaxkgt) then
          if( self%NBND(Ltt+1) .ne. self%NBND(Ltt)) lprint=1
          
          do i=1,4
             if(self%nk(i,Ltt+1) .ne. self%nk(i,Ltt) .or. self%sk(i,Ltt+1) .ne. self%sk(i,Ltt)) lprint=1
          enddo
          
          LL = Ltt
          
          if(Ltt .eq. Lmaxkgt-1 .and. lprint .eq. 0) then
             lprint = 1
             LL = Lmaxkgt
          endif
       else
          lprint = 1
       endif

       if(lprint .eq. 1) then
          write(nfile,'(2I5,4(I5,F7.2))') LL, self%NBND(Ltt), (self%nk(i,LL),self%sk(i,LL), i=1,mint)
       endif
    enddo

  end subroutine write_inputdata

 
  subroutine convert_input_file( self, basis_1e,basis_2e,sub_basis_2e )
    implicit none
    type(input), intent(inout) :: self
    type(basis_input),intent(in) :: basis_1e,basis_2e,sub_basis_2e

    integer :: nunit
    character(:), allocatable :: filename, tempchar, tempchar2
    integer :: tempi1, tempi2, l, i
    logical :: ex, templ, star_prev, use_second_basis

    filename='data.in_new'

    inquire(file=filename,exist=ex)
    if(ex) error stop '*** ERROR converting old data.in to new format: tried to make file ''data.in_new'' but it already exists'

    open(newunit=nunit,file=filename,action='write',status='new')

    write(nunit,'(A)') '#MCCC-FN'
    write(nunit,'(A)')
    write(nunit,'(A)') 'SCATTERING_SYSTEM'
    
      if(self%Z1==1 .and. self%Z2==1) then
        if(self%Zasym == 1 ) then !H2+
          tempchar = 'H2+'
        elseif(self%Zasym == 0) then !H2
          tempchar = 'H2'
        else
          error stop '*** ERROR converting old data.in to new format: couldn''t determine target'
        endif
      elseif( (self%Z1+self%Z2==3) .and. (self%Z1*self%Z2==2) ) then
        if(self%Zasym == 1 ) then !HeH+
          tempchar = 'HeH+'
        else
          error stop '*** ERROR converting old data.in to new format: couldn''t determine target'
        endif
      else
        error stop '*** ERROR converting old data.in to new format: couldn''t determine target'
      endif
      
      write(nunit,'("  TARGET                    : ",A)') tempchar
      write(nunit,'("  INTERNUCLEAR_SEPARATION   : ",F8.5)') self%Rd
      write(nunit,'("  ORIGIN                    : ",F7.5)') self%origin
      if(self%Zproj == -1) then
        tempchar = 'ELECTRON'
      elseif(self%Zproj == 1) then
        tempchar = 'POSITRON'
      else
        error stop '*** ERROR converting old data.in to new format: couldn''t determine projectile'
      endif
      write(nunit,'("  PROJECTILE                : ",A)') tempchar
      write(nunit,'("  ENERGY                    : ",F7.3," eV")') self%energy*self%eV
      write(nunit,'("  M_TOT                     : "I0,", ",I0)') self%Mtot_start, self%Mtot_stop
      if(.not.self%good_parity) then
        tempchar = '0'
      else
        if(self%ipar == 0) then
          tempchar = '-1, 1'
        elseif(self%ipar == -1) then
          tempchar = '-1'
        elseif(self%ipar == 1) then
          tempchar = '1'
        endif
      endif
      write(nunit,'("  PARITY                    : ",A)') tempchar
      write(nunit,'("  NUMBER_ENTRANCE_STATES    : ",I0)') self%nent
      write(nunit,'("  PROJECTILE_LMAX           : ",I0)') self%Lpmax
    write(nunit,'(A)') 'END'
    write(nunit,'(A)')
    write(nunit,'(A)') 'CALCULATION_MODES'
      if(self%calculation_type == 0) then
        tempchar = 'SPHERICAL'
      elseif(self%calculation_type == 2) then
        tempchar = 'SPHEROIDAL'
      else
        error stop '*** ERROR converting old data.in to new format: new input only accepts calculation_mode = 0 or 2' 
      endif
      write(nunit,'("  COORDINATE_SYSTEM         : ",A)') tempchar
      if(self%calculation_mode > 0) then
        tempchar = 'ON'
      else
        tempchar = 'OFF'
      endif
      write(nunit,'("  SCATTERING                : ",A)') tempchar
    write(nunit,'(A)') 'END'
    write(nunit,'(A)')
    write(nunit,'(A)') 'OUTPUT'
      if(self%iosc > 0) then
        tempchar = 'ON'
      else
        tempchar = 'OFF'
      endif
      write(nunit,'("  OSCILLATOR_STRENGTHS      : ",A)') tempchar
      write(nunit,'(A)') '  DIPOLE_MOMENTS            : OFF'
      if(self%calculation_mode == 1) then
        write(nunit,'(A)') '  DCS                       : OFF'
        write(nunit,'(A)') '  TV_FILES                  : OFF'
      elseif(self%calculation_mode == 2) then
        write(nunit,'(A)') '  DCS                       : ON'
        write(nunit,'(A)') '  TV_FILES                  : ON'
      elseif(self%calculation_mode == 3) then
        write(nunit,'(A)') '  DCS                       : OFF'
        write(nunit,'(A)') '  TV_FILES                  : ON'
      endif
      write(nunit,'(A)') '  ONE_EL_BASIS_LIST         : OFF'
      write(nunit,'(A)') '  TWO_EL_CONFIG_LIST        : OFF'
      write(nunit,'(A)') '  CI_COEFFICIENTS           : OFF' 
    write(nunit,'(A)') 'END'
    write(nunit,'(A)')
    write(nunit,'(A)') 'FIRST_DIAGONALISATION_BASIS'
      if(basis_1e%labot /= 0) error stop '*** ERROR converting old data.in to new format: new input doesn''t support one-el labot > 0' 
      write(nunit,'("  LMAX                      : ",I0)') basis_1e%latop
      tempi1 = basis_1e%nps(0)
      templ=.false.
      if(basis_1e%latop > 0) then
        templ = .true.
        do l=1, basis_1e%latop
          if(basis_1e%nps(l) /= basis_1e%nps(l-1)-1) templ = .false.
        enddo
      endif
      tempchar='   '
      write(tempchar,'(I0)') basis_1e%latop
      if(templ) then
        write(nunit,'("  BASIS_SIZE                : ",I0,"-l")') tempi1
      else
        if(basis_1e%latop == 0) then
          write(nunit,'("  BASIS_SIZE                : ",I0)') basis_1e%nps
        else
          write(nunit,'("  BASIS_SIZE                : ",'//tempchar//'(I0,", "),I0)') basis_1e%nps
        endif
      endif
      if(basis_1e%latop == 0) then
        write(nunit,'("  EXPONENTIAL_FALLOFF       : ",F6.3)') basis_1e%alpha
      else
        write(nunit,'("  EXPONENTIAL_FALLOFF       : ",'//tempchar//'(F6.3,","),F6.3)') basis_1e%alpha
      endif
    write(nunit,'(A)') 'END'
    write(nunit,*)
    write(nunit,'(A)') 'ONE_ELECTRON_STATES'
      if(self%Mt_min /= 0) error stop '*** ERROR converting old data.in to new format: new input doesn''t support one-el Mt_min > 0' 
      write(nunit,'("  M_MAX                     : ",I0)') self%Mt_max
      write(tempchar,'(I0)') self%Mt_max+1
      if(.not.data_in%good_parity) then
        self%nst(:,1) = self%nst(:,0)
        self%nst(:,-1) = 0
      endif
      write(nunit,'("  NST_POS_PAR               : ",'//tempchar//'(I0,", "),I0)') self%nst(:,1)
      write(nunit,'("  NST_NEG_PAR               : ",'//tempchar//'(I0,", "),I0)') self%nst(:,-1)
    write(nunit,'(A)') 'END'
    write(nunit,'(A)')
   

    if(self%calculation_type == 0 .and. dataMSC%MSC_nconfig > 0 .or. self%calculation_type == 2 .and. self%use_MSC >= 0) then
      tempchar2 = ''
      use_second_basis = .true.
    else
      !tempchar2 = '#'
      tempchar2 = ''
      use_second_basis = .false.
    endif
    write(nunit,'(A)')tempchar2//'SECOND_DIAGONALISATION_BASIS'
      if(use_second_basis) then
        write(nunit,'(A)') tempchar2//'  USE_SECOND_BASIS          : ON'
      else
        write(nunit,'(A)') tempchar2//'  USE_SECOND_BASIS          : OFF'
      endif
      
      if(basis_2e%labot /= 0) error stop '*** ERROR converting old data.in to new format: new input doesn''t support two-el labot > 0' 
      write(nunit,'("'//tempchar2//'  LMAX                      : ",I0)') basis_2e%latop
      tempi1 = basis_2e%nps(0)
      templ=.false.
      if(basis_2e%latop > 0) then
        templ = .true.
        do l=1, basis_2e%latop
          if(basis_2e%nps(l) /= basis_2e%nps(l-1)-1) templ = .false.
        enddo
      endif
      tempchar='   '
      write(tempchar,'(I0)') basis_2e%latop
      if(templ) then
        write(nunit,'("'//tempchar2//'  BASIS_SIZE                : ",I0,"-l")') tempi1
      else
        if(basis_2e%latop==0) then
          write(nunit,'("'//tempchar2//'  BASIS_SIZE                : ",I0)') basis_2e%nps
        else
          write(nunit,'("'//tempchar2//'  BASIS_SIZE                : ",'//tempchar//'(I0,", "),I0)') basis_2e%nps
        endif
      endif
      if(basis_2e%latop == 0) then
        write(nunit,'("'//tempchar2//'  EXPONENTIAL_FALLOFF       : ",F6.3)') basis_2e%alpha
      else
        write(nunit,'("'//tempchar2//'  EXPONENTIAL_FALLOFF       : ",'//tempchar//'(F6.3,","),F6.3)') basis_2e%alpha
      endif
      write(nunit,'(A)') tempchar2
      if(basis_1e%latop == 0 .and. basis_1e%nps(0) == 1 .and. basis_1e%alpha(0) == basis_2e%alpha(0)) then
        tempi1 = 0
      elseif(self%calculation_type == 0) then
        tempi1 = dataMSC%MSC_nconfig
      elseif(self%calculation_type == 2) then
        tempi1 = self%use_MSC
      endif
      if(.not.use_second_basis) tempi1 = 0
      write(nunit,'("'//tempchar2//'  INSERT_ONE_EL_STATES      : ",I0)') tempi1
      write(nunit,'(A)') tempchar2
      write(nunit,'(A)') tempchar2//'  INSERT_NATURAL_ORBITALS   : 0'
      write(nunit,'(A)') tempchar2//'  ONLY_NATURAL_ORBITALS     : OFF'
      write(nunit,'(A)') tempchar2//'  NATORB_GS_M               : 0'
      write(nunit,'(A)') tempchar2//'  NATORB_GS_PAR             : 1'
      write(nunit,'(A)') tempchar2//'  NATORB_GS_SPIN            : 0'
      write(nunit,'(A)') tempchar2
      if(self%use_sub > 0) then
        tempchar = 'ON'
      else
        tempchar = 'OFF'
      endif
      write(nunit,'("'//tempchar2//'  USE_SUB_BASIS             : ",A)') tempchar
      if(sub_basis_2e%labot /= 0) error stop '*** ERROR converting old data.in to new format: new input doesn''t support sub-basis labot > 0' 
      write(nunit,'("'//tempchar2//'  SUB_LMAX                  : ",I0)') sub_basis_2e%latop
      tempi1 = sub_basis_2e%nps(0)
      templ=.false.
      if(sub_basis_2e%latop > 0) then
        templ = .true.
        do l=1, sub_basis_2e%latop
          if(sub_basis_2e%nps(l) /= sub_basis_2e%nps(l-1)-1) templ = .false.
        enddo
      endif
      tempchar='   '
      write(tempchar,'(I0)') sub_basis_2e%latop+1
      print*, 'latop:', sub_basis_2e%latop
      print*, 'tempchar:', tempchar
      if(templ) then
        write(nunit,'("'//tempchar2//'  SUB_BASIS_SIZE            : ",I0,"-l")') tempi1
      else
        if(sub_basis_2e%latop == 0) then
          write(nunit,'("'//tempchar2//'  SUB_BASIS_SIZE            : ",I0)') sub_basis_2e%nps
        else
          write(nunit,'("'//tempchar2//'  SUB_BASIS_SIZE            : ",'//tempchar//'(I0,", "),I0)') sub_basis_2e%nps
        endif
      endif
      if(sub_basis_2e%latop == 0) then
        write(nunit,'("'//tempchar2  //'  SUB_BASIS_EXP_FALLOFF     : ",F6.3)') sub_basis_2e%alpha
      else
        write(nunit,'("'//tempchar2  //'  SUB_BASIS_EXP_FALLOFF     : ",'//tempchar//'(F6.3,","),F6.3)') sub_basis_2e%alpha
      endif

    write(nunit,'(A)') tempchar2//'END'
    write(nunit,*)

    if(.not.self%hlike) then
      inquire(file='F5',exist=ex)
      if(.not.ex) error stop '*** ERROR converting old data.in to new format: missing F5 file'
      open(unit=3,file='F5',action='read')
      read(3,*) !target
      read(3,*) self%Mmax2el
      allocate(self%nst_2e(0:self%Mmax2el,-1:1,0:1))
      read(3,*) (self%nst_2e(l,1,0), l=0,self%Mmax2el) ! positive parity, singlet
      read(3,*) (self%nst_2e(l,-1,0), l=0,self%Mmax2el) ! negative parity, singlet
      read(3,*) (self%nst_2e(l,1,1), l=0,self%Mmax2el) ! positive parity, triplet
      read(3,*) (self%nst_2e(l,-1,1), l=0,self%Mmax2el) ! negative parity, triplet
      read(3,*) !inc
      read(3,*) self%l_ion_core
      allocate(self%n_ion_core(0:self%l_ion_core))
      backspace(3)
      read(3,*) self%l_ion_core, (self%n_ion_core(l), l=0, self%l_ion_core)
      read(3,*) self%l12max, self%M12max
      if(self%l12max < basis_2e%latop) error stop '*** ERROR converting old data.in to new format: cannot handle l12max < basis_2e%latop'
      self%l12max = min(self%l12max,basis_2e%latop)
      allocate(self%nkin(0:self%l12max), self%nkout(0:self%l12max))
      read(3,*) (self%nkout(l), l=0,self%l12max)
      read(3,*) (self%nkin(l), l=0,self%l12max)    
      close(3)
      write(nunit,'(A)') 'TWO_ELECTRON_CONFIGURATIONS'
        write(nunit,'("  L_ION_CORE                : ",I0)') self%l_ion_core
        if(self%l_ion_core == 0) then
          write(tempchar,'(I0)') self%l_ion_core
          write(nunit,'("  N_ION_CORE                : ",I0)') self%n_ion_core
        else
          write(nunit,'("  N_ION_CORE                : ",'//tempchar//'(I0,", "),I0)') self%n_ion_core
        endif
        !write(nunit,'("  L12MAX                    : ",I0)') self%l12max
        write(nunit,'("  M12MAX                    : ",I0)') self%M12max
        write(tempchar,'(I0)') self%l12max
        write(nunit,'("  NKIN                      : ",'//tempchar//'(I2.1,", "),I2.1)') self%nkin
        !write(nunit,'("  NKOUT                     : ",'//tempchar//'(I2.1,", "),I2.1)') self%nkout
        write(nunit,'(A)') 'END'
        write(nunit,*)
        write(nunit,'(A)') 'TWO_ELECTRON_STATES'
        write(nunit,'("  M_MAX                     : ",I0)') self%Mmax2el
        if(self%Mmax2el == 0) then
          tempchar='I0'
        else
          write(tempchar,'(I0)') self%Mmax2el
          tempchar = tempchar//'(I0,", "),I0'
        endif
        write(nunit,'("  NST_POS_PAR_SINGLET       : ",'//tempchar//')') self%nst_2e(:,1,0)
        write(nunit,'("  NST_NEG_PAR_SINGLET       : ",'//tempchar//')') self%nst_2e(:,-1,0)
        write(nunit,'("  NST_POS_PAR_TRIPLET       : ",'//tempchar//')') self%nst_2e(:,1,1)
        write(nunit,'("  NST_NEG_PAR_TRIPLET       : ",'//tempchar//')') self%nst_2e(:,-1,1)
      write(nunit,'(A)') 'END' 
      write(nunit,*)
    endif

    write(nunit,'(A)') 'RADIAL_GRID'
      write(nunit,'("  RMAX                      : ",F5.1)') self%rmax
      write(nunit,'("  QMAX                      : ",F5.1)') self%qmax
      write(nunit,'("  NDOUBLE                   : ",I0)') self%ndouble
      write(nunit,'("  NPDBL                     : ",I0)') self%npdbl
      write(nunit,'("  NPWAVE                    : ",I0)') self%npwave
      write(nunit,'("  FORMCUT                   : ",ES7.1)') self%formcut
      write(nunit,'("  REGCUT                    : ",ES7.1)') self%regcut
      write(nunit,'("  EXPCUT                    : ",ES7.1)') self%expcut
    write(nunit,'(A)') 'END'
    write(nunit,'(A)')

    write(nunit,'(A)') 'KGRID'
      star_prev = .false.
      l=0
      do while (l <= self%Lpmax)
        templ = (l == self%Lpmax)
        if(l < self%Lpmax) then
          if(self%NBND(l) /= self%NBND(l+1) .or. any(self%nk(:,l) /= self%nk(:,l+1)) .or. any(self%sk(:,l) /= self%sk(:,l+1))) then
            templ = .true.
          endif
        endif
        if(star_prev) templ = .true.
        if(templ) then
          write(nunit,'(2X,2(I2,", "), 3(I2,",",F4.2,", "), I3,",",F4.2)') l,self%NBND(l), (self%nk(i,l),self%sk(i,l), i=1,4)
          star_prev = .false.
        else
          do while (l < self%Lpmax .and. self%NBND(l) == self%NBND(l+1) .and. all(self%nk(:,l) == self%nk(:,l+1)) .and. all(self%sk(:,l) == self%sk(:,l+1)))
            l=l+1
          enddo
          !l=l-1
          write(nunit,'(2X," *, ",I2,", ", 3(I2,",",F4.2,", "), I3,",",F4.2)') self%NBND(l), (self%nk(i,l),self%sk(i,l), i=1,4)
          star_prev = .true.
        endif
        l=l+1
      enddo
    write(nunit,'(A)') 'END'
    write(nunit,*)
    write(nunit,'(A)') 'CALCULATION_OPTIONS'
      write(nunit,'("  LTMAX                     : ",I0)') self%ltmax
      if(self%iborn == 1) then
        tempchar = 'ON'
      else
        tempchar = 'OFF'
      endif
      write(nunit,'("  BORN_COMPLETION           : ",A)') tempchar                 
      if(self%iAnalitBorb == 1) then
        tempchar = 'ON'
      else
        tempchar = 'OFF'
      endif
!      write(nunit,'("  ONLY_BORN                 : ",A)') tempchar                
      write(nunit,'("  THETA                     : ",F3.1)') self%theta               
      if(self%ifirst ==1) then
        tempchar = 'ON'
      else
        tempchar = 'OFF'
      endif
      write(nunit,'("  EXCHANGE                  : ",A)') tempchar              
      if(self%ifirst ==-1) then
        tempchar = 'ON'
      else
        tempchar = 'OFF'
      endif
      write(nunit,'("  UBA                       : ",A)') tempchar             
      if(self%inc == 1) then
        tempchar = 'ON'
      else
        tempchar = 'OFF'
      endif
      write(nunit,'("  REARRANGE                 : ",A)') tempchar
      if(self%idistpot == 1) then
        tempchar = 'ON'
      elseif(self%idistpot > 1) then
        error stop '*** ERROR converting old data.in to new format: new input cannot handle idistpot > 1'
      else
        tempchar = 'OFF'
        self%ldw = 0
      endif
      write(nunit,'("  DISTORTING_POTENTIAL      : ",A)') tempchar
      write(nunit,'("  DISTPOT_LMAX              : ",I0)') self%ldw
      write(nunit,'("  PROJ_BOUNDSTATE_NMAX      : ",I0)') self%ndw
      write(nunit,'("  PROJ_BOUNDSTATE_ALPHA     : ",F3.1)') self%aldw
      if(self%l_ps < 0) then
        tempchar = 'OFF'
      else
        tempchar = 'ON'
      endif
      write(nunit,'("  POSITRONIUM_FORMATION     : ",A)') tempchar
      write(nunit,'("  NO_SECOND_SPIN_CHANNEL    : ",A)') 'OFF'
      write(nunit,'("  OPTICAL_POTENTIAL         : ",A)') 'OFF'
      write(nunit,'("  WEAK_COUPLING             : ",A)') 'OFF'
      write(nunit,'("  P_SPACE                   : ",A)') 'ALL'
    write(nunit,'(A)') 'END'

    close(nunit)


    write(*,*)
    write(*,*) ' Finished converting data.in to new format - created data.in_new'
    write(*,*)

  end subroutine convert_input_file






!-=-= SWITCH =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine switch(var,w,block)
    use input_parser
    implicit none

    !-=-= INPUT
    logical, intent(out) :: var
    character(len=*), intent(in) :: w,block

    !-=-= LOCAL VARIABLES
    character(len=255) :: w2
                  
    !if(nitems == 1) then
    !  var = .true.
    !else
    if(nitems == 3 ) then
      call readu(w2)
      select case(w2)
        case('ON')
          var = .true.
        case('OFF')
          var = .false.
        case default
          call invalid_argument(w,block,w2)
      end select
    else
      call nitem_error(w,block)
    endif
  end subroutine switch
  

!-=-= HEADING_ERROR =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine heading_error(w)
    implicit none
    character(len=*), intent(in) :: w
  
    write(*,'(A)') '*** ERROR reading input: trailing text after '//trim(adjustl(w))//' heading.'
    error stop 
  end subroutine heading_error

!-=-= READ_LABEL =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine read_label(label,block)
    use input_parser
    implicit none
    character(len=*), intent(out) :: label
    character(len=*), intent(in) :: block
    character(len=255) :: colon
    logical :: error
    integer :: i
     
    error = .false.

    call readu(label)
    if(trim(adjustl(label)) == 'END') return
    if(nitems == 1) error = .true.
    if(.not. error) then
      call readu(colon)
      if(trim(adjustl(colon)) /= ':') error = .true.
    endif
    if(error) then
      write(*,*) '***ERROR: expected ":" after <'//trim(adjustl(label))//'> keyword'
      error stop
    endif
    
    do i=1, size(labels_read)
      if(trim(adjustl(labels_read(i))) == trim(adjustl(label))) then
        write(*,*) '*** ERROR reading input: Keyword <'//trim(adjustl(label))//'> appears more than once in the <'//&
          &trim(adjustl(block))//'> block'
        error stop
      endif
      if(trim(adjustl(labels_read(i))) == '') exit
    enddo
    labels_read(i) = trim(adjustl(label))

  end subroutine read_label
  
!-=-= NITEM_ERROR =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
subroutine nitem_error(w,block)
    implicit none
    character(len=*), intent(in) :: w, block
  
    write(*,'(A)') '*** ERROR reading input: incorrect number of arguments to the <'//trim(adjustl(w))//'> keyword &
      &in the <'//trim(adjustl(block))//'> block'
    error stop 
  end subroutine nitem_error
  
  
!-=-= EOF_DURING_BLOCK =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine eof_during_block(w)
    implicit none
    character(len=*), intent(in) :: w
  
    write(*,'(A)') '*** ERROR reading input: end of file reached while reading '//trim(adjustl(w))//' options.' 
    error stop 
  end subroutine eof_during_block
  

!-=-= INVALID_label_IN_BLOCK =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine invalid_label_in_block(w1, w2)
    implicit none
    character(len=*), intent(in) :: w1, w2
  
    write(*,'(A)') '*** ERROR reading input: invalid keyword '//trim(adjustl(w1))//' in the '//trim(adjustl(w2))//' block.'
    error stop
  end subroutine invalid_label_in_block
  

!-=-= INVALID_ARGUMENT =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine invalid_argument_i(label, block, i)
    implicit none
    character(len=*), intent(in) :: label, block
    integer, intent(in) :: i
    character(len=100) :: arg

    write(arg,*) i
    call invalid_argument_c(label,block,arg)
  
  end subroutine invalid_argument_i
  
  subroutine invalid_argument_r(label, block, r)
    implicit none
    character(len=*), intent(in) :: label, block
    real(dpf), intent(in) :: r
    character(len=100) :: arg

    write(arg,*) r
    call invalid_argument_c(label,block,arg)

  end subroutine invalid_argument_r
  
  subroutine invalid_argument_c(label, block, w)
    implicit none
    character(len=*), intent(in) :: label, block
    character(len=*), intent(in) :: w
    character(len=100) :: arg

    write(*,'(A)') '*** ERROR reading input: invalid argument "'//trim(adjustl(w))//'" to the <'//trim(adjustl(label))//'> &
      &keyword in the <'//trim(adjustl(block))//'> block'
    error stop
  end subroutine invalid_argument_c

!-=-= READ_BLOCK =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
subroutine read_block(block,blockname,core_input)
    use input_parser
    implicit none
    character(len=*), intent(out) :: block, blockname
    logical, intent(in) :: core_input
    integer :: i

    call readu(block)
    if(core_input) then 
      blockname = trim(adjustl(block))//'(CORE)'
    else
      blockname = trim(adjustl(block))
    endif
        
    do i=1, size(blocks_read)
      if(trim(adjustl(blocks_read(i))) == trim(adjustl(blockname))) then
        print*, 'BLOCKS_READ(i): |'//trim(adjustl(blocks_read(i)))//'|'
        write(*,*) '*** ERROR reading input: Block <'//trim(adjustl(blockname))//'> appears more than once'
        error stop
      endif
      if(trim(adjustl(blocks_read(i))) == '') exit
    enddo
    blocks_read(i) = trim(adjustl(blockname))

    labels_read = '' !at start of new block reset the labels_read array

  end subroutine read_block

  subroutine read_basis_size(latop, nps, block)
    use input_parser
    implicit none
    integer, intent(in) :: latop
    integer, dimension(0:latop), intent(out) :: nps
    character(len=*), intent(in) :: block
    character(len=100) :: tempc
    character(len=:), allocatable :: tempc1
    integer :: len_tempc1, len_tempc1_minus1, l

    call reada(tempc)
    tempc1 = trim(adjustl(tempc))
    len_tempc1 = len(tempc1)
    len_tempc1_minus1 = max(len_tempc1-1,1)
    if(tempc1(len_tempc1_minus1:len_tempc1) == '-l') then !N_l = N - l rule
      if(nitems > 3) then
        write(*,*) '*** ERROR reading input: too many arguments for the (N_l = N - l) rule in the <'//trim(adjustl(block))//'> block'
        error stop
      endif
      read(tempc1(1:len_tempc1-2),*) nps(0)
      do l=1, latop
        nps(l) = nps(0) - l
      enddo
    elseif(nitems == 3) then !Not N_l = N - l rule, and only one N supplied
      call reread(-1)
      call readi(nps(0))
      nps(1:latop) = nps(0)
    elseif(nitems == 2 + latop+1) then !specify nps per l
      call reread(-1)
      do l=0, latop
        call readi(nps(l))
      enddo
    else
      write(*,*) '*** ERROR reading input: wrong number of arguments to the <BASIS_SIZE> keyword in the <'//trim(adjustl(block))//'> block'
      error stop
    endif

    if(any(nps < 0)) then 
      write(*,*) '*** ERROR reading input: negative basis size in the <'//trim(adjustl(block))//'> block'
      write(*,*) ' BASIS_SIZE: ', nps
      error stop
    endif

  end subroutine read_basis_size
  
  subroutine read_alpha(latop, alpha, block)
    use input_parser
    implicit none
    integer, intent(in) :: latop
    real(dpf), dimension(0:latop), intent(out) :: alpha
    character(len=*), intent(in) :: block
    character(len=100) :: tempc
    character(len=:), allocatable :: tempc1
    integer :: len_tempc1, l

    if(nitems == 3) then !only one alpha supplied
      call readf(alpha(0))
      alpha(1:latop) = alpha(0)
    elseif(nitems == 2 + latop+1) then !specify alpha per l
      do l=0, latop
        call readf(alpha(l))
      enddo
    else
      write(*,*) '*** ERROR reading input: wrong number of arguments to the <EXPONENTIAL_FALLOFF> keyword in the <'//trim(adjustl(block))//'> block'
      error stop
    endif
    
    if(any(alpha < 0.0d0)) then 
      write(*,*) '*** ERROR reading input: negative alpha in the <'//trim(adjustl(block))//'> block'
      write(*,*) ' ALPHA: ', alpha
      error stop
    endif
  end subroutine read_alpha
  
  subroutine read_ALL_or_N(N,label,block)
    use input_parser
    implicit none
    integer, intent(out) :: N
    character(len=*), intent(in) :: label, block
    character(len=255) :: w

    call readu(w)
    if(trim(adjustl(w)) == 'ALL') then
      N = -1
    else
      call reread(-1) 
      call readi(N)
      if(N < 0) call invalid_argument(label,block,N)
    endif

  end subroutine read_ALL_or_N

end module input_data



!            select case (label)
!              case('MAX_L_PER_KGRID')
!                allocate(L_k_input(0:nitems-2))
!                L_k_input(0) = 0
!                do i=1, nitems-2
!                  call readi(L_k_input(i))
!                  if(L_k_input(i) < L_k_input(i-1)) call invalid_argument(label,block,L_k_input(i))
!                enddo
!              case('NUM_POINTS_FIRST_INT')
!                if(nitems-2 /= size(L_k_input)) call nitem_error(label,block)
!                do i=1, nitems-2
!                  if(L_k_input(i-1) >= self%Lpmax .and. i>1) exit
!                  call readi(tempi)
!                  if(tempi <= 0) call invalid_argument(label,block,tempi)
!                  self%nk(1,L_k_input(i-1):min(L_k_input(i),self%Lpmax)) = tempi
!                enddo
!              case('NUM_POINTS_SECOND_INT')
!                if(nitems-2 /= size(L_k_input)) call nitem_error(label,block)
!                do i=1, nitems-2
!                  if(L_k_input(i-1) >= self%Lpmax .and. i>1) exit
!                  call readi(tempi)
!                  if(tempi <= 0) call invalid_argument(label,block,tempi)
!                  self%nk(2,L_k_input(i-1):min(L_k_input(i),self%Lpmax)) = tempi
!                enddo
!              case('NUM_POINTS_THIRD_INT')
!                if(nitems-2 /= size(L_k_input)) call nitem_error(label,block)
!                do i=1, nitems-2
!                  if(L_k_input(i-1) >= self%Lpmax .and. i>1) exit
!                  call readi(tempi)
!                  if(tempi <= 0) call invalid_argument(label,block,tempi)
!                  self%nk(3,L_k_input(i-1):min(L_k_input(i),self%Lpmax)) = tempi
!                enddo
!              case('NUM_POINTS_AROUND_SINGULARITY')
!                if(nitems-2 /= size(L_k_input)) call nitem_error(label,block)
!                do i=1, nitems-2
!                  if(L_k_input(i-1) >= self%Lpmax .and. i>1) exit
!                  call readi(tempi)
!                  if(tempi < 0) call invalid_argument(label,block,tempi)
!                  if(tempi == 0) tempi = -10 !User requested zero points -> code expects negative number to ignore this
!                  self%nk(4,L_k_input(i-1):min(L_k_input(i),self%Lpmax)) = tempi
!                enddo
!              case('SECOND_INT_START')
!                if(nitems-2 /= size(L_k_input)) call nitem_error(label,block)
!                do i=1, nitems-2
!                  if(L_k_input(i-1) >= self%Lpmax .and. i>1) exit
!                  call readf(tempr)
!                  if(tempr <= 0) call invalid_argument(label,block,tempr)
!                  self%sk(1,L_k_input(i-1):min(L_k_input(i),self%Lpmax)) = tempr
!                enddo
!              case('THIRD_INT_START')
!                if(nitems-2 /= size(L_k_input)) call nitem_error(label,block)
!                do i=1, nitems-2
!                  if(L_k_input(i-1) >= self%Lpmax .and. i>1) exit
!                  call readf(tempr)
!                  if(tempr <= self%sk(2,L_k_input(i-1))) call invalid_argument(label,block,tempr)
!                  self%sk(2,L_k_input(i-1):min(L_k_input(i),self%Lpmax)) = tempr
!                enddo
!              case('FALLOFF_POWER')
!                if(nitems-2 /= size(L_k_input)) call nitem_error(label,block)
!                do i=1, nitems-2
!                  if(L_k_input(i-1) >= self%Lpmax .and. i>1) exit
!                  call readf(tempr)
!                  if(tempr <= 0.0d0) call invalid_argument(label,block,tempr)
!                  self%sk(3,L_k_input(i-1):min(L_k_input(i),self%Lpmax)) = tempr
!                enddo
!              case('WIDTH_AROUND_SINGULARITY')
!                if(nitems-2 /= size(L_k_input)) call nitem_error(label,block)
!                do i=1, nitems-2
!                  if(L_k_input(i-1) >= self%Lpmax .and. i>1) exit
!                  call readf(tempr)
!                  if(tempr <= 0.0d0) call invalid_argument(label,block,tempr)
!                  self%sk(4,L_k_input(i-1):min(L_k_input(i),self%Lpmax)) = tempr
!                enddo
!              case('END')
!                exit
!              case default
!                call invalid_label_in_block(label,block)
!            end select 



