module target_states

  use state_class

  implicit none

!$$********   This is declaration of an  object of the type basis_state. 


  type(basis_state), target, public :: TargetStates, TargetStates1el, TargetStates2el
  type(basis_state), pointer, public:: Target_Basis  ! TargetStates or TargetStates2el depending upon 1 or 2 electron target 

  type(basis_state) :: TargetStates_unrearranged !Liam added: save original description of one-electron states in terms of  
                                                 !unrearranged Laguerre basis (need for natural orbitals)

  contains

  subroutine multiply_by_scalar(self, s)
    !Liam: this subroutine multiplies the CI coefficients of a state by a scalar s
    implicit none

    type(state), intent(inout) :: self
    real*8, intent(in) :: s

    self%CI = self%CI * s
  end subroutine multiply_by_scalar

  subroutine sum_into_1el(state_l, state_r)
    use one_electron_func
    use ovlpste1me
    implicit none

    type(state), intent(inout) :: state_l
    type(state), intent(in) :: state_r
    type(sturmian_nr) :: temp_sturm
    integer :: i, nl, nr, ll, lr, j, jR, jL, ml, mr, Nmax, indA, indB
    logical :: found
    real*8 :: CI_l, CI_r, CIi, CIj,Bele, Hele
    type(basis_sturmian_nr) :: temp_bst
    type(state), pointer :: state_jR

    if(state_l%nam /= state_r%nam) then
      print*, bst_nr%b(state_l%na)%l
      print*, bst_nr%b(state_r%na)%l
      !error stop '*** ERROR in sum_into_1el: l%nam /= r%nam'
    endif

    !if(any(state_l%ma /= state_r%ma)) error stop '*** ERROR in sum_into_1el: l%ma /= r%ma'

    ml = state_l%ma(1)
    mr = state_r%ma(1)
    if(any(state_l%ma /= ml)) error stop '*** different m in state_l'
    if(any(state_r%ma /= ml)) error stop '*** wrong m in state_r'
    
    do i=1, state_r%nam
      nr = state_r%na(i)
      lr = bst_nr%b(nr)%l
      CI_r = state_r%CI(i)
      found = .false.
      do j=1, state_l%nam
        nl = state_l%na(j)
        ll = bst_nr%b(nl)%l
        if(lr == ll) then
          found = .true.
          CI_l = state_l%CI(j)
          call copy(temp_sturm, bst_nr%b(nr))
          call multiply_byscalar(CI_r, temp_sturm)
          call multiply_byscalar(CI_l, bst_nr%b(nl))
          call sum_into(bst_nr%b(nl), temp_sturm)
          bst_nr%ortint(nl,:) = CI_l*bst_nr%ortint(nl,:)
          bst_nr%ortint(:,nl) = CI_l*bst_nr%ortint(:,nl)
          bst_nr%ortint(nl,:) = bst_nr%ortint(nl,:) + CI_r*bst_nr%ortint(nr,:)
          bst_nr%ortint(:,nl) = bst_nr%ortint(:,nl) + CI_r*bst_nr%ortint(:,nr)
          bst_nr%ham1el(nl,:) = CI_l*bst_nr%ham1el(nl,:)
          bst_nr%ham1el(:,nl) = CI_l*bst_nr%ham1el(:,nl)
          bst_nr%ham1el(nl,:) = bst_nr%ham1el(nl,:) + CI_r*bst_nr%ham1el(nr,:)
          bst_nr%ham1el(:,nl) = bst_nr%ham1el(:,nl) + CI_r*bst_nr%ham1el(:,nr)
        
          state_l%CI(j) = 1.0d0
         
          exit
        endif !lr == ll
      enddo !j
      if(.not.found) then
        !state_l does not have a function of this l
        Nmax = bst_nr%N
        call new_basis(temp_bst,Nmax)
        call copy_basis(temp_bst,bst_nr)
        call destruct(bst_nr)
        Nmax = Nmax + 1
        call new_basis(bst_nr,Nmax)
        do j=1, bst_nr%N-1
          call copy(bst_nr%b(j), temp_bst%b(j))
        enddo
        call copy(bst_nr%b(Nmax), bst_nr%b(nr))
        call multiply_byscalar(CI_r, bst_nr%b(Nmax))
        bst_nr%ortint(:,Nmax) = CI_r*bst_nr%ortint(:,nr)
        bst_nr%ortint(Nmax,:) = CI_r*bst_nr%ortint(nr,:)
        bst_nr%ham1el(:,Nmax) = CI_r*bst_nr%ham1el(:,nr)
        bst_nr%ham1el(Nmax,:) = CI_r*bst_nr%ham1el(nr,:)
        
      endif
    enddo !i
    
    !Now update the ovlpst and e1me matrices

    !First find the index of state_l in the TargetStates%b array
    do jL = 1, TargetStates%Nmax
      if(TargetStates%b(jL)%M == state_l%M .and. TargetStates%b(jL)%parity == state_l%parity .and. &
        &TargetStates%b(jL)%spin == state_l%spin .and. TargetStates%b(jL)%inum == state_l%inum) exit
    enddo
    if(jL > TargetStates%Nmax) error stop '*** ERROR in sum_into_1el: could not find index of state_l'

    !Now update the matrices - this code copied from spheroidal_12.f90
    do jR=1, TargetStates%Nmax
      state_jR => TargetStates%b(jR)
      bele = 0d0; Hele = 0d0
      do i = 1, get_nam(state_l)
         indA = get_na(state_l,i)
         CIi = get_CI(state_l,i)

         do j = 1, get_nam(state_jR)
            indB = get_na(state_jR,j)
            CIj = get_CI(state_jR,j)
            bele = bele + CIi*CIj*bst_nr%ortint(indA,indB)
            Hele = Hele + CIi*CIj*bst_nr%Ham1el(indA,indB)
         end do
      end do

      ovlpst(jL,jR) = bele; e1me(jL,jR) = Hele
      ovlpst(jR,jL) = bele; e1me(jR,jL) = Hele
    enddo


    
  end subroutine sum_into_1el
  
  subroutine sum_into_2el(state_l, state_r)
    use one_electron_func
    use input_data
    !Liam: this state sums two states state_l and state_r and stores the result in state_l
    !      (Hence state_r is summed INTO state_l)
    !
    !      state_l = state_l + state_r
    !
    !      basis should be the underlying basis in which the states are expanded
    !
    !      states_basis should be the basis containing the two states to be summed
    !
    implicit none

    type(state), intent(inout) :: state_l
    type(state), intent(in) :: state_r

    integer :: cl, cr !configuration indices for left and right states
    integer :: nconfigs, al, bl, ar, br, nco
    type(state) :: tmp_state
    type(state), pointer :: orb_l, orb_r
    logical :: found
    real*8 :: CI_l, CI_r
    logical, dimension(TargetStates%Nmax) :: core

    core = .false.
    core(1:TargetStates2el%nicm) = .true.

    if(data_in%theta /= 0.0d0) error stop '*** ERROR: the sum_into_2el subroutine does not yet work with theta /= 0.0'

    nconfigs = state_l%nam
    if(nconfigs /= state_r%nam) error stop '*** ERROR in sum_into_st: states do not have the same number of configurations'

    do cl=1, nconfigs
      
      al = state_l%na(cl); bl = state_l%nb(cl)

      if(al > bl) then          !any configurations with al > bl will reference the same 
        state_l%CI(cl) = 1.0d0  !underlying functions as the corresponding al < bl configuration
        cycle                   !so they should not be modified here except to change the CI coefficient to 1.0
      endif                     !to match what was done for the al < bl configuration.

      !Search for corresponding configuration in state_r
      found=.false.
      do cr=1, nconfigs
        ar = state_r%na(cr); br = state_r%nb(cr)
        !All matching configurations will have ar == al
        if(ar == al) then
          !If it is a (core)(core) configuration we expect to see an identical configuration
          !in state_r
          if(core(bl) .and. br == bl) then
            if(found) error stop '*** ERROR: multiple configuration matches'
            found = .true.
            !We simply sum the CI coefficients
            state_l%CI(cl) = state_l%CI(cl) + state_r%CI(cr)
         
            !write(*,'(2("(",I2,")")," matched to ",2("(",I2,")"))') al,bl,ar,br
            !If it is a (core)(not core) configuration then we need to find a (core)(not core) configuration
            !in state_r with the same A orbital and a different B orbital
          elseif(.not.core(bl) .and. br /= bl .and..not.core(br)) then
            if(found) error stop '*** ERROR: multiple configuration matches'
            found = .true.
            orb_l => TargetStates%b(bl)
            orb_r => TargetStates%b(br)
            CI_l = state_l%CI(cl); CI_r = state_r%CI(cr)
            call copy(tmp_state, orb_r)
            call multiply_by_scalar(tmp_state, CI_r)
            call multiply_by_scalar(orb_l, CI_l)
            call sum_into_1el(orb_l, tmp_state)
            state_l%CI(cl) = 1.0d0
            !write(*,'(2("(",I2,")")," matched to ",2("(",I2,")"))') al,bl,ar,br
          endif

        endif !ar == al

      enddo !cr

      if(.not.found) then
        print*, '*** ERROR in sum_into_2el: unmatched configuration'
        print*, 'state_l:', state_l%label
        print*, 'config, na, nb:', cl, al, bl
      endif !found
      
      !Once the corresponding configuration has been found in state_r and they have been summed, the
      !CI coefficient in state_l for this configuration will have been absorbed into the one-electron 
      !orbital so now we set it to 1.0
    enddo !cl

  end subroutine sum_into_2el
  
  subroutine sum_into_2el_b4_rearrange(state_l, state_r)
    use one_electron_func
    use input_data
    !Same as sum_into_2el but for the state representations before rearrange12() is called
    implicit none

    type(state), intent(inout) :: state_l
    type(state), intent(in) :: state_r

    integer :: cl, cr !configuration indices for left and right states
    integer :: al, bl, ar, br, nam
    type(state) :: tmp_state
    type(state), pointer :: orb_l, orb_r
    logical :: found
    real*8 :: CI_l, CI_r
    real*8, dimension(:), allocatable :: CI
    integer, dimension(:), allocatable :: na, nb, ma, mb, phase

    !print*, 'BEFORE SUM:'
    !print*, 'STATE_R%na:', state_r%na
    !print*, 'STATE_R%nb:', state_r%nb
    !print*, 'STATE_R%CI:', state_r%CI
    !print*, 'STATE_L%na:', state_l%na
    !print*, 'STATE_L%nb:', state_l%nb
    !print*, 'STATE_L%CI:', state_l%CI
    !print*, ''


    do cr=1, state_r%nam
      ar = state_r%na(cr); br = state_r%nb(cr)
     
      !print*, 'ar, br:', ar, br

      !Search for corresponding configuration in state_l
      do cl=1, state_l%nam
        al = state_l%na(cl); bl = state_l%nb(cl)
        !print*, '   -> al, bl:', al, bl
        if(al == ar .and. bl == br) exit
      enddo
      if(cl == state_l%nam+1) then !didn't find corresponding configuration, add it:
        nam = state_l%nam
        allocate(CI(nam+1),na(nam+1),nb(nam+1),ma(nam+1),mb(nam+1),phase(nam+1))
        CI(1:nam) = state_l%CI; CI(nam+1) = 0.0d0
        na(1:nam) = state_l%na; na(nam+1) = ar
        nb(1:nam) = state_l%nb; nb(nam+1) = br
        ma(1:nam) = state_l%ma; ma(nam+1) = state_r%ma(cr) 
        mb(1:nam) = state_l%mb; mb(nam+1) = state_r%mb(cr)
        phase = (-1)**state_l%spin
        !print*, 'ADDING CONFIG', ar, br
        !print*, '  -> na:', na
        !print*, '     nb:', nb
        call construct_st_(tmp_state, .false., state_l%m, state_l%parity, state_l%spin, state_l%energy, state_l%inum, &
          &nam+1, CI, na, ma, nb, mb,phase)
        !print*, '  -> na:', tmp_state%na
        !print*, '     nb:', tmp_state%nb
        deallocate(CI,na,nb,ma,mb,phase)
        call destruct_st(state_l)
        call copy_st(state_l,tmp_state)
        call destruct_st(tmp_state)
        !print*, '  -> na:', state_l%na
        !print*, '     nb:', state_l%nb
      endif
      !at this point cl is the index of the corresponding configuration in state_l
      state_l%CI(cl) = state_l%CI(cl) + state_r%CI(cr)
      !print*, ''

    enddo
    
    !print*, 'AFTER SUM:'
    !print*, 'STATE_L%na:', state_l%na
    !print*, 'STATE_L%nb:', state_l%nb
    !print*, 'STATE_L%CI:', state_l%CI

  end subroutine sum_into_2el_b4_rearrange
    
  subroutine test_summing_routines
    use one_electron_func
    use grid_radial
    use sturmian_class
    use input_data
    implicit none
    type(state), pointer :: state_l, state_r
    type(sturmian_nr), pointer :: sturm_l, sturm_r, pn_i, pn_f
    type(sturmian_nr) :: temp_sturm
    type(state) :: temp_state
    type(state), pointer ::stateR, stateL, gs, onestate_i, onestate_f, st
    integer :: l, n, i1, i2, M, P, S, i, a, b, nicm
    real*8 :: result_l, result_v, rl, rv, rl1, rv1, CIi, CIj, bele, Hele, dE, sum_new, sum_two, rCI_i, rCI_f
    integer :: jL, jR, nFns, indA, indB, j, nsti, ncm_i, msti, icl_i, licl_i, nicl_i, ncm_f, icl_f, licl_f, nicl_f, n_i, lfn_i,lfn_f
    integer :: Nmax1el, Nmax2el, numcoreMO, numcoreLag, numFCMO, numFCLag, dM, itest
    character*6 :: coretype
    character*8 :: orbtype
    integer, dimension(3,-1:1, 100) :: orbs_for_test

!    !First test sturmiam sum
!    sturm_l => bst_nr%b(1)
!
!    l = sturm_l%l
!    do i=2, bst_nr%n
!      if (bst_nr%b(i)%l == l) then
!        sturm_r => bst_nr%b(i)
!        exit
!      endif
!    enddo
!    !now sturm_r has the same l as sturm_l

!    call copy(temp_sturm, sturm_l)

!    call sum_into(sturm_l, sturm_r)
!
!    open(unit=600,file='sturm_sum.test',action='write',status='replace')
!    do i=sturm_l%minf, sturm_l%maxf
!      write(600,*) grid%gridr(i), temp_sturm%f(i), sturm_r%f(i), sturm_l%f(i)
!    enddo 
!    close(600)
!    !!! SUM_INTO IS WORKING

    Nmax1el = TargetStates%Nmax
    Nmax2el = TargetStates2el%Nmax
    nicm = TargetStates2el%nicm

!    !Check there are enough orbitals of each type to run test
!    !First print them out
!    print*, '===== ONE ELECTRON ORBITALS ====='
!    do i=1, Nmax1el
!      st => TargetStates%b(i)
!      if(i <= nicm) then
!        coretype = '(CORE)'
!      else
!        coretype = '( FC )'
!      endif
!      if(st%energy == 0.0d0) then
!        orbtype = 'Laguerre'
!      else
!        orbtype = '   MO   '
!      endif
!      print*, i, coretype//' '//orbtype//' ', nint(st%M), st%parity
!    enddo
!    print*, ''
!  
!    orbs_for_test = 0 
!    do i=1, nicm
!      if(TargetStates%b(i)%energy /= 0.0d0) exit
!    enddo
!    if(i>nicm) error stop '*** did not find a core MO'
!    print*, 'For testing coreMO+coreMO sum, initial state will be',i
!    orbs_for_test(1,:,1) = i
!    gs => TargetStates%b(i)
!    M = gs%M
!    P = gs%parity
!    n = 0
!    print*, 'Looking for 2 M=', M, 'core orbitals'
!    do i=1, nicm
!      if(n==2) exit
!      if(TargetStates%b(i)%energy == 0.0d0) cycle
!      st => TargetStates%b(i)
!      if(st%M == M .and. st%parity == -P) then
!        n = n + 1
!        print*, '  > found state', i
!        orbs_for_test(n+1,0,1) = i
!      endif
!    enddo
!    if(n < 2) error stop '*** not enough ∆M=0 core orbitals for testing' 
!    n = 0
!    print*, 'Looking for 2 M=', M+1, 'core orbitals'
!    do i=1, nicm
!      if(n==2) exit
!      if(TargetStates%b(i)%energy == 0.0d0) cycle
!      st => TargetStates%b(i)
!      if(st%M == M+1 .and. st%parity == -P) then
!        n = n + 1
!        print*, '  > found state', i
!        orbs_for_test(n+1,1,1) = i
!      endif
!    enddo
!    if(n < 2) error stop '*** not enough ∆M=1 core orbitals for testing' 
!    n = 0
!    print*, 'Looking for 2 M=', M-1, 'core orbitals'
!    do i=1, nicm
!      if(n==2) exit
!      if(TargetStates%b(i)%energy == 0.0d0) cycle
!      st => TargetStates%b(i)
!      if(st%M == M-1 .and. st%parity == -P) then
!        n = n + 1
!        print*, '  > found state', i
!        orbs_for_test(n+1,-1,1) = i
!      endif
!    enddo
!    if(n < 2) error stop '*** not enough ∆M=-1 core orbitals for testing' 
!    
!    do i=1, nicm
!      if(TargetStates%b(i)%energy == 0.0d0) exit
!    enddo
!    if(i>nicm) error stop '*** did not find a core Laguerre'
!    print*, 'For testing coreLag+coreLag sum, initial state will be',i
!    orbs_for_test(1,:,2) = i
!    gs => TargetStates%b(i)
!    M = gs%M
!    P = gs%parity
!    n = 0
!    print*, 'Looking for 2 M=', M, 'core Laguerre orbitals'
!    do i=1, nicm
!      if(n==2) exit
!      if(TargetStates%b(i)%energy /= 0.0d0) cycle
!      st => TargetStates%b(i)
!      if(st%M == M .and. st%parity == -P) then
!        n = n + 1
!        print*, '  > found state', i
!        orbs_for_test(n+1,0,2) = i
!      endif
!    enddo
!    if(n < 2) error stop '*** not enough ∆M=0 core Laguerre orbitals for testing' 
!    n = 0
!    print*, 'Looking for 2 M=', M+1, 'core Laguerre orbitals'
!    do i=1, nicm
!      if(n==2) exit
!      if(TargetStates%b(i)%energy /= 0.0d0) cycle
!      st => TargetStates%b(i)
!      if(st%M == M+1 .and. st%parity == -P) then
!        n = n + 1
!        print*, '  > found state', i
!        orbs_for_test(n+1,1,2) = i
!      endif
!    enddo
!    if(n < 2) error stop '*** not enough ∆M=1 core Laguerre orbitals for testing' 
!    n = 0
!    print*, 'Looking for 2 M=', M-1, 'core Laguerre orbitals'
!    do i=1, nicm
!      if(n==2) exit
!      if(TargetStates%b(i)%energy /= 0.0d0) cycle
!      st => TargetStates%b(i)
!      if(st%M == M-1 .and. st%parity == -P) then
!        n = n + 1
!        print*, '  > found state', i
!        orbs_for_test(n+1,-1,2) = i
!      endif
!    enddo
!    if(n < 2) error stop '*** not enough ∆M=-1 core Laguerre orbitals for testing' 
!    
!    do i=nicm+1, Nmax1el
!      if(TargetStates%b(i)%energy == 0.0d0) exit
!    enddo
!    if(i>Nmax1el) error stop '*** did not find a FC Laguerre'
!    print*, 'For testing FCLag+FCLag sum, initial state will be',i
!    orbs_for_test(1,:,3) = i
!    gs => TargetStates%b(i)
!    M = gs%M
!    P = gs%parity
!    n = 0
!    print*, 'Looking for 2 M=', M, 'FC Laguerre orbitals'
!    do i=nicm, Nmax1el
!      if(n==2) exit
!      if(TargetStates%b(i)%energy /= 0.0d0) cycle
!      st => TargetStates%b(i)
!      if(st%M == M .and. st%parity == -P) then
!        n = n + 1
!        print*, '  > found state', i
!        orbs_for_test(n+1,0,3) = i
!      endif
!    enddo
!    if(n < 2) error stop '*** not enough ∆M=0 FC Laguerre orbitals for testing' 
!    n = 0
!    print*, 'Looking for 2 M=', M+1, 'FC Laguerre orbitals'
!    do i=nicm, Nmax1el
!      if(n==2) exit
!      if(TargetStates%b(i)%energy /= 0.0d0) cycle
!      st => TargetStates%b(i)
!      if(st%M == M+1 .and. st%parity == -P) then
!        n = n + 1
!        print*, '  > found state', i
!        orbs_for_test(n+1,1,3) = i
!      endif
!    enddo
!    if(n < 2) error stop '*** not enough ∆M=1 FC Laguerre orbitals for testing' 
!    n = 0
!    print*, 'Looking for 2 M=', M-1, 'FC Laguerre orbitals'
!    do i=nicm, Nmax1el
!      if(n==2) exit
!      if(TargetStates%b(i)%energy /= 0.0d0) cycle
!      st => TargetStates%b(i)
!      if(st%M == M-1 .and. st%parity == -P) then
!        n = n + 1
!        print*, '  > found state', i
!        orbs_for_test(n+1,-1,3) = i
!      endif
!    enddo
!    if(n < 2) error stop '*** not enough ∆M=-1 FC Laguerre orbitals for testing' 
!    
!    do i=1, nicm
!      if(TargetStates%b(i)%energy /= 0.0d0) exit
!    enddo
!    if(i>Nmax1el) error stop '*** did not find a core MO'
!    print*, 'For testing core MO+FCLag sum, initial state will be',i
!    orbs_for_test(1,:,4) = i
!    gs => TargetStates%b(i)
!    M = gs%M
!    P = gs%parity
!    n = 0
!    print*, 'Looking for 2 M=', M, 'FC Laguerre orbitals'
!    do i=nicm, Nmax1el
!      if(n==2) exit
!      if(TargetStates%b(i)%energy /= 0.0d0) cycle
!      st => TargetStates%b(i)
!      if(st%M == M .and. st%parity == -P) then
!        n = n + 1
!        print*, '  > found state', i
!        orbs_for_test(n+1,0,4) = i
!      endif
!    enddo
!    if(n < 2) error stop '*** not enough ∆M=0 FC Laguerre orbitals for testing' 
!    n = 0
!    print*, 'Looking for 2 M=', M+1, 'FC Laguerre orbitals'
!    do i=nicm, Nmax1el
!      if(n==2) exit
!      if(TargetStates%b(i)%energy /= 0.0d0) cycle
!      st => TargetStates%b(i)
!      if(st%M == M+1 .and. st%parity == -P) then
!        n = n + 1
!        print*, '  > found state', i
!        orbs_for_test(n+1,1,4) = i
!      endif
!    enddo
!    if(n < 2) error stop '*** not enough ∆M=1 FC Laguerre orbitals for testing' 
!    n = 0
!    print*, 'Looking for 2 M=', M-1, 'FC Laguerre orbitals'
!    do i=nicm, Nmax1el
!      if(n==2) exit
!      if(TargetStates%b(i)%energy /= 0.0d0) cycle
!      st => TargetStates%b(i)
!      if(st%M == M-1 .and. st%parity == -P) then
!        n = n + 1
!        print*, '  > found state', i
!        orbs_for_test(n+1,-1,4) = i
!      endif
!    enddo
!    if(n < 2) error stop '*** not enough ∆M=-1 FC Laguerre orbitals for testing' 
!
!
!    do itest=1, 4
!      print*, ''
!      if(itest == 1) then
!        print*, '---- TESTING coreMO+coreMO sum ----'
!      elseif(itest == 2) then
!        print*, '---- TESTING coreLag+coreLag sum ----'
!      elseif(itest == 3) then
!        print*, '---- TESTING FCLag+FCLag sum ----'
!      endif
!
!      do dM=-1, 1 
!        print*, '-> dM = ', dM
!        i = orbs_for_test(1,dM,itest)
!        i1 = orbs_for_test(2,dM,itest)
!        i2 = orbs_for_test(3,dM,itest)
!        i1 = orbs_for_test(3,dM,itest)
!        i2 = orbs_for_test(2,dM,itest)
!        gs => TargetStates%b(i)
!        state_l => TargetStates%b(i1)
!        state_r => TargetStates%b(i2)
!   
!        call multiply_by_scalar(state_l, 2.0d0)
!
!        dE = 1.0d0
!        call oscstr1_st(0,i1,i,state_l,gs,dE,result_l,result_v,rl1,rv1)
!        print*, 'DIPOLE matrix element between state ',i,' and ', i1, ':', rl1, rv1
!        
!        call oscstr1_st(0,i2,i,state_r,gs,dE,result_l,result_v,rl,rv)
!        print*, 'DIPOLE matrix element between state ',i,' and ', i2, ':', rl, rv
!        rl1 = rl1 + rl
!        rv1 = rv1 + rv
!        print*, 'SUMMED DIPOLE matrix elements:                                ', rl1, rv1
!        
!        print*, 'Now state',i2,'summed into state',i1
!        call sum_into_1el(state_l, state_r)
!        
!        call oscstr1_st(0,i1,i,state_l,gs,dE,result_l,result_v,rl,rv)
!        print*, 'DIPOLE matrix element between state ',i,' and ', i1, ':', rl, rv
!        if(abs((rl-rl1)/rl1) < 1e-5) then
!          print*, '****** TEST PASSED! *******'
!        else
!          error stop '****** TEST FAILED!!!!!! *******'
!        endif
!        print*, ''
!
!      enddo !dM
!    enddo !itest

!! SUM_INTO_1EL IS WORKING !




  gs => TargetStates2el%b(1)
  state_l => TargetStates2el%b(2)
  state_r => TargetStates2el%b(3)

!  print*, 'BEFORE SUM:'
!  print*, 'gs:'
!  do i=1, gs%nam
!    a = gs%na(i)
!    b = gs%nb(i)
!    print*, 'CONFIG, A, B:', i, a, b
!    print*, 'BST A:', TargetStates%b(a)%na
!    print*, 'SUM f:', (sum(bst_nr%b(TargetStates%b(a)%na(j))%f), j=1, TargetStates%b(a)%nam)
!    print*, 'BST B:', TargetStates%b(b)%na
!    print*, 'SUM f:', (sum(bst_nr%b(TargetStates%b(b)%na(j))%f), j=1, TargetStates%b(b)%nam)
!    print*, ''
!  enddo
!  print*, 'state_l:'
!  do i=1, state_l%nam
!    a = state_l%na(i)
!    b = state_l%nb(i)
!    print*, 'CONFIG, A, B:', i, a, b
!    print*, 'BST A:', TargetStates%b(a)%na
!    print*, 'SUM f:', (sum(bst_nr%b(TargetStates%b(a)%na(j))%f), j=1, TargetStates%b(a)%nam)
!    print*, 'BST B:', TargetStates%b(b)%na
!    print*, 'SUM f:', (sum(bst_nr%b(TargetStates%b(b)%na(j))%f), j=1, TargetStates%b(b)%nam)
!    print*, ''
!  enddo
!  print*, 'state_r:'
!  do i=1, state_r%nam
!    a = state_r%na(i)
!    b = state_r%nb(i)
!    print*, 'CONFIG, A, B:', i, a, b
!    print*, 'BST A:', TargetStates%b(a)%na
!    print*, 'SUM f:', (sum(bst_nr%b(TargetStates%b(a)%na(j))%f), j=1, TargetStates%b(a)%nam)
!    print*, 'BST B:', TargetStates%b(b)%na
!    print*, 'SUM f:', (sum(bst_nr%b(TargetStates%b(b)%na(j))%f), j=1, TargetStates%b(b)%nam)
!    print*, ''
!  enddo
!
!!  do i=1, bst_nr%n
!!    write(*,'(1000(ES8.1,X))') bst_nr%ortint(i,:)
!!  enddo
!!  print*, ''
!!  call sum_into_2el(state_l, state_r)
!!  
!!  print*, 'AFTER SUM:'
!!  do i=1, bst_nr%n
!!    write(*,'(1000(ES8.1,X))') bst_nr%ortint(i,:)
!!  enddo
!!  print*, ''
!  do i=1, state_l%nam
!    a = state_l%na(i)
!    b = state_l%nb(i)
!    print*, 'CONFIG, A, B:', i, a, b
!    print*, 'BST A:', TargetStates%b(a)%na
!    print*, 'SUM f:', (sum(bst_nr%b(TargetStates%b(a)%na(j))%f), j=1, TargetStates%b(a)%nam)
!    print*, 'BST B:', TargetStates%b(b)%na
!    print*, 'SUM f:', (sum(bst_nr%b(TargetStates%b(b)%na(j))%f), j=1, TargetStates%b(b)%nam)
!    print*, ''
!  enddo
!



    !Now test two-electron sum
    M = TargetStates2el%b(1)%M
    P = TargetStates2el%b(1)%parity
    S = nint(2*TargetStates2el%b(1)%spin)
    do i1=2, TargetStates2el%Nmax
      if(abs(TargetStates2el%b(i1)%M - M) <= 1 .and. TargetStates2el%b(i1)%parity == -P .and. nint(TargetStates2el%b(i1)%spin) == S) then
        state_l => TargetStates2el%b(i1)
        exit
      endif
    enddo
    do i2=i1+1, TargetStates2el%Nmax
      if(TargetStates2el%b(i2)%M == TargetStates2el%b(i1)%M .and. TargetStates2el%b(i2)%parity == -P .and. nint(TargetStates2el%b(i2)%spin) == S) then
        state_r => TargetStates2el%b(i2)
        exit
      endif
    enddo

    print*, 'state_l:', state_l%label, i1
    print*, 'state_r:', state_r%label, i2
  
  !  state_l => TargetStates2el%b(1)

  print*, gs%label
    do i=1, gs%nam
      a = gs%na(i)
      b = gs%nb(i)
      nicm = TargetStates2el%nicm
      if(a <= nicm .and. b <= nicm) print*, i, '(CORE)(CORE)', a,b, '|'//TargetStates%b(a)%label//'|', '|'//TargetStates%b(b)%label//'|'
      if(a <= nicm .and. b > nicm) print*, i, '(CORE)(NOT CORE)', a,b, '|'//TargetStates%b(a)%label//'|', '|'//TargetStates%b(b)%label//'|'
      if(b <= nicm .and. a > nicm) print*, i, '(NOT CORE)(CORE)', a,b, '|'//TargetStates%b(a)%label//'|', '|'//TargetStates%b(b)%label//'|'
     ! if(a > nicm .or. b > nicm) gs%CI(i) = 0.0d0
    enddo
  print*, state_l%label
    do i=1, state_l%nam
      a = state_l%na(i)
      b = state_l%nb(i)
      nicm = TargetStates2el%nicm
      if(a <= nicm .and. b <= nicm) print*, i, '(CORE)(CORE)', a,b, '|'//TargetStates%b(a)%label//'|', '|'//TargetStates%b(b)%label//'|'
      if(a <= nicm .and. b > nicm) print*, i, '(CORE)(NOT CORE)', a,b, '|'//TargetStates%b(a)%label//'|', '|'//TargetStates%b(b)%label//'|'
      if(b <= nicm .and. a > nicm) print*, i, '(NOT CORE)(CORE)', a,b, '|'//TargetStates%b(a)%label//'|', '|'//TargetStates%b(b)%label//'|'
     ! if(a > nicm .or. b > nicm) state_l%CI(i) = 0.0d0
    enddo
    print*, state_r%label
    do i=1, state_r%nam
      a = state_r%na(i)
      b = state_r%nb(i)
      nicm = TargetStates2el%nicm
      if(a <= nicm .and. b <= nicm) print*, i, '(CORE)(CORE)', a,b, '|'//TargetStates%b(a)%label//'|', '|'//TargetStates%b(b)%label//'|'
      if(a <= nicm .and. b > nicm) print*, i, '(CORE)(NOT CORE)', a,b, '|'//TargetStates%b(a)%label//'|', '|'//TargetStates%b(b)%label//'|'
      if(b <= nicm .and. a > nicm) print*, i, '(NOT CORE)(CORE)', a,b, '|'//TargetStates%b(a)%label//'|', '|'//TargetStates%b(b)%label//'|'
     ! if(a > nicm .or. b > nicm) state_r%CI(i) = 0.0d0
    enddo
    
    print*, 'state_l%nb:', state_l%nb


    if(maxval(state_l%na) > TargetStates%Nmax) then
      print*, 'label:', state_l%label
      print*, 'nam:', state_l%nam
      print*, 'na:', state_l%na
      print*, 'maxval:', maxval(state_l%na)
      print*, 'NMAX:', targetstates%nmax
      print*, 'nusemax:', state_l%nusemax
      print*, 'nuse:', state_l%nuse
      error stop
    elseif(maxval(state_r%na) > TargetStates%Nmax) then
      print*, 'label:', state_r%label
      print*, 'nam:', state_r%nam
      print*, 'na:', state_r%na
      print*, 'maxval:', maxval(state_r%na)
      print*, 'NMAX:', targetstates%nmax
      print*, 'nusemax:', state_r%nusemax
      print*, 'nuse:', state_r%nuse
      error stop
    endif

    call oscstrength_2e_MOrep(1,i1,1.0d0,result_l,result_v,rl1,rv1)
    print*, 'DIPOLE matrix element between state 1 and ', i1, ':', rl1, rv1
    dE = state_r%energy - TargetStates%b(1)%energy
    dE = 1.0d0
    call oscstrength_2e_MOrep(1,i2,1.0d0,result_l,result_v,rl,rv)
    print*, 'DIPOLE matrix element between state 1 and ', i2, ':', rl, rv
    print*, 'SUMMED DIPOLE matrix elements:                                ', rl1+rl, rv1+rv
    call sum_into_2el(state_l, state_r)

    print*, 'Now state',i2,'summed into state',i1
    call oscstrength_2e_MOrep(1,i1,1.0d0,result_l,result_v,rl,rv)
    print*, 'DIPOLE matrix element between state 1 and ', i1, ':', rl, rv

    print*, ''
    print*, 'NOW A SECOND TIME:'
    call oscstrength_2e_MOrep(1,i1,1.0d0,result_l,result_v,rl1,rv1)
    print*, 'DIPOLE matrix element between state 1 and ', i1, ':', rl1, rv1
    dE = state_r%energy - TargetStates%b(1)%energy
    dE = 1.0d0
    call oscstrength_2e_MOrep(1,i2,1.0d0,result_l,result_v,rl,rv)
    print*, 'DIPOLE matrix element between state 1 and ', i2, ':', rl, rv
    print*, 'SUMMED DIPOLE matrix elements:                                ', rl1+rl, rv1+rv

    call sum_into_2el(state_l, state_r)

    print*, 'Now state',i2,'summed into state',i1
    call oscstrength_2e_MOrep(1,i1,1.0d0,result_l,result_v,rl,rv)
    print*, 'DIPOLE matrix element between state 1 and ', i1, ':', rl, rv
   
    print*, ''
    print*, 'Now we multiply state', i1, ' by 3'
    print*, 'and multiply state', i2, ' by 7'
    call multiply_by_scalar(TargetStates2el%b(i1), 3.0d0)
    call multiply_by_scalar(TargetStates2el%b(i2), 7.0d0)
    call oscstrength_2e_MOrep(1,i1,1.0d0,result_l,result_v,rl1,rv1)
    print*, 'DIPOLE matrix element between state 1 and ', i1, ':', rl1, rv1
    dE = state_r%energy - TargetStates%b(1)%energy
    dE = 1.0d0
    call oscstrength_2e_MOrep(1,i2,1.0d0,result_l,result_v,rl,rv)
    print*, 'DIPOLE matrix element between state 1 and ', i2, ':', rl, rv
    print*, 'SUMMED DIPOLE matrix elements:                                ', rl1+rl, rv1+rv
    call sum_into_2el(state_l, state_r)

    print*, 'Now state',i2,'summed into state',i1
    call oscstrength_2e_MOrep(1,i1,1.0d0,result_l,result_v,rl,rv)
    print*, 'DIPOLE matrix element between state 1 and ', i1, ':', rl, rv

    stop


      
  end subroutine test_summing_routines
  
  subroutine  construct_st_(self,hlike,m,parity,spin,energy,inum,ncm,CI,no1,mo1,no2,mo2,phase)

    implicit none

    type(state), intent(inout):: self
    logical, intent(in):: hlike
    real*8, intent(in):: m  
    integer, intent(in):: parity 
    real*8, intent(in):: spin
    real*8, intent(in):: energy
    integer, intent(in):: inum
    integer, intent(in):: ncm
    real*8, dimension(ncm), intent(in) :: CI
    integer, dimension(ncm), intent(in) :: no1, no2, phase, mo1, mo2
    optional:: no2, phase, mo2
    integer:: i

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
          if ( self%m /= self%ma(i) ) then
             print*,"Error construct_st_: angular projection for target state m and mo1 different"
             print*, 'self%m, self%ma(i) =', self%m, self%ma(i)
             stop
          end if
       end do
       self%CI(1:ncm) = CI(1:ncm)
    else
       self%label = "  -  "
       allocate(self%CI(ncm),self%na(ncm),self%ma(ncm),self%nb(ncm),self%mb(ncm))
       self%nam=ncm
       self%na(1:ncm) = no1(1:ncm)
       self%ma(1:ncm) = mo1(1:ncm)
       self%nb(1:ncm) = no2(1:ncm)
       self%mb(1:ncm) = mo2(1:ncm)
       do i = 1, ncm
          if ( self%m /= self%ma(i)+self%mb(i) ) then
             print*,"Error construct_st_: angular projection for target state m and mo1+mo2 different"
             print*, 'self%m, self%ma(i), self%mb(i) =', self%m, self%ma(i), self%mb(i)
             stop
          end if
       end do
       self%CI(1:ncm) = CI(1:ncm)
    endif

!    print*,'created a state with energy:', energy
   
  end   subroutine  construct_st_


end module target_states
