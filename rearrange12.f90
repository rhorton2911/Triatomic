subroutine rearrange12_st(TargetStates1el, TargetStates2el,nicm,Nmax,Nmax1el)

  use sturmian_class
  use state_class
  use one_electron_func
  use grid_radial
  use ovlpste1me
  use MPI_module
  use input_data 
  
  
  implicit none
  
  type(basis_state):: TargetStates1el, TargetStates2el 
  integer, intent(in):: nicm,Nmax,Nmax1el
  
  integer:: nst, ma, inc, n1, n2, i, istcore, nstcore, nam, nst1el, ipar, j, nd, ncm, is, nstp, nspbst
  real*8:: energy, rma, tmp
!!$ representation of f.c. 1el states via orginal one-electron states
  integer, dimension(:,:), allocatable:: igiven
  real*8, dimension(:,:), allocatable:: CIigiven
  integer, dimension(:), allocatable:: igiven_max, iTmp
  
!!$ representation of 2el states via f.c. one-electron states
  integer, dimension(:), allocatable:: manst1el, iparnst1el 
  integer, dimension(Nmax):: nstton1
  integer, dimension(Nmax,nicm):: nstton2
  integer, dimension(Nmax,nicm):: nstA
  real*8, dimension(Nmax,nicm):: CInstA
  logical iscore, islarger, itmplog
  type(basis_state):: tmp1elst
  
  integer, dimension(:), allocatable:: no1,no2,mo1,mo2, phase
  real*8, dimension(2*Nmax1el):: CI
  integer, dimension(:), allocatable:: nobst, mobst
  real*8, dimension(:), allocatable:: CInobst, CITmp
  integer:: itmpmax
!!$ Mark: Addition solving non-uniqueness
  integer, dimension(:,:), allocatable:: nst_coreorb_cfig 
  integer, dimension(nicm):: ind_coreorb
  integer:: n_coreorb_cfigs, n_new_cfigs
  
  interface
    subroutine rearrange(bst_nr,LMAX,       TargetStates,use_lag_ham1el,lag_ham1el_m)
      use input_data
      use grid_radial
      use sturmian_class
      use state_class
      use MPI_module
      type(basis_sturmian_nr), intent(inout) :: bst_nr 
      integer, intent(in):: Lmax   ! max L valus in one-electron basis bst
      type(basis_state), intent(inout):: TargetStates
      logical, intent(in):: use_lag_ham1el ! switch use lag_ham1el_m - Laguerre basis <Lp|H|L> 
      real*8, dimension(bst_nr%n,bst_nr%n,-Lmax:Lmax), intent(in):: lag_ham1el_m  ! 
      optional:: lag_ham1el_m
    end subroutine
  end interface
  
  nstton1(:) = 0 
  nstton2(:,:) = 0
  CI(:) = 0d0
  
  CInstA(:,:) = 0d0
  nstA(:,:) = 0
!!$ Mark: Addition solving non-uniqueness
  n_coreorb_cfigs = 0
  ind_coreorb(:) = TargetStates2el%ncore(:)
 if ( data_in%non_uniq) then
!!$ Addition to not include core orbitals into FCO. 
!!$ nst_coreorb_cfig = number of symmetric configurations that have 
!!$ "inner" and "outer" core orbital that aren't symmetric
!!$ for each 2e-state and inner core orbital
     allocate(nst_coreorb_cfig(Nmax,nicm) )
     nst_coreorb_cfig(:,:) = 0
  end if
!!$>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!!$  This part of the code deal with building frozen-core representation using target states of H2+
!!$      
!!$   \Psi(1,2)  =  \sum_{icore} \psi^{core}_{icore}(1) \psi^{fccore}_{icore}(2)
!!$
!!$  first determine the number of frozen-core (f.c.) orbitals (H2+ states actually)
  
!!$ nst1el is the number of frozen core (molecular) orbitals + core (molecular) orbitals
  nst1el = nicm
  
  do nst=1, Nmax  ! nst - loop over number of states
     
     do istcore=1,nicm  ! go through the core one-electron  molecular orbitals
        
        nstcore = TargetStates2el%ncore(istcore) ! index to core (Laguerre or molecular) orbital/1el state
        
        nam = get_nam(TargetStates2el%b(nst)) ! Number of configurations for 2el states nst
        
        i = 0 ! start the counter for f.c. 1el. state representation here
        j = 0 ! to deal with symetric configuration
        do inc = 1, nam
           n1 = get_na(TargetStates2el%b(nst),inc,1)  ! index to one-electron TargetStates
           if(n1 .ne. nstcore) cycle   ! Build frozen core orbital for each core orbital
           
           n2 = get_na(TargetStates2el%b(nst),inc,2) 
!!$ At this point n2 is index to outer excited electron orbital or
!!$ another core obital 
           
!!$ We only want to deal with symmetric configurations |ab>
!!$ then antisymmetrise later |ba>.
!!$ Say we have symmetric config |ab>, then |ba> is a possibility.
!!$ We know it is |ba> if n2 index comes before n1 index and 
!!$ n2 index is a core orbital
!!$ i.e. need to cycle if (n1 > n2 .AND n2 == any(TargetStates2el%ncore(:)) )
           
!!$ itmplog =.true. if n2 index to a core (Laguerre or molecular) orbital/1el state 
!!$ and n2 > n1, note n1 must be a core orbital 
!!$ .OR. n2 not index to a core orbital 
           itmplog = islarger(n2,istcore)
           
           if( n1 .eq. n2) then ! symmetric configuration
              j = n1
!!$              print*,"n1, n2",n1,n2," Sym"
              cycle
           elseif(itmplog) then 
!!$ Mark: Addition solving non-uniqueness
!!$ All core orbitals (type 2 and 1) should not be included inside the
!!$ rearranged frozen core orbital.
!!$ Frozen core orbitals not created for n2 core orbital 
              if ( data_in%non_uniq ) then 
                 if ( any(n2 == ind_coreorb(:)) ) then 
!!$ Don't create FCO for this configuration n1 and n2 are core orbitals
!!$ Note: Symmetric configurations (n1=n2) won't be here
                    n_coreorb_cfigs = n_coreorb_cfigs + 1
!!$                    print*,"n1, n2",n1,n2," Core"
                    cycle 
                 end if
              end if ! Non-uniqueness being solved
              ! contniue to build new f.c. state
!!$              print*,"n1, n2",n1,n2," FCO"
           else
              cycle ! to avoid double counting of configurations 
           endif
           
           ! f.c. orbital representation is done below
           ! only states n2 with the same value of magnetic quantum number  and parity will be here           
           ! need only to count the first occurance to add a new f.c. 1el state
           if(i .eq. 0 ) then
              nst1el = nst1el + 1  
           endif
           
           i = i + 1
           
        enddo ! nam - loop over 2e-state configurations     
     enddo ! istcore - loop over core one-electron orbitals
  enddo ! nst - loop over number of states
  
  allocate(igiven_max(nst1el),CIigiven(nst1el,Nmax1el),igiven(nst1el,Nmax1el))
  igiven_max(:) = 0
  igiven(:,:) = 0
  CIigiven(:,:) = 0d0
  allocate(manst1el(nst1el), iparnst1el(nst1el))
  manst1el(:) = 0
  iparnst1el(:) = 0
!!$ Mark changes for non-uniqueness two-electron configs with two core orbitals
!  allocate(no1(nst1el),no2(nst1el),mo1(nst1el),mo2(nst1el), phase(nst1el))
  n_new_cfigs = nst1el + n_coreorb_cfigs
  allocate(no1(n_new_cfigs),no2(n_new_cfigs),mo1(n_new_cfigs),mo2(n_new_cfigs), phase(n_new_cfigs))
  no1(:) = 0
  no2(:) = 0
  mo1(:) = 0
  mo2(:) = 0
  phase(:) = 0
  
  if (myid==0) print '("rearrange12(): Nmax,Nmax1el,nicm,nst1el: ",4i7)',   Nmax, Nmax1el, nicm, nst1el
  
  do istcore=1,nicm  ! go through the core one-electron states
     
     nstcore = TargetStates2el%ncore(istcore)
     manst1el(istcore) = get_ang_mom_proj(TargetStates1el%b(nstcore))
     iparnst1el(istcore) = get_par_st(TargetStates1el%b(nstcore))
     
  enddo
  
  nst1el = nicm
  
  do nst = 1, Nmax
     
     do istcore = 1, nicm  ! go through the core one-electron (molecular) orbitals 
!!$       if (myid==0)  print*
        nstcore = TargetStates2el%ncore(istcore)
        
        nam = get_nam(TargetStates2el%b(nst))
        
        i = 0 ! start the counter for f.c. 1el. state representation here
        j = 0 ! to deal with symetric configuration
        do inc = 1, nam
           n1 = get_na(TargetStates2el%b(nst),inc,1)  ! index to one-electron TargetStates
           if(n1 .ne. nstcore) cycle
           
           n2 = get_na(TargetStates2el%b(nst),inc,2)  
           
!!$           if (myid==0) print*, '--', n1, n2
           
           itmplog = islarger(n2,istcore)
           
           if( n1 .eq. n2) then
              j = n1
              nstA(nst,istcore) = istcore
              CInstA(nst,istcore) = get_CI(TargetStates2el%b(nst),inc)
!!$              write(*,'("Sym core build nst, nst1el, n1, n2:       ",4I4)') nst, nst1el, n1,n2
              cycle
           elseif(itmplog) then 
!!$ Mark: Addition solving non-uniqueness
!!$ All core orbitals (type 2 and 1) should not be included inside the
!!$ rearranged frozen core orbital.
!!$ Keeping these core orbitals stored in the original form
!!$ Note: Symmetric configurations (n1=n2) won't be here
              if ( data_in%non_uniq) then
                 if (any(n2 == ind_coreorb(:)) ) then
!!$ nst_coreorb_cfig keeps track of two-electron states which
!!$ we want to have the original configurations for core orbitals n1 and n2 
                    nst_coreorb_cfig(nst,istcore) = nst_coreorb_cfig(nst,istcore) + 1
!              write(*,'("MCore build nst,build nst, nst1el, n1, n2:",4I4,F10.5)') nst, nst1el, n1,n2,get_CI(TargetStates2el%b(nst),inc)
                    cycle
                 end if
              end if ! Non-uniqueness being solved
!!$ End addition
              ! contniue to build new f.c. state
           else
              cycle ! to avoid double counting of configurations 
           endif
           
           ! f.c. orbital representation is done below
           ! only states n2 with the same value of magnetic quantum number  and parity will be here           
           ! need only to count the first occurance to add a new f.c. 1el state

           !!$ Mark Addition. Some FCO have zero two-electron CI function. Cycle over these
           if ( get_CI(TargetStates2el%b(nst),inc) == 0d0 ) cycle     

           if(i .eq. 0 ) then
              nst1el = nst1el + 1  
              nstton2(nst,istcore) = nst1el
              manst1el(nst1el) = get_ang_mom_proj(TargetStates1el%b(n2))
              iparnst1el(nst1el) =   get_par_st(TargetStates1el%b(n2))
           endif

           i = i + 1
!!$ igiven: indexes of n2 for FCO nst1el
           igiven(nst1el,i) = n2
!!$              write(*,'("FCO build nst, nst1el, n1, n2:            ",4I4)') nst, nst1el, istcore, n2
!!$ CIigiven:  CI coefficients of |n1 n2> that make up the FCO nst1el
           CIigiven(nst1el,i) = get_CI(TargetStates2el%b(nst),inc)
!!$           if (myid==0) print'(10i5)', nst, n1,nst1el,n2, i
!!$ igiven_max: Number of n2 molecular orbitals which make up the FCO nst1el           
           igiven_max(nst1el) = i
           
        enddo ! inc -- configurations making up this two-electron state
        

!!$ Mark: If state has configs 1s1s, 1s2s, .., 1snl, where nl not a core orbital
!!$ The symmetric config gets sucked into the FCO when more than just
!!$ a symmetric config exists. 
!!$ This needs to be avoided when solving non-uniqueness)
        ! Add symmetric config to the list iff there are nonsymmetric already.
        if (.NOT. data_in%non_uniq) then ! Mark: Addition
           if( i .gt. 0 .and. j .ne. 0) then
!!$ if j = 0 then there was no symmetric configuration
!!$ if i = 0 then there was no f.c. orbital for this core orbital 
              i = i + 1
              igiven_max(nst1el) = i ! Number of orbitals to comprise the FC orbital.
              igiven(nst1el,i) = j
              CIigiven(nst1el,i) = CInstA(nst,istcore)/2d0   ! /sqrt(2(1+delta))
              nstA(nst,istcore) = 0  ! set to zero to avoid including it in two places
!!$           if (myid==0) print*,'adding symmetric configuration:', j, j
!!$           if (myid==0) print*, 'igiven_max(nst1el)=', igiven_max(nst1el)
           end if
        end if ! End addition


     enddo ! istcore -- all core orbitals

  enddo ! nst -- two-electron state
!!$  if (myid==0) print*

  call new_basis_st(tmp1elst,nst1el,.true.,0)

!!$  if (myid==0) print*, 'nicm = ', nicm
  do i=1,nicm  ! go through the core one-electron states
     
     nstcore = TargetStates2el%ncore(i)
     call copy_st(tmp1elst%b(i),TargetStates1el%b(nstcore))
     tmp1elst%Nstates = i
     
  enddo
  
  
  nspbst = basis_size(bst_nr)
!!$  if (myid==0) print*, 'nspbst=',nspbst
  allocate(nobst(nspbst),mobst(nspbst),CInobst(nspbst))
  nobst(:) = 0 
  mobst(:) = 0
  CInobst(:) = 0d0
  
  
!!$ Finish dealing with frozen-core representation via target states of H2+
!!$>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  

 
  do i=nicm+1,nst1el
     
!!$  need to transform from representation of frozen-core orbitals in terms of one-electron 
!!$  target states TargetStates1el to reprresentation via underlying one-electron functions with 
!!$  fixed orbital angular momentum given by  bst_nr
     nd = igiven_max(i)   ! Number of orbitals to comprise the FC orbital.
     allocate( iTmp(1:nd), CITmp(1:nd) )
     iTmp(:) = igiven(i,1:nd) ; CITmp(:) = CIigiven(i,1:nd)
    
     !TODO: for some F5 files this will not populate na array properly (try 1*1Sg, 2*2Su)
     call make_nomo(Targetstates1el,bst_nr,nd,CITmp,iTmp,manst1el(i),nspbst,nobst,mobst,CInobst,itmpmax)
     deallocate( iTmp, CITmp )
!!$--------
     ipar = iparnst1el(i)
     rma = manst1el(i)
     energy = 0d0
     j = 0
     nd = itmpmax
     tmp1elst%Nstates = i
     call  construct_st(tmp1elst%b(i),.true.,rma,ipar,0.5d0,energy,i,nd,CInobst(1:nd),nobst(1:nd),mobst(1:nd))
  enddo
  
!!$  if (myid==0) print*, 'size of tmp1elst:',  basis_size_st(tmp1elst)
  call rearrange(bst_nr,get_max_L(bst_nr),tmp1elst,.FALSE.)
 
  do nst=1, Nmax
     
!!$ determine how many antisymmetric configurations are used for this state
!!$ it is equal to the number of core orbital used in description of this state
     is = nint(TargetStates2el%b(nst)%spin)
     i = 0
     phase(:) = (-1)**(is)
     do istcore=1,nicm
        
        if(nstA(nst,istcore) .ne. 0) then   ! symmetric configuration
!!$           print*,"Symm config nst,istcore",nst,istcore
           i = i + 1 
           no1(i) = istcore 
           no2(i) =  istcore 
           nstcore = TargetStates2el%ncore(istcore)
           mo1(i) = get_ang_mom_proj(TargetStates1el%b(nstcore))
           mo2(i) = mo1(i)
           CI(i) = CInstA(nst,istcore)           
        endif
        
!!$ MARK MODIFICATION
!!$ Solving non-uniqueness
!!$ All core orbitals (type 2 and 1) should not be included inside the
!!$ rearranged frozen core orbital.
!!$ Frozen core orbitals not created for n2 core orbital
        if ( data_in%non_uniq) then
           if ( nst_coreorb_cfig(nst,istcore) /= 0  ) then 
!!$ nst_coreorb_cfig(nst,istcore) /= 0
!!$ then this state 'nst' with "inner" core orbital 'istcore' has configuration 
!!$ with an "outer" core orbital
              nstcore = TargetStates2el%ncore(istcore) ! index to core orbital

!!$ Find these configurations with two core orbitals
              do inc = 1, get_nam(TargetStates2el%b(nst)) 
                 n1 = get_na(TargetStates2el%b(nst),inc,1)  ! index to one-electron 'TargetStates1el' 
!!$ make sure this two-electron
!!$ configuration has the same "inner" core orbital as the two-electron
!!$ configuration that has two core orbitals
                 if (n1 .ne. nstcore) cycle 
                 n2 = get_na(TargetStates2el%b(nst),inc,2)
                 itmplog = islarger(n2,istcore)
                 
                 if ( n1 .eq. n2) cycle ! symmetric configuration. Included above
                 if (itmplog .AND. any(n2 == ind_coreorb(:)) ) then
                    ! continue                        
                 else 
                    cycle ! Wrong configuration, should be included in FCO
                 end if
                 
                 i = i + 1
!!$ Core-orbital one-electron target states were copied to
!!$ 'tmp1elst%b' in pisitions 1 to nicm. Accounting for the new indexes.
                 no1(i) = istcore
                 do j = 1, nicm
                    if (n2 == TargetStates2el%ncore(j)) exit
                 end do
                 no2(i) = j
                 mo1(i) = get_ang_mom_proj(TargetStates1el%b(nstcore))
                 ! M = m_1 + m_2, m_2 = M - m_1
                 mo2(i) = get_ang_mom_proj(TargetStates2el%b(nst)) - mo1(i) 
                 if(no1(i) .eq. no2(i)) then   ! This condition is never met.
                    CI(i) = CIigiven(istcore,istcore)
                 else
                ! CI(i) = sqrt(2d0)  ! as per  setupCI() 1d0
                    CI(i) = sqrt(2d0) * get_CI(TargetStates2el%b(nst),inc) 
                 end if ! config type CI
              end do ! inc -two-electron configs for state
           end if ! Config with 2 core orbitals
        end if ! Non-uniqueness being solved. end addition
        
        if(nstton2(nst,istcore) .ne. 0) then  ! f.c. configuration
           i = i + 1 
           no1(i) = istcore 
           nst1el = nstton2(nst,istcore)
           no2(i) = nst1el
           nstcore = TargetStates2el%ncore(istcore)
           mo1(i) = get_ang_mom_proj(TargetStates1el%b(nstcore))
           mo2(i) = manst1el(nst1el)
           if(no1(i) .eq. no2(i)) then   ! This condition is never met  
              CI(i) = CIigiven(istcore,istcore)
           else
              CI(i) = sqrt(2d0)  ! as per  setupCI() 1d0
           endif
        endif
     enddo
     ncm = i

     call setupCI(TargetStates2el%b(nst),ncm,CI,no1,mo1,no2,mo2,phase)        
     
  enddo ! nst - 2el target states

  call destruct_basis_st(TargetStates1el)
  call new_basis_st(TargetStates1el,tmp1elst%Nstates,.true.,0)
  if (myid==0) print*, 'tmp1elst%Nstates:', tmp1elst%Nstates
  TargetStates1el%Nmax = tmp1elst%Nmax
  TargetStates1el%Nstates = tmp1elst%Nstates
  do i=1,tmp1elst%Nstates
     call copy_st(TargetStates1el%b(i), tmp1elst%b(i))

  enddo

!!$ Note  arrays ovlpst(:,:) and e1me(:,:)  might need redefinition (change of diminsions) as number of one-elctron target states is changed and can be larger than previously defined.
!!$! 1el
!!$!  if (myid==0) print*, '1el'
!!$  do nst=1,TargetStates1el%Nmax
!!$     do nstp=1,nst
!!$        tmp = ovlp_st(TargetStates1el%b(nst),TargetStates1el%b(nstp))
!!$        ovlpst(nst,nstp) = tmp
!!$        ovlpst(nstp,nst) = tmp
!!$!        if (myid==0) print*, '3. ', nst, nstp, tmp
!!$     enddo
!!$  enddo

!  call print_orbitals !Liam added this to test rearrange

  return
end subroutine rearrange12_st

subroutine print_orbitals
  use sturmian_class
  use state_class
  use one_electron_func
  use target_states
  implicit none

  integer :: st, orba, orbb, sturm
  type(state), pointer :: pstate, porba, porbb
  type(sturmian_nr), pointer :: psturm
  integer :: l, m, n, nicm, config, a, b
  character(len=:), allocatable :: statelabel, orblabel
  real*8 :: CI
  logical, dimension(TargetStates%Nmax) :: core
  character*13 :: corelabel

  print*, ''; print*, ''
  print*, '<<<<<< PRINTING ORBITALS >>>>>>>>>'
  print*, ''

  nicm = TargetStates2el%nicm
  do orba=1, TargetStates%Nmax
    core(orba) = (orba <= nicm)
  enddo
  
  do st=1, TargetStates2el%Nmax
    pstate => TargetStates2el%b(st)
    statelabel = trim(adjustl(pstate%label))
    write(*,'(A)') '-----------| TWO-ELECTRON STATE: '//statelabel//' |-----------'
    write(*,'(A)') ' Config   na   nb  type(core/fc)'
!                   '   100   100  100  (CORE)*( FC )
    do config=1, pstate%nam
      a = pstate%na(config)
      b = pstate%nb(config)
      CI = pstate%CI(config)
      if(core(a) .and. core(b)) then
        corelabel = '(CORE)*(CORE)'
      elseif(core(a) .and. .not.core(b)) then
        corelabel = '(CORE)*( FC )'
      elseif(.not.core(a) .and. core(b)) then
        corelabel = '( FC )*(CORE)'
      endif
      write(*,'(3X,I3.1,3X,I3.1,2X,I3.1,2X,A13)') config, a, b, corelabel
    enddo !config

  enddo !state


end subroutine print_orbitals











subroutine rearrange12_new(orbitals_st, states_2e)
  !
  ! Rearranges the two-electron target states from (antisymmetrised) CI form
  !
  !   Phi_n(R1,R2) = SUM_ab C_ab^(n) phi_a(R1) phi_b(R2)
  !
  ! to one with Core Orbitals (CO) and Frozen Core (FC) orbitals
  !
  !   Phi_n(R1,R2) = SUM_a phi_a^(CO)(R1) phi_a^(n)(R2)
  !
  ! where the frozen core orbital is
  !
  !   phi_a^(n)(R2) = SUM_a'b delta_a'a C_ab^(n) phi_b(R2)
  !                 = SUM_l C_l^(n) X_l^(n)(r2) Y_l^m(theta2,phi2)
  !
  ! which is itself rearranged to a familiar one-electron state form.
  !
  !
  ! Coded June 2015 by JS.
  !
  use one_electron_func
  use state_class
  implicit none

  type(basis_state), intent(inout) :: orbitals_st, states_2e

  type(basis_state) :: orbitals_tmp
  type(state), pointer :: state_2e
  integer :: j2e,n2e, jCO,nCO, nFC
  integer, dimension(states_2e%nicm) :: vCO

  interface
    subroutine frozen_core_orbitals( first, state_2e, nCO,vCO, nFC, orbitals_tmp )
      use state_class
      logical :: first
      type(state), intent(inout) :: state_2e
      integer, intent(in) :: nCO
      integer, dimension(nCO), intent(in) :: vCO
      integer, intent(inout) :: nFC
      type(basis_state), intent(inout), optional :: orbitals_tmp
    end subroutine
    subroutine rearrange(bst_nr,LMAX,       TargetStates,use_lag_ham1el,lag_ham1el_m)
      use input_data
      use grid_radial
      use sturmian_class
      use state_class
      use MPI_module
      type(basis_sturmian_nr), intent(inout) :: bst_nr 
      integer, intent(in):: Lmax   ! max L valus in one-electron basis bst
      type(basis_state), intent(inout):: TargetStates
      logical, intent(in):: use_lag_ham1el ! switch use lag_ham1el_m - Laguerre basis <Lp|H|L> 
      real*8, dimension(bst_nr%n,bst_nr%n,-Lmax:Lmax), intent(in):: lag_ham1el_m  ! 
      optional:: lag_ham1el_m
    end subroutine
  end interface
  
  n2e = basis_size(states_2e)
  nCO = states_2e%nicm; vCO(:) = states_2e%ncore(:)

  ! The first step is to determine the number of orbitals (states) we need.
  nFC = 0   ! Tally of frozen core orbitals.
  do j2e = 1, n2e
     state_2e => states_2e%b(j2e)
     call frozen_core_orbitals(.true., state_2e, nCO,vCO, nFC)
  end do
  write(*,'(/A,I4,A,I4,A/)') 'rearrange12: ', nCO, ' core orbitals and ', nFC, ' frozen core orbitals.'

! Now create the (state) basis to hold all the (CO+FC) orbitals.
  call new_basis( orbitals_tmp, nCO+nFC, .true., states_2e%basis_type )


  ! The second step is to fill the orbital basis, starting with core orbitals.
  do jCO = 1, nCO
     call copy( orbitals_tmp%b(jCO), orbitals_st%b(vCO(jCO)) )
     states_2e%ncore(jCO) = jCO   ! Core orbitals now occupy first nCO spaces.
  end do

  ! And now make the frozen core orbitals.
  nFC = 0
  do j2e = 1, n2e
     state_2e => states_2e%b(j2e)
     call frozen_core_orbitals( .false., state_2e, nCO,vCO, nFC, orbitals_tmp )
  end do


  ! Finally, we re-make the orbital basis with the CO+FC states.
  call destruct_basis_st(orbitals_st)
  call new_basis( orbitals_st, nFC+nCO, .true., orbitals_tmp%basis_type )
  do jCO = 1, nFC+nCO
     call copy( orbitals_st%b(jCO), orbitals_tmp%b(jCO) )
  end do
  call rearrange( bst_nr, get_max_L(bst_nr), orbitals_st, .FALSE. )


end subroutine rearrange12_new


subroutine frozen_core_orbitals( first, state_2e, nCO,vCO, nFC, orbitals_tmp )
  !
  use state_class
  implicit none

  logical :: first   ! True for first run through -- only calculate array sizes.
  type(state), intent(inout) :: state_2e   ! The 2e state to rearrange.
  integer, intent(in) :: nCO   ! The number of core orbitals.
  integer, dimension(nCO), intent(in) :: vCO   ! Indices of the core orbitals.
  integer, intent(inout) :: nFC
  type(basis_state), intent(inout), optional :: orbitals_tmp

  type(state) :: orbital_FC
  logical :: nonsymmetric
  integer :: jCO, jConf,nConf, nOrb,nConfNew, indCO,ind1,ind2
  real*8 :: mCO
  integer, dimension(:), allocatable :: vOrb
  real*8, dimension(:), allocatable :: CIOrb
  integer, dimension(nCO) :: ind1Vec,ind2Vec, m1Vec,m2Vec, phaseVec
  real*8, dimension(nCO) :: CIVec

  nConf = get_nam(state_2e)   ! Number of configurations that build this state.
  nConfNew = 0   ! Number of configurations to be used after the rearrangement.
  allocate( vOrb(nConf), CIOrb(nConf) )   ! Max configurations it could be.

  do jCO = 1, nCO   ! Loop through the core orbitals of all 2e states.
     indCO = vCO(jCO)   ! Index of this core orbital: a.

     nOrb = 0   ! The number of orbitals that will make up the FCO.
     vOrb(:) = 0   ! Indices of the orbitals that will make the FC orbital.

     ! We will only make a FCO if there is at least one nonsymmetric config.
     ! If there is only one, symmetric config |aa>, we link both to the same CO.
     nonsymmetric = .false.

     do jConf = 1, nConf   ! Loop through the configurations of this 2e state.
        ind1 = get_na(state_2e,jConf,1)   ! Index of the r1 orbital: a'
        ind2 = get_na(state_2e,jConf,2)   ! Index of the r2 orbital: b
        if (ind1 /= indCO) cycle   ! delta_a'a

        ! The configurations have been antisymmetrised by construct_st/setupCI,
        !   but we want to ignore |ba> for now and antisymmetrise again later.
        ! We know it is |ba> if the r2 index comes before the r1 index,
        !   and if the second index is one of the core orbital indices.
        if ( ind1>ind2 .and. any(ind2==vCO(:)) ) cycle
        if (ind1 /= ind2) nonsymmetric = .true.   ! Need to make a FCO.

        ! Add to the list of orbitals to make a FCO (unless only symmetric).
        nOrb = nOrb + 1
        if (.not.first) then
           vOrb(nOrb) = ind2
           if (ind1 == ind2) then
              CIOrb(nOrb) = get_CI(state_2e,jConf) / 2d0
           else
              CIOrb(nOrb) = get_CI(state_2e,jConf)
           end if
        end if

     end do

     if (nOrb == 0) cycle
     if (.not.first) mCO = get_ang_mom_proj( orbitals_tmp%b(jCO) )
     nConfNew = nConfNew + 1

     if (nonsymmetric) then   ! At least one nonsymmetric configuration.
        nFC = nFC + 1   ! We need to make a FCO for this state & CO.

        if (.not.first) then   ! Combine all r2 orbitals into one FCO.
           call make_FCO( nOrb,vOrb(1:nOrb),CIOrb(1:nOrb), orbital_FC )
           orbitals_tmp%b(nCO+nFC) = orbital_FC

           CIVec(nConfNew) = dsqrt(2d0)   ! Cancels with 1/sqrt(2) in setupCI.
           ind1Vec(nConfNew) = jCO; ind2Vec(nConfNew) = nCO+nFC
           m1Vec(nConfNew) = mCO; m2Vec(nConfNew) = get_ang_mom_proj(orbital_FC)
        end if

     elseif (.not.first) then   ! Symmetric-only configuration.
        CIVec(nConfNew) = CIOrb(nOrb) * 2d0   ! To cancel out 1/2 from above.
        ind1Vec(nConfNew) = jCO; ind2Vec(nConfNew) = jCO
        m1Vec(nConfNew) = mCO; m2Vec(nConfNew) = mCO
     end if

  end do


  if (.not.first) then
     phaseVec(1:nConfNew) = (-1)**get_spin(state_2e)
     call setupCI( state_2e, nConfNew,CIVec(1:nConfNew),ind1Vec(1:nConfNew),m1Vec(1:nConfNew),ind2Vec(1:nConfNew),m2Vec(1:nConfNew),phaseVec(1:nConfNew) )
  end if


end subroutine frozen_core_orbitals


subroutine make_FCO( nOrb,vOrb,CIOrb, orbital_FC )
  !
  ! Combines a list of orbital state indices (vOrb) and CI coefficients (CIOrb)
  !   into a frozen core orbital thru the underlying single particle functions.
  !
  ! Essentially a re-coded make_nomo().
  !
  ! NOTE that the orbitals states are from TargetStates (module target_states),
  !   while the s.p. functions are from bst_nr (module one_electron_func).
  !
  ! Coded June 2015 by JS.
  !
  use grid_radial
  use input_data
  use one_electron_func
  use target_states
  implicit none

  integer, intent(in) :: nOrb
  integer, dimension(nOrb), intent(in) :: vOrb
  real*8, dimension(nOrb), intent(in) :: CIOrb
  type(state), intent(out) :: orbital_FC

  type(state), pointer :: orbital_st
  integer :: nBst, jOrb,ind2, jSP,indSP,nSP, par
  integer, dimension(:), allocatable :: vSP, mVec
  real*8 :: CIConf,CISP, m
  real*8, dimension(:), allocatable :: CIVec


  nBst = basis_size(bst_nr)
  allocate( vSP(nBst), CIVec(nBst) )
  vSP(:) = 0; CIVec(:) = 0d0

  ! For each orbital, find its underlying s.p. functions contribution to the FC.
  do jOrb = 1, nOrb
     ind2 = vOrb(jOrb)   ! Index of the orbital state.
     CIConf = CIOrb(jOrb)   ! CI coefficient of the two-electron configuration.

     orbital_st => TargetStates%b(ind2)
     do jSP = 1, get_nam(orbital_st)
        indSP = get_na(orbital_st,jSP)   ! S.p. function index within the basis.
        CISP = get_CI(orbital_st,jSP)   ! CI coefficient of this s.p. function.

        ! The cumulative CI of every s.p. function in the basis for the FC.
        CIVec(indSP) = CIVec(indSP) + CIConf*CISP
     end do
  end do

  ! Make a shorter list of s.p. functions, with only the ones that contribute.
  nSP = 0   ! Tally of s.p. functions that will contribute to the FC.
  do indSP = 1, nBst   ! Loop over the whole s.p. basis.
     if (CIVec(indSP) /= 0d0) then
        nSP = nSP + 1
        vSP(nSP) = indSP   ! Vector of contributing s.p. function indices.
        CIVec(nSP) = CIVec(indSP)   ! Collapse, because nSP <= indSP always.
     end if
  end do


  ! Gather all the pieces and build the FC state for this target state and CO.
  ! All contributing orbital states will have the same m, parity, and spin.
  orbital_st => TargetStates%b(vOrb(1))
  m = get_ang_mom_proj(orbital_st)
  par = get_par(orbital_st)
  allocate( mVec(nSP) ); mVec(:) = m

  call construct_st( orbital_FC,.true., m,par,0.5d0,0d0,0, nSP,CIVec(1:nSP),vSP(1:nSP),mVec(1:nSP) )


end subroutine make_FCO


!!$
!!$ Check if this 1el state n is a core 1el state
function iscore(n)

  use state_class
  use target_states

  implicit none

  integer, intent(in):: n
  logical:: iscore
  integer:: ic, ncst

  iscore = .false.
  do ic=1,TargetStates2el%nicm
     ncst = TargetStates2el%ncore(ic)
     if(n .eq. ncst) then
        iscore = .true.
        exit
     endif
  enddo

  return
end function iscore
!!$
!!$ Check if this 1el state n is a core 1el state and if its  number is larger that core state number k
function islarger(n,k)

  use state_class
  use target_states
  use MPI_module

  implicit none
  
  logical:: islarger
  integer, intent(in):: n,k
  logical  iscore
  integer:: ic, ncst

  iscore = .false.
  islarger = .false.
!!$ iscore=.true. if n = index to a core orbital.
  do ic=1,TargetStates2el%nicm        ! go through the core one-electron (molecular) orbitals
     ncst = TargetStates2el%ncore(ic) ! index to core (Laguerre or molecular) orbital/1el state
     if(n .eq. ncst) then
        iscore = .true.                 
        exit
     endif
  enddo

!!$ islarger = .true. if n /= index of a core orbital
!!$ OR n = index to a core orbital .AND. ic  > k
  if(iscore) then 
     if(ic .gt. k) then
        islarger = .true.
!        if (myid==0) print*, '!!', n, ic,k
     endif
  else
     islarger = .true.
!     if (myid==0) print*, '>>', n, ic,k
  endif

  return
end function islarger

!!$-----------------------------------------------------------------------------------------------------------------
!!$  transform from repreesentation of frozen-core orbitals in terms of one-electron 
!!$  target states TargetStates1el to reprresentation via underlying one-electron functions with 
!!$  fixed orbital ongular momentum given by  bst
subroutine  make_nomo(TargetStates1el,bst,nd,CIigiven,igiven,manst1el,nspbst,nobst,mobst,CInobst,itmpmax)

  use sturmian_class
  use state_class
  use ovlpste1me
  use MPI_module

  implicit none
  
  type(basis_state), intent(in):: TargetStates1el
  type(basis_sturmian_nr), intent(in):: bst   ! has fixed  angular momentum 
  integer, intent(in):: nd
  real*8, dimension(nd), intent(in):: CIigiven
  integer, dimension(nd), intent(in):: igiven
  integer, intent(in)::   manst1el, nspbst
  integer, dimension(nspbst), intent(out):: nobst, mobst
  real*8, dimension(nspbst), intent(out):: CInobst
  integer, intent(out):: itmpmax

  integer:: i, n,nst, mst, nam, j, np, nstp
  real*8:: CI1el, CI2el, tmp, CI2elp
  real*8, dimension(nspbst):: CItmp

  nobst(:) = 0 
  mobst(:) = 0
  CInobst(:) = 0d0
  itmpmax = 0

  CItmp(:) = 0d0

  mobst(:) = manst1el  !! they have all the same value m

  mst = manst1el
  
  do n=1,nd

     nst = igiven(n)  ! index to one-electron target states  TargetStates1el
!     if (myid==0) print*, 'n,nst=', n,nst
     CI2el = CIigiven(n)  ! this is coef. with which those one-electron states contribute this the current FC orbital
     
     nam = get_nam(TargetStates1el%b(nst))  ! representation of target state nst  via one-electron  bst  basis
     do i=1,nam
        
        j = get_na(TargetStates1el%b(nst),i)   ! index to one-electron basis function
        CI1el = get_CI(TargetStates1el%b(nst),i)

        CInobst(j) = CInobst(j) + CI1el*CI2el
!        if (myid==0) write(*,'(4i5,2E15.5)') n,nst,i,j, CI1el,CI2el
     enddo

  enddo


  i = 0
  do j=1,nspbst

     if(CInobst(j) .eq. 0 ) cycle
     
     i = i + 1
     
     CItmp(i) = CInobst(j)
     nobst(i) = j
     
!     if (myid==0) print*,'make_nomo: ***', j, get_ang_mom(bst%b(j)), CItmp(i)

  enddo
  
!  if (myid==0) print*, '>>>', i, nspbst

  CInobst(:) = 0d0
  CInobst(1:i) = CItmp(1:i)

  itmpmax = i

  
  
!!$!test
!!$
!!$  tmp = 0d0
!!$
!!$ do n=1,nd
!!$
!!$     nst = igiven(n)  ! index to one-electron target states  TargetStates1el
!!$     CI2el = CIigiven(n)  ! this is coef. with which those one-electron states contribute this the current FC orbital
!!$
!!$     
!!$     do np=1,nd
!!$        
!!$        nstp = igiven(np)  ! index to one-electron target states  TargetStates1el
!!$        CI2elp = CIigiven(np)  ! this is coef. with which those one-electron states contribute this the current FC orbital
!!$        
!!$        tmp = tmp + CI2el*CI2elp*ovlpst(nst,nstp)
!!$
!!$     enddo
!!$
!!$  enddo
!!$
!!$  print*, 'tmp=',tmp
!!$
!!$
!!$
!!$
!!$  tmp = 0d0
!!$
!!$ do n=1, itmpmax
!!$
!!$     nst = nobst(n)
!!$     CI2el = CInobst(n)  ! this is coef. with which those one-electron states contribute this the current FC orbital
!!$
!!$     
!!$     do np=1, itmpmax
!!$        
!!$        nstp = nobst(np)
!!$        CI2elp = CInobst(np)  ! this is coef. with which those one-electron states contribute this the current FC orbital
!!$        
!!$        tmp = tmp + CI2el*CI2elp * bst%ortint(nst,nstp)
!!$
!!$     enddo
!!$
!!$  enddo
!!$
!!$  print*, 'tmp=',tmp

end subroutine make_nomo
