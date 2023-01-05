!!$------------------------------------------------------------------------

!To use natural orbitals:
!Create a file 'natorb' with a number in it specifiying the number of natural orbitals to include.
! Three modes:
! - Positive number N: Hybrid basis formed by inserting N natural orbitals into the original two-electron basis
! - N = 0: Basis formed from all natural orbitals with no Laguerre functions
! - Negative number N: basis formed from |N| number of natural orbitals with no Laguerre functions

module natural_orbitals_module

  real*8 :: fullCI_energy

  contains
subroutine natorbitals(TargetGS) ! TargetStates2el%b(1)

  use sturmian_class
  use state_class
  use one_electron_func
  use target_states
  use input_data
  use ovlpste1me

  implicit none  

! passed via "use target_states":  type(basis_state), intent(in):: TargetStates ! - basis of one-electron target states 
  type(state), intent(in):: TargetGS  ! ground state of H2
  type(state) :: testGS, temp2el, temp1el
  type(state), pointer :: orb, orbA, orbB, orb_prev, orb_next
  real*8, dimension(:,:), allocatable :: CI
  integer:: nusemax_full !  number of orbitals used to describe this state (useful for two-electron states only)
  integer :: nusemax ! number m>=0 orbitals
  integer, dimension(:), allocatable:: nuse, nuse_full   ! idex to orbitals used to describe the state
  integer:: nam  ! size of array na(:)
  integer, dimension(:), allocatable:: nuse_inv, nuse_negm
  integer:: n, nc, i, ii, j, jj, ni, nj, norbmax, n_prev, n_next

  real*8, dimension(:,:), allocatable :: z, bmat
  real*8, dimension(:), allocatable:: w, ww
  integer:: matz, ierr, itmp
  real*8:: tmp
  integer, dimension(:), allocatable:: Num

  logical :: eigen_sort

  type(basis_state) :: natorbs
  integer :: num_natorb, num_natorb_incl
  real*8 :: m, wtemp
  integer :: par, nci, ma, pa, mb, pb, m_prev, m_next, m_need
  integer, dimension(:,:), allocatable :: inum
  integer, dimension(:), allocatable :: no1, mo1, no2, mo2, phase
  real*8, dimension(:), allocatable :: civec
  character*1 :: pmchar1, pmchar2
  character*5 :: labelA, labelB, label, label_prev, label_next

  type(basis_sturmian_nr) :: bst_data, bst_MSC, bst_sub, bst_temp
  integer :: nBasis_MSC, nBasis_data, max_latop, lorb, nst_basis_MSC, max_m, im
  real*8, dimension(:,:,:), allocatable:: lag_ham1el_m
  
  character*20, parameter :: symcharsUP = 'u.gSPDFGHIKLMNOQRTUV'
  character*20, parameter :: symcharsLOW = 'u.gspdfghiklmnoqrtuv'
  
  logical :: spheroidal
  
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
  
  !type(basis_sturmian) :: TargetStates2el_fullCI, TargetStates_temp

  fullCI_energy = TargetGS%energy
  !call copy_basis_st(TargetStates2el_fullCI, TargetStates2el)

  spheroidal = data_in%calculation_type == 2 .or. data_in%calculation_type == 3

  !Write labels for the hybrid basis
  do i=1, TargetStates%Nmax
    orb => TargetStates%b(i)
    
    if(orb%energy == 0) then !Laguerre function
      write(orb%label,'(X,I1,A1,A1,X)') orb%n, symcharsLOW(orb%l+4:orb%l+4), symcharsUP(nint(abs(orb%m))+4:nint(abs(orb%m))+4)
    else !MO
      write(orb%label,'(X,I1,A1,A1,A1)') orb%n, symcharsLOW(orb%l+4:orb%l+4), symcharsUP(nint(abs(orb%m))+4:nint(abs(orb%m))+4), &
        &symcharsLOW(orb%parity+2:orb%parity+2)
    endif
  enddo

  write(*,*)
  write(*,'(" =*=*=*=*=*=*=*=*= MAKING NATURAL ORBITALS =*=*=*=*=*=*=*=*=*=*")')
  write(*,*)

  if(get_ang_mom_proj(TargetGS) /= 0) error stop '*** ERROR in natorbitals: M /= 0'

!!$ Fist need to transfom the CI array in TargetGS (one-dimensional, CI(nc) - per configuration)  into a two-demensional array C(i,j) - per one-electron orbitals (one-electron H2+ type target states)  
!!$ NOTE that indexes "i,j" in C(i,j) are not to the onelectron orbitals but to  pointers to them, this is because
!!$ not all orbitals might be used for given diagonalization to get CI(nc)

  !number of one-electron orbitals used for this state
  nusemax_full = TargetGS%nusemax
  
  !number of configurations
  nam = TargetGS%nam

!  print*, ''
!  write(*,'(A)') '========== EXPANSION OF GROUND STATE IN TargetStates BASIS =========='
!  write(*,'("config   orbA   (mA,pA)   orbB   (mB,pB)     CI")')
!  do i=1, nam
!    ni = TargetGS%na(i)
!    nj = TargetGS%nb(i)
!    orbA => TargetStates%b(ni)
!    orbB => TargetStates%b(nj)
!    mA = get_ang_mom_proj(orbA)
!    pA = get_par(orbA)
!    mB = get_ang_mom_proj(orbB)
!    pB = get_par(orbb)
!    labelA = trim(adjustl(get_label(orbA)))
!    labelB = trim(adjustl(get_label(orbB)))
!    write(*,'(I4,I5,"(",A4,") (",I2,",",I2,")",I3,"(",A4,") (",I2,",",I2,")  ",ES12.5)') &
!      &i, ni, labelA, mA, pA, nj, labelB, mB, pB, TargetGS%CI(i)
!  enddo
!  write(*,'(A)') '======================================================='
!  print*, ''
  
  !number of one-electron orbitals not including M<0 orbitals
  nusemax = 0
  do i=1, nusemax_full
    orb => TargetStates%b(i)
    if(get_ang_mom_proj(orb) >= 0) nusemax = nusemax + 1
  enddo

  !Array containing indices of one-electron orbitals used for this state
  allocate(nuse_full(nusemax_full))
  nuse_full(:) = TargetGS%nuse(:)
  
  !Array containing indices of M>=0 one-electron orbitals used for this state
  !write(*,'(A)') '===== M >=0 orbitals used to represent this state ===='
  !write(*,'("index  orb    ( m, p)"   )')
  allocate(nuse(nusemax))
  j=0
  do i=1, nusemax_full
    ni = nuse_full(i)
    if(get_ang_mom_proj(TargetStates%b(ni)) < 0) cycle
    j = j + 1
    nuse(j) = ni
    orb => TargetStates%b(ni)
    !write(*,'(I3,I5,"(",A4,")",2X,"(",I2,",",I2,")")') j, ni, trim(adjustl(get_label(orb))), get_ang_mom_proj(orb), get_par(orb)
  enddo
  !write(*,'(A)') '======================================================'
  !print*, ''
  
  !write(*,'(A)') '===== All orbitals used to represent this state ===='
  !write(*,'("index  orb    ( m, p)"   )')
  j=0
  do i=1, nusemax_full
    ni = nuse_full(i)
    j = j + 1
    orb => TargetStates%b(ni)
    !write(*,'(I3,I5,"(",A4,")",2X,"(",I2,",",I2,")")') j, ni, trim(adjustl(get_label(orb))), get_ang_mom_proj(orb), get_par(orb)
  enddo
  !write(*,'(A)') '======================================================'
  !print*, ''

  allocate(CI(nusemax,nusemax))
  CI(:,:) = 0d0
  allocate(bmat(nusemax,nusemax))
  bmat(:,:) = 0d0
  do i=1,nusemax
     bmat(i,i) = 1d0
  enddo

  allocate(w(nusemax))
  allocate(z(nusemax,nusemax))

  allocate(Num(nusemax))
  
  norbmax = MAXVAL(nuse_full)
  allocate(nuse_inv(norbmax), nuse_negm(norbmax))
  nuse_inv(:) = 0

  n_prev = -1
  m_prev = -1
  label_prev = 'xxx'
  do i=1,nusemax_full
    n = nuse_full(i)
    orb => TargetStates%b(n)
    m = get_ang_mom_proj(orb)
    label = trim(adjustl(get_label(orb)))

    do j=1, nusemax
      orbA => TargetStates%b( nuse(j) )
      if(orbA%n == orb%n .and. orbA%l == orb%l .and. orbA%parity == orb%parity .and. orbA%spin == orb%spin &
        &.and. orbA%inum == orb%inum .and. orbA%m == abs(m)) exit
    enddo !j
    if(j > nusemax) then
      print*, 'i, n, label:', i, n, label
      error stop '*** ERROR in natorbs: nuse_inv'
    endif
    nuse_inv(n) = j
    
    do j=1, nusemax_full
      orbA => TargetStates%b( nuse_full(j) )
      if(orbA%n == orb%n .and. orbA%l == orb%l .and. orbA%parity == orb%parity .and. orbA%spin == orb%spin &
        &.and. orbA%inum == orb%inum .and. orbA%m == -m) exit
    enddo !j
    if(j > nusemax_full) then
      print*, 'i, n, label:', i, n, label
      error stop '*** ERROR in natorbs: nuse_negm'
    endif
    nuse_negm(n) = nuse_full(j)

  enddo

  tmp = 0d0
  do nc=1,nam

     ni = TargetGS%na(nc)
     nj = TargetGS%nb(nc)
     
     i = nuse_inv(ni)
     j = nuse_inv(nj)

     labelA = TargetStates%b(ni)%label
     labelB = TargetStates%b(nj)%label
     
    ! write(*,'("nc, ni, nj, i, j, labeli, labelj, CI:",5i5,2A6,ES15.5)',advance='no') nc, ni, nj, i, j, labelA, labelB, TargetGS%CI(nc)
     if(get_ang_mom_proj(TargetStates%b(ni)) < 0) then
    !   write(*,*) '   (SKIP)'
       cycle
     else
    !   write(*,*) '   (KEEP)'
     endif

     CI(i,j) = TargetGS%CI(nc)

     tmp = tmp + CI(i,j)*CI(i,j)

  enddo
!  print*, 'tmp=', tmp

  print*
  if(nusemax <= 10 .and. .false.) then
    print*,'matrix: note i,j are not indexes to the one-electron orbitals, to get it: n = nuse(i) '
    print'(13X,50("(",A4,")",9X))', (trim(adjustl(TargetStates%b(nuse(i))%label)), i=1, nusemax)
    do i=1,nusemax
      print'("(",A4,")",2X,50ES15.5)',trim(adjustl(TargetStates%b(nuse(i))%label)),(CI(i,j), j=1,nusemax)
    enddo
  endif

  !
  !      do n=1,nspm
  !         write(*,'(50F10.5)') (CI(n,m), m=1,nspm)         
  !      enddo

  ! Assume that s.p.orbitals come from H2+ diagonalization and form an orthonomal basis
  ! CI array is for H2 ground state: it is a symmetric array and can be diagonalized.
  ! the dimension of the CI array: nusemax x nusemax
  
  matz=2

  call  rsg(nusemax,nusemax,CI,bmat,w,matz,z,ierr)
  
  if(ierr /= 0) then
    print*, 'natorb: ierr =', ierr
    error stop
  endif

  eigen_sort = .true.
! order natural orbitals by the magnitude of their eigenvalues
! This array will be sorted by insertion sort algorithm.
! sorting by nat.orb. weights.
  if(eigen_sort) then
    do n=1,nusemax
       Num(n) = n
    end do
    do n=2,nusemax
       itmp = Num(n)
       j = n
       do while(j.ge.2)
          if(abs(w(itmp)).gt.abs(w(Num(j-1)))) then
             Num(j) = Num(j-1)
             j = j - 1
          else
             exit
          endif
       end do     
       Num(j) = itmp
    end do
 
    w(1:nusemax) = w(Num(1:nusemax))
    z(1:nusemax,1:nusemax) = z(1:nusemax,Num(1:nusemax))
  endif !eigen_sort

  ! print for a natural orbital its representation through original orbitals

!!!Create basis to store the natural orbitals
  !count number of natural orbitals (non-zero eigenvalues)
  num_natorb=0
  do jj=1, nusemax
    if(w(jj) == 0.0d0) cycle
    num_natorb = num_natorb + 1
    !check all orbitals making up this natural orbital have the same m and parity
    i=0
    do ii=1, nusemax
      if(z(ii,jj) == 0.0d0) cycle
      i=i+1
      n = nuse(ii)
      if(i==1) then
        m = get_ang_mom_proj( TargetStates%b(n) )
        par = get_par( TargetStates%b(n) )
      endif

      if(get_ang_mom_proj( TargetStates%b(n) ) /= m .or. m < 0) error stop '*** NATORB: M ERROR'
      if(data_in%good_parity .and. get_par( TargetStates%b(n) ) /= par ) then
        print*, get_par( TargetStates%b(n) ), par 
        error stop '*** NATORB: PARITY ERROR'
      endif

    enddo !ii
    
    !Make an additional natural orbital with -m if m > 0
    if(m > 0) num_natorb = num_natorb + 1

  enddo !jj

  print*, 'Creating basis for', num_natorb, 'natural orbitals'
  
  call new_basis(natorbs,num_natorb,.true.,data_in%calculation_type)
  allocate(ww(num_natorb))
  max_m = 0
  !Now go through the same loops and make the natural orbitals
  j=0
  do jj=1, nusemax
    if(w(jj) == 0.0d0) cycle
    j = j + 1
    ww(j) = w(jj)
    nci = count(z(:,jj) /= 0.0d0)
    allocate(no1(nci), mo1(nci), civec(nci))
    i=0
    do ii=1, nusemax
      if(z(ii,jj) == 0.0d0) cycle
      i = i + 1
      no1(i) = nuse(ii)
      mo1(i) = get_ang_mom_proj(TargetStates%b(no1(i)))
      civec(i) = z(ii,jj)
    enddo !ii
    m = get_ang_mom_proj(TargetStates%b(no1(1)))
    par = get_par(TargetStates%b(no1(1)))
    
    call construct_st_(natorbs%b(j), .true., dble(m), par, 0.5d0, 0.0d0, 1, nci, civec, no1, mo1)
    
    ii = maxloc(abs(civec),1)
    ni = no1(ii)
    ni = no1(1)
    orbA => TargetStates%b(ni)
    natorbs%b(j)%n = orbA%l+1 !orbA%n
    natorbs%b(j)%l = orbA%l
    
    max_m = max(nint(abs(m)),max_m)

    if(m>0) then
      j = j + 1
      ww(j) = w(jj)
      m = -m
      mo1 = -mo1
      do ii=1, nci
        if(get_ang_mom_proj( TargetStates%b(no1(ii)) ) > 0) no1(ii) = nuse_negm(no1(ii))
      enddo
      call construct_st_(natorbs%b(j), .true., dble(m), par, 0.5d0, 0.0d0, 1, nci, civec, no1, mo1)
      natorbs%b(j)%n = orbA%n
      natorbs%b(j)%l = orbA%l
    endif

   
    deallocate(no1,mo1,civec)

  enddo !jj
  
  do j=1, num_natorb
    orb => natorbs%b(j)
    
    m = orb%m
    par = orb%parity

    if(m>=0) then
      write(orb%label,'("N",I1,"-",A1,A1)') j, symcharsUP(abs(nint(m))+4:abs(nint(m))+4), symcharsLOW(par+2:par+2)
    else
      write(orb%label,'("n",I1,"-",A1,A1)') j, symcharsUP(abs(nint(m))+4:abs(nint(m))+4), symcharsLOW(par+2:par+2)
    endif

    orb%energy = 1.0d0
  enddo
  print*, ''
  print*, 'EIGENVALUES SORTED BY MAGNITUDE:'
  print*, ww(1:num_natorb)
  print*, ''

  !  print*, '---------------------------------------'
  !  do j=1, num_natorb
  !    orb => natorbs%b(j)
  !    write(*,'(" Natural orbital #",I0,", eigenvalue = ",ES12.5)') j, ww(j)
  !    write(*,'(" na         ma  CI")')
  !    do i=1, orb%nam
  !      ii = orb%na(i)
  !      m = orb%ma(i)
  !      write(*,'(I3,"(",A4,")",I5,1X,ES12.5)') ii, trim(adjustl(TargetStates%b(ii)%label)), nint(m), orb%CI(i)
  !    enddo
  !    write(*,*)
  !  enddo
  print*, ''
  write(*,'("--NATURAL ORBITALS -------------------------")')
  print*, 'i, eigenvalue,  label, ( M, P,S)'
  do j=1, min(num_natorb,10)
    orb => natorbs%b(j)
    write(*,'(I2,X,ES12.5,2X,A5,2X,"(",I2,",",I2,",",I1,")")') j, ww(j), orb%label, nint(orb%m), orb%parity, nint(orb%spin)
  enddo
  write(*,'("--------------------------------------------")')
  print*, ''
 
  !reconstruct ground state
  do j=1, num_natorb
    orbA => natorbs%b(j)
    m = get_ang_mom_proj(orbA)
    if(m == 0) then
      orbB => natorbs%b(j)
    elseif(m > 0) then
      orbB => natorbs%b(j+1)
    else
      orbB => natorbs%b(j-1)
    endif

    nci = orbA%nam**2
    allocate(no1(nci), mo1(nci), no2(nci), mo2(nci), civec(nci), phase(nci))
    i=0
    do ii=1, orbA%nam
      do jj=1, orbB%nam
        i=i+1
        civec(i) = orbA%ci(ii)*orbB%ci(jj)
        no1(i) = orbA%na(ii)
        no2(i) = orbB%na(jj)
        mo1(i) = orbA%ma(ii)
        mo2(i) = orbB%ma(jj)
      enddo
    enddo
    civec = civec * ww(j)
    phase(:) = (-1)**(targetGS%spin)
   
    call construct_st_(temp2el, .false., targetGS%m, targetGS%parity, targetGS%spin, 0.0d0, 1, nci, civec, no1, mo1, no2, mo2,phase)
    deallocate(no1, mo1, no2, mo2, civec, phase)

    if(j==1) then
      call copy_st(testGS, temp2el)
    else
      call sum_into_2el_b4_rearrange(testGS, temp2el)
    endif

    call destruct_st(temp2el)
    
  enddo
  
!  write(*,'(A)') '========== EXPANSION OF TESTGS IN TargetStates BASIS =========='
!  write(*,'("config   orbA   (mA,pA)   orbB   (mB,pB)     CI")')
!  do i=1, nam
!    ni = TargetGS%na(i)
!    nj = TargetGS%nb(i)
!    orbA => TargetStates%b(ni)
!    orbB => TargetStates%b(nj)
!    mA = get_ang_mom_proj(orbA)
!    pA = get_par(orbA)
!    mB = get_ang_mom_proj(orbB)
!    pB = get_par(orbb)
!    labelA = trim(adjustl(get_label(orbA)))
!    labelB = trim(adjustl(get_label(orbB)))
!    write(*,'(I4,I5,"(",A4,") (",I2,",",I2,")",I3,"(",A4,") (",I2,",",I2,")  ",ES12.5)') &
!      &i, ni, labelA, mA, pA, nj, labelB, mB, pB, TargetGS%CI(i)
!  enddo
!  write(*,'(A)') '======================================================='
!  print*, ''

  deallocate(nuse,nuse_inv,w,z,bmat,Num) 

!Below was to sort the NOs by symmetry but I decided against it
!  !Sort 
!  jj=0
!  do i=0, max_m
!    do par = 1, -1, -2
!      do ii=1, -1, -2
!        if(i==0 .and. ii==-1) cycle
!        im = i*ii
!        jj = jj + 1 !jj = symmetry counter
!        !Find NO in this symmetry with largest eigenvalue
!        wtemp = 0.0d0
!        do j=1, num_natorb
!          orb => natorbs%b(j)
!          if(nint(orb%m) == im .and. orb%parity == par) then
!            wtemp = max(wtemp,abs(ww(j)))
!          endif
!        enddo
!        do j=1, num_natorb
!          orb => natorbs%b(j)
!          if(nint(orb%m) == im .and. orb%parity == par .and. abs(ww(j)) == wtemp) then
!            call copy_st(temp1el,natorbs%b(jj))
!            call copy_st(natorbs%b(jj),orb)
!            call copy_st(orb,temp1el)
!            wtemp = ww(jj)
!            ww(jj) = ww(j)
!            ww(j) = wtemp
!            exit
!          endif
!        enddo
!      enddo
!    enddo
!  enddo

  !Set inum of natural orbitals now that they are in final order
  allocate(inum(-get_max_L(bst_nr):get_max_L(bst_nr),-1:1))
  inum = 0
  do j=1, num_natorb
    orb => natorbs%b(j)
    m = orb%m
    par = orb%parity
    inum(nint(m),par) = inum(nint(m),par) + 1
    orb%inum = inum(nint(m),par)
    orb%n = orb%l + orb%inum
  enddo
  deallocate(inum)
  
  num_natorb_incl = data_in%num_nat_orb
  if(num_natorb_incl == -1) then
    num_natorb_incl = 0
  elseif(data_in%only_nat_orbs) then
    num_natorb_incl = -num_natorb_incl
  endif

!  call new_basis_st(TargetStates_temp, TargetStates%Nmax+abs(num_natorb_incl), .true., TargetStates%basis_type)
!
!  do j=1, TargetStates%Nmax
!    call copy_st(TargetStates_temp%b(j), TargetStates%b(j))
!  enddo
!
!  do j=1, num_natorb_incl
!    call copy_st(TargetStates_temp%b(TargetStates_backup%Nmax+j), natorbs%b(j))
!  enddo
!
!  call copy_basis_st(TargetStates, TargetStates_temp)
!
!  call destruct_basis_st(TargetStates_temp)
!
!  deallocate(TargetStates2el%b(1)%CI, TargetStates2el%b(1)%na, TargetStates2el%b(1)%nb, &
!    &TargetStates2el%b(1)%ma, TargetStates2el%b(1)%mb, TargetStates2el%b(1)%nuse)
!
!  TargetStates2el%b(1)%nusemax =abs(num_natorb_incl)
!  allocate(TargetStates2el%nuse(abs(num_natorb_incl)))  
!  allocate(TargetStates2el%na(abs(num_natorb_incl)))  
!  allocate(TargetStates2el%nb(abs(num_natorb_incl)))  
!  allocate(TargetStates2el%ma(abs(num_natorb_incl)))  
!  allocate(TargetStates2el%mb(abs(num_natorb_incl)))  
!
!  TargetStates2el%b(1)%nam = abs(num_natorb_incl)
!  do j=1, abs(num_natorb_incl)
!    TargetStates2el%b(1)%nuse(j) = TargetStates%Nmax-abs(num_natorb_incl)-1+j
!    TargetStates2el%b(1)%na(j) = TargetStates%Nmax-abs(num_natorb_incl)-1+j
!    TargetStates2el%b(1)%nb(j) = TargetStates%Nmax-abs(num_natorb_incl)-1+j
!    TargetStates2el%b(1)%ma(j) = 
!  enddo
!


  !Unrearrange required even if rearrange is not called - to represent natorbs back in terms of bst_nr rather than TargetStates
  call unrearrange(natorbs,TargetStates_unrearranged)

  call destruct(bst_nr)
  
  if(spheroidal) then
    call construct_spheroidal(bst_data, basis_1e)  
    call expand_basis_m(bst_data)

    ! bst_sub: short-ranged Laguerre functions to replace orbitals.
    if (data_in%use_sub <= 0) then !no sub basis
      call construct_spheroidal(bst_MSC, basis_2e)  
      call expand_basis_m(bst_MSC)
    else
      call construct_spheroidal(bst_temp, basis_2e)  
      call expand_basis_m(bst_temp)
      call construct_spheroidal(bst_sub, sub_basis_2e)  
      call expand_basis_m(bst_sub)
      call combine_basis_nr(bst_MSC,bst_temp,bst_sub,natorbs%basis_type)
    end if
  else
    call construct(bst_data, data_in)  ! Construct temporary basis from data.in inputs called bst_data
    call construct(bst_MSC, dataMSC, dataMSC%alpha_nl )   ! Construct temporary basis from dataMSC inputs called bst_MSC
  endif
  
  call destruct_basis_st(TargetStates)
  call hybrid_basis_natorbs(bst_nr, bst_MSC, bst_data, natorbs, TargetStates, num_natorb_incl, lag_ham1el_m)
  
 print*, '   FINISHED MAKING NATURAL ORBITALS'

  return
end subroutine natorbitals
  

subroutine unrearrange(natorbs,TS)
  !convert natorbs from representation in terms of TS basis to representation in terms of original bst_nr basis
  use state_class
  implicit none
  type(basis_state), intent(inout) :: natorbs
  type(basis_state), intent(in) :: TS
  integer :: num_natorb, Nmax, j, i, ii, jj, nam, nam_new, nam_TS, numLagMax
  type(state), pointer :: natorb, TSorb
  integer :: n
  integer, dimension(:), allocatable :: na, ma
  real*8, dimension(:), allocatable :: CI
  
  num_natorb = natorbs%Nmax
  Nmax  = TS%Nmax

  numLagMax = 0
  do i=1, Nmax
    numLagMax = max(numLagMax, maxval(TS%b(i)%na))
  enddo

  allocate( na(numLagMax), ma(numLagMax), CI(numLagMax) )
  do j=1, num_natorb
    natorb => natorbs%b(j)
    nam = natorb%nam

    na = 0
    ma = 0
    nam_new = 0
    CI = 0.0d0

    do i=1, nam
      TSorb => TS%b( natorb%na(i) ) !one-electron state 
      nam_TS = TSorb%nam

      do ii=1, nam_TS !loop through indices to Laguerre functions underlying the TargetStates orbital
        
        n = TSorb%na(ii) !index to underlying Laguerre function

        if(.not.any(na(1:nam_new) == n)) then !We have not already added this Laguerre function 
          nam_new = nam_new + 1
          na(nam_new) = n
          ma(nam_new) = TSorb%ma(ii)
          jj = nam_new
        else
          do jj=1, nam_new
            if(na(jj) == n) exit
          enddo
        endif

        CI(jj) = CI(jj) + TSorb%CI(ii)*natorb%CI(i) !update the CI coefficient for this Laguerre function

      enddo ! Laguerre functions

    enddo ! One-electron orbitals

    deallocate(natorb%CI, natorb%na, natorb%ma)
    natorb%nam = nam_new
    allocate(natorb%CI(nam_new), natorb%na(nam_new), natorb%ma(nam_new))
    natorb%CI = CI(1:nam_new)
    natorb%na = na(1:nam_new)
    natorb%ma = ma(1:nam_new)
 
  enddo ! Natural orbitals

end subroutine unrearrange

subroutine hybrid_basis_natorbs(bst, bst1, bst2, natorbs, TargetStates, num_natorb_incl, lag_ham1el_m)
  !TargetStates will be comprised of bst_1, plus a number of natural orbitals which are represented in terms of both bst_1 and bst_2
  use sturmian_class
  use state_class
  use ovlpste1me
  use input_data
  implicit none
  type(basis_sturmian_nr), intent(in) :: bst1, bst2
  type(basis_sturmian_nr), intent(inout) :: bst
  type(basis_state), intent(in) :: natorbs
  type(basis_state), intent(out) :: TargetStates
  integer, intent(in) :: num_natorb_incl
  real*8, dimension(:,:,:), allocatable, intent(out) :: lag_ham1el_m

  integer :: n_bst, n_bst1, n_bst2, Nmax, indi, indj, ni, nj
  real*8 :: CIi, CIj
  integer :: basis_type, max_latop, lorb, n, k, l, par, m, n_natorb, nst1, nst2, m1, m2, ist, i, j, n_start, n_stop
  integer :: im, im_min, im_max
  
  type(sturmian_nr), pointer :: sturm, sturmi, sturmj
  type(state), pointer :: natorb, state1, state2

  logical :: hlike, state_copied, spheroidal

  real*8, dimension(1) :: CI
  integer, dimension(1) :: no, mo
  real*8 :: energy, overlap, Ham_1e
  integer, dimension(:,:), allocatable :: inum
  real*8 :: Hyl_Ham_1e, Hyl_overlap
  
  character*20, parameter :: symcharsUP = 'u.gSPDFGHIKLMNOQRTUV'
  character*20, parameter :: symcharsLOW = 'u.gspdfghiklmnoqrtuv'
  
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

  n_bst1 = bst1%n !second basis
  n_bst2 = bst2%n !first basis
  n_bst  = n_bst1 + n_bst2

  spheroidal = data_in%calculation_type == 2 .or. data_in%calculation_type == 3

  basis_type = natorbs%basis_type
    
  if(spheroidal) then !spheroidal code: the first basis goes first
    call combine_basis_nr(bst,bst2,bst1,basis_type)
    !combine_basis_nr calculated ortint and ham1el for spherical case - replace now for spheroidal case
    do i = 1, n_bst
       sturmi => bst%b(i)
       do j = i, n_bst
          sturmj => bst%b(j)
          overlap = Hyl_overlap(sturmi,sturmj)
          Ham_1e = Hyl_Ham_1e(sturmi,sturmj)
          bst%ortint(i,j) = overlap; bst%ham1el(i,j) = Ham_1e
          bst%ortint(j,i) = overlap; bst%ham1el(j,i) = Ham_1e
       end do
    end do

    n_start = n_bst2+1
    n_stop = n_bst

!    print*, 'n_bst1, n_bst2, n_bst:', n_bst1, n_bst2, n_bst
!    print*, 'n_start, n_stop:', n_start, n_stop
!    do i=1, n_bst
!      sturmi => bst%b(i)
!      if(i==n_start) print*, '===== N_START ====='
!      print*, 'i, k, l, al, m:', i, sturmi%k, sturmi%l, sturmi%alpha, sturmi%m
!      if(i==n_stop) print*, '===== N_STOP ====='
!    enddo

  else !Spherical code: the first basis goes second
    call combine_basis_nr(bst,bst1,bst2,basis_type)

    n_start = 1
    n_stop = n_bst1

  endif
 
  max_latop = get_max_L(bst) ! Largest l of combined Laguerre basis

  allocate(inum(-max_latop:max_latop,-1:1))
  inum = 0

  allocate(lag_ham1el_m(n_bst,n_bst,-max_latop:max_latop))
  lag_ham1el_m(:,:,:) = 0d0
  
  if(spheroidal) then
    Nmax = n_bst1
  else
    ! Maximum number of Molecular states possible from bst1 basis 
    Nmax = 0
    do n = 1, n_bst1
      lorb = get_ang_mom(bst1%b(n))
      Nmax = Nmax + (2*lorb + 1)
    end do
  endif

  if(num_natorb_incl == 0) Nmax = natorbs%Nmax
  if(num_natorb_incl < 0) Nmax = -num_natorb_incl

  !Nmax = num_natorb_incl

  print*, 'Nmax:', nmax

  if(allocated(e1me)) deallocate(e1me,ovlpst)
  allocate(e1me(Nmax, Nmax))   ! H2+ Ham matrix element = sum C <f(1)|K_1+V_1|i(1)>
  allocate(ovlpst(Nmax, Nmax)) ! Overlap one- or two-electron functions/states 
  e1me(:,:) = 0d0
  ovlpst(:,:) = 0d0 

  hlike=.true.
  call new_basis_st(TargetStates,Nmax,hlike,basis_type)
  
  if(num_natorb_incl <= 0) then
    
    do ist=1, Nmax
      call copy_st(TargetStates%b(ist), natorbs%b(ist))
    enddo

  else
        
    CI(1) = 1.0d0
    energy = 0.0d0
    ist = 0
    do n=n_start, n_stop
      sturm => bst%b(n)
      k = get_k(sturm)
      l = get_ang_mom(sturm)
      par = (-1)**l

      no(1) = n

      !different treatment if spherical or spheroidal:
      !spherical assumes bst contains only m=0 orbitals and explicitly loops over m = -l -> l here
      !spheroidal assumed bst already contains m = -1 -> l functions
      if(spheroidal) then
        im_min = 0
        im_max = 0
      else
        im_min = -l
        im_max = l
      endif

      do im=im_min, im_max

        if(spheroidal) then
          m = get_ang_mom_proj(sturm)
        else
          m = im
        endif

        ist = ist + 1
        state1 => TargetStates%b(ist)
        mo(1) = m
        inum(m,par) = inum(m,par) + 1

        state_copied = .false.
        do n_natorb=1, num_natorb_incl
          natorb => natorbs%b(n_natorb)
          if(natorb%inum == inum(m,par) .and. nint(natorb%m) == m .and. natorb%parity == par) then
            if(state_copied) error stop '*** ERROR in hybrid_basis_natorbs: state already copied'
            call copy_st(state1,natorb)
            state1%inum = k
            state1%l = l
            state_copied = .true.
            exit
          endif
        enddo

        if(state_copied) cycle
        call construct_st(state1, hlike, dble(m), par, 0.5d0, energy, k, 1, CI, no, mo)
        call set_n_majconf(state1,k+l)
        call set_l_majconf(state1,l)
      
        write(state1%label,'(X,I1,A1,A1,X)') k+l, symcharsLOW(l+4:l+4), symcharsUP(abs(m)+4:abs(m)+4)
      
      enddo
    enddo

  endif
 
!  print*, '===== HYBRID BASIS BEFORE REARRANGE ====='
!  do i=1, TargetStates%Nmax
!    state1 => TargetStates%b(i)
!    print*, 'i, label, n, l, m, p, inum, nam:', i, state1%label, state1%n, state1%l, state1%m, state1%parity, state1%inum, state1%nam
!  enddo

  if(spheroidal) then
    
    if(data_in%inc == 1) call rearrange(bst,get_max_L(bst),TargetStates,.FALSE.)

    ! Calculate matrix elements between orbital states.
    do nst1 = 1, Nmax
       state1 => TargetStates%b(nst1)

       do nst2 = nst1, Nmax
          state2 => TargetStates%b(nst2)

          overlap = 0d0; Ham_1e = 0d0

          do indi = 1, get_nam(state1)
             ni = get_na(state1,indi); CIi = get_CI(state1,indi)
             do indj = 1, get_nam(state2)
                nj = get_na(state2,indj); CIj = get_CI(state2,indj)
                overlap = overlap + CIi*CIj*bst%ortint(ni,nj)
                Ham_1e = Ham_1e + CIi*CIj*bst%ham1el(ni,nj)
             end do
          end do

          ovlpst(nst1,nst2) = overlap; e1me(nst1,nst2) = Ham_1e
          ovlpst(nst2,nst1) = overlap; e1me(nst2,nst1) = Ham_1e
       end do
    end do
  else

    call Hybrid_H1el_st(Nmax,lag_ham1el_m, basis_size_nr(bst), get_max_L(bst))
  
    do nst1 = 1, Nmax
       state1 => TargetStates%b(nst1)
       m1 = get_ang_mom_proj(state1)
  
       do nst2 = 1, nst1
          state2 => TargetStates%b(nst2)
          m2 = get_ang_mom_proj(state2)
  
          if ( m1 /= m2 .OR. get_par(state1) /= get_par(state2) ) cycle
  
          ovlpst(nst1,nst2) = ovlp_st(state1,state2,bst)
          ovlpst(nst2,nst1) =  ovlpst(nst1,nst2)
  
       end do
    end do
  
    if(data_in%inc == 1) call rearrange(bst,   get_max_L(bst),   TargetStates,.TRUE.,lag_ham1el_m)

  endif !spheroidal
    
  print*, '===== FINAL HYBRID BASIS ====='
  do i=1, TargetStates%Nmax
    state1 => TargetStates%b(i)
    print*, 'i, label, n, l, m, p, inum:', i, state1%label, state1%n, state1%l, state1%m, state1%parity, state1%inum
  enddo
  
end subroutine hybrid_basis_natorbs

end module natural_orbitals_module
  


  
