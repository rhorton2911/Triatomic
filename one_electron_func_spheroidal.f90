subroutine construct_1el_basis_spheroidal(nStates)
  !
  ! JS
  !
  use input_data
  use grid_radial
  use MPI_module
  use one_electron_func
  use ovlpste1me
  use spheroidal_class
  use sturmian_class
  use natural_orbitals_module
  use target_states

  implicit none

  integer, intent(in) :: nStates

  type(basis_sturmian_nr) :: bst_1e,bst_2e,bst_sub
  type(state), pointer :: statei,statej
  type(sturmian_nr), pointer :: sturmi,sturmj
  character*9 :: fileName
  integer :: jFunc,nFunc, nst, n1,n2,n3, i,j,jState, k,l, par,m,mMin,mMax, ERR, nDebug
  integer :: indi,indj, ni,nj
  integer, dimension(data_debug%n_sturm_1e) :: vDebug
  integer, dimension(:), allocatable :: nArray, mArray
  real*8 :: R,Z1,Z2, CIi,CIj, overlap,Ham_1e, Hyl_overlap,Hyl_Ham_1e
  real*8, dimension(:), pointer :: rhoVec, sturmVec
  real*8, dimension(:), allocatable :: energies
  real*8, dimension(:,:), allocatable :: overlapMatrix, HMatrix, CIMatrix, overlapMatrix2
  logical :: hlike, isDebug=.false.
  integer :: jj
  real*8 :: norm

  !FOR THE DSYGVX SUBROUTINE
  real*8, dimension(:), allocatable :: WORK
  integer, dimension(:), allocatable :: IWORK, IFAIL
  integer :: LWORK, NFOUND
  real*8, external :: DLAMCH


  hlike = .true.   ! True because we create the one electron states first.
  rhoVec => grid%gridr
  R = data_in%Rd
  Z1 = data_in%Z1
  Z2 = data_in%Z2

  ! Create space for the target states.
  call new_basis(TargetStates, nStates, hlike, data_in%calculation_type)
  jState = 0
  if (myid==0) write(*,'(A,I4,A/)') 'Created space for', nStates, ' one-electron target states.'

  if (R /= 0d0) then
     TargetStates%en_ion = Z1*Z2/R
  else ! Atomic Case
     TargetStates%en_ion = 0d0 
  end if

  ! Construct the radial basis functions.
  call construct_spheroidal(bst_nr, basis_1e)
  nFunc = basis_size(bst_nr)
  do i = 1, nFunc
     sturmi => bst_nr%b(i)
     do j = 1, i
        sturmj => bst_nr%b(j)
        overlap = Hyl_overlap(sturmi, sturmj)
        Ham_1e = Hyl_Ham_1e(sturmi, sturmj)
        
        bst_nr%ortint(i,j) = overlap; bst_nr%ham1el(i,j) = Ham_1e
        bst_nr%ortint(j,i) = overlap; bst_nr%ham1el(j,i) = Ham_1e
     end do
  end do


  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
  ! Debugging prep.
  if ( data_debug%list_sturm_1e > 0 ) then
     isDebug = .true.
     open(11, file='1e/sturm')
     write(11,'(4A4,A8,A12)') '# n', 'k','l','m', 'alpha', 'first value'
  end if
  nDebug = min( data_debug%n_sturm_1e, nFunc )
  vDebug(:) = 0
  if ( nDebug > 0 ) then
     vDebug(1:nDebug) = abs( data_debug%v_sturm_1e(:) )
     open(13, file='1e/matrices')
     write(13,'(A/A,99(S,I0.2,SP,2I2,2X))')'# Overlap matrix', '# ', (/( get_k(bst_nr%b(vDebug(i))),get_ang_mom(bst_nr%b(vDebug(i))),get_ang_mom_proj(bst_nr%b(vDebug(i))), i=1,nDebug )/)
  end if

  ! Debugging.
  do i = 1, nFunc
     sturmi => bst_nr%b(i)
     k = get_k(sturmi)
     l = get_ang_mom(sturmi)
     m = get_ang_mom_proj(sturmi)

     if (isDebug) write(11,'(4I4,F8.4,F12.8)') i, k,l,m, get_alpha(sturmi), value(sturmi,0)
     
     if (nDebug > 0) then
        if ( any(i==data_debug%v_sturm_1e(:)) ) then
           sturmVec => fpointer(sturmi)
           write(fileName,'(A,I0.2,SP,2I2)') '1e/', k,l,m
           open(12, file=fileName)
           do j = get_minf(sturmi), get_maxf(sturmi)
              write(12,'(2Es20.12)') rhoVec(j), sturmVec(j)
           end do
           close(12) ; print*, 'Made ', fileName
        end if
     end if

     if ( any(i==vDebug(:)) ) write(13,'(99F8.4)') (/( bst_nr%ortint(i,vDebug(j)), j=1,nDebug )/)
  end do
  if (isDebug) then
     close(11) ; print*, 'Made 1e/sturm'
  end if
  if ( nDebug > 0 ) then
     write(13,'(/A/A,99(S,I0.2,SP,2I2,2X))')'# One-electron Hamiltonian matrix', '# ', (/( get_k(bst_nr%b(vDebug(j))),get_ang_mom(bst_nr%b(vDebug(j))),get_ang_mom_proj(bst_nr%b(vDebug(j))), j=1,nDebug )/)
     do i = 1, nFunc
        if ( any(i==vDebug(:)) ) write(13,'(99F8.4)') (/( bst_nr%ham1el(i,vDebug(j)), j=1,nDebug )/)
     end do
     close(13) ; print*, 'Made 1e/matrices'
  end if
  ! End debugging.


  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
  ! Loop through each (m,par) symmetry and make one-electron states.
  allocate( nArray(nFunc), mArray(nFunc) )
  do m = data_in%Mt_min, data_in%Mt_max   ! As specified by the input file.
     do par = 1,-1,-1   ! 1,-1 if parity is a good quantum number, 0 if not.
        if (data_in%good_parity) then
           if (par == 0) cycle
           nst = data_in%nst(m,par)
        else
           if (par /= 0) cycle
           nst = sum( data_in%nst(m,:) )
        end if
        if (myid==0) write(*,'(3(A,I2))') 'Symmetry (m,par) = (',m,',',par,')'

        if(nst > 0) then   ! States will be constructed.

           !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
           ! This block of code determines which of the Sturmian functions
           ! will be used in constructing the particular (m,par) state.

           nArray(:) = -1   ! The indices of the extended basis functions.
           mArray(:) = -1   ! The corresponding values of m.
           
           jFunc = 0
           do j = 1, nFunc   ! Loop through the functions.
              if (.not.data_in%good_parity .or. (data_in%good_parity .and. (-1)**get_ang_mom(bst_nr%b(j))==par)) then
                 if (get_ang_mom_proj(bst_nr%b(j)) == m) then
                    jFunc = jFunc + 1
                    nArray(jFunc) = j
                    mArray(jFunc) = m
                 end if
              endif
           enddo ! func
           if (myid==0) write(*,'(3(A,I3))') 'Making ', nst, ' states from ', jFunc,' basis functions.'

           !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
           ! This code builds the Hamiltonian matrix with the chosen functions.
           allocate( overlapMatrix(jFunc,jFunc), overlapMatrix2(jFunc,jFunc), HMatrix(jFunc,jFunc), CIMatrix(jFunc,jFunc), energies(jFunc) )

           overlapMatrix(:,:) &
                = bst_nr%ortint( nArray(1:jFunc), nArray(1:jFunc) )
           HMatrix(:,:) &
                = bst_nr%ham1el( nArray(1:jFunc), nArray(1:jFunc) )
           if (R > 0d0) HMatrix(:,:) = HMatrix(:,:) + Z1*Z2/R*overlapMatrix(:,:)

           overlapMatrix2 = overlapMatrix

           !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
           ! This block of code solves the generalised eigenvalue problem Hv=abv
           ! for CI coefficient vectors v and corresponding eigenenergies a.
           !call rsg( jFunc, jFunc, HMatrix, overlapMatrix, energies, 1, CIMatrix, ERR )
           
           allocate(IFAIL(jFunc), IWORK(5*jFunc))
           allocate(WORK(1))
           LWORK = -1
           call dsygvx( 1, 'V', 'I', 'U', jFunc, HMatrix, jFunc, overlapMatrix, jFunc, 0.0d0,0.0d0, 1,jFunc, 2*DLAMCH('S'), &
             &NFOUND, energies, CIMatrix, jFunc, WORK, LWORK, IWORK, IFAIL, ERR)
           
           LWORK = WORK(1)
           deallocate(WORK)
           allocate(WORK(LWORK))
           call dsygvx( 1, 'V', 'I', 'U', jFunc, HMatrix, jFunc, overlapMatrix, jFunc, 0.0d0,0.0d0, 1,jFunc, 2*DLAMCH('S'), &
             &NFOUND, energies, CIMatrix, jFunc, WORK, LWORK, IWORK, IFAIL, ERR)

           deallocate(WORK, IFAIL, IWORK)

           if(myid==0) then             
              write(*,'(/A)') 'Energies in a.u. :'
              write(*,'(5F20.5)') energies
              write(*,*) 'Energies in eV   :'
              write(*,'(5F20.5)') energies * data_in%eV  !  27.211
           endif
             
           !Liam added - remove very small CI coefficients
           where(abs(CIMatrix) < data_in%CI_min) CIMatrix = 0.0d0

           ! Build the target states requested in data.in.
           do j = 1, min(jFunc, nst)
              jState = jState + 1
              if( sum(CIMatrix(:,j)) < 0d0 ) CIMatrix(:,j) = -CIMatrix(:,j)
            
              !Liam added - renormalise CIMatrix after we set some elements to zero
              norm=0.0d0
              do i=1,jFunc
                do jj=1, jFunc
                  norm = norm + CIMatrix(i,j)*CIMatrix(jj,j)*overlapMatrix2(i,jj)
                enddo
              enddo
              !renormalise CIVec
              CIMatrix(:,j) = CIMatrix(:,j) / sqrt(norm)

              call construct_st(TargetStates%b(jState), hlike, dble(m),par, 0.5d0,energies(j), j,jFunc, CIMatrix(:,j), nArray(1:jFunc),mArray(1:jFunc) )

              if (m > 0) then   ! We have degeneracy.
                 jState = jState + 1
                 call construct_st(TargetStates%b(jState), hlike, dble(-m),par, 0.5d0,energies(j), j,jFunc, CIMatrix(:,j), nArray(1:jFunc),-mArray(1:jFunc) )
              endif
           enddo

           deallocate( overlapMatrix, overlapMatrix2, HMatrix, CIMatrix, energies )

        endif
     enddo ! par
  enddo ! m

  if (myid==0) write(*,'(//2(A,I4)//)') 'Number of states created:', jState, ' and requested:', nStates

  if (jState /= nStates) stop

  call sort_by_energy_basis_st(TargetStates)
  call calc_spectro_factors(TargetStates, bst_nr, basis_1e%labot,basis_1e%latop)
  if (myid==0) call print_energy_basis_st(TargetStates)

  if (data_in%hlike) then
     if (data_in%inc==1) then
        call rearrange(bst_nr,get_max_L(bst_nr), TargetStates,.FALSE.)
        call spheroidal_separate(TargetStates, bst_oid)
     end if
     return   ! The end for hydrogen-like targets.
  end if

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
  ! Preparing for two-electron targets: setting up the second basis.

  ! bst_1e: old basis functions if we want 1e states, otherwise empty.
  if (data_in%use_MSC == 0) then   ! We discard the one-electron target states.
     call new_basis(bst_1e,0)   ! Dummy basis.
  else
     nFunc = basis_size(bst_nr) - get_num_m0(bst_nr)   ! Count of m=0 functions.
     call expand_basis_m(bst_nr)   ! Expand m -> -m functions.
     call new_basis(bst_1e, basis_size(bst_nr))
     call copy_basis(bst_1e, bst_nr)

     ! Shift the state links so that -m states point to -m functions.
     do jState = 1, nStates
        statej => TargetStates%b(jState)
        if (get_ang_mom_proj(statej) < 0) then
           do j = 1, get_nam(statej)
              statej%na(j) = statej%na(j) + nFunc
           end do
        end if
     end do
        
  end if
  n1 = basis_size(bst_1e)
  call destruct(bst_nr)

  ! bst_2e: Laguerre functions to diagonalise the two-electron Hamiltonian.
  if (data_in%use_MSC < 0) then
     if (myid==0) write(*,'(/A)') 'Using H2+ states to diagonalise H2'
     call new_basis(bst_2e,0)   ! Dummy basis.
  else   ! This option is for using a basis of mainly Laguerre functions.
     if (myid==0) write(*,'(/A)') 'Constructing main MSC basis.'
     call construct_spheroidal(bst_2e, basis_2e)   ! Low alpha "main" basis.
     call expand_basis_m(bst_2e)
  end if
  n2 = basis_size(bst_2e)

  ! bst_sub: short-ranged Laguerre functions to replace orbitals.
  if (data_in%use_sub <= 0) then
     call new_basis(bst_sub, 0)   ! Dummy basis.
  else
     if (myid==0) write(*,'(/A)') 'Constructing sub MSC basis.'
     call construct_spheroidal(bst_sub, sub_basis_2e)   ! High alpha.
     call expand_basis_m(bst_sub)
  end if
  n3 = basis_size(bst_sub)
  nFunc = n1 + n2 + n3


  ! Combine all the bases into one.
  call new_basis(bst_nr, nFunc)
  j = 0
  do i = 1, basis_size(bst_1e)   ! Basis for one-electron targets.
     j = j + 1
     call copy( bst_nr%b(j), bst_1e%b(i) )
  end do
  do i = 1, basis_size(bst_2e)   ! "Main" two-electron basis.
     j = j + 1
     call copy( bst_nr%b(j), bst_2e%b(i) )
  end do
  do i = 1, basis_size(bst_sub)   ! "Sub" two-electron basis.
     j = j + 1
     call copy( bst_nr%b(j), bst_sub%b(i) )
  end do

!!$  call write_basis(bst_1e, 'bst_1e')
!!$  call write_basis(bst_2e, 'bst_2e')
!!$  call write_basis(bst_sub,'bst_sub')
!!$  call write_basis(bst_nr, 'bst_nr')
  call destruct(bst_1e)
  call destruct(bst_2e)
  call destruct(bst_sub)

  ! Calculate matrix elements between basis functions.
  do i = 1, nFunc
     sturmi => bst_nr%b(i)
     do j = i, nFunc
        sturmj => bst_nr%b(j)
        overlap = Hyl_overlap(sturmi,sturmj)
        Ham_1e = Hyl_Ham_1e(sturmi,sturmj)
        bst_nr%ortint(i,j) = overlap; bst_nr%ham1el(i,j) = Ham_1e
        bst_nr%ortint(j,i) = overlap; bst_nr%ham1el(j,i) = Ham_1e
     end do
  end do


  ! Make the orbital states from the basis depending on the representation.
  if (data_in%use_MSC < 0) then   ! One-electron state representation.
     if (myid==0) write(*,'(/A)') 'Using the above H2+ states to diagonalise H2'
     if (data_in%use_sub <= 0) then
        j = nFunc + 1   ! Make sure no replacements will be made.
        if (myid==0) write(*,'(A)') 'No orbitals to be replaced'
     else
        j = n1 + n2 + data_in%use_sub   ! First of the substitute functions.
        if (myid==0) write(*,'(A,I2,A)') 'Replacing ', nFunc-j+1, ' states with short-ranged Laguerre functions'
     end if
     call make_state_basis(bst_nr,TargetStates, j)

  else   ! Laguerre/Hylleraas function representation.
     if (myid==0) write(*,'(/A)') 'Using a new Laguerre basis to diagonalise H2'
     j = data_in%use_MSC   ! Number of 1e states to insert into the basis.
     if (myid==0) write(*,'(A,I2,A)') 'Replacing ', j, ' Laguerres with H2+ states'
     call make_Hyl_basis(bst_nr,TargetStates, n1,j,n3)

  end if
  call copy(TargetStates_unrearranged,TargetStates)

  if (data_in%inc==1) call rearrange(bst_nr,get_max_L(bst_nr), TargetStates, .false.)
!  call write_basis(bst_nr, 'bst_2_rearr_1e')

  ! Calculate matrix elements between orbital states.
  jState = basis_size(TargetStates)
  allocate( ovlpst(jState,jState), e1me(jState,jState) )
  do i = 1, jState
     statei => TargetStates%b(i)

     norm = 0.0d0
     do j = i, jState
        statej => TargetStates%b(j)

        overlap = 0d0; Ham_1e = 0d0
        do indi = 1, get_nam(statei)
           ni = get_na(statei,indi); CIi = get_CI(statei,indi)
           do indj = 1, get_nam(statej)
              nj = get_na(statej,indj); CIj = get_CI(statej,indj)
              overlap = overlap + CIi*CIj*bst_nr%ortint(ni,nj)
              Ham_1e = Ham_1e + CIi*CIj*bst_nr%ham1el(ni,nj)
           end do
        end do

        ovlpst(i,j) = overlap; e1me(i,j) = Ham_1e
        ovlpst(j,i) = overlap; e1me(j,i) = Ham_1e
     end do
  end do


end subroutine construct_1el_basis_spheroidal

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine write_states(states, nStates)
  use input_data
  use one_electron_func
  use state_class
  use sturmian_class
  implicit none

  type(basis_state), intent(in) :: states
  integer, intent(in) :: nStates

  type(state), pointer :: statei
  type(sturmian_nr), pointer :: sturmj
  integer :: matLen,matWid, i,j, l,lMax
  integer, dimension(:), allocatable :: mVec
  integer, dimension(:,:), allocatable :: lMat
  real*8, dimension(:,:), allocatable :: coeffMat


  lMax = -1
  matLen = min( nStates, states%Nmax )
  matWid = 0

  do i = 1, matLen
     statei => states%b(i)
     matWid = max( matWid, get_nam(statei) )
  end do
  allocate( lMat(matLen,matWid), mVec(matLen), coeffMat(matLen,matWid) )
  lMat(:,:) = 0
  mVec(:) = 0
  coeffMat(:,:) = 0d0

  do i = 1, matLen
     statei => states%b(i)
     mVec(i) = get_ang_mom_proj(statei)

     do j = 1, get_nam(statei)
        sturmj => bst_nr%b( get_na(statei,j) )
        l = get_ang_mom(sturmj)
        lMat(i,j) = l
        if ( l > lMax ) lMax = l
        coeffMat(i,j) = get_CI(statei,j)
     end do
  end do

  call write_harmonics( data_in%calculation_type,data_in%Rd, lMax,maxval(mVec(:)), matLen,matWid, lMat,mVec,coeffMat )


end subroutine write_states

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

real*8 function Hyl_overlap(sturmp, sturm)
  !
  use grid_radial
  use input_data
  use spheroidal_class
  use sturmian_class
  implicit none

  type(sturmian_nr), intent(in) :: sturmp, sturm

  type(spheroidal_fn) :: fnp,fn
  integer :: dk
  real*8 :: R, kp,k, lp,l, mp,m, lamp,lam, tmp,Yint
  

  Hyl_overlap = 0d0
  R = data_in%Rd
  lamp = 2d0*get_alpha(sturmp)
  lam = 2d0*get_alpha(sturm)

  ! Determine the overlap analytically or numerically.
  if (data_in%calculation_type==2 .and. lamp==lam) then   ! Can do analytically.
     lp = dble( get_ang_mom(sturmp) ); l = dble( get_ang_mom(sturm) )
     mp = dble( get_ang_mom_proj(sturmp) ); m = dble( get_ang_mom_proj(sturm) )
     if (mp /= m) return
     
     kp = dble( get_k(sturmp) ) ; k = dble( get_k(sturm) )
     dk = abs( nint(kp-k) )
     mp = abs(mp); m = abs(m)
  
     if (lp == l) then
        k = min(kp,k)
        if (dk == 0) then
           Hyl_overlap = (6d0*k*k - 6d0*k + 6d0*k*m - 3d0*m + m*m + 2d0)/lam/lam + (2d0*k - 1d0 + m)*R/lam + R*R/6d0
        elseif (dk == 1) then
           Hyl_overlap = -(4d0*k + 2d0*m + lam*R) * dsqrt(k*(k+m)) /lam/lam
        elseif (dk == 2) then
           Hyl_overlap = dsqrt(k * (k+m) * (k+1) * (k+1+m)) /lam/lam
        end if
     end if

     tmp = Yint(lp,mp, 2d0,0d0, l,m)   ! (eta^2-1/3) term.
     if ( dk==0 .and. tmp/=0d0 ) Hyl_overlap = Hyl_overlap - R*R/6d0*tmp

  else   ! Must calculate numerically.
     call convert_from_sturm_to_oid(sturmp, fnp)
     call convert_from_sturm_to_oid(sturm, fn)
     Hyl_overlap = oid_overlap(fnp,fn)

  end if


end function Hyl_overlap

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

real*8 function Hyl_Ham_1e(sturmp, sturm)
  !
  use grid_radial
  use input_data
  use sturmian_class
  implicit none

  type(sturmian_nr), intent(in) :: sturmp, sturm

  integer :: calc_type, i1,i2, dk
  real*8 :: al,R, Zplus,Zminus, kp,lp,mp, k,l,m, tmp,Yint
  real*8, dimension(:), pointer :: rhoVec,weightVec, fpVec,fVec,gVec
  real*8, dimension(:), allocatable :: HVec, RVec


  Hyl_Ham_1e = 0d0
  tmp = 0d0
  calc_type = data_in%calculation_type
  R = data_in%Rd
  Zplus = data_in%Z1 + data_in%Z2
  Zminus = abs(data_in%Z1 - data_in%Z2)

  kp = dble( get_k(sturmp) );            k = dble( get_k(sturm) )
  lp = dble( get_ang_mom(sturmp) );      l = dble( get_ang_mom(sturm) )
  mp = dble( get_ang_mom_proj(sturmp) ); m = dble( get_ang_mom_proj(sturm) )
  al = get_alpha(sturm)
  if (mp /= m) return   ! Conservation of angular momentum projection always.
  if (lp/=l .and. data_in%good_parity) return   ! Conservation of angular momentum for homogeneous diatomics.
  m = abs(m)

  rhoVec => grid%gridr
  weightVec => grid%weight
  fpVec => fpointer(sturmp)
  fVec => fpointer(sturm)
  i1 = max( get_minf(sturmp), get_minf(sturm) )
  i2 = min( get_maxf(sturmp), get_maxf(sturm) )
  allocate( HVec(i1:i2), RVec(i1:i2) )
  HVec(:) = 0d0
  RVec(:) = 1d0 / (rhoVec(i1:i2)+R)

  ! Solve analytically or numerically depending on the functions coming in.
  if (lp == l) then
     if (calc_type==2 .and. get_alpha(sturmp)==al) then   ! Analytical.
        dk = nint( abs(kp-k) )
        k = min(kp,k)
        
        if (dk == 0) then
           tmp = (k*k - k + k*m - m/2d0 + 1d0 - 2d0*Zplus*R + 2d0*l*(l+1d0) + (2d0*k-1d0+m)*(al*R-2*Zplus/al)) /4d0
        elseif (dk == 1) then
           tmp = dsqrt(k*(k+m)) * (al*R+2d0*Zplus/al) /4d0
        elseif (dk == 2) then
           tmp = -dsqrt( k*(k+m)*(k+1d0)*(k+1d0+m) ) /8d0
        end if
        if (m > 0d0) HVec(:) = -m*m*R/8d0 * fVec(i1:i2)*RVec(:)
        
     elseif (calc_type == 2) then
        HVec(:) = ( -al*al/2d0*rhoVec(i1:i2)*rhoVec(i1:i2) + al*(k-0.5d0+m/2d0-al*R/2d0-Zplus/al)*rhoVec(i1:i2) + (l*(l+1d0) - m*m/4d0 + al*R*(2d0*k-1d0+m-Zplus/al))/2d0 - m*m*R/8d0*RVec(:) ) * fVec(i1:i2)
        gVec => gpointer(sturm)
        HVec(:) = HVec(:) - rhoVec(i1:i2)/2d0*gVec(i1:i2)
        
     elseif (calc_type == 3) then
        stop 'Needs amending due to new definition of gVec'
        HVec(:) = ( -al*al/2d0*rhoVec(i1:i2) + al*(k-1d0+m/2d0) - Zplus + (k-1d0 + l*(l+1d0) - m/2d0*(m/2d0-1d0) + R*(al+Zplus))/2d0*RVec(:) - (m*m/4d0-1d0)*R/2d0*RVec(:)*RVec(:) ) * fVec(i1:i2)
        if ( k > 1d0 ) then
           gVec => gpointer(sturm)
           HVec(:) = HVec(:) - sqrt((k-1d0)*(k+m-1d0))/2d0*RVec(:)*gVec(i1:i2)
        end if
        HVec(:) = HVec(:) * RVec(:)   ! 1/(rho+R) from fpVec.
     
     end if
  else   ! lp/=l ==> Z1/=Z2. Calculate the heterogeneous term.
!DF 21-03-2019 Check the sign of the below term , it i s"+" in Yevgeniy code: Horbitals.f90 H1el_Ham 
     if (kp==k) tmp = tmp - R/2d0*Zminus*Yint(lp,m, 1d0,0d0, l,m)
  end if
  Hyl_Ham_1e = tmp + sum(fpVec(i1:i2) * HVec(:) * weightVec(i1:i2))


end function Hyl_Ham_1e

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine spheroidal_separate(states, statesOid)
  !
  use grid_radial
  use input_data
  use one_electron_func
  use state_class
  use sturmian_class
  use target_data
  implicit none

  type(basis_state), intent(in) :: states
  type(spheroidal_basis), intent(out) :: statesOid

  type(state), pointer :: tmpState
  type(sturmian_nr), pointer :: tmpSturm
  type(spheroidal_fn), pointer :: tmpFn
  integer :: jState,nStates, nTerms, n,lambda, i1,i2, j
  integer, dimension(:), allocatable :: naVec, lVec
  real*8 :: R,c2,A,stdev, l,m
  real*8, dimension(:), pointer :: vec1,vec2, weightVec
  real*8, dimension(:), allocatable :: AVec, DVec
  real*8, dimension(:,:), allocatable :: ratioMat


  R = data_in%Rd
  weightVec => grid%weight
  nStates = basis_size(states)
  call create_spheroidal_basis(statesOid, nStates)

  write(*,'(//A/)') 'Separating one-electron target states into spheroidal functions'
  do jState = 1, nStates
     tmpState => states%b(jState)
     n = get_n_majconf(tmpState); lambda = get_l_majconf(tmpState)
     m = dble( get_ang_mom_proj(tmpState) )
     c2 = get_elec_energy(states,jState) * R*R/2d0

     write(*,'(A,I2,3A,F12.8,A)',advance='no') 'State #',jState, ' (',get_label(tmpState),') with c^2 = ',c2, ' and A = '
     call get_target_sepconst( data_target, n,lambda,nint(m), A,j )
     if( j /= 0 ) then
        write(*,'(F16.12,A,I1,A)') A, ' [', j, ']'
     else
        print*, 'NO DATA FOUND'
     endif

     nTerms = get_nam(tmpState)
     allocate( naVec(nTerms), DVec(nTerms), lVec(nTerms) )
     naVec(:) = (/( get_na(tmpState,j), j=1,nTerms )/)

     tmpSturm => bst_nr%b(naVec(1)) ; vec2 => fpointer(tmpSturm)
     i1 = get_minf(tmpSturm); i2 = get_maxf(tmpSturm)
     do j = 2, nTerms
        tmpSturm => bst_nr%b(naVec(j))
        i1 = max( i1, get_minf(tmpSturm) ) ; i2 = min( i2, get_maxf(tmpSturm) )
     end do
     allocate( ratioMat(nTerms-1,i1:i2), AVec(i1:i2) )
     
     do j = 1, nTerms-1
        l = get_ang_mom( bst_nr%b(naVec(j)) ); lVec(j) = l
        vec1 => vec2; vec2 => fpointer( bst_nr%b(naVec(j+1)) )
        DVec(j) = sum( vec1(i1:i2) * weightVec(i1:i2) )
        ratioMat(j,:) = vec1(i1:i2) / vec2(i1:i2)
        
        AVec(:) = -dble(l*(l+1)) + c2*( 2d0*(l*l+l+m*m-1d0)/(2d0*l-1d0)/(2d0*l+3d0) - dsqrt((l+1d0-m)*(l+1d0+m)*(l+2d0-m)*(l+2d0+m)/(2d0*l+1d0)/(2d0*l+5d0))/(2d0*l+3d0)/ratioMat(j,:) )
        if (j > 1) AVec(:) = Avec(:) - c2*dsqrt((l-1d0-m)*(l-1d0+m)*(l-m)*(l+m)/(2d0*l-3d0)/(2d0*l+1d0))/(2d0*l-1d0)*ratioMat(j-1,:)

        A = sum(AVec(:)) / dble(i2-i1+1)
        stdev = dsqrt( sum((AVec(:)-A)**2) )
        write(*,'(8X,A,I2,8X,A,F12.8,8X,A,F16.12,A,F16.12)') 'l =',nint(l), 'D = ',DVec(j), 'A = ',A, ' +/-', stdev
     end do

     j = nTerms
     l = get_ang_mom( bst_nr%b(naVec(j)) ); lVec(j) = l
     DVec(j) = sum( vec2(i1:i2) * weightVec(i1:i2) )
     stdev = dsqrt( sum(DVec(:)*DVec(:)) )
     write(*,'(8X,A,I2,8X,A,F12.8,8X,A,F12.8/)') 'l =',nint(l), 'D = ',DVec(j), 'normalised by dividing by', stdev
     DVec(:) = DVec(:) / stdev

     j = maxloc( abs(DVec(:)), 1 )
     tmpSturm => bst_nr%b(naVec(j))
     vec1 => fpointer(tmpSturm)
     i1 = get_minf(tmpSturm); i2 = get_maxf(tmpSturm)
     AVec(:) = vec1(i1:i2) / DVec(j)
     tmpFn => get_spheroidal_fn( statesOid, jState )
     call create_spheroidal_fn( tmpFn, n,get_energy(tmpState),lambda,nint(m), i1,i2,grid%nr,AVec, nTerms,lVec,DVec )

     deallocate( naVec, DVec, lVec, ratioMat, AVec )

  end do


end subroutine spheroidal_separate

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine make_state_basis(bst,states, jSub)
  !
  use MPI_module
  use state_class
  use sturmian_class
  implicit none

  type(basis_sturmian_nr), intent(in) :: bst
  type(basis_state), intent(inout) :: states
  integer, intent(in) :: jSub   ! Index of the first "sub" function in bst.

  type(state), pointer :: statei
  type(sturmian_nr), pointer :: sturmj
  logical :: found, hlike
  integer :: i,j, k,l,m, par, ncm,nVec(1),mVec(1)
  real*8 :: spin,E, CIVec(1)


  ! Dummy variables for any orbitals created from "sub" basis functions.
  hlike = .true.
  spin = 0.5d0
  E = 0d0
  ncm = 1
  CIVec(:) = 1d0

  ! Loop through each one-electron state to see if it should be replaced.
  do i = 1, basis_size(states)
     statei => states%b(i)
     k = get_inum_st(statei)
     l = get_l_majconf(statei)
     m = get_ang_mom_proj(statei)

     ! Look through the "sub" functions (start with jSub) for a matching k,l,m.
     j = jSub-1; found = .false.
     do while (j<basis_size(bst) .and. .not.found)
        j = j + 1
        sturmj => bst%b(j)
        if (get_k(sturmj)==k .and. get_ang_mom(sturmj)==l .and. get_ang_mom_proj(sturmj)==m) found = .true.
     end do

     if (found) then   ! Remake the orbital from the found function.
        if (myid==0) write(*,'(4(A,I2))') 'State number ',i, ' replaced by Laguerre function k =',k, ', l =',l, ', m =',m
        call destruct_st(statei)
        par = (-1)**l
        nVec(:) = j   ! Link to the "sub" Laguerre function.
        mVec(:) = m
        call construct_st(statei,hlike,dble(m),par,spin,E,k,ncm,CIVec,nVec,mVec)
        call set_n_majconf(statei, k+l)
        call set_l_majconf(statei, l)
     end if
  end do


end subroutine make_state_basis


subroutine make_Hyl_basis(bst,states, nBst,nStates,nSub)
  !
  use MPI_module
  use state_class
  use sturmian_class
  implicit none

  type(basis_sturmian_nr), intent(in) :: bst   ! New expanded & combined basis.
  type(basis_state), intent(inout) :: states   ! In: 1e states/Out: 2e orbitals.
  integer, intent(in) :: nBst      ! Number of functions used for 1e states.
  integer, intent(in) :: nStates   ! Number of 1e states to substitute in.
  integer, intent(in) :: nSub      ! Number of Laguerres to substitute in.

  type(basis_state) :: statesTmp
  type(state), pointer :: stateNew,stateOld,statei,statej
  type(sturmian_nr), pointer :: sturmi,sturmj
  logical :: hlike, found
  integer :: nOld,nNew, i,j,n, k,l,m,par, basis_type,ncm,nVec(1),mVec(1)
  real*8 :: spin,E, CIVec(1)
  integer, dimension(:,:), allocatable :: inum
  logical :: dont_substitute, sub_by_l, sub_by_label, sub_by_dom
  character*8, parameter :: lchars='spdfghiklmnoqrtuv'
  character*8, parameter :: Mchars='SPDFGHIKLMNOQRTUV'

  real*8, dimension(:), allocatable :: ovlps
  integer :: indi, indj, ni, nj
  real*8 :: CIi, CIj
    
  if(nStates > states%Nmax) then
    if(myid==0) write(*,*) '*** ERROR: requested more 1-el states in the second basis than kept from the&
      & first diagonalisation'
    error stop
  endif


  !Liam: was playing around with different ways to substitute 1el states into basis
  !- turning these all off for now to use default approach
  dont_substitute = .false.
  sub_by_dom = .false.
  sub_by_l = .false.
  sub_by_label = .false.
  if(sub_by_dom) dont_substitute = .true.

  hlike = .true.
  basis_type = states%basis_type
  spin = 0.5d0
  E = 0d0
  ncm = 1
  CIVec(:) = 1d0

  nOld = basis_size(bst)   ! No. of functions in the combined basis.
  nNew = nOld - nBst - nSub   ! No. of functions in the main two-electron basis.

  if(dont_substitute) nNew = nNew + nStates

  call new_basis(statesTmp,nStates, hlike,basis_type)
  do i = 1, nStates
     call copy( statesTmp%b(i), states%b(i) )
  end do
  call destruct_basis_st(states)
  call new_basis(states,nNew, hlike,basis_type)

  allocate(inum(-get_max_L(bst):get_max_L(bst),-1:1))
  inum = 0
    
  n = 0

  if(dont_substitute) then
    do i=1, nStates
      n = n + 1
        call copy( states%b(i), statesTmp%b(i) )
    enddo
  endif

! These lines, and the end-if's below, order the orbitals similarly to in the spherical case but unfortunately don't solve the mysterious degeneracy problem.
!  do l = 0, get_max_L(bst)
!  do m = -l, l
  do i = 1, nOld - nBst - nSub
     sturmi => bst%b(nBst+i)
     k = get_k(sturmi)
     l = get_ang_mom(sturmi)
     m = get_ang_mom_proj(sturmi)
     par = (-1)**l
!     if (get_ang_mom(sturmi) /= l) cycle
!     if (get_ang_mom_proj(sturmi) /= m) cycle

    inum(m,par) = inum(m,par) + 1

     n = n+1
     stateNew => states%b(n)

!     print*, 'LAGUERRE i, m, par, inum:', i, m, par, inum(m,par)

     ! First, see if it can replaced by a one-electron target state.
     j = 0; found = .false.
     do while ( j<nStates .and. .not.found .and. .not.dont_substitute)
        j = j + 1
        stateOld => statesTmp%b(j)
        if(sub_by_l) then
          if ( get_inum_st(stateOld)==k .and. get_l_majconf(stateOld)==l .and. get_ang_mom_proj(stateOld)==m ) found = .true.
        elseif(sub_by_label) then
          if ( get_n_majconf(stateOld)==k+l .and. get_l_majconf(stateOld)==l .and. get_ang_mom_proj(stateOld)==m ) found = .true.
        else
          if ( get_inum_st(stateOld)==inum(m,par) .and. get_par(stateOld)==par .and. get_ang_mom_proj(stateOld)==m ) found = .true.
        endif
        !if ( k==1 .and. get_l_majconf(stateOld)==l .and. get_par(stateOld)==par .and. get_ang_mom_proj(stateOld)==m ) found = .true.
!        print*, '   -> H2+ j, label, m, par, l, inum:', j, stateOld%label, stateOld%m, stateOld%parity, stateOld%l, stateOld%inum
!        if(found) print*, '        :: found!'
     end do

     ! Simply copy the target state into its corresponding orbital spot.
     if (found) then
        if (myid==0) write(*,'(4(A,I2))') 'Laguerre function k =',k, ', l =',l, ', m =',m, ' replaced by one-electron state ',j
        call copy(stateNew, stateOld)
        stateNew%inum = k
        stateNew%l = l
        stateNew%n = k+l
        write(stateNew%label,'(i1,2A)') k+l, lchars(l+1:l+1), Mchars(abs(m)+1:abs(m)+1)
        cycle   ! Don't do any more checks or convert a sturmian into a state.
     end if

     ! Next, see if it can be replaced by a short-ranged "sub" function.
     j = nOld - nSub; found = .false.
     do while ( j<nOld .and. .not.found )
        j = j + 1
        sturmj => bst%b(j)
        if ( get_k(sturmj)==k .and. get_ang_mom(sturmj)==l .and. get_ang_mom_proj(sturmj)==m ) found = .true.
     end do

     ! Point to this "sub" function instead of the "main" one.
     if (found) then
        if (myid==0) write(*,'(3(A,I2),2(A,F5.2))') 'Laguerre function k =',k, ', l =',l, ', m =',m, ' alpha =',get_alpha(sturmi), ' --->', get_alpha(sturmj)
        nVec(:) = j
     else
        nVec(:) = nBst + i
     end if

     ! Make an orbital out of the Laguerre function.
     par = (-1)**l
     mVec(:) = m
     call construct_st(stateNew,hlike,dble(m),par,spin,E, k,ncm,CIVec,nVec,mVec)
     call set_n_majconf(stateNew, k+l)
     call set_l_majconf(stateNew, l)
     write(stateNew%label,'(i1,2A)') k+l, lchars(l+1:l+1), Mchars(abs(m)+1:abs(m)+1)
  end do
!  end do
!  end do
  
  call destruct_basis_st(statesTmp)


  if(sub_by_dom) then
    !at this point 'states' should be H2+ states plus full Laguerre basis - we just determine which Laguerres to remove

    allocate(ovlps(nNew))

    do i=1, nStates
      statei => states%b(i)
    
      ovlps = 0.0d0
      
      do j=nStates+1, nNew
        statej => states%b(j)
        ovlps(j) = 0.0d0
        do indi = 1, get_nam(statei)
           ni = get_na(statei,indi); CIi = get_CI(statei,indi)
           do indj = 1, get_nam(statej)
              nj = get_na(statej,indj); CIj = get_CI(statej,indj)
              ovlps(j) = ovlps(j) + CIi*CIj*bst%ortint(ni,nj)
           end do !indj
        end do !indi
      enddo !j

      n = maxloc(abs(ovlps),1) !index of state to remove
  
      call new_basis(statesTmp,nNew, hlike,basis_type)
      do j = 1, nNew
        call copy( statesTmp%b(j), states%b(j) )
      end do
      call destruct_basis_st(states)
      call new_basis(states,nNew-1, hlike,basis_type)
      do j=1, n-1
        call copy( states%b(j), statesTmp%b(j))
      enddo
      do j=n+1, nNew
        call copy( states%b(j-1), statesTmp%b(j))
      enddo
      nNew = nNew - 1
      call destruct_basis_st(statesTmp)

      statei => states%b(i)
      print*, 'FOR STATE '//statei%label//': removed Laguerre '//states%b(n)%label

    enddo !i
    
    do i=1, nStates
      statei => states%b(i)
    
      ovlps = 0.0d0
      
      do j=nStates+1, nNew
        statej => states%b(j)
        ovlps(j) = 0.0d0
        do indi = 1, get_nam(statei)
           ni = get_na(statei,indi); CIi = get_CI(statei,indi)
           do indj = 1, get_nam(statej)
              nj = get_na(statej,indj); CIj = get_CI(statej,indj)
              ovlps(j) = ovlps(j) + CIi*CIj*bst%ortint(ni,nj)
           end do !indj
        end do !indi
      enddo !j

      n = maxloc(abs(ovlps),1) !index of state to remove
  
      call new_basis(statesTmp,nNew, hlike,basis_type)
      do j = 1, nNew
        call copy( statesTmp%b(j), states%b(j) )
      end do
      call destruct_basis_st(states)
      call new_basis(states,nNew-1, hlike,basis_type)
      do j=1, n-1
        call copy( states%b(j), statesTmp%b(j))
      enddo
      do j=n+1, nNew
        call copy( states%b(j-1), statesTmp%b(j))
      enddo
      nNew = nNew - 1
      call destruct_basis_st(statesTmp)

      statei => states%b(i)
      print*, 'FOR STATE '//statei%label//': removed Laguerre '//states%b(n)%label

    enddo !i

  endif


end subroutine make_Hyl_basis
