!Module: one_electron_func_group
!Purpose: contains functions for setting up the one-electron 
!         diagonalisation of the target and including the results in
!         a basis for the two electron structure problem.

module one_electron_func_group



contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine: construct_1el_basis_nr_group
!Purpose: constructs basis for two electron structure as a mixture of
!         one-electron laguerre functions and molecular orbitals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine construct_1el_basis_nr_group(Number_one_electron_func,basis)
    use input_data       !data_in, dataMSC
    use grid_radial
    use sturmian_class
    !use vnc_module
    use one_electron_func  !Contains bst_nr variable
    use target_states    !TargetStates, TargetStates1el, TargetStates2el
		use basismodule
    use ovlpste1me       !ovlpst, e1me
    !use natural_orbitals_module
    use MPI_module       !myid

    implicit none
!
    integer:: ii, jj, kk
    !type(smallinput):: indata
		type(basis_sturmian_nr):: basis !Formerly bst_data
		!Indices for first laguerre basis
    integer, dimension(:), allocatable:: sturm_ind_list_one
		integer, dimension(:), allocatable:: k_list_one, l_list_one, m_list_one
		integer:: num_func, rad_func
		integer:: k,l,m
    real(dpf):: Rij, cosij
		logical:: found
		!Indices for second laguerre basis
    integer, dimension(:), allocatable:: sturm_ind_list_two
		integer, dimension(:), allocatable:: k_list_two, l_list_two, m_list_two
		integer:: num_func_two, rad_func_two
		integer:: basissize
		!Indices for combined basis
    integer, dimension(:), allocatable:: sturm_ind_list_nr
		integer, dimension(:), allocatable:: k_list_nr, l_list_nr, m_list_nr
		integer:: totfuncs

		integer:: lagnmax, basissize_nr
		integer:: hamlimmin, hamlimmax, maxL

    integer, intent(in):: Number_one_electron_func
    real*8, dimension(:,:), allocatable:: H, b, b2, CI
    real*8, dimension(:), allocatable:: w
    integer:: i, j, n, nd, N_oneel, jstart,jstop, ipar, ipar_n, l_n, li, lj, Nmax, lorb, itmp
    real*8:: al, Rd, Z1, Z2

    integer:: matz, irerr
    integer:: ma, ni, nj
    real*8:: tmp, res, tmp1, rval,  tmp_sign_CI
    real*8:: energy
    integer, dimension(:), allocatable:: no, mo
    logical:: hlike
    type(sturmian_nr), pointer:: pi, pj
    ! MSC Method Declarations
    type(basis_state):: tmp1elst
    type(basis_sturmian_nr):: bst_MSC, bst_data
    integer:: nBasis_data, nBasis_MSC, Nmax1el
    integer:: nst_Basis_MSC
    integer:: m_n, nst, itmpmax, nam, nspbst
    real*8:: rma
    real*8, dimension(:,:), allocatable:: CIigiven
    real*8, dimension(:), allocatable:: CInobst 
    integer, dimension(:,:), allocatable:: igiven
    integer, dimension(:), allocatable:: nobst, mobst, manst1el
    real*8, dimension(:,:,:), allocatable:: lag_ham1el_m
    integer:: max_latop
    real*8 :: norm

    basis_type = data_in%calculation_type
    hlike = .true.

		!------------------------------------Set up first Laguerre basis and perform 1e diagonalisation------------------------------!

    !Initialise sturmian data types, requires input data type defined in the MCCC code.

    call construct(basis, data_in)
		call countConfigs(basis, data_in, k_list_one, l_list_one, m_list_one, sturm_ind_list_one, num_func, rad_func)

    !Create space for basis of Mo's
    call new_basis_st(TargetStates,num_func,hlike,basis_type)

    if (myid == 0) write(*,'("allocated space for one-electron target states Nmax=",I5)') Number_one_electron_func
		!Store energy of ionised system, in this case, ionised H3++ (just three protons)
    if ( sum(data_in%Rvec(:)) .gt. 0.0_dpf ) then
       TargetStates%en_ion = 0.0_dpf
       !Loop over nuclei
       do ii=1, 3
          do jj = ii+1, 3
             !Use law of cosines to compute distance between nuclei
	           cosij = cos(data_in%thetavec(ii))*cos(data_in%thetavec(jj)) + &
	                   sin(data_in%thetavec(ii))*sin(data_in%thetavec(jj))*cos(data_in%phivec(ii)-data_in%phivec(jj))
             Rij = sqrt(data_in%Rvec(ii)**2 + data_in%Rvec(jj)**2 - &
	           2*data_in%Rvec(ii)*data_in%Rvec(jj)*cosij)
	           !Account for degenerate case where nuclei coincide
             if (Rij .gt. 0.0_dpf) then
                TargetStates%en_ion = TargetStates%en_ion + dble(data_in%charge(ii)*data_in%charge(jj))/Rij
             end if
          end do
       end do
    else ! Atomic Case
       TargetStates%en_ion = 0d0 
    end if      
          
    !Perform one-electron structure calculation and store in TargetStates 
		print*, "1e STRUCTURE"
    call one_electron_structure_group(TargetStates, basis, num_func, sturm_ind_list_one, k_list_one, l_list_one, m_list_one)

		print*, "SORTING DONE"
    call sort_by_energy_basis_st(TargetStates)


		!Spectro factors used to assign principal klm to a state, shown in state.core_parts
		print*, "SPECTRO FACTORS"
    call calc_spectro_factors_group(TargetStates, basis, data_in%labot, data_in%latop, sturm_ind_list_one, m_list_one)
    !if (myid == 0) then
    !   call print_energy_basis_st(TargetStates)
    !end if

    !Populate arary e1me to be used later in H12 diagonalization. Off diagonal calculated later.
    Nmax = basis_size_st(TargetStates)
    if(allocated(e1me)) deallocate(e1me)
    allocate(e1me(Nmax, Nmax))  
    allocate(ovlpst(Nmax, Nmax))  
    e1me(:,:) = 0d0
    ovlpst(:,:) = 0d0
    do j=1,Nmax          
       e1me(j,j) = get_energy_st(TargetStates%b(j))
       ovlpst(j,j) = 1d0       
    enddo

	  !----------------------------------- Construct Second Laguerre Basis and Merge -------------------------------------!
    !Construct second laguerre basis, called bst_MSC
    call construct(bst_MSC, dataMSC, dataMSC%alpha_nl)   

		!Construct second set of k_list, l_list and m_list variables 
		call countConfigs(bst_MSC, dataMSC, k_list_two, l_list_two, m_list_two, sturm_ind_list_two, num_func_two, rad_func_two)
         
    !Number of Sturmian functions in Laguerre basis and second laguerre basis
    nBasis_MSC = num_func_two     !basis_size(bst_MSC)
    nBasis_data = num_func        !basis_size(basis)
    
    ! Combine both laguerre bases, store in bst_nr
    ! Calculates <Lfp|Lf> for e-H2 exchange matrix elements.
    ! rearrange.f90 will carry through rearrangement
		print*, "COMBINE BASIS"
    call combine_basis_nr(bst_nr,bst_MSC,basis,basis_type)

		!Merge k_list, m_list, l_list, etc arrays for bst_nr
		totfuncs = nBasis_MSC + nBasis_data
		allocate(k_list_nr(totfuncs), l_list_nr(totfuncs), m_list_nr(totfuncs), sturm_ind_list_nr(totfuncs))

		do ii = 1, nBasis_MSC
			 k_list_nr(ii) = k_list_two(ii)
			 l_list_nr(ii) = l_list_two(ii)
			 m_list_nr(ii) = m_list_two(ii)
			 sturm_ind_list_nr(ii) = sturm_ind_list_two(ii)
	  end do
		do ii = 1, nBasis_data
			 jj = ii + nBasis_MSC
			 k_list_nr(jj) = k_list_one(ii)
			 l_list_nr(jj) = l_list_one(ii)
			 m_list_nr(jj) = m_list_one(ii)
			 !sturm_ind_list_nr(jj) = sturm_ind_list_one(ii) + basis%n
			 sturm_ind_list_nr(jj) = sturm_ind_list_one(ii) + bst_MSC%n
			 !print*, bst_MSC%n, sturm_ind_list_nr(jj), basis%n
	  end do

    !max_latop = get_max_L(bst_nr) ! Largest l of laguerre basis
		max_latop = 0      !Basis now includes m

		!lag_ham1el_m stores matrix elements for combined basis
		!lagnmax = basis_size(bst_nr)   !old size, with fixed m basis

		lagnmax = totfuncs
    if (allocated(lag_ham1el_m)) deallocate(lag_ham1el_m)
    allocate(lag_ham1el_m(totfuncs,totfuncs,-max_latop:max_latop))
    lag_ham1el_m(:,:,:) = 0d0
    
    if(allocated(e1me)) then
       deallocate(e1me) !
       deallocate(ovlpst) ! Basis functions overlaps
    end if
    
    !Count maximum number of Molecular states possible from second laguerre basis (2nd diagonalisation)
    itmp = 0
    do n = 1, nBasis_MSC
       !lorb = get_ang_mom(bst_MSC%b(n))
       itmp = itmp + 1   ! (2*lorb + 1) !Basis now includes m
    end do
    nst_Basis_MSC = itmp   

    allocate(e1me(nst_Basis_MSC, nst_Basis_MSC))   !Stores 1e matrix elements
    allocate(ovlpst(nst_Basis_MSC, nst_Basis_MSC)) !Overlap one- or two-electron functions/states 
    e1me(:,:) = 0d0
    ovlpst(:,:) = 0d0 
    
	  !Overwrite molecular orbitals with second laguerre basis and populate e1me and ovlpst
    !Hybrid calculates <Lfp|H|Lf> for e-H2 exchange matrix elements, which will carry through rearrangement      
		!lagnmax = basis_size(bst_nr)   !old size, with fixed m basis
		print*, "HYBRID MSC BASIS"
    call Hybrid_MSCbasis(bst_nr,TargetStates,nst_Basis_MSC,nBasis_MSC,lag_ham1el_m,lagnmax,max_latop,&
			                   m_list_nr, sturm_ind_list_nr)

    !Liam added: save the original description of the one-electron states before rearrange is called - needed for natural orbitals
    call copy(TargetStates_unrearranged,TargetStates)

    !Destruct temporary bases    
    call destruct(bst_MSC)
    call destruct(basis)

 	  !!Diagonalise 1e hamiltonian a second time in new sturmian basis, as CI coefficients are needed later
    !!Molecular states from second diagonalisation using second sturmian basis stored in TargetStates1el
    !call new_basis_st(TargetStates1el,num_func_two,hlike,basis_type)

	  !Save new states to TargetStates1el
    !call one_electron_structure_group(TargetStates1el, bst_MSC, indata, num_func_two, sturm_ind_list_two, k_list_two, l_list_two, m_list_two)

    !Call rearrange to represent 1e states in a one-electron basis with well defined (lm)
		if (.not. data_in%good_m) then
		   hamlimmin = 0
		   hamlimmax = 0
	  else
			 hamlimmin = -data_in%latop
			 hamlimmax = data_in%latop
	  end if
		maxL = get_max_L(bst_nr)
		basissize_nr = totfuncs  !lag_ham1el_m 

		print*, "REARRANGE"
    call rearrange(bst_nr,maxL,TargetStates,.true.,basissize_nr,hamlimmin,hamlimmax,lag_ham1el_m, totfuncs, m_list_nr, sturm_ind_list_nr)

!!$!------------------------------   

!!$! Use the same basis for H2+ 1s configuration and other configurations
   ! if (data_in%inc == 1 .AND. dataMSC%MSC_nconfig == 0) then
   !    if (myid == 0) then
   !       print*,'Calling rearrange without molecular orbital basis.'
   !       print*,"inc, MSC_nconfig", data_in%inc, dataMSC%MSC_nconfig
   !    end if

!!$! Calculating <Lfp|H|Lf> for e-H2 exchange matrix elements.
!!$! will carry through rearrangement  
   !    if (allocated(lag_ham1el_m))  deallocate(lag_ham1el_m)
   !    allocate(lag_ham1el_m(basis_size(bst_nr),basis_size(bst_nr),-data_in%latop:data_in%latop))
   !    lag_ham1el_m(:,:,:) = 0d0
   !    call Hybrid_H1el_st(Nmax,lag_ham1el_m,basis_size(bst_nr),get_max_L(bst_nr))
   ! 
   !    !Liam added: save the original description of the one-electron states before rearrange is called - needed for natural orbitals
   !    call copy(TargetStates_unrearranged,TargetStates)
   !    call rearrange(bst_nr,get_max_L(bst_nr),TargetStates,.TRUE.,lag_ham1el_m)
   !         
!!$! Use the different bases for H2+ 1sSg and 2pSu states and other configurations (2s, 3s, 3p, ..)       
   ! else if ( ABS(dataMSC%MSC_nconfig) /= 0  ) then  
   !    
   !    call destruct(bst_nr)
   !    
   !    call construct(bst_data, data_in)  ! Construct temporary basis from data.in inputs called bst_data
!  !     call construct(bst_MSC, dataMSC)   ! Construct temporary basis from dataMSC inputs called bst_MSC


   !    if (data_in%inc == 1) then !!$ inc must = 1 for two-electron targets. Checked in H12
   !       call rearrange(bst_nr,get_max_L(bst_nr),TargetStates,.TRUE.,lag_ham1el_m)
   !    end if
   ! 
!  !     print*
!  !     print*,'Write state 1 in a file'
!  !     call write_1elState(bst_nr, TargetStates%b(7))
!  !     print*
   ! 

!!$! Diagonalise H2 with 1sSg MO and AO 
   !    if ( dataMSC%MSC_nconfig >= 1 .AND. (.NOT. data_in%hlike) ) return

!!$! Above used for one and two-electron targets. Above can be used in H12.f90 now.      

   !       ! Molecular states from second diagonalisation using MSC basis stored in TargetStates1el
   !       call new_basis_st(TargetStates1el,Number_one_electron_func,hlike,basis_type)
   !       
!!$! Start second diagonalisation using MSC basis
   !       if(allocated(no) .OR. allocated (mo) ) then
   !          deallocate(no,mo)
   !       end if
   !       allocate(no(nBasis_MSC),mo(nBasis_MSC)) 
   !       no (:) = 0
   !       mo (:) = 0
   !       
   !       N_oneel = 0 
!!$! Partial Waves serperated by Angular Projection M and Parity 
   !       do ma = -data_in%Mt_max, data_in%Mt_max

   !          do ipar = 1, -1, -1

   !             if(data_in%good_parity) then
   !               if(ipar == 0) cycle
   !             elseif(.not.data_in%good_parity) then
   !               if(ipar /= 0) cycle
   !             endif
   !               
   !             nst = data_in%nst(ABS(ma),ipar)

   !             if(nst == 0) cycle

   !             if(nst < 0) nst = -1

   !             ! Determine number of configurations nd per partial wave
   !             i = 0
   !             do n = 1, nst_Basis_MSC
   !                m_n = get_ang_mom_proj( TargetStates%b(n) )
   !                ipar_n = get_par( TargetStates%b(n) )
   !                if ( m_n == ma .and. (ipar == ipar_n .or. .not.data_in%good_parity)) then
   !                   i = i + 1
   !                   no(i) = n ! Index to configuration for particular partial wave
   !                end if
   !             end do
   !             nd = i 

   !             if(nst == -1) nst = nd
   !             
   !             if(nd == 0) then
   !                if (myid == 0) print*,'Error: Molecular State Configuration. nd=0'
   !                if (myid == 0) write(*, '("Symmetry and size: ma, ipar, nd = ",2I3,I5)') ma, ipar, nd 
   !                stop
   !             endif
   !             
   !             mo(1:nd) = ma
   !             
   !             ! Temporary arrays
   !             allocate(H(nd,nd))
   !             allocate(b(nd,nd))
   !             allocate(CI(nd,nd))
   !             allocate(w(nd))
   !             
   !             H(:,:) = 0d0 
   !             b(:,:) = 0d0
!!$! Populate H = <f|H_T|i> and b = <f|i> using e1me and ovlpst.
   !             do ni = 1, nd
   !                i = no(ni)
   !                do nj = 1, ni
   !                   j = no(nj)
   !                   
   !                   H(ni,nj) = e1me(i,j)
   !                   b(ni,nj) = ovlpst(i,j)
   !                   
   !                   if ( Rd /= 0d0) then
   !                      H(ni,nj) = H(ni,nj) + Z1 * Z2 / Rd * b(ni,nj)
   !                   end if
   !                   
   !                   H(nj,ni) = H(ni,nj)
   !                   b(nj,ni) = b(ni,nj)
   !                   
   !                end do ! nj
   !             end do ! ni
   !       
   !             if(data_in%N_core_el > 0) H = H + core_energy*b
   !             print*, '>> add core_en', sum(H), sum(b)
   !             
   !             if (myid == 0) print*
   !             if (myid == 0) write(*, '("Symmetry and size: ma, ipar, nd = ",2I3,I5)') ma, ipar, nd   
   !             ! Diagonalise
   !             !matz=2
   !             !call  rsg(nd,nd,H,b,w,matz,CI,irerr)
   !             allocate(IFAIL(nd), IWORK(5*nd))
   !             allocate(WORK(1))
   !             LWORK = -1
   !             call dsygvx( 1, 'V', 'I', 'U', nd, H, nd, b, nd, 0.0d0,0.0d0, 1,nd, 2*DLAMCH('S'), &
   !               &NFOUND, w, CI, nd, WORK, LWORK, IWORK, IFAIL, irerr)
   !             
   !             LWORK = WORK(1)
   !             deallocate(WORK)
   !             allocate(WORK(LWORK))
   !             call dsygvx( 1, 'V', 'I', 'U', nd, H, nd, b, nd, 0.0d0,0.0d0, 1,nd, 2*DLAMCH('S'), &
   !               &NFOUND, w, CI, nd, WORK, LWORK, IWORK, IFAIL, irerr)
   !             deallocate(IFAIL, IWORK, WORK)
   !             
   !                
   !             if (myid == 0) then
   !                write(*,'("irerr =",I3)') irerr
   !                if(data_in%N_core_el > 0) then
   !                  print*, ' Energies in a.u. (including core energy)'
   !                else
   !                  print*, ' Energies in a.u.'
   !                endif
   !                write(*,'(5F15.5)') (real(w(i)), i=1,nd)
   !                print*
   !                if(data_in%N_core_el > 0) then
   !                  print*, ' Energies in eV (including core energy)'
   !                else
   !                  print*, ' Energies in eV'
   !                endif
   !                write(*,'(5F15.5)') (real(data_in%eV*w(i)), i=1,nd)
   !             end if
   !             
!!$! Construct states from MSC diagonalisation and store them in TargetStates1el        
   !             jstart = 1
   !             if(data_in%N_core_el > 0) jstart = 2
   !             jstop = min(nd,jstart+nst-1)
   !             do j = jstart, jstop
   !                energy = w(j)
   !                N_oneel = N_oneel + 1
   !                tmp_sign_CI = SUM(CI(1:nd,j))
   !                if (  tmp_sign_CI < 0d0 ) then
   !                   CI(1:nd,j) = -CI(1:nd,j)
   !                end if
   !                TargetStates1el%Nstates = N_oneel
   !                call  construct_st(TargetStates1el%b(N_oneel),hlike,dble(ma),ipar,0.5d0,energy,j,nd,CI(1:nd,j),no(1:nd),mo(1:nd))
   !             end do
   !             
   !             
   !             deallocate(H)
   !             deallocate(b)
   !             deallocate(CI)
   !             deallocate(w)
   !             
   !          end do ! parity
   !       end do ! ma

   !       
!!$! First diagonalisation stored in TargetStates used Laguerre Basis function configuratiosn from data_in
!!$! Second diagonalisation stored in TargetStates1el. Used Target States of First diagonalisation
!!$! Save CI coefficients from Second diagonalisation in CIigiven. 
   !       
   !       Nmax1el = TargetStates1el%Nstates ! Second Diagonalisation number of states
   !       Nmax1el = min(Nmax1el,Number_one_electron_func) ! Number of states requested
   !       TargetStates1el%Nmax = Nmax1el
   !       
!!$! nBasis_MSC = Second Diagonalisation Laguerre basis size. Each target state will be constructed from MO/target states for that symmetry.
!!$! CIigiven = CI coefficients to configurations per state
!!$! igiven_max = Number of one-electron configurations per state
!!$! igiven =  Index to one-electron target state configurations
!!$! manst1el = Angular projection of configuration
   !       allocate( CIigiven(Nmax1el,nBasis_MSC), igiven(Nmax1el,nBasis_MSC), manst1el(Nmax1el))  
!!$! igiven_max(Nmax1el)    igiven_max(:) = 0
   !       igiven(:,:) = 0
   !       CIigiven(:,:) = 0d0         
   !       manst1el(:) = 0
   !       
   !       do nst = 1, Nmax1el
   !          manst1el(nst) = get_ang_mom_proj(TargetStates1el%b(nst)) ! Assign angular projection
   !          nam = get_nam(TargetStates1el%b(nst)) ! Number of one-electron configurations per state
   !          !                igiven_max(nst) = nam
   !          do j = 1, nam
   !             igiven(nst,j) = get_na(TargetStates1el%b(nst),j)
   !             CIigiven(nst,j) = get_CI(TargetStates1el%b(nst),j)
   !          end do
   !       end do
   !       
!!$! nsbst = size of combined basis (bst_nr). nobst pointer to all basis functions which make state
!!$! nobst = index to combined Laguerre Basis or rearranged Basis functions
   !       
   !       nspbst = basis_size(bst_nr) 
   !       allocate(nobst(nspbst),mobst(nspbst),CInobst(nspbst))
   !       nobst(:) = 0 
   !       mobst(:) = 0 
   !       CInobst(:) = 0d0
   !       
!!$! Multiply CI coefficients from first and second diagonalisation. This allows us to use the rearrange routine.
!!$! After combining CI coefficients construct new states.
   !       call new_basis_st(tmp1elst,Nmax1el,hlike,basis_type)            
   !       do nst = 1, Nmax1el
   !          nd = get_nam(TargetStates1el%b(nst)) ! Size of second diagonalisation
   !          
   !          call make_nomo(Targetstates,bst_nr,nd,CIigiven(nst,1:nd),igiven(nst,1:nd),manst1el(nst),nspbst,nobst,mobst,CInobst,itmpmax)               

   !          ipar = get_par(TargetStates1el%b(nst))
   !          rma = manst1el(nst)
   !          energy = get_energy_st(TargetStates1el%b(nst))
   !          nd = itmpmax
   !          tmp1elst%Nstates = nst

   !          call  construct_st(tmp1elst%b(nst),.true.,rma,ipar,0.5d0,energy,get_inum_st(TargetStates1el%b(nst)),nd,CInobst(1:nd),nobst(1:nd),mobst(1:nd))
   !       end do


!!$! Replace molecular orbitals with short-ranged functions
!!$! e1me and ovlpst are calculated here for H2+ and H2
!  !           call replace_molecular_orbitals(tmp1elst)
   !       
   !       if (data_in%inc == 1) then
   !          call rearrange(bst_nr,get_max_L(bst_nr),tmp1elst,.FALSE.)
   !       end if

   !       call destruct_basis_st(TargetStates)
   !       call new_basis_st(TargetStates,tmp1elst%Nstates,.true.,0)
   !       TargetStates%Nmax = tmp1elst%Nmax
   !       TargetStates%Nstates = tmp1elst%Nstates
   !       if ( Rd /= 0d0 ) then
   !          TargetStates%en_ion = Z1 * Z2 / Rd
   !       else ! Atomic Case
   !          TargetStates%en_ion = 0d0
   !       end if

   !       if(data_in%N_core_el > 0) then
   !         TargetStates%en_ion = TargetStates%en_ion + core_energy !core_energy in vnc_module - does not contain nuclear repulsion
   !       endif

   !       do nst = 1, tmp1elst%Nstates
   !          call copy_st(TargetStates%b(nst), tmp1elst%b(nst))
   !       end do
   !       call destruct_basis_st(tmp1elst)

   !       call sort_by_energy_basis_st(TargetStates)
!!$! calc_spectro_factors sets target states: l_major, n - principal quantum
!nu!mber and correspond Laguerre basis k stored in inum
   !       call calc_spectro_factors(TargetStates, bst_nr, min(dataMSC%labot,data_in%labot), max(dataMSC%latop,data_in%latop))
   !       if (myid == 0) then
   !          call print_energy_basis_st(TargetStates)
   !       end if
   !    
   ! else
   ! 
   !    if (myid==0) print*,"rearrange and Molecular State Configuration basis were not called"

   ! end if ! inc rearrange and MSC

		deallocate(k_list_one,l_list_one,m_list_one,sturm_ind_list_one)
		deallocate(k_list_two,l_list_two,m_list_two,sturm_ind_list_two)


    
end subroutine construct_1el_basis_nr_group





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine: one_electron_structure_group
!Purpose: performs a one-electron structure calculation for in general
!         a non-linear molecule, using the input laguerre basis. Possibly
!         subject to symmetry restrictions.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine one_electron_structure_group(oneestatebasis, basis, num_func, sturm_ind_list,k_list, l_list, m_list)  
    use grid_radial
    use sturmian_class
		use input_data
    use state_class
    use basismodule
    use numbers
    use ieee_arithmetic
		use vnc_module
    implicit none
		!Nuclear coordinates and additional options
    !type(smallinput):: indata
		!Given sturmian basis
    type(basis_sturmian_nr)::basis
    !State basis type for storing one electron states
    type(basis_state):: oneestatebasis
    integer:: nr  
		!Matrix elements, energies and CI coeficients to be calculated.
    integer:: ii, jj, kk, n, m
    real(dpf), dimension(:), allocatable:: w
    real(dpf), dimension(:,:), allocatable:: z
    complex(dpf), dimension(:,:), allocatable:: V, H, KMat, B, VPot, VPotTemp
    logical, dimension(:), allocatable:: use_list
    complex(dpf), dimension(:,:,:), allocatable:: VRadMatEl
	  !Basis function index arrays and matrix elements
    real(dpf), dimension(:,:), allocatable:: realH, realB, realK, BMat
    real(dpf), dimension(:,:,:), allocatable:: angular
    integer, dimension(:):: k_list, m_list, l_list, sturm_ind_list
    real(dpf), dimension(:), allocatable::  energies
    integer:: num_func     !, l, m, k, rad_func
		integer:: num_lambda
		real(dpf):: largest, largestZ
    !Arrays for testing case
    !integer, dimension(20):: testm, testl, testk
		!For writing one electron states to TargetStates
    integer:: numfound, u1, u2
    complex(dpf), dimension(:,:), allocatable:: tempK, tempB
    integer:: inum, ncm, counter
    integer, dimension(:), allocatable:: no1, no2, mo1, mo2, phase
		integer:: mmaj, largestind
	  !For renormalising CI coeffs
    real(dpf):: norm
		!Write wave functions to file
    real(dpf), dimension(:), allocatable:: func
    character(len=40)::filename
    logical:: writeWaveFunc = .false.
    !Data required for calling lapack routine dsygv
    real(dpf), dimension(:), allocatable:: work
    integer:: lda, ldb
    integer:: lwork, info
    logical:: uselapack
		!For calling rsg instead
		integer:: i1, i2, ier
		!For calculating nuclear-nuclear interaction
		real(dpf):: Rij, cosij
		!For writing energies to file
		integer:: nstates
 
    nr = grid%nr
		uselapack = .true.

		!Now set up in main and stored in vnc_module's "vnc" variable
    !!Set up expansion of potential, use formula: SUM_l=0^l=L (2l+1) = (L+1)^2
    !num_lambda = (indata%lambdamax+1)**2
 
    !allocate(VPot(nr,num_lambda))
    !VPot(:,:) = 0.0_dpf
    !allocate(VPotTemp(nr, num_lambda))
    !do ii = 1, 3
    !   call getVPotNuc(grid, VPotTemp, indata%R(ii), indata%theta(ii), &
 	  !            indata%phi(ii), indata%charge(ii), indata)
    !   VPot(:,:) = VPot(:,:) + VPotTemp(:,:)
    !end do
    !deallocate(VPotTemp)

    allocate(VPot(size(vnc(:,1)),size(vnc(1,:))))
		VPot(:,:) = vnc(:,:)	
    
    allocate(H(num_func,num_func),V(num_func,num_func))
    H(:,:) = 0.0_dpf
    V(:,:) = 0.0_dpf
 
    allocate(realB(num_func,num_func))
		realB(:,:) = 0.0_dpf
    do ii = 1, num_func
			 do jj = 1, num_func
					if (m_list(ii) .eq. m_list(jj)) then
						 n = sturm_ind_list(ii)
						 m = sturm_ind_list(jj)
             realB(ii,jj) = basis%ortint(n,m)
				  end if
			 end do
	  end do
    call getKMatM(realK,realB,basis,num_func,sturm_ind_list)   !Assumes B loops over klm
    allocate(B(num_func,num_func), KMat(num_func,num_func))
    B(:,:) = realB(:,:)
    KMat(:,:) = realK(:,:)
    deallocate(realB, realK)
 
    !!Defined arrays for test case
    !testm(:) = indata%testm(:) 
    !testl(:) = indata%testl(:)
    !testk(:) = indata%testk(:) 

		!!Set up restricted basis for test case with special small basis
    allocate(use_list(num_func))
    use_list(:) = .true.
    !if ((indata%isobasisop .eq. 1) .and. (indata%isoscop .eq. 1)) then
    !   use_list(:) = .false.
    !   numfound = 0
    !   do ii = 1, 20
    !      print*, testk(ii), testl(ii), testm(ii)
    !      do jj = 1, num_func
    !   	    if (((k_list(jj) .eq. testk(ii)) .and. &
    !   	    (l_list(jj) .eq. testl(ii)))  .and. &
    !         (m_list(jj) .eq. testm(ii)))  then
    !   
    !   	       use_list(jj) = .true.
    !   	       numfound = numfound + 1
    !         end if
    !   	 end do
    !   end do
    !   if (numfound .ne. 20) then
    !   	 print*, "ERROR: basis functions with required symmetries for &
    !   	 &test case not found, require at least lmax=3 and N=5. &
    !   	 & Stopping."
    !   	 error stop
    !   end if

    !   !Produce K and B matrix for restricted basis, just copy over the
    !   !right elements from the full K and B matrices
    !   allocate(tempK(num_func,num_func), tempB(num_func, num_func))
    !   tempK(:,:) = KMat(:,:)
    !   tempB(:,:) = B(:,:)
    !   deallocate(KMat, B)
    !   allocate(KMat(20,20), B(20,20))
    !   u1 = 1
    !   do ii = 1, num_func
    !      u2 = 1
    !      do jj= 1, num_func
    !         if (use_list(ii) .and. use_list(jj)) then
    !            KMat(u1,u2) = tempK(ii,jj)
    !            B(u1,u2) = tempB(ii,jj)
    !            u2 = u2 + 1
    !         end if
    !      end do
    !      if (use_list(ii)) then
    !         u1 = u1 + 1
    !      end if
    !   end do
    !   deallocate(tempK, tempB, H, V)	
	  ! 		    
    !   num_func = 20
    !   allocate(H(num_func,num_func), V(num_func,num_func))
    !   H(:,:) = 0.0_dpf
    !   V(:,:) = 0.0_dpf
    !end if

    print*, "GET ANGULAR MATRIX ELEMENTS"
		num_lambda = (data_in%ltmax+1)**2
    !!Precalculate angular integrals appearing in V-matrix elements
    allocate(angular(num_lambda,num_func,num_func))
    call getAngular(num_func, angular, l_list, m_list)
 
    !Precalculate radial matrix elements
    allocate(VRadMatEl(num_lambda,basis%n,basis%n))
    VRadMatEl(:,:,:) = 0.0_dpf
    print*, "GET RADIAL MATRIX ELEMENTS"
    call getRadMatEl(basis, VPot, grid, VRadMatEl)
 
    !If real spherical harmonics are used, V matrix elements will have complex part zero.
    call getVMatEl(sturm_ind_list, V, VRadMatEl, num_func, angular, use_list)
    print*, "V Matrix Elements Computed"

 
    deallocate(angular)
    deallocate(VRadMatEl)
    H(:,:) = KMat(:,:) + V(:,:)
 
    allocate(realH(num_func,num_func), realB(num_func,num_func))
    realH(:,:) = real(H(:,:))
    realB(:,:) = real(B(:,:)) 

    !Nuclear-nuclear matrix elements, being multiples of the overlap matrix, do not affect electronic wave functions.
    print*, "ADD NUCLEAR INTERACTION ENERGY"
    do ii=1, 3
       do jj = ii+1, 3
          !Use law of cosines to compute distance between nuclei
 	       cosij = cos(data_in%thetavec(ii))*cos(data_in%thetavec(jj)) + &
 	               sin(data_in%thetavec(ii))*sin(data_in%thetavec(jj))*cos(data_in%phivec(ii)-data_in%phivec(jj))
          Rij = sqrt(data_in%Rvec(ii)**2 + data_in%Rvec(jj)**2 - &
 	       2*data_in%Rvec(ii)*data_in%Rvec(jj)*cosij)
 	       !Account for degenerate case where nuclei coincide
          if (Rij .gt. 0.0_dpf) then
             realH(:,:) = realH(:,:) + realB(:,:)*dble(data_in%charge(ii)*data_in%charge(jj))/Rij
          end if
       end do
    end do
    print*, "NUCLEAR INTERACTION CALCULATED"
 

    !Check for infinite values to be safe. Sometimes results from use of older versions of Yint 
	  !subroutine within and openMp region, as older versions were not thread safe.
    largest = 0.0_dpf
    do ii = 1, num_func
       do jj = 1, num_func
          if (isnan(realH(ii,jj))) then
             print*, "realH: ", ii,jj
             print*, realH(ii,jj)
             error stop
          end if
          if (isnan(realB(ii,jj))) then
             print*, "realB: ", ii,jj
             error stop
          end if
          if (.not. ieee_is_finite(realH(ii,jj))) then
             print*, "realH infinite: ", ii,jj
             print*, H(ii,jj)
             print*, realH(ii,jj)
             error stop
          end if
          if (.not. ieee_is_finite(realB(ii,jj))) then
             print*, "realB infinite: ", ii,jj
             error stop
          end if
    
          if ((abs(realH(ii,jj)) .gt. largest) .and. (abs(realH(ii,jj)) .gt. 0.0_dpf)) then
             largest = realH(ii,jj)
          end if
       end do
    end do
    do ii = 1, num_func
       do jj = 1, num_func
          if (abs(realH(ii,jj)) .lt. 1E-8*abs(largest)) then
             realH(ii,jj) = 0.0_dpf
          end if
       end do
    end do
    print*, "Removed small matrix elements"

    allocate(w(num_func),z(num_func,num_func))

    allocate(BMat(num_func,num_func))
    BMat(:,:) = realB(:,:)

    !Diagonalise one-electron hamiltonian
    if (uselapack .eqv. .false.) then
       print*, "CALL RSG"
       call rsg(num_func,num_func,realH,realB,w,1,z,ier)
    else
       print*, "CALL DSYGV"
       lda = num_func
       ldb = num_func
       allocate(work(1))
       call DSYGV(1, 'V', 'U', num_func, realH, lda, realB, ldb, w, work, -1, info) !-1 -> workspace query, get best lwork
       lwork = int(work(1))
       deallocate(work)
       allocate(work(lwork))
       call DSYGV(1, 'V', 'U', num_func, realH, lda, realB, ldb, w, work, lwork, info) 
       z(:,:) = realH(:,:) !On exit, dsygv overwrites matrix with eigenvectors
       print*, "1e diagonalisation, info= ",  info
       if (info .ne. 0) then
          print*, "ERROR in DSYGV, info !=0"
          error stop
       end if

       deallocate(work)
    end if


    !Find largest expansion coefficient, ignore those below a certain
    !magnitude for stability/to get rid of underflow 
	  print*, "REMOVING SMALL CI COEFFICIENTS"
    largestZ = 0.0_dpf
    do ii = 1, num_func
       do jj = 1, num_func
          if (abs(z(ii,jj)) .gt. abs(largestZ)) then
             largestZ = z(ii,jj)
          end if
       end do
    end do
    do ii = 1, num_func
       do jj = 1, num_func
          if ((abs(z(ii,jj)) .lt. 1E-10*abs(largestZ)) .and. (abs(z(ii,jj)) .gt. 0.0_dpf)) then
             z(ii,jj) = 0.0_dpf
          end if
       end do
    end do

	  print*, "RENORMALISING CI COEFFICIENTS"
    !Renormalise coefficients so norm is 1
	  !Added openMP to speed this up for large basis sizes
    !$OMP PARALLEL DO DEFAULT(SHARED) & 
    !$OMP PRIVATE(ii, jj, norm) &
    !$OMP SCHEDULE(DYNAMIC)
    do n = 1, num_func         
       norm = 0.0_dpf
       do ii=1, num_func
          do jj = 1, num_func
             norm = norm + z(ii,n)*z(jj,n)*BMat(ii,jj)
          end do
       end do
       z(:,n) = z(:,n)/sqrt(norm)
    end do
	  !$OMP END PARALLEL DO

	  print*, "FIXING SIGN ISSUE"
    !Fix CI coefficient sign issue
    do ii = 1, num_func
       if (sum(z(:,ii)) .lt. 0.0_dpf) then
          z(:,ii) = -z(:,ii)
       end if
    end do

    deallocate(realH,realB,BMat)

    !Create basis of one-electron states, store in state_basis data type
    call new_basis_st(oneestatebasis,num_func, .true. , 0)
    do ii =1, oneestatebasis%Nmax
       inum = ii        !Index of the state for given symmetry
       ncm = num_func   !Number of CI coefficients in expansion
       allocate(no1(ncm), no2(ncm), mo1(ncm), mo2(ncm), phase(ncm))
       !no1, mo1 store superindex ii and magnetic number m of basis funcs
       !used in expansion of the state.
       !Useful if we use only a subset of the full basis for symmetry reasons.
       do jj =1, num_func
          no1(jj)=jj
          mo1(jj)=m_list(jj)
       end do
       no2(:)=no1(:)
       mo2(:)=mo1(:)
       phase(:)=1.0_dpf   !Phase of CI coefficients           
			 !Find largest CI coefficient
       largestZ = 0.0_dpf
       do jj = 1, num_func
          if (abs(z(jj,ii)) .gt. abs(largestZ)) then
             largestZ = z(jj,ii)
			       largestind=jj
          end if
       end do
			 mmaj = m_list(largestind)
			 !print*, ii, "k,l,M: ", k_list(largestind), l_list(largestind), mmaj
       call construct_st(oneestatebasis%b(ii),.true.,dble(mmaj),0,0.5_dpf,w(ii),inum,num_func,z(:,ii),no1,mo1,no2,mo2,phase)
       call set_l_majconf(oneestatebasis%b(ii), l_list(largestind))
			 call set_n_majconf(oneestatebasis%b(ii), k_list(largestind)+l_list(largestind))  !n=k+l
       deallocate(no1,no2,mo1,mo2,phase)
    end do

    if((.not. (sum(data_in%Rvec(:)) .gt. 0.0_dpf)) .and. (writeWaveFunc)) then
       write(filename,'(A13,I0,A4)') 'oneeradfuncN=', num_func, ".txt"
       filename=TRIM(filename)
       open(77,file=filename) 
       write(filename,'(A8,I0,A4)') 'oneeCIN=', num_func, ".txt"
       filename=TRIM(filename)
       open(80,file=filename)

       write(77,*) "One electron ground state radial wave function, hydrogen like case"
       write(77,*) "r (a.u)   psi(r)"

       allocate(func(grid%nr))
       func(:) = 0.0_dpf
       do ii = 1, num_func
          i1 = basis%b(ii)%minf
          i2 = basis%b(ii)%maxf
          func(i1:i2) = func(i1:i2) + z(ii,1)*basis%b(ii)%f(i1:i2)
       end do

       do ii = 1, grid%nr
          write(77,*) grid%gridr(ii), func(ii)
       end do

       do ii = 1, num_func
          write(80,*) ii, z(ii,1)
       end do

       deallocate(func)
       close(77)
       close(80)
    end if

    !Number of energies calculated will be equal to the size of the basis used.
    !Allocate array to store energy of each state
    nstates = num_func
    allocate(energies(nstates))
    energies(:) = 0.0_dpf
    nstates=0
    nstates=num_func
    energies(:) = w(:)

    !deallocate(k_list, l_list, m_list)
    deallocate(H,B,V,KMat)
    deallocate(use_list)
    deallocate(w,z)    

    print*, "WRITE ENERGIES TO FILE"
    !Write energes of states to file
    open(80,file="1eenergies.txt")
    write(80,*) "State Energies (Ha)"
    write(80,*) "R1=", data_in%Rvec(1), " R2=", data_in%Rvec(2), "R3=", data_in%Rvec(3)
    write(80,*) "N (index), E(Ha)"
    do ii = 1, min(nstates, 100)
       write(80,*) ii, energies(ii)
    end do
    close(80)

 end subroutine one_electron_structure_group


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !Subroutine: Hybrid_MSCbasis
 !Purpose: merges the one-electron molecular orbitals and the second laguerre basis
 !         into a single basis for use in the two electron-structure calculation
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine Hybrid_MSCbasis(bst,TargetStates,Nmax,NBasisTwo,lag_ham1el_m,Lag_nmax,maxL,&
			                      m_list_nr, sturm_ind_list_nr)
    use input_data
		use basismodule
    use sturmian_class
    use state_class
    use ovlpste1me
    use MPI_module
    use vnc_module
    
    implicit none

		!type(smallinput):: indata
    type(basis_sturmian_nr), intent(in) :: bst   ! this is Sturmian basis 
    type(basis_state), intent(inout) :: TargetStates
    integer, intent(inout):: Nmax, NBasisTwo
    integer, intent(in):: Lag_nmax, maxL
    real*8, dimension(Lag_nmax,Lag_nmax,-maxL:maxL), intent(inout):: lag_ham1el_m
 
    type(basis_state):: IonStates !
    type(state), pointer :: state1, state2, MolOrb
    integer::  basis_type, n, nd, lst, ma, ipar, k, nst1, nst2, Ntmp, ist, m1,m2
    logical:: hlike
    real*8, dimension(1):: CI
    integer, dimension(1):: no, mo
    real*8:: energy
    real*8:: a, b
    integer :: kMO, lMO, mMO, j
    character(len=:), allocatable :: labelMO
    logical :: replaced_with_MO

		real*8:: cosij, Rij
		integer:: ii, jj
	  integer, dimension(:):: m_list_nr, sturm_ind_list_nr
		integer:: mmin, mmax, ind

    if ( dataMSC%MSC_nconfig == 0 ) stop ! Should not be here!!!
    ! Number of Sturmians in second basis
    if (myid == 0) print*, 'Size of MSC bst basis: ',NBasisTwo

    if(dataMSC%MSC_nconfig > TargetStates%Nmax) then
      if(myid==0) write(*,*) '*** ERROR: requested more 1-el states in the second basis than kept from the&
        & first diagonalisation'
      error stop
    endif


    !Record ground state and vary index to data.in Laguerre functions
    basis_type = data_in%calculation_type
    hlike = .true.
    call new_basis_st(IonStates,ABS(dataMSC%MSC_nconfig),hlike,basis_type)
    do n = 1, ABS(dataMSC%MSC_nconfig)
       !call copy_H2plus_st(IonStates%b(n),TargetStates%b(n),NBasis_F5)
       call copy_H2plus_st(IonStates%b(n),TargetStates%b(n),NBasisTwo)       !<<----- shifts na array of MO's by NBasisTwo to match new bst_nr basis 
    end do
    call destruct_basis_st(TargetStates)

    !!$ create target state basis of size Nmax
    call new_basis_st(TargetStates,Nmax,hlike,basis_type)
    if ( sum(data_in%Rvec(:)) .gt. 0.0_dpf ) then
       TargetStates%en_ion = 0.0_dpf
       !Loop over nuclei
       do ii=1, 3
          do jj = ii+1, 3
             !Use law of cosines to compute distance between nuclei
	           cosij = cos(data_in%thetavec(ii))*cos(data_in%thetavec(jj)) + &
	                   sin(data_in%thetavec(ii))*sin(data_in%thetavec(jj))*cos(data_in%phivec(ii)-data_in%phivec(jj))
             Rij = sqrt(data_in%Rvec(ii)**2 + data_in%Rvec(jj)**2 - &
	           2*data_in%Rvec(ii)*data_in%Rvec(jj)*cosij)
	           !Account for degenerate case where nuclei coincide
             if (Rij .gt. 0.0_dpf) then
                TargetStates%en_ion = TargetStates%en_ion + dble(data_in%charge(ii)*data_in%charge(jj))/Rij
             end if
          end do
       end do
    else ! Atomic Case
       TargetStates%en_ion = 0d0 
    end if      

    if(data_in%N_core_el > 0) then
      TargetStates%en_ion = TargetStates%en_ion + core_energy !core_energy in vnc_module - does not contain nuclear repulsion
    endif
  
    !These parameters used for Sturmians being stored in State objects
    energy = 0d0
    CI(1) = 1d0
    nd = 1
    
    ist = 0
    do ii=1, NBasisTwo
			 n = sturm_ind_list_nr(ii)
       k = get_k(bst%b(n))
       lst = get_ang_mom(bst%b(n))
       ipar = (-1)**lst

			 mmin = -lst
			 mmax = lst
			 if (.not. data_in%good_m) then
					mmin = 1
					mmax = mmin
			 end if

       do ind = mmin, mmax
					if (data_in%good_m) then
						 ma = ind
				  else 
						 ma = m_list_nr(ii)
				  end if
          ist = ist + 1
          TargetStates%Nstates = ist
          mo(1) = ma
					if (data_in%good_m) then
             no(1) = n
				  else
						 !No is index of basis function, not equal to sturmian index when m is included in basis
						 no(1) = ii
				  end if
    
          replaced_with_MO = .false.
          !Liam: added this to allow including an arbitrary number of MOs in the hybrid basis
					!Diatomic: Loop through MO's to find one with principal contribution (kl) matching sturmian (kl), subject to m symmetry
					!Non-linear: Loop through MO's find one with principal (klm) matching sturmian (klm), subject to group symmetry.
          do j=1, dataMSC%MSC_nconfig
             MolOrb => IonStates%b(j)
             lMO = get_l_majconf(MolOrb)
             kMO = get_n_majconf(MolOrb) - lMO 
					   mMO = get_ang_mom_proj(MolOrb)    !Gives mmajconf in non-diatomic case
             labelMO = trim(adjustl(get_label(MolOrb)))

					 	!Require sturmian replaced to be assigned carefully to avoid overcompleteness
					 	!Scheme: replace sturmians with orbitals with principal parts same as the sturmians' (klm)
             if((kMO == k .and. lMO == lst .and. mMO == ma)) then
                call copy_st(TargetStates%b(ist),MolOrb)
                replaced_with_MO = .true.
                if(myid == 0) write(*, &
									'(" Replaced Sturmian k=",I0,", l=",I0,", m=",I0," with one-electron state ",I0," ( ",A," )")') &
                  k, lst, ma, j, labelMO
             endif
          enddo !j
             
          if(.not.replaced_with_MO) then  !This state will be a Sturmian
             call construct_st(TargetStates%b(ist),hlike,dble(ma),ipar,0.5d0,energy,k,nd,CI,no,mo)
             call set_n_majconf(TargetStates%b(ist),k+lst)
             call set_l_majconf(TargetStates%b(ist),lst)
             call set_label(TargetStates%b(ist),k+lst,lst,ma)
          endif
          
         !This is the old code for just adding the 2pSu state 
         !
         ! if ( ABS(dataMSC%MSC_nconfig) == 2 .AND. ma == 0 .AND. k == 1 .AND. lst == 1 ) then
         !   call copy_st(TargetStates%b(ist),IonStates%b(2))
         ! else if ( ABS(dataMSC%MSC_nconfig) <= 2) then 
         !   call construct_st(TargetStates%b(ist),hlike,dble(ma),ipar,0.5d0,energy,k,nd,CI,no,mo)
         !   call set_n_majconf(TargetStates%b(ist),k)
         !   call set_l_majconf(TargetStates%b(ist),lst)
         ! else
         !   print*," Have not coded for lagrer MSC basis. Only 1s and 2p are done"
         !   stop
         ! end if
       end do

    end do 

    !For testing
		!do ii = 1, TargetStates%Nstates
		!	 print*, TargetStates%b(ii)%n, TargetStates%b(ii)%l, TargetStates%b(ii)%M, TargetStates%b(ii)%label
	  !end do
		!stop

    call destruct_basis_st(IonStates)

    if ( ist /= Nmax ) then
       if (myid == 0) print*, 'one_electron_func.f90: Hybrid_MSCbasis(): ist .ne. Nmax :', ist, Nmax
       stop
    end if

    !!$ Calculates e1me and lag_ham1el_m
		!NOTE: must pass maxL argument in order for array indexing to be correct, maxL changes in H3+ mode
    call Hybrid_H1el_st_group(Nmax,lag_ham1el_m, Lag_nmax, maxL, m_list_nr, sturm_ind_list_nr)
    !call Hybrid_H1el_st_group(Nmax,lag_ham1el_m, Lag_nmax, get_max_L(bst), m_list_nr, sturm_ind_list_nr)

    do nst1 = 1, Nmax
       state1 => TargetStates%b(nst1)
       m1 = get_ang_mom_proj(state1)

       do nst2 = 1, nst1
          state2 => TargetStates%b(nst2)
          m2 = get_ang_mom_proj(state2)

          if ( (data_in%good_m .and. m1 /= m2) .OR. (data_in%good_parity .and. get_par(state1) /= get_par(state2)) ) cycle

          ovlpst(nst1,nst2) = ovlp_st_group(state1,state2,bst,m_list_nr,sturm_ind_list_nr)
          ovlpst(nst2,nst1) =  ovlpst(nst1,nst2)

          !if (myid == 0) print*, nst1, nst2,  e1me(nst1,nst2), ovlpst(nst1,nst2)

       end do
    end do

!    a =  H1el_st(TargetStates%b(2),TargetStates%b(3))    
!    call Hlagorb(bst,2,3,0,b)
!    if (myid == 0) print*, '>>> a,b:',a,b
!
!    a = bst%ortint(2,3)
!    b =  ovlp_st(TargetStates%b(2),TargetStates%b(3))
!    if (myid == 0) print*, '!!! a,b:',a,b
!

  end subroutine Hybrid_MSCbasis

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Subroutine: Hybrid_H1el_st_group
  !Purpose: calculates one-electron hamiltonian matrix elements and overlap
  !         matrix in hybrid basis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Hybrid_H1el_st_group(Nmax,lag_ham1el_m,nmax_lag, max_latop, m_list_nr, sturm_ind_list_nr)
  !!$ Calculates the one-electron molecular Hamiltonian matrix elements !$ lag_ham1el_m = <Lp|H|L> for Laguerre basis functions for particular m.
  !!$ e1me = <np|H|n> for target states with underlying laguerre basis functions 
  
      use ovlpste1me
      use target_states
      use MPI_module
      use state_class
      use input_data
			use basismodule
  
      implicit none
      integer, intent(in):: Nmax ! Number of states      
      integer, intent(in):: nmax_lag, max_latop ! Number of Lag function and max l
      real*8, dimension(nmax_lag,nmax_lag,-max_latop:max_latop), intent(inout):: lag_ham1el_m
  
      type(state), pointer:: state_j, state_i
      real*8:: Ci, Cj, Ham
      real*8:: H1el_Lag
      integer:: nsti, nstj, i, j, n, ni, nj             ! Indexes
      integer:: ma

			integer:: mi, mj, sturmi, sturmj
			integer, dimension(:):: m_list_nr   !Combined m_list for bst_nr = (bst_MSC, bst_data)
			integer, dimension(:):: sturm_ind_list_nr   !Combined m_list for bst_nr = (bst_MSC, bst_data)



      !$OMP PARALLEL DO DEFAULT(SHARED) & 
      !$OMP PRIVATE(nstj,i,j,ma,Ham,ni,Ci,nj,Cj,mi,mj,sturmi,sturmj) &
      !$OMP SCHEDULE(DYNAMIC)
      do nsti = 1, Nmax
         !state_i => TargetStates%b(nsti)
         !ma = state_i%m
				 ma = TargetStates%b(nsti)%m
				 if (.not. data_in%good_m) then
						ma = 0   !Basis contains m
			   end if 
  
         do nstj = 1, nsti
            !state_j => TargetStates%b(nstj)
  
            Ham = 0d0
  
            !if (data_in%good_m .and. (state_i%m /= state_j%m)) cycle
            !if (data_in%good_parity .and. (state_i%parity /= state_j%parity)) cycle
            if (data_in%good_m .and. (TargetStates%b(nsti)%m /= TargetStates%b(nstj)%m)) cycle
            if (data_in%good_parity .and. (TargetStates%b(nsti)%parity /= TargetStates%b(nstj)%parity)) cycle

            do i = 1, TargetStates%b(nsti)%nam
							 !If state is an MO, ni refers to combined bst_nr, while ci refers original basis with sturmians stored in bst_data
               !ni = state_i%na(i)
               !Ci = get_CI(state_i,i)
							 ni = TargetStates%b(nsti)%na(i)
							 Ci = TargetStates%b(nsti)%CI(i)
  
               do j = 1, TargetStates%b(nstj)%nam
                  !nj = state_j%na(j)
                  !Cj = get_CI(state_j,j)
									nj = TargetStates%b(nstj)%na(j)
									Cj = TargetStates%b(nstj)%CI(j)

									mi = m_list_nr(ni)
									mj = m_list_nr(nj) 
									sturmi = sturm_ind_list_nr(ni)
									sturmj = sturm_ind_list_nr(nj)
 
                  !One-electron Hamiltonian for each Laguerre function 
                  lag_ham1el_m(ni,nj,ma) = H1el_Lag(sturmi,sturmj,mi,mj)
                  lag_ham1el_m(nj,ni,ma) = lag_ham1el_m(ni,nj,ma)
                  Ham = Ham + lag_ham1el_m(ni,nj,ma)  * Ci * Cj 
               end do
            end do
  
            e1me(nsti,nstj) = Ham
            e1me(nstj,nsti) = e1me(nsti,nstj)
  
         end do ! nstj
      end do ! nsti    
      !$OMP END PARALLEL DO


  
  end subroutine Hybrid_H1el_st_group





  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Subroutine: countConfigs
	!Purpose: counts the number of one-electon configurations for a given
	!         sturmian basis resolved in (kl) only. Fills arrays l_list,
	!         m_list, etc, with the (klm) indices of the fully (klm) resolved 
	!         basis. Uses data_in_sturm values to determine number of functions per l
	!Note: only called in non-linear molecule mode.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine countConfigs(basis, data_in_sturm, k_list, l_list, m_list, sturm_ind_list, num_func, rad_func)
		  use input_data  
			use sturmian_class
		  implicit none
			type(input):: data_in_sturm
			type(basis_sturmian_nr):: basis
			integer, dimension(:), allocatable:: k_list, l_list, m_list
			integer, dimension(:), allocatable:: sturm_ind_list
			integer, intent(inout):: num_func, rad_func
      integer:: ii, jj, kk    
			integer:: k, l, m, n    !Sturmian indices
			logical:: found
 
      !Triatomic molecules do not have good quantum numbers such as angular momentum and parity (like H2+).
      !Instead, the clamped nuclei electronic V-matrix is diagonal in different irreducible representations of a 
      !given point group. Since we are not currently using a symmetry adapted basis, we can only loop over lm
      num_func=0
      do l = data_in_sturm%labot, data_in_sturm%latop
         !Now need to loop over l and m, m is no longer conserved.
         do m = -l, l
            num_func = num_func + data_in_sturm%nps(l) !Array nps indexed from zero
         end do
      end do
      if (num_func .eq. 0) then
         print*, "ERROR: num_func=0: no basis functions. Stopping"
         error stop
      end if

      !Define arrays to keep track of k,l,m values for each index
      !Indexing scheme: specify (l,m), then k goes from 1 to N_l for each such pair 
      !i.e l -> m -> k_lm
      allocate(k_list(num_func), l_list(num_func), m_list(num_func), sturm_ind_list(num_func))
      ii = 0
      rad_func = 0
      kk = 0
      do l = data_in_sturm%labot, data_in_sturm%latop
         do k = 1, data_in_sturm%nps(l)
            do m = -l, l
               !jj  = kk 

               ii = ii + 1
               k_list(ii) = k
               l_list(ii) = l
  	           m_list(ii) = m
  	           !print*, "L,M,k,ii: ", l, m, k, ii
  
  	           !Indexing changed for non-diatomics, indices for basis now include l, k AND m. For a given l, the same radial basis
  	           !functions phi_{kl} are used for each value of m from -l to l. Introduce an array to track which radial basis 
  	           !function belongs to which index: sturm_ind_list
  	           !jj = jj + 1
  
  	           found = .false.
  	           n = 1
  	           do while (.not. found)
                  !print*, basis%b(n)%l, basis%b(n)%m, basis%b(n)%k, "L,K: ", l, k
  	              if ((basis%b(n)%l .eq. l) .and. (basis%b(n)%k .eq. k)) then
                     found = .true.
  	              else
                     n = n + 1
  	              end if
  	           end do
  
  	           sturm_ind_list(ii) = n
  	           !print*, jj
            end do
         end do
         !The radial part of the basis is identical for any given (l,k) pair, regardless of m
         !do k = 1, data_in_sturm%nps(l)
         !   rad_func = rad_func + 1
         !end do
				 rad_func = rad_func + data_in_sturm%nps(l)
  
         kk = kk + data_in_sturm%nps(l)
      end do
	
  end subroutine countConfigs


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!Subroutine: countConfigs
	!Purpose: counts the number of one-electon configurations for a given
	!         sturmian basis resolved in (kl) only. Fills arrays l_list,
	!         m_list, etc, with the (klm) indices of the fully (klm) resolved 
	!         basis. Uses data_in_sturm values to determine number of functions per l
	!Note: only called in non-linear molecule mode.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine countConfigs_old(basis, data_in_sturm, k_list, l_list, m_list, sturm_ind_list, num_func, rad_func)
		  use input_data  
			use sturmian_class
		  implicit none
			type(input):: data_in_sturm
			type(basis_sturmian_nr):: basis
			integer, dimension(:), allocatable:: k_list, l_list, m_list
			integer, dimension(:), allocatable:: sturm_ind_list
			integer, intent(inout):: num_func, rad_func
      integer:: ii, jj, kk    
			integer:: k, l, m, n    !Sturmian indices
			logical:: found
 
      !Triatomic molecules do not have good quantum numbers such as angular momentum and parity (like H2+).
      !Instead, the clamped nuclei electronic V-matrix is diagonal in different irreducible representations of a 
      !given point group. Since we are not currently using a symmetry adapted basis, we can only loop over lm
      num_func=0
      do l = data_in_sturm%labot, data_in_sturm%latop
         !Now need to loop over l and m, m is no longer conserved.
         do m = -l, l
            num_func = num_func + data_in_sturm%nps(l) !Array nps indexed from zero
         end do
      end do
      if (num_func .eq. 0) then
         print*, "ERROR: num_func=0: no basis functions. Stopping"
         error stop
      end if

      !Define arrays to keep track of k,l,m values for each index
      !Indexing scheme: specify (l,m), then k goes from 1 to N_l for each such pair 
      !i.e l -> m -> k_lm
      allocate(k_list(num_func), l_list(num_func), m_list(num_func), sturm_ind_list(num_func))
      ii = 0
      rad_func = 0
      kk = 0
      do l = data_in_sturm%labot, data_in_sturm%latop
         do m = -l, l
            jj  = kk 
            do k = 1, data_in_sturm%nps(l)
               ii = ii + 1
               k_list(ii) = k
               l_list(ii) = l
  	           m_list(ii) = m
  	           !print*, "L,M,k,ii: ", l, m, k, ii
  
  	           !Indexing changed for non-diatomics, indices for basis now include l, k AND m. For a given l, the same radial basis
  	           !functions phi_{kl} are used for each value of m from -l to l. Introduce an array to track which radial basis 
  	           !function belongs to which index: sturm_ind_list
  	           jj = jj + 1
  
  	           found = .false.
  	           n = 1
  	           do while (.not. found)
                  !print*, basis%b(n)%l, basis%b(n)%m, basis%b(n)%k, "L,K: ", l, k
  	              if ((basis%b(n)%l .eq. l) .and. (basis%b(n)%k .eq. k)) then
                     found = .true.
  	              else
                     n = n + 1
  	              end if
  	           end do
  
  	           sturm_ind_list(ii) = n
  	           !print*, jj
            end do
         end do
         !The radial part of the basis is identical for any given (l,k) pair, regardless of m
         do k = 1, data_in_sturm%nps(l)
            rad_func = rad_func + 1
         end do
  
         kk = kk + data_in_sturm%nps(l)
      end do
	
  end subroutine countConfigs_old


end module one_electron_func_group



