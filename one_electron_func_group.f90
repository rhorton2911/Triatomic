!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine: construct_1el_basis_nr_group
!Purpose: constructs basis for two electron structure as a mixture of
!         one-electron laguerre functions and molecular orbitals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine construct_1el_basis_nr(Number_one_electron_func)
    use input_data
    use grid_radial
    use sturmian_class
    !use vnc_module
    use one_electron_func
    !use target_states
    !use ovlpste1me
    !use natural_orbitals_module
    !use MPI_module

    implicit none

    integer, intent(in):: Number_one_electron_func
!
    integer:: ii, jj




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

    !Create space for basis of Mo's
    call new_basis_st(TargetStates,Number_one_electron_func,hlike,basis_type)

    if (myid == 0) write(*,'("allocated space for one-electron target states Nmax=",I5)') Number_one_electron_func
		!Store energy of ionised system, in this case, ionised H3++ (just three protons)
    if ( sum(indata%R(:)) .gt. 0.0_dpf ) then
       TargetStates%en_ion = Z1 * Z2 / Rd
       !Loop over nuclei
       do ii=1, 3
          do jj = ii+1, 3
             !Use law of cosines to compute distance between nuclei
	           cosij = cos(indata%theta(ii))*cos(indata%theta(jj)) + &
	                   sin(indata%theta(ii))*sin(indata%theta(jj))*cos(indata%phi(ii)-indata%phi(jj))
             Rij = sqrt(indata%R(ii)**2 + indata%R(jj)**2 - &
	           2*indata%R(ii)*indata%R(jj)*cosij)
	           !Account for degenerate case where nuclei coincide
             if (Rij .gt. 0.0_dpf) then
                TargetStates%en_ion = TargetStates%en_ion + dble(indata%charge(ii)*indata%charge(jj))/Rij
             end if
          end do
       end do
			 TargetStates%en_ion = 
    else ! Atomic Case
       TargetStates%en_ion = 0d0 
    end if      
    N_oneel = 0
          

    !Store one-electron molecular orbitals in TargetStates





    call sort_by_energy_basis_st(TargetStates)

    !call calc_spectro_factors(TargetStates, bst_nr, data_in%labot, data_in%latop)
    if (myid == 0) then
       call print_energy_basis_st(TargetStates)
    end if

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

    !Construct second laguerre basis, called bst_MSC
    call construct(bst_MSC, dataMSC, dataMSC%alpha_nl)   
         
    !Number of Sturmian functions in Laguerre basis and second laguerre basis
    nBasis_MSC =  basis_size(bst_MSC)
    nBasis_data = basis_size(bst_data)
    
    ! Combine both laguerre bases, store in bst_nr
    ! Calculates <Lfp|Lf> for e-H2 exchange matrix elements.
    ! rearrange.f90 will carry through rearrangement
    call combine_basis_nr(bst_nr,bst_MSC,bst_data,basis_type)
    max_latop = get_max_L(bst_nr) ! Largets l of laguerre basis

    if (allocated(lag_ham1el_m)) deallocate(lag_ham1el_m)
    allocate(lag_ham1el_m(basis_size(bst_nr),basis_size(bst_nr),-max_latop:max_latop))
    lag_ham1el_m(:,:,:) = 0d0
    
    if(allocated(e1me)) then
       deallocate(e1me) !
       deallocate(ovlpst) ! Basis functions overlaps
    end if
    
    !Count maximum number of Molecular states possible from second laguerre basis (2nd diagonalisation)
    itmp = 0
    do n = 1, nBasis_MSC
       lorb = get_ang_mom(bst_MSC%b(n))
       itmp = itmp + (2*lorb + 1)
    end do
    nst_Basis_MSC = itmp   
    allocate(e1me(nst_Basis_MSC, nst_Basis_MSC))   !Stores 1e matrix elements
    allocate(ovlpst(nst_Basis_MSC, nst_Basis_MSC)) !Overlap one- or two-electron functions/states 
    e1me(:,:) = 0d0
    ovlpst(:,:) = 0d0 
    
		!Overwrite molecular orbitals with second laguerre basis and populate e1me and ovlpst
    !Hybrid calculates <Lfp|H|Lf> for e-H2 exchange matrix elements, which will carry through rearrangement      
    call Hybrid_MSCbasis(bst_nr,TargetStates,nst_Basis_MSC,nBasis_MSC,lag_ham1el_m,basis_size(bst_nr),max_latop)
    
    !Liam added: save the original description of the one-electron states before rearrange is called - needed for natural orbitals
    call copy(TargetStates_unrearranged,TargetStates)

    !!$ Destruct temporary bases    
    call destruct(bst_MSC)
    call destruct(bst_data)

	  !Diagonalise 1e hamiltonian a second time in new sturmian basis, as CI coefficients are needed later
    !Molecular states from second diagonalisation using second sturmian basis stored in TargetStates1el
    call new_basis_st(TargetStates1el,Number_one_electron_func,hlike,basis_type)


	  !Save new states to TargetStates1el


	  !Call rearrange to express the basis in terms of functions with will defined (lm) 














!!$------------------------------   

!!$ Use the same basis for H2+ 1s configuration and other configurations
    if (data_in%inc == 1 .AND. dataMSC%MSC_nconfig == 0) then
       if (myid == 0) then
          print*,'Calling rearrange without molecular orbital basis.'
          print*,"inc, MSC_nconfig", data_in%inc, dataMSC%MSC_nconfig
       end if

!!$ Calculating <Lfp|H|Lf> for e-H2 exchange matrix elements.
!!$ will carry through rearrangement  
       if (allocated(lag_ham1el_m))  deallocate(lag_ham1el_m)
       allocate(lag_ham1el_m(basis_size(bst_nr),basis_size(bst_nr),-data_in%latop:data_in%latop))
       lag_ham1el_m(:,:,:) = 0d0
       call Hybrid_H1el_st(Nmax,lag_ham1el_m,basis_size(bst_nr),get_max_L(bst_nr))
    
       !Liam added: save the original description of the one-electron states before rearrange is called - needed for natural orbitals
       call copy(TargetStates_unrearranged,TargetStates)
       call rearrange(bst_nr,get_max_L(bst_nr),TargetStates,.TRUE.,lag_ham1el_m)
            
!!$ Use the different bases for H2+ 1sSg and 2pSu states and other configurations (2s, 3s, 3p, ..)       
    else if ( ABS(dataMSC%MSC_nconfig) /= 0  ) then  
       
       call destruct(bst_nr)
       
       call construct(bst_data, data_in)  ! Construct temporary basis from data.in inputs called bst_data
!       call construct(bst_MSC, dataMSC)   ! Construct temporary basis from dataMSC inputs called bst_MSC


       if (data_in%inc == 1) then !!$ inc must = 1 for two-electron targets. Checked in H12
          call rearrange(bst_nr,get_max_L(bst_nr),TargetStates,.TRUE.,lag_ham1el_m)
       end if
    
!       print*
!       print*,'Write state 1 in a file'
!       call write_1elState(bst_nr, TargetStates%b(7))
!       print*
    

!!$ Diagonalise H2 with 1sSg MO and AO 
       if ( dataMSC%MSC_nconfig >= 1 .AND. (.NOT. data_in%hlike) ) return

!!$ Above used for one and two-electron targets. Above can be used in H12.f90 now.      

          ! Molecular states from second diagonalisation using MSC basis stored in TargetStates1el
          call new_basis_st(TargetStates1el,Number_one_electron_func,hlike,basis_type)
          
!!$ Start second diagonalisation using MSC basis
          if(allocated(no) .OR. allocated (mo) ) then
             deallocate(no,mo)
          end if
          allocate(no(nBasis_MSC),mo(nBasis_MSC)) 
          no (:) = 0
          mo (:) = 0
          
          N_oneel = 0 
!!$ Partial Waves serperated by Angular Projection M and Parity 
          do ma = -data_in%Mt_max, data_in%Mt_max

             do ipar = 1, -1, -1

                if(data_in%good_parity) then
                  if(ipar == 0) cycle
                elseif(.not.data_in%good_parity) then
                  if(ipar /= 0) cycle
                endif
                  
                nst = data_in%nst(ABS(ma),ipar)

                if(nst == 0) cycle

                if(nst < 0) nst = -1

                ! Determine number of configurations nd per partial wave
                i = 0
                do n = 1, nst_Basis_MSC
                   m_n = get_ang_mom_proj( TargetStates%b(n) )
                   ipar_n = get_par( TargetStates%b(n) )
                   if ( m_n == ma .and. (ipar == ipar_n .or. .not.data_in%good_parity)) then
                      i = i + 1
                      no(i) = n ! Index to configuration for particular partial wave
                   end if
                end do
                nd = i 

                if(nst == -1) nst = nd
                
                if(nd == 0) then
                   if (myid == 0) print*,'Error: Molecular State Configuration. nd=0'
                   if (myid == 0) write(*, '("Symmetry and size: ma, ipar, nd = ",2I3,I5)') ma, ipar, nd 
                   stop
                endif
                
                mo(1:nd) = ma
                
                ! Temporary arrays
                allocate(H(nd,nd))
                allocate(b(nd,nd))
                allocate(CI(nd,nd))
                allocate(w(nd))
                
                H(:,:) = 0d0 
                b(:,:) = 0d0
!!$ Populate H = <f|H_T|i> and b = <f|i> using e1me and ovlpst.
                do ni = 1, nd
                   i = no(ni)
                   do nj = 1, ni
                      j = no(nj)
                      
                      H(ni,nj) = e1me(i,j)
                      b(ni,nj) = ovlpst(i,j)
                      
                      if ( Rd /= 0d0) then
                         H(ni,nj) = H(ni,nj) + Z1 * Z2 / Rd * b(ni,nj)
                      end if
                      
                      H(nj,ni) = H(ni,nj)
                      b(nj,ni) = b(ni,nj)
                      
                   end do ! nj
                end do ! ni
          
                if(data_in%N_core_el > 0) H = H + core_energy*b
                print*, '>> add core_en', sum(H), sum(b)
                
                if (myid == 0) print*
                if (myid == 0) write(*, '("Symmetry and size: ma, ipar, nd = ",2I3,I5)') ma, ipar, nd   
                ! Diagonalise
                !matz=2
                !call  rsg(nd,nd,H,b,w,matz,CI,irerr)
                allocate(IFAIL(nd), IWORK(5*nd))
                allocate(WORK(1))
                LWORK = -1
                call dsygvx( 1, 'V', 'I', 'U', nd, H, nd, b, nd, 0.0d0,0.0d0, 1,nd, 2*DLAMCH('S'), &
                  &NFOUND, w, CI, nd, WORK, LWORK, IWORK, IFAIL, irerr)
                
                LWORK = WORK(1)
                deallocate(WORK)
                allocate(WORK(LWORK))
                call dsygvx( 1, 'V', 'I', 'U', nd, H, nd, b, nd, 0.0d0,0.0d0, 1,nd, 2*DLAMCH('S'), &
                  &NFOUND, w, CI, nd, WORK, LWORK, IWORK, IFAIL, irerr)
                deallocate(IFAIL, IWORK, WORK)
                
                   
                if (myid == 0) then
                   write(*,'("irerr =",I3)') irerr
                   if(data_in%N_core_el > 0) then
                     print*, ' Energies in a.u. (including core energy)'
                   else
                     print*, ' Energies in a.u.'
                   endif
                   write(*,'(5F15.5)') (real(w(i)), i=1,nd)
                   print*
                   if(data_in%N_core_el > 0) then
                     print*, ' Energies in eV (including core energy)'
                   else
                     print*, ' Energies in eV'
                   endif
                   write(*,'(5F15.5)') (real(data_in%eV*w(i)), i=1,nd)
                end if
                
!!$ Construct states from MSC diagonalisation and store them in TargetStates1el        
                jstart = 1
                if(data_in%N_core_el > 0) jstart = 2
                jstop = min(nd,jstart+nst-1)
                do j = jstart, jstop
                   energy = w(j)
                   N_oneel = N_oneel + 1
                   tmp_sign_CI = SUM(CI(1:nd,j))
                   if (  tmp_sign_CI < 0d0 ) then
                      CI(1:nd,j) = -CI(1:nd,j)
                   end if
                   TargetStates1el%Nstates = N_oneel
                   call  construct_st(TargetStates1el%b(N_oneel),hlike,dble(ma),ipar,0.5d0,energy,j,nd,CI(1:nd,j),no(1:nd),mo(1:nd))
                end do
                
                
                deallocate(H)
                deallocate(b)
                deallocate(CI)
                deallocate(w)
                
             end do ! parity
          end do ! ma

          
!!$ First diagonalisation stored in TargetStates used Laguerre Basis function configuratiosn from data_in
!!$ Second diagonalisation stored in TargetStates1el. Used Target States of First diagonalisation
!!$ Save CI coefficients from Second diagonalisation in CIigiven. 
          
          Nmax1el = TargetStates1el%Nstates ! Second Diagonalisation number of states
          Nmax1el = min(Nmax1el,Number_one_electron_func) ! Number of states requested
          TargetStates1el%Nmax = Nmax1el
          
!!$ nBasis_MSC = Second Diagonalisation Laguerre basis size. Each target state will be constructed from MO/target states for that symmetry.
!!$ CIigiven = CI coefficients to configurations per state
!!$ igiven_max = Number of one-electron configurations per state
!!$ igiven =  Index to one-electron target state configurations
!!$ manst1el = Angular projection of configuration
          allocate( CIigiven(Nmax1el,nBasis_MSC), igiven(Nmax1el,nBasis_MSC), manst1el(Nmax1el))  
!!$ igiven_max(Nmax1el)    igiven_max(:) = 0
          igiven(:,:) = 0
          CIigiven(:,:) = 0d0         
          manst1el(:) = 0
          
          do nst = 1, Nmax1el
             manst1el(nst) = get_ang_mom_proj(TargetStates1el%b(nst)) ! Assign angular projection
             nam = get_nam(TargetStates1el%b(nst)) ! Number of one-electron configurations per state
             !                igiven_max(nst) = nam
             do j = 1, nam
                igiven(nst,j) = get_na(TargetStates1el%b(nst),j)
                CIigiven(nst,j) = get_CI(TargetStates1el%b(nst),j)
             end do
          end do
          
!!$ nsbst = size of combined basis (bst_nr). nobst pointer to all basis functions which make state
!!$ nobst = index to combined Laguerre Basis or rearranged Basis functions
          
          nspbst = basis_size(bst_nr) 
          allocate(nobst(nspbst),mobst(nspbst),CInobst(nspbst))
          nobst(:) = 0 
          mobst(:) = 0 
          CInobst(:) = 0d0
          
!!$ Multiply CI coefficients from first and second diagonalisation. This allows us to use the rearrange routine.
!!$ After combining CI coefficients construct new states.
          call new_basis_st(tmp1elst,Nmax1el,hlike,basis_type)            
          do nst = 1, Nmax1el
             nd = get_nam(TargetStates1el%b(nst)) ! Size of second diagonalisation
             
             call make_nomo(Targetstates,bst_nr,nd,CIigiven(nst,1:nd),igiven(nst,1:nd),manst1el(nst),nspbst,nobst,mobst,CInobst,itmpmax)               

             ipar = get_par(TargetStates1el%b(nst))
             rma = manst1el(nst)
             energy = get_energy_st(TargetStates1el%b(nst))
             nd = itmpmax
             tmp1elst%Nstates = nst

             call  construct_st(tmp1elst%b(nst),.true.,rma,ipar,0.5d0,energy,get_inum_st(TargetStates1el%b(nst)),nd,CInobst(1:nd),nobst(1:nd),mobst(1:nd))
          end do


!!$ Replace molecular orbitals with short-ranged functions
!!$ e1me and ovlpst are calculated here for H2+ and H2
!             call replace_molecular_orbitals(tmp1elst)
          
          if (data_in%inc == 1) then
             call rearrange(bst_nr,get_max_L(bst_nr),tmp1elst,.FALSE.)
          end if

          call destruct_basis_st(TargetStates)
          call new_basis_st(TargetStates,tmp1elst%Nstates,.true.,0)
          TargetStates%Nmax = tmp1elst%Nmax
          TargetStates%Nstates = tmp1elst%Nstates
          if ( Rd /= 0d0 ) then
             TargetStates%en_ion = Z1 * Z2 / Rd
          else ! Atomic Case
             TargetStates%en_ion = 0d0
          end if

          if(data_in%N_core_el > 0) then
            TargetStates%en_ion = TargetStates%en_ion + core_energy !core_energy in vnc_module - does not contain nuclear repulsion
          endif

          do nst = 1, tmp1elst%Nstates
             call copy_st(TargetStates%b(nst), tmp1elst%b(nst))
          end do
          call destruct_basis_st(tmp1elst)

          call sort_by_energy_basis_st(TargetStates)
!!$ calc_spectro_factors sets target states: l_major, n - principal quantum
!number and correspond Laguerre basis k stored in inum
          call calc_spectro_factors(TargetStates, bst_nr, min(dataMSC%labot,data_in%labot), max(dataMSC%latop,data_in%latop))
          if (myid == 0) then
             call print_energy_basis_st(TargetStates)
          end if
       
    else
    
       if (myid==0) print*,"rearrange and Molecular State Configuration basis were not called"

    end if ! inc rearrange and MSC



    
end subroutine construct_1el_basis_nr


