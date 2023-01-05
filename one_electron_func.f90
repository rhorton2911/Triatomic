subroutine construct_1el_basis_nr(Number_one_electron_func)
    use input_data
    use grid_radial
    use sturmian_class
    use vnc_module
    use one_electron_func
    use target_states
    use ovlpste1me
    use natural_orbitals_module
    use MPI_module

    implicit none

    integer, intent(in):: Number_one_electron_func
!
    real*8, dimension(:,:), allocatable:: H, b, b2, CI
    real*8, dimension(:), allocatable:: w
    integer:: i, j, n, nd, N_oneel, jstart,jstop, ipar, ipar_n, l_n, li, lj, Nmax, lorb, itmp
    real*8:: al, Rd, Z1, Z2
!
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

    !FOR THE DSYGVX SUBROUTINE
    real*8, dimension(:), allocatable :: WORK
    integer, dimension(:), allocatable :: IWORK, IFAIL
    integer :: LWORK, NFOUND
    real*8, external :: DLAMCH

    !For OMP around the H matrix
    integer :: numpairs, pair
    integer, dimension(:,:), allocatable :: pairlist
    integer, external :: OMP_GET_THREAD_NUM
  
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



    basis_type = data_in%calculation_type
    hlike = .true.
    Rd = data_in%Rd
    Z1 = data_in%Z1
    Z2 = data_in%Z2

! Create space for basis of  1-electron pseudostates
    call new_basis_st(TargetStates,Number_one_electron_func,hlike,basis_type)
    if (myid == 0) write(*,'("allocated space for one-electron target states Nmax=",I5)') Number_one_electron_func
!!$ Ionisation of H2+, we remove 1 electron. The potentials that are left are the repulsion between the two protons. 
!!$ Hence I(H2+)=E(H2+)-E(H2++)=(K_1+V_1+1/R)-(1/R) In fixed nuclei approximation
    if ( Rd /= 0d0 ) then
       TargetStates%en_ion = Z1 * Z2 / Rd
    else ! Atomic Case
       TargetStates%en_ion = 0d0 
    end if      
    if(data_in%N_core_el > 0) then
      TargetStates%en_ion = TargetStates%en_ion + core_energy !core_energy in vnc_module - does not contain nuclear repulsion
    endif
    N_oneel = 0
                
! For all given (la,nd,al) construct Sturmian basis: is kept in module one_electron_func
    if (myid == 0) print*, 'Start making nonrrel  Sturmian basis'
    call construct(bst_nr,data_in)
    
   
    !if(data_in%N_core_el > 0) call remove_core_from_bst(bst_nr)

    nspm =  basis_size(bst_nr) ! number of sturmian finctions
    if (myid == 0) print*, 'Size of Sturmian basis: nspm=', nspm
    allocate(no(nspm),mo(nspm)) 
    no (:) = 0
    mo (:) = 0
    
    do ma = data_in%Mt_min, data_in%Mt_max
       do ipar = 1,-1,-1

         if(data_in%good_parity) then 
           if(ipar == 0) cycle
           nst = data_in%nst(ma,ipar)
         elseif(.not. data_in%good_parity) then
           if(ipar /= 0) cycle
           nst = data_in%nst(ma,0)
         endif
         
         if(nst == 0) cycle

         if(nst < 0) nst = -1
         
          
!!$   form set of (one-electron) configurations to be used to digonalise the molecular ion (H2+, HeH++, etc)
          i = 0
          do n=1,nspm
             l_n = get_ang_mom(bst_nr%b(n))
             ipar_n = (-1)**(l_n)
             if( l_n .ge. ma .and. (ipar .eq. ipar_n .or. .not.data_in%good_parity)) then
                i = i + 1
                no(i) =n
             endif
          end do
          
          nd = i    ! number of sturmian orbitals to be used in diagonlization
          if(nd .eq. 0) then
             if (myid == 0) print*,'***ERROR: nd = 0'
             if (myid == 0) print*, 'm, par:', ma, ipar
             stop
          endif
          mo(1:nd) = ma

          !Liam added: let the user request ALL states
          if(nst == -1) nst = nd
          
          if (myid == 0) write(*, '(/"Symmetry and size: ma, ipar, nd = ",2I3,I5)') ma, ipar, nd

          ! Temporary arrays
          allocate(H(nd,nd))
          allocate(b(nd,nd))
          allocate(b2(nd,nd))
          allocate(CI(nd,nd))
          allocate(w(nd))
          
          H(:,:) = 0d0
          b(:,:) = 0d0
         
!!$ Calculate H matrix
!!$ Here we rely on the special form of Laguerre functions of order (l+1)
!!$ in order to  calculate overlaps and 1/r integrals analytically
          b(1:nd,1:nd) = bst_nr%ortint(no(1:nd),no(1:nd))
        
          numpairs = (nd*(nd+1))/2
          allocate(pairlist(numpairs,2))
          pair=0
          do ni=1, nd
            do nj=1, ni
              pair = pair + 1
              pairlist(pair,1) = ni
              pairlist(pair,2) = nj
            enddo
          enddo

!          !$OMP PARALLEL DO &
!          !$OMP DEFAULT(SHARED) PRIVATE(ni, nj, i, pi, j, pj, tmp) 
          do pair=1, numpairs
            ni = pairlist(pair,1)
            nj = pairlist(pair,2)
            
              i = no(ni)
              pi => bst_nr%b(i)
              j = no(nj)
              pj => bst_nr%b(j)

!!$ Nuclei repulsion term for diatomic molecules:
!!$ Note: H_2 Hamiltonian H = H2+ + H2+ + V12.
!!$ Have implemented a -1/R factor in H2 (H12.f90) structure code to compensate for the 2/R.
              if ( Rd /= 0d0) then
                 tmp = Z1*Z2/Rd * b(ni,nj)
                 H(ni,nj) = H(ni,nj) + tmp 
              end if

!!$ One-electron configurations matrix elements
!!$ H_T= K + V + Z^2/R 
              call Hlagorb(bst_nr,i,j,ma,tmp)


              H(ni,nj) = H(ni,nj) + tmp 
              H(nj,ni) = H(ni,nj)
             
!!$           if (myid == 0) print*, i, j, ni, nj, H(ni,nj), tmp
                
          end do
!          !$OMP END PARALLEL DO

          deallocate(pairlist)


          if(any(H /= H)) error stop '*** ERROR: H /= H'


          b2 = b !save overlap matrix for renormalising CI matrix after diagonalisation
          if(data_in%N_core_el > 0) H = H + core_energy*b
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

          !Liam added to remove noise in CI coefficients
          !Set any CI coefficients < CI_min to zero
          where(abs(CI) < data_in%CI_min) CI = 0.0d0
          
          !Renormalise CI
          do n=1, nd
            norm=0.0d0
            do i=1,nd
              do j=1, nd
                norm = norm + CI(i,n)*CI(j,n)*b2(i,j)
              enddo!j
            enddo !i
            CI(:,n) = CI(:,n) / sqrt(norm)
          enddo !n

             
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

!!$ Create basis of 1-electron pseudostates from Laguere basis and CI coef.
!!$          if(la .le. data_in%la_core) then
!!$             jstart = data_in%npr_core(la) - la + 1
!!$          else 
!!$             jstart = 1
!!$          endif

          jstart = 1
          if(data_in%N_core_el > 0 .and. .not. data_in%pseudo_pot) jstart = 2
          jstop = min(nd,jstart+nst-1)
          do j=jstart,jstop
             energy = w(j)

             N_oneel = N_oneel + 1
             tmp_sign_CI = SUM(CI(1:nd,j))
             if(  tmp_sign_CI .lt. 0d0 ) then
                CI(1:nd,j) = -CI(1:nd,j)
             endif
             TargetStates%Nstates = N_oneel
             call  construct_st(TargetStates%b(N_oneel),hlike,dble(ma),ipar,0.5d0,energy,j,nd,CI(1:nd,j),no(1:nd),mo(1:nd))
!!$ NOTE here we might need to acount for degeneracy of the molecule energy levels for ma .ne. 0
!!$ this could lead to introdicing addtonal state with -ma 
             if( ma .ne. 0 ) then
                N_oneel = N_oneel + 1
                TargetStates%Nstates = N_oneel
                call  construct_st(TargetStates%b(N_oneel),hlike,-dble(ma),ipar,0.5d0,energy,j,nd,CI(1:nd,j),no(1:nd),-1*mo(1:nd))
             endif
          end do
             
          deallocate(H)
          deallocate(b)
          deallocate(b2)
          deallocate(CI)
          deallocate(w)
       end do ! parity
    end do ! ma
    
    if( TargetStates%Nstates .ne. TargetStates%Nmax ) then
       if (myid == 0) print*, 'one_electron_func.f90: TargetStates%Nstates .ne. TargetStates%Nmax:', TargetStates%Nstates, TargetStates%Nmax
       stop
    endif
          
    if (myid == 0) print*
    call sort_by_energy_basis_st(TargetStates)

    call calc_spectro_factors(TargetStates, bst_nr, data_in%labot, data_in%latop)
    if (myid == 0) then
       call print_energy_basis_st(TargetStates)
    end if

    if(data_in%N_core_el > 0) return !TODO: REMOVE


!!$  populate arary e1me to be used later in H12 diagonalization
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
!!$ This constructs a Laguerre basis with different values of alpha
       call construct(bst_MSC, dataMSC, dataMSC%alpha_nl )   ! Construct temporary basis from dataMSC inputs called bst_MSC
            
       ! Number of Sturmian functions in Laguerre basis and MSC basis (2nd diagonalisation) 
       nBasis_MSC =  basis_size(bst_MSC)
       nBasis_data = basis_size(bst_data)
       
!!$ Combine Bases
!!$ Calculates <Lfp|Lf> for e-H2 exchange matrix elements.
!!$ rearrange.f90 will carry through rearrangement
       call combine_basis_nr(bst_nr,bst_MSC,bst_data,basis_type)
       max_latop = get_max_L(bst_nr) ! Largets l of laguerre basis

       if (allocated(lag_ham1el_m)) deallocate(lag_ham1el_m)
       allocate(lag_ham1el_m(basis_size(bst_nr),basis_size(bst_nr),-max_latop:max_latop))
       lag_ham1el_m(:,:,:) = 0d0
       
       if(allocated(e1me)) then
          deallocate(e1me) !
          deallocate(ovlpst) ! Basis functions overlaps
       end if
       
       ! Maximum number of Molecular states possible from MSC basis (2nd diagonalisation)
       itmp = 0
       do n = 1, nBasis_MSC
          lorb = get_ang_mom(bst_MSC%b(n))
          itmp = itmp + (2*lorb + 1)
       end do
       nst_Basis_MSC = itmp   ! Maximum number of Molecular states from MSC basis 
       allocate(e1me(nst_Basis_MSC, nst_Basis_MSC))   ! H2+ Ham matrix element = sum C <f(1)|K_1+V_1|i(1)>
       allocate(ovlpst(nst_Basis_MSC, nst_Basis_MSC)) ! Overlap one- or two-electron functions/states 
       e1me(:,:) = 0d0
       ovlpst(:,:) = 0d0 
       
!!$ Copy H2+ ground state to 1s configuration and 2pPu state to 2p m=0 configuration.
!!$ All other configurations built from bst_MSC Laguerre basis 
!!$ functions. Configurations are stored in TargetStates. Populate e1me and ovlpst
!!$ Hybrid calculates <Lfp|H|Lf> for e-H2 exchange matrix elements.
!!$ will carry through rearrangement      
       call Hybrid_MSCbasis(bst_nr,TargetStates,nst_Basis_MSC,nBasis_MSC,lag_ham1el_m,basis_size(bst_nr),max_latop)
       
       !Liam added: save the original description of the one-electron states before rearrange is called - needed for natural orbitals
       call copy(TargetStates_unrearranged,TargetStates)

       if (data_in%inc == 1) then !!$ inc must = 1 for two-electron targets. Checked in H12
          call rearrange(bst_nr,get_max_L(bst_nr),TargetStates,.TRUE.,lag_ham1el_m)
       end if
    
!       print*
!       print*,'Write state 1 in a file'
!       call write_1elState(bst_nr, TargetStates%b(7))
!       print*
    
!!$ Destruct temporary bases       

       call destruct(bst_MSC)
       call destruct(bst_data)

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


 subroutine Hybrid_MSCbasis(bst,TargetStates,Nmax,NBasis_F5,lag_ham1el_m,Lag_nmax,maxL)
!
! Mark: This routine obtains H2+ ground state from Laguerre functions specified in data.in.
! and specifies the excited A.O. Laguerre functions from F5.    
! 
    use input_data
    use sturmian_class
    use state_class
    use ovlpste1me
    use MPI_module
    use vnc_module
    
    implicit none

    type(basis_sturmian_nr), intent(in) :: bst   ! this is Sturmian basis 
    type(basis_state), intent(inout) :: TargetStates
    integer, intent(inout):: Nmax, NBasis_F5
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

    if ( dataMSC%MSC_nconfig == 0 ) stop ! Should not be here!!!
    ! Number of Sturmians from F5 file
    if (myid == 0) print*, 'Size of MSC bst basis: ',NBasis_F5

    if(dataMSC%MSC_nconfig > TargetStates%Nmax) then
      if(myid==0) write(*,*) '*** ERROR: requested more 1-el states in the second basis than kept from the&
        & first diagonalisation'
      error stop
    endif


!!$ Record ground state and vary index to data.in Laguerre functions
    basis_type = data_in%calculation_type
    hlike = .true.
    call new_basis_st(IonStates,ABS(dataMSC%MSC_nconfig),hlike,basis_type)
    do n = 1, ABS(dataMSC%MSC_nconfig)
       call copy_H2plus_st(IonStates%b(n),TargetStates%b(n),NBasis_F5)
    end do
    call destruct_basis_st(TargetStates)

!!$ create target state basis of size Nmax
    call new_basis_st(TargetStates,Nmax,hlike,basis_type)
    if ( data_in%Rd /= 0d0 ) then
       TargetStates%en_ion = data_in%Z1 * data_in%Z2 / data_in%Rd
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
    do n=1, NBasis_F5
       k = get_k(bst%b(n))
       lst = get_ang_mom(bst%b(n))
       ipar = (-1)**lst

       do ma = -lst,lst
          ist = ist + 1
          TargetStates%Nstates = ist
          mo(1) = ma
          no(1) = n
    
          replaced_with_MO = .false.
          !Liam: added this to allow including an arbitrary number of MOs in the hybrid basis
          do j=1, dataMSC%MSC_nconfig
            MolOrb => IonStates%b(j)
            lMO = get_l_majconf(MolOrb)
            kMO = get_n_majconf(MolOrb) - lMO 
            mMO = get_ang_mom_proj(MolOrb)
            labelMO = trim(adjustl(get_label(MolOrb)))

            if(kMO == k .and. lMO == lst .and. mMO == ma) then
              call copy_st(TargetStates%b(ist),MolOrb)
              replaced_with_MO = .true.
              if(myid == 0) write(*,'(" Replaced Sturmian k=",I0,", l=",I0,", m=",I0," with one-electron state ",I0," ( ",A," )")') &
                & k, lst, ma, j, labelMO
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
    call destruct_basis_st(IonStates)

    if ( ist /= Nmax ) then
       if (myid == 0) print*, 'one_electron_func.f90: Hybrid_MSCbasis(): ist .ne. Nmax :', ist, Nmax
       stop
    end if

!!$ Calculates e1me and lag_ham1el_m
    call Hybrid_H1el_st(Nmax,lag_ham1el_m, basis_size_nr(bst), get_max_L(bst))

    do nst1 = 1, Nmax
       state1 => TargetStates%b(nst1)
       m1 = get_ang_mom_proj(state1)

       do nst2 = 1, nst1
          state2 => TargetStates%b(nst2)
          m2 = get_ang_mom_proj(state2)

          if ( m1 /= m2 .OR. (data_in%good_parity .and. get_par(state1) /= get_par(state2)) ) cycle

          ovlpst(nst1,nst2) = ovlp_st(state1,state2,bst)
          ovlpst(nst2,nst1) =  ovlpst(nst1,nst2)

!!$         if (myid == 0) print*, nst1, nst2,  e1me(nst1,nst2), ovlpst(nst1,nst2)

       end do
    end do

!!$    a =  H1el_st(TargetStates%b(2),TargetStates%b(3))    
!!$    call Hlagorb(bst,2,3,0,b)
!!$    if (myid == 0) print*, '>>> a,b:',a,b
!!$
!!$    a = bst%ortint(2,3)
!!$    b =  ovlp_st(TargetStates%b(2),TargetStates%b(3))
!!$    if (myid == 0) print*, '!!! a,b:',a,b


  end subroutine Hybrid_MSCbasis

subroutine Hybrid_H1el_st(Nmax,lag_ham1el_m,nmax_lag, max_latop)
!!$ Calculates the one-electron molecular Hamiltonian matrix elements 
!!$ lag_ham1el_m = <Lp|H|L> for Laguerre basis functions for particular m.
!!$ e1me = <np|H|n> for target states with underlying laguerre basis functions 

  use ovlpste1me
  use target_states
  use MPI_module
  use state_class
  use input_data

  implicit none
  integer, intent(in):: Nmax ! Number of states      
  integer, intent(in):: nmax_lag, max_latop ! Number of Lag function and max l
  real*8, dimension(nmax_lag,nmax_lag,-max_latop:max_latop), intent(inout):: lag_ham1el_m

  type(state), pointer:: state_j, state_i
  real*8:: Ci, Cj, Ham
  real*8:: H1el_Lag
  integer:: nsti, nstj, i, j, n, ni, nj             ! Indexes
  integer:: ma


  do nsti = 1, Nmax
   state_i => TargetStates%b(nsti)
   ma = state_i%m

   do nstj = 1, nsti
   state_j => TargetStates%b(nstj)

    Ham = 0d0

     if (state_i%m /= state_j%m) cycle
     if (data_in%good_parity .and. (state_i%parity /= state_j%parity)) cycle
     do i = 1, state_i%nam
        ni = state_i%na(i)
        Ci = get_CI(state_i,i)

           do j = 1, state_j%nam
              nj = state_j%na(j)
              Cj = get_CI(state_j,j)

!!$ One-electron Hamiltonian for each Laguerre function 
              lag_ham1el_m(ni,nj,ma) = H1el_Lag(ni,nj,ma)
              lag_ham1el_m(nj,ni,ma) = lag_ham1el_m(ni,nj,ma)
              Ham = Ham + lag_ham1el_m(ni,nj,ma)  * Ci * Cj

           end do
        end do

        e1me(nsti,nstj) = Ham
        e1me(nstj,nsti) = e1me(nsti,nstj)

     end do ! nstj
  end do ! nsti    

end subroutine Hybrid_H1el_st


subroutine replace_molecular_orbitals(tmp1elst)

!!$ This subroutine  reverts one-electron molecular orbitals stored in  S' back to short-ranged atomic orbitals stored in S.
!!$ This new basis is stored back in S' and one-electron matrix elements are calculated for e1me and ovlpst, which allows 
!!$ calculations of H2 structure.
!!$ Both the data structures S' = tmp1elst and S = TargetStates point to the same one-electron functions stored
!!$ in bst_nr. At present THIS SUBROUTINE NEEDS TO BE BEFORE THE FINAL REARRANGE OF S', otherwise bst_nr will change for S.
! S = one-electron molecular orbital basis resulting from second and/or first diagonalisation = TargetStates 
! S' = one-electron molecular orbital basis resulting from second diagonalisation = tmp1elst


  use input_data
  use ovlpste1me
  use target_states
  use MPI_module
  use state_class
  use sturmian_class
  use one_electron_func

  implicit none
  type(basis_state), intent(inout):: tmp1elst
  integer, dimension(TargetStates%Nmax):: S_SRF_index
  integer, dimension(tmp1elst%Nmax):: Sp_SRF_index
  integer:: Nmax, k, l_majconf, S_SRF_Nmax, Sp_SRF_Nmax
  integer:: j, i, n, nst, jst, ind_j, ind_n, mnst
  real*8:: CIj, CIn
  real*8, dimension(tmp1elst%Nmax,tmp1elst%Nmax):: ovlpst_temp, e1me_temp
  integer:: z
  real*8:: temp_ovlp, ovlp2

! In case there is not short-ranged functions and initialising 
Sp_SRF_index(:) = -1

  Nmax = basis_size_st(tmp1elst)
if ( data_in%use_sub > 0 .AND. (.NOT. data_in%hlike) )  then

!!$ Loop through S basis = TargetStates
!!$ and find short ranged functions
   j = 0
do i = 1, basis_size_st(TargetStates)
   l_majconf = get_l_majconf(TargetStates%b(i))
   if ( l_majconf < dataMSC%labot_diff .OR. l_majconf > dataMSC%latop_diff ) cycle
   if ( get_inum_st(TargetStates%b(i))  > dataMSC%nps_diff(l_majconf)  ) cycle
   if ( dataMSC%nps_diff(l_majconf) <= 0 ) cycle
   j = j + 1
   S_SRF_index(j) = i  ! index to short ranged functions in S
!   print*,j,i,get_inum_st(TargetStates%b(i)) 
end do
   S_SRF_Nmax = j ! Number of orbitals/target state in S built from short ranged functions
!!$ Loop through S' basis = tmp1elst
!!$ and find index to MO that should be replaced with short ranged functions



   i = 0
do nst = 1, Nmax
   k = get_inum_st(tmp1elst%b(nst))
   l_majconf = get_l_majconf(tmp1elst%b(nst))
   mnst = get_ang_mom_proj(tmp1elst%b(nst))
   do j = 1, S_SRF_Nmax
! Match short-ranged function with function from S' basis 
      if ( k == get_inum_st(TargetStates%b(S_SRF_index(j))) .AND. l_majconf == get_l_majconf(TargetStates%b(S_SRF_index(j))) .AND. mnst == get_ang_mom_proj(TargetStates%b(S_SRF_index(j))) ) then
         i = i + 1
         Sp_SRF_index(i) = nst  ! Index to orbitals/target states in S' that need to be replaced with a SRF
!         print*,"k,  l, m",k,l_majconf,mnst
!         print*,"i, S', S",i,nst, S_SRF_index(j)

!!$ TESTING
!if ( k == 2 .AND. l_majconf == 0 .AND. mnst == 0 ) then
!print*,""
!        print*,"2s function expansion with MO basis"
!        if( 1 < get_nam(TargetStates%b(S_SRF_index(j)))  ) print*,"issues"
!
!!!$ May need to replace the 2 index to the 2s function        ind_n = get_na(TargetStates%b(S_SRF_index(2)),1)
!         ovlp2 = 0d0
!        do jst = 1, Nmax
!          temp_ovlp = 0d0
!        do n = 1, get_nam(TargetStates%b(S_SRF_index(j)))
!         ind_n = get_na(TargetStates%b(S_SRF_index(j)),n)
!         CIn = get_CI(TargetStates%b(S_SRF_index(j)),n)
!         do z = 1, get_nam(tmp1elst%b(jst))
!            ind_j = get_na(tmp1elst%b(jst),z) ! Point to functions in bst_nr, needs to be done before final rearrange
!            CIj = get_CI(tmp1elst%b(jst),z)
!               temp_ovlp = temp_ovlp + bst_nr%ortint(ind_j,ind_n) * CIj * CIn
! 	     	end do ! z 
!            end do ! n 
!           ovlp2 = ovlp2 + temp_ovlp * temp_ovlp
!          print*,"nst ovlp",jst, temp_ovlp
!         end do ! jst
!	print*,"<2s|2s>",ovlp2
!print*,""
!end if
!!$*********************



!!$ Copying orbital/state from S to S' basis
         call copy_st(tmp1elst%b(nst),TargetStates%b(S_SRF_index(j)))
         exit
      end if ! SRF
   end do ! j
end do ! nst
Sp_SRF_Nmax = i ! Number of orbitals/target state in S' that have been replaced by SRF 

  ovlpst_temp(:,:) = 0d0
  e1me_temp(:,:) = 0d0

!!$ Loop through replaced states/orbitals of S' == tmp1elst
!!$ Calculating these states/orbitals of S' ovlp and e1me 
!!$ that have atleast one-short ranged function
do i = 1, Sp_SRF_Nmax
   nst = Sp_SRF_index(i)

   do jst = 1, Nmax 
         do j = 1, get_nam(tmp1elst%b(jst))
            ind_j = get_na(tmp1elst%b(jst),j) ! Point to functions in bst_nr, needs to be done before final rearrange
            CIj = get_CI(tmp1elst%b(jst),j)
            do n = 1, get_nam(tmp1elst%b(nst))
               ind_n = get_na(tmp1elst%b(nst),n) ! Point to functions in bst_nr, needs to be done before final rearrange
               CIn = get_CI(tmp1elst%b(nst),n)
               ovlpst_temp(nst,jst) = ovlpst_temp(nst,jst) + bst_nr%ortint(ind_n,ind_j) * CIj * CIn
               e1me_temp(nst,jst) = e1me_temp(nst,jst) + bst_nr%ham1el(ind_n,ind_j) * CIj * CIn
            end do
         end do
      if ( ANY(Sp_SRF_index == jst) ) then ! < phi | phi >
! continue
      else ! < phi | phi' >  
        ovlpst_temp(jst,nst) = ovlpst_temp(nst,jst)
        e1me_temp(jst,nst) = e1me_temp(nst,jst)
      end if ! S or S' function overlap
   end do ! jst loop over MOorbitals/states of S' molecular orbital 
end do ! i loop / nst loop over  MOorbitals/states of S' molecular orbital that were replaced

!!$ Calculate the rest of the S' e1me and ovlpst matrix
!!$ for states/orbitals that aren't SRF i.e. can use properties of target states

end if ! hlike and SRF
          if(allocated(e1me)) deallocate(e1me)
          if(allocated(ovlpst)) deallocate(ovlpst)
          allocate(e1me(Nmax, Nmax))
          allocate(ovlpst(Nmax, Nmax))
          e1me(:,:) = 0d0
          ovlpst(:,:) = 0d0
          do jst = 1, Nmax
             if ( ANY(Sp_SRF_index == jst)   ) then
             e1me(jst,:) = e1me_temp(jst,:)
             e1me(:,jst) = e1me_temp(:,jst)    
             ovlpst(jst,:) = ovlpst_temp(jst,:)
             ovlpst(:,jst) = ovlpst_temp(:,jst)

             else ! not being replaced with a short-ranged function (SRF)
                ovlpst(jst,jst) = 1d0
                if ( data_in%Rd /= 0d0 ) then
                   e1me(jst,jst) = get_energy_st(tmp1elst%b(jst)) - data_in%Z1 * data_in%Z2 / data_in%Rd
                else ! Atomic Case
                   e1me(jst,jst) = get_energy_st(tmp1elst%b(jst))
                end if ! atomic
             end if ! SRF
          end do ! j 

!      do nst = 1, Nmax
! do jst = 1, Nmax
!if ( get_ang_mom_proj(tmp1elst%b(nst)) == get_ang_mom_proj(tmp1elst%b(jst)) ) then
!        write(*,'("nst ",I2," l ",I2," m ",I2," jst ",I2," l ",I2," m ",I2)') nst, get_l_majconf(tmp1elst%b(nst)), get_ang_mom_proj(tmp1elst%b(nst)), jst, get_l_majconf(tmp1elst%b(jst)), get_ang_mom_proj(tmp1elst%b(jst))
!
!        print*,"e1me,ovlpst",e1me(nst,jst),ovlpst(nst,jst)
!else if ( get_ang_mom_proj(tmp1elst%b(nst)) == get_ang_mom_proj(tmp1elst%b(jst)) .AND. e1me(nst,jst) /= 0d0 ) then
!print*,""
!print*,"**********************************"
!
!        write(*,'("nst ",I2," l ",I2," m ",I2," jst ",I2," l ",I2," m ",I2)') nst, get_l_majconf(tmp1elst%b(nst)), get_ang_mom_proj(tmp1elst%b(nst)), jst, get_l_majconf(tmp1elst%b(jst)), get_ang_mom_proj(tmp1elst%b(jst))
!        print*,"e1me,ovlpst",e1me(nst,jst),ovlpst(nst,jst)
!print*,"**********************************"
!print*,""
!end if
! end do
!      end do




end subroutine replace_molecular_orbitals

subroutine remove_core_from_bst(bst)
  use sturmian_class
  use vnc_module
  use target_states
  implicit none

  type(basis_sturmian_nr), intent(inout) :: bst
  integer :: i, n, l, j, jj
  type(state), pointer :: orb
  type(sturmian_nr), pointer :: sturm

  do i=1, CoreOrbitals%Nmax
    orb => CoreOrbitals%b(i)
    n = orb%n
    l = orb%l
    do j=1, bst%n
      sturm => bst%b(j)
      if(sturm%k == n-l .and. sturm%l == l) exit
    enddo !j
    if(j <= bst%n) then !found a sturmian to remove
      do jj=j+1, bst%n
        call copy(bst%b(jj-1),bst%b(jj))
      enddo !jj
      call destruct(bst%b(bst%n))
      bst%n = bst%n-1
    endif

  enddo !i
end subroutine remove_core_from_bst
