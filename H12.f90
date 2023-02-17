subroutine structure12

  use input_data
  use sturmian_class
  use state_class
  use target_states
  use ovlpste1me
  use one_electron_func
  use MPI_module
  use natural_orbitals_module
!!$ Mark added
  use vmat_exch_module

  implicit none  

  integer:: lorbmax
  parameter( lorbmax = 20 )
  
  character(LEN=20):: target
  integer:: Mmax
  integer, dimension(:,:,:), allocatable:: nstate  

  integer:: ma, ip, is, itmp, inc, labot, latop, nsp, ncm, nc, ncp, l,m, Nmax, Nmax1el, nst, nstp
  real*8:: rma, ris,  en_ion,   tmp_sign_CI
  integer, dimension(:), allocatable:: phase
  logical:: hlike
  real*8, allocatable,dimension(:):: alpha
  integer, allocatable,dimension(:):: nps 
  type(input):: dataF5  
  integer, allocatable,dimension(:):: no1,no2, mo1, mo2
  real*8:: tmp, resultb, resultH, Etot

  integer:: matz, ioerr, i, j
  real*8, dimension(:,:), allocatable:: H, b, CI
  real*8, dimension(:), allocatable:: w
  integer:: nicm, ni, icheck
  integer, allocatable,dimension(:):: ncore, mncore
  integer:: nsp1,nsp2, nsp1p, nsp2p, mnsp1, mnsp2, mnsp1p, mnsp2p, Nmax_bound, n
  integer:: l12max, ltmp, M12max
  integer, dimension(0:lorbmax):: nk1, nk2, nk1_saved
  integer:: l_ion_core
  integer, dimension(0:lorbmax):: n_ion_core
! Mark
  integer:: Ntmp, lorb, NBasis_F5, NBasis_data, N1el
  type(basis_sturmian_nr):: bst_F5,bst_data
!!$ Non-uniqueness
  integer:: k_core, l_core, l_majcore
  integer::  nspm_temp, nst_core, nam, n_core
!
  common /icheck_potl_Nmax_open/ icheck_Nmax_open 
  integer:: icheck_Nmax_open

  integer :: mode !Liam added for automatically including all states
  integer :: nstates

  !Liam added for parallising H matrix calculation
  integer, dimension(:), allocatable :: nclist, ncplist
  integer :: jj, num_loop

  !FOR THE DSYGVX SUBROUTINE
  real*8, dimension(:), allocatable :: WORK
  integer, dimension(:), allocatable :: IWORK, IFAIL
  integer :: LWORK, NFOUND
  real*8, external :: DLAMCH

  logical :: do_natorb, done_natorb, read_CI
  integer :: natorb_GS_M, natorb_GS_par, natorb_GS_spin, ncm_in
  real*8, dimension(:), allocatable :: natorb_CI

  character*4 :: labelA, labelB
  type(state), pointer :: orbA, orbB
  integer :: mAA, pA, mB, pB, nj
  character(len=6) :: configLabel
  character*17, parameter :: lchars='spdfghiklmnoqrtuv'
  character*20, parameter :: symchars = 'u gSPDFGHIKLMNOQRTUV'
  type(state), pointer :: stateL,stateR
  integer :: n1, n2, l1, l2, nst1, nst2
  character(len=:), allocatable :: stlabel
  integer :: domConfig
  logical :: F5_ex

  type(state), pointer :: tempstate

  hlike = .false.

  done_natorb = .false.
  do_natorb = (data_in%num_nat_orb > 0 .or. data_in%num_nat_orb == -1)
  if(do_natorb) then
    natorb_GS_M = data_in%natorb_GS_M
    natorb_GS_par = data_in%natorb_GS_par
    natorb_GS_spin = data_in%natorb_GS_spin
  endif
  1000 continue

  nk1(:) = 0
  nk2(:) = 0
  n_ion_core(:) = 0

  inquire(file='F5',exist=F5_ex)
  if(data_in%new_input_format .and. F5_ex) then
    if(myid==0) write(*,*) '*** ERROR: F5 file is not used with new data.in format - please delete it'
    error stop
  elseif(.not.data_in%new_input_format .and. .not.F5_ex) then
    if(myid==0) write(*,*) '*** ERROR: F5 file missing'
    error stop
  endif


  if (myid==0) then
     print*
     print*,'Start CI routines for quasi two-electron di-atomic molecule'
     print*
  end if

  if(.not.data_in%new_input_format) then
   
    open(3,file='F5')
    read(3,*) target
    if (myid==0) write(*,*) 'target'
    if (myid==0) write(*,*) target
  
    read(3,*) Mmax

    allocate(nstate(0:Mmax,-1:1,0:1))
    nstate(:,:,:) = 0

    read(3,*) (nstate(l,1,0), l=0,Mmax) ! positive parity, singlet
    read(3,*) (nstate(l,-1,0), l=0,Mmax) ! negative parity, singlet
    read(3,*) (nstate(l,1,1), l=0,Mmax) ! positive parity, triplet
    read(3,*) (nstate(l,-1,1), l=0,Mmax) ! negative parity, triplet
    if(.not.data_in%good_parity) then
      nState(:,0,:) = nState(:,1,:) + nState(:,-1,:)
      where(nState < 0) nState = -1 !To keep functionality that "-1" = "all states"
    endif

  else !new data.in format - no F5 file

    Mmax = data_in%Mmax2el
    allocate(nstate(0:Mmax,-1:1,0:1))
    nstate(:,:,:) = data_in%nst_2e

  endif

  if (myid==0) write(*,*) 'Mmax abs value'
  if (myid==0) write(*,*) Mmax
  
  if (data_in%good_parity) then
    if(myid==0) then
      write(*,*) 'nState(M,+,0) -- singlet, positive parity :'
      write(*,*) '     ', ( nState(l,1,0), l=0,MMax )
      write(*,*) 'nState(M,-,0) -- singlet, negative parity :'
      write(*,*) '     ', ( nState(l,-1,0), l=0,MMax )
      write(*,*) 'nState(M,+,1) -- triplet, positive parity :'
      write(*,*) '     ', ( nState(l,1,1), l=0,MMax )
      write(*,*) 'nState(M,-,1) -- triplet, negative parity :'
      write(*,*) '     ', ( nState(l,-1,1), l=0,MMax )
    endif
  else
    if(myid==0) then
      write(*,*) 'nState(M,0,0) -- singlet, no parity:'
      write(*,*) '     ', ( nState(l,0,0), l=0,MMax )
      write(*,*) 'nState(M,0,1) -- triplet, no parity:'
      write(*,*) '     ', ( nState(l,0,1), l=0,MMax )
    endif
  end if
    

  basis_type = data_in%calculation_type
  hlike = .false.

!!$  Target states were created for one-electron target: one-electon ion, get the ground state energy of the ion and then destroy the old target states
  en_ion = get_energy_st(TargetStates%b(1))  
  if (myid==0) print*, 'en_ion=', en_ion
!  call new_basis_st(TargetStates2el,Nmax,hlike,basis_type)
  TargetStates2el%en_ion = en_ion
  
!$$   make s.p. Laguerre basis
  dataF5%calculation_type = basis_type

!!$ Use Laguerre basis or MSC basis from data.in file
!!$ Laguerre basis used with MSC=1 and having first diagonalisation with just the 1s function
  if ( ABS(dataMSC%MSC_nconfig)  <= 2 .AND. dataMSC%MSC_nconfig /= 0 ) then  ! Molecular State Configuration Basis
     if (myid==0) print*,"Using Molecular orbital/atomic orbital Basis", dataMSC%MSC_nconfig
     dataF5%latop = dataMSC%latop
     dataF5%labot = dataMSC%labot
     allocate(dataF5%alpha(0:dataMSC%latop))
     allocate(dataF5%nps(0:dataMSC%latop))
     dataF5%nps(:) = dataMSC%nps(:)
     dataF5%alpha(:) = dataMSC%alpha(:)
  else ! Full Molecular State Basis
     if (myid==0) print*,"Using H2+ molecular orbitals from first basis in data.in", dataMSC%MSC_nconfig
     dataF5%latop = data_in%latop
     dataF5%labot = data_in%labot
     allocate(dataF5%alpha(0:data_in%latop))
     allocate(dataF5%nps(0:data_in%latop))
     dataF5%nps(:) = data_in%nps(:)
     dataF5%alpha(:) = data_in%alpha(:)
  end if

  if (data_in%inc /= 1 ) then
     if (myid==0) print*,"Error: MSC methods requires the use of rearrange for two-electron targets. Set inc = 1 in data.in"
     stop
  end if


  latop = dataF5%latop 
  labot = dataF5%labot 
  if (myid==0) write(*,'("labot,latop: ",2I5)') labot, latop
  if ( allocated(nps) ) then
     deallocate(alpha,nps)
  end if
  allocate(alpha(0:latop))
  allocate(nps(0:latop))
  nps(:) = dataF5%nps(:) 
  alpha(:) = dataF5%alpha(:) 
  if (myid==0) write(*,'("nps(l),alpha(l): ",20(I5,F10.4))') (nps(l),alpha(l), l=labot,latop)

  if (myid==0) print*
  
  if(.not.data_in%new_input_format) then
    read(3,*) inc   ! inc=1 for call rearange
  else
    inc = 1
  endif
  if (myid==0) write(*,*) inc, '            */ inc'

  if ( ABS(dataMSC%MSC_nconfig) >= 3 .OR. dataMSC%MSC_nconfig == 0 ) then 
     if (myid==0) print*
     if (myid==0) print*, 'Use all target states from one-electron diagonalisation to diagonalise two-electron H2'
     if (myid==0) print*, 'MSC basis in data.in is ignored'

     N1el = basis_size_st(TargetStates)
     if( allocated(e1me) ) deallocate(e1me,ovlpst)
     allocate( e1me(N1el,N1el), ovlpst(N1el,N1el) )
     e1me(:,:) = 0d0
     ovlpst(:,:) = 0d0     
!!$ One-electron target states have + Z^2/R. Need to remove from one-electron Ham.
!!$ Re-added in two-electron matrix elements
     do n=1,N1el
        if ( data_in%Rd /= 0d0 ) then
           e1me(n,n) = get_energy_st(TargetStates%b(n)) - data_in%Z1 * data_in%Z2 / data_in%Rd
        else ! Atomic Case
           e1me(n,n) = get_energy_st(TargetStates%b(n))
        end if
        ovlpst(n,n) = 1d0       
     end do
     !     end if
     if (myid==0) print*
  else if ( dataMSC%MSC_nconfig > 0 ) then 
     if (myid==0) print*
     if (myid==0) print*," MO Hybrid bases for H2+ and H2. "
     if (myid==0) print*," basis = some MOs with the rest as Laguerre atomic orbitals"
     if (myid==0) print*
  else if ( dataMSC%MSC_nconfig < 0 .AND. dataMSC%MSC_nconfig >= -2 ) then
     if (myid==0) print*
     if (myid==0) print*," MO Hybrid bases for H2+ and H2. "
     if (myid==0) print*,"basis = accurate 1sSg and maybe 2pSu MO with the rest as H2+ MO/target states"
     if (myid==0) print*
  else
     if (myid==0) print*,"MSC_nconfig defined incorrectly",  dataMSC%MSC_nconfig 
     stop   
  end if
!!$----------------------
  ! Maximum number of one-electron orbitals for core
  itmp = 0 
  do n = 0, latop
     itmp = itmp + nps(n) * (2*n + 1) 
  end do
  Ntmp = itmp
  ncm = Ntmp   
  allocate(ncore(ncm),mncore(ncm))   

  if(.not.data_in%new_input_format) then
    read(3,*) l_ion_core,(n_ion_core(l), l=0,l_ion_core)
  else
    l_ion_core = data_in%l_ion_core
    n_ion_core(0:l_ion_core) = data_in%n_ion_core
  endif

  if (myid==0) print*,  l_ion_core, (n_ion_core(l), l=0,l_ion_core)
  if(l_ion_core .gt. 20) then
     if (myid==0) print*, 'increase size of array n_ion_core(20) to at least l_ion_core=', l_ion_core
     stop
  endif

  if (myid==0) print*

!!$ read CI model only once, might want to change later to read per symmetry as in atomic code.
!!$ nk2 refers to the maximum allowed number of basis functions per l of outer electron.
  if(.not.data_in%new_input_format) then
    read(3,*) l12max, M12max   ! 
    l12max = min(l12max, latop)
    read(3,*) (nk2(i), i=0,l12max)   ! 
    read(3,*) (nk1(i), i=0,l12max)   ! 
    close(3)
  else
    l12max = data_in%l12max
    M12max = data_in%M12max
    nk2(0:l12max) = data_in%nkout(0:l12max)
    nk1(0:l12max) = data_in%nkin(0:l12max)
  endif

  if (myid==0) write(*,*) l12max, M12max, '           */ loutmax, M12max'
  if (myid==0) write(*,*) (nk2(i), i=0,l12max), '           */ nkout'
  if(do_natorb .and. .not. done_natorb) then
    nk1 = nk2
    nstate = 0
    nstate(natorb_GS_M,natorb_GS_par,natorb_GS_spin) = 1
  endif

  if (myid==0) write(*,*) (nk1(i), i=0,l12max), '           */ nkin'
!!$
  if(l_ion_core .gt. l12max) then
     if (myid==0) print*,'l_ion_core .gt. l12max', l_ion_core, l12max
     l_ion_core = l12max
     if (myid==0) print*, 'redefine l_ion_core to l12max value'
  endif
  nk1_saved(:) = nk1(:)


  nicm = 0
  ncm = -1

!!$ MARK: Make orthogonal one-electron basis for solving non-uniqueness
!!$ for electron scattering. Need to modify e1me and ovlpst.
!!$ Using Gram-Schmidt Orhtogonalisation because molecular orbital
!!$ and Laguerre basis functions with different alpha will not be
!!$ orthogonal even if we us an orthongal Laguerre basis of (2l+2).
 if(myid == 0) print*, '---->>>>  data_in%non_uniq =',  data_in%non_uniq
 if ( data_in%non_uniq) then
     call Orthonormalise_one_electron_basis(TargetStates,bst_nr)
  end if

  Nmax = 0
  
  if(myid==0 .and. data_in%print_2el_config) open(unit=777,file='two-el-config.out',action='write',status='replace')

  do mode=0, 1
  
    if(mode==1) call new_basis_st(TargetStates2el,Nmax,hlike,basis_type)

    do ma = -Mmax, Mmax
  
  !!$ For Positron Scattering
  !!$ Use FC for states with ma larger than M12max. Greatly reduces Number of states
       if ( abs(ma) > M12max ) then
          nk1(:) = 0
          nk1(0) = 1
       else
          nk1(:) = nk1_saved(:)
       end if
       
       rma = ma
       do ip=1,-1,-1    ! parity
        if (data_in%good_parity) then
           if (ip == 0) cycle
        else
           if (ip /= 0) cycle
        end if

          do is = 0,1   ! spin
             ris = is

             nstates = nstate(abs(ma),ip,is)

             
             if(nstates <= 0 .and. nstates /= -1) cycle
             
             if (myid==0 .and. mode==1) print*
             if (myid==0 .and. mode==1) write(*,'("Symmetry (M,par,spin) : ",3I3)') ma,ip,is

             if(nstates > 0) then
               if(mode==0) then
                 Nmax = Nmax + nstates
                 cycle
               else
                 if(myid==0) write(*,*) 'Number of requested states = ', nstates
               endif
             elseif(nstates == -1) then
               if(myid==0 .and. mode==1) write(*,*) 'Number of requested states = ALL'
             endif

             
  !!$ Set up list of configurations, the CI model in F5 file is read in config12()
  !!$    first find number of configurations
             ncm = 0  ! number of configuration
             if(allocated(no1)) then
                deallocate(no1,mo1)
             endif
             allocate(no1(TargetStates%Nmax*TargetStates%Nmax),mo1(TargetStates%Nmax*TargetStates%Nmax))
             no1 (:) = 0
             mo1 (:) = 0
             
             call config12_tmp_st(TargetStates,ma,ip,is,l12max,nk1(0:l12max),nk2(0:l12max),l_ion_core,n_ion_core(0:l12max),ncm,no1,mo1)  ! call first time to find ncm
  !!$ Determening core-orbitals
  !!$ This has been moved to here because for example type 1 core orbitals for 2p are not all included  
  !!$ for example: m=0, pi=1 state: configs 2p(m=-1) 2p(m=1), 2p(m=1) 2p(m=-1) are kept here.
  !!$ in config12_st the anti symmetric confgis are excluded i.e. 2p(m=1) 2p(m=-1)
  !!$ and the 2p(m=1) orbital is not saved as a core orbital!
             do nc=1,ncm
  !              if (myid==0) print*, no1(nc), no2(nc), mo1(nc), mo2(nc)
                icheck = 0  ! checking if this orbital wa sa;ready included
                do ni=1,nicm
                   if(ncore(ni).eq.no1(nc) .and. mncore(ni).eq.mo1(nc)) then
                      icheck = 1
                      exit
                   end if
                end do
                if(icheck.eq.0) then
                   nicm = nicm + 1
                   ncore(nicm) = no1(nc)
                   mncore(nicm) = mo1(nc)
  !                 if (myid==0) print*, 'nicm, ncore(nicm)', nicm,ncore(nicm)
  !                 if (myid==0) print*, 'k, l, m', get_inum_st(TargetStates%b(ncore(nicm))), get_l_majconf(TargetStates%b(ncore(nicm))), mncore(nicm)
                end if
             enddo
             
             if(allocated(no1)) deallocate(no1,mo1)
             if(allocated(no2)) deallocate(no2,mo2) !coded like this for a reason
             allocate(no1(ncm),mo1(ncm),no2(ncm),mo2(ncm)) 
             no1 (:) = 0
             mo1 (:) = 0
             no2 (:) = 0
             mo2 (:) = 0
             call config12_st(TargetStates,ma,ip,is,l12max,nk1(0:l12max),nk2(0:l12max),l_ion_core,n_ion_core(0:l12max),ncm,no1,mo1,no2,mo2)!,mode) ! call second time to populate arrays           
             
             if (myid==0 .and. mode==1) write(*,*) 'Number of configurations = ', ncm

             if(nstates == -1) then
               nstates = ncm
               nstate(abs(ma),ip,is) = nstates
             endif
             if(mode==0) then
               Nmax = Nmax + nstates
               cycle
             endif
           
             if (ncm < nStates) stop 'ERROR : Not enough configurations to make up the requested number of states.'
             
             if(nstates == 0) then
               if(myid==0) print*, 'MPS ', ma, ip, is, 'has no configs'
               cycle
             endif


  !!$ END Set up list of configurations           
             
  
  !!$  ------------------------------------------------                  
             
  !!$   Temporary arrays
             if(allocated(H)) deallocate(H,b)
             if(allocated(w)) deallocate(w,CI)
             allocate(w(ncm),CI(ncm,ncm))
           
             if(do_natorb .and. .not.done_natorb) then
               inquire(file='natorb_CI',exist=read_CI)
               if(read_CI) then
                 open(unit=111,file='natorb_CI',action='read')
                 read(111,*) ncm_in, w(1)
                 if(ncm /= ncm) then
                   write(*,'("***ERROR reading the natorb_CI file: number of configurations (",I0,") in file /= ",I0)') ncm_in, ncm
                   error stop
                 endif
                 do jj=1, ncm
                   read(111,*) CI(jj,1)
                 enddo !jj
                 close(111)
                 print*, '  GROUND-STATE CI EXPANSION READ FROM FILE TO GO INTO NATURAL ORBITALS'
                 goto 2000
               endif !read_CI
             endif !do_natorb
  
  !!$ Form H matrix
             allocate(H(ncm,ncm),b(ncm,ncm))
             num_loop = ncm*(ncm+1)/2
             allocate(nclist(num_loop), ncplist(num_loop))
             num_loop = 0
             do nc = 1, ncm
               do ncp = nc, ncm
                 num_loop = num_loop + 1
                 nclist(num_loop) = nc
                 ncplist(num_loop) = ncp
               enddo !ncp
             enddo !nc
             if(num_loop /= ncm*(ncm+1)/2) error stop 'ERROR in spheroidal_12.f90: incorrect num_loop'
             
             Nmax1el = basis_size_st(TargetStates)
              !Liam added OMP here - speeds things up a bit
              !$OMP PARALLEL DO DEFAULT(SHARED) & 
              !$OMP PRIVATE(nc,ncp,nsp1,mnsp1,nsp2,mnsp2,nsp1p,mnsp1p,nsp2p,mnsp2p,resultH,resultb) &
              !$OMP SCHEDULE(DYNAMIC)
              do jj = 1, num_loop
                nc = nclist(jj)
                ncp = ncplist(jj)
  
                nsp1 = no1(nc) 
                mnsp1 = mo1(nc) 
                nsp2 = no2(nc) 
                mnsp2 = mo2(nc) 
  
                nsp1p = no1(ncp) 
                mnsp1p = mo1(ncp) 
                nsp2p = no2(ncp) 
                mnsp2p = mo2(ncp) 
                
                b(nc,ncp) = 0d0
                H(nc,ncp) = 0d0
                
                call H12me_st_notortog(data_in%Rd,is,TargetStates,Nmax1el,e1me,ovlpst,bst_nr,nsp1,nsp2,nsp1p,nsp2p,mnsp1,mnsp2,mnsp1p,mnsp2p,resultH,resultb)
  !              call H12me_st_notortog_test(data_in%Z,data_in%Rd,is,TargetStates,Nmax1el,e1me,ovlpst,bst_nr,nsp1,nsp2,nsp1p,nsp2p,mnsp1,mnsp2,mnsp1p,mnsp2p,resultH,resultb)
  !              print*,no2(nc),no2(ncp),resultH,resultb 
  
                
                b(nc,ncp) = resultb
                b(ncp,nc) = resultb
                H(nc,ncp) = resultH 
                H(ncp,nc) = resultH 
                  
             enddo !jj
             !$OMP END PARALLEL DO

             deallocate(nclist,ncplist)
  
             if(.false.) then
                print*, 'H matrix'
                do nc=1,ncm
                   write(*,'(100Es12.4)') (real(H(nc,i)), i=1,ncm)
                enddo
                print*, 'b matrix'
                do nc=1,ncm
                   write(*,'(100Es12.4)') (real(b(nc,i)), i=1,ncm)
                enddo
             endif
  
             matz=2
             !call rsg(ncm,ncm,H,b,w,matz,CI,ioerr)
             
             allocate(IFAIL(ncm), IWORK(5*ncm))
             allocate(WORK(1))
             LWORK = -1

!             print*, 'BST_NR:'
!             do jj=1, bst_nr%n
!               print*, 'jj, k, l, m, alpha:', jj, bst_nr%b(jj)%k, bst_nr%b(jj)%l, bst_nr%b(jj)%m, bst_nr%b(jj)%alpha
!             enddo
!             print*, ''
!
!
!
!              print*, 'ONE_ELECTRON_STATES:'
!              do jj=1, TargetStates%Nmax
!                print*, 'jj = ', jj
!                print*, 'na:', TargetStates%b(jj)%na
!                print*, 'CI:', TargetStates%b(jj)%CI
!                print*, ''
!              enddo
!             
!              print*, 'OVLPST:'
!             do jj=1, TargetStates%Nmax
!               print*, ovlpst(:,jj)
!             enddo
!             print*, ''
!              
!             print*, 'E1Me:'
!             do jj=1, TargetStates%Nmax
!               print*, e1me(:,jj)
!             enddo
!             print*, ''
                
             call dsygvx( 1, 'V', 'I', 'U', ncm, H, ncm, b, ncm, 0.0d0,0.0d0, 1,ncm, 2*DLAMCH('S'), &
               &NFOUND, w, CI, ncm, WORK, LWORK, IWORK, IFAIL, ioerr)
           
             LWORK = WORK(1)
             deallocate(WORK)
             allocate(WORK(LWORK))
             
             call dsygvx( 1, 'V', 'I', 'U', ncm, H, ncm, b, ncm, 0.0d0,0.0d0, 1,ncm, 2*DLAMCH('S'), &
               &NFOUND, w, CI, ncm, WORK, LWORK, IWORK, IFAIL, ioerr)
          
             !Liam added to remove noise in CI coefficients
              !Set any CI coefficients < CI_min to zero
            where(abs(CI) < data_in%CI_min) CI = 0.0d0
           
             deallocate(WORK, IFAIL, IWORK)

             if (myid==0) then
                write(*,'("ioerr =",I3)') ioerr
                if(ioerr /=0) error stop
                print*, ' Energies in a.u.'
                write(*,'(5F15.5)') (real(w(i)), i=1,ncm)
                print*, ' Energies in eV - en_ion'
                write(*,'(5F15.5)') (real(data_in%eV*(w(i)-en_ion)), i=1,ncm)
             end if
           
             2000 continue !continue here from goto statement before diagonalisation if ground-state CI already read in for natorbs
             
             if(nstate(ABS(ma),ip,is) .gt. ncm) then
                if (myid==0) print*,'H12.f: nstate(ma,ip,is) .gt. ncm :', nstate(ABS(ma),ip,is), ncm
                if (myid==0) print*, 'increase number of s.p. orbitals or decrease number of states in F5'
                stop
             endif
             
             if(allocated(phase)) deallocate(phase)
             allocate(phase(ncm))
             phase(:) = (-1)**(is)
             do nc=1,nstate(ABS(ma),ip,is)
                
                ! Mark: Fix sign of coefficients
                tmp_sign_CI = SUM(CI(1:ncm,nc))
  !!$              print*,"nc,SUM_CI",nc,tmp_sign_CI
                if(  tmp_sign_CI .lt. 0d0 ) then
                   CI(1:ncm,nc) = -CI(1:ncm,nc)
                endif
  
                TargetStates2el%Nstates = TargetStates2el%Nstates + 1
              
                if(do_natorb .and. .not.done_natorb .and. .not.read_CI) then !write ground-state CI to file for next time
                  if(TargetStates2el%Nstates > 1) error stop '*** ERROR writing natorb_CI: Nstates > 1 should not be possible'
                  open(unit=111,file='natorb_CI',action='write',status='new')
                  write(111,*) ncm, w(1)
                  do jj=1, ncm
                    write(111,*) CI(jj,1)
                  enddo
                  close(111)
                endif
                
                tempstate => TargetStates2el%b(TargetStates2el%Nstates)
                call construct_st(tempstate,hlike,rma,ip,ris,w(nc),nc,ncm,CI(1:ncm,nc),no1,mo1,no2,mo2,phase)

  !!$ Positive and negative m's are calculated independently.
  !!$              if(ma .ne. 0) then
  !!$                 TargetStates2el%Nstates = TargetStates2el%Nstates + 1
  !!$
  !!$!                call construct_st(TargetStates2el%b(TargetStates2el%Nstates),hlike,-rma,ip,ris,w(nc),nc,ncm,CI(1:ncm,nc),no1,mo1,no2,mo2,phase)
  !!$                 call construct_st(TargetStates2el%b(TargetStates2el%Nstates),hlike,-rma,ip,ris,w(nc),nc,ncm,CI(1:ncm,nc),no1,-mo1,no2,-mo2,phase)
  !!$              endif
  
                !if (1.3<=data_in%Rd .and. data_in%Rd<=1.5) call set_label(TargetStates2el%b(TargetStates2el%Nstates))
                !Liam: I commented out the above and replace with:
                call set_label(TargetStates2el%b(TargetStates2el%Nstates))
                ! - the state labels are assigned according to the state index within each symmetry, and as a consequence of
                !   the von Neumann-Wigner non-crossing rule (avoided crossings) the state order within each symmetry must
                !   be independent of R.
               
               stateL => TargetStates2el%b(TargetStates2el%Nstates) 
               stLabel = trim(adjustl(stateL%label))
               if(TargetStates2el%b(TargetStates2el%Nstates)%m < 0) stLabel = stLabel//'-'
               if(myid==0 .and. data_in%print_CI) then
                 open(unit=7777,file='CI-'//stLabel,action='write',status='replace')
                 domConfig = maxloc(abs(CI(:,nc)),1)
                 do j=1, ncm
                   nst1 = no1(j)
                   nst2 = no2(j)
                   n1 = get_n_majconf(TargetStates%b(nst1))
                   n2 = get_n_majconf(TargetStates%b(nst2))
                   l1 = get_l_majconf(TargetStates%b(nst1))
                   l2 = get_l_majconf(TargetStates%b(nst2))
  
                   write(configLabel,'(2(I0,A1))') n1, lchars(l1+1:l1+1), n2, lchars(l2+1:l2+1)
  
                   write(7777,'(A6,2X,ES12.5)',advance='no') configLabel, CI(j,nc)
                   if(j == domConfig) then
                     write(7777,*) '  *** DOMINANT CONFIGURATION'
                     stateL%domconfig = trim(adjustl(configLabel))
                   else
                     write(7777,*)
                   endif
  
                 enddo
                 close(7777)
               endif !print_CI
              

             enddo
                
          enddo ! end is loop
       enddo ! end ip loop
    enddo    ! end ma loop

  enddo !end mode loop
  if(myid==0 .and. data_in%print_2el_config) close(777)

  if (myid==0) print*
  if(TargetStates2el%Nstates .ne. Nmax) then
    print*, targetstates2el%nstates, nmax
     print*,'Number of calculated states is less than the number of ordered states',myid
     print*,'Check the CI configuration list, most likely ncm=0 for some symmetry',myid
     stop
  endif

  TargetStates2el%nicm = nicm   ! nicm was determined in above CI loops
  allocate(TargetStates2el%ncore(1:nicm))
  allocate(TargetStates2el%mncore(1:nicm))
  TargetStates2el%ncore(1:nicm) = ncore(1:nicm)
  TargetStates2el%mncore(1:nicm) = mncore(1:nicm)
  if (myid==0) print*,'nicm=',nicm
!  do l=1, TargetStates2el%nicm
!     if (myid==0) print*,l, TargetStates2el%ncore(l), TargetStates2el%mncore(l)
!  enddo

  if (myid==0) print*, ''
  if (myid==0) print*, 'Printing core-orbitals'
  do ni = 1, nicm
     if (myid==0) print*, 'ni, ncore(ni)', ni,ncore(ni)
     if (myid==0) print*, 'k, l, m', get_inum_st(TargetStates%b(ncore(ni))), get_l_majconf(TargetStates%b(ncore(ni))), mncore(ni)
  enddo
  if (myid==0) print*, ''
 
  
  if (myid==0) print*, " Sort states"
  call sort_by_energy_basis_st(TargetStates2el)
!!$ find Nmax_open to be compared with Nmax_open in check_potl_ordering()
  Etot = get_energy_st(TargetStates2el%b(1)) + data_in%energy
  do n=1,Nmax
     tmp = Etot - get_energy_st(TargetStates2el%b(n))
     if(tmp .gt. 0.0 ) then  ! open state
        TargetStates2el%Nmax_open = n
     else
        exit
     endif
  enddo

  icheck_Nmax_open = 0
  if(data_in%SF_tol == 1.0d0) call check_potl_ordering(TargetStates2el,icheck_Nmax_open)   !! in 
  !call write_basis(bst_nr, 'bst_2e_0')
  !call write_configs(bst_nr, TargetStates, TargetStates2el, .false.)

  if (myid==0) then
     print*, " Print energies"
     call print_energy_basis_st(TargetStates2el)
  end if
  Nmax_bound = 0

  do n=1,Nmax
     tmp = get_energy_st(TargetStates2el%b(n))
     if(tmp .lt. 0) then
        Nmax_bound = n
        else
           exit
     endif
  enddo
  TargetStates2el%Nmax_bound = Nmax_bound
  
  if(do_natorb.and..not.done_natorb) then
    call natorbitals(TargetStates2el%b(1))
    done_natorb = .true.
    deallocate(nstate,dataF5%nps,dataF5%alpha,alpha,nps,ncore,mncore)
    call destruct_basis_st(TargetStates2el)
    goto 1000
  endif
  

!!$  if (myid==0) print*, ' CI representation of two-electron target states'
!!$  do n=1,Nmax
!!$     do i=1,get_nam(TargetStates2el%b(n))
!!$        if (myid==0) write(*,'("n,i, CI(i):",2i5,2X,e15.5)') n,i, get_CI(TargetStates2el%b(n),i)       
!!$     enddo
!!$  enddo
 
!!$ Mark: Testing configurations
!  do nst = 1,TargetStates2el%Nmax 
!     do nstp = 1, TargetStates2el%b(nst)%nam   
!         write(*,'("nst,n1,n2,CI,"3I3,F10.5)'),nst,get_na(TargetStates2el%b(nst),nstp,1),get_na(TargetStates2el%b(nst),nstp,2),get_CI(TargetStates2el%b(nst),nstp)
!     end do   ! nstp
!  end do  ! nst

!  do nst = 1, TargetStates2el%Nmax
!     do nstp = nst, TargetStates2el%Nmax
!        call Test_Overlap_st(nstp,nst)
!     end do
!  end do
 
  if(inc == 1) then
     call rearrange12_st(TargetStates,TargetStates2el,TargetStates2el%nicm,TargetStates2el%Nmax,TargetStates%Nmax)
!     call rearrange12_new(TargetStates, TargetStates2el)
  end if
  !call write_basis(bst_nr,'bst_2e_rearr_0')
  !call write_configs(bst_nr, TargetStates, TargetStates2el, .true.)
!  call orthonormal12(TargetStates, TargetStates2el)

  ! Re-calculate the one-electron orbital states matrix elements.
  deallocate(ovlpst, e1me)
  Nmax1el = basis_size(TargetStates)
  allocate( ovlpst(Nmax1el,Nmax1el), e1me(Nmax1el,Nmax1el) )
  
  !Liam added OMP to speed this up
  !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(Nmax1el, TargetStates, bst_nr, ovlpst, e1me)
  do nstp = 1, Nmax1el
     do nst = nstp, Nmax1el
        
        resultb = 0d0; resultH = 0d0
        do ip = 1, get_nam(TargetStates%b(nstp))
           nsp1p = get_na(TargetStates%b(nstp),ip)
           do i = 1, get_nam(TargetStates%b(nst))
              nsp1 = get_na(TargetStates%b(nst),i)
              
              tmp = get_CI(TargetStates%b(nstp),ip)*get_CI(TargetStates%b(nst),i)
              resultb = resultb + tmp*bst_nr%ortint(nsp1p,nsp1)
              resultH = resultH + tmp*bst_nr%Ham1el(nsp1p,nsp1)
           end do
        end do

        ovlpst(nstp,nst) = resultb; e1me(nstp,nst) = resultH
        ovlpst(nst,nstp) = resultb; e1me(nst,nstp) = resultH
     end do
  end do
  !$OMP END PARALLEL DO

  deallocate(H,b,CI,w,phase)
 
!!$ Mark: for debugging.
!!$ Check that <np|H|n> and <np|n> work by recalculating
!!$ two-electron state energy with the arrays 
!!$ used in the exchange-matrix elements

! call H2st_test()

!!$ n_ion_core = all configurations for these orbitals i.e. 
!!$ usually config up to 1s,nl... n_ion_core,nl
!!$ n_ion_core orbitals have max l l_ion_core
!!$ Solving non-uniqueness
  if ( data_in%non_uniq ) then
     allocate(is_core_POMO(TargState_Orthonorm%Nmax),is_core_MO(TargetStates%Nmax))
!!$ is_core_MO = 0 1s1s,1s2s,1s2p,..,1snl <=> nl orbital type 0 
!!$ is_core_MO = 1 2s2s <=> 2s orbital type 1
!!$ is_core_MO = 2 1s1s,1s2s,1s2p,.., 1snl <=> 1s orbital type 2 
!!$ is_core_MO - to determine which orthonormal molecular orbitals
!!$ are core orbitals for the projection operator


     if ( data_in%non_uniq .AND. Mmax > M12max ) then
        print*,"2e-target states Mmax > M12max"
        print*,"Non-uniqueness won't be solved properly with this input"
        stop
     end if

     is_core_POMO(:) = 0
     is_core_MO(:) = 0

     do n_core = 1, TargetStates2el%nicm
!!$ If we used rearrange all core-orbital one-electron target states
!!$ are the first nicm states in the array

!!$ is_core_POMO is for the orthogonal one-electron basis before H2 diagonalisation 
        if (inc == 1) then
!           nst_core = n_core
!!$ nst_core = n_core will not work here because the Orthnormal basis 
!!$ is in a different order compared to rearrange TargetStates basis
           nst_core = TargetStates2el%ncore(n_core)
        else 
           nst_core = TargetStates2el%ncore(n_core)
        end if 

        k_core = get_inum_st(TargState_Orthonorm%b(nst_core))
        l_majcore = get_l_majconf(TargState_Orthonorm%b(nst_core))
        if (l_ion_core >= 0) then
           if ( k_core <= n_ion_core(l_majcore) ) then
              is_core_POMO(nst_core) = 2
!!$              print*,"PO: st_ind, k, l", nst_core, k_core, l_majcore,"type 2"
           else if ( k_core <= nk1(l_majcore) ) then
              is_core_POMO(nst_core) = 1
!!$              print*,"PO: st_ind, k, l", nst_core, k_core, l_majcore,"type 1"
           else
              print*,"issues"
              print*,"PO: st_ind, k, l_majcore", nst_core, k_core, l_majcore
              print*,"nk1(l_majcore)", nk1(l_majcore)
              stop
           end if
        end if

        if (inc == 1) then
           nst_core = n_core
        else 
           nst_core = TargetStates2el%ncore(n_core)
        end if 
!!$ is_core_MO is for the H2 Target State basis after constructing states
        k_core = get_inum_st(TargetStates%b(nst_core))
        l_majcore = get_l_majconf(TargetStates%b(nst_core))
        if (l_ion_core >= 0) then
           if ( k_core <= n_ion_core(l_majcore) ) then
              is_core_MO(nst_core) = 2
!!$              print*,"MO: st_ind, k, l", nst_core, k_core, l_majcore,"type 2"
           else if ( k_core <= nk1(l_majcore) ) then
              is_core_MO(nst_core) = 1
!!$              print*,"MO: st_ind, k, l", nst_core, k_core, l_majcore,"type 1"
           else
              print*,"issues"
              print*,"MO: st_ind, k, l", nst_core, k_core, l_majcore
              stop
           end if
        end if

     end do ! n_core - loop over core obitals
  end if ! non-uniqueness?


end subroutine structure12





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine: structure12group
!Purpose: custom version of the structure12 subroutine for use in H3+ mode.
!         Creates and diagonalises 2e hamiltonian for non-linear molecule,
!         using the hybrid basis produced in construct_1el_basis_nr_group
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine structure12group(basis,oneestates,num_states,indata)
    use numbers
    use basismodule
    use input_data      !data_in
    use sturmian_class  
    use state_class     !states data type
    use target_states   !Contains TargetStates2el, state_class_module.f90
		use ovlpste1me      !Contains overlap matrix of hybrid basis 
		use one_electron_func    !Contains e1me, 1e matrix elements in hybrid basis, bst_nr
    use target_states   !Contains global variables for storing 1 and 2e target states
    use MPI_module      !Contains myid global variable
    implicit none
    type(basis_sturmian_nr)::basis   !Sturmian basis used for 1e diagonalisation
    type(basis_state):: oneestates   !One electron spatial molecular orbitals
    type(smallinput):: indata
    integer:: ii, jj, counter, kk
    integer:: num_states !Size of one electron basis
    logical:: hlike
    integer:: maxstates
    !Variables controlling the configurations used
    integer:: l_ion_core
    integer, dimension(:), allocatable:: n_ion_core, nk1, nk2
    !List of configurations
    integer:: is !Spin
    integer:: numcon
    integer:: par
    real(dpf):: mval
    integer, dimension(:), allocatable:: phase
    integer:: mode     ! basis_type
    type(state),pointer :: tempstate
    !repo1, repo2 replace mo1, mo2 as m proj goes to group irrep in
    !non-linear case
    integer, dimension(:), allocatable:: no1, repo1, no2, repo2
    integer:: l12max
    !Configuration indices and matrix elements
    real(dpf), dimension(:,:), allocatable:: H, b
    real(dpf):: Helement, belement
    integer:: nc, ncp
    integer:: nsp1, nsp2, nsp1p, nsp2p 
    integer:: repnsp1, repnsp2,repnsp1p,repnsp2p
    integer:: num_4econfigs, numloops
    integer, dimension(:), allocatable:: nclist, ncplist
    integer:: rep, val
    !Required inputs for lapack routine dsygvx
    real*8:: abstol
    real*8, dimension(:), allocatable:: w, work
    real*8, dimension(:,:), allocatable:: ci
    integer, dimension(:), allocatable:: iwork, ifail
    integer:: info, lwork, nfound
    real*8, external :: DLAMCH
    !File IO variables
    character(len=100):: line
    !For checking values
    integer:: bfind 
    real*8:: tmpCI2
    type(sturmian_nr), pointer:: p2 
		integer:: Nmax1el
			
    !do ii = 1, oneestates%Nmax
    !   print*, oneestates%b(ii)%energy
    !end do

    !counter = 1
    !do ii = 1, oneestates%Nmax
		!   print*, oneestates%b(ii)%energy
    !  ! print*, "NEW BASIS SIZE: ", oneestates%b(ii)%nam
    !  ! do jj = 1, oneestates%b(ii)%nam
    !  !    print*, basis%b(counter)%k, basis%b(counter)%l, basis%b(counter)%m
    ! 	!	 counter = counter + 1
    ! 	!end do 
    !end do
    !stop
    
    !l12max = data_in%l12max

    maxstates=0 
    do mode = 0, 1  !mode=0 (count states), mode=1 (diagonisalise H12)
       if (mode .eq. 1) then
          !maxstates counted when mode=0
          hlike = .false.
          call new_basis_st(TargetStates2el,maxstates,hlike,basis_type)
       end if
    
       do is = 0, 1 !Spin loop
          !Call config12 to set up list of configurations to include: no1, mo1, no2, mo2 arrays
          numcon = oneestates%Nmax**2
          allocate(no1(numcon), repo1(numcon))
          allocate(n_ion_core(numcon))
          no1(:) = 0
          repo1(:) = 0
          numcon = 0

          l_ion_core = 0         
          n_ion_core(:) = 0

          !Call to find number of configurations to use
          call config12_tmp_st_group(oneestates, is, l12max, l_ion_core, n_ion_core, numcon, no1, repo1)
          if (mode .eq. 0) then !Count basis size
             maxstates=maxstates + numcon
          end if

          if(allocated(no1)) deallocate(no1,repo1)
          if(allocated(no2)) deallocate(no2,repo2) !coded like this for a reason

          if (mode .eq. 1) then !Perform structure
             allocate(no1(numcon),repo1(numcon),no2(numcon),repo2(numcon))
             no1(:) = 0
             repo1(:) = 0 
             no2(:) = 0
             repo2(:) = 0
             !Create list of configurations 
             rep = 3  !Index of irrep of symmetry group 

             !Produces lists no1, no2 of indices of one-electron states to use
             !to build two electron configs based on symmetry rules, etc.
             call config12_st_group(oneestates, rep,  is, l12max, l_ion_core, n_ion_core, numcon, no1, repo1, no2, repo2)

             !Fill arrays with indices of pairs of 2e configs
             num_4econfigs = numcon*(numcon+1)/2
             allocate(nclist(num_4econfigs), ncplist(num_4econfigs))
             nclist(:) = 0
             ncplist(:) = 0

             !Setup nc and ncp lists to index the set of PAIRS of 2e configurations
             numloops = 0
             do ii=1, numcon
                do jj = ii, numcon
                   numloops = numloops + 1
                   nclist(numloops) = ii
                   ncplist(numloops) = jj
                end do
             end do

             if (myid==0) write(*,*) 'Spin = ', is
             if (myid==0) write(*,*) 'Number of 2e configurations = ', numcon

             !Declare H and b matrices
             allocate(H(numcon,numcon), b(numcon,numcon))
             
             H(:,:) = 0.0_dpf
             b(:,:) = 0.0_dpf
             !do ii = 1, oneestates%Nmax
             !   nsp1 = get_nam(oneestates%b(ii))
             !   print*, oneestates%b(ii)%energy, nsp1
             !   do jj = 1, nsp1
             ! 	    val = get_na(oneestates%b(ii),jj,1)
             !      print*, val, basis%b(val)%l, basis%b(val)%m
             !   end do
             !end do
             !stop
             
             !Calculate H12 interaction matrix elements <pn|H|n>, 
             !|n> = |e1> * |e2> and |np> = |e1p> * |e2p>
             !$OMP PARALLEL DO DEFAULT(SHARED) & 
             !$OMP PRIVATE(nc,ncp,nsp1,repnsp1,nsp2,repnsp2,nsp1p,repnsp1p,nsp2p,repnsp2p,Helement,belement) &
             !$OMP SCHEDULE(DYNAMIC)
             do jj = 1, num_4econfigs
                nc = nclist(jj)
                ncp = ncplist(jj)

                !Indices of 1e molcule orbitals used in 2e configurations |n> and |np> for a given spin
                nsp1 = no1(nc) 
                nsp2 = no2(nc) 
                nsp1p = no1(ncp) 
                nsp2p = no2(ncp) 

                repnsp1 = repo1(nc) 
                repnsp2 = repo2(nc) 
                repnsp1p = repo1(ncp) 
                repnsp2p = repo2(ncp) 

                !call H12me_st_group(is, indata, oneestates, basis, nsp1, nsp2, nsp1p, nsp2p, repnsp1,repnsp2,repnsp1p,repnsp2p, Helement, belement)
								Nmax1el = oneestates%Nmax
                !call H12me_st_group_notortog(is, indata, oneestates, Nmax1el,e1me,ovlpst,basis, nsp1, nsp2, nsp1p, nsp2p, repnsp1,repnsp2,repnsp1p,repnsp2p, Helement, belement)
                call H12me_st_group_notortog(is, indata, oneestates, Nmax1el,e1me,ovlpst,bst_nr, nsp1, nsp2, nsp1p, nsp2p, repnsp1,repnsp2,repnsp1p,repnsp2p, Helement, belement)

                H(nc, ncp) = Helement
                H(ncp, nc) = Helement
                b(nc, ncp) = belement
                b(ncp, nc) = belement
                !print*, nc, ncp, Helement, belement
             end do
             !$OMP END PARALLEL DO

             !if ((is .eq. 0) .and. (numcon .lt. 5)) then
             !   open(90,file='HMatReese')
             !   write(90,*) "HMat"
             !   do kk = 1, numcon
             !      write(90,*) H(kk,:)
             !   end do
             !   write(90,*) "bmat"
             !   do kk = 1, numcon
             !      write(90,*) b(kk,:)
             !   end do
             !   close(90)
             !end if 

						 !print*, "H Matrix"
						 !do ii = 1, numcon
						 ! 	print*, H(ii,:)
						 !end do
						 !print*, "b Matrix"
						 !do ii = 1, numcon
						 ! 	print*, b(ii,:)
						 !end do

             !Diagonalise two electron hamiltonian to obtain H3+ electronic states 
             allocate(work(1))
             allocate(w(numcon))    !Stores eigenvalues
             allocate(iwork(5*numcon))
             allocate(ifail(numcon))
             allocate(ci(numcon,numcon))
             lwork = -1
             abstol = 2.0_dpf*DLAMCH('S') !twice underflow threshold for doubles
             info = 0
             nfound = 0
              		 
             !Workspace query with lwork=-1
             call DSYGVX(1, 'V', 'I', 'U', numcon, H, numcon, b, numcon,    &
             0.0_dpf, 0.0_dpf, 1, numcon, abstol, nfound, w, ci, numcon, work, lwork, &
             iwork, ifail, info)

             lwork = int(work(1))
             deallocate(work)
             allocate(work(lwork))

             !Perform diagonialisation
             call DSYGVX(1, 'V', 'I', 'U', numcon, H, numcon, b, numcon,    &
             0.0_dpf, 0.0_dpf, 1, numcon, abstol, nfound, w, ci, numcon, work, lwork, &
             iwork, ifail, info) 
             deallocate(work, iwork, ifail)

             !Save results of diagonalisation in TargetStates2el
             if (allocated(phase)) then
                deallocate(phase)
             end if
             allocate(phase(numcon))
             phase(:) = (-1)**is
             do ii = 1, numcon
                !Resolved sign problem
                if ( sum(ci(:,ii)) .lt. 0.0_dpf) then
                   ci(:,ii) = -ci(:,ii)
                end if

                TargetStates2el%Nstates = TargetStates2el%Nstates + 1
                tempstate => TargetStates2el%b(TargetStates2el%Nstates)
                mval =  0.0_dpf  !m is not a good quantum number, set to zero
                par =   0  !Neither is parity 
                call construct_st(tempstate,hlike,mval,par,dble(is),w(ii),ii,numcon,ci(:,ii),no1,repo1,no2,repo2,phase)
                tempstate%label = '-' 
             end do
             deallocate(phase)
             deallocate(no1, repo1, no2, repo2)
             deallocate(nclist, ncplist, H, b)
             deallocate(w,ci)
          end if
          deallocate(n_ion_core)
       end do !is loop
    end do  !mode loop

    call sort_by_energy_basis_st(TargetStates2el)

    open(81,file='2eenergies.txt')
    write(81,*) "  Two Electron States, Target Name: H3+"
    write(81,*) "  Nuclear Geometry:                   (R1, R2, R3)= ", indata%R(1), indata%R(2), indata%R(3)
    write(81,*) "                          (theta1, theta2, theta3)= ", indata%theta(1), indata%theta(2), indata%theta(3)
    write(81,*) "                                (phi1, phi2, phi3)= ", indata%phi(1), indata%phi(2), indata%phi(3)
    write(81,*) "  N     Label       S     E(a.u)"
    do ii = 1, min(TargetStates2el%Nstates,200)
       write(line,"(2X,I6,6X,A6,5X,I1,3X,f18.10)") ii, TargetStates2el%b(ii)%label, int(TargetStates2el%b(ii)%spin), TargetStates2el%b(ii)%energy
       write(81,*) ADJUSTL(TRIM(line))
    end do

    close(81)

    
end subroutine structure12group




!-----------------------------------------------------------------
subroutine config12(bstF5,ma,ip,is,l12max,nk1,nk2,l_ion_core,n_ion_core,ncm,no1,mo1,no2,mo2)

  use sturmian_class
  use MPI_module 

  implicit none

  type(basis_sturmian_nr), intent(in):: bstF5
  integer, intent(in):: ma, ip, is
  integer, intent(in):: l12max
  integer, dimension(0:l12max), intent(in):: nk1, nk2
  integer, intent(in):: l_ion_core
  integer, dimension(0:l12max), intent(in):: n_ion_core
  integer, intent(inout):: ncm
  integer, dimension(ncm),intent(inout):: no1, mo1, no2, mo2

  integer ico, nsp1, nsp2, nspm, lnsp1, mnsp1, lnsp2, mnsp2, icheck
  integer:: k1,k2
 
  nspm =  basis_size(bstF5) ! number of sturmian finctions

  ico = 0
  do nsp1=1,nspm
     lnsp1 =  get_ang_mom(bstF5%b(nsp1))
     do mnsp1 = -lnsp1,lnsp1
        
        do nsp2=nsp1,nspm
           lnsp2 =  get_ang_mom(bstF5%b(nsp2))
           if((-1)**(lnsp1+lnsp2) .ne. ip) cycle
           
           do mnsp2 = -lnsp2,lnsp2

!!$--------------------- logic block: include or not config ------------
              if(mnsp1+mnsp2 .ne. ma) cycle
              if(nsp1 .eq. nsp2) then
                 if(mnsp1 .gt. mnsp2) cycle  ! to avoid the same config. counted twice
                 if(mnsp1 .eq. mnsp2) then
                    if(is .eq. 1) cycle  ! symmetry condition
                 endif
              endif
              
              k1 = get_k(bstF5%b(nsp1))   ! inner orbital
              k2 = get_k(bstF5%b(nsp2))
              if(k1 .gt. nk1(lnsp1) ) cycle
              if(k2 .gt. nk2(lnsp2) ) cycle

              if(lnsp1 .gt. l12max) cycle
              if(lnsp2 .gt. l12max) cycle

              if(l_ion_core .ge. 0) then
                 if( k2 .le. nk1(lnsp2)  .or. k1 .le. n_ion_core(lnsp1)) then
                    continue ! include this configuration
                 else
                    cycle   ! ! exclude this configuration
                 endif
              endif


              if(ncm .gt. 0) then
                 icheck = 1
                 call testsameconfig(ncm,no1,no2,mo1,mo2,ico,nsp1,mnsp1,nsp2,mnsp2,icheck)
                 if(icheck .eq. 0) cycle 
              endif

!              if (myid==0) print*, nk1(lnsp1), k1,  lnsp1, k2, lnsp2

!!$--------------------- end logic block: include or not config ------------
              
              ico = ico + 1
              
              if(ncm .gt. 0) then
                 no1(ico) = nsp1
                 mo1(ico) = mnsp1
                 no2(ico) = nsp2
                 mo2(ico) = mnsp2
!                 if (myid==0) print*, nsp1, nsp2, mnsp1, mnsp2
              endif
              
           enddo
        enddo
     enddo
  enddo
  
  ncm = ico
     
return
end subroutine config12
!
!
subroutine config12_tmp(bstF5,ma,ip,is,l12max,nk1,nk2,l_ion_core,n_ion_core,ncm)

  use sturmian_class
  use MPI_module
  
  implicit none

  type(basis_sturmian_nr), intent(in):: bstF5
  integer, intent(in):: ma, ip, is
  integer, intent(in):: l12max! l12max
  integer, dimension(0:l12max), intent(in):: nk1, nk2
  integer, intent(in):: l_ion_core
  integer, dimension(0:l12max), intent(in):: n_ion_core
  integer, intent(inout):: ncm

  integer ico, nsp1, nsp2, nspm, lnsp1, mnsp1, lnsp2, mnsp2, icheck
  integer:: k1,k2
 
  nspm =  basis_size(bstF5) ! number of sturmian finctions

  ico = 0
  do nsp1=1,nspm
     lnsp1 =  get_ang_mom(bstF5%b(nsp1))
     do mnsp1 = -lnsp1,lnsp1
        
        do nsp2=nsp1,nspm
           lnsp2 =  get_ang_mom(bstF5%b(nsp2))
           if((-1)**(lnsp1+lnsp2) .ne. ip) cycle
           
           do mnsp2 = -lnsp2,lnsp2

!!$--------------------- logic block: include or not config ------------
              if(mnsp1+mnsp2 .ne. ma) cycle
              if(nsp1 .eq. nsp2) then
                 if(mnsp1 .gt. mnsp2) cycle  ! to avoid the same config. counted twice
                 if(mnsp1 .eq. mnsp2) then
                    if(is .eq. 1) cycle  ! symmetry condition
                 endif
              endif
              
              k1 = get_k(bstF5%b(nsp1))  ! inner orbital
              k2 = get_k(bstF5%b(nsp2))
              if(k1 .gt. nk1(lnsp1) ) cycle
              if(k2 .gt. nk2(lnsp2) ) cycle

              if(lnsp1 .gt. l12max) cycle
              if(lnsp2 .gt. l12max) cycle

              if(l_ion_core .ge. 0) then
                 if( k2 .le. nk1(lnsp2)  .or. k1 .le. n_ion_core(lnsp1)) then
                    continue ! include this configuration
                 else
                    cycle   ! ! exclude this configuration
                 endif
              endif
              
!!$--------------------- end logic block: include or not config ------------
              
              ico = ico + 1
              
           enddo
        enddo
     enddo
  enddo
  
  ncm = ico
     
return
end subroutine config12_tmp
!
!
subroutine config12_tmp_st(TargetStates1el,ma,ip,is,l12max,nk1,nk2,l_ion_core,n_ion_core,ncm,no1,mo1)
  !
  ! Estimate the number of configurations (ncm) for this symmetry (ip,ma,is).
  ! Not all combinations of orbitals are used! They are selected based on:
  !   - l_ion_core, n_ion_core: usually [0 1].
  !   - nk1: usually larger, e.g. [3 2 1].
  !
  ! Configurations can be grouped into two types:
  !   - Ionic: For a few number of r1 orbitals, use ALL r2 orbitals.
  !   - CI: For a larger number of r1 orbitals, use a few r2 orbitals.
  !
  ! Typically, we have l_ion_core=0 and n_ion_core=[1] which gives us
  !   - 1s+1s, 1s+2s, ... , 1s+Ns   with N=nk2(0)
  !   - 1s+2p, 1s+3p, ... , 1s+Np   with N=nk2(1)
  !   -   :      :     :      :
  !   - etc., subject to (-1)^(l1+l2)=ip and m1+m2=ma selection rules.
  ! Then with e.g. nk1=[3 2 1] we also have the full CI configurations of
  !   - 2s+2s, 2s+3s, 2s+2p, 2s+3p, 2s+3d
  !   -        3s+3s, 3s+2p, 3s+3p, 3s+3d
  !   -               2p+2p, 2p+3p, 2p+3d
  !   -                      3p+3p, 3p+3d
  !   -                             3d+3d
  !   - again subject to (-1)^(l1+l2)=ip and m1+m2=ma selection rules.
  !
  use state_class
  use MPI_module

  implicit none

  type(basis_state), intent(in):: TargetStates1el
  integer, intent(in):: ma, ip, is
  integer, intent(in):: l12max
  integer, dimension(0:l12max), intent(in):: nk1, nk2
  integer, intent(in):: l_ion_core
  integer, dimension(0:l12max), intent(in):: n_ion_core
  integer, intent(inout):: ncm
  integer, dimension(TargetStates1el%Nmax*TargetStates1el%Nmax),intent(inout):: no1, mo1

  integer, dimension(TargetStates1el%Nmax*TargetStates1el%Nmax):: no2, mo2
  integer:: Nmax, ico, nst1, nst2, ipar1, ipar2, mnst1, mnst2, icheck
  integer:: k1,k2, l1_majconf, l2_majconf, n1_majconf, n2_majconf

  Nmax =  basis_size_st(TargetStates1el) ! number of  one-el. target states

  ico = 0
  do nst1=1,Nmax
     ipar1 = get_par_st(TargetStates1el%b(nst1))     
     mnst1 = get_ang_mom_proj(TargetStates1el%b(nst1))
     if(abs(mnst1) .gt. l12max) cycle
        
     l1_majconf = get_l_majconf(TargetStates1el%b(nst1))
     n1_majconf = get_n_majconf(TargetStates1el%b(nst1))

     do nst2=1,Nmax
        ipar2 = get_par_st(TargetStates1el%b(nst2))
        mnst2 = get_ang_mom_proj(TargetStates1el%b(nst2))

        if(abs(mnst2) .gt. l12max) cycle
        ! If homogeneous (ip=+/-1) config parity must equal symmetry parity.
        ! If heterogeneous (ip=0) there is no parity.
        if(ip/=0 .and. ipar1*ipar2/=ip) cycle

!!$--------------------- logic block: include or not config ------------
        if(mnst1+mnst2 .ne. ma) cycle
        if(nst1 .eq. nst2) then
           if(is .eq. 1) cycle  ! symmetry condition
        endif
        
        l2_majconf = get_l_majconf(TargetStates1el%b(nst2))
        n2_majconf = get_n_majconf(TargetStates1el%b(nst2))

        k1 = n1_majconf - l1_majconf
        k2 = n2_majconf - l2_majconf
        !k1 = get_inum_st(TargetStates1el%b(nst1))
        !k2 = get_inum_st(TargetStates1el%b(nst2))

        if(k1 .gt. nk1(abs(mnst1)) ) cycle  ! inner orbital
        if(k2 .gt. nk2(abs(mnst2)) ) cycle

        if(l1_majconf .gt. l12max) cycle
        if ( k1 > nk1(l1_majconf) ) cycle   ! Restrict r1 angular momentum.
        if(l2_majconf .gt. l12max) cycle

        if(l_ion_core .ge. 0) then
           ! First test if an ionic configuration.
           if ( k1 > n_ion_core(l1_majconf) ) then   ! Not ionic.
              ! Second test if a full CI configuration.
              if ( k2 > nk1(l2_majconf) ) cycle   ! Not full CI either.
           endif
        endif

        ! Remove any (antisymmetric) duplicates.
        if(ico .gt. 0) then
           icheck = 1
           call testsameconfig(ico,no1,no2,mo1,mo2,ico,nst1,mnst1,nst2,mnst2,icheck)
           if(icheck .eq. 0) cycle
        endif        
!!$--------------------- end logic block: include or not config ------------

        ico = ico + 1
!!$ This is specified to count the number of core orbitals
        no1(ico) = nst1
        mo1(ico) = mnst1
        no2(ico) = nst2
        mo2(ico) = mnst2
     enddo
  enddo

  ncm = ico

  return
end subroutine config12_tmp_st
!
!
subroutine config12_st(TargetStates1el,ma,ip,is,l12max,nk1,nk2,l_ion_core,n_ion_core,ncm,no1,mo1,no2,mo2)
  !
  ! See the previous subroutine for information on these configuration.  !
  use state_class
  use MPI_module
  use input_data

  implicit none

  type(basis_state), intent(in):: TargetStates1el
  integer, intent(in):: ma, ip, is
  integer, intent(in):: l12max
  integer, dimension(0:l12max), intent(in):: nk1, nk2
  integer, intent(in):: l_ion_core
  integer, dimension(0:l12max), intent(in):: n_ion_core
  integer, intent(inout):: ncm
  integer, dimension(ncm),intent(inout):: no1, mo1, no2, mo2
  character*8, parameter :: lchars='spdfghiklmnoqrtuv'
  character*8, parameter :: Mchars='SPDFGHIKLMNOQRTUV'
  character*3, parameter :: pchars='u g'
  character*8 :: orb_type

  logical :: ex

  integer:: ico, Nmax, nst1, nst2, ipar1, ipar2, mnst1, mnst2, icheck
  integer:: k1,k2,  l1_majconf, l2_majconf, n1_majconf, n2_majconf

  Nmax =  basis_size_st(TargetStates1el) ! number of one-el. target states

  if(myid==0 .and. data_in%print_1el_basis) then
    open(unit=666,file='one-el-basis.out',action='write',status='replace')
    write(666,'(A)') 'index, nlm, m, p, orbital type'
  endif

  if(myid==0 .and. data_in%print_2el_config.and. ma >=0 ) then
    write(777,'("===       SYMMETRY: ",I0,2A,"      ===")')  2*is+1, Mchars(abs(ma)+1:abs(ma)+1), pchars(ip+2:ip+2)
    write(777,'(A)') 'index, nlm(1), nlm(2),    m1, m2'
  endif

  ico = 0
  do nst1=1,Nmax
     ipar1 =  get_par_st(TargetStates1el%b(nst1))
     mnst1 = get_ang_mom_proj(TargetStates1el%b(nst1))
     l1_majconf = get_l_majconf(TargetStates1el%b(nst1))
     n1_majconf = get_n_majconf(TargetStates1el%b(nst1))
     if(myid==0 .and. data_in%print_1el_basis) then
       if(get_energy(TargetStates1el%b(nst1)) == 0.0d0) then
         orb_type = 'Laguerre'
         write(666,'(I5.1,2X,I1,2A, 2X, I2.1, 2X, I2.1, 2X, A8)') nst1, n1_majconf, lchars(l1_majconf+1:l1_majconf+1), Mchars(abs(mnst1)+1:abs(mnst1)+1), mnst1, get_par(TargetStates1el%b(nst1)), orb_type
         !write(666,'(I5.1,2X,A4 X, I2.1, 2X, A8)') nst1, trim(adjustl(get_label(TargetStates1el%b(nst1)))), mnst1, orb_type
       else
         orb_type = '   MO   '
         write(666,'(I5.1,2X,A4, X, I2.1, 2X, I2.1, 2X, A8)') nst1, trim(adjustl(get_label(TargetStates1el%b(nst1)))), mnst1, get_par(TargetStates1el%b(nst1)), orb_type
       endif
     endif

     if(abs(mnst1) .gt. l12max) cycle
     do nst2=1,Nmax
        ipar2 =  get_par_st(TargetStates1el%b(nst2))
        mnst2 = get_ang_mom_proj(TargetStates1el%b(nst2))

        if(abs(mnst2) .gt. l12max) cycle
        ! If homogeneous (ip=+/-1) config parity must equal symmetry parity.
        ! If heterogeneous (ip=0) there is no parity.
        if(ip/=0 .and. ipar1*ipar2/=ip) cycle

!!$--------------------- logic block: include or not config ------------
        if(mnst1+mnst2 .ne. ma) cycle
        if(nst1 .eq. nst2) then
           if(is .eq. 1) cycle  ! symmetry condition
        endif
        
        l2_majconf = get_l_majconf(TargetStates1el%b(nst2))
        n2_majconf = get_n_majconf(TargetStates1el%b(nst2))
        
        k1 = n1_majconf - l1_majconf
        k2 = n2_majconf - l2_majconf

        !k1 = get_inum_st(TargetStates1el%b(nst1))
        !k2 = get_inum_st(TargetStates1el%b(nst2))

        if(k1 .gt. nk1(abs(mnst1)) ) cycle   !  inner orbital
        if(k2 .gt. nk2(abs(mnst2)) ) cycle
        
        if(l1_majconf .gt. l12max) cycle
        if ( k1 > nk1(l1_majconf) ) cycle
        if(l2_majconf .gt. l12max) cycle
        if ( k2 > nk2(l2_majconf) ) cycle

        if(l_ion_core .ge. 0) then
           if ( k1 > n_ion_core(l1_majconf) ) then
              if ( k2 > nk1(l2_majconf) ) cycle
           endif
        endif

        ! Remove any (antisymmetric) duplicates.
        if(ncm .gt. 0) then
           icheck = 1
           call testsameconfig(ncm,no1,no2,mo1,mo2,ico,nst1,mnst1,nst2,mnst2,icheck)
           if(icheck .eq. 0) cycle
        endif
!!$--------------------- end logic block: include or not config ------------

        ico = ico + 1
  
        if(myid==0 .and. data_in%print_2el_config .and. ma >= 0) write(777,'(I5,2X,2(I1,2A,5X),3X,2(I2.1,2X), 2(A6))') ico, &
          &n1_majconf, lchars(l1_majconf+1:l1_majconf+1), Mchars(abs(mnst1)+1:abs(mnst1)+1), &
          &n2_majconf, lchars(l2_majconf+1:l2_majconf+1), Mchars(abs(mnst2)+1:abs(mnst2)+1), mnst1, mnst2, &
          &TargetStates1el%b(nst1)%label, TargetStates1el%b(nst2)%label

        if(ncm .gt. 0) then
           no1(ico) = nst1
           mo1(ico) = mnst1
           no2(ico) = nst2
           mo2(ico) = mnst2
!          if (myid==0) print'("st: ",8i5)', nst1, nst2, k1, k2, mnst1, mnst2, get_l_majconf(TargetStates1el%b(nst1)), get_l_majconf(TargetStates1el%b(nst2))
! write(*,'(3(3I3,3X))') ico,nst1,nst2, k1,l1_majconf,mnst1, k2,l2_majconf,mnst2
        endif

     enddo
  enddo

  ncm = ico

  if(myid==0 .and. data_in%print_1el_basis) close(666)
  if(myid==0 .and. data_in%print_2el_config .and. ma >= 0) write(777,*)

return
end subroutine config12_st
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine: config12_tmp_st_group
!Purpose: calculates the number of configurations that will be used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine config12_tmp_st_group(TargetStates1el,is,l12max,l_ion_core,n_ion_core,ncm,no1,repo1)
  !
	! Estimate the number of configurations (ncm) to include for this spin (is)
  ! Not all combinations of orbitals are used! They are selected based on:
  !   - l_ion_core, n_ion_core: usually [0 1].
  !   - nk1: usually larger, e.g. [3 2 1].
  !
  ! Configurations can be grouped into two types:
  !   - Ionic: For a few number of r1 orbitals, use ALL r2 orbitals.
  !   - CI: For a larger number of r1 orbitals, use a few r2 orbitals.
	!   - Frozen core: for a single (ground state) r1 orbital, use many r2
	!                  orbitals (special case of ionic)
  !
  use state_class
  use MPI_module

  implicit none

  type(basis_state), intent(in):: TargetStates1el
  integer, intent(in):: is
  integer, intent(in):: l12max
  integer, intent(in):: l_ion_core
  integer, dimension(0:l12max), intent(in):: n_ion_core
  integer, intent(inout):: ncm
  integer,dimension(TargetStates1el%Nmax*TargetStates1el%Nmax),intent(inout):: no1, repo1

  integer, dimension(TargetStates1el%Nmax*TargetStates1el%Nmax):: no2, repo2
  integer:: Nmax, ico, nst1, nst2, ipar1, ipar2, repnst1, repnst2, icheck
  integer:: k1,k2, l1_majconf, l2_majconf, n1_majconf, n2_majconf, m1_majconf, m2_majconf

  Nmax =  basis_size_st(TargetStates1el) ! number of  one-el. target states

  ico = 0
  do nst1=1,Nmax
     !mnst1 = get_ang_mom_proj(TargetStates1el%b(nst1))
	   repnst1 = -1 
	   m1_majconf = get_ang_mom_proj(TargetStates1el%b(nst1))
     l1_majconf = get_l_majconf(TargetStates1el%b(nst1))
     n1_majconf = get_n_majconf(TargetStates1el%b(nst1)) 

     do nst2=1,Nmax

        !mnst2 = get_ang_mom_proj(TargetStates1el%b(nst2))
!!$--------------------- logic block: include or not config ------------

				!For now just use frozen core approximation
				if (nst1 .gt. 1) then
				   cycle
			  end if

        !if(mnst1+mnst2 .ne. ma) cycle  !Replace with group representation condition
        if(nst1 .eq. nst2) then
           if(is .eq. 1) cycle  ! spin symmetry condition
        endif
        
				repnst2 = -1
	      m2_majconf = get_ang_mom_proj(TargetStates1el%b(nst2))
        l2_majconf = get_l_majconf(TargetStates1el%b(nst2))
        n2_majconf = get_n_majconf(TargetStates1el%b(nst2))

				!k1 and k2 are "approximate" quantum numbers used to enumerate one
				!electron m.o's in the diatomic case.
        !k1 = n1_majconf - l1_majconf
        !k2 = n2_majconf - l2_majconf
        !if(k1 .gt. nk1(abs(mnst1)) ) cycle  ! inner orbital
        !if(k2 .gt. nk2(abs(mnst2)) ) cycle
        !if(l1_majconf .gt. l12max) cycle
        !if(l2_majconf .gt. l12max) cycle

        ! Remove any (antisymmetric) duplicates.
        if(ico .gt. 0) then
           icheck = 1
           call testsameconfiggroup(ico,no1,no2,repo1,repo2,ico,nst1,repnst1,nst2,repnst2,icheck)
           if(icheck .eq. 0) cycle
        endif        
!!$--------------------- end logic block: include or not config ------------

        ico = ico + 1
!!$ This is specified to count the number of core orbitals
        no1(ico) = nst1
        repo1(ico) = repnst1
        no2(ico) = nst2
        repo2(ico) = repnst2
     enddo
  enddo

  ncm = ico

  return
end subroutine config12_tmp_st_group
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Subroutine: config12_st_group
!Purpose: produces a list of configurations to use 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine config12_st_group(TargetStates1el,rep,is,l12max,l_ion_core,n_ion_core,ncm,no1,repo1,no2,repo2)
  !
  ! See the previous subroutine for information on these configuration.  !
  use state_class
  use MPI_module
  use input_data

  implicit none

  type(basis_state), intent(in):: TargetStates1el
  integer, intent(in):: rep, is
  integer, intent(in):: l12max
  integer, intent(in):: l_ion_core
  integer, dimension(0:l12max), intent(in):: n_ion_core
  integer, intent(inout):: ncm
  integer, dimension(ncm),intent(inout):: no1, repo1, no2, repo2
  character*8, parameter :: lchars='spdfghiklmnoqrtuv'
  character*8, parameter :: Mchars='SPDFGHIKLMNOQRTUV'
  character*4, parameter, dimension(6) :: replabels=(/'A1p ', 'A2p ','Ep  ', 'A1pp', 'A2pp', 'Epp '/)
  character*8 :: orb_type

  logical :: ex

  integer:: ico, Nmax, nst1, nst2, repnst1, repnst2, icheck
  integer:: l1_majconf, l2_majconf, n1_majconf, n2_majconf, m1_majconf, m2_majconf
	integer:: mval, m1val, m2val

  Nmax =  basis_size_st(TargetStates1el) ! number of one-el. target states

  if(myid==0 .and. data_in%print_1el_basis) then
    open(unit=666,file='one-el-basis.out',action='write',status='replace')
    write(666,'(A)') 'index, nlm, m, p, orbital type'
  endif

  if(myid==0 .and. data_in%print_2el_config.and. rep >=0 ) then
    write(777,'("===       SYMMETRY: ",I0,A,"      ===")')  2*is+1, replabels(rep)
    write(777,'(A)') 'index, nlm(1), nlm(2),    rep1,  rep2'
  endif

  ico = 0
  do nst1=1,Nmax
     repnst1 = -1
     m1_majconf = get_ang_mom_proj(TargetStates1el%b(nst1))
     l1_majconf = get_l_majconf(TargetStates1el%b(nst1))
     n1_majconf = get_n_majconf(TargetStates1el%b(nst1))

     if(myid==0 .and. data_in%print_1el_basis) then
        if(get_energy(TargetStates1el%b(nst1)) == 0.0d0) then
           orb_type = 'Laguerre'
	         mval = m1_majconf
           write(666,'(I5.1,2X,I1,2A, 2X, I2.1, 2X, I2.1, 2X, A8)') nst1, n1_majconf, lchars(l1_majconf+1:l1_majconf+1), &
	               Mchars(abs(mval)+1:abs(mval)+1), mval, get_par(TargetStates1el%b(nst1)), orb_type 
           !write(666,'(I5.1,2X,A4 X, I2.1, 2X, A8)') nst1, trim(adjustl(get_label(TargetStates1el%b(nst1)))), mnst1, orb_type
        else
           orb_type = '   MO   '
           write(666,'(I5.1,2X,A4, X, I2.1, 2X, I2.1, 2X, A8)') nst1, trim(adjustl(get_label(TargetStates1el%b(nst1)))), m1_majconf, get_par(TargetStates1el%b(nst1)), orb_type
        endif
     endif

     do nst2=1,Nmax
   !!$--------------------- logic block: include or not config ------------
      	 !For now just use frozen core approximation
        if (nst1 .gt. 1) then
           cycle
        end if
       
        !if(mnst1+mnst2 .ne. ma) cycle  !Replace with group representation condition
        if(nst1 .eq. nst2) then
           if(is .eq. 1) cycle  !spin symmetry condition
        endif
           
   	    repnst2 = -1
   	    m2_majconf = get_ang_mom_proj(TargetStates1el%b(nst2))
        l2_majconf = get_l_majconf(TargetStates1el%b(nst2))
        n2_majconf = get_n_majconf(TargetStates1el%b(nst2))
        
        ! Remove any (antisymmetric) duplicates.
        if(ncm .gt. 0) then
           icheck = 1
           call testsameconfiggroup(ncm,no1,no2,repo1,repo2,ico,nst1,repnst1,nst2,repnst2,icheck)
           if(icheck .eq. 0) cycle
        endif
   !!$------------------ end logic block: include or not config ------------
   
        ico = ico + 1
     
   	    m1val = m1_majconf 
   	    m2val = m2_majconf 
        if(myid==0 .and. data_in%print_2el_config .and. rep >= 0) write(777,'(I5,2X,2(I1,2A,5X),3X,2(I2.1,2X), 2(A6))') ico, &
          &n1_majconf, lchars(l1_majconf+1:l1_majconf+1), Mchars(abs(m1val)+1:abs(m1val)+1), &
          &n2_majconf, lchars(l2_majconf+1:l2_majconf+1), Mchars(abs(m2val)+1:abs(m2val)+1), replabels(rep), replabels(rep), &
          &TargetStates1el%b(nst1)%label, TargetStates1el%b(nst2)%label
      
        if(ncm .gt. 0) then
           no1(ico) = nst1
           repo1(ico) = repnst1
           no2(ico) = nst2
           repo2(ico) = repnst2
        endif
     enddo
  enddo

  if (ico .ne. ncm) then
     print*, "ERROR: number of configurations chosen does not match earlier tally in config12_tmp_st_group, stopping"
	   stop
  end if

  ncm = ico

  if(myid==0 .and. data_in%print_1el_basis) close(666)
  if(myid==0 .and. data_in%print_2el_config .and. rep >= 0) write(777,*)

return
end subroutine config12_st_group
!
!
subroutine config12_tmp_st_old(TargetStates1el,ma,ip,is,l12max,nk1,nk2,l_ion_core,n_ion_core,ncm)

  use state_class
  use MPI_module
  
  implicit none

  type(basis_state), intent(in):: TargetStates1el
  integer, intent(in):: ma, ip, is
  integer, intent(in):: l12max
  integer, dimension(0:l12max), intent(in):: nk1, nk2
  integer, intent(in):: l_ion_core
  integer, dimension(0:l12max), intent(in):: n_ion_core
  integer, intent(inout):: ncm

  integer Nmax, ico, nst1, nst2, ipar1, ipar2, mnst1, mnst2, icheck
  integer:: k1,k2, l1_majconf, l2_majconf
 
  Nmax =  basis_size_st(TargetStates1el) ! number of  one-el. target states
  
  ico = 0
  do nst1=1,Nmax
     ipar1 =  get_par_st(TargetStates1el%b(nst1))
     mnst1 = get_ang_mom_proj(TargetStates1el%b(nst1))
     if(abs(mnst1) .gt. l12max) cycle     

!     do nst2=1,Nmax
     do nst2=nst1,Nmax
        ipar2 =  get_par_st(TargetStates1el%b(nst2))
        mnst2 = get_ang_mom_proj(TargetStates1el%b(nst2))
        if(abs(mnst2) .gt. l12max) cycle
        
        if(ipar1*ipar2 .ne. ip) cycle               

!!$--------------------- logic block: include or not config ------------
        if(mnst1+mnst2 .ne. ma) cycle
        if(nst1 .eq. nst2) then
           if(is .eq. 1) cycle  ! symmetry condition
        endif
        
        k1 = get_inum_st(TargetStates1el%b(nst1))
        k2 = get_inum_st(TargetStates1el%b(nst2))
        if(k1 .gt. nk1(abs(mnst1)) ) cycle  ! inner orbital
        if(k2 .gt. nk2(abs(mnst2)) ) cycle
       
        l1_majconf = get_l_majconf(TargetStates1el%b(nst1))
        l2_majconf = get_l_majconf(TargetStates1el%b(nst2))
        if(l1_majconf .gt. l12max) cycle
        if(l2_majconf .gt. l12max) cycle
        
        if(l_ion_core .ge. 0) then
           !           if (myid==0) print*, l1_majconf, k1, n1(l1_majconf), n_ion_core(l1_majconf)
           if( k2 .le. nk1(l2_majconf)  .or. k1 .le. n_ion_core(l1_majconf)) then
              continue ! include this configuration
           else
              cycle   ! ! exclude this configuration
           endif
        endif
        
!!$--------------------- end logic block: include or not config ------------
        
        ico = ico + 1
        
     enddo
  enddo
  
  ncm = ico
  
  return
end subroutine config12_tmp_st_old
!
!
subroutine config12_st_old(TargetStates1el,ma,ip,is,l12max,nk1,nk2,l_ion_core,n_ion_core,ncm,no1,mo1,no2,mo2)

  use state_class
  use MPI_module

  implicit none

  type(basis_state), intent(in):: TargetStates1el
  integer, intent(in):: ma, ip, is
  integer, intent(in):: l12max
  integer, dimension(0:l12max), intent(in):: nk1, nk2
  integer, intent(in):: l_ion_core
  integer, dimension(0:l12max), intent(in):: n_ion_core
  integer, intent(inout):: ncm
  integer, dimension(ncm),intent(inout):: no1, mo1, no2, mo2

  integer ico, Nmax, nst1, nst2, ipar1, ipar2, mnst1, mnst2, icheck
  integer:: k1,k2,  l1_majconf, l2_majconf

  Nmax =  basis_size_st(TargetStates1el) ! number of one-el. target states

  ico = 0
  do nst1=1,Nmax
     ipar1 =  get_par_st(TargetStates1el%b(nst1))
     mnst1 = get_ang_mom_proj(TargetStates1el%b(nst1))
     if(abs(mnst1) .gt. l12max) cycle
!     do nst2=1,Nmax
     do nst2=nst1,Nmax
        ipar2 =  get_par_st(TargetStates1el%b(nst2))
        mnst2 = get_ang_mom_proj(TargetStates1el%b(nst2))
        if(abs(mnst2) .gt. l12max) cycle

             if(ipar1*ipar2 .ne. ip) cycle               

!!$--------------------- logic block: include or not config ------------
        if(mnst1+mnst2 .ne. ma) cycle
        if(nst1 .eq. nst2) then
           if(is .eq. 1) cycle  ! symmetry condition
        endif
        
        k1 = get_inum_st(TargetStates1el%b(nst1))
        k2 = get_inum_st(TargetStates1el%b(nst2))
        if(k1 .gt. nk1(abs(mnst1)) ) cycle   !  inner orbital
        if(k2 .gt. nk2(abs(mnst2)) ) cycle

        l1_majconf = get_l_majconf(TargetStates1el%b(nst1))
        l2_majconf = get_l_majconf(TargetStates1el%b(nst2))
        if(l1_majconf .gt. l12max) cycle
        if(l2_majconf .gt. l12max) cycle

        if(l_ion_core .ge. 0) then
           if( k2 .le. nk1(l2_majconf)  .or. k1 .le. n_ion_core(l1_majconf)) then
              continue ! include this configuration
           else
              cycle   ! ! exclude this configuration
           endif
        endif        
        
        if(ncm .gt. 0) then
           icheck = 1
           call testsameconfig(ncm,no1,no2,mo1,mo2,ico,nst1,mnst1,nst2,mnst2,icheck)
           if(icheck .eq. 0) cycle 
        endif
!!$--------------------- end logic block: include or not config ------------
        
        ico = ico + 1

        if(ncm .gt. 0) then
           no1(ico) = nst1
           mo1(ico) = mnst1
           no2(ico) = nst2
           mo2(ico) = mnst2
!          if (myid==0) print'("st: ",8i5)', nst1, nst2, k1, k2, mnst1, mnst2, get_l_majconf(TargetStates1el%b(nst1)), get_l_majconf(TargetStates1el%b(nst2))
!          write(*,'(3(3I3,3X))') ico,nst1,nst2, k1,l1_majconf,mnst1, k2,l2_majconf,mnst2
        endif
        
     enddo
  enddo
  
  ncm = ico
  
return
end subroutine config12_st_old
!
!
subroutine testsameconfig(ncm,no1,no2,mo1,mo2,ico,nsp1,mnsp1,nsp2,mnsp2,icheck)

	integer:: ncm, n
  integer, dimension(ncm),intent(in):: no1, mo1, no2, mo2
  integer, intent(in):: nsp1,mnsp1,nsp2,mnsp2, ico
  integer, intent(inout):: icheck

  do n=1,ico

     if( (nsp1 .eq. no1(n) .and. mnsp1 .eq. mo1(n)) .and. (nsp2 .eq. no2(n) .and. mnsp2 .eq. mo2(n)) ) then
        icheck = 0
        exit
     endif

     if( (nsp2 .eq. no1(n) .and. mnsp2 .eq. mo1(n)) .and. (nsp1 .eq. no2(n) .and. mnsp1 .eq. mo2(n)) ) then
        icheck = 0
        exit
     endif

  enddo


  return
end subroutine testsameconfig



subroutine testsameconfiggroup(ncm,no1,no2,repo1,repo2,ico,nsp1,repnsp1,nsp2,repnsp2,icheck)

	integer:: ncm, n
  integer, dimension(ncm),intent(in):: no1, repo1, no2, repo2
  integer, intent(in):: nsp1,repnsp1,nsp2,repnsp2, ico
  integer, intent(inout):: icheck

  do n=1,ico
     if( (nsp1 .eq. no1(n) .and. repnsp1 .eq. repo1(n)) .and. (nsp2 .eq. no2(n) .and. repnsp2 .eq. repo2(n)) ) then
        icheck = 0
        exit
     endif

     if( (nsp2 .eq. no1(n) .and. repnsp2 .eq. repo1(n)) .and. (nsp1 .eq. no2(n) .and. repnsp1 .eq. repo2(n)) ) then
        icheck = 0
        exit
     endif

  enddo


  return
end subroutine testsameconfiggroup



!$**************************************************************************************************
!$$   This subroutine is to calculate matrix elements of H for H_2 molecule for 4 one-electron functions 
!$$   with given orbital angular momentum and its projection.
!$$   antisymmetric configurations 

subroutine H12me(Znuc,Rd,is,bstF5,nsp1,nsp2,nsp1p,nsp2p,mnsp1,mnsp2,mnsp1p,mnsp2p,resultH,resultb)
  
  use sturmian_class
  use MPI_module
  use input_data

  implicit none

  real*8, intent(in):: Znuc, Rd
  integer, intent(in):: is
  type(basis_sturmian_nr), intent(in):: bstF5
  integer, intent(in):: nsp1,nsp2,nsp1p,nsp2p
  integer, intent(in):: mnsp1,mnsp2,mnsp1p,mnsp2p
  real*8, intent(out):: resultH,resultb

  type(sturmian_nr), pointer:: p1, p2, p1p, p2p
  real*8:: oneelme, result,  twoelme, tmp
  real*8 :: Z1, Z2

  
  resultH = 0d0
  resultb = 0d0

  oneelme = 0d0
  twoelme = 0d0
  tmp = 0d0

  p1 => bstF5%b(nsp1)
  p2 => bstF5%b(nsp2)
  p1p => bstF5%b(nsp1p)
  p2p => bstF5%b(nsp2p)

  Z1 = data_in%Z1
  Z2 = data_in%Z2
  
  
!!$   overlap and one-electron operator matrix              
  if(mnsp1 .eq. mnsp1p .and. mnsp2 .eq. mnsp2p) then
     resultb =  bstF5%ortint(nsp1,nsp1p) * bstF5%ortint(nsp2,nsp2p)
     
     call Hlagorb(bstF5,nsp1,nsp1p,mnsp1,result)
     oneelme =  result * bstF5%ortint(nsp2,nsp2p)
     call Hlagorb(bstF5,nsp2,nsp2p,mnsp2,result)
     oneelme =  oneelme + result * bstF5%ortint(nsp1,nsp1p)
     
  endif  
  
  if((nsp1 .eq. nsp2 .and. mnsp1 .eq. mnsp2) .and. (nsp1p .ne. nsp2p .or. mnsp1p .ne. mnsp2p)) then
     ! <aa | ..| b1 b2>
     resultb = sqrt(2d0) * resultb 
     oneelme = sqrt(2d0) * oneelme
  elseif((nsp1 .ne. nsp2 .or. mnsp1 .ne. mnsp2) .and. (nsp1p .eq. nsp2p .and. mnsp1p .eq. mnsp2p)) then
     ! <a1a2 | ..| bb>
     resultb = sqrt(2d0) * resultb 
     oneelme = sqrt(2d0) * oneelme
  elseif((nsp1 .eq. nsp2 .and. mnsp1 .eq. mnsp2) .and. (nsp1p .eq. nsp2p .and. mnsp1p .eq. mnsp2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...  
  else
     !<a1 a2 | ..| b1 b2>
     if(mnsp1 .eq. mnsp2p .and. mnsp2 .eq. mnsp1p) then
        resultb = resultb + (-1)**(is) * bstF5%ortint(nsp1,nsp2p) * bstF5%ortint(nsp2,nsp1p)
        
        call Hlagorb(bstF5,nsp1,nsp2p,mnsp1,result)
        oneelme =  oneelme + (-1)**(is) * result * bstF5%ortint(nsp2,nsp1p)
        
        call Hlagorb(bstF5,nsp2,nsp1p,mnsp2,result)
        oneelme =  oneelme + (-1)**(is) * result * bstF5%ortint(nsp1,nsp2p)                       
     endif
  endif
  
  
!!$   two-electron operator  matrix
  call V12me(p1,p2,p1p,p2p,mnsp1,mnsp2,mnsp1p,mnsp2p,result)
  twoelme  = result
  !                if (myid==0)  print*, 'result', result
  if((nsp1 .eq. nsp2 .and. mnsp1 .eq. mnsp2) .and. (nsp1p .ne. nsp2p .or. mnsp1p .ne. mnsp2p)) then
     ! <aa | ..| b1 b2>
     twoelme  =  sqrt(2d0) * twoelme 
  elseif((nsp1 .ne. nsp2 .or. mnsp1 .ne. mnsp2) .and. (nsp1p .eq. nsp2p .and. mnsp1p .eq. mnsp2p)) then
     ! <a1a2 | ..| bb>
     twoelme  =  sqrt(2d0) * twoelme 
  elseif((nsp1 .eq. nsp2 .and. mnsp1 .eq. mnsp2) .and. (nsp1p .eq. nsp2p .and. mnsp1p .eq. mnsp2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...  
  else
     !<a1 a2 | ..| b1 b2>
     call V12me(p1,p2,p2p,p1p,mnsp1,mnsp2,mnsp2p,mnsp1p,result)
     twoelme  = twoelme + (-1)**(is) * result
  endif
    
!!$ Nuclei repulsion term for diatomic molecules:
  tmp = 0d0
  if( Rd .ne. 0) then
     tmp = Z1*Z2/Rd  * resultb
  endif
  
  resultH = oneelme + tmp + twoelme
  
  
  return
end subroutine H12me
!$**************************************************************************************************
!
!$$  two-electron configurations are made from one-electron target states, 
!$$  these states are not orthogonal
!$$  the idea here is to use H2+ 1s orbital together with Laguerre expansion (large exp. fall-off) 
!$$  to capture electron-electron corellations for low lying states.

!$$   This subroutine is to calculate matrix elements of H for H_2 molecule for 4 one-electron target states
!$$   with given projection of orbital angular momentum and parity
!$$   antisymmetric configurations 

subroutine H12me_st_notortog(Rd,is,TargetStates1el,Nmax1el,e1me,ovlpst,bst,nst1,nst2,nst1p,nst2p,mnst1,mnst2,mnst1p,mnst2p,resultH,resultb)
  
  use sturmian_class
  use state_class
  use MPI_module
  use input_data

  implicit none

  real*8, intent(in):: Rd
  integer, intent(in):: is
  type(basis_state), intent(in):: TargetStates1el
  integer, intent(in):: Nmax1el
  real*8, dimension(Nmax1el,Nmax1el), intent(in):: e1me, ovlpst
  type(basis_sturmian_nr), intent(in):: bst
  integer, intent(in):: nst1,nst2,nst1p,nst2p
  integer, intent(in):: mnst1,mnst2,mnst1p,mnst2p
  real*8, intent(out):: resultH,resultb

  type(sturmian_nr), pointer:: p1, p2, p1p, p2p
  integer:: nsp1,nsp2,nsp1p,nsp2p, i1,i2,i1p,i2p
  real*8:: tmp, tmpCI1, tmpCI1p, tmpCI2, tmpCI2p, ttt, result, oneelME, twoelME
  integer:: i1max, i2max, i1pmax, i2pmax
  real*8:: a1,a2, a3,a4, a5, a6, Z1, Z2
  
  resultH = 0d0
  resultb = 0d0

  oneelme = 0d0
  twoelme = 0d0
  tmp = 0d0
  
  Z1 = data_in%Z1
  Z2 = data_in%Z2

  i1max = get_nam(TargetStates1el%b(nst1))
  i2max = get_nam(TargetStates1el%b(nst2))
  i1pmax = get_nam(TargetStates1el%b(nst1p))
  i2pmax = get_nam(TargetStates1el%b(nst2p))


!!$   overlap and one-electron operator matrix              
  if(mnst1 .eq. mnst1p .and. mnst2 .eq. mnst2p) then

     resultb = ovlpst(nst1,nst1p) * ovlpst(nst2,nst2p)

     oneelme =   e1me(nst1,nst1p) * ovlpst(nst2,nst2p)

     oneelme =  oneelme + e1me(nst2,nst2p) * ovlpst(nst1,nst1p)

  endif  

  if((nst1 .eq. nst2) .and. (nst1p .ne. nst2p)) then
     ! <aa | ..| b1 b2>
     resultb = sqrt(2d0) * resultb 
     oneelme = sqrt(2d0) * oneelme
  elseif((nst1 .ne. nst2) .and. (nst1p .eq. nst2p)) then
     ! <a1a2 | ..| bb>
     resultb = sqrt(2d0) * resultb 
     oneelme = sqrt(2d0) * oneelme
  elseif((nst1 .eq. nst2) .and. (nst1p .eq. nst2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...  
  else
     !<a1 a2 | ..| b1 b2>
     if(mnst1 .eq. mnst2p .and. mnst2 .eq. mnst1p) then

        resultb = resultb + (-1)**(is) * ovlpst(nst1,nst2p) * ovlpst(nst2,nst1p)
!!$        resultb = resultb + (-1)**(is) * ovlp_st(TargetStates1el%b(nst1),TargetStates1el%b(nst2p)) * ovlp_st(TargetStates1el%b(nst2),TargetStates1el%b(nst1p))
        
        oneelme =  oneelme + (-1)**(is) * e1me(nst1,nst2p) * ovlpst(nst2,nst1p) ! ovlp_st(TargetStates1el%b(nst2),TargetStates1el%b(nst1p)) *  H1el_st(TargetStates1el%b(nst1),TargetStates1el%b(nst2p))
        
        oneelme =  oneelme + (-1)**(is) * e1me(nst2,nst1p) *  ovlpst(nst1,nst2p) ! ovlp_st(TargetStates1el%b(nst1),TargetStates1el%b(nst2p)) *  H1el_st(TargetStates1el%b(nst2),TargetStates1el%b(nst1p))

     endif
  endif
  
  
!!$   two-electron operator  matrix
  ttt = 0d0
  do i1=1,i1max
     nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
     tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
     p1 => bst%b(nsp1)
     do i2=1,i2max
        nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
        tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
        p2 => bst%b(nsp2)
        
        do i1p=1,i1pmax
           nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
           tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)     
           p1p => bst%b(nsp1p)
           do i2p=1,i2pmax
              nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
              tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
              p2p => bst%b(nsp2p)
              
              tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p

              call V12me(p1,p2,p1p,p2p,mnst1,mnst2,mnst1p,mnst2p,result)              
              
              ttt = ttt +  tmp * result

           enddo
        enddo
     enddo
  enddo

  twoelME = ttt

  !if (myid==0)  print*, 'result', result
  if((nst1 .eq. nst2) .and. (nst1p .ne. nst2p)) then
     ! <aa | ..| b1 b2>
     twoelme  =  sqrt(2d0) * twoelme 
  elseif((nst1 .ne. nst2) .and. (nst1p .eq. nst2p)) then
     ! <a1a2 | ..| bb>
     twoelme  =  sqrt(2d0) * twoelme 
  elseif((nst1 .eq. nst2) .and. (nst1p .eq. nst2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...  
  else
     !<a1 a2 | ..| b1 b2>
     ttt = 0d0
     do i1=1,i1max
        nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
        tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
        p1 => bst%b(nsp1)
        do i2=1,i2max
           nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
           tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
           p2 => bst%b(nsp2)
           
           do i1p=1,i1pmax
              nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
              tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)
              p1p => bst%b(nsp1p)     
              do i2p=1,i2pmax
                 nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
                 tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
                 p2p => bst%b(nsp2p)
                 
                 tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p
                 
                 call V12me(p1,p2,p2p,p1p,mnst1,mnst2,mnst2p,mnst1p,result)              
                 
                 ttt = ttt + tmp * result
                 
              enddo
           enddo
        enddo
     enddo
     twoelme  = twoelme + (-1)**(is) * ttt



  endif
    
!!$ Nuclei repulsion term for diatomic molecules:
  tmp = 0d0
  if( Rd .ne. 0) then
     tmp = Z1*Z2/Rd  * resultb
  endif

  resultH = oneelme + tmp + twoelme
  
  return
end subroutine H12me_st_notortog


subroutine H12me_st_notortog_test(Znuc,Rd,is,TargetStates1el,Nmax1el,e1me,ovlpst,bst,nst1,nst2,nst1p,nst2p,mnst1,mnst2,mnst1p,mnst2p,resultH,resultb)

  use sturmian_class
  use state_class
  use MPI_module
  use input_data

  implicit none

  real*8, intent(in):: Znuc, Rd
  integer, intent(in):: is
  type(basis_state), intent(in):: TargetStates1el
  integer, intent(in):: Nmax1el
  real*8, dimension(Nmax1el,Nmax1el), intent(in):: e1me, ovlpst
  type(basis_sturmian_nr), intent(in):: bst
  integer, intent(in):: nst1,nst2,nst1p,nst2p
  integer, intent(in):: mnst1,mnst2,mnst1p,mnst2p
  real*8, intent(out):: resultH,resultb

  type(sturmian_nr), pointer:: p1, p2, p1p, p2p
  integer:: nsp1,nsp2,nsp1p,nsp2p, i1,i2,i1p,i2p
  real*8:: tmp, tmpCI1, tmpCI1p, tmpCI2, tmpCI2p, ttt, result, oneelME, twoelME
  integer:: i1max, i2max, i1pmax, i2pmax
  real*8:: a1,a2, a3,a4, a5, a6, Z1, Z2

  resultH = 0d0
  resultb = 0d0

  oneelme = 0d0
  twoelme = 0d0
  tmp = 0d0
  
  Z1 = data_in%Z1
  Z2 = data_in%Z2

  i1max = get_nam(TargetStates1el%b(nst1))
  i2max = get_nam(TargetStates1el%b(nst2))
  i1pmax = get_nam(TargetStates1el%b(nst1p))
  i2pmax = get_nam(TargetStates1el%b(nst2p))


!!$   overlap and one-electron operator matrix              
  if(mnst1 .eq. mnst1p .and. mnst2 .eq. mnst2p) then

     resultb = ovlpst(nst1,nst1p) * ovlpst(nst2,nst2p)

!     oneelme =   e1me(nst1,nst1p) * ovlpst(nst2,nst2p)

     oneelme =  oneelme + e1me(nst2,nst2p) * ovlpst(nst1,nst1p)

  endif

  if((nst1 .eq. nst2) .and. (nst1p .ne. nst2p)) then
     ! <aa | ..| b1 b2>
!     resultb = sqrt(2d0) * resultb
!     oneelme = sqrt(2d0) * oneelme
     resultb =  resultb
     oneelme =  oneelme

  elseif((nst1 .ne. nst2) .and. (nst1p .eq. nst2p)) then
     ! <a1a2 | ..| bb>
!     resultb = sqrt(2d0) * resultb
!     oneelme = sqrt(2d0) * oneelme
     resultb =  resultb
     oneelme =  oneelme

  elseif((nst1 .eq. nst2) .and. (nst1p .eq. nst2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...  
  else
     !<a1 a2 | ..| b1 b2>
     if(mnst1 .eq. mnst2p .and. mnst2 .eq. mnst1p) then

        resultb = resultb ! + (-1)**(is) * ovlpst(nst1,nst2p) * ovlpst(nst2,nst1p)
!!$        resultb = resultb + (-1)**(is) * ovlp_st(TargetStates1el%b(nst1),TargetStates1el%b(nst2p)) * ovlp_st(TargetStates1el%b(nst2),TargetStates1el%b(nst1p))

        oneelme =  oneelme !+ (-1)**(is) * e1me(nst1,nst2p) * ovlpst(nst2,nst1p) ! ovlp_st(TargetStates1el%b(nst2),TargetStates1el%b(nst1p)) *  H1el_st(TargetStates1el%b(nst1),TargetStates1el%b(nst2p))

        oneelme =  oneelme !+ (-1)**(is) * e1me(nst2,nst1p) *  ovlpst(nst1,nst2p) ! ovlp_st(TargetStates1el%b(nst1),TargetStates1el%b(nst2p)) *  H1el_st(TargetStates1el%b(nst2),TargetStates1el%b(nst1p))

     endif
  endif


!!$   two-electron operator  matrix
  ttt = 0d0
  do i1=1,i1max
     nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
     tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
     p1 => bst%b(nsp1)
     do i2=1,i2max
        nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
        tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
        p2 => bst%b(nsp2)

        do i1p=1,i1pmax
           nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
           tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)
           p1p => bst%b(nsp1p)
           do i2p=1,i2pmax
              nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
              tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
              p2p => bst%b(nsp2p)

              tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p

              call V12me(p1,p2,p1p,p2p,mnst1,mnst2,mnst1p,mnst2p,result)

              ttt = ttt +  tmp * result

           enddo
        enddo
     enddo
  enddo

  twoelME = ttt

  !if (myid==0)  print*, 'result', result
  if((nst1 .eq. nst2) .and. (nst1p .ne. nst2p)) then
     ! <aa | ..| b1 b2>
     twoelme  =  sqrt(2d0) * twoelme
  elseif((nst1 .ne. nst2) .and. (nst1p .eq. nst2p)) then
     ! <a1a2 | ..| bb>
     twoelme  =  sqrt(2d0) * twoelme
  elseif((nst1 .eq. nst2) .and. (nst1p .eq. nst2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...  
  else
     !<a1 a2 | ..| b1 b2>
     ttt = 0d0
     do i1=1,i1max
        nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
        tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
        p1 => bst%b(nsp1)
        do i2=1,i2max
           nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
           tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
           p2 => bst%b(nsp2)

           do i1p=1,i1pmax
              nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
              tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)
              p1p => bst%b(nsp1p)
              do i2p=1,i2pmax
                 nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
                 tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
                 p2p => bst%b(nsp2p)

                 tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p

                 call V12me(p1,p2,p2p,p1p,mnst1,mnst2,mnst2p,mnst1p,result)

                 ttt = ttt + tmp * result

              enddo
           enddo
        enddo
     enddo
     twoelme  = twoelme + (-1)**(is) * ttt



  endif

!!$ Nuclei repulsion term for diatomic molecules:
  tmp = 0d0
  if( Rd .ne. 0) then
     tmp = Z1*Z2/Rd  * resultb
  endif

!  resultH = oneelme + tmp + twoelme
  resultH = oneelme + tmp ! + twoelme


  return
end subroutine H12me_st_notortog_test

!$**************************************************************************************************
!$$   two-electron configurations are made from one-electron target states, 
!$$   This subroutine is to calculate matrix elements of H for H_2 molecule for 4 one-electron target states
!$$   with given projection of orbital angular momentum and parity
!$$   antisymmetric configurations 

subroutine H12me_st(Znuc,Rd,is,TargetStates1el,bst,nst1,nst2,nst1p,nst2p,mnst1,mnst2,mnst1p,mnst2p,resultH,resultb)
  
  use sturmian_class
  use state_class
  use MPI_module
  use input_data

  implicit none

  real*8, intent(in):: Znuc, Rd
  integer, intent(in):: is
  type(basis_state), intent(in):: TargetStates1el
  type(basis_sturmian_nr), intent(in):: bst
  integer, intent(in):: nst1,nst2,nst1p,nst2p
  integer, intent(in):: mnst1,mnst2,mnst1p,mnst2p
  real*8, intent(out):: resultH,resultb

  type(sturmian_nr), pointer:: p1, p2, p1p, p2p
  integer:: nsp1,nsp2,nsp1p,nsp2p, i1,i2,i1p,i2p
  real*8:: tmp, tmpCI1, tmpCI1p, tmpCI2, tmpCI2p, ttt, result, oneelME, twoelME
  integer:: i1max, i2max, i1pmax, i2pmax
  real*8 :: Z1, Z2
  
  resultH = 0d0
  resultb = 0d0
  
  Z1 = data_in%Z1
  Z2 = data_in%Z2

  i1max = get_nam(TargetStates1el%b(nst1))
  i2max = get_nam(TargetStates1el%b(nst2))
  i1pmax = get_nam(TargetStates1el%b(nst1p))
  i2pmax = get_nam(TargetStates1el%b(nst2p))
!  if (myid==0) print*, 'i1max=', i1max
!  if (myid==0) print*, 'CI=', get_CI(TargetStates1el%b(nst1),1)
!  if (myid==0) print*, 'nst1 =',  get_na(TargetStates1el%b(nst1),1,1)


!!$  deal with one-electron ME: diagonal for one-electron target states

  oneelME = 0d0
  
  if(nst1 .eq. nst1p .and. nst2 .eq. nst2p) then

     resultb = 1d0

     oneelME = get_energy_st(TargetStates1el%b(nst1)) + get_energy_st(TargetStates1el%b(nst2))
     
  endif

!!$   two-electron operator  matrix
  twoelME = 0d0

  ttt = 0d0
  do i1=1,i1max
     nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
     tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
     p1 => bst%b(nsp1)
     do i2=1,i2max
        nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
        tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
        p2 => bst%b(nsp2)
        
        do i1p=1,i1pmax
           nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
           tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)     
           p1p => bst%b(nsp1p)
           do i2p=1,i2pmax
              nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
              tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
              p2p => bst%b(nsp2p)
              
              tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p

              call V12me(p1,p2,p1p,p2p,mnst1,mnst2,mnst1p,mnst2p,result)              
              
              ttt = ttt +  tmp * result

           enddo
        enddo
     enddo
  enddo

  twoelME = ttt

  !                 if (myid==0) print*, 'result', result
  if((nst1 .eq. nst2) .and. (nst1p .ne. nst2p)) then
     ! <aa | ..| b1 b2>
     twoelme  =  sqrt(2d0) * twoelme 
  elseif((nst1 .ne. nst2) .and. (nst1p .eq. nst2p)) then
     ! <a1a2 | ..| bb>
     twoelme  =  sqrt(2d0) * twoelme 
  elseif((nst1 .eq. nst2) .and. (nst1p .eq. nst2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...  
  else
     !<a1 a2 | ..| b1 b2>
     ttt = 0d0
     do i1=1,i1max
        nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
        tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
        p1 => bst%b(nsp1)
        do i2=1,i2max
           nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
           tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
           p2 => bst%b(nsp2)
           
           do i1p=1,i1pmax
              nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
              tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)
              p1p => bst%b(nsp1p)     
              do i2p=1,i2pmax
                 nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
                 tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
                 p2p => bst%b(nsp2p)
                 
                 tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p
                 
                 call V12me(p1,p2,p2p,p1p,mnst1,mnst2,mnst2p,mnst1p,result)              
                 
                 ttt = ttt + tmp * result
                 
              enddo
           enddo
        enddo
     enddo
     twoelme  = twoelme + (-1)**(is) * ttt
  endif
  
!!$ Nuclei repulsion term for diatomic molecules:
!!$ Note: H_2 Hamiltonian H = H2+ + H2+ + V12.
!!$ H2+ code has Znuc*Znuc/R. Need to compensate for the 2*Znuc*Znuc/R.
!!$ Implemented a -Znuc*Znuc/R factor in H2 structure code.
  tmp = 0d0
  if( Rd .ne. 0) then
!     tmp = Znuc*Znuc/Rd  * resultb
     tmp = - Z1*Z2/Rd  * resultb
  endif

  resultH = oneelme + tmp +  twoelme

  return
end subroutine H12me_st


!$**************************************************************************************************
!$$ Subroutine: H12me_st_group
!$$ Purpose:  two-electron configurations are made from one-electron target states, 
!$$   This subroutine is to calculate matrix elements of H for H_3^+ molecule for 4 one-electron target states
!$$   with given symmetry group representation as a label
!$$   antisymmetric configurations 

subroutine H12me_st_group(is,indata,TargetStates1el,bst,nst1,nst2,nst1p,nst2p,repnsp1,repnsp2,repnsp1p,repnsp2p,resultH,resultb)
  
  use basismodule 
  use sturmian_class
  use state_class
  use MPI_module
  use input_data

  implicit none

  integer, intent(in):: is
  type(smallinput):: indata
  type(basis_state), intent(in):: TargetStates1el
  type(basis_sturmian_nr), intent(in):: bst
  integer, intent(in):: nst1,nst2,nst1p,nst2p
  integer, intent(in):: repnsp1,repnsp2,repnsp1p,repnsp2p   !Group irreps of the states
  real*8, intent(out):: resultH,resultb

  type(sturmian_nr), pointer:: p1, p2, p1p, p2p
  integer:: nsp1,nsp2,nsp1p,nsp2p, i1,i2,i1p,i2p
  integer:: mnst1,mnst2,mnst1p,mnst2p
  real*8:: tmp, tmpCI1, tmpCI1p, tmpCI2, tmpCI2p, ttt, result, oneelME, twoelME
  integer:: i1max, i2max, i1pmax, i2pmax
  integer:: ii, jj
  real*8:: Rij, cosij
  
  resultH = 0d0
  resultb = 0d0

  i1max = get_nam(TargetStates1el%b(nst1))
  i2max = get_nam(TargetStates1el%b(nst2))
  i1pmax = get_nam(TargetStates1el%b(nst1p))
  i2pmax = get_nam(TargetStates1el%b(nst2p))
!  if (myid==0) print*, 'i1max=', i1max
!  if (myid==0) print*, 'CI=', get_CI(TargetStates1el%b(nst1),1)
!  if (myid==0) print*, 'nst1 =',  get_na(TargetStates1el%b(nst1),1,1)


!!$  deal with one-electron ME: diagonal for one-electron target states

  oneelME = 0d0
  
  if(nst1 .eq. nst1p .and. nst2 .eq. nst2p) then

     resultb = 1d0

     oneelME = get_energy_st(TargetStates1el%b(nst1)) + get_energy_st(TargetStates1el%b(nst2))
  endif

!!$   two-electron operator  matrix
  twoelME = 0d0

  ttt = 0d0
  do i1=1,i1max
     nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
     tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
     p1 => bst%b(nsp1)

     do i2=1,i2max
        nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
        tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
        p2 => bst%b(nsp2)
        
        do i1p=1,i1pmax
           nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
           tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)     
           p1p => bst%b(nsp1p)
           do i2p=1,i2pmax
              nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
              tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
              p2p => bst%b(nsp2p)

              tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p

	      mnst1 = get_ang_mom_proj_nr(p1)
	      mnst2 = get_ang_mom_proj_nr(p2)
	      mnst2p = get_ang_mom_proj_nr(p2p)
	      mnst1p = get_ang_mom_proj_nr(p1p)

              call V12me_group(indata,p1,p2,p1p,p2p,mnst1,mnst2,mnst1p,mnst2p,result)              
              
              ttt = ttt +  tmp * result

           enddo
        enddo
     enddo
  enddo

  twoelME = ttt

  ! if (myid==0) print*, 'result', result
  if((nst1 .eq. nst2) .and. (nst1p .ne. nst2p)) then
     ! <aa | ..| b1 b2>
     twoelme  =  sqrt(2d0) * twoelme 
  elseif((nst1 .ne. nst2) .and. (nst1p .eq. nst2p)) then
     ! <a1a2 | ..| bb>
     twoelme  =  sqrt(2d0) * twoelme 
  elseif((nst1 .eq. nst2) .and. (nst1p .eq. nst2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...  
  else
     !<a1 a2 | ..| b1 b2>
     ttt = 0d0
     do i1=1,i1max
        nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
        tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
        p1 => bst%b(nsp1)
        do i2=1,i2max
           nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
           tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
           p2 => bst%b(nsp2)
           
           do i1p=1,i1pmax
              nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
              tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)
              p1p => bst%b(nsp1p)     
              do i2p=1,i2pmax
                 nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
                 tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
                 p2p => bst%b(nsp2p)
                 
                 tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p

                 mnst1 = get_ang_mom_proj_nr(p1)
                 mnst2 = get_ang_mom_proj_nr(p2)
                 mnst2p = get_ang_mom_proj_nr(p2p)
                 mnst1p = get_ang_mom_proj_nr(p1p)

                 
                 !Calculates V12 matrix elements between functions with well defined lm
                 call V12me_group(indata,p1,p2,p2p,p1p,mnst1,mnst2,mnst2p,mnst1p,result)              
                 
                 ttt = ttt + tmp * result
                 
              enddo
           enddo
        enddo
     enddo
     twoelme  = twoelme + (-1)**(is) * ttt
  endif
  
!!$ Nuclei repulsion term for triatomic molecules (really more general):
!!$ Note: H_3+ Hamiltonian H = H3++ + H3++ + V12 + Hnuclei-nuclei
!!$ Where H3++_no is the H3++ hamiltonian without nuclear-nuclear interaction
   tmp = 0d0
   do ii=1, 3
      do jj = ii+1, 3
         !Use law of cosines to compute distance between nuclei
	       cosij = cos(indata%theta(ii))*cos(indata%theta(jj)) + &
	         sin(indata%theta(ii))*sin(indata%theta(jj))*cos(indata%phi(ii)-indata%phi(jj))
         Rij = sqrt(indata%R(ii)**2 + indata%R(jj)**2 - &
	       2*indata%R(ii)*indata%R(jj)*cosij)
	       !Account for degenerate case where nuclei coincide
         if (Rij .gt. 0.0_dpf) then
            tmp = tmp + (dble(indata%charge(ii)*indata%charge(jj))/Rij) * resultb
         end if
      end do
   end do

!Old comment kept for posterity.
!!$ H2+ code has Znuc*Znuc/R. Need to compensate for the 2*Znuc*Znuc/R.
!!$ Implemented a -Znuc*Znuc/R factor in H2 structure code.
!  if( Rd .ne. 0) then
!!     tmp = Znuc*Znuc/Rd  * resultb
!     tmp = - Z1*Z2/Rd  * resultb
!  endif

  !H3++ code includes nuclear interaction, need to compensate 
  !for 2*sum_ij z_i*z_j/R_ij introduced by using two sets of 1e states, subtract tmp
  resultH = oneelme - tmp + twoelme

  return
end subroutine H12me_st_group



!$**************************************************************************************************
!$$ Subroutine: H12me_st_group_notortog
!$$ Purpose:  two-electron configurations are made from one-electron target states, 
!$$   This subroutine is to calculate matrix elements of H for H_3^+ molecule for 4 one-electron target states
!$$   with given symmetry group representation as a label
!$$   antisymmetric configurations 
!$$
!$$   Uses a mixture of laguerre function with large alpha and 1e molecular orbitals to correctly desribe
!$$   excited states and improve convergence.
subroutine H12me_st_group_notortog(is,indata,TargetStates1el,Nmax1el,e1me,ovlpst,bst,nst1,nst2,nst1p,nst2p,repnsp1,repnsp2,repnsp1p,repnsp2p,resultH,resultb)
  
  use basismodule 
  use sturmian_class
  use state_class
  use MPI_module
  use input_data

  implicit none

  integer, intent(in):: is
  type(smallinput):: indata
  type(basis_state), intent(in):: TargetStates1el
	integer:: Nmax1el
	real*8, dimension(Nmax1el,Nmax1el):: e1me, ovlpst
  type(basis_sturmian_nr), intent(in):: bst
  integer, intent(in):: nst1,nst2,nst1p,nst2p
  integer, intent(in):: repnsp1,repnsp2,repnsp1p,repnsp2p   !Group irreps of the states
  real*8, intent(out):: resultH,resultb

  type(sturmian_nr), pointer:: p1, p2, p1p, p2p
  integer:: nsp1,nsp2,nsp1p,nsp2p, i1,i2,i1p,i2p
  integer:: mnst1,mnst2,mnst1p,mnst2p
  real*8:: tmp, tmpCI1, tmpCI1p, tmpCI2, tmpCI2p, ttt, result, oneelME, twoelME
  integer:: i1max, i2max, i1pmax, i2pmax
  integer:: ii, jj
  real*8:: Rij, cosij
  
  resultH = 0d0
  resultb = 0d0

  i1max = get_nam(TargetStates1el%b(nst1))
  i2max = get_nam(TargetStates1el%b(nst2))
  i1pmax = get_nam(TargetStates1el%b(nst1p))
  i2pmax = get_nam(TargetStates1el%b(nst2p))
!  if (myid==0) print*, 'i1max=', i1max
!  if (myid==0) print*, 'CI=', get_CI(TargetStates1el%b(nst1),1)
!  if (myid==0) print*, 'nst1 =',  get_na(TargetStates1el%b(nst1),1,1)

!!$  deal with one-electron ME: sum of two terms H_1*S_2 + H_2*S_1

  oneelME = 0d0
  
	!Basis functions from different irrps will be orthogonal
  if(repnsp1 .eq. repnsp1p .and. repnsp2 .eq. repnsp2p) then

     !resultb = 1d0
     !oneelME = get_energy_st(TargetStates1el%b(nst1)) + get_energy_st(TargetStates1el%b(nst2))

     resultb = ovlpst(nst1,nst1p) * ovlpst(nst2,nst2p)

     oneelme =   e1me(nst1,nst1p) * ovlpst(nst2,nst2p)

     oneelme =  oneelme + e1me(nst2,nst2p) * ovlpst(nst1,nst1p)

		 !if (nst2 .eq. 4) then
		 ! 	print*, "NST, NST2P", nst2, nst2p
		 ! 	print*, ovlpst(nst1,nst1p), ovlpst(nst2,nst2p), resultb
	   !end if

  endif

  if((nst1 .eq. nst2) .and. (nst1p .ne. nst2p)) then
     ! <aa | ..| b1 b2>
     resultb = sqrt(2d0) * resultb 
     oneelme = sqrt(2d0) * oneelme
  elseif((nst1 .ne. nst2) .and. (nst1p .eq. nst2p)) then
     ! <a1a2 | ..| bb>
     resultb = sqrt(2d0) * resultb 
     oneelme = sqrt(2d0) * oneelme
  elseif((nst1 .eq. nst2) .and. (nst1p .eq. nst2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...  
  else
     !<a1 a2 | ..| b1 b2>
     if(repnsp1 .eq. repnsp2p .and. repnsp2 .eq. repnsp1p) then

        resultb = resultb + (-1)**(is) * ovlpst(nst1,nst2p) * ovlpst(nst2,nst1p)
!!$     resultb = resultb + (-1)**(is) * ovlp_st(TargetStates1el%b(nst1),TargetStates1el%b(nst2p)) * ovlp_st(TargetStates1el%b(nst2),TargetStates1el%b(nst1p))
        
        oneelme =  oneelme + (-1)**(is) * e1me(nst1,nst2p) * ovlpst(nst2,nst1p) ! ovlp_st(TargetStates1el%b(nst2),TargetStates1el%b(nst1p)) *  H1el_st(TargetStates1el%b(nst1),TargetStates1el%b(nst2p))
        
        oneelme =  oneelme + (-1)**(is) * e1me(nst2,nst1p) *  ovlpst(nst1,nst2p) ! ovlp_st(TargetStates1el%b(nst1),TargetStates1el%b(nst2p)) *  H1el_st(TargetStates1el%b(nst2),TargetStates1el%b(nst1p))


		 !if (nst2 .eq. 4) then
		 ! 	print*, "PART 2 NST, NST2P", nst2, nst2p
		 ! 	print*, ovlpst(nst1,nst1p), ovlpst(nst2,nst2p), resultb
	   !end if



     endif
  endif

!!$   two-electron operator  matrix
  twoelME = 0d0

  ttt = 0d0
  do i1=1,i1max
     nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
     tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
     p1 => bst%b(nsp1)

     do i2=1,i2max
        nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
        tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
        p2 => bst%b(nsp2)
        
        do i1p=1,i1pmax
           nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
           tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)     
           p1p => bst%b(nsp1p)
           do i2p=1,i2pmax
              nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
              tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
              p2p => bst%b(nsp2p)

              tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p

	            mnst1 = get_ang_mom_proj_nr(p1)
	            mnst2 = get_ang_mom_proj_nr(p2)
	            mnst2p = get_ang_mom_proj_nr(p2p)
	            mnst1p = get_ang_mom_proj_nr(p1p)

              call V12me_group(indata,p1,p2,p1p,p2p,mnst1,mnst2,mnst1p,mnst2p,result)              
              
              ttt = ttt +  tmp * result

           enddo
        enddo
     enddo
  enddo

  twoelME = ttt

  ! if (myid==0) print*, 'result', result
  if((nst1 .eq. nst2) .and. (nst1p .ne. nst2p)) then
     ! <aa | ..| b1 b2>
     twoelme  =  sqrt(2d0) * twoelme 
  elseif((nst1 .ne. nst2) .and. (nst1p .eq. nst2p)) then
     ! <a1a2 | ..| bb>
     twoelme  =  sqrt(2d0) * twoelme 
  elseif((nst1 .eq. nst2) .and. (nst1p .eq. nst2p)) then
     ! orbitals are the same  <aa|...|bb>
     ! all is already done...  
  else
     !<a1 a2 | ..| b1 b2>
     ttt = 0d0
     do i1=1,i1max
        nsp1 = get_na(TargetStates1el%b(nst1),i1,1)
        tmpCI1 = get_CI(TargetStates1el%b(nst1),i1)
        p1 => bst%b(nsp1)
        do i2=1,i2max
           nsp2 = get_na(TargetStates1el%b(nst2),i2,1)
           tmpCI2 = get_CI(TargetStates1el%b(nst2),i2)
           p2 => bst%b(nsp2)
           
           do i1p=1,i1pmax
              nsp1p = get_na(TargetStates1el%b(nst1p),i1p,1)
              tmpCI1p = get_CI(TargetStates1el%b(nst1p),i1p)
              p1p => bst%b(nsp1p)     
              do i2p=1,i2pmax
                 nsp2p = get_na(TargetStates1el%b(nst2p),i2p,1)
                 tmpCI2p = get_CI(TargetStates1el%b(nst2p),i2p)
                 p2p => bst%b(nsp2p)
                 
                 tmp = tmpCI1 * tmpCI1p * tmpCI2 * tmpCI2p

                 mnst1 = get_ang_mom_proj_nr(p1)
                 mnst2 = get_ang_mom_proj_nr(p2)
                 mnst2p = get_ang_mom_proj_nr(p2p)
                 mnst1p = get_ang_mom_proj_nr(p1p)

                 
                 !Calculates V12 matrix elements between functions with well defined lm
                 call V12me_group(indata,p1,p2,p2p,p1p,mnst1,mnst2,mnst2p,mnst1p,result)              
                 
                 ttt = ttt + tmp * result
                 
              enddo
           enddo
        enddo
     enddo
     twoelme  = twoelme + (-1)**(is) * ttt
  endif
  
!!$ Nuclei repulsion term for triatomic molecules (really more general):
!!$ Note: H_3+ Hamiltonian H = H3++ + H3++ + V12 + Hnuclei-nuclei
!!$ Where H3++_no is the H3++ hamiltonian without nuclear-nuclear interaction
   tmp = 0d0
   do ii=1, 3
      do jj = ii+1, 3
         !Use law of cosines to compute distance between nuclei
	       cosij = cos(indata%theta(ii))*cos(indata%theta(jj)) + &
	         sin(indata%theta(ii))*sin(indata%theta(jj))*cos(indata%phi(ii)-indata%phi(jj))
         Rij = sqrt(indata%R(ii)**2 + indata%R(jj)**2 - &
	       2*indata%R(ii)*indata%R(jj)*cosij)
	       !Account for degenerate case where nuclei coincide
         if (Rij .gt. 0.0_dpf) then
            tmp = tmp + (dble(indata%charge(ii)*indata%charge(jj))/Rij) * resultb
         end if
      end do
   end do

!Old comment kept for posterity.
!!$ H2+ code has Znuc*Znuc/R. Need to compensate for the 2*Znuc*Znuc/R.
!!$ Implemented a -Znuc*Znuc/R factor in H2 structure code.
!  if( Rd .ne. 0) then
!!     tmp = Znuc*Znuc/Rd  * resultb
!     tmp = - Z1*Z2/Rd  * resultb
!  endif

  !H3++ code includes nuclear interaction, need to compensate 
  !for 2*sum_ij z_i*z_j/R_ij introduced by using two sets of 1e states, subtract tmp
  resultH = oneelme + twoelme + tmp

  return
end subroutine H12me_st_group_notortog









!$**************************************************************************************************
!
! This is matrix element of V(1,2) electron-electron potential for 4 functions: 
! <n1 n2 | V(1,2) | n1p n2p >, where V(1,2) = sum_{lam} V_{lam}(1,2)
! Note: potential V(12) can not change the total ang.mom. projection M of a configuration
! This subroutine uses either real or complex spherical harmonics in the
! expansion of V(1,2), with corresonding gaunt coefficients for each
! case.
subroutine V12me_group(indata,pn1,pn2,pn1p,pn2p,m1,m2,m1p,m2p,result)

  use basismodule
  use sturmian_class
  use grid_radial 
  use MPI_module

  implicit none
  
  type(sturmian_nr), intent(in):: pn1,pn2,pn1p,pn2p 
  integer, intent(in):: m1, m1p, m2, m2p
  real*8, intent(out):: result
  
  real*8:: Yint
  

  integer:: l1, l2, l1p, l2p
  integer:: maxfm, minfm, i1, i2, minfun, maxfun, minfun1, maxfun1, lam, lammin, lammax
  real*8:: rlam, rq, reslam
  real*8:: tmp, tmp1, tmp2, sum1, sum2
  integer:: minf1,maxf1, minf1p,maxf1p,minf2,maxf2, minf2p,maxf2p
  real*8, pointer, dimension(:):: f1, f1p, f2, f2p 
  real*8, dimension(grid%nr):: temp, fun, fun1, fun11
  real*8:: rl1,rl2,rl1p,rl2p,rm1, rm2,rm1p,rm2p
  real*8:: factor, qmin, qmax, q, pi
  type(smallinput):: indata

  pi = 4.0d0*atan(1.0d0)

  result = 0d0

  l1 = get_ang_mom(pn1) 
  l2 = get_ang_mom(pn2)   
  l1p = get_ang_mom(pn1p)
  l2p = get_ang_mom(pn2p)
  rl1 = l1
  rl2 = l2
  rl1p = l1p
  rl2p = l2p

  rm1 = m1
  rm2 = m2
  rm1p = m1p
  rm2p = m2p

  !print*, "L2p, m2p: ", l2p, m2p, "L2, m2: ", l2, m2
  
  f1 => fpointer(pn1)
  f1p => fpointer(pn1p)
  
  f2 => fpointer(pn2)
  f2p => fpointer(pn2p)
  
  minf1 = get_minf(pn1)
  maxf1 = get_maxf(pn1)
  minf1p = get_minf(pn1p)
  maxf1p = get_maxf(pn1p)
  
  minf2 = get_minf(pn2)
  maxf2 = get_maxf(pn2)
  minf2p = get_minf(pn2p)
  maxf2p = get_maxf(pn2p)

  maxfm = min(maxf1,maxf1p)
  minfm = max(minf1,minf1p)

  fun1(minfm:maxfm) = f1(minfm:maxfm)*f1p(minfm:maxfm)*grid%weight(minfm:maxfm)
  
  minfun = max(minf2,minf2p)
  maxfun = min(maxf2,maxf2p)
  
  !Liam: removed weights to use accurate form 
  !fun(minfun:maxfun) = f2(minfun:maxfun)*f2p(minfun:maxfun) * grid%weight(minfun:maxfun)
  fun(minfun:maxfun) = f2(minfun:maxfun)*f2p(minfun:maxfun)

  lammin=max(abs(l1-l1p),abs(l2-l2p))
  lammax=min(l1+l1p,l2+l2p)

  if (indata%harmop .eq. 1) then
     lammin = 0
  end if

  do lam=lammin,lammax
     if (indata%harmop .eq. 0) then
	!Complex harmonics, only one q gives non-zero gaunt coeffs, no loop
	qmin = lam
	qmax = lam
     else if (indata%harmop .eq. 1) then
        !Real harmonics, multiple q give non-zero gaunt coeffs, loop
	qmin = -lam
	qmax = lam
     else
        print*, "ERROR: neither real nor complex harmonics selected, stopping. V12me_group"
        error stop
     end if

	
     do q = qmin, qmax
        rlam = lam

        if (indata%harmop .eq. 0) then
	   !Regular gaunt coeffs vanish for q != m1p - m1
           rq = m1p - m1
           tmp1 = (-1)**(nint(rq))*Yint(rl1,rm1,rlam,-rq,rl1p,rm1p)
           tmp2 = Yint(rl2,rm2,rlam,rq,rl2p,rm2p)
        else if (indata%harmop .eq. 1) then
	   !Real gaunt coeffs accept four different values of q
           rq = q
           tmp1 = Xint(rl1,rm1,rlam,rq,rl1p,rm1p)
           tmp2 = Xint(rl2,rm2,rlam,rq,rl2p,rm2p)
           factor = 4*pi/dble(2*lam+1) !4pi/(2l+1) not included in Xint
           tmp1 = tmp1*factor
	end if

        tmp = tmp1 * tmp2

!        if (myid==0) print*, 'lam, q, tmp', lam, rq,tmp1,tmp2

        if( tmp .eq. 0d0) then
           cycle
        endif 

        !call form(lam,fun,minfun,maxfun,maxfm,temp,i1,i2)
        call form_accurate(lam,fun,minfun,maxfun,maxfm,temp,i1,i2)
        
        minfun1 = max(i1,minfm)
        maxfun1 = min(i2,maxfm)
        
        fun11(minfun1:maxfun1) = fun1(minfun1:maxfun1) * temp(minfun1:maxfun1)
        
        reslam = SUM(fun11(minfun1:maxfun1))
    
        result = result + reslam * tmp
     end do
     
  enddo
  

end subroutine V12me_group




!$**************************************************************************************************
!
! This is matrix element of V(1,2) electron-electron potential for 4 functions: 
! <n1 n2 | V(1,2) | n1p n2p >, where V(1,2) = sum_{lam} V_{lam}(1,2)
! Note: potential V(12) can not change the total ang.mom. projection M of a configuration
subroutine V12me(pn1,pn2,pn1p,pn2p,m1,m2,m1p,m2p,result)

  use sturmian_class
  use grid_radial 
  use MPI_module

  implicit none
  
  type(sturmian_nr), intent(in):: pn1,pn2,pn1p,pn2p 
  integer, intent(in):: m1, m1p, m2, m2p
  real*8, intent(out):: result
  
  real*8:: Yint
  

  integer:: l1, l2, l1p, l2p
  integer:: maxfm, minfm, i1, i2, minfun, maxfun, minfun1, maxfun1, lam, lammin, lammax
  real*8:: rlam, rq, reslam
  real*8:: tmp, tmp1, tmp2, sum1, sum2
  integer:: minf1,maxf1, minf1p,maxf1p,minf2,maxf2, minf2p,maxf2p
  real*8, pointer, dimension(:):: f1, f1p, f2, f2p 
  real*8, dimension(grid%nr):: temp, fun, fun1, fun11
  real*8:: rl1,rl2,rl1p,rl2p,rm1, rm2,rm1p,rm2p


  result = 0d0

  l1 = get_ang_mom(pn1) 
  l2 = get_ang_mom(pn2)   
  l1p = get_ang_mom(pn1p)
  l2p = get_ang_mom(pn2p)
  rl1 = l1
  rl2 = l2
  rl1p = l1p
  rl2p = l2p

  rm1 = m1
  rm2 = m2
  rm1p = m1p
  rm2p = m2p
  
  f1 => fpointer(pn1)
  f1p => fpointer(pn1p)
  
  f2 => fpointer(pn2)
  f2p => fpointer(pn2p)
  
  minf1 = get_minf(pn1)
  maxf1 = get_maxf(pn1)
  minf1p = get_minf(pn1p)
  maxf1p = get_maxf(pn1p)
  
  minf2 = get_minf(pn2)
  maxf2 = get_maxf(pn2)
  minf2p = get_minf(pn2p)
  maxf2p = get_maxf(pn2p)


  maxfm = min(maxf1,maxf1p)
  minfm = max(minf1,minf1p)

  fun1(minfm:maxfm) = f1(minfm:maxfm)*f1p(minfm:maxfm)*grid%weight(minfm:maxfm)
  
  minfun = max(minf2,minf2p)
  maxfun = min(maxf2,maxf2p)
  
  !Liam: removed weights to use accurate form 
  !fun(minfun:maxfun) = f2(minfun:maxfun)*f2p(minfun:maxfun) * grid%weight(minfun:maxfun)
  fun(minfun:maxfun) = f2(minfun:maxfun)*f2p(minfun:maxfun)


  lammin=max(abs(l1-l1p),abs(l2-l2p))
  lammax=min(l1+l1p,l2+l2p)

  do lam=lammin,lammax

     rlam = lam
     rq = m1p - m1
     tmp1 = (-1)**(nint(rq))*Yint(rl1,rm1,rlam,-rq,rl1p,rm1p)
     tmp2 = Yint(rl2,rm2,rlam,rq,rl2p,rm2p)
     tmp = tmp1 * tmp2

!     if (myid==0) print*, 'lam, q, tmp', lam, rq,tmp1,tmp2

     if( tmp .eq. 0d0) then
        cycle
     endif
     

     !call form(lam,fun,minfun,maxfun,maxfm,temp,i1,i2)
     call form_accurate(lam,fun,minfun,maxfun,maxfm,temp,i1,i2)
     
     minfun1 = max(i1,minfm)
     maxfun1 = min(i2,maxfm)
     
     fun11(minfun1:maxfun1) = fun1(minfun1:maxfun1) * temp(minfun1:maxfun1)
     
     reslam = SUM(fun11(minfun1:maxfun1))
    
     result = result + reslam * tmp
     
  enddo
  

end subroutine V12me




!$*************************************************************************************************


subroutine convert_from_st_to_sp(bst,TargetStates1el,TargetStates2el,nspm,Nmax1el,Nmax)
  
  use sturmian_class
  use state_class
  use grid_radial
  use MPI_module

  implicit none
  
  type(basis_sturmian_nr), intent(in):: bst
  type(basis_state), intent(in):: TargetStates1el
  type(basis_state):: TargetStates2el 
  integer, intent(in):: nspm,Nmax,Nmax1el

  integer:: nst, namst, i, ist, n1st, n2st, nam1, nam2, i1, i2, n1, n2, ncm
  real*8:: CIst, CI1, CI2
!  real*8, dimension(nspm,nspm):: CIno
!  real*8, dimension(nspm*nspm):: CItmp
  real*8, dimension(:,:), allocatable:: CIno
  real*8, dimension(:), allocatable:: CItmp
  integer, dimension(:), allocatable:: no1, no2, mo1, mo2, phase
  integer, dimension(nspm):: nspar
  integer:: newnspm, j, jn1, jn2
  integer, dimension(:), allocatable:: arnsp

  if (myid==0) print*, 'start: convert_from_st_to_sp()'
  if (myid==0) print*, 'nspm,Nmax1el,Nmax:', nspm,Nmax1el,Nmax




!!!$ Mark Additions
!  if ( non_uniq_log ) then
!     allocate(is_core_MO(TargetStates1el%Nmax))
!!!$ is_core_MO = 0 1s1s,1s2s,1s2p,..,1snl <=> nl orbital type 0 
!!!$ is_core_MO = 1 2s2s <=> 2s orbital type 1
!!!$ is_core_MO = 2 1s1s,1s2s,1s2p,.., 1snl <=> 1s orbital type 2 
!!!$ is_core_MO - to determine which orthonormal molecular orbitals
!!!$ are core orbitals for the projection operator
!     is_core_MO(:) = 0
!
!     do n_core = 1, TargetStates2el%nicm
!!!$ If we used rearrange all core-orbital one-electron target states
!!!$ are the first nicm states in the array
!        if (inc == 1) then
!           nst_core = TargetStates2el%ncore(n_core)
!        else 
!           nst_core = TargetStates2el%ncore(n_core)
!        end if 
!
!        k_core = get_inum_st(TargetStates1el%b(nst_core))
!        l_majcore = get_l_majconf(TargetStates1el%b(nst_core))
!        if (l_ion_core >= 0) then
!           if ( k_core <= nk1(l_majcore) ) then
!                 is_core_MO(nst_core) = 1
!              print*,"st_ind, k, l", nst_core, k_core, l_majcore
!           else if ( k_core <= n_ion_core(l_majcore) ) then
!                 is_core_MO(nst_core) = 2
!              print*,"st_ind, k, l", nst_core, k_core, l_majcore
!           else
!              print*,"issues"
!              print*,"st_ind, k, l", nst_core, k_core, l_majcore
!              stop
!           end if
!        end if
!     end do ! n_core - loop over core obitals
!  end if
!!!$ end Additions




!!$ this loop can be paralalized
  do nst=1,Nmax
     
!     if (myid==0) print*, '>> nst =', nst
     
     namst = get_nam(TargetStates2el%b(nst))
!     if (myid==0) print*, 'namst=', namst
     

!!$ as nspm can be very large we need to avoid dealing with all one-electron orbitals
!!$ and use only those orbitals that are used for given state i 
!!$  First create two arrays nspar(1:nspm) and arnsp(1:newnspm) 
!!$
!!$ create array nspar(1:nspm) where we record 1 if the orbital is used in the description of the given two-electron state ist

     nspar(:) = 0
     
     do ist=1,namst
        
        n1st = get_na(TargetStates2el%b(nst),ist,1)
        nam1 = get_nam(TargetStates1el%b(n1st))
        do i1=1,nam1           
           n1 = get_na(TargetStates1el%b(n1st),i1)
           nspar(n1) = 1
        enddo
        
        n2st = get_na(TargetStates2el%b(nst),ist,2)
        nam2 = get_nam(TargetStates1el%b(n2st))
        do i2=1,nam2           
           n2 = get_na(TargetStates1el%b(n2st),i2)
           nspar(n2) = 1
        enddo
     enddo
!!$ check how many orbitals are used for decription of the state ist
     newnspm = 0
     do i=1,nspm
        if(nspar(i) .ne. 0) then
           newnspm = newnspm + 1        
        endif
     enddo
     allocate(arnsp(newnspm))
     arnsp(:) = 0
!     if (myid==0) print*, 'newnspm=', newnspm
     
     j = 0
     do i=1,nspm
        if(nspar(i) .ne. 0) then
           j = j + 1        
           arnsp(j) = i
           nspar(i) = j
        endif
     enddo
!!$ at this stage arrays nspar() and arnsp() allow to move forward and backwards between old enumeration (1 to nspm) 
!!$ and new enumeration (1 to newnspm)
!!$ Hopefully newnspm << nspm
     
     allocate(CIno(newnspm,newnspm))
     allocate(CItmp(newnspm*newnspm),no1(newnspm*newnspm),no2(newnspm*newnspm),mo1(newnspm*newnspm),mo2(newnspm*newnspm),phase(newnspm*newnspm))
     CIno(:,:) = 0d0
     
     do ist=1,namst
        
        n1st = get_na(TargetStates2el%b(nst),ist,1)
        n2st = get_na(TargetStates2el%b(nst),ist,2)
        
        CIst = get_CI(TargetStates2el%b(nst),ist)
        
        nam1 = get_nam(TargetStates1el%b(n1st))
        nam2 = get_nam(TargetStates1el%b(n2st))
!!$        if (myid==0) print*, 'n1st, n2st=', n1st, n2st
        do i1 = 1, nam1
           
           n1 = get_na(TargetStates1el%b(n1st),i1)
           CI1 = get_CI(TargetStates1el%b(n1st),i1)
!           mo1(n1) = get_angmom_proj(TargetStates1el%b(n1st))
           
           jn1 = nspar(n1)
           
           do i2 = 1, nam2
              
              n2 = get_na(TargetStates1el%b(n2st),i2)
              CI2 = get_CI(TargetStates1el%b(n2st),i2)
!              mo2(n2) = get_angmom_proj(TargetStates1el%b(n2st))
              jn2 = nspar(n2)

              CIno(jn1,jn2) = CIno(jn1,jn2) + CIst * CI1 * CI2
!!$              if (myid==0) print*, n1, n2, jn1, jn2,  CIno(jn1,jn2)
           end do ! i2
           
        end do ! i1

     end do ! ist


     i = 0
     do jn1 = 1, newnspm
        do jn2 = jn1, newnspm
           
           if(CIno(jn1,jn2) .eq. 0 ) cycle
        
           i = i + 1
           
           n1 = arnsp(jn1)
           no1(i) = n1
           n2 = arnsp(jn2)
           no2(i) = n2
           mo1(i) = get_ang_mom_proj(bst%b(n1))
           mo2(i) = get_ang_mom_proj(bst%b(n2))
           if(n1 .eq. n2)  then
              CItmp(i) = CIno(jn1,jn2)
              phase(i) = 1d0     
           else
              CItmp(i) = CIno(jn1,jn2) * sqrt(2d0)
              phase(i) = nint(CIno(jn1,jn2)/CIno(jn2,jn1))
           endif
!!$           if (myid==0) print*, i, n1, n2, CItmp(i), phase(i)
        end do ! jn2
     end do ! jn1
     
     ncm = i

!!!$ Mark: Additions
!!!$ Defining core orbitals as a function of the 
!!!$ two-electron state, two-electron config number and electrons 1 or 2 
!!!$ THIS WILL NOT WORK IF WE DON'T USE REARRANGE!!!
!     if ( non_uniq_log ) then
!        i = 0
!        do jn1 = 1, newnspm
!           do jn2 = jn1, newnspm
!              
!              if(CIno(jn1,jn2) .eq. 0 ) cycle
!              i = i + 1
!              
!              do ist=1,namst
!                 
!                 n1st = get_na(TargetStates2el%b(nst),ist,1)
!                 n2st = get_na(TargetStates2el%b(nst),ist,2)
!                 
!                 nam1 = get_nam(TargetStates1el%b(n1st))
!                 nam2 = get_nam(TargetStates1el%b(n2st))
!                 
!                 do i1=1,nam1
!                    n1 = get_na(TargetStates1el%b(n1st),i1)
!                    
!!!$ nspar contains index to j1 and j2 
!                    if ( jn1 /= nspar(n1) ) cycle
!                    
!                    do i2=1,nam2
!                       
!                       n2 = get_na(TargetStates1el%b(n2st),i2)
!                       
!                       if ( jn2 /= nspar(n2) ) cycle
!                       
!                       is_core_orb(nst,i,1) = is_core_MO(n1st) 
!                       is_core_orb(nst,i,2) = is_core_MO(n2st)  
!                       
!                    end do ! i2 - atomic orbital loop "outer" electron
!                 end do ! i1 - atomic orbital loop "inner" electron
!              end do ! ist - loop over two-electron configurations 
!           end do ! jn1
!        end do ! jn2
!     end if  ! non-uniqueness solved ?
!
!!!$ End of additions

                     
     call setupCI(TargetStates2el%b(nst),ncm,CItmp(1:ncm),no1(1:ncm),mo1(1:ncm),no2(1:ncm),mo2(1:ncm),phase(1:ncm))    
     deallocate(arnsp,CIno, CItmp,no1,no2,mo1,mo2,phase)


  enddo  ! nst - two-electron state
  
  if (myid==0) print*, 'finish convert_from_st_to_sp()'

end subroutine convert_from_st_to_sp


subroutine H2st_test()
!!$
  use input_data
  use grid_radial
  use sturmian_class
  use vnc_module
  use target_states
  use one_electron_func
  
  implicit none
  
  integer:: nstf, nsti        ! State index 
  
  type(sturmian_nr), pointer:: tn1, tnp1, tn2, tnp2     ! One-electron orbitals
  integer:: ind1, ind2, indp2, indp1                      ! index toone-electron orbital
  integer::  ne_con,  nep_con, Ncon, Npcon, n_orb, np_orb, n_orb_max, np_orb_max 
  integer:: nuse_norb, nuse_nporb 
  real*8, pointer, dimension(:):: fp1, f1, fp2, f2
  real*8, dimension(grid%nr)::  fun1, epot_i
  real*8:: H2_energy, HT_ham1el
  real*8:: temp_overlap2, overlap2, temp
  
  real*8:: Yint, ang1, ang2
  real*8:: CIp, CI  
  integer:: lambda, lambda_min, lambda_max, mu
  integer:: Lap1,Mp1,La1,M1,Lap2,Mp2,La2,M2    ! A.O. QN and Molecular State M 
  integer:: Mp12, M12, Spin, Spinp
  integer:: minf1, maxf1, minfp1, maxfp1, minf2, maxf2, minfp2, maxfp2
  integer:: ir1, ir2, jr1, jr2, or1, or2
  real*8, pointer, dimension(:):: weight, gridr

  if ( data_in%non_uniq ) then
     print*,"Subroutine H2st_test is not coded for the two-electron"
     print*," state being represented by 1e molecular states."       
     stop   
  end if        
  
  weight => grid%weight
  gridr => grid%gridr
  
 print*,"Checking target state reproduce H_T energy CIp*CI*<np1,np2|H_T|n1,n2>"
 print*,"nstf, nsti, CIp*CI*<np1,np2|H_T|n1,n2> using bst_nr%ham1el, two-electron energy"

  do nstf = 1, TargetStates2el%Nmax 
  do nsti = 1, nstf

  HT_ham1el = 0d0

  ! Final State number nstf
  Npcon =  TargetStates2el%b(nstf)%nam          ! Number of A.O. and Molecular ion configurations.        
  Mp12 = NINT(TargetStates2el%b(nstf)%M )        ! Molecular state Angular Projection
  Spinp = NINT(TargetStates2el%b(nstf)%spin)    ! 2e Molecular State Spin
  np_orb_max = TargetStates2el%b(nstf)%nusemax  ! Number of orbitals used to describe this state
  
  ! Initial State number nsti
  Ncon = TargetStates2el%b(nsti)%nam         
  M12 = NINT(TargetStates2el%b(nsti)%M)       
  Spin = NINT(TargetStates2el%b(nsti)%spin)
  n_orb_max = TargetStates2el%b(nsti)%nusemax         
  
  if (Spin == Spinp .AND. Mp12 == M12)  then      
     
     ! FINAL State COORDINATE 1
     do np_orb = 1, np_orb_max
        
        nuse_nporb = TargetStates2el%b(nstf)%nuse(np_orb)
        
        ! INITIAL State COORDINATE 0 
        do n_orb = 1, n_orb_max
           
           nuse_norb = TargetStates2el%b(nsti)%nuse(n_orb)   
           
           ! Below sums all the overlaps for COORDINATE 2  same
           ! configuratins(1s,2s,..) in COORDINATE 1
           overlap2 = 0d0
           
           ! Looping over FINAL Molecular State orbitals. COORDINATE 2
           do nep_con = 1, Npcon          
              
              indp1 = TargetStates2el%b(nstf)%na(nep_con)  ! Final state number np. nep A.O.
              
              if ( nuse_nporb /= indp1 ) cycle
              
              ! Quantum numbers and functions for BETA
              indp2 = TargetStates2el%b(nstf)%nb(nep_con)                
              tnp2 => bst_nr%b(indp2)                                         
              fp2 => fpointer(tnp2)                                
              Lap2 = get_ang_mom(tnp2)                             
              Mp2 = get_ang_mom_proj(tnp2)          
              CIp = get_CI(TargetStates2el%b(nstf),nep_con)  
              
              ! Looping over INITIAL Molecular State orbitals.  COORDINATE 2
              do ne_con = 1, Ncon       
                 
                 ind1 = TargetStates2el%b(nsti)%na(ne_con)  
                 
                 if ( nuse_norb /= ind1 ) cycle            
                 
                 ! Quantum numbers and functions for DELTA
                 ind2 = TargetStates2el%b(nsti)%nb(ne_con) 
                 tn2 => bst_nr%b(ind2)                                         
                 f2 => fpointer(tn2)                               
                 La2 = get_ang_mom(tn2)                            
                 M2 = get_ang_mom_proj(tn2)            
                 CI = get_CI(TargetStates2el%b(nsti),ne_con)  
                 
                 ! Selections Rules  
                 if ( Lap2 /= La2 .OR. Mp2 /= M2  ) cycle
                 
                 temp_overlap2 = 0d0
                 ! DO OVERLAP OF COORDINATE SPACE 2  
                 minf2 = get_minf(tn2)
                 maxf2 = get_maxf(tn2)
                 minfp2 = get_minf(tnp2)
                 maxfp2 = get_maxf(tnp2)      
                 ir1 = max(minf2,minfp2)
                 ir2 = min(maxf2,maxfp2)  
!                 temp_overlap2 = SUM(fp2(ir1:ir2) * f2(ir1:ir2) * weight(ir1:ir2))
                 temp_overlap2 = bst_nr%ortint(indp2,ind2)
                 
                 overlap2 = overlap2 + CIp * CI * temp_overlap2
                 
              end do  ! INITIAL STATE COORDINATE 2
           end do    ! FINAL STATE COORDINATE 2       
           
           if ( overlap2 /= 0d0 )  then
              
              ! COORDINATE 1 RADIAL INTEGRALS
              ! Quantum numbers and functions for ALPHA
              indp1 = nuse_nporb                                  ! Final state number np. nep A.O.
              tnp1 => bst_nr%b(indp1)                             !           
              fp1 => fpointer(tnp1)                               ! One electron functions
              Lap1 = get_ang_mom(tnp1)                            ! Gets Angular momentum A.O.
              Mp1 = get_ang_mom_proj(tnp1)                        ! Get angular projection of A.O. 
              
              ! Quantum numbers and functions for GAMMA
              ind1 = nuse_norb               
              tn1 => bst_nr%b(ind1)                                          
              f1 => fpointer(tn1)                               
              La1 = get_ang_mom(tn1)                            
              M1 = get_ang_mom_proj(tn1)  
              

              HT_ham1el = HT_ham1el + 2d0 * bst_nr%ham1el(indp1,ind1) * overlap2 
           end if ! overlap2
        end do ! INITIAL STATE COORDINATE 1    
     end do  ! FINAL STATE COORDINATE 1

!!$ Calculating -Cp*C*<np1|k1><kp0,np2|V_02|n0,n2>
!!$ Looping over FINAL Molecular State orbitals
  do nep_con =  1, Npcon 
        ! FINAL STATE COORDINATE 1 
        ! Quantum numbers and functions for ALPHA
        indp1 = TargetStates2el%b(nstf)%na(nep_con)
        tnp1 => bst_nr%b(indp1)       !           
        fp1 => fpointer(tnp1)         ! One electron functions
        Lap1 = get_ang_mom(tnp1)      ! Gets Angular momentum A.O.
        Mp1 = get_ang_mom_proj(tnp1)  ! Get angular projection of A.O. 
        minfp1 = get_minf(tnp1)
        maxfp1 = get_maxf(tnp1)
        
        ! FINAL STATE COORDINATE 2 
        ! Quantum numbers and functions for BETA
        indp2 = TargetStates2el%b(nstf)%nb(nep_con)
        tnp2 => bst_nr%b(indp2)
        fp2 => fpointer(tnp2)
        Lap2 = get_ang_mom(tnp2)
        Mp2 = get_ang_mom_proj(tnp2)
        CIp = get_CI(TargetStates2el%b(nstf),nep_con)
        minfp2 = get_minf(tnp2)
        maxfp2 = get_maxf(tnp2) 


     ! Looping over INITIAL Molecular State orbitals
     do ne_con = 1, Ncon       
        
        ! INITIAL STATE COORDINATE 0 
        ! Quantum numbers and functions for GAMMA
        ind1 = TargetStates2el%b(nsti)%na(ne_con) 
        tn1 => bst_nr%b(ind1)
        f1 => fpointer(tn1)
        La1 = get_ang_mom(tn1)
        M1 = get_ang_mom_proj(tn1)
        minf1 = get_minf(tn1)
        maxf1 = get_maxf(tn1)
        
        ! INITIAL STATE COORDINATE 2 
        ! Quantum numbers and functions for DELTA
        ind2 = TargetStates2el%b(nsti)%nb(ne_con) 
        tn2 => bst_nr%b(ind2)                                         
        f2 => fpointer(tn2)                               
        La2 = get_ang_mom(tn2)                            
        M2 = get_ang_mom_proj(tn2)            
        CI = get_CI(TargetStates2el%b(nsti),ne_con)  
        minf2 = get_minf(tn2)
        maxf2 = get_maxf(tn2)

        or1 = max(minfp1,minf1)     
        or2 = min(maxfp1,maxf1)
        
        !Liam: removed weights to use accurate form
        !fun1(or1:or2) = fp1(or1:or2) * f1(or1:or2) * weight(or1:or2)
        fun1(or1:or2) = fp1(or1:or2) * f1(or1:or2) !* weight(or1:or2)
        
        lambda_min = max(ABS(Lap1 - La1), ABS(Lap2 - La2))
        lambda_max = min(Lap1 + La1, Lap2 + La2)
        mu = Mp2 - M2
        do lambda = lambda_min, lambda_max
           ang1 = 0d0
           ang2 = 0d0
           ang1 = dble((-1)**(mu))*Yint(dble(Lap1),dble(Mp1),dble(lambda),dble(-mu),dble(La1),dble(M1))
           ang2 = Yint(dble(Lap2),dble(Mp2),dble(lambda),dble(mu),dble(La2),dble(M2))
           if ( ang1 * ang2 == 0d0 ) cycle

           ! Integrate over coordinate 1 
           !call form(lambda,fun1,or1,or2,grid%nr,epot_i,jr1,jr2)
           call form_accurate(lambda,fun1,or1,or2,grid%nr,epot_i,jr1,jr2)

              ir1 = max(minfp2,minf2)
              ir2 = min(maxfp2,maxf2)
              ir1 = max(jr1,ir1)
              ir2 = min(jr2,ir2) 
        
              temp = SUM(fp2(ir1:ir2) * epot_i(ir1:ir2) * f2(ir1:ir2) *weight(ir1:ir2))       
              temp = CI * CIp * temp * ang1 * ang2

              HT_ham1el = HT_ham1el + temp
        end do ! lambda


     end do ! INITIAL STATE 
  end do ! FINAL STATE
  end if ! Spin
  
  H2_energy = 0d0
  if ( nsti == nstf ) then
  if ( data_in%Rd /= 0d0 ) then
        HT_ham1el = HT_ham1el + 1d0 / data_in%Rd 
  end if
  H2_energy = get_energy(TargetStates2el%b(nsti))        
  end if      

 write(*,'(2I4,5X,2(F12.5,3X))') nstf,nsti,HT_ham1el,H2_energy
 end do ! nsti 
 end do ! nstf  
  
end subroutine H2st_test

subroutine Orthonormalise_one_electron_basis(TargetStates1el,bst)
!!$ MARK: This subroutine passes through the one-electron target state basis
!!$ 'TargetStates1el' and underlying one-electron Sturmian basis 'bst'. 
!!$ The one-electron target state basis 'TargetStates1el' is orthonomalised
!!$ using the Gram-Schmidt procedure.
!!$ Matrix elements e1me and ovlpst are also modified so that we can
!!$ perform two-electron diagonalisation
!!$ and construct two-electron target state with an orthonormal basis.
!!$ The goal is to end up with one-electron orbitals that are orthonormal,
!!$ an essential property for solving non-uniqueness.

  use input_data
  use grid_radial
  use sturmian_class
  use state_class
  use ovlpste1me
  use vmat_exch_module
  use MPI_module

  implicit none
! Will hold orthonormal one-electron (molecular orbital) basis
  type(basis_state), intent(inout):: TargetStates1el
! Will hold orthonormal one-electron basis underlying basis functions
  type(basis_sturmian_nr), intent(inout):: bst

  integer:: Nmax_1el_basis    ! One-electron basis size
  integer:: Nmax_SturmBasis     ! Underlying one-electron basis size
  integer:: basis_type          ! Spherical or Spheroidal
  integer:: nst_j, nst_k, nst_n ! One-electron basis indexes 
  integer:: Nj_1el, Nk_1el, Nn_1el ! Number of underlying atomic orbitals per one-electron (molecular orbital) basis
  integer:: nk, nn
  integer:: ind_k, ind_n  ! Index to underlying atomic orbitals basis 'bst'
  real*8:: CI_k, CI_n
  real*8:: sum1
  real*8, allocatable, dimension(:):: CI_new
  integer, allocatable, dimension(:):: na_new
!!$ For testing
  logical:: test_log = .FALSE.
  type(basis_state):: tmp1el_basis    ! One-electron basis
  type(basis_sturmian_nr):: temp_bst  ! Underlying one-electron basis
  real*8, allocatable, dimension(:,:):: tmp_ovlpst, tmp_e1me 

  if (myid == 0 ) write(*,'("")')
  if (myid == 0 ) write(*,'(" Gram-Schmidt orthogonalisation")')
  if (myid == 0 ) write(*,'("of one-electron (molecular orbital) target state basis")')
  if (myid == 0 ) write(*,'("")')

! Size of one-electron molecular orbital basis
  Nmax_1el_basis = basis_size_st(TargetStates1el)

!!$ Testing this subroutine requires the below before modifying structures 
  if (test_log ) then
!!$ Create new bases and copy one-electron basis and underlying functions
     call new_basis_st(tmp1el_basis,Nmax_1el_basis,.true.,data_in%calculation_type)
     do nst_j = 1, Nmax_1el_basis
        call copy_st(tmp1el_basis%b(nst_j),TargetStates1el%b(nst_j))
     end do
     Nmax_SturmBasis = basis_size(bst)
     call new_basis(temp_bst,Nmax_SturmBasis)
     call copy_basis_nr(temp_bst,bst)
     allocate(tmp_ovlpst(Nmax_1el_basis,Nmax_1el_basis), tmp_e1me(Nmax_1el_basis,Nmax_1el_basis))
     tmp_ovlpst(:,:) = ovlpst(:,:)
     tmp_e1me(:,:) = e1me(:,:)
  end if

!!$ u_j - old nonorthogonal molecular orbital set
!!$ Overlaps <u_j|u_k> in ovlpst(j,k)
!!$  <u_j|H_T|u_k> in e1me(j,k)

!!$ form overlap array <u_j|v_k>
!!$ u_j - old nonorthogonal molecular orbital set,
!!$ v_k - new orthogonal molecular orbital set but not yet normalised
!!$ Only elements with  j >= k  required.
!!$ This algorithm utilises the property <v_k|v_k> = <v_k|u_k>
!!$ |v_k> = |u_k> - sum^{k-1}_{n=1} <v_n|u_k> |v_n> / <v_n|v_n>
!!$ <u_j||v_k> = <u_j||u_k> - sum^{k-1}_{n=1} <v_n|u_k> <u_j||v_n> / <v_n|v_n>
!!$ hence <u_j||v_k> = <u_j||u_k> - sum^{k-1}_{n=1} <v_n|u_k> <u_j||v_n> / <v_n|u_n>
  do nst_j = 1, Nmax_1el_basis
     do nst_k = 1, nst_j
        sum1 = 0d0
        do nst_n = 1, nst_k - 1
           sum1 = sum1 + ovlpst(nst_k,nst_n) * ovlpst(nst_j,nst_n) / ovlpst(nst_n,nst_n)
        end do ! nst_n
     ovlpst(nst_j,nst_k) = ovlpst(nst_j,nst_k) - sum1
!     if (myid == 0 .AND. test_log) write(*,'("j,k =",2I3,", <u_j|v_k> =",F10.5)') nst_j, nst_k, real(ovlpst(nst_j,nst_k))
     end do  ! nst_k 
  end do  ! nst_j
!  if (myid == 0 .AND. test_log) write(*,'("")')

!!$ Calculating <v_j|H_T|v_k>
!!$ First calculating <u_j|H_T|v_k> for j >= k 
!!$ Same ideas as calculating overlap elements <u_j|v_k>
  do nst_j = 1, Nmax_1el_basis
     do nst_k = 1, nst_j
        sum1 = 0d0
        do nst_n = 1, nst_k - 1
           sum1 = sum1 + ovlpst(nst_k,nst_n) * e1me(nst_j,nst_n) / ovlpst(nst_n,nst_n)
        end do ! nst_n
     e1me(nst_j,nst_k) = e1me(nst_j,nst_k) - sum1
!     if (myid == 0 .AND. test_log) write(*,'("j,k =",2I3,", <u_j|H_T|v_k> =",F10.5)') nst_j, nst_k, real(e1me(nst_j,nst_k))
     end do  ! nst_k 
  end do  ! nst_j
!  if (myid == 0 .AND. test_log) write(*,'("")')
!!$ Now <v_j|H_T|v_k> 
!!$ <v_j| = <u_j| - sum^{j-1}_{n=1} <u_j|v_n> <v_n| / <v_n|v_n>
!!$ <v_j| = <u_j|H_T|v_k> - sum^{j-1}_{n=1} <u_j|v_n> <v_n|H_T|v_k> / <v_n|v_n>
!!$ hence <v_j|H_T|v_k> = <u_j|H_T|v_k> - sum^{j-1}_{n=1} <u_j|v_n> <v_n|H_T|v_k> / <v_n|u_n>
  do nst_j = 1, Nmax_1el_basis
     do nst_k = 1, nst_j
        sum1 = 0d0
        do nst_n = 1, nst_j - 1
           sum1 = sum1 + ovlpst(nst_j,nst_n) * e1me(nst_n,nst_k) / ovlpst(nst_n,nst_n)
        end do ! nst_n
     e1me(nst_j,nst_k) = e1me(nst_j,nst_k) - sum1
     e1me(nst_k,nst_j) = e1me(nst_j,nst_k)
     if (myid == 0 .AND. test_log) write(*,'("j,k =",2I3,", <v_j|H_T|v_k> =",F10.5)') nst_j, nst_k, real(e1me(nst_j,nst_k))
     end do  ! nst_k 
  end do  ! nst_j
  if (myid == 0 .AND. test_log) write(*,'("")')

!!$ Storing e1me = <vb_j|H_T|vb_k> and ovlpst = <vb_j|vb_k>
!!$ vb_j - new set of orthonormal functions

!!$ Calculating <vb_j|H_T|vb_k> = N_j.N_k.<v_j|H_T|v_k>
!!$ = <v_j|H_T|v_k> / sqrt(<v_j|v_j> <v_k|v_k>)
!!$ = <v_j|H_T|v_k> / sqrt(<v_j|u_j> <v_k|u_k>)
  do nst_j = 1, Nmax_1el_basis
     do nst_k = 1, nst_j
        e1me(nst_j,nst_k) = e1me(nst_j,nst_k) / dsqrt(ovlpst(nst_j,nst_j)*ovlpst(nst_k,nst_k))
        e1me(nst_k,nst_j) = e1me(nst_j,nst_k)
        if (myid == 0 .AND. test_log) write(*,'("j,k =",2I3,", <vb_j|H_T|vb_k> =",F10.5)') nst_j, nst_k,real(e1me(nst_j,nst_k))
     end do  ! nst_k 
  end do  ! nst_j
  if (myid == 0 .AND. test_log) write(*,'("")')

!!$ Must have ovlpst(j,k) = <v_j|u_k> to
!!$ construct orthonormal one-electron (molecular orbital) basis functions
  call Orthonormalise_one_electron_basis_func(TargetStates1el,bst)

!!$ Create new and copy one-electron basis for the Projection Operator I_0 
  call new_basis_st(TargState_Orthonorm,Nmax_1el_basis,.true.,data_in%calculation_type)
  do nst_j = 1, Nmax_1el_basis
     call copy_st(TargState_Orthonorm%b(nst_j),TargetStates1el%b(nst_j))
  end do
  Nmax_SturmBasis = basis_size(bst)
  call new_basis(bst_orth,Nmax_SturmBasis)
  call copy_basis_nr(bst_orth,bst)


!  Nmax_SturmBasis = basis_size(bst)
!  call new_basis(bst_orth,Nmax_SturmBasis)
!  call copy_basis_nr(bst_orth,bst)

!!$ Testing <u_j|v_k>, <v_j|v_k>, <u_j|H_T|v_k>, <vb_j|vb_k>
!!$ <v_j|H_T|v_k>, <vb_j|H_T|vb_k> 
!!$ Testing construction of orthonormal functions
  if (test_log ) then
!!$ Orthnormal function |vb>
!!$ <vb_j| = N_k.<v_j|, N_k = 1 / sqrt(<v_k|v_k>)
!!$ hence <vb_j| = N_j.<u_j| - N_j.sum^{j-1}_{n=1} <u_j|v_n> <vb_n| / sqrt(<v_n|v_n>)
!!$ <vb_j|vb_k> = N_j.N_k.<u_j|v_k> - N_j.sum^{k-1}_{n=1} <v_n|u_k> <vb_n|vb_k> / sqrt(<v_n|v_n>)
     do nst_j = 1, Nmax_1el_basis
        do nst_k = 1, nst_j
           sum1 = 0d0
           do nst_n = 1, nst_j - 1
              sum1 = sum1 + ovlpst(nst_j,nst_n) * ovlpst(nst_n,nst_k) / dsqrt(ovlpst(nst_n,nst_n))
           end do ! nst_n
           ovlpst(nst_j,nst_k) = (ovlpst(nst_j,nst_k) / dsqrt(ovlpst(nst_k,nst_k)) - sum1) / dsqrt(ovlpst(nst_j,nst_j))
           ovlpst(nst_k,nst_j) = ovlpst(nst_j,nst_k)
           if (myid == 0 .AND. test_log) write(*,'("j,k =",2I3,", <vb_j|vb_k>=",F10.5)') nst_j, nst_k,real(ovlpst(nst_j,nst_k))
        end do  ! nst_k 
     end do  ! nst_j
     call Test_orthonormalise_one_electron_basis(tmp1el_basis,temp_bst,tmp_ovlpst,tmp_e1me)
  end if

!!$ Forming overlaps <vb_j|vb_k> in ovlpst
  ovlpst(:,:) = 0d0      
  do nst_j = 1, Nmax_1el_basis
     ovlpst(nst_j,nst_j) = 1d0   
  end do

end subroutine Orthonormalise_one_electron_basis

subroutine Orthonormalise_one_electron_basis_func(TargetStates1el,bst)
!!$ Must have ovlpst(j,k) = <u_j|v_k>, 
!!$ where u_j - original molecular orbital
!!$ v_k - new orthogonal molecular orbital


  use input_data
  use grid_radial
  use sturmian_class
  use state_class
  use ovlpste1me

  implicit none
  type(basis_state), intent(inout):: TargetStates1el ! One-electron basis
  type(basis_sturmian_nr), intent(inout):: bst ! Underlying one-electron basis

  integer:: Nmax_1el_basis    ! One-electron basis size
  integer:: Nmax_SturmBasis     ! Underlying one-electron basis size
  integer:: nst_k, nst_n ! One-electron basis indexes 
  integer:: Nk_1el, Nn_1el ! Number of underlying atomic orbitals per one-electron (molecular orbital) basis
  integer:: nk, nn
  integer:: ind_k, ind_n  ! Index to underlying atomic orbitals basis 'bst'
  real*8:: CI_k, CI_n
  real*8, allocatable, dimension(:):: CI_new
  integer, allocatable, dimension(:):: na_new, ma_new
  integer:: l_n
  real*8:: temp

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

! Size of one-electron molecular orbital basis
  Nmax_1el_basis = basis_size_st(TargetStates1el)
! Size of underlying atomic orbital basis
  Nmax_SturmBasis = basis_size(bst)

  allocate(CI_new(1:Nmax_SturmBasis),na_new(1:Nmax_SturmBasis), ma_new(1:Nmax_SturmBasis))

!!$ This may be able to be done better
!!$ Construct indexes to underlying atomic orbital basis
  do nn = 1, Nmax_SturmBasis
     na_new(nn) = nn
  end do

!!$ Orthnormal function |vb>
!!$ |vb_k> = N_k.|v_k>, N_k = 1 / sqrt(<v_k|v_k>) = 1 / sqrt(<v_k|u_k>)
!!$ |vb_k> = N_k.|u_k> - N_k.sum^{k-1}_{n=1} <v_n|u_k> |vb_n> / N_n.<v_n|v_n> 
!!$ |vb_k> = N_k.|u_k> - N_k.sum^{k-1}_{n=1} <v_n|u_k> |vb_n> / sqrt(<v_n|v_n>) 
!!$ |vb_k> = N_k.|u_k> - N_k.sum^{k-1}_{n=1} <v_n|u_k> |vb_n> / sqrt(<v_n|u_n>) 
  do nst_k = 1, Nmax_1el_basis
     Nk_1el = get_nam(TargetStates1el%b(nst_k))
     CI_new(:) = 0d0

     do nk = 1, Nk_1el
        ind_k = get_na(TargetStates1el%b(nst_k),nk,1)
        CI_k = get_CI(TargetStates1el%b(nst_k),nk)
        if ( CI_k == 0d0 ) cycle
        CI_new(ind_k) = CI_k / dsqrt(ovlpst(nst_k,nst_k))
     end do ! nk

     do nst_n = 1, nst_k - 1
        Nn_1el = get_nam(TargetStates1el%b(nst_n))
        temp = ovlpst(nst_k,nst_n) / dsqrt(ovlpst(nst_n,nst_n))
        if ( temp == 0d0 ) cycle

        do nn = 1, Nn_1el ! Loop over atomic orbitals of molecular orbital nst_n
           ind_n = get_na(TargetStates1el%b(nst_n),nn,1)
           CI_n = get_CI(TargetStates1el%b(nst_n),nn)

           if ( CI_n == 0d0 ) cycle
           CI_new(ind_n) = CI_new(ind_n) - temp * CI_n /dsqrt(ovlpst(nst_k,nst_k))
        end do ! nn

     end do ! nst_n

     ma_new = get_ang_mom_proj(TargetStates1el%b(nst_k))

     call modify_CI_state(TargetStates1el%b(nst_k),Nmax_SturmBasis,na_new,CI_new,ma_new)

  end do  ! nst_k

  if ( data_in%inc == 1 ) then 
     call rearrange(bst,get_max_L(bst),TargetStates1el,.FALSE.)
  end if


end subroutine Orthonormalise_one_electron_basis_func



subroutine Test_orthonormalise_one_electron_basis(TargetStates1el,bst,tmp_ovlpst_uu,tmp_e1me_uu )
!!$ MARK: This subroutine passes through the one-electron target state basis
!!$ 'TargetStates1el' and underlying one-electron Sturmian basis 'bst'. 
!!$ The one-electron target state basis 'TargetStates1el' is orthonomalised
!!$ using the Gram-Schmidt procedure.
!!$ Matrix elements e1me and ovlpst are also modified so that we can
!!$ perform two-electron diagonalisation
!!$ and construct two-electron target state with an orthonormal basis.
!!$ The goal is to end up with one-electron orbitals that are orthonormal,
!!$ an essential property for solving non-uniqueness.

  use input_data
  use grid_radial
  use sturmian_class
  use state_class
  use MPI_module

  implicit none
  type(basis_state), intent(inout):: TargetStates1el
  type(basis_sturmian_nr), intent(inout):: bst
  real*8, dimension(TargetStates1el%Nmax,TargetStates1el%Nmax), intent(in)::  tmp_ovlpst_uu, tmp_e1me_uu

  type(basis_state):: tmp1el_basis    ! One-electron basis
  type(basis_sturmian_nr):: bst_orth  ! Underlying one-electron basis  
  integer:: Nmax_1el_basis    ! One-electron basis size
  integer:: Nmax_SturmBasis     ! Underlying one-electron basis size
  integer:: basis_type          ! Spherical or Spheroidal
  integer:: nst_j, nst_k, nst_n ! One-electron basis indexes 
  integer:: Nj_1el, Nk_1el, Nn_1el ! Number of underlying atomic orbitals per one-electron (molecular orbital) basis
  integer:: nj, nk, nn             ! Scroll through atomic orbitals 
  integer:: ind_j, ind_k, ind_n  ! Index to underlying atomic orbitals basis 'bst'
  real*8:: CI_j, CI_k, CI_n         
  real*8, allocatable, dimension(:,:)::  e1me_uu, e1me_uv, e1me_vv  !<..|H_T|..>
  real*8, allocatable, dimension(:,:)::  ovlp_uu, ovlp_uv, ovlp_vv  !<u|u>, <u|v>, <v|v>
  real*8, allocatable, dimension(:):: CI_new
  integer, allocatable, dimension(:):: na_new, ma_new
  integer:: l_n
  real*8:: sum1, sum2, temp 
  
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

!!$ Create new and copy one-electron basis 
!!$ Probably not neeeded
  basis_type = data_in%calculation_type   ! Spherical or spheroidal
  Nmax_1el_basis = basis_size_st(TargetStates1el)
  call new_basis_st(tmp1el_basis,Nmax_1el_basis,.true.,basis_type)
  do nst_j = 1, Nmax_1el_basis
       call copy_st(tmp1el_basis%b(nst_j),TargetStates1el%b(nst_j))
  end do

!!$ Create new and copy underlying one-electron functions
  Nmax_SturmBasis = basis_size(bst)
  call new_basis(bst_orth,Nmax_SturmBasis)
  call copy_basis_nr(bst_orth,bst)

!!$ Creates orthonormal one-electron (TargetStates) basis via Gram-Schmidt
!!$ orthogonalisation.
  if (myid == 0) write(*,'("********************************")') 
  if (myid == 0) write(*,'("     Gram-Schmidt Testing")') 
  if (myid == 0) write(*,'("********************************")') 

  allocate(ovlp_uu(Nmax_1el_basis,Nmax_1el_basis), ovlp_uv(Nmax_1el_basis,Nmax_1el_basis), ovlp_vv(Nmax_1el_basis,Nmax_1el_basis))
  allocate(e1me_uu(Nmax_1el_basis,Nmax_1el_basis), e1me_uv(Nmax_1el_basis,Nmax_1el_basis), e1me_vv(Nmax_1el_basis,Nmax_1el_basis))

!!$ Form overlaps <u_j|u_k>   and <u_j|H_T|u_k> 
!!$ u_j - old nonorthogonal molecular orbital set
!!$ Underlying atomic orbital basis 'bst' overlaps in 'bst%ortint(:,:)'
!!$ and target Hamiltonian matrix elements in bst%ham1el(ind_j,ind_k)  
  do nst_j = 1, Nmax_1el_basis
     Nj_1el = get_nam(TargetStates1el%b(nst_j))
     do nst_k = 1, nst_j
        Nk_1el = get_nam(TargetStates1el%b(nst_k))
        sum1 = 0d0
        sum2 = 0d0
        do nj = 1, Nj_1el
           ind_j = get_na(TargetStates1el%b(nst_j),nj,1)
           CI_j = get_CI(TargetStates1el%b(nst_j),nj)
           do nk =1, Nk_1el
                ind_k = get_na(TargetStates1el%b(nst_k),nk,1)
                CI_k = get_CI(TargetStates1el%b(nst_k),nk)
                sum1 = sum1 + CI_j * CI_k * bst%ortint(ind_j,ind_k)
                sum2 = sum2 + CI_j * CI_k * bst%ham1el(ind_j,ind_k)
           end do ! nk
        end do ! nj
        ovlp_uu(nst_j,nst_k) = sum1
        ovlp_uu(nst_k,nst_j) = ovlp_uu(nst_j,nst_k)  
        e1me_uu(nst_j,nst_k) = sum2
        e1me_uu(nst_k,nst_j) = e1me_uu(nst_j,nst_k)
!!$ Checking that the inital basis underlying atomic orbitals overlaps agree with molecular orbital overlaps
!        if (myid == 0) write(*,'("j,k =",2I3,", <u_j|H_T|u_k> =",F10.5,", e1me_uu(j,k) =",F10.5)')  nst_j, nst_k, real(e1me_uu(nst_j,nst_k)), tmp_e1me_uu(nst_j,nst_k) 
!        if (myid == 0) write(*,'("j,k =",2I3,", <u_j|u_k> =",F10.5,", ovlpst_uu(j,k) =",F10.5)')  nst_j, nst_k, real(ovlp_uu(nst_j,nst_k)), tmp_ovlpst_uu(nst_j,nst_k)
     end do  ! nst_k 
  end do  ! nst_j
!  if (myid == 0) print*,""
!!$ Form overlaps <u_j|u_k>  and <u_j|H_T|u_k>
!!$ u_j - old nonorthogonal molecular orbital set
  e1me_uu(:,:) = tmp_e1me_uu(:,:)
  ovlp_uu(:,:) = tmp_ovlpst_uu(:,:)

!!$ form overlap array <u_j|v_k>
!!$ u_j - old nonorthogonal molecular orbital set,
!!$ v_k - new orthogonal molecular orbital set but not yet normalised
!!$ Only elements with  j >= k  required.
!!$ This algorithm utilises the property <v_k|v_k> = <v_k|u_k>
!!$ |v_k> = |u_k> - sum^{k-1}_{n=1} <v_n|u_k> |v_n> / <v_n|v_n>
!!$ <u_j||v_k> = <u_j||u_k> - sum^{k-1}_{n=1} <v_n|u_k> <u_j||v_n> / <v_n|v_n>
!!$ hence <u_j||v_k> = <u_j||u_k> - sum^{k-1}_{n=1} <v_n|u_k> <u_j||v_n> / <v_n|u_n>
  ovlp_uv(:,:) = ovlp_uu(:,:)
  do nst_j = 1, Nmax_1el_basis
     do nst_k = 1, nst_j
        sum1 = 0d0
        do nst_n = 1, nst_k - 1
           sum1 = sum1 + ovlp_uv(nst_k,nst_n) * ovlp_uv(nst_j,nst_n) / ovlp_uv(nst_n,nst_n) 
        end do ! nst_n
        ovlp_uv(nst_j,nst_k) = ovlp_uu(nst_j,nst_k) - sum1 
!        if (myid == 0) write(*,'("j,k =",2I3,", <u_j|v_k> =",F10.5)') nst_j, nst_k, real(ovlp_uv(nst_j,nst_k))
     end do  ! nst_k 
  end do  ! nst_j
!  if (myid == 0) print*,""

!!$ Calculating <v_j|H_T|v_k>
!!$ First calculating <u_j|H_T|v_k> for j >= k 
!!$ Same ideas as calculating overlap elements <u_j|v_k>
  e1me_uv(:,:) = e1me_uu(:,:)
  do nst_j = 1, Nmax_1el_basis
     do nst_k = 1, nst_j
        sum1 = 0d0
        do nst_n = 1, nst_k - 1
           sum1 = sum1 + ovlp_uv(nst_k,nst_n) * e1me_uv(nst_j,nst_n) / ovlp_uv(nst_n,nst_n)
        end do ! nst_n
        e1me_uv(nst_j,nst_k) = e1me_uu(nst_j,nst_k) - sum1
!        if (myid == 0) write(*,'("j,k =",2I3,", <u_j|H_T|v_k> =",F10.5)') nst_j, nst_k, real(e1me_uv(nst_j,nst_k))
     end do  ! nst_k 
  end do  ! nst_j
!  if (myid == 0) print*,""
!!$ Now <v_j|H_T|v_k> 
!!$ <v_j| = <u_j| - sum^{j-1}_{n=1} <u_j|v_n> <v_n| / <v_n|v_n>
!!$ <v_j|H_T|v_k> = <u_j|H_T|v_k> - sum^{j-1}_{n=1} <u_j|v_n> <v_n|H_T|v_k> / <v_n|v_n>
!!$ hence <v_j|H_T|v_k> = <u_j|H_T|v_k> - sum^{j-1}_{n=1} <u_j|v_n> <v_n|H_T|v_k> / <v_n|u_n>
  e1me_vv(:,:) = e1me_uv(:,:)
  do nst_j = 1, Nmax_1el_basis
     do nst_k = 1, nst_j
        sum1 = 0d0
        do nst_n = 1, nst_j - 1
           sum1 = sum1 + ovlp_uv(nst_j,nst_n) * e1me_vv(nst_n,nst_k) / ovlp_uv(nst_n,nst_n)
        end do ! nst_n
        e1me_vv(nst_j,nst_k) = e1me_uv(nst_j,nst_k) - sum1
        e1me_vv(nst_k,nst_j) = e1me_vv(nst_j,nst_k)
        if (myid == 0) write(*,'("j,k =",2I3,", <v_j|H_T|v_k> =",F10.5)') nst_j, nst_k, real(e1me_vv(nst_j,nst_k))
     end do  ! nst_k 
  end do  ! nst_j
  if (myid == 0) print*,""


!!$ Calculating <v_j|v_k>
!!$ Same algorithm as <v_j|H_T|v_k>
  ovlp_vv(:,:) = ovlp_uv(:,:)
  do nst_j = 1, Nmax_1el_basis
     do nst_k = 1, nst_j
        sum1 = 0d0
        do nst_n = 1, nst_j - 1
           sum1 = sum1 + ovlp_uv(nst_j,nst_n) * ovlp_vv(nst_n,nst_k) / ovlp_uv(nst_n,nst_n)
        end do ! nst_n
        ovlp_vv(nst_j,nst_k) = ovlp_uv(nst_j,nst_k) - sum1
        ovlp_vv(nst_k,nst_j) = ovlp_vv(nst_j,nst_k)
        if (myid == 0 ) write(*,'("j,k =",2I3,", <v_j|v_k> =",F10.5)') nst_j, nst_k, real(ovlp_vv(nst_j,nst_k))
     end do  ! nst_k 
  end do  ! nst_j
  if (myid == 0 ) print*,""

!!$ vb_j - new set of orthonormal functions
!!$ Calculating <vb_j|H_T|vb_k> = N_j.N_k.<v_j|H_T|v_k>
!!$ = <v_j|H_T|v_k> / sqrt(<v_j|v_j> <v_k|v_k>)
  do nst_j = 1, Nmax_1el_basis
     do nst_k = 1, nst_j
        e1me_vv(nst_j,nst_k) = e1me_vv(nst_j,nst_k) / dsqrt(ovlp_vv(nst_j,nst_j)*ovlp_vv(nst_k,nst_k))
        e1me_vv(nst_k,nst_j) = e1me_vv(nst_j,nst_k)
        if (myid == 0 ) write(*,'("j,k =",2I3,", <vb_j|H_T|vb_k> =",F10.5)') nst_j, nst_k,real(e1me_vv(nst_j,nst_k))
     end do  ! nst_k 
  end do  ! nst_j

!!$ Construct orthonormal molecular orbital set vb_k with
!!$ |v_k> = |u_k> - sum^{k-1}_{n=1} <v_n|u_k> |v_n> / <v_n|v_n>
!!$ |vb_k> = |v_k> / sqrt(<v_k|v_k>)

  allocate(CI_new(1:Nmax_SturmBasis),na_new(1:Nmax_SturmBasis),ma_new(1:Nmax_SturmBasis))

!!$ Constructs an orthnormalised one-electron (molecular orbital) basis

!!$ Construct indexes to underlying atomic orbital basis
  do nn = 1, Nmax_SturmBasis
     na_new(nn) = nn
  end do 

!!$ Orthomal function |vb>
!!$ |vb_k> = N_k.|v_k>, N_k = 1 / sqrt(<v_k|v_k>)
!!$ |vb_k> = N_k.|u_k> - N_k.sum^{k-1}_{n=1} <v_n|u_k> |vb_n> / N_n.<v_n|v_n> 
!!$ |vb_k> = N_k.|u_k> - N_k.sum^{k-1}_{n=1} <v_n|u_k> |vb_n> / sqrt(<v_n|v_n>) 
  do nst_k = 1, Nmax_1el_basis 
     Nk_1el = get_nam(TargetStates1el%b(nst_k))
     CI_new(:) = 0d0

     do nk = 1, Nk_1el
        ind_k = get_na(TargetStates1el%b(nst_k),nk,1)
        CI_k = get_CI(TargetStates1el%b(nst_k),nk)
        if ( CI_k == 0d0 ) cycle
        CI_new(ind_k) = CI_k / dsqrt(ovlp_vv(nst_k,nst_k))
     end do ! nk
        
     do nst_n = 1, nst_k - 1
        Nn_1el = get_nam(tmp1el_basis%b(nst_n))
        temp = ovlp_uv(nst_k,nst_n) / dsqrt(ovlp_vv(nst_n,nst_n))
        if ( temp == 0d0 ) cycle

        do nn = 1, Nn_1el ! Loop over atomic orbitals of molecular orbital nst_n
           ind_n = get_na(tmp1el_basis%b(nst_n),nn,1)
           CI_n = get_CI(tmp1el_basis%b(nst_n),nn)

           if ( CI_n == 0d0 ) cycle
           CI_new(ind_n) = CI_new(ind_n) - temp * CI_n / dsqrt(ovlp_vv(nst_k,nst_k)) 
        end do ! nn
     end do ! nst_n
     ma_new = get_ang_mom_proj(tmp1el_basis%b(nst_k))
     call modify_CI_state(tmp1el_basis%b(nst_k),Nmax_SturmBasis,na_new,CI_new,ma_new)
  end do  ! nst_k

  if ( data_in%inc == 1 ) then 
     call rearrange(bst,get_max_L(bst),tmp1el_basis,.FALSE.)
  end if

  if (myid == 0 ) print*,"********************"
  if (myid == 0 ) print*,"      REARRANGE"
  if (myid == 0 ) print*,"********************"
!!$ Calculate <vb_j|H_T|vb_k> with the orthonormalised one-electron (molecular orbital)
!!$ underlying atomic orbital basis. If this matches with <vb_j|H_T|vb_k>
!!$ printed above, this means orthonormal molecular orbital basis correctly stored
!!$ and calculated.
  do nst_j = 1, Nmax_1el_basis
     Nj_1el = get_nam(tmp1el_basis%b(nst_j))
     do nst_k = 1, nst_j
        Nk_1el = get_nam(tmp1el_basis%b(nst_k))
        sum1 = 0d0
        sum2 = 0d0
        do nj = 1, Nj_1el
           ind_j = get_na(tmp1el_basis%b(nst_j),nj,1)
           CI_j = get_CI(tmp1el_basis%b(nst_j),nj)
           do nk =1, Nk_1el
                ind_k = get_na(tmp1el_basis%b(nst_k),nk,1)
                CI_k = get_CI(tmp1el_basis%b(nst_k),nk)
                sum1 = sum1 + CI_j * CI_k * bst%ortint(ind_j,ind_k)
                sum2 = sum2 + CI_j * CI_k * bst%ham1el(ind_j,ind_k)
           end do ! nk
        end do ! nj
        ovlp_vv(nst_j,nst_k) = sum1
        ovlp_vv(nst_k,nst_j) = ovlp_vv(nst_j,nst_k)
        e1me_vv(nst_j,nst_k) = sum2
        e1me_vv(nst_k,nst_j) = e1me_vv(nst_j,nst_k)
        if (myid == 0) write(*,'("j,k =",2I3,", <vb_j|H_T|vb_k> =",F10.5)') nst_j, nst_k, real(e1me_vv(nst_j,nst_k))
        if (myid == 0) write(*,'("j,k =",2I3,", <vb_j|vb_k> =",F10.5)') nst_j, nst_k, real(ovlp_vv(nst_j,nst_k))
     end do  ! nst_k 
  end do  ! nst_j
  if (myid == 0) print*,""
!  stop

end subroutine Test_orthonormalise_one_electron_basis 

subroutine Test_Overlap_st(nstp,nst)

  use sturmian_class
  use target_states
  use grid_radial
  use target_states
  use one_electron_func
  use input_data

  implicit none 
  integer, intent(in)::nstp, nst

!!$ TESTING
  integer:: nconfp, nconf
  integer:: n1,n2,np1,np2, i1, i2, ip1, ip2
  integer:: ind1, ind2, indp1, indp2
  type(sturmian_nr), pointer::tn1, tn2, tnp1, tnp2  ! pointers to one-electron atomic orbitals
  real*8:: CI_stp, CI_st, CI_1, CI_2, CI_1p, CI_2p ! 2e configuration and 1e molecular orbial coefficients
  real*8, pointer, dimension(:)::  fp1, f1, fp2, f2  ! underlying atomic orbital radial functions
  integer:: minf1, maxf1, minfp1, maxfp1
  integer:: minf2, maxf2, minfp2, maxfp2
  integer:: ir1, ir2
  real*8:: overlap1, overlap2, overlap12
  real*8:: r_ovlp1,r_ovlp12
  integer:: Lap1,Mp1,Lap2,Mp2,La1,M1,La2,M2,Mp_st,M_st  ! A.O. QN and Molecular State M 
  integer:: Spin, Spinp

  M_st = NINT(TargetStates2el%b(nst)%M)       
  Mp_st = NINT(TargetStates2el%b(nstp)%M )      ! Molecular state Angular Projection
  Spinp = NINT(TargetStates2el%b(nstp)%spin) ! 2e Molecular State Spin
  Spin = NINT(TargetStates2el%b(nst)%spin)

  if ( Mp_st /= M_st .OR. Spinp /= Spin ) return

  overlap12 = 0d0 
  r_ovlp12 = 0d0 

!!$     <ip1,ip2|i1,i2>    
  print*,""
  print*,"Overlap contributions <nstp,nst>",nstp,nst
  do nconfp = 1 , TargetStates2el%b(nstp)%nam
     CI_stp = get_CI(TargetStates2el%b(nstp),nconfp)
     np1 = get_na(TargetStates2el%b(nstp),nconfp,1)
     np2 = get_na(TargetStates2el%b(nstp),nconfp,2)
     Mp1 = get_ang_mom_proj(TargetStates%b(np1))      ! Angular momentum projection 1e molecular orbital
     Mp2 = get_ang_mom_proj(TargetStates%b(np2))      ! Angular momentum projection 1e molecular orbital

     do nconf = 1, TargetStates2el%b(nst)%nam
        CI_st = get_CI(TargetStates2el%b(nst),nconf)
        n1 = get_na(TargetStates2el%b(nst),nconf,1)
        n2 = get_na(TargetStates2el%b(nst),nconf,2)
        M1 = get_ang_mom_proj(TargetStates%b(n1))      ! Angular momentum projection 1e molecular orbital
        M2 = get_ang_mom_proj(TargetStates%b(n2)) 

        if ( Mp1 /= M1 .OR. Mp2 /= M2 ) cycle

        do ip1 = 1, get_nam(TargetStates%b(np1))
           indp1 = get_na(TargetStates%b(np1),ip1)
           CI_1p = get_CI(TargetStates%b(np1),ip1)
           tnp1 => bst_nr%b(indp1)   !           
           fp1 => fpointer(tnp1)     ! One electron functions
           Lap1 = get_ang_mom(tnp1)
           minfp1 = get_minf(tnp1)
           maxfp1 = get_maxf(tnp1)

           do i1 = 1, get_nam(TargetStates%b(n1))
              ind1 = get_na(TargetStates%b(n1),i1)
              CI_1 = get_CI(TargetStates%b(n1),i1)
              tn1 => bst_nr%b(ind1)
              f1 => fpointer(tn1)     ! One electron functions
              La1 = get_ang_mom(tn1)
              minf1 = get_minf(tn1)
              maxf1 = get_maxf(tn1)
              ir1 = max(minf1,minfp1)
              ir2 = min(maxf1,maxfp1)
              if ( Lap1 /= La1 ) cycle
              overlap1 = CI_1p * CI_1 * SUM( fp1(ir1:ir2) *  f1(ir1:ir2) * grid%weight(ir1:ir2))
              r_ovlp1 = CI_1p * CI_1 * SUM( fp1(ir1:ir2) * grid%gridr(ir1:ir2) * f1(ir1:ir2) * grid%weight(ir1:ir2))  
              if ( n1==1 .AND. np1 ==1 ) then
                   write(*,'("np1,n1,np2,n2,<np1|n1>,          "4I3,F12.7)') np1,n1,np2,n2, overlap1
                   write(*,'("np1,n1,np2,n2,<np1|r|n1>,        "4I3,F12.7)') np1,n1,np2,n2, r_ovlp1
              end if
              if (  n1 /= 1 ) then
                if (np1 == 1) then
                   write(*,'("np1,n1,np2,n2,Cst.<np1|n1>,      "4I3,F12.7)') np1,n1,np2,n2,CI_st*overlap1
                   write(*,'("np1,n1,np2,n2,Cst.<np1|r|n1>,    "4I3,F12.7)') np1,n1,np2,n2,CI_st*r_ovlp1
                else if ( np1 /= 2 ) then
                   write(*,'("np1,n1,np2,n2,Cst.Cstp<np1|n1>,  "4I3,F12.7)') np1,n1,np2,n2,CI_st*CI_stp*overlap1
                   write(*,'("np1,n1,np2,n2,Cst.Cstp<np1|r|n1>,"4I3,F12.7)') np1,n1,np2,n2,CI_st*CI_stp*r_ovlp1
                end if
              end if

              do ip2 = 1, get_nam(TargetStates%b(np2))
                 indp2 = get_na(TargetStates%b(np2),ip2) ! Index to underlying atomic orbitals 
                 CI_2p = get_CI(TargetStates%b(np2),ip2)
                 tnp2 => bst_nr%b(indp2)
                 fp2 => fpointer(tnp2)
                 Lap2 = get_ang_mom(tnp2)
                 minfp2 = get_minf(tnp2)
                 maxfp2 = get_maxf(tnp2)

                 do i2 = 1, get_nam(TargetStates%b(n2))
                    ind2 = get_na(TargetStates%b(n2),i2)
                    CI_2 = get_CI(TargetStates%b(n2),i2)
                    tn2 => bst_nr%b(ind2)
                    f2 => fpointer(tn2)     ! One electron functions
                    La2 = get_ang_mom(tn2)
                    minf2 = get_minf(tn2)
                    maxf2 = get_maxf(tn2)
                    ir1 = max(minf2,minfp2)
                    ir2 = min(maxf2,maxfp2)
                    if ( Lap2 /= La2 ) cycle
                    overlap2 = CI_2p * CI_2 * SUM( fp2(ir1:ir2) *  f2(ir1:ir2) * grid%weight(ir1:ir2))
                    overlap12 = overlap12 + CI_st * CI_stp * overlap1 * overlap2
                    r_ovlp12 = r_ovlp12 +  CI_st * CI_stp * r_ovlp1 * overlap2
                 end do ! i2
              end do ! ip2  
           end do ! i1 
        end do ! ip1
     end do ! nst
  end do ! nstp
  write(*,'("nstp,nst, <nstp|nst>,  "2I3,F12.5)') nstp,nst, overlap12
  write(*,'("nstp,nst, <nstp|r|nst>,"2I3,F12.5)') nstp,nst, r_ovlp12


end subroutine Test_Overlap_st
