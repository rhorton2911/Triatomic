
!-=-= OPTICAL_POTENTIAL MODULE =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!     AUTHOR: LIAM SCARLETT
!
!     DESCRIPTION: SUBROUTINES FOR CALCULATING OPTICAL POTENTIALS

      module optical_potential
      
        !PUBLIC VARIABLES
          logical, public :: complex_optical
          real*8, dimension(:,:), allocatable :: WMat
          complex*16, dimension(:,:), allocatable :: CWMat
          integer :: nst_P, nchm_P, nchm_open_P, nkgmax_P

        !PUBLIC SUBROUTINES
          public construct_optical_potential

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

contains

!-=-= CONSTRUCT_OPTICAL_POTENTIAL =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !
  ! This subroutine takes a pre-calculated V-matrix and calculated the corresponding optical potential matrix 
  ! for a specified number of channels <nchm_P> in P space
  !
  !The weak-coupling optical potential is calculated directly:
  !
  ! <f|W|i> = ∑_m <f|V|m> w_m <m|V|i>
  !
  ! where f and i are states/channels/k points in P space
  ! and m is the same in Q space
  ! w_m are the Green's function & integration weights.
  ! 
  !
  ! The full optical potential is calculated by solving the
  ! matrix equation:
  !       P        Q
  !    -------------------- -------     -------
  !    |     |            | |     |     |     |
  !  P |  I  |      -WV   | |  W  |     | WWC |
  !    |     |            | |     |     |     |
  !    -------------------- -------     -------
  !    |     |            | |     |     |     |
  !    |     |            | |     |  =  |     |
  !    |     |            | |     |     |     |
  !  Q |  0  |   I - WV   | |  X  |     | XWC |
  !    |     |            | |     |     |     |
  !    |     |            | |     |     |     |
  !    -------------------- -------     -------
  !
  ! The output matrix WMAT contains the optical potential PLUS P-space V-matrix
  !
  subroutine construct_optical_potential(nchm, nchm_open, nchm_P, nkgmax, nkgmax_P, gridk_onopen, gfw_in, npk)
    use input_data
    use vmat_module 
    implicit none

    !-=-= INPUT/OUTPUT
    integer,            intent(in) :: nchm      !total number of channels (P+Q)
    integer,            intent(in) :: nchm_open !number of open channels (P+Q)
    integer,            intent(in) :: nchm_P    !number of channels to keep in P space
    integer,            intent(in) :: nkgmax    !number of total k-grid points (P+Q)
    integer,            intent(in) :: nkgmax_P  !number of k-grid points in P space
    integer,     dimension(nchm+1),             intent(in)              :: npk 
    real*8,    dimension(nkgmax),             intent(in)              :: gfw_in
    real*8,    dimension(nchm_open),             intent(in)              :: gridk_onopen

    !-=-= LOCAL VARIABLES
    integer :: nch, nchf, nchi, nkf, nki, nkm
    real*8, dimension(nkgmax) :: gfw
    real*8, dimension(:,:), allocatable :: Xmat
    complex*16, dimension(:,:), allocatable :: CXmat
    integer, dimension(:), allocatable :: ipiv  
    integer :: info, nunit, nkgmax_X
    character(len=:), allocatable :: filename
    character(len=3) :: nchflabel, nchilabel
    real*8 :: temp, tt
    complex*16 :: tempc
    logical :: complex_input
    real*8, dimension(nkgmax_P,nkgmax_P) :: vmat_P
    complex*16, dimension(nkgmax_P,nkgmax_P) :: Cvmat_P
    real*8, parameter :: pi = acos(-1.0d0)
    logical :: weak_coupling

    print*, '  complex_optical:', complex_optical
    print*, '>>> optical:', vmat(npk(1),npk(2)),vmat(npk(2),npk(1))

    weak_coupling = data_in%weak_coupling

    if(weak_coupling) then !Xmat will just be WWC 
      nkgmax_X = nkgmax_P
    else                   !Xmat will be WC + XWC
      nkgmax_X = nkgmax 
    
      !Alocate matrix for full optical potential matrix
      if(complex_optical) then
        allocate(CWmat(nkgmax_P,nkgmax_P))
      else
        allocate(Wmat(nkgmax_P,nkgmax_P))
      endif

    endif

    !Save the P-space portion of the V matrix to add to the optical potential at the end
    if(complex_optical) then
      CVmat_P = CVmat(1:nkgmax_P,1:nkgmax_P)
    else
      Vmat_P = Vmat(1:nkgmax_P,1:nkgmax_P)
    endif

    !Allocate the "X" matrix (RHS) of matrix equation
    if(complex_optical) then
      allocate(CXmat(nkgmax_X,nkgmax_P))
      CXmat = 0.0d0
    else
      allocate(Xmat(nkgmax_X,nkgmax_P))
      Xmat = 0.0d0
    endif

    !GFW array is equal to gfw_in but for open onshell points we add the on-shell subtraction weight
    gfw = gfw_in
    do nch=1, nchm
      if(nch<=nchm_open) gfw(npk(nch)) = -pi/gridk_onopen(nch)
    enddo
  
!    !calculate weak-coupling optical potential
!    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(temp,tempc,nki,nkf,nkm)
    do nki=1, nkgmax_P
      do nkf=1, nkgmax_X
       
        !if(nkf <= nkgmax_P .and. nkf > nki) cycle
        temp = 0.0d0
        tempc = 0.0d0
        do nch=nchm_P+1, nchm
          nkm = npk(nch)
          if(nch<=nchm_open) then !also means complex_optical = .true.
            !on-shell subtraction if there are open channels in Q space
            !Factor of (0,1)=sqrt(-1) is here because we define Green's function as real
            if(weak_coupling) then
              tempc = tempc +  Vmat(nkm,nkf)*gfw(nkm)* Vmat(nkm,nki)*(0.0d0,1.0d0)
            else
              tempc = tempc + CVmat(nkm,nkf)*gfw(nkm)*CVmat(nkm,nki)*(0.0d0,1.0d0)
            endif
          endif
          do nkm=npk(nch)+1, npk(nch+1)-1 
            
            if(complex_optical .and. .not.weak_coupling) then
              tempc = tempc + CVmat(nkm,nkf)*gfw(nkm)*CVmat(nkm,nki)
            else
              temp  = temp  +  Vmat(nkm,nkf)*gfw(nkm)* Vmat(nkm,nki)
            endif

          enddo !nkm
        enddo !nch
        
        if(complex_optical) then
          CXmat(nkf,nki) = tempc
!          if(nkf<=nkgmax_P) CXmat(nki,nkf) = tempc
        else
          Xmat(nkf,nki) = temp
!          if(nkf<=nkgmax_P) Xmat(nki,nkf) = temp
        endif
      
      enddo !nkf
    enddo !nki
!    !$OMP END PARALLEL DO
        
    !Write the weak-coupling matrix to file
    if(data_in%print_vmatrix) then
      if(complex_optical) then
        call write_matrix_comp(CXmat, nst_P, nchm_P, 'VW-WC')
      else
        call write_matrix(Xmat, nst_P, nchm_P, 'VW-WC')
      endif
    endif

    if(weak_coupling) then
      deallocate(Vmat) !Vmat always allocated for weak-coupling
      if(complex_optical) then
        allocate(CVmat(nkgmax_P,nkgmax_P))
        CVmat = CVmat_P + CXmat
      else
        allocate(Vmat(nkgmax_P,nkgmax_P))
        Vmat = Vmat_P + Xmat
      endif

      print*, 'VP:', Vmat_P(1,1)
      print*, 'X:', Xmat(1,1)
      print*, 'V:', Vmat(1,1)



      return

    endif


    print*, 'MAKING FULL OPTICAL POTENTIAL'


  !Prepare Vmat and Xmat for calculating full optical potential
  do nkf=1, nkgmax
    do nki=1, nkgmax
      if(complex_optical) then
        CVmat(nkf,nki) = CVmat(nkf,nki) * sqrt(abs(gfw(nkf)*gfw(nki)))
        if(nki<=nkgmax_P) CXmat(nkf,nki) = CXmat(nkf,nki) * sqrt(abs(gfw(nkf)))
      else
        Vmat(nkf,nki) = Vmat(nkf,nki) * sqrt(abs(gfw(nkf)*gfw(nki)))
        if(nki<=nkgmax_P) Xmat(nkf,nki) = Xmat(nkf,nki) * sqrt(abs(gfw(nkf)))
      endif
    enddo
  enddo
  
  if(complex_optical) then 
      
    !absorb factor of i=sqrt(-1) from green's function into Vmat
    do nch=nchm_P+1,nchm_open
      nki = npk(nch)
      do nkf=1,nkgmax 
        CVmat(nkf,nki) = cmplx(0.0d0,1.0d0)*CVmat(nkf,nki)
      enddo
    enddo
      
    CVmat = -CVmat
    do nkf=1, nkgmax
      do nki=1,nkgmax_P
        CVmat(nkf,nki) = 0.0d0
      enddo
      CVmat(nkf,nkf) = sign(1.0d0,gfw(nkf)) + CVmat(nkf,nkf)
    enddo
    
  else
    Vmat = -Vmat
    do nkf=1, nkgmax
      do nki=1,nkgmax_P
        Vmat(nkf,nki) = 0.0d0
      enddo
      Vmat(nkf,nkf) = sign(1.0d0,gfw(nkf)) + Vmat(nkf,nkf)
    enddo

  endif !complex optical

  allocate(ipiv(nkgmax))
  ipiv = 0 !initialise to zero to avoid DIVZERO error when compiling with -ftrapuv
  info = 0

  if(complex_optical) then
    call ZGESV(nkgmax,nkgmax_P,CVmat,nkgmax,ipiv,CXmat,nkgmax,info)
  else
    call DGESV(nkgmax,nkgmax_P,Vmat,nkgmax,ipiv,Xmat,nkgmax,info)
  endif
  deallocate(ipiv)

  do nki=1,nkgmax_P 
    do nkf=1, nkgmax_P
      if(gfw(nkf)/=0.0d0) then
        if(complex_optical) then
          CWmat(nkf,nki) = CXmat(nkf,nki)/sqrt(abs(gfw(nkf)))*sign(1.0d0,gfw(nkf))
        else
          Wmat(nkf,nki) = Xmat(nkf,nki)/sqrt(abs(gfw(nkf)))*sign(1.0d0,gfw(nkf))
        endif
      else
        if(complex_optical) then
          CWmat(nkf,nki) = 0.0d0
        else
          Wmat(nkf,nki) = 0.0d0
        endif
      endif
    enddo
  enddo

  if(data_in%print_vmatrix) then
    if(complex_optical) then
      call write_matrix_comp(CWmat, nst_P, nchm_P, 'VW-FULL')
    else
      call write_matrix(Wmat, nst_P, nchm_P, 'VW-FULL')
    endif
  endif

  !Add the optical matrix to the P-space V matrix
  if(complex_optical) then
    deallocate(CVmat)
    allocate(CVmat(nkgmax_P,nkgmax_P))
    CVmat = CWmat + CVmat_P
  else
    deallocate(Vmat)
    allocate(Vmat(nkgmax_P,nkgmax_P))
    Vmat = Wmat + Vmat_P
  endif
   
  deallocate(Wmat)

end subroutine construct_optical_potential

subroutine write_matrix(matrix, nst, nch, label)
  use channels
  use target_states
  use kgrid
  implicit none
  
  real*8, dimension(:,:), intent(in) :: matrix
  integer, intent(in) :: nst, nch
  character(*), intent(in) :: label
  character*20 :: filename
  integer, dimension(nst) :: sst, mst, parst
  integer :: nkmax, j, nkf, nki
  real*8, dimension(nst) :: en_est
  character(len=:), allocatable :: nstchar, nchchar, nkchar
 
  nkmax = npk(nch+1)-1

  nstchar = '0000000000'
  write(nstchar,'(I0)') nst
  nstchar = trim(adjustl(nstchar))
  nchchar = '0000000000'
  write(nchchar,'(I0)') nch
  nchchar = trim(adjustl(nchchar))
  nkchar = '0000000000'
  write(nkchar,'(I0)') nkmax
  nkchar = trim(adjustl(nkchar))


  do j=1, nst
    mst(j) = get_ang_mom_proj(TargetStates2el%b(j))
    parst(j) = get_par(TargetStates2el%b(j))
    sst(j) = 2*get_spin(TargetStates2el%b(j))+1
    en_est(j) = get_energy(TargetStates2el%b(j))
  enddo

  write(filename,'(A,"-",A3)') trim(adjustl(label)), symlabel
  open(unit=605,file=trim(adjustl(filename)),action='write',status='replace')
  write(605,*) nkmax, nch, symlabel
  write(605,'('//nstchar//'(I0,X))') (mst(j), j=1, nst)
  write(605,'('//nstchar//'(I0,X))') (parst(j), j=1, nst)
  write(605,'('//nstchar//'(I0,X))') (sst(j), j=1, nst)
  write(605,'('//nchchar//'(I0,X))') (st_ch(j),j=1,nch)
  write(605,'('//nchchar//'(I0,X))') (Lp_ch(j),j=1,nch)
  write(605,'('//nchchar//'(I0,X))') (Mp_ch(j),j=1,nch)
  write(605,'('//nstchar//'(ES15.8,X))') (en_est(j),j=1,nst)
  write(605,'('//nkchar//'(ES15.8,X))') (gridk(1:npk(j+1)-npk(j),j), j=1,nch)
  write(605,'('//nkchar//'(ES15.8,X))') (weightk(1:npk(j+1)-npk(j),j), j=1,nch)
  write(605,'('//nchchar//'(I0,X),I0)') npk(1:nch+1)
  do nkf=1,nkmax
    write(605,'('//nkchar//'("(",ES15.8,",",ES15.8,")",X))') (matrix(nkf,nki), 0.0d0, nki=1,nkmax)
  enddo
  close(605)
end subroutine write_matrix

subroutine write_matrix_comp(matrix, nst, nch, label)
  use channels
  use target_states
  use kgrid
  implicit none
  complex*16, dimension(:,:), intent(in) :: matrix
  integer, intent(in) :: nst, nch
  character*2, intent(in) :: label
  character*20 :: filename
  integer, dimension(nst) :: sst, mst, parst
  real*8, dimension(nst) :: en_est

  integer :: nkmax, j, nkf, nki
  character(len=:), allocatable :: nstchar, nchchar, nkchar
  
  nkmax = npk(nch+1)-1

  nstchar = '0000000000'
  write(nstchar,'(I0)') nst
  nstchar = trim(adjustl(nstchar))
  nchchar = '0000000000'
  write(nchchar,'(I0)') nch
  nchchar = trim(adjustl(nchchar))
  nkchar = '0000000000'
  write(nkchar,'(I0)') nkmax
  nkchar = trim(adjustl(nkchar))

  do j=1, nst
    mst(j) = get_ang_mom_proj(TargetStates2el%b(j))
    parst(j) = get_par(TargetStates2el%b(j))
    sst(j) = 2*get_spin(TargetStates2el%b(j))+1
    en_est(j) = get_energy(TargetStates2el%b(j))
  enddo

  write(filename,'("V-",A2,"-P",I0,"-",A3)') label, nst, symlabel
  open(unit=605,file=trim(adjustl(filename)),action='write',status='replace')
  write(605,*) nkmax, nch, symlabel
  write(605,'('//nstchar//'(I0,X))') (mst(j), j=1, nst)
  write(605,'('//nstchar//'(I0,X))') (parst(j), j=1, nst)
  write(605,'('//nstchar//'(I0,X))') (sst(j), j=1, nst)
  write(605,'('//nchchar//'(I0,X))') (st_ch(j),j=1,nch)
  write(605,'('//nchchar//'(I0,X))') (Lp_ch(j),j=1,nch)
  write(605,'('//nchchar//'(I0,X))') (Mp_ch(j),j=1,nch)
  write(605,'('//nstchar//'(ES15.8,X))') (en_est(j),j=1,nst)
  write(605,'('//nkchar//'(ES15.8,X))') (gridk(1:npk(j+1)-npk(j),j), j=1,nch)
  write(605,'('//nkchar//'(ES15.8,X))') (weightk(1:npk(j+1)-npk(j),j), j=1,nch)
  write(605,'('//nchchar//'(I0,X),I0)') npk(1:nch+1)
  do nkf=1,nkmax
    write(605,'('//nkchar//'("(",ES15.8,",",ES15.8,")",X))') (matrix(nkf,nki), nki=1,nkmax)
  enddo
  close(605)
end subroutine write_matrix_comp

end module optical_potential



!! BELOW IS THE PREVIOUS OPTICAL SUBROUTINE FOR POSTERITY

!! module optical_potential
!!   !This module stores the optical potential for use later in scat.f90
!!   !along with some info about number of states/channels in P space
!! 
!!   real*8, allocatable :: Wmat(:,:) !The optical potential
!!   complex*16, allocatable :: CWmat(:,:) !The optical potential
!!   integer :: nst_P    !number of states in P space (remaining states are in Q space)
!!   integer :: nchm_P   !number of channels in P space
!!   integer :: nchm_open_P
!!   integer :: nkgmax_P !number of k-grid points in P space
!!   logical :: optical
!!   logical :: weak_coupling
!!   logical :: complex_optical
!! contains
!! 
!! subroutine construct_optical_potential(Etot, Mtot, par_global, rspin_tot)
!!   !This subroutine constructs the optical potential 
!!   !
!!   !The weak-coupling optical potential is calculated directly:
!!   !
!!   ! <f|W|i> = ∑_m <f|V|m> w_m <m|V|i>
!!   !
!!   ! where f and i are states/channels/k points in P space
!!   ! and m is the same in Q space
!!   ! w_m are the Green's function & integration weights.
!!   ! 
!!   !
!!   !The full optical potential is calculated by solving the
!!   ! matrix equation:
!!   !       P        Q
!!   !    -------------------- -------     -------
!!   !    |     |            | |     |     |     |
!!   !  P |  I  |      -WV   | |  W  |     | WWC |
!!   !    |     |            | |     |     |     |
!!   !    -------------------- -------     -------
!!   !    |     |            | |     |     |     |
!!   !    |     |            | |     |  =  |     |
!!   !    |     |            | |     |     |     |
!!   !  Q |  0  |   I - WV   | |  X  |     | XWC |
!!   !    |     |            | |     |     |     |
!!   !    |     |            | |     |     |     |
!!   !    -------------------- -------     -------
!! 
!!   use channels
!!   use kgrid
!!   use input_data
!!   use MPI_module
!!  
!!   !Need the below for calculating V-matrix elements
!!   use contwaves
!!   use one_electron_func
!!   use sturmian_class
!!   use grid_radial
!!   use target_states
!!   
!!   implicit none
!!   
!!   !-INPUT DATA----------------------------------------------------------------!
!!   real*8, intent(in) :: Etot 
!!   integer, intent(in) :: Mtot
!!   real*8, intent(in) :: rspin_tot
!!   integer, intent(in) :: par_global
!! 
!!   !-LOOP INDICES--------------------------------------------------------------!
!!   integer :: nchf, nchi, nch
!!   integer :: nkf, nki, nqf, nqi, nq, kq, nkm, kqf, kqi
!!   integer :: nchloop, nchloopmax, numdone
!!   integer, allocatable, dimension(:) :: nchfloop, nchiloop, nchdone
!!   
!!   !-CALCULATING V MATRICES----------------------------------------------------!
!!   real*8, allocatable, dimension(:,:) :: Vmat, Vmatt, Vmatt_ex, VmatP
!!   complex*16, allocatable, dimension(:,:) :: CVmat !complex V matrix for when optical potential is complex
!!   integer :: Nmax_1e_orb, Nmax, nr, nqmf, nqmi, ifirst, Zasym, iborn, nstf, nsti
!!   real*8 :: Rd, Znuc, Zproj, theta
!!   real*8, dimension(grid%nr) :: grid_r, weight
!!   real*8, allocatable :: ortchil(:,:), flchil(:,:)
!!   logical :: spheroidal
!!   real*8 :: tolerance
!!   real*8 :: maxV
!!   
!!   !-CALCULATING W MATRIX------------------------------------------------------!
!!   real*8, allocatable, dimension(:,:) :: WmatWC !weak-coupling optical potential
!!   complex*16, allocatable, dimension(:,:) :: CWmatWC !complex weak-coupling optical potential
!!   real*8, dimension(nkgmax) :: gfw, Cgfw !Green's function and k integration weights (real, complex)
!!   real*8, allocatable, dimension(:,:) :: Xmat
!!   complex*16, allocatable, dimension(:,:) :: CXmat !for when W is complex
!!   integer :: nsize, io
!!   integer :: nchm_open_Q
!!   complex*16, allocatable, dimension(:,:) :: CMat
!! 
!!   !-REQUIRED BY DGESV SUBROUTINE----------------------------------------------!
!!   integer, allocatable, dimension(:) :: ipiv
!!   integer :: info
!!   
!!   real*8, parameter :: pi = acos(-1.0d0)
!!   real*8 :: t1, t2, temp, time
!!   complex*16 :: tempc
!!   integer :: hours, mins, secs
!!   integer :: percent, percent_prev
!!   real*8, external :: omp_get_wtime
!!   character*15 :: filename
!!   real*8, dimension(nchm,nchm) :: strengths
!!   real*8 :: spinf, spini, mem_vmat
!!   integer :: parstf, parsti, mstf, msti
!!   character*3 :: stflabel, stilabel
!! 
!!   !-RECONSTRUCTING LOCAL POTENTIAL--------------------------------------------!
!!   real*8, dimension(grid%nr,2) :: W_r, WC_r !second index angle = 0, pi/2
!!   real*8, dimension(grid%nr) :: temp_r, WCtemp_r
!!   type(sturmian_nr), pointer :: tnf, tni
!!   real*8, dimension(:), pointer :: wavef, wavei
!!   integer :: minf, maxf, Lf, Li, Mf, Mi, ir
!!   complex*16, external :: YLM
!!   
!!   !-SUBROUTINE START----------------------------------------------------------!
!!       
!!   spheroidal = data_in%calculation_type==2 .or. data_in%calculation_type==3
!!   tolerance = grid%expcut
!!   inquire(file='weak_coupling',exist=weak_coupling)
!! 
!!   nr = grid%nr
!!   grid_r = grid%gridr
!!   
!!   if(spheroidal) then
!!     !This can be changed when confirmed that spherical form works with bool weights
!!     weight = grid%bool
!!   else
!!     weight = grid%weight
!!   endif
!!   
!!   Rd = data_in%Rd
!!   Znuc = data_in%Z
!!   Zproj = data_in%Zproj
!!   Zasym = data_in%Zasym
!!   theta = data_in%theta
!!   ifirst = data_in%ifirst
!! 
!!   if(theta/=0.0d0) then
!!     error stop 'STOPPING in optical.f90: optical potential requires theta=0.0'
!!   elseif(ifirst==-1) then
!!     error stop 'STOPPING in optical.f90: cannot calculate optical potential in Born mode'
!!   endif
!! 
!!   !Determine number of channels and number of open channels in P space
!!   do nchm_P=1, nchm
!!     if(st_ch(nchm_P) > nst_P) exit
!!   enddo
!!   nchm_P = nchm_P - 1
!!   nchm_open_P = min(nchm_open,nchm_P)
!!   
!!   !If there are closed channels in P space then the optical potential will be complex
!!   complex_optical = nchm_open > nchm_P
!!   if(complex_optical) then
!!     print*, 'nchm_open > nchm_P: optical potential is complex'
!!   else
!!     print*, 'nchm_open <= nchm_P: optical potential is real'
!!   endif
!!   
!!   !set nkgmax_P to k-grid point at end of P space
!!   nkgmax_P = npk(nchm_P+1) - 1
!! 
!!   !Wmat or CWmat will be used depending on whether or not the optical potential is complex
!!   if(allocated(Wmat)) deallocate(Wmat)
!!   if(allocated(CWmat)) deallocate(CWmat)
!!   if(complex_optical) then
!!     allocate(CWmat(nkgmax_P,nkgmax_P),CWmatWC(nkgmax_P,nkgmax_P))
!!   else
!!     allocate(Wmat(nkgmax_P,nkgmax_P),WmatWC(nkgmax_P,nkgmax_P))
!!   endif
!! 
!!   !These will be required by exchange matrix subroutines
!!   Nmax_1e_orb = basis_size(bst_nr)
!!   if(spheroidal) then
!!     allocate(ortchil(Nmax_1e_orb,nkgmax), flchil(Nmax_1e_orb,nkgmax))
!!     call make_spheroidal_overlaps(Nmax_1e_orb, ortchil,flchil)
!!   else !spherical
!!     allocate(ortchil(Nmax_1e_orb,nkgmax),flchil(Nmax_1e_orb,nkgmax) )
!!     call make_overlaps_H12orb(ortchil,flchil,Nmax_1e_orb)
!!   endif
!! 
!!   !Remove small numbers to improve stability
!!   where(abs(ortchil)<data_in%expcut) ortchil = 0.0d0
!!   where(abs(flchil)<data_in%expcut) flchil = 0.0d0
!! 
!!   if(weak_coupling) then
!!     !weak-coupling optical potential only requires V-matrix elements between
!!     !P and Q space
!!     mem_vmat = dble(nkgmax) * dble(nkgmax_P) * 8.0d0 / (10.0d0**9.0d0)
!!     if(complex_optical) then
!!       mem_vmat = mem_vmat * 2
!!       print*,"V matrix memory(GB)", mem_vmat
!!       print*, 'nkgmax:', nkgmax
!!       allocate(CVmat(nkgmax_P,nkgmax))
!!     else
!!       print*,"V matrix memory(GB)", mem_vmat
!!       print*, 'nkgmax:', nkgmax
!!       allocate(Vmat(nkgmax_P,nkgmax))
!!     endif
!!   else
!!     !full optical potential requires full V matrix
!!     mem_vmat = dble(nkgmax) * dble(nkgmax) * 8.0d0 / (10.0d0**9.0d0)
!!     if(complex_optical) then
!!       mem_vmat = mem_vmat * 2
!!       print*,"V matrix memory(GB)", mem_vmat
!!       print*, 'nkgmax:', nkgmax
!!       allocate(CVmat(nkgmax,nkgmax))
!!       allocate(CXmat(nkgmax,nkgmax_P))
!!     else
!!       print*,"V matrix memory(GB)", mem_vmat
!!       print*, 'nkgmax:', nkgmax
!!       allocate(Vmat(nkgmax,nkgmax))
!!       allocate(Xmat(nkgmax,nkgmax_P))
!!     endif
!!   endif
!! 
!!   !Set arrays nchfloop and nchiloop to turn nested iteration over matrix into one loop for better OMP performance
!!   nchloopmax=nchm*nchm
!!   allocate(nchfloop(nchloopmax),nchiloop(nchloopmax),nchdone(nchloopmax))
!!   nchdone = 0
!!   nch = 0 
!!   do nchi=1, nchm
!!     do nchf=1, nchi
!!       if(weak_coupling .and. nchf > nchm_P) cycle
!!       nch = nch+ 1
!!       nchfloop(nch) = nchf; nchiloop(nch) = nchi
!!     enddo !nchf
!!   enddo !nchi
!!   nchloopmax = nch
!! 
!!   print*, 'MAKING VMATRIX FOR OPTICAL POTENTIAL'
!!   percent = 0; percent_prev = 0
!!   t1=omp_get_wtime()
!! 
!!   if(complex_optical) then
!!     CVmat = 0.0d0           !set entire V matrix to zero because the below loops may not cover
!!   else                      !all channel blocks
!!     Vmat = 0.0d0 
!!   endif
!! 
!!   !$OMP PARALLEL DO &
!!   !$OMP DEFAULT(PRIVATE) &
!!   !$OMP SHARED(nchfloop, nchiloop ,nchloop, nchloopmax, nchm_P, npk, nchm, Vmat, CVmat, rspin_tot, Lp_ch) &
!!   !$OMP SHARED(cont_oid, nr, grid_r, weight, Nmax_1e_orb, Nmax, ortchil, flchil, Etot, st_ch, st1_ch) &
!!   !$OMP SHARED(Rd, Znuc, Zproj, Zasym, theta, ifirst, chil, Mtot, ncwaves, spheroidal, weak_coupling, tolerance, nchdone, complex_optical)
!!   do nchloop=1,nchloopmax
!!     nchf=nchfloop(nchloop)
!!     nchi=nchiloop(nchloop)
!!       
!!     !Stuff for printing percent completed
!!     numdone = sum(nchdone)
!!     percent = (10*numdone)/nchloopmax * 10
!!     percent_prev = (10*(numdone-1))/nchloopmax * 10
!!     if(percent > percent_prev) write(*,'( "(",I3,"%)" )') percent
!! 
!!     nqmf = npk(nchf+1)-npk(nchf)
!!     nqmi = npk(nchi+1)-npk(nchi)
!!     nstf = st1_ch(nchf)
!!     nsti = st1_ch(nchi)
!!     allocate(Vmatt(nqmf,nqmi), Vmatt_ex(nqmf,nqmi))
!! 
!!     vmatt = 0.0d0
!!     vmatt_ex = 0.0d0
!! 
!!     !Calculate direct V-matrix elements
!!     if(spheroidal) then 
!!       call vdirect_2e_oid(rspin_tot, nchf,nchi, nqmf,nqmi, npk(nchf),npk(nchi), cont_oid, nr, grid_r, weight, Vmatt)
!!     else !spherical  
!!       call vdirect_2e(Mtot,rspin_tot, nchf,nchi, nstf,nsti, nqmf, nqmi, npk(nchf),npk(nchi) ,chil,nr,weight,grid_r,Zproj,Zasym,iBorn,vmatt)
!!     endif
!!     
!!     !Calculate exchange V-matrix elements
!!     if(ifirst == 1) then
!!       if(spheroidal) then
!!         call vexch_2e_oid(nchf,nchi, nqmf,nqmi, npk(nchf),npk(nchi), cont_oid,Nmax_1e_orb,ncwaves,ortchil,flchil, nr,grid_r,weight, Etot,Rd,Znuc,Zproj,rspin_tot,theta, Vmatt_ex)
!!       else
!!         call vexch_2e(Mtot,nchf,nchi,st_ch(nchf),st_ch(nchi),nqmf,nqmi,npk(nchf),npk(nchi),chil,Nmax_1e_orb,ncwaves,ortchil,flchil,nr,weight,grid_r,Zproj,Etot,theta,Znuc,Rd,rspin_tot,Vmatt_ex)
!!       endif !spheroidal
!!     endif !ifirst
!!         
!!     if(maxval(abs(Vmatt+vmatt_ex))==0.0d0) then
!!       print*, '>>> FULLY ZERO MATRIX: stf,sti, Lf,Li:', st_ch(nchf),st_ch(nchi),Lp_ch(nchf),Lp_ch(nchi)
!!     endif
!!   
!!     !Set small numbers to zero to reduce instability
!!     where(abs(Vmatt) < data_in%formcut) Vmatt = 0.0d0
!!     where(abs(Vmatt_ex) < data_in%formcut) Vmatt_ex = 0.0d0
!! 
!!     !Set exchange matrix elements to zero for first and last k
!!     vmatt_ex(2,:) = 0.0d0
!!     vmatt_ex(:,nqmi) = 0.0d0
!! 
!!     !Where exchange matrix element is insignificant compared to direct set to zero
!!     where(abs((vmatt-vmatt_ex)) < data_in%formcut * abs(vmatt)) vmatt_ex = 0.0d0
!! 
!!     !Copy matrix elements into main Vmat (or CVmat if complex) array
!!     do nkf=1,nqmf
!!       do nki=1,nqmi
!!         nqf = npk(nchf)+nkf-1
!!         nqi = npk(nchi)+nki-1
!!         if(complex_optical) then
!!           CVmat(nqf,nqi) = cmplx(Vmatt(nkf,nki)+Vmatt_ex(nkf,nki), 0.0d0)
!!           if(nqf/=nqi .and. .not.(weak_coupling.and.nchi>nchm_P)) CVmat(nqi,nqf) = CVmat(nqf,nqi)
!!         else
!!           Vmat(nqf,nqi) = Vmatt(nkf,nki) + Vmatt_ex(nkf,nki)
!!           if(nqf/=nqi .and. .not.(weak_coupling.and.nchi>nchm_P)) Vmat(nqi,nqf) = Vmat(nqf,nqi)
!!         endif
!!       enddo !nki
!!     enddo !nkf
!! 
!!     deallocate(Vmatt)
!!     deallocate(Vmatt_ex)
!! 
!!     nchdone(nchloop) = 1
!! 
!!   enddo !nchloop
!!   !$OMP END PARALLEL DO
!!   
!!   t2=omp_get_wtime()
!!   write(*,*)
!!   time = t2 - t1
!!   hours = time/3600
!!   mins = (time-3600*hours)/60
!!   secs = (time-3600*hours-60*mins)
!!   write(*,'("Time calculating V matrix: ",I2.2, " hours, ", I2.2, " mins, ", I2.2, " secs")') hours, mins, secs
!! 
!!   deallocate(ortchil)
!!   deallocate(flchil)
!!   deallocate(nchfloop)
!!   deallocate(nchiloop)
!!   deallocate(nchdone)
!! 
!!   !Now Vmat (or CVmat if complex) contains the V matrix and we continue on to calculate the optical potential
!! 
!!   !setup gfw array with Green's function & integration weights
!!   !integration is done in Q space only
!!   do nch = 1, nchm
!!     if(nch<=nchm_open) then
!!       gfw(npk(nch)) = -pi/gridk(1,nch) !Open channel 
!!     else
!!       gfw(npk(nch)) = 0.0d0
!!     endif
!!     do nq = npk(nch)+1, npk(nch+1) - 1
!!       kq = nq - npk(nch) + 1
!!       gfw(nq) = weightk(kq,nch) * grfn(kq,nch)
!!     end do
!!   end do
!! 
!!   print*, 'CALCULATING OPTICAL POTENTIAL:'
!!   t1=omp_get_wtime()
!! 
!!   !Calculate weak-coupling optical potential
!!   !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(temp,tempc,nki,nkf,nkm)
!!   do nki=1, nkgmax_P
!!     do nkf=1, nkgmax
!!       if(weak_coupling .and. nkf > nkgmax_P) cycle
!!       if(nkf <= nkgmax_P .and. nkf > nki) cycle
!!       temp = 0.0d0
!!       tempc = 0.0d0
!!       do nch=nchm_P+1, nchm
!!         nkm = npk(nch)
!!         if(nch<=nchm_open) then !also means complex_optical == .true.
!!           !on-shell subtraction if there are open channels in Q space
!!           !Factor of (0,1)=sqrt(-1) is here because we define Green's function as real
!!           tempc = tempc + CVmat(nkf,nkm)*gfw(nkm)*CVmat(nki,nkm)*(0.0d0,1.0d0)
!!         endif
!!         do nkm=npk(nch)+1, npk(nch+1)-1
!!           if(complex_optical) then
!!             tempc = tempc + CVmat(nkf,nkm)*gfw(nkm)*CVmat(nki,nkm)
!!           else
!!             temp = temp + Vmat(nkf,nkm)*gfw(nkm)*Vmat(nki,nkm)
!!           endif
!!         enddo !nkm
!!       enddo !nch
!!       
!!       if(complex_optical) then
!!         if(nkf<=nkgmax_P) CWmatWC(nkf,nki) = tempc
!!         if(nkf<=nkgmax_P) CWmatWC(nki,nkf) = tempc
!!       else
!!         if(nkf<=nkgmax_P) WmatWC(nkf,nki) = temp
!!         if(nkf<=nkgmax_P) WmatWC(nki,nkf) = temp
!!       endif
!!       
!!       if(.not.weak_coupling .and. complex_optical) then
!!         CXmat(nkf,nki) = tempc
!!         if(nkf<=nkgmax_P) CXmat(nki,nkf) = tempc
!!       elseif(.not.weak_coupling) then
!!         Xmat(nkf,nki) = temp
!!         if(nkf<=nkgmax_P) Xmat(nki,nkf) = temp
!!       endif
!!     enddo !nkf
!!   enddo !nki
!!   !$OMP END PARALLEL DO
!! 
!!   allocate(VmatP(1:nkgmax_P,1:nkgmax_P))
!!   do nki=1,nkgmax_P
!!     do nkf=1,nkgmax_P
!!       if(complex_optical) then
!!         VmatP(nkf,nki) = dble(CVmat(nkf,nki))
!!       else
!!         VmatP(nkf,nki) = Vmat(nkf,nki)
!!       endif
!!     enddo
!!   enddo
!! 
!!   !write V matrix to file
!!   if(complex_optical) then
!!     call write_matrix_comp(CVmat, nst_P, nchm_P, 'VV')
!!     call write_matrix_comp(CWmatWC, nst_P, nchm_P, 'WC')
!!   else
!!     call write_matrix(Vmat, nst_P, nchm_P, 'VV')
!!     call write_matrix(WmatWC, nst_P, nchm_P, 'WC')
!!   endif
!!   
!!   if(weak_coupling) then
!!     if(complex_optical) then
!!       CWmat = CWmatWC
!!       call write_matrix_comp(CWmatWC, nst_P, nchm_P, 'WW')
!!       call write_matrix_comp(CWmatWC+VmatP, nst_P, nchm_P, 'VW')
!!     else
!!       Wmat = WmatWC
!!       call write_matrix(WmatWC, nst_P, nchm_P, 'WW')
!!       call write_matrix(WmatWC+VmatP, nst_P, nchm_P, 'VW')
!!     endif
!!     return
!!   endif
!!   
!!   !prepare Vmat and Xmat for calculating full optical potential
!!   do nkf=1, nkgmax
!!     do nki=1, nkgmax
!!       if(complex_optical) then
!!         CVmat(nkf,nki) = CVmat(nkf,nki) * sqrt(abs(gfw(nkf)*gfw(nki)))
!!         if(nki<=nkgmax_P) CXmat(nkf,nki) = CXmat(nkf,nki) * sqrt(abs(gfw(nkf)))
!!       else
!!         Vmat(nkf,nki) = Vmat(nkf,nki) * sqrt(abs(gfw(nkf)*gfw(nki)))
!!         if(nki<=nkgmax_P) Xmat(nkf,nki) = Xmat(nkf,nki) * sqrt(abs(gfw(nkf)))
!!       endif
!!     enddo
!!   enddo
!!   
!!   if(complex_optical) then 
!!     
!!     !absorb factor of i=sqrt(-1) from green's function into Vmat
!!     do nch=nchm_P+1,nchm_open
!!       nki = npk(nch)
!!       do nkf=1,nkgmax 
!!         CVmat(nkf,nki) = cmplx(0.0d0,1.0d0)*CVmat(nkf,nki)
!!       enddo
!!     enddo
!!     
!!     CVmat = -CVmat
!!     do nkf=1, nkgmax
!!       do nki=1,nkgmax_P
!!         CVmat(nkf,nki) = 0.0d0
!!       enddo
!!       CVmat(nkf,nkf) = sign(1.0d0,gfw(nkf)) + CVmat(nkf,nkf)
!!     enddo
!!   
!!   else
!!     Vmat = -Vmat
!!     do nkf=1, nkgmax
!!       do nki=1,nkgmax_P
!!         Vmat(nkf,nki) = 0.0d0
!!       enddo
!!       Vmat(nkf,nkf) = sign(1.0d0,gfw(nkf)) + Vmat(nkf,nkf)
!!     enddo
!! 
!!   endif !complex optical
!! 
!!   nsize = nkgmax 
!! 
!!   allocate(ipiv(nsize))
!!   if(complex_optical) then
!!     call ZGESV(nsize,nkgmax_P,CVmat,nkgmax,ipiv,CXmat,nkgmax,info)
!!   else
!!     call DGESV(nsize,nkgmax_P,Vmat,nkgmax,ipiv,Xmat,nkgmax,info)
!!   endif
!!   deallocate(ipiv)
!!   
!!   if(complex_optical) then 
!!     deallocate(CVmat)
!!   else
!!     deallocate(Vmat)
!!   endif
!! 
!!   do nki=1,nkgmax_P 
!!     do nkf=1, nkgmax_P
!!       if(gfw(nkf)/=0.0d0) then
!!         if(complex_optical) then
!!           CWmat(nkf,nki) = CXmat(nkf,nki)/sqrt(abs(gfw(nkf)))*sign(1.0d0,gfw(nkf))
!!         else
!!           Wmat(nkf,nki) = Xmat(nkf,nki)/sqrt(abs(gfw(nkf)))*sign(1.0d0,gfw(nkf))
!!         endif
!!       else
!!         if(complex_optical) then
!!           CWmat(nkf,nki) = 0.0d0
!!         else
!!           Wmat(nkf,nki) = 0.0d0
!!         endif
!!       endif
!!     enddo
!!   enddo
!! 
!!   if(complex_optical) then
!!     deallocate(CXmat)
!!   else
!!     deallocate(Xmat)
!!   endif
!!   
!!   if(complex_optical) then
!!     call write_matrix_comp(CWmat, nst_P, nchm_P, 'WW')
!!     call write_matrix_comp(CWmat-CWmatWC, nst_P, nchm_P, 'WB')
!!     call write_matrix_comp(CWmat+VmatP, nst_P, nchm_P, 'VW')
!!     deallocate(CWmatWC)
!!   else
!!     call write_matrix(Wmat, nst_P, nchm_P, 'WW')
!!     call write_matrix(Wmat-WmatWC, nst_P, nchm_P, 'WB')
!!     call write_matrix(Wmat+VmatP, nst_P, nchm_P, 'VW')
!!     deallocate(WmatWC)
!!   endif
!!   deallocate(VmatP)
!! 
!!   t2=omp_get_wtime()
!!   time = t2 - t1
!!   hours = time/3600
!!   mins = (time-3600*hours)/60
!!   secs = (time-3600*hours-60*mins)
!!   write(*,'("Time calculating W matrix: ",I2.2, " hours, ", I2.2, " mins, ", I2.2, " secs")') hours, mins, secs
!! 
!! end subroutine construct_optical_potential
!!     
!! subroutine write_matrix(matrix, nst, nch, label)
!!   use channels
!!   use target_states
!!   use kgrid
!!   implicit none
!!   
!!   real*8, dimension(:,:), intent(in) :: matrix
!!   integer, intent(in) :: nst, nch
!!   character*2, intent(in) :: label
!!   character*20 :: filename
!!   integer, dimension(nst) :: sst, mst, parst
!!   integer :: nkmax, j, nkf, nki
!! 
!!   nkmax = npk(nch+1)-1
!! 
!!   do j=1, nst
!!     mst(j) = get_ang_mom_proj(TargetStates2el%b(j))
!!     parst(j) = get_par(TargetStates2el%b(j))
!!     sst(j) = 2*get_spin(TargetStates2el%b(j))+1
!!   enddo
!! 
!!   write(filename,'("V-",A2,"-P",I0,"-",A3)') label, nst, symlabel
!!   open(unit=605,file=trim(adjustl(filename)),action='write',status='replace')
!!   write(605,*) nkmax, nch, symlabel
!!   write(605,'('//nstchar//'(I0,X))') (mst(j), j=1, nst)
!!   write(605,'('//nstchar//'(I0,X))') (parst(j), j=1, nst)
!!   write(605,'('//nstchar//'(I0,X))') (sst(j), j=1, nst)
!!   write(605,'('//nchchar//'(I0,X))') (st_ch(j),j=1,nch)
!!   write(605,'('//nchchar//'(I0,X))') (Lp_ch(j),j=1,nch)
!!   write(605,'('//nchchar//'(I0,X))') (Mp_ch(j),j=1,nch)
!!   write(605,'('//nkchar//'(ES15.8,X))') (gridk(1:npk(j+1)-npk(j),j), j=1,nch)
!!   write(605,'('//nkchar//'(ES15.8,X))') (weightk(1:npk(j+1)-npk(j),j), j=1,nch)
!!   write(605,'(<nch+1>(I0,X))') npk(1:nch+1)
!!   do nkf=1,nkmax
!!     write(605,'('//nkchar//'("(",ES15.8,",",ES15.8,")",X))') (matrix(nkf,nki), 0.0d0, nki=1,nkmax)
!!   enddo
!!   close(605)
!! end subroutine write_matrix
!! 
!! subroutine write_matrix_comp(matrix, nst, nch, label)
!!   use channels
!!   use target_states
!!   use kgrid
!!   implicit none
!!   complex*16, dimension(:,:), intent(in) :: matrix
!!   integer, intent(in) :: nst, nch
!!   character*2, intent(in) :: label
!!   character*20 :: filename
!!   integer, dimension(nst) :: sst, mst, parst
!! 
!!   integer :: nkmax, j, nkf, nki
!! 
!!   nkmax = npk(nch+1)-1
!!   
!!   do j=1, nst
!!     mst(j) = get_ang_mom_proj(TargetStates2el%b(j))
!!     parst(j) = get_par(TargetStates2el%b(j))
!!     sst(j) = 2*get_spin(TargetStates2el%b(j))+1
!!   enddo
!! 
!!   write(filename,'("V-",A2,"-P",I0,"-",A3)') label, nst, symlabel
!!   open(unit=605,file=trim(adjustl(filename)),action='write',status='replace')
!!   write(605,*) nkmax, nch, symlabel
!!   write(605,'('//nstchar//'(I0,X))') (mst(j), j=1, nst)
!!   write(605,'('//nstchar//'(I0,X))') (parst(j), j=1, nst)
!!   write(605,'('//nstchar//'(I0,X))') (sst(j), j=1, nst)
!!   write(605,'('//nchchar//'(I0,X))') (st_ch(j),j=1,nch)
!!   write(605,'('//nchchar//'(I0,X))') (Lp_ch(j),j=1,nch)
!!   write(605,'('//nchchar//'(I0,X))') (Mp_ch(j),j=1,nch)
!!   write(605,'('//nkchar//'(ES15.8,X))') (gridk(1:npk(j+1)-npk(j),j), j=1,nch)
!!   write(605,'('//nkchar//'(ES15.8,X))') (weightk(1:npk(j+1)-npk(j),j), j=1,nch)
!!   write(605,'(<nch+1>(I0,X))') npk(1:nch+1)
!!   do nkf=1,nkmax
!!     write(605,'('//nkchar//'("(",ES15.8,",",ES15.8,")",X))') (matrix(nkf,nki), nki=1,nkmax)
!!   enddo
!!   close(605)
!! end subroutine write_matrix_comp
!! 
!! end module optical_potential
















