subroutine add_phenom_pol(spin_mult, Mtot, par_global, nchm, vmat_in, chilphase, tdist, ifirstorder)
  !Author:          Liam Scarlett
  !Date created:    29-Aug-2018
  !Date modified:   30-Aug-2018
  !
  !This subroutine adds a phenomenological polarisation potential of the form
  !
  !   V_pol(r0) = -alpha0 * v(r0) - alpha2 * v(r0) * P2[cos(theta)],
  !
  !where
  !
  !   v(r0) = 1/(2*r0**4) * {1 - exp[-(r_0/a)^b]},
  !
  !P2 is the Legendre polynomial of order 2, and alpha0, alpha2 are the 
  !spherical and non-spherical polarisabilities, respectively.

  use MPI_module
  use kgrid
  use channels, only: Mp_ch, Lp_ch, st_ch, nent_ch, nchm_open
  use optical_potential, only: write_matrix
  implicit none

  !-INPUT-----------------------------------------------------------------------

    !Quantum numbers
    integer, intent(in) :: spin_mult, Mtot, par_global
    
    !number of channels
    integer, intent(in) :: nchm 
    
    !Input V-matrix - the polarisation potential matrix elements
    ! will be added to this
    real*8, dimension(nkgmax,nkgmax), intent(inout) :: vmat_in

    real*8, dimension(nkgmax,nkgmax) :: vmat, vmat_pol

    !distorting phase shifts
    real*8, dimension(nkgmax), intent(in) :: chilphase

    !Analytical term <k_f| U_0 |q_i> of physical T-matrix
    complex*16,dimension(nchm), intent(in):: tdist

    !ifirstorder = 1 if using Born approx, otherwise = 0
    integer, intent(in):: ifirstorder

    !If user requests, optimise cutoff radius to yield cross section 
    ! equal to ICS_target (read from file)
    real*8 :: ICS_target
  
  !-LOCAL-STUFF-----------------------------------------------------------------
    !on-shell T-matrix elements
    complex*16, dimension(nchm_open,nent_ch) :: ton, von, von_target
    complex*16, allocatable, dimension(:,:) :: ton_target
  
    !Need these things to call tmatsolver
    complex*16, dimension(nent_ch,nkgmax) :: TMat
    real*8,dimension(nchm_open,nchm_open) :: kon
    real*8,dimension(nkgmax-nchm,nchm_open) :: vdr
    real*8, dimension(nkgmax) :: gfw
    integer :: irec = 0

    !various loop indices
    integer:: nop,nkop, nch,nq,kq, nchf, nchi, nqf, nqch, nqchf, mkop, n,m, i,j
    integer :: nchm_potl, nent_potl, kf, ki
    real*8, parameter :: pi = acos(-1d0)
    real*8 :: tmp
    real*8 :: a1, a2, a1_start, a2_start, bbb, spher_pol, nonspher_pol, ICS1, ICS2, ICS_final, a_final
    real*8, allocatable :: aaa(:), ICS(:)
    integer :: nfile, n_iter, io, max_iter
    character*9 :: iterchar
    integer :: temp_S, temp_M, temp_par
    integer :: switch, n_steps, nchm_elastic
    character*3 :: sym
    character*20, parameter :: symchars = 'u gSPDFGHIKLMNOQRTUV'
    character*1 :: x
    real*8 :: da, da_temp
    logical :: found_root
    character*10 :: tempchar
    complex*16 :: tempcomp
    real*8 :: ICS_tol, a_tol
    real*8, allocatable :: thetaVec(:), DCS(:,:)
    real*8 :: scaling_coef
    logical :: scaling
    real*8 :: fixed_a
    logical :: potl_found
    real*8 :: ICS_temp

    nchm_elastic = 1
    do while(st_ch(nchm_elastic) == 1)
      nchm_elastic = nchm_elastic + 1
    enddo
    nchm_elastic = nchm_elastic - 1

    ICS_tol = 0.00001d0
    a_tol = ICS_tol !TODO: CHECK?
    ICS_target = -1.0d0

    open(newunit=nfile,file='pol.in',action='read')
    read(nfile,*) a1_start, a2_start, bbb, fixed_a
    read(nfile,*) spher_pol, nonspher_pol
    ICS_target = -1.0d0
    da = 1.0d0

    if(fixed_a > 0.0d0) then !we want to fix the cutoff radius and search for a scaling factor
      scaling = .true.
      if (a1_start == a2_start) stop 'STOP: need different a1_start, a2_start for scaling factor search.'
    else
      scaling_coef = 1.0d0
      scaling = .false.
    endif

    write(sym,'(I1,A1,A1)') spin_mult, symchars(Mtot+4:Mtot+4),symchars(par_global+2:par_global+2)
      
    print*, ''
    print*, '--------------------------------------------------'
    print*, '| Adding phenomenological polarisation potential |'
    print*, '--------------------------------------------------'
    print*, '  Exponent = ', bbb
    print*, '  Spherical polarisability = ', spher_pol
    print*, '  Nonspherical polarisability = ', nonspher_pol
    
    !determine search interval [a1_start,a2_start] and bisection width da for this symmetry
    do
      read(nfile,*,iostat=io) temp_M, temp_par, temp_S, a1, a2, da_temp, ICS_temp
      if(io /= 0) exit !reached end of file

      if(temp_M == Mtot .and. temp_par == par_global .and. temp_S == spin_mult  ) then !found correct symmetry
        a1_start = a1
        a2_start = a2
        da = da_temp
        ICS_target = ICS_temp
        exit
      endif
    enddo
    close(nfile)
    
    !If search a1 = a2 just add polarisation potential - don't optimise cutoff radius
    if(a1_start == a2_start) then 
      print*, '  Cutoff radius = ', a1_start
      call vpol_2e(nkgmax,npk,a1_start,bbb,spher_pol,nonspher_pol,vmat_pol)
      vmat_in = vmat_in + vmat_pol
      call write_matrix(Vmat_in, 1, nchm, 'VW')
      print*, 'WROTE MATRIX'
      return
    endif

    !read potl.in file to determine target ICS
!    potl_found = .false.
!    open(newunit=nfile,file='potl.in',action='read')
!    call read_potlhead(nfile) !read potl header
!
!    do
!      read(nfile,'("M,parity,spin,nchm_open,nent_ch: ",5I5)',err=6, end=5) temp_M, temp_par, temp_S, nchm_potl, nent_potl
!      write(*,'("M,parity,spin,nchm_open,nent_ch: ",5I5)') M, temp_par, temp_S, nchm_potl, nent_potl
!       
!      allocate(ton_target(nchm_potl,nent_potl))
!      do i=1,nent_potl
!        do n=1,nchm_potl
!          read(nfile,'("(",E15.8,",",E15.8,") ")',err=6,end=6,ADVANCE='No') tempcomp
!          if(i<=nchm .and. n <= nchm) then
!            von(n,i) = tempcomp
!          endif
!        enddo
!        read(nfile,'("")')
!      enddo
!      do i=1,nent_potl
!        do n=1,nchm_potl
!          read(nfile,'("(",E15.8,",",E15.8,") ")',err=6,end=6,ADVANCE='No') ton_target(n,i)
!        enddo
!        read(nfile,'("")')
!      enddo
!
!      if(temp_M == Mtot .and. temp_par == par_global .and. nint(2.0d0*(temp_S+0.5d0)+1.0d0) == spin_mult) then
!        print*, 'TON_TARGET FOUND'
!        potl_found = .true.
!        exit
!      endif
!
!      deallocate(ton_target)
!    enddo
!    !now ton_target and von_target contain T and V matrix elements for the
!    !present symmetry read from the potl.in file
!        
!    5 continue
!    close(nfile)
!    if(.not.potl_found) stop 'STOPPING in phenom_pol.f90 : could not find right symmetry in potl.in'
!
!    !set target ICS
!    ICS_target = 4.0d0*pi**3 * sum(abs(ton_target)**2)

    if(ICS_target<0.0d0) error stop 'STOPPING in phenom_pol.f90 : did not set ICS_target'


    if ((a1_start < 0.0d0 .or. a2_start < a1_start) .and. .not. scaling) error stop ' ERROR: invalid search interval for cutoff radius'
    
    if(.not.scaling) then
      print*, '  Optimising cutoff radius to yield ICS', ICS_target
      print*, ''
      print*, '  ----------------------------------------------'
      print*, '    Iteration | Cutoff radius    | ICS'
      print*, '    ------------------------------------------'
    else
      print*, '  Cutoff radius = ', fixed_a
      print*, '  Optimising scaling factor to yield ICS', ICS_target
      print*, ''
      print*, '  ----------------------------------------------'
      print*, '    Iteration | Scaling factor   | ICS'
      print*, '    ------------------------------------------'
    endif

    !setup gfw array for calls to tmatsolver
    do nch = 1, nchm
      do nq = npk(nch), npk(nch+1) - 1
        kq = nq - npk(nch) + 1
        gfw(nq) = weightk(kq,nch) * grfn(kq,nch)
      end do
    end do
   
    !first we start at a2_start and step down in increments of da until a root
    !is found, then the bisection method narows down on the root
    max_iter = ceiling((log(a2_start-a1_start)-log(1.0D-10))/log(2.0d0)) !max iterations from bisection method
    n_steps = ceiling((a2_start-a1_start)/da) + 1 !max steps in original search
    max_iter = max_iter + n_steps
    allocate(aaa(max_iter), ICS(max_iter))

    a1=a2_start+da
    found_root = .false.
    
    print*, 'MAX_ITER:', max_iter
    print*, 'a2_start, a1_start:', a2_start, a1_start
    print*, 'da:', da
    print*, 'size(aaa):', size(aaa)

    do n_iter=1, n_steps
      write(iterchar,'(I9.0)') n_iter
      iterchar = trim(adjustl(iterchar))//'         '
      a1 = a1 - da
      if (a1 < 0.0d0 .and. .not. scaling) a1 = 0.0d0
      
      if(scaling) then
        call vpol_2e(nkgmax,npk,fixed_a,bbb,spher_pol,nonspher_pol,vmat_pol)
        vmat = vmat_in + a1*vmat_pol
      else
        call vpol_2e(nkgmax,npk,a1,bbb,spher_pol,nonspher_pol,vmat_pol)
        vmat = vmat_in + vmat_pol
      endif
      
      call tmatsolver(spin_mult,nqmax,nchm,st_ch,npk,nkgmax,vmat,gfw,chilphase,tdist,gridk,nchm,nchm,ifirstorder,ton,irec,TMat,.false.)
      
      PRINT*, 'NCHM_ELASTIC:', NCHM_ELASTIC
      PRINT*, 'SIZE(TON):', SIZE(TON)
      ICS1 = 4.0d0 * pi**3 * sum(abs(ton(1:nchm_elastic,1:nchm_elastic))**2)


      print*, 'n_iter:', n_iter
      print*, 'iterchar:', iterchar
      print*, 'size(aaa), size(ICS):', size(aaa), size(ICS)
      print*, 'aaa(n_iter):', aaa(n_iter)
      print*, 'ICS(n_iter):', ICS(n_iter)
      aaa(n_iter) = a1
      ICS(n_iter) = ICS1
      write(*,'(5X, A9, 2X, F16.13, 3X, ES12.5, 3X, ES12.5)', advance='no') iterchar, aaa(n_iter), ICS(n_iter)

      if(n_iter > 1 .and. sign(1.0d0,ICS1-ICS_target) /= sign(1.0d0,ICS2-ICS_target) .and. .not. found_root) then
        write(*,'("  --  Found root")')
        a1_start = a1
        a2_start = a1+da
        found_root = .true.
        exit 
      else
        write(*,*)
      endif
      ICS2 = ICS1
    enddo

    !if the above found a root then a1_start and a2_start have been modified to
    !be one da step on either side of the root and the bisection search will
    !begin with those endpoints.
    !If not, then the bisection search will begin with the original search
    !interval
      
    a1 = a1_start
    a2 = a2_start
      
    if(scaling) then
      call vpol_2e(nkgmax,npk,fixed_a,bbb,spher_pol, nonspher_pol, vmat_pol)
      vmat = vmat_in + a1*vmat_pol
    else
      call vpol_2e(nkgmax,npk,a1,bbb,spher_pol, nonspher_pol, vmat_pol)
      vmat = vmat_in + vmat_pol
    endif
    call tmatsolver(spin_mult,nqmax,nchm,st_ch,npk,nkgmax,vmat,gfw,chilphase,tdist,gridk,nchm,nchm,ifirstorder,ton,irec,TMat,.false.)
    ICS1 = 4.0d0 * pi**3 * sum(abs(ton(1:nchm_elastic,1:nchm_elastic))**2)
      
    if(scaling) then
      call vpol_2e(nkgmax,npk,fixed_a,bbb,spher_pol, nonspher_pol, vmat_pol)
      vmat = vmat_in + a2*vmat_pol
    else
      call vpol_2e(nkgmax,npk,a2,bbb,spher_pol, nonspher_pol, vmat_pol)
      vmat = vmat_in + vmat_pol
    endif
    call tmatsolver(spin_mult,nqmax,nchm,st_ch,npk,nkgmax,vmat,gfw,chilphase,tdist,gridk,nchm,nchm,ifirstorder,ton,irec,TMat,.false.)
    ICS2 = 4.0d0 * pi**3 * sum(abs(ton(1:nchm_elastic,1:nchm_elastic))**2)
      
    do 
      n_iter = n_iter + 1
      write(iterchar,'(I9.0)') n_iter
      iterchar = trim(adjustl(iterchar))//'         '
     
      aaa(n_iter) = (a1+a2)/2.0d0
        
      if(scaling) then
        call vpol_2e(nkgmax,npk,fixed_a,bbb,spher_pol, nonspher_pol, vmat_pol)
        vmat = vmat_in + aaa(n_iter)*vmat_pol
      else
        call vpol_2e(nkgmax,npk,aaa(n_iter),bbb,spher_pol, nonspher_pol, vmat_pol)
        vmat = vmat_in + vmat_pol
      endif
      call tmatsolver(spin_mult,nqmax,nchm,st_ch,npk,nkgmax,vmat,gfw,chilphase,tdist,gridk,nchm,nchm,ifirstorder,ton,irec,TMat,.false.)
      ICS(n_iter) = 4.0d0 * pi**3 * sum(abs(ton(1:nchm_elastic,1:nchm_elastic))**2)
    
      write(*,'(5X, A9, 2X, F16.13, 3X, ES12.5)') iterchar, aaa(n_iter), ICS(n_iter)
      
      print*, 'ICS(n_iter):', ICS(n_iter)
      print*, 'ICS_target:', ICS_target
      print*, 'ICS_tol:', ICS_tol
      print*, 'a2:', a2
      print*, 'a1:', a1
      print*, 'a_tol:', a_tol
      print*, 'abs(ICS(n_iter)-ICS_target)/ICS_target:', abs(ICS(n_iter)-ICS_target)/ICS_target
      
      if(abs(ICS(n_iter)-ICS_target)/ICS_target <= ICS_tol .or. (a2-a1)<= a_tol) exit
      if(n_iter>3) then
        if(abs(ICS(n_iter)-ICS(n_iter-1)) <= ICS_tol .and. abs(ICS(n_iter-1)-ICS(n_iter-2)) <= ICS_tol) exit
      endif
       
      if(sign(1.0d0,ICS1-ICS_target) == sign(1.0d0,ICS2-ICS_target)) then 
        !(if endpoints of the current interval are both on the same side of the target ICS)

        if(abs(ICS1-ICS_target) .gt. abs(ICS2-ICS_target)) then
          !if the right endpoint is closer to the target ICS, move the left endpoint in
          a1 = aaa(n_iter)
          ICS1 = ICS(n_iter)
        else
          !if the left endpoint is closer to the target ICS, move the right endpoint in
          a2 = aaa(n_iter)
          ICS2 = ICS(n_iter)
        endif
      !else: endpoints are on different sides of the target ICS
      elseif(sign(1.0d0,ICS(n_iter)-ICS_target) == sign(1.0d0,ICS1-ICS_target)) then
        !if middle point is on the same side as the left endpoint, move the
        !left endpoint in
        a1 = aaa(n_iter)
        ICS1 = ICS(n_iter)
      else
        !otherwise, move the right endpoint in
        a2 = aaa(n_iter)
        ICS2 = ICS(n_iter)
      endif
    enddo
    print*, '  ----------------------------------------------'
        
    ICS_final = ICS(n_iter)
    a_final = aaa(n_iter)

    call sort_two_arrays(aaa(1:n_iter), ICS(1:n_iter), n_iter)
    open(newunit=nfile,file=sym//'.log',action='write',status='replace')
    write(nfile,'("Cut-off radius     ",A3," pw ICS", 5X, "TARGET ICS")') sym
    do i=1, n_iter
      write(nfile,'(F16.13,4(3X,ES11.5))') aaa(i), ICS(i), ICS_target
    enddo
    close(nfile)

    !check for situations where ICS oscillates around ICS_target
    switch = 0
    if(abs(ICS1-ICS_target) < abs(ICS_final-ICS_target)) then
      ICS_final = ICS1
      a_final = a1
      switch = 1
    endif
    if(abs(ICS2-ICS_target) < abs(ICS_final-ICS_target)) then
      ICS_final = ICS2
      a_final = a2
      switch = 1
    endif

    if(abs(ICS1-ICS_target)/ICS_target > 0.1d0) then
      if(scaling) then
        print*, '  Could not find scaling factor in given search interval'
      else
        print*, '  Could not find cutoff radius in given search interval'
      endif
      error stop 'stopping...'
    endif
  
    if(scaling) then
      write(*,'("  OPTIMISED SCALING FACTOR FOR (2S+1,M,par)=(",I0,",",I0,",",I0,"): ",F16.13)') spin_mult, Mtot, par_global, a_final
    else
      write(*,'("  OPTIMISED CUTOFF RADIUS FOR (2S+1,M,par)=(",I0,",",I0,",",I0,"): ",F16.13)') spin_mult, Mtot, par_global, a_final
    endif
    print*, ''
  
    if(switch == 1) then
      if(scaling) then
        call vpol_2e(nkgmax,npk,fixed_a,bbb,spher_pol, nonspher_pol, vmat_pol)
      else
        call vpol_2e(nkgmax,npk,a_final,bbb,spher_pol, nonspher_pol, vmat_pol)
      endif
    endif
    if(scaling) then    
      vmat_in = vmat_in + a_final*vmat_pol
    else
      vmat_in = vmat_in + vmat_pol
    endif

    call write_matrix(Vmat_in, 1, nchm, 'VW')

    return

    6 error stop ' ERROR: problem reading potl.in file'
end subroutine add_phenom_pol

subroutine eps_from_ton(ton, nchm, eps)
  use kgrid
  implicit none
  integer, intent(in) :: nchm
  complex*16, dimension(nchm,nchm), intent(in) :: ton
  real*8, intent(out) :: eps

  complex*16, dimension(nchm,nchm) :: tmpA, tmpB
  real*8, dimension(nchm,nchm) :: kon, z, ID
  real*8, dimension(nchm) :: eigenphases
  integer :: ierr, i
  integer, parameter :: matz = 0
  real*8, parameter :: pi = acos(-1.0d0)

  eps = 0.0d0
  return

  tmpA = - (0.0d0,1.0d0)*pi*ton * gridk(1,1)
  do i=1,nchm
    tmpA(i,i) = 1.0d0 + tmpA(i,i)
  enddo
  tmpB = ton * gridk(1,1) **2
  call MATINV_comp(tmpA,nchm,nchm,tmpB)
  kon = dble(tmpB)

  ID = 0.0d0
  do i=1,nchm
    ID(i,i) = 1.0d0
  enddo
  !call rsg(nchm,nchm,kon,ID,eigenphases,matz,z,ierr)
  eigenphases = atan(eigenphases)
  do i=1,nchm
    if(eigenphases(i)<0.0d0) eigenphases(i) = eigenphases(i) + pi
  enddo
  eps = sum(eigenphases)
end subroutine eps_from_ton

subroutine eps_from_kon(kon, nchm, eps)
  use kgrid
  implicit none
  integer, intent(in) :: nchm
  real*8, dimension(nchm,nchm), intent(in) :: kon
  real*8, intent(out) :: eps

  real*8, dimension(nchm,nchm) :: ID, tmpA, z
  real*8, dimension(nchm) :: eigenphases
  integer :: ierr, i
  integer, parameter :: matz = 0
  real*8, parameter :: pi = acos(-1.0d0)
  
  eps = 0.0d0
  return

  tmpA = kon * gridk(1,1) **2

  ID = 0.0d0
  do i=1,nchm
    ID(i,i) = 1.0d0
  enddo
  !call rsg(nchm,nchm,tmpA,ID,eigenphases,matz,z,ierr)
  eigenphases = atan(eigenphases)
  do i=1,nchm
    if(eigenphases(i)<0.0d0) eigenphases(i) = eigenphases(i) + pi
  enddo
  eps = sum(eigenphases)
end subroutine eps_from_kon

subroutine vpol_2e_oid(nkgmax,npk,aaa,bbb,spher_pol,nonspher_pol,vmatt)
  use contwaves
  use channels
  use spheroidal_class
  use input_data
  use grid_radial
  use kgrid, only: gridk
  
  use one_electron_func
  use ovlpste1me
  use sturmian_class
  use target_states

  implicit none
  
  integer, intent(in) :: nkgmax
  integer, intent(in):: npk(nchm+1) ! On-shell k-grid point for particular channel
  real*8, intent(in) :: aaa, bbb, spher_pol, nonspher_pol
  real*8, dimension(nkgmax,nkgmax), intent(inout):: vmatt

  type(spheroidal_fn), pointer :: wavef, wavei, wavefp, waveip
  real*8, pointer, dimension(:) :: radvecf, radveci
  integer :: lamf, lami, m, lf, li, jlf, jli, j, kf, ki, neta, lmin, lmax, mmin, mmax, i1eta, i2eta, i1min, i2max
  integer :: nchf, nchi, i1, i2, irho, ieta, nchfp, nchip, mf
  real*8 :: tmp, angcoeff
  real*8, allocatable :: angintegrand(:), ttt(:), ttti(:), angint(:,:,:,:) !angint(lf,li,m,rho)
  real*8, allocatable :: etaweights(:), eta(:)
  real*8 :: R, rho
  integer :: lambdaf, lambdai
  real*8, external :: PLM !associated Legendre polynomial
  real*8, external :: PL
  real*8, external :: vfunc_oid, Yint
  real*8 :: vmat_converted(nkgmax,nkgmax)
  character*20 :: filename
  real*8 :: deta, costheta
  real*8,allocatable :: temperoni(:)
  real*8, parameter :: pi = acos(-1.0d0)
  real*8 :: t1, t2, lmin_temp
  real*8, external :: omp_get_wtime
  real*8 :: weight(grid%nr), rhogrid(grid%nr)
  character*5 :: aaachar
  neta = 101
  deta = 1.0d0/dble(neta-1)

  R = data_in%Rd

  weight = grid%weight
  rhogrid = grid%gridr

  if(mod(2,neta)==0) neta = neta + 1
  allocate(angintegrand(neta), eta(neta), etaweights(neta), ttt(grid%nr), ttti(grid%nr))
  do j=1,neta
    eta(j) = 0.0d0 + dble(j-1)*deta
    etaweights(j) = dble(mod(j+1,2)*2+2)
  enddo
  etaweights(1) = 1.0d0
  etaweights(neta) = 1.0d0
  etaweights = etaweights * deta / 3.0d0

  !find the largest number of l terms, and min/max bounds for rho integration
  lmax = 0
  mmax = 0
  i1min = grid%nr
  i2max = 0
  do kf=1, nkgmax
    wavef => get_spheroidal_fn(cont_oid, kf)

    do j=get_num_terms(wavef), 1, -1
      if(get_term_D(wavef,j)>grid%expcut) exit
    enddo
    j=get_num_terms(wavef)
    lf = get_term_l(wavef, j)
    mf = get_m(wavef)
    if (lf > lmax) lmax = lf
    if (mf > mmax) mmax = mf
    j = get_rad_min(wavef)
    if ( j < i1min ) i1min = j
    j = get_rad_max(wavef)
    if ( j > i2max ) i2max = j
  enddo
  lmin = lmax
  mmin = mmax
  do kf=1, nkgmax
    wavef => get_spheroidal_fn(cont_oid, kf)
    lf = get_term_l(wavef, 1)
    mf = get_m(wavef)
    if (lf < lmin) lmin = lf
    if (mf < mmin) mmin = mf
  enddo

  t1 = omp_get_wtime()
  if(allocated(angint)) deallocate(angint)
  allocate(angint(lmin:lmax,lmin:lmax,mmin:mmax,i2max))
  angint = 0.0d0
  angintegrand = 0.0d0
  !calculate angular integral for each (lf,li,m,rho)
  !$OMP PARALLEL DO &
  !$OMP DEFAULT(private) &
  !$OMP SHARED(angint, deta, R, aaa, bbb, etaweights, eta, spher_pol, nonspher_pol) &
  !$OMP SHARED(i1min, i2max, lmin, lmax, mmin, mmax, nchm, neta, rhogrid, Mp_ch) &
  !$OMP SCHEDULE(guided) 
  !!$OMP COLLAPSE(2)
  do lf=lmin, lmax, 2
    do li=lf, lmax, 2
      !conservation of parity:
      if(mod(li,2) /= mod(lf,2)) cycle
      !conservation of l when not including quadrupole polarisability:
      !if(nonspher_pol == 0.0d0 .and. lf /= li) cycle
      !selection rule for lf,li from CGC when including quadrupole polarisability:
      !if(li /=lf .and. li < abs(lf-2) .or. li > lf + 2) cycle
      !print*, 'lf,li,mpch:', lf, li, Mp_ch(nchm)
      do m=minval(abs(Mp_ch(1:nchm))), min(min(lf,li), abs(Mp_ch(nchm)))
        tmp = 0.0d0
        !write(filename,'("ang-",I0,".",I0,".",I0)') lf, li, m
        !open(unit=111,file=filename,action='write')

        !-- ANALYTICAL EVALUATION OF ANGULAR INTEGRAL ---------------------

          !add spherical polarisability term, includes delta(lf,li)
          if(lf == li) tmp = spher_pol

          !add quadrupole polarisability term
          if(li >= abs(lf-2) .and. li <= lf + 2 .and. nonspher_pol > 0.0d0) then
            tmp = tmp + nonspher_pol * Yint(dble(lf),dble(m),2.0d0,0.0d0,dble(li),dble(m))
          endif

          angint(lf,li,m,i1min:i2max) = -1.0d0 * tmp * (1.0d0-exp(-(rhogrid(i1min:i2max)/aaa)**bbb)) / rhogrid(i1min:i2max)**2 / 2.0d0
        !------------------------------------------------------------------

        !-- NUMERICAL EVALUATION OF ANGULAR INTEGRAL ----------------------
              
          !evaluate coefficients arising from writing spherical harmonics in
          !terms of Legendre polynomials
          angcoeff = 1.0d0
          do j=li-m+1,li+m
            angcoeff = angcoeff * dble(j)
          enddo
          do j=lf-m+1,lf+m
            angcoeff = angcoeff * dble(j)
          enddo
          angcoeff = sqrt(dble(2*lf+1)*dble(2*li+1)/angcoeff)
            
          do irho = i1min, i2max
            rho = rhogrid(irho)
            do ieta = 1, neta

              !evaluate spherical cos(theta_0) for the P_2(cos theta) in quadrupole term
              costheta = (2.0d0*rho/R+1.0d0)*eta(ieta) / sqrt((2.0d0*rho/R+1.0d0)**2 + eta(ieta)**2-1.0d0)
              if(ieta == neta) costheta = 1.0d0

              angintegrand(ieta) = PLM(lf,m,eta(ieta))*vfunc_oid(rho,eta(ieta),aaa,bbb)*PLM(li,m,eta(ieta)) &
                & * (spher_pol + nonspher_pol * PL(2,costheta)) * ((rho+R/2.0d0)**2 + (eta(ieta)*R/2.0d0)**2)
            enddo !eta

            !integral:
            call minmaxi(angintegrand, neta, i1eta, i2eta)
            tmp = angcoeff * sum(angintegrand(i1eta:i2eta) * etaweights(i1eta:i2eta))

            !if numerical integral has converged to the analytical form then
            !exit the rho loop and keep the analytical form for the remaining
            !(larger) rho's
            if(R/rho <= 0.1d0 .and. abs(angint(lf,li,m,irho)-tmp) <= 1.0D-3) exit 

            !otherwise, replace the analytical form with the numerical result
            angint(lf,li,m,irho) = tmp

          enddo !rho
        !--------------------------------------------------------------------

        if(li /= lf) angint(li,lf,m,:) = angint(lf,li,m,:)
      enddo !mi
    enddo !li
  enddo !lf
  !$OMP END PARALLEL DO
  t2 = omp_get_wtime()
  !print*, ':::::::: time spent on angular integrals:', t2-t1

  !V-matrix calculation
  vmatt = 0.0d0
  do nchf = 1, nchm
    wavef => get_spheroidal_fn(cont_oid, npk(nchf))
    lambdaf = get_lambda(wavef)
    m = get_m(wavef)
    do nchi = 1, nchm
      wavei => get_spheroidal_fn(cont_oid, npk(nchi))
      lambdai = get_lambda(wavei)

      if (mod(lambdaf,2) /= mod(lambdai,2)) cycle !all matrix elements conserve pseudo parity

      if (m /= get_m(wavei)) cycle !all matrix elements conserve m

      if(st_ch(nchf) /= st_ch(nchi)) cycle !matrix elements only between equal states

      t1 = omp_get_wtime()
      !$OMP PARALLEL DO &
      !$OMP DEFAULT(private) &
      !$OMP SHARED(vmatt, npk, nchf, nchi, m, cont_oid, weight, angint, lmin, lmax, grid) &
      !$OMP SCHEDULE(guided)
      !!$OMP COLLAPSE(2)
      kfloop: do kf = npk(nchf), npk(nchf+1)-1
        wavef => get_spheroidal_fn(cont_oid, kf)
        kiloop: do ki =  max(kf,npk(nchi)), npk(nchi+1)-1
          !if(kf > ki) cycle !V matrix is symmetric, only calculate half
          wavei => get_spheroidal_fn(cont_oid, ki)

          i1 = max(get_rad_min(wavef), get_rad_min(wavei))
          i2 = min(get_rad_max(wavef), get_rad_max(wavei))
          if(i1>=i2) cycle

          ttt = 0.0d0
          !summation over spherical angular momentum in expansion of spheroidal harmonics
          do jlf = 1, get_num_terms(wavef)
            lf = get_term_l(wavef,jlf)
            if(lf>lmax) cycle
            if(lf<lmin) cycle
            ttti = 0.0d0
            do jli = 1, get_num_terms(wavei)
              li = get_term_l(wavei,jli)
              if(li>lmax) cycle
              if(li<lmin) cycle
        
              if(get_term_D(wavei,jli) < grid%expcut) cycle
              ttti(i1:i2) = ttti(i1:i2) + get_term_D(wavei,jli) * angint(lf,li,abs(m),i1:i2)
            enddo !li
            
            if(get_term_D(wavef,jlf) < grid%expcut) cycle
            ttt(i1:i2) = ttt(i1:i2) + get_term_D(wavef,jlf) * ttti(i1:i2)
          enddo !lf

          radvecf => get_rad_vec(wavef)

          radveci => get_rad_vec(wavei)

          ttt(i1:i2) = (2.0d0/pi) * radvecf(i1:i2) * ttt(i1:i2) * radveci(i1:i2)* weight(i1:i2)

          call minmaxi(ttt, i2, i1,i2)

          !print*, 'kf,ki:', kf,ki
          vmatt(kf,ki) = sum(ttt(i1:i2))
        
        enddo kiloop
      enddo kfloop
      !$OMP END PARALLEL DO
      t2 = omp_get_wtime()

      !copy lower triangle of V-matrix to upper triangle
      do kf=1, nkgmax
        do ki=1, kf-1
          vmatt(kf,ki) = vmatt(ki,kf)
        enddo
      enddo

    enddo 
  enddo 

return

  vmat_converted = 0.0d0
  !V_{lfmf,limi} = sum_{lambdaf,lambdai} ...
  
  do nchf=1, nchm !loop over lfmf
    !get nchf wave with on-shell k just to set lf and mf
    wavef => get_spheroidal_fn(cont_oid, npk(nchf))
    m = get_m(wavef)
    lf = get_lambda(wavef)
    do nchi=1, nchm !loop over limi
      wavei => get_spheroidal_fn(cont_oid,npk(nchi))
      !m still conserved here
      if(get_m(wavei) /= m) cycle
      !l conserved in spherical matrix elements (without nonspher pol)
      li = get_lambda(wavei)
      if(lf /= li) cycle

      !now we calculate spherical V_{lfmf,limi} from
      !spheroidal matrix elements for each kf,ki
      do kf=npk(nchf), npk(nchf+1)-1
        do ki=npk(nchi), npk(nchi+1)-1

          !now sum over lambdaf,lambdai.
          !sum over all channels and pick ones with correct m
          ! - equivalent to summing over lambda
          do nchfp=1,nchm
            wavef => get_spheroidal_fn(cont_oid,npk(nchfp))
            if(get_m(wavef) /= m) cycle
            lambdaf = get_lambda(wavef)
            do nchip=1,nchm
              wavei => get_spheroidal_fn(cont_oid,npk(nchip))
              if(get_m(wavei) /= m) cycle
              lambdai = get_lambda(wavei)

              !here we want to add the lambdaf, lambdai spheroidal V matrix
              !and associated D terms to the spherical V matrix element.
              !This corresponds to the nchf,nchi term of the spheroidal
              !V matrix, now we need to find the lf,li terms of the 
              !lambdaf and lambdai D terms:
              do jlf=1,get_num_terms(wavef)
                if(get_term_l(wavef,jlf) /= lf) cycle
                do jli=1,get_num_terms(wavei)
                  if(get_term_l(wavei,jli) /= li) cycle
                  vmat_converted(kf,ki) = vmat_converted(kf,ki) &
                    & + (0.0d0,1.0d0)**(lambdai-li+lf-lambdaf) &
                    & * get_term_D(wavef,jlf) * get_term_D(wavei,jli) &
                    & * vmatt(kf,ki)
                  if(kf==1 .and. ki==1) then
                    print*, 'ADDING TERM TO lf,li=', lf, li
                    print*, 'lambdaf, lambdai=', lambdaf, lambdai
                    print*, '  term_Df, termDi = ', get_term_D(wavef,jlf),get_term_D(wavei, jli)
                    print*, '  num terms = ', get_num_terms(wavef),get_num_terms(wavei)
                    print*, ''
                  endif
                enddo !jli
              enddo !jlf
            enddo !nchip
          enddo !nchfp
        enddo !ki
      enddo !kf
      
      write(filename,'("vmat-",I0,".",I0)') nchf, nchi
      open(unit=111,file=trim(adjustl(filename))//'converted', action='write', status='replace')
      do kf = npk(nchf), npk(nchf+1)-1
        write(111,'(10000(E12.5,3X))') gridk(kf-npk(nchf)+1,nchf), vmat_converted(kf,npk(nchi):npk(nchi+1)-1)
      enddo
      close(111)
    
    enddo !nchi
  enddo !nchf

end subroutine vpol_2e_oid

function vfunc_oid(rho, eta, aaa, bbb) result(V)
  use input_data
  implicit none
  real*8, intent(in) :: rho, eta, aaa, bbb
  real*8 :: V
  real*8 :: r

  r = data_in%Rd/2.0d0 * sqrt((2.0d0*rho/data_in%Rd + 1.0d0)**2 + eta**2 -1)

  V = -1.0d0/(2.0d0*r**4) * (1.0d0 - exp(-(r/aaa)**bbb))

end function vfunc_oid
          
      
subroutine vpol_2e(nkgmax,npk,aaa,bbb,spher_pol, nonspher_pol, vmatt)
  use input_data
  use channels

  implicit none
  integer, intent(in) :: nkgmax
  integer, intent(in):: npk(nchm+1) 
  real*8, intent(in) :: aaa, bbb, spher_pol, nonspher_pol
  real*8, dimension(nkgmax,nkgmax), intent(inout):: vmatt

  if(data_in%calculation_type == 0 .or. data_in%calculation_type == 1) then
    !spherical
    call vpol_2e_spherical(nkgmax,npk,aaa,bbb,spher_pol, nonspher_pol, vmatt)
  else
    !spheroidal
    call vpol_2e_oid(nkgmax,npk,aaa,bbb,spher_pol, nonspher_pol, vmatt)
  endif
end subroutine vpol_2e

subroutine vpol_2e_spherical(nkgmax,npk,aaa,bbb,spher_pol, nonspher_pol, vmatt)
  use sturmian_class
  use target_states
  use one_electron_func
  use channels
  use contwaves
  use input_data
  use distpot
  use vnc_module
  use grid_radial
  use kgrid, only: gridk

  implicit none
  integer, intent(in) :: nkgmax
  integer, intent(in):: npk(nchm+1) ! On-shell k-grid point for particular channel
  real*8, intent(in) :: aaa, bbb, spher_pol, nonspher_pol
  real*8, dimension(nkgmax,nkgmax), intent(inout):: vmatt

  integer:: ki, kf, i, ir
  integer :: nchf, nchi
  real*8, dimension(grid%nr):: ttt, temp, fun0, fun1
  type(sturmian_nr), pointer:: tn0,tnp0,tn1, tn2, tnp1, tnp2     ! One-electron orbitals
  real*8, pointer, dimension(:)::  fp0, f0, fp1, f1, fp2, f2
  ! Radial integrations
  integer::  minf, maxf, minfp, maxfp,  ir1, ir2, irt1,irt2   
  integer::  ir2_nuc, ir2_e_e, ir_all
  real*8:: sum_rad
  real*8, parameter :: pi = acos(-1d0)
  character*20 :: filename
  
  real*8:: Yint
  integer:: Lap0, La0, Mp0, M0! Projectile Quantum Numbers (QN)  
  
  vmatt = 0.0d0

  do nchf = 1, nchm
    do nchi = 1, nchm
      ! Projectile quantum numbers
      Lap0 = Lp_ch(nchf)  
      Mp0 = Mp_ch(nchf)
      La0 = Lp_ch(nchi)  
      M0 = Mp_ch(nchi)

      if(Mp0 /= M0) cycle

      ttt = 0.0d0
      if(Lap0 == La0) then
        !NIC: modify these two radial parts
        ttt = ttt - spher_pol/(2.0d0*grid%gridr**4)*(1.0d0-exp(-(grid%gridr/aaa)**bbb))
      endif
      if(La0>=abs(Lap0-2) .and. La0<=Lap0+2 .and. nonspher_pol > 0.0d0) then
        ttt = ttt - nonspher_pol/(2.0d0*grid%gridr**4)*(1.0d0-exp(-(grid%gridr/aaa)**bbb)) &
          & * Yint(dble(Lap0),dble(Mp0),2.0d0,0.0d0,dble(La0),dble(M0))
      endif
      
      call minmaxi(ttt,grid%nr,ir1,ir2)
      
      if (ir1 >= ir2) cycle

      do ki=npk(nchi), npk(nchi+1)-1
        tn0 => chil%b(ki)
        f0 => fpointer(tn0)
        minf = get_minf(tn0)
        maxf = get_maxf(tn0)
    
        ! $OMP PARALLEL DO SHARED(chil, vmatt, minf, maxf, ki, npk, f0) PRIVATE(tnp0, fp0, minfp, maxfp, ir2, ir1, fun0, kf)
        do kf=npk(nchf), npk(nchf+1)-1
            tnp0 => chil%b(kf)
            fp0 => fpointer(tnp0)
            minfp = get_minf(tnp0)
            maxfp = get_maxf(tnp0)
        
            ir2 = min(maxf,maxfp,ir2) 
            ir1 = max(minf,minfp,ir1)

            fun0(ir1:ir2) = (f0(ir1:ir2) * fp0(ir1:ir2)) * grid%weight(ir1:ir2) * ttt(ir1:ir2)  

            vmatt(kf,ki) = sum(fun0(ir1:ir2)) * (2d0/pi)
        end do !kf
        ! $OMP END PARALLEL DO
      end do !ki
      write(filename,'("vmat-",I0,".",I0)') nchf, nchi
      open(unit=111,file=trim(adjustl(filename)), action='write', status='replace')
      do kf = npk(nchf), npk(nchf+1)-1
        write(111,'(10000(E12.5,3X))') gridk(kf-npk(nchf)+1,nchf),vmatt(kf,npk(nchi):npk(nchi+1)-2)
      enddo
      close(111)


    enddo !nchi
  enddo !nchf

end subroutine vpol_2e_spherical

subroutine sort_two_arrays(X, Y, N)
  implicit none
  integer, intent(in) :: N
  real*8, intent(inout) :: X(0:N-1), Y(0:N-1)
  integer :: i, j
  real*8 :: a, b

  do i=1, N-1
    a=X(i)
    b=Y(i)
    j=i-1
    do while (j >= 0)
      if(X(j) <= a) exit
      X(j+1) = X(j)
      Y(j+1) = Y(j)
      j=j-1
    enddo !j
    X(j+1)=a
    Y(j+1)=b
  enddo !N
end subroutine sort_two_arrays
  
subroutine sort_three_arrays(X, Y, Z, N)
    implicit none
    real*8, intent(inout) :: X(0:N-1), Y(0:N-1), Z(0:N-1)
    integer, intent(in) :: N
    integer :: i, j
    real*8 :: a, b, c

    do i=1, N-1
      a=X(i)
      b=Y(i)
      c=Z(i)
      j=i-1
      do while (j >= 0)
        if(X(j) <= a) exit
        X(j+1) = X(j)
        Y(j+1) = Y(j)
        Z(j+1) = Z(j)
        j=j-1
      enddo !j
      X(j+1)=a
      Y(j+1)=b
      Z(j+1)=c
    enddo !N
  
end subroutine sort_three_arrays

