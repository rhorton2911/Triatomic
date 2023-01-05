module spheroidal_class

  public :: get_index, get_lambda, get_m, get_momentum, get_energy, &
       get_rad_min, get_rad_max, get_rad_val, get_rad_vec, &
       get_num_terms, get_term_l, get_term_D, get_l_vec, get_D_vec, &
       create_spheroidal_fn, destroy_spheroidal_fn, get_spheroidal_fn, &
       copy_spheroidal_fn, convert_from_sturm_to_oid, &
       create_spheroidal_basis, destroy_spheroidal_basis, get_num_fns, &
       oid_overlap, angular_coeffs_2, angular_coeffs_3


  type, public :: spheroidal_fn
     !private   !LHS removed this to make it easier to broadcast continuum waves across MPI tasks
     integer :: j   ! An index or principle quantum number for this function.
     real*8 :: E   ! Energy.
     integer :: lambda   ! The "pseudo-angular-momentum" or spheroidal l.
     integer :: m   ! The angular momentum projection of this function.
     integer :: nTerms   ! The size of the spherical harmonic expansion.
     integer :: radMin,radMax   ! Minimum & maximum indices of radFn.
     real*8, pointer, dimension(:) :: radFn   ! The (separated) radial part.
     real*8, pointer, dimension(:) :: DVec   ! Coefficients for the Ylm's.
     integer, pointer, dimension(:) :: lVec   ! The l's of these Ylm's.

  end type spheroidal_fn


  type, public :: spheroidal_basis
     type(spheroidal_fn), pointer, dimension(:) :: b   ! Vector of functions.
     integer :: numFns   ! Number of spheroidal functions.
  end type spheroidal_basis


  interface get_energy
     module procedure get_energy
  end interface


contains !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

  function get_index(self)
    implicit none
    integer :: get_index
    type(spheroidal_fn), intent(in) :: self

    get_index = self%j
  end function get_index

  function get_lambda(self)
    implicit none
    integer :: get_lambda
    type(spheroidal_fn), intent(in) :: self

    get_lambda = self%lambda
  end function get_lambda

  function get_m(self)
    implicit none
    integer :: get_m
    type(spheroidal_fn), intent(in) :: self

    get_m = self%m
  end function get_m

  subroutine set_m( self, mIn )
    implicit none
    type(spheroidal_fn), intent(inout) :: self
    integer, intent(in) :: mIn

    self%m = mIn
  end subroutine set_m

  function get_num_terms(self)
    implicit none
    integer :: get_num_terms
    type(spheroidal_fn), intent(in) :: self

    get_num_terms = self%nTerms
  end function get_num_terms

  function get_rad_min( self )
    implicit none
    integer :: get_rad_min
    type(spheroidal_fn), intent(in) :: self

    get_rad_min = self%radMin
  end function get_rad_min

  function get_rad_max( self )
    implicit none
    integer :: get_rad_max
    type(spheroidal_fn), intent(in) :: self

    get_rad_max = self%radMax
  end function get_rad_max

  function get_momentum( self )
    implicit none
    real*8 :: get_momentum
    type(spheroidal_fn), intent(in) :: self

    if ( self%E < 0d0 ) then
       stop 'negative energy!'
    else
       get_momentum = dsqrt(2d0*self%E)
    end if
  end function get_momentum

  function get_energy( self )
    implicit none
    real*8 :: get_energy
    type(spheroidal_fn), intent(in) :: self

    get_energy = self%E
  end function get_energy

  subroutine set_energy( self, EIn )
    implicit none
    type(spheroidal_fn), intent(inout) :: self
    real*8, intent(in) :: EIn

    self%E = EIn
  end subroutine set_energy

  function get_num_fns( self )
    implicit none
    integer :: get_num_fns
    type(spheroidal_basis), intent(in) :: self

    get_num_fns = self%numFns
  end function get_num_fns

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

  function get_rad_val( self, i )
    implicit none
    real*8 :: get_rad_val
    type(spheroidal_fn), intent(in) :: self
    integer, intent(in) :: i
    integer :: i1,i2

    i1 = self%radMin ; i2 = self%radMax
    if ( i == 0 ) then
       get_rad_val = self%radFn(i1)
    elseif ( i1<=i .and. i<=i2 ) then
       get_rad_val = self%radFn(i)
    else
       get_rad_val = 0d0
    end if
  end function get_rad_val

  function get_rad_vec( self )
    implicit none
    real*8, pointer, dimension(:) :: get_rad_vec
    type(spheroidal_fn), intent(in) :: self

    get_rad_vec => self%radFn(:)
  end function get_rad_vec

  function get_l_vec( self )
    implicit none
    integer, pointer, dimension(:) :: get_l_vec
    type(spheroidal_fn), intent(in) :: self

    get_l_vec => self%lVec(:)
  end function get_l_vec

  function get_D_vec( self )
    implicit none
    real*8, pointer, dimension(:) :: get_D_vec
    type(spheroidal_fn), intent(in) :: self

    get_D_vec => self%DVec(:)
  end function get_D_vec

  function get_term_l( self, termNum )
    implicit none
    integer :: get_term_l
    type(spheroidal_fn), intent(in) :: self
    integer, intent(in) :: termNum
    
    get_term_l = self%lVec(termNum)
  end function get_term_l  

  function get_term_D( self, termNum )
    implicit none
    real*8 :: get_term_D
    type(spheroidal_fn), intent(in) :: self
    integer, intent(in) :: termNum

    get_term_D = self%DVec(termNum)
  end function get_term_D

  function get_spheroidal_fn( self, fnNum )
    implicit none
    type(spheroidal_fn), pointer :: get_spheroidal_fn
    type(spheroidal_basis), intent(in) :: self
    integer, intent(in) :: fnNum

    get_spheroidal_fn => self%b(fnNum)
  end function get_spheroidal_fn

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

  subroutine init_spheroidal_fn(self)
    implicit none
    type(spheroidal_fn), intent(inout) :: self

    self%j = 0
    self%E = 0d0
    self%lambda = -1
    self%m = -1
    self%radMin = 0
    self%radMax = 0
    self%nTerms = 0
  end subroutine init_spheroidal_fn
 
  subroutine create_spheroidal_fn(self, j,E,lambda,m, radMin,radMax,gridMax,radFn, nTerms,lVec,DVec)
    implicit none

    type(spheroidal_fn), intent(inout) :: self
    real*8, intent(in) :: E
    integer, intent(in) :: j,lambda,m, radMin,radMax,gridMax, nTerms
    real*8, dimension(1:gridMax), intent(in) :: radFn
    real*8, dimension(nTerms), intent(in) :: DVec
    integer, dimension(nTerms), intent(in) :: lVec

    self%j = j
    self%E = E
    self%lambda = lambda
    self%m = m
    self%radMin = radMin
    self%radMax = radMax
    self%nTerms = nTerms

    allocate( self%radFn(1:radMax), self%lVec(1:nTerms), self%DVec(1:nTerms) )
!    self%radFn(:) = 0d0
    self%radFn(radMin:radMax) = radFn(radMin:radMax)
    self%lVec(1:nTerms) = lVec(1:nTerms)
    self%DVec(1:nTerms) = DVec(1:nTerms)

  end subroutine create_spheroidal_fn

  subroutine copy_spheroidal_fn( newFn, oldFn )
    implicit none

    type(spheroidal_fn), intent(in) :: oldFn
    type(spheroidal_fn), intent(out) :: newFn

    integer, dimension(oldFn%nTerms) :: lVec
    real*8, dimension(oldFn%nTerms) :: DVec

    lVec(:) = oldFn%lVec(:)
    DVec(:) = oldFn%DVec(:)
    call create_spheroidal_fn( newFn, oldFn%j,oldFn%E,oldFn%lambda,oldFn%m, oldFn%radMin,oldFn%radMax,oldFn%radMax,oldFn%radFn,oldFn%nTerms,lVec(:),DVec(:) )

  end subroutine copy_spheroidal_fn

  subroutine destroy_spheroidal_fn( self )
    implicit none
    type(spheroidal_fn), intent(inout) :: self

    deallocate( self%radFn, self%lVec,self%DVec )
    call init_spheroidal_fn(self)

  end subroutine destroy_spheroidal_fn

  subroutine convert_from_sturm_to_oid(sturm, oid)
    use sturmian_class
    implicit none
    type(sturmian_nr), intent(in) :: sturm
    type(spheroidal_fn), intent(out) :: oid

    integer :: j, lambda,m, radMin,radMax, nTerms
    real*8 :: E
    real*8, dimension(:), pointer :: sturmVec
    real*8, dimension(:), allocatable :: radVec
    integer, dimension(1) :: lVec
    real*8, dimension(1) :: DVec

    j = get_k(sturm)
    E = get_energy(sturm)
    lambda = get_ang_mom(sturm); m = get_ang_mom_proj(sturm)
    nTerms = 1; lVec(1) = lambda; DVec(1) = 1d0

    radMin = get_minf(sturm); radMax = get_maxf(sturm)
    allocate( radVec(1:radMax) )
    sturmVec => fpointer(sturm)

!    call create_spheroidal_fn(oid, j,E,lambda,m, radMin,radMax,radMax,radVec, nTerms,lVec,DVec)
    oid%j = j
    oid%E = E
    oid%lambda = lambda
    oid%m = m
    oid%radMin = radMin
    oid%radMax = radMax
    oid%nTerms = nTerms
    
    !allocate( oid%radFn(1:radMax), oid%lVec(1:nTerms), oid%DVec(1:nTerms) )
    allocate( oid%lVec(1:nTerms), oid%DVec(1:nTerms) )
    oid%radFn(radMin:radMax) => fpointer(sturm) !sturmVec(radMin:radMax)
    oid%lVec(1:nTerms) = lVec(1:nTerms)
    oid%DVec(1:nTerms) = DVec(1:nTerms)

  end subroutine convert_from_sturm_to_oid


  subroutine create_spheroidal_basis(self, numFnsIn)
    implicit none
    type(spheroidal_basis), intent(inout) :: self
    integer, intent(in) :: numFnsIn
    integer :: fnNo

    self%numFns = numFnsIn
    if (numFnsIn > 0) then
       allocate( self%b(numFnsIn) )
       do fnNo = 1, numFnsIn
          call init_spheroidal_fn( self%b(fnNo) )
       end do
    end if

  end subroutine create_spheroidal_basis

  subroutine destroy_spheroidal_basis( self )
    implicit none
    type(spheroidal_basis), intent(inout) :: self
    integer :: fnNo

    do fnNo = 1, self%numFns
       call destroy_spheroidal_fn( self%b(fnNo) )
    end do
    self%numFns = 0
    if ( associated(self%b) ) deallocate( self%b )

  end subroutine destroy_spheroidal_basis

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

  function oid_overlap(fnp,fn, i1,i2,inVec)
    !
    use grid_radial
    use input_data
    implicit none
    
    type(spheroidal_fn), intent(in) :: fnp,fn
    integer, intent(in), optional :: i1,i2
    real*8, dimension(1:grid%nr), optional :: inVec
    real*8 :: oid_overlap   ! Output.
    
    integer :: j1,j2
    real*8, dimension(:), pointer :: fpVec,fVec
    real*8, dimension(1:grid%nr) :: jVec
    
    
    oid_overlap = 0d0
    
    j1 = 1; j2 = grid%nr

!OLD METHODS
!    jVec(:) = grid%gridr(:)   ! Pass in as rho, comes out as angular quadratic.
    call angular_quadratic_2(fnp,fn, j1,j2,jVec)
    if (j1 >= j2) return   ! Angular quadratic is zero => overlap is zero.

!NEW METHODS
!!$    ! Do the angular and Jacobian parts here numerically.
!!$    mp = dble(fnp%m); m = dble(fn%m)
!!$    do jp = 1, fnp%nTerms
!!$       lp = dble( fnp%lVec(jp) ); Dp = fnp%DVec(jp)
!!$       do j = 1, fn%nTerms
!!$          l = dble( fn%lVec(j) ); D = fn%DVec(j)
!!$          call angular_coeffs_2(lp,mp, l,m, rhoTmp,etaTmp)
!!$          rhoCoeff = rhoCoeff + Dp*D*rhoTmp
!!$          etaCoeff = etaCoeff + Dp*D*etaTmp
!!$       end do
!!$    end do
!!$    if (rhoCoeff==0d0 .and. etaCoeff==0d0) return   ! Overlap=0 due to ang.mom.
!!$
!!$    ! Make jVec = rhoCoeff*rho(rho+R) + etaCoeff*R^2.
!!$    jVec(:) = rhoCoeff*grid%gridr(:)*(grid%gridr(:)+data_in%Rd) + etaCoeff*R*R

    ! If overlapping with a radial function (e.g. distpot) insert it here.
    if (present(inVec)) then
       j1 = i1; j2 = i2
       jVec(j1:j2) = jVec(j1:j2) * inVec(j1:j2)
    end if

    ! Do the radial integration numerically.
    fpVec => get_rad_vec(fnp); fVec => get_rad_vec(fn)
    j1 = max(j1,get_rad_min(fnp)); j1 = max(j1,get_rad_min(fn))
    j2 = min(j2,get_rad_max(fnp)); j2 = min(j2,get_rad_max(fn))
    oid_overlap = sum(fpVec(j1:j2)*fVec(j1:j2) * jVec(j1:j2)*grid%weight(j1:j2))


  end function oid_overlap

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

  subroutine angular_quadratic_2(fnp,fn, j1,j2,jVec)
    !
    use input_data
    use grid_radial
    implicit none

    type(spheroidal_fn), intent(in) :: fnp,fn
    integer, intent(inout) :: j1,j2
    real*8, dimension(j1:j2), intent(out) :: jVec
    
    integer :: j2tmp, indp,ind
    real*8 :: R, lp,mp, l,m, rhoTerm,etaTerm, Dp,DpD, tmp,Yint
!    real*8, dimension(j1:j2) :: rhoVec
    real*8, dimension(:), pointer :: rhoVec


    j2tmp = j2; j2 = j1   ! This flags that the angular quadratic is zero.
    mp = dble(fnp%m); m = dble(fn%m)
    if (mp /= m) return   ! Conservation of orbital angular momentum projection.

    R = data_in%Rd
!    rhoVec(:) = jVec(:); jVec(:) = 0d0
    rhoVec => grid%gridr

    rhoTerm = 0d0; etaTerm = 0d0
    do indp = 1, fnp%nTerms
       lp = dble(fnp%lVec(indp)); Dp = fnp%DVec(indp)
       do ind = 1, fn%nTerms
          l = dble(fn%lVec(ind)); DpD = Dp*fn%DVec(ind)
          if (lp == l) rhoTerm = rhoTerm + DpD

          tmp = Yint(lp,mp, 2d0,0d0, l,m)
          if (tmp /= 0d0) etaTerm = etaTerm + DpD*tmp
       end do
    end do
    if (rhoTerm==0d0 .and. etaTerm==0d0) return
    
    j2 = j2tmp   ! Angular quadratic is non-zero.
    jVec(:) = (rhoVec(:)*(rhoVec(:)+R)+R*R/6d0)*rhoTerm - R*R/6d0*etaTerm
    if (data_in%calculation_type == 3) jVec(:) = jVec(:) /(rhoVec(:)+R)/(rhoVec(:)+R)


  end subroutine angular_quadratic_2

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

  subroutine angular_coeffs_2(lp,mp, l,m, rhoCoeff,etaCoeff)
    !
    implicit none

    real*8, intent(in) :: lp,mp, l,m
    real*8, intent(out) :: rhoCoeff,etaCoeff

    real*8 :: Yint

    rhoCoeff = 0d0; etaCoeff = 0d0
    if (mp /= m) return   ! Conservation or orbital angular momentum projection.

    if (lp == l) then
       rhoCoeff = 1d0
       etaCoeff = 1d0/6d0
    end if
    
    etaCoeff = etaCoeff - Yint(lp,mp, 2d0,0d0, l,m)/6d0

  end subroutine angular_coeffs_2


  subroutine angular_coeffs_3(lp,mp, lam,mu, l,m, rhoCoeff,etaCoeff)
    !
    implicit none

    real*8, intent(in) :: lp,mp, lam,mu, l,m
    real*8, intent(out) :: rhoCoeff,etaCoeff

    real*8 :: tmp, Yint, F3Y
    

    rhoCoeff = 0d0; etaCoeff = 0d0
    if ((-1)**nint(lp+lam+l) /= 1) return   ! Conservation of parity.
    if (mp /= mu + m) return   ! Conservation of angular momentum projection.

    tmp = Yint(lp,mp, lam-2d0,mu, l,m)
    if (tmp /= 0d0) etaCoeff = -1d0/4d0 * dsqrt((lam-1d0)**2-mu**2)*dsqrt(lam**2-mu**2) /(2d0*lam-1d0)/(2d0*lam+1d0) * tmp

    tmp = Yint(lp,mp, lam,mu, l,m)
    if (tmp /= 0d0) then
       rhoCoeff = tmp
       etaCoeff = etaCoeff + 1d0/2d0 * (lam*(lam+1d0)+mu*mu-1d0) /(2d0*lam-1d0)/(2d0*lam+3d0) * tmp
    end if

    tmp = Yint(lp,mp, lam+2d0,mu, l,m)
    if (tmp /= 0d0) etaCoeff = etaCoeff - 1d0/4d0 * dsqrt((lam+1d0)**2-mu**2)*dsqrt((lam+2d0)**2-mu**2) /(2d0*lam+1d0)/(2d0*lam+3d0) * tmp
    

  end subroutine angular_coeffs_3


end module spheroidal_class

