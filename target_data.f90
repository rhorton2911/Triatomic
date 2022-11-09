module target_data
  !
  ! Module contains target data in the structure "data_target" (which works similarly to data_in).
  ! Arrays are defined ONCE, early in input_data.f90 once the target name is read from input file.
  ! Their dimensions should resize to be as big as the calculation.
  ! Data for other targets can be added easily in a separate subroutine.
  ! 
  ! Modified 12/7/12 by J.S. from M.Z.'s original code.
  !
  use numbers
  public :: get_target_energy, get_target_oscstr
  !set_data_H2I


  type target_params
     integer :: nMax, lMax, mMax

     real(dp), dimension(:,:,:), pointer :: en_exact, sep_exact
     real(dp), dimension(:,:,:,:,:,:), pointer :: osc_exact
     integer, dimension(:,:,:), pointer :: en_ref, sep_ref
     integer, dimension(:,:,:,:,:,:), pointer :: osc_ref
  end type target_params


  type(target_params) :: data_target


contains

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine get_target_energy( self, n, l, mDeg, en_exact, en_ref )
  
  implicit none
  type(target_params), intent(in) :: self
  integer, intent(in) :: n, l, mDeg
  real(dp), intent(out) :: en_exact
  integer, intent(out) :: en_ref

  integer :: m


  m = abs(mDeg)

  if( n <= self%nMax .and. l <= self%lMax .and. m <= self%mMax ) then
     en_exact = self%en_exact(n,l,m)
     en_ref = self%en_ref(n,l,m) 
  else
     en_exact = 0d0
     en_ref = 0
  end if


end subroutine get_target_energy

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine get_target_oscstr( self, nL,lL,mLDeg, nR,lR,mRDeg, osc_exact, osc_ref )

  implicit none
  type(target_params), intent(in) :: self
  integer, intent(in):: nL,lL,mLDeg, nR,lR,mRDeg
  real(dp), intent(out):: osc_exact
  integer, intent(out):: osc_ref

  integer :: mL, mR


  mL = abs(mLDeg)
  mR = abs(mRDeg)

  if( nL <= self%nMax .and. lL <= self%lMax .and. mL <= self%mMax .and. nR &
     <= self%nMax .and. lR <= self%lMax .and. mR <= self%mMax ) then
     osc_exact = self%osc_exact( nL,lL,mL, nR,lR,mR )
     osc_ref = self%osc_ref( nL,lL,mL, nR,lR,mR )
  else
     osc_exact = 0d0
     osc_ref = 0
  end if


end subroutine get_target_oscstr

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine get_target_sepconst( self, n, l, mDeg, sep_exact, sep_ref )
  
  implicit none
  type(target_params), intent(in) :: self
  integer, intent(in) :: n, l, mDeg
  real(dp), intent(out) :: sep_exact
  integer, intent(out) :: sep_ref

  integer :: m


  m = abs(mDeg)

  if( n <= self%nMax .and. l <= self%lMax .and. m <= self%mMax ) then
     sep_exact = self%sep_exact(n,l,m)
     sep_ref = self%sep_ref(n,l,m) 
  else
     sep_exact = 0d0
     sep_ref = 0
  end if


end subroutine get_target_sepconst

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine set_data_H2plus( self, nMax, lMax, mMax )
  implicit none

  type(target_params) :: self
  integer, intent(in) :: nMax, lMax, mMax

  integer :: n, l, m


  n = max( nMax, 6 )   ! Ensure we have enough space
  l = max( lMax, 5 )   !  to cover both the exact data
  m = max( mMax, 2 )   !  as well as the requested states.

  self%nMax = n
  self%lMax = l
  self%mMax = m
  allocate( self%en_exact(n,0:l,-m:m), self%en_ref(n,0:l,-m:m) )
  allocate( self%osc_exact(n,0:l,-m:m, n,0:l,-m:m), self%osc_ref(n,0:l,-m:m, n,0:l,-m:m) )
  allocate( self%sep_exact(n,0:l,-m:m), self%sep_ref(n,0:l,-m:m) )

  
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
  ! Energies taken from papers, where  R = 2.0 a.u. 
  
!  self%en_exact(:,:,:) = 0.5  
  self%en_exact(:,:,:) = 0d0
  self%en_exact(1,0,0) = -6.026342144949D-1
  self%en_exact(2,0,0) =  1.391351246617D-1
  self%en_exact(3,0,0) =  3.223189549591D-1
  self%en_exact(2,1,0) = -1.675343922024D-1
  self%en_exact(3,1,0) =  2.445868349143D-1
  self%en_exact(4,1,0) =  3.626870757439D-1
  self%en_exact(3,2,0) =  2.642223711745D-1
  self%en_exact(4,2,0) =  3.692081223681D-1
  self%en_exact(4,3,0) =  3.733561298504D-1
  self%en_exact(5,3,0) =  4.191557039293D-1
  self%en_exact(5,4,0) =  4.196265315607D-1
  self%en_exact(6,5,0) =  4.443288012689D-1
  self%en_exact(2,1,1) =  7.122818010413D-2
  self%en_exact(3,2,1) =  2.733003733563D-1
  self%en_exact(3,1,1) =  2.991351700867D-1
  self%en_exact(4,2,1) =  3.732898693347D-1
  self%en_exact(4,3,1) =  3.738010794706D-1
  self%en_exact(5,4,1) =  4.196839242591D-1
  self%en_exact(3,2,2) =  2.872673181892D-1
  self%en_exact(4,3,2) =  3.750374570632D-1

  self%en_exact(:,:,:) = self%en_exact(:,:,:) - 0.5d0   ! Converts molecular energy to electronic.

  
  ! Refernces for exact Energies
  !  Ref[1]: Madsen and Peek 1971
  
  self%en_ref(:,:,:) = 0
  self%en_ref(1,0,0) = 1
  self%en_ref(2,0,0) = 1
  self%en_ref(3,0,0) = 1
  self%en_ref(2,1,0) = 1
  self%en_ref(3,1,0) = 1
  self%en_ref(4,1,0) = 1
  self%en_ref(3,2,0) = 1
  self%en_ref(4,2,0) = 1
  self%en_ref(4,3,0) = 1
  self%en_ref(5,3,0) = 1
  self%en_ref(5,4,0) = 1
  self%en_ref(6,5,0) = 1
  self%en_ref(2,1,1) = 1
  self%en_ref(3,2,1) = 1
  self%en_ref(3,1,1) = 1
  self%en_ref(4,2,1) = 1
  self%en_ref(4,3,1) = 1
  self%en_ref(5,4,1) = 1
  self%en_ref(3,2,2) = 1
  self%en_ref(4,3,2) = 1
  
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
  ! Oscillator Strengths n,l,m -> np,lp,mp stored as (n,l,m,np,lp,mp)
  
  self%osc_exact(:,:,:,:,:,:) = 0d0
  self%osc_exact(1,0,0,2,1,0) = 0.319    ! 1s(m=0)-2p(m=0)
  self%osc_exact(1,0,0,2,1,-1) = 0.46     ! 1s(m=0)-2p(m=-1)  
  self%osc_exact(1,0,0,2,1,1) = 0.46     ! 1s(m=0)-2p(m=1)
  self%osc_exact(1,0,0,3,1,0) = 8.24E-4  ! 1s(m=0)-3p(m=0) 
  self%osc_exact(2,1,0,2,0,0) = 1.36E-1  ! 2p(m=0)-2s(m=0)
  self%osc_exact(2,1,0,3,2,-1) = 0.268    ! 2p(m=0)-3d(m=-1)
  self%osc_exact(2,1,0,3,2,1) = 0.268    ! 2p(m=0)-3d(m=1)
  self%osc_exact(2,1,0,3,2,0) = 2.18E-1  ! 2p(m=0)-3d(m=0)  
  self%osc_exact(2,1,-1,3,2,-1) = 0.275    ! 2p(m=-1)-3d(m=-1)
  self%osc_exact(2,1,1,3,2,1) = 0.275    ! 2p(m=1)-3d(m=1)
  self%osc_exact(1,0,0,4,1,0) = 5.54E-5  ! 1s(m=0)-4p(m=0)             
  self%osc_exact(1,0,0,4,3,0) = 4.07E-5  ! 1s(m=0)-4f(m=0)  
  self%osc_exact(2,1,0,3,0,0) = 1.46E-2  ! 2p(m=0)-3s(m=0)
  
  ! Oscillator Strengths References
  ! Ref[0]: None, Ref[1]: D.R. Bates 1954, Ref[2]: D.R. Bates 1953, Ref[3]: D.R. Bates 1951
  
  
  self%osc_ref(:,:,:,:,:,:) = 0
  self%osc_ref(1,0,0,2,1,0) = 3      ! 1s(m=0)-2p(m=0)
  self%osc_ref(1,0,0,2,1,-1) = 2      ! 1s(m=0)-2p(m=-1)  
  self%osc_ref(1,0,0,2,1,1) = 2      ! 1s(m=0)-2p(m=1)
  self%osc_ref(1,0,0,3,1,0) = 1      ! 1s(m=0)-3p(m=0) 
  self%osc_ref(2,1,0,2,0,0) = 1      ! 2p(m=0)-2s(m=0)
  self%osc_ref(2,1,0,3,2,-1) = 2      ! 2p(m=0)-3d(m=-1)
  self%osc_ref(2,1,0,3,2,1) = 2      ! 2p(m=0)-3d(m=1)
  self%osc_ref(2,1,0,3,2,0) = 1      ! 2p(m=0)-3d(m=0)  
  self%osc_ref(2,1,-1,3,2,-1) = 2      ! 2p(m=-1)-3d(m=-1)
  self%osc_ref(2,1,1,3,2,1) = 2      ! 2p(m=1)-3d(m=1)
  self%osc_ref(1,0,0,4,1,0) = 1  ! 1s(m=0)-4p(m=0)             
  self%osc_ref(1,0,0,4,3,0) = 1  ! 1s(m=0)-4f(m=0)  
  self%osc_ref(2,1,0,3,0,0) = 1  ! 2p(m=0)-3s(m=0)

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
  ! Separation constants data.

  self%sep_exact(:,:,:) = 0d0
  self%sep_exact(1,0,0) = -1.393538844365D0
  self%sep_exact(2,0,0) = -4.732635792362D-1
  self%sep_exact(3,0,0) = -2.350163373006D-1
  self%sep_exact(2,1,0) = -2.521958176764D0
  self%sep_exact(3,1,0) = -2.202549502969D0
  self%sep_exact(4,1,0) = -2.109334455147D0
  self%sep_exact(3,2,0) = -6.226856306735D0
  self%sep_exact(4,2,0) = -6.125266676162D0
  self%sep_exact(4,3,0) = -1.212403991946D1
  self%sep_exact(5,3,0) = -1.207913362274D1
  self%sep_exact(5,4,0) = -2.007937550198D1
  self%sep_exact(6,5,0) = -3.005520938239D1
  self%sep_exact(2,1,1) = -2.682595167359D0
  self%sep_exact(3,2,1) = -6.258284857989D0
  self%sep_exact(3,1,1) = -2.320638001893D0
  self%sep_exact(4,2,1) = -6.144561697965D0
  self%sep_exact(4,3,1) = -1.213470103338D1
  self%sep_exact(5,4,1) = -2.008347603906D1
  self%sep_exact(3,2,2) = -6.364329962716D0
  self%sep_exact(4,3,2) = -1.216647631675D1

  ! Separation constants references.
  ! 0: None. 1: Madsen and Peek 1971
  self%sep_ref(:,:,:) = 0
  self%sep_ref(1,0,0) = 1
  self%sep_ref(2,0,0) = 1
  self%sep_ref(3,0,0) = 1
  self%sep_ref(2,1,0) = 1
  self%sep_ref(3,1,0) = 1
  self%sep_ref(4,1,0) = 1
  self%sep_ref(3,2,0) = 1
  self%sep_ref(4,2,0) = 1
  self%sep_ref(4,3,0) = 1
  self%sep_ref(5,3,0) = 1
  self%sep_ref(5,4,0) = 1
  self%sep_ref(6,5,0) = 1
  self%sep_ref(2,1,1) = 1
  self%sep_ref(3,2,1) = 1
  self%sep_ref(3,1,1) = 1
  self%sep_ref(4,2,1) = 1
  self%sep_ref(4,3,1) = 1
  self%sep_ref(5,4,1) = 1
  self%sep_ref(3,2,2) = 1
  self%sep_ref(4,3,2) = 1


end subroutine set_data_H2plus


! The below is yet to be properly implemented for H2.
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!
!!$subroutine set_data_H2(self, nMax,lMax,mMax)
!!$  implicit none
!!$
!!$  type(target_params) :: self
!!$  integer, intent(in) :: nMax, lMax, mMax
!!$
!!$  integer :: n, l, m
!!$
!!$
!!$  n = max( nMax, 0 )   ! Ensure we have enough space
!!$  l = max( lMax, 0 )   !  to cover both the exact data
!!$  m = max( mMax, 0 )   !  as well as the requested states.
!!$
!!$  self%nMax = n
!!$  self%lMax = l
!!$  self%mMax = m
!!$  allocate( self%en_exact(n,0:l,-m:m), self%en_ref(n,0:l,-m:m) )
!!$  allocate( self%osc_exact(n,0:l,-m:m, n,0:l,-m:m), self%osc_ref(n,0:l,-m:m, n,0:l,-m:m) )
!!$  allocate( self%sep_exact(n,0:l,-m:m), self%sep_ref(n,0:l,-m:m) )
!!$
!!$end subroutine set_data_H2



end module target_data
