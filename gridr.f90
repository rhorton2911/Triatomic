module grid_radial
 use numbers
 implicit none
 public

! This is declaration of the type rgrid.
 type rgrid
    integer:: npwave                      ! the number of points per oscillation
    integer:: npdbl                       ! the number of points with the same dx per interval
    real(dpf):: rmax                         ! maximum value of r
    real(dpf):: qmax                         ! maximum value of q for which integration is accurate
    integer:: ndouble                     ! number of doubling
    real(dpf):: formcut
    real(dpf):: regcut
    real(dpf):: expcut
! the above parameters are used to set up r-grid decribed below
    real(dpf), dimension(:), pointer:: gridr   ! r-grid values
    real(dpf), dimension(:), pointer:: weight  ! r-grid weights for Simpson rule
    real(dpf), dimension(:), pointer:: dr      ! dr(i) = grid(i) - grid(i-1)
    real(dpf), dimension(:), pointer:: dr_on_three, dr_three_on_eight
    real(dpf), dimension(:), pointer:: bool    ! r-grid weights for Bool rule
    integer:: nr                          ! number of points in r-grid
    integer, dimension(:), pointer:: jdouble ! j points where dx doubles
 
! rpow_f(r,l) = r**l, rpow_g(r,l) = 1/r**(l+1)
! i2_f = nr,  i1_g = 1
    integer:: ltmax                       ! used to set up arrays power_f(i,l) and power_g(i,l) for e-e Coulomb potential integration
    real(dpf), dimension(:,:), allocatable :: rpow_f, rpow_g
    real(dpf), dimension(:,:), pointer :: LegendreP, LegendreQ
    integer, dimension(:), allocatable:: i1_f,i2_g, i1_P,i2_Q
! note that i1_g = 1, and i2_f = nr
    
end type rgrid

!****  This is declaration of an object of the type rgrid.

 type(rgrid):: grid

!****

! This is public method for type rgrid 
 public setgrids

contains

!======================================================================      
!  This subroutine sets up radial grid:  grid. It is used to perform integrations
!  from zero to infinity using Simpson's rule. The structure is such
!  that a wave with momentum QMAX will have NPWAVE points per oscillation. This
!  determines HMAX, the largest step. At some point, XDBLE, the intervals,
!  dX, are progressively halved. This is done to accomodate irregular
!  solutions.
!  INPUT:
!    npwave  - number of points per oscillation
!    qmax  - the biggest momentum (a.u.) that can be calculated by this mesh.
!    ndoble - the number of doubling of intervals is NDOUBLE
!    rmax  - the largest "r=x" in the meshes.
!  OUTPUT:
!    gridr -  R grid
!    nr    - Number of X points
!    jdouble - j points where dx doubles
!    weight -  weights for Simpson rule
!    bool   -  weights for Bool rule
!
! It is assumed that wave function is exact zero at r=0 and at r=rmax (or at 
! corrsponding maxf, minf points). Therefore the first (r=0) and last (r=rmax)
! terms in Simpson composed rule are exact zero and r grid starts from 
! the second term and finishes at the one before last.
!======================================================================      
subroutine setgrids(self)
   use input_data
   use MPI_module
   implicit none
   type(rgrid) self

   integer:: npwave, npdbl, nleft, j, nj, nd, max_nj, jb, i1,i2
   real(dpf):: hmax, hmin, rdble, rleft, r, dr
   integer:: lna, nlimit
!!$ Resetting data_in%Rd at an exact grid point 
   integer:: point_count, interv_n, interv_n_min, interv_n_max
   real(dpf):: diff_i, diff_f, temp_n
   real(dpf):: rdble_interv


! This is r-grid parameters used in setgrids(grid)
   self%npwave = data_in%npwave
   self%npdbl = data_in%npdbl
   self%qmax = data_in%qmax
   self%ndouble = data_in%ndouble
   self%rmax = data_in%rmax
   self%formcut = data_in%formcut
   self%regcut = data_in%regcut
   self%expcut = data_in%expcut
   self%ltmax = data_in%ltmax



!  jdouble stores the J index for GRIDR where DR doubles, with the first
!  point set to 1 and the last to NR. Used in the numerov integrations.

!  Set up GRIDR

!  Make sure NPDBL is a multiple of 4 (required for Boole's rule - Liam changed from multiple of 2)
   self%npdbl=(self%npdbl/4) * 4
   if (self%npdbl.lt.4) then
      print*,'Need to have at least 4 equally spaced points in GRIDR'
      stop 'Need to have at least 4 equally spaced points in GRIDR'
   end if
  
   hmax = 3.14d0/dble(self%npwave)/self%qmax
   hmin = hmax/dble(2.0d0**self%ndouble)
!  The value of the R point from which dR is constant, = HMAX, is RDBLE
!  RDBLE = NPDBL * hmin * (2**NDOUBLE-1)
   rdble = dble(self%npdbl) * hmin * (2**self%ndouble-1)
!  The remaining part from rdble to rmax is:
   rleft = self%rmax - rdble
!  nleft = int(rleft / hmax) / 2 * 2
   nleft = int(rleft / hmax) / 4 * 4
!  The total number of points:
   self%nr = nleft + self%npdbl * self%ndouble

!!$ Increase the number of points in the rdble region
!!$ until an even index point is close to R/2
   if ( data_in%Rd == 0d0 .OR. data_in%calculation_type==2 .OR. data_in%calculation_type==3 ) then
      interv_n_min = self%ndouble + 2
      interv_n_max = interv_n_min   ! Not actually used in this case.
   else
      if ( data_in%Rd/2d0 >= rdble ) then
!!$ R/2 in the interval where dr=hmax
         j = INT((data_in%Rd/2d0-rdble)/hmax)
         if (mod(j,2)==1) j = j - 1 
         diff_i = ABS(data_in%Rd/2d0-(rdble+dble(j)*hmax))
      else
!!$ R/2 in intervals where dr doubles every npdbl points
         diff_f = dble(self%npdbl) * hmin
         interv_n = 1
!!$ Check which interval has R/2
         do while ( diff_f < data_in%Rd/2d0 )
            interv_n = interv_n + 1
            diff_f = diff_f + 2d0**(interv_n-1) * dble(self%npdbl) * hmin
         end do
!!$ interv_n = interval 1 dr=hmin, 2 dr=2*hmin etc  
!!$ i.e. Points around R/2 have dr=2d0**(interv_n-1)*hmin
!!$ diff_f last r value of the interval that contains R/2
!!$ temp_n =  r value at the beginning of the interval
!!$ that contains R/2
         temp_n = diff_f - 2d0**(interv_n-1) * self%npdbl * hmin
         j = INT( (data_in%Rd/2d0 - temp_n) / (2d0**(interv_n-1) * hmin ))
         if (mod(j,2)==1) j = j - 1
!!$ j is closest even index point below R/2
         diff_i = data_in%Rd/2d0 - temp_n - (2d0**(interv_n-1) * hmin * j)
      end if

!!$ diff_i = r difference between R/2 and the even index r grid point below R/2
      diff_f = diff_i
      interv_n = 1
!!$ interv_n_max and interv_n_min 
      interv_n_max = self%ndouble + 3 - interv_n
      interv_n_min = interv_n_max

!!$ Find the largest dr interval which we can add 2 points
      do while ( diff_i < 2d0 * hmax / (2d0**(interv_n-1)) .AND. interv_n <= self%ndouble+1 )
         interv_n =  interv_n + 1
         interv_n_max = self%ndouble + 3 - interv_n
      end do

!!$ Find smallest dr inerval in sequence with the largest that can add 2 points 
      do while ( diff_f > 2d0 * hmax / (2d0**(interv_n-1)) .AND. interv_n <= self%ndouble+1)
         diff_f = diff_f - 2d0 * hmax / (2d0**(interv_n-1))
         interv_n_min = self%ndouble + 3 - interv_n
         interv_n =  interv_n + 1
      end do

      point_count = 0
      rdble_interv = 0d0
      do j = interv_n_min, interv_n_max 
         rdble_interv = rdble_interv + 2d0 * hmax / (2d0**(self%ndouble+2-j)) 
         point_count = point_count + 2
      end do

!  RDBLE = NPDBL * hmin * (2**NDOUBLE-1)
      rdble = float(self%npdbl) * hmin * (2**self%ndouble-1) + rdble_interv
!  The remaining part from rdble to rmax is:
      rleft = self%rmax - rdble
!  nleft = int(rleft / hmax) / 2 * 2
      nleft = int(rleft / hmax) / 4 * 4!  The total number of points:
      self%nr = nleft + self%npdbl * self%ndouble + point_count 
   end if ! Rd=0 or Spheroidal

   if (myid == 0) then
      print*,' SETTING UP RADIAL GRID:'
      print*,'   > NDOUBLE:',self%ndouble
      print*,'   > HMIN:',hmin
      print*,'   > HMAX:',hmax
      print*,'   > NPDBL:',self%npdbl
      print*,'   > RDBLE:',rdble, self%nr - nleft
      print*,'   > NR:', self%nr
      write(*,*)
   endif

  allocate(self%gridr(self%nr))
  allocate(self%weight(self%nr))
  allocate(self%dr(self%nr), self%dr_on_three(self%nr), self%dr_three_on_eight(self%nr))
  allocate(self%bool(self%nr))
  allocate(self%jdouble(self%ndouble+2))

 !!$ Mark
  point_count = 0   
!!$
  self%jdouble(1) = 1
  do nd = 2, self%ndouble + 1
     self%jdouble(nd) = self%npdbl * (nd - 1)
! Mark: changed
     if (data_in%Rd /= 0d0 .AND. nd >= interv_n_min .AND. data_in%calculation_type /= 2 .AND. data_in%calculation_type /= 3) then
     if (interv_n_max >=  nd) point_count = point_count + 2
        self%jdouble(nd) = self%jdouble(nd) + point_count 
     end if
  end do
  self%jdouble(self%ndouble+2) = self%nr

  dr = hmin
  r = 0.0
  j = 0
  jb = 0
  do nd = 1, self%ndouble+1
! For all intervals except for the last one max_nj should be equal to npdbl, for first interval it does not give corect result and is corrected in the line below, for the last interval it should give nleft.
     max_nj = self%jdouble(nd+1) - self%jdouble(nd) 
     if(nd .eq. 1) max_nj = max_nj + 1
     do nj = 1, max_nj
        j = j + 1
        self%gridr(j) = r + float(nj) * dr
        self%dr(j) = dr
!  Simpson's rule weights
        self%weight(j) = float(mod(j,2) * 2 + 2) * dr / 3.0d0
     end do
     self%weight(j) = dr  !  this accounts for change of integration step at the boundary: (dr + 2*dr)/3 = dr
     r = self%gridr(j)
!
!  Bool's rule weights
     do nj=1,max_nj,4
        self%bool(jb+nj) = dr * 32.0d0
        self%bool(jb+nj+1) = dr * 12.0d0
        self%bool(jb+nj+2) = dr * 32.0d0
        self%bool(jb+nj+3) = dr * 14.0d0
     end do
     jb = self%jdouble(nd+1) 
     self%bool(jb) = dr * 21.0d0
!
     dr = dr * 2.0d0
  end do

  self%dr_on_three = self%dr/3.0d0
  self%dr_three_on_eight = self%dr*0.375d0

  self%weight(j) = hmax/3.0d0  ! for j=nr
  self%bool(self%nr) = 7.0d0
  self%bool(1:self%nr) = self%bool(1:self%nr) * 2.0d0 / 45.0d0

  !Liam added: use boole's rule if requested in data.in
  if(data_in%iweight == 1) self%weight = self%bool

  if (myid == 0)  print*,'   > Last R and NR:', self%gridr(self%nr), self%nr
  

!  call  rpow_construct(self%ltmax,self%nr,self%gridr,self%regcut,self%expcut)


  if( data_in%calculation_type == 0 ) then   ! Simple powers for the spherical case.

     if (myid == 0) print*, '    > Set rpow, ltmax =', self%ltmax

     allocate(self%rpow_f(1:self%nr,0:self%ltmax))
     allocate(self%rpow_g(1:self%nr,0:self%ltmax))
     allocate(self%i1_f(0:self%ltmax))
     allocate(self%i2_g(0:self%ltmax))

     self%i1_f(0) = 1
     self%i2_g(0) = self%nr
     self%rpow_f(1:self%nr,0)=1.0
     self%rpow_g(1:self%nr,0)=1.0/self%gridr(1:self%nr)

     do lna=1,self%ltmax

        self%i1_f(lna) = self%i1_f(lna-1)
        self%i2_g(lna) = self%i2_g(lna-1)
        self%rpow_f(1:self%nr,lna)=self%rpow_f(1:self%nr,lna-1)*self%gridr(1:self%nr)
        self%rpow_g(1:self%nr,lna)=self%rpow_g(1:self%nr,lna-1)/self%gridr(1:self%nr)
        
        do while (self%rpow_f(self%i1_f(lna),lna) .lt. self%regcut)
           self%i1_f(lna) = self%i1_f(lna)+1
        end do
        do while (self%rpow_g(self%i2_g(lna),lna) .lt. self%expcut*self%expcut)
           self%i2_g(lna) = self%i2_g(lna)-1
        end do

     end do

!!!$ Mark: Temporary changes to fix data_in%Rd onto a rgrid point that is
!!!$ at the beginning of a simpson integration (i.e. even index)  
!!!$ mod(odd,2) = 1, mod(even,2) = 0
     if ( data_in%Rd /= 0d0 ) then
        j = 4  
!!$ Diatomics Z/|r-R/2| point at r=R/2
        diff_i = ABS(self%gridr(j-2) - data_in%Rd/2d0 ) 
        diff_f = ABS(self%gridr(j) - data_in%Rd/2d0 ) 
        do while ( diff_f <= diff_i )
           diff_i = diff_f
           diff_f = ABS(self%gridr(j) - data_in%Rd/2d0 )
!!$ Keeping points even 
           j = j + 2 
        end do
        data_in%Rd = 2d0 * self%gridr(j-4)
        if (myid == 0) then
           print*,"    > Inter-nuclear distance R redefined",data_in%Rd 
        endif
     end if ! Rd/=0
     
  elseif( data_in%calculation_type==2 .or. data_in%calculation_type==3 ) then   ! Legendre polynomials for spheroidal case.
     nlimit = self%ltmax * (self%ltmax+3) / 2
     allocate( self%LegendreP(1:self%nr,0:nlimit), self%LegendreQ(1:self%nr,0:nlimit), self%i1_P(0:nlimit), self%i2_Q(0:nlimit) )
     
     call Legendre_construct(self, data_in%Rd)

     self%i1_P(:) = 1
     self%i2_Q(:) = self%nr

     do lna = 0, nlimit

        i1 = 1
        do while (i1<self%nr .and. abs(self%LegendreP(i1,lna))<self%regcut)
           i1 = i1 + 1
        end do
!        self%i1_P(lna) = i1

        i2 = self%nr
        do while (i2>1 .and. abs(self%LegendreQ(i2,lna))<self%expcut**2)
           i2 = i2 - 1
        end do
!        self%i2_Q(lna) = i2

     end do

  endif

  if(myid == 0) write(*,*)
end subroutine setgrids

subroutine destruct_gridr(self)

  implicit none
  type(rgrid) self

  if( allocated(self%rpow_f) ) then
     deallocate( self%rpow_f, self%rpow_g )
     deallocate(self%i1_f)
     deallocate(self%i2_g)
  elseif( allocated(self%i1_P) ) then
     deallocate(self%LegendreP,self%LegendreQ, self%i1_P,self%i2_Q)
  endif
  deallocate(self%gridr,self%weight,self%bool)

end subroutine destruct_gridr

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine Legendre_construct(self, Rd)
  !
  ! Calculates associated Legendre polynomials of argument greater than 1 for
  ! use in the electron-electron repulsion spheroidal partial wave expansion.
  ! Functions are stored as LegendreP(l,m,rho) and LegendreQ(l,m,rho) for the
  ! regular and irregular functions respectively.
  !
  ! Coded 12/10/12 by J.S.
  !
  use Associated_Legendre_Functions   ! The main module for calculating these polynomials.
  use Special_Functions   ! Variables controlling the calculations are stored in this module.
  use MPI_module
  implicit none

  type(rgrid), intent(inout) :: self
  integer :: n
  real(dpf), intent(in) :: Rd
  real(dpf) :: Ron2, tmp


  l_max = self%ltmax
  m_max = self%ltmax
  n_points = self%nr

  self%LegendreP(:,:) = 0d0
  self%LegendreQ(:,:) = 0d0
  
  if (Rd == 0d0) then
     n = 0
     self%LegendreP(:,0) = 1d0
     self%LegendreQ(:,0) = 1d0 / self%gridr(:)

     do l = 1, l_max
        n = n + l
        self%LegendreP(:,n)= self%LegendreP(:,n-l)*dble(2*l-1)*self%gridr(:)
        self%LegendreQ(:,n)= self%LegendreQ(:,n-l)/dble(2*l+1)/self%gridr(:)

        do m = 1, l
           self%LegendreP(:,n+m) = self%LegendreP(:,n)
           self%LegendreQ(:,n+m) = (-1)**m*self%LegendreQ(:,n)
        end do
     end do

  else
     Directive = 'both'
     recur = 'Wronskian'   ! Better for when both regular and irregular functions are calculated.
     Leg%D%B%Dir=recur
     
     allocate( Factor(0:l_max+m_max), x(n_points), Leg%R_LM%F(0:l_max,0:m_max), Leg%I_LM%F(0:l_max,0:m_max),&
	       Leg%R_L%F(0:l_max), Leg%R_LM%P(0:l_max,0:m_max,n_points), Leg%I_LM%Q(0:l_max,0:m_max,n_points) )
     x = 2d0/Rd * self%gridr(:) + 1d0
  
     if (myid==0) then
       write(*,*) " Calculating Legendre polynomials:"
       write(*,3) l_max, m_max, n_points, Directive, Control, normalize, eps, recur, test_wron
     endif
3    Format(/,10x,'Maximum L                             = ',i6,1x, &
          'Maximum M                      = ',i6, &
          /,10x,'Number of Points                      = ',i6,1x, &
          'Type Functions Computed        = ',a10,1x, &
          /,10x,'Type Calculation               = ',a24,1x, &
          'Normalize Regular Functions on (-1,1) = ',l1,1x, &
          /,10x,'Continued Fraction Convergence = ',1pe15.8,1x, &
          'Backward Recurrence Method            = ',a24,1x, &
          /,10x,'Test Wronskian                 = ',l1) 
     
     call Factorials
     call Legendre( R_LM=Leg%R_LM, I_LM=Leg%I_LM )

     ! Re-normalise the polynomials to avoid R dependence and factorials.
     n = 0
     Ron2 = 1d0
     do l = 0, l_max
        n = n + l
        tmp = Ron2   ! lambda! (R/2)^lambda
        do m = 0, l
!           self%LegendreP(l,m,:) = Leg%R_LM%P(l,m,:) * tmp
           self%LegendreP(:,n+m) = Leg%R_LM%P(l,m,:) * tmp
           if ((l-m) .eq. 0 ) cycle
           tmp = tmp / dble(l-m)   ! (lambda-mu)! (R/2)^lambda
        end do

        Ron2 = Ron2 * Rd/2d0; tmp = Ron2   ! lambda! (R/2)^(lambda+1)
        do m = 0, l
!           self%LegendreQ(l,m,:) = Leg%I_LM%Q(l,m,:) / tmp
           self%LegendreQ(:,n+m) = Leg%I_LM%Q(l,m,:) / tmp
           tmp = tmp * dble(l+m+1)   ! (lambda+mu)! (R/2)^(lambda+1)
        end do

        Ron2 = Ron2 * dble(l+1)   ! Puts the lambda! in lambda! (R/2)^lambda
     end do

     deallocate( Factor, x, Leg%R_LM%F, Leg%I_LM%F, Leg%R_L%F, Leg%R_LM%P, Leg%I_LM%Q )

  end if

!!$  do l = 0, l_max
!!$     do m = 0, m_max
!!$        open(100+10*l+m)
!!$        do i = 1, n_points
!!$           write(100+10*l+m,'(3Es20.12)') self%gridr(i), self%LegendreP(l,m,i), self%LegendreQ(l,m,i)
!!$        end do
!!$        close(100+10*l+m)
!!$     end do
!!$  end do
!!$stop!!!

end subroutine Legendre_construct


end module grid_radial

