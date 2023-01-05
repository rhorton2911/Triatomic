module spheroidal_waves_module
  contains
subroutine spheroidal_waves(mode, lambda,m,A,R,k, Zasym, Zplus,     Zminus, nTerms, lVec,DVec, jstart,jstop,reg, pshift)
  !
  use grid_radial
  use MPI_module
  implicit none

  integer, intent(in) :: mode, lambda,m
  real*8, intent(out) :: A
  real*8, intent(in) :: R,k, Zasym,Zplus,Zminus
  integer, intent(inout) :: nTerms
  integer, dimension(1:nTerms), intent(out), optional :: lVec
  real*8, dimension(1:nTerms), intent(out), optional  :: DVec
  integer, intent(out), optional :: jstart,jstop
  real*8, dimension(1:grid%nr), intent(out), optional :: reg
  real*8, intent(out), optional :: pshift

  integer :: par,par01, lMin,lMax, ind,increment
  real*8 :: c

  if (k < 0d0) then
     nTerms = 1
     if (abs(mode) >= 1) then
        lVec(:) = 0
        DVec(:) = 0d0
     end if
     if (abs(mode) >= 2) then
        jstart = 1
        jstop = 1
        reg(:) = 0d0
        pshift = 0d0
     end if
     return
  end if
  c = k*R/2d0

  ! The l's in lVec should now be >=m, have the same parity as l,
  lMin = abs(m)
  lMax = lambda + ceiling(2**(log10(abs(c))+3.0)) !This was Jeremy's rule - works fine for FN scattering

  !lMax = lambda + 20 !Liam added: want much more accurate V-matrix elements for VCC - overkilling the spherical l doesn't really
                     !  affect computation time so why not...
 
  if (Zminus == 0d0) then   ! Homogenous diatomic.
     par = (-1)**lambda   ! +1 for positive parity and -1 for negative parity.
     if ((-1)**lMin /= par) lMin = lMin + 1   ! lMin has the same parity as lam.
     if ((-1)**lMax /= par) lMax = lMax - 1   ! lMax has the same parity as lam.
     increment = 2   ! Skip every other l due to parity.
  else   ! Heterogenous diatomic -- no parity.
     increment = 1   ! Do not skip any l's.
  end if

  ! Number of terms in both the lVec and the DVec.
  if (mode == 0) nTerms = (lMax-lMin)/increment + 1

  ! Solve the angular differential equation with a basis of spherical harmonics.
  if (abs(mode) >= 1) then
     lVec(:) = (/( lMin + (ind-1)*increment, ind=1,nTerms )/)
     call spheroidal_angular(lambda,m,A,R,abs(k),Zminus, nTerms,lVec,DVec)
  end if

  ! Use A to solve the radial differential equation numerically.
  if (abs(mode) >= 2) call spheroidal_radial(mode, lambda,abs(m),A,R,abs(k),Zasym, nTerms,lVec,DVec, jstart,jstop,reg, pshift)


end subroutine spheroidal_waves
!
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
!
subroutine spheroidal_angular(lIn,mIn,A,R,k,Zminus, nTerms,lVec,DVec)
  !
  use MPI_module
  implicit none

  integer, intent(in) :: lIn,mIn, nTerms
  real*8, intent(out) :: A
  real*8, intent(in) :: R,k, Zminus
  integer, dimension(1:nTerms), intent(in) :: lVec
  real*8, dimension(1:nTerms), intent(out) :: DVec

  real*8, parameter :: acc=1E-6

  integer :: ind,indl
  real*8 :: cSquared, l,m, check
  real*8, dimension(nTerms) :: AVec
  real*8, dimension(nTerms,nTerms) :: recurMat

  ! Variables for the LAPACK routine DGEEV.
  character*1 :: JOBVL,JOBVR
  integer :: N, LDA, LDVL,LDVR, LWORK, INFO
  real*8, dimension(nTerms) :: WR,WI
  real*8, dimension(40*nTerms) :: WORK
  real*8, dimension(nTerms,nTerms) :: VL,VR


  A = 0d0
  cSquared = k*k * R*R / 4d0
  m = dble(mIn)
  DVec(:) = 0d0

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
  ! Build the coefficient matrix.
  recurMat(:,:) = 0d0; indl = 0
  do ind = 1, nTerms
     l = dble(lVec(ind))
     if (nint(l) == lIn) indl = ind

     if (ind > 1) recurMat(ind,ind-1) = cSquared * dsqrt( (l-1d0-m)*(l-1d0+m)*(l-m)*(l+m) /(2d0*l-3d0)/(2d0*l+1d0) ) /(2d0*l-1d0)
     recurMat(ind,ind) = l*(l+1d0) - 2d0*cSquared*(l*(l+1d0)+m*m-1d0)/(2d0*l-1d0)/(2d0*l+3d0)
     if (ind < nTerms) recurMat(ind,ind+1) = cSquared * dsqrt( (l+1d0-m)*(l+1d0+m)*(l+2d0-m)*(l+2d0+m) /(2d0*l+1d0)/(2d0*l+5d0) ) /(2d0*l+3d0)
  end do
  if (indl == 0) stop 'spheroidal_angular(): l not contained in lVec'

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
  ! Solve the recursion matrix for the eigenvalues.

  ! Set up DGEEV subroutine.
  JOBVL = 'N'   ! 1.(in) : Left eigenvectors of the matrix are not computed.
  JOBVR = 'V'   ! 2.(in) : Right eigenvectors of the matrix are computed.
  N = nTerms   ! 3.(in) : The order of the matrix.
!  A(:,:) = recurMat(:,:)   ! 4.(inout) : The matrix to be solved. (LDA,N)
  LDA = nTerms   ! 5.(in) : The leading dimension of the matrix.
  WR(:) = 0d0   ! 6.(out) : The real parts of the computed eigenvalues. (N)
  WI(:) = 0d0   ! 7.(out) : The imaginary parts of the computed eigenvalues. (N)
  VL(:,:) = 0d0   ! 8.(out) : If JOBVL='N', VL is not referenced. (LDVL,N)
  LDVL = nTerms   ! 9.(in) : The leading dimension of the array VL. [sic]
  VR(:,:) = 0d0   ! 10.(out) : The right eigenvectors are stored in the columns of VR. (LDVR,N)
  LDVR = nTerms   ! 11.(in) : The leading dimension of the matrix VR.
  WORK(:) = 0d0   ! 12.(workspace) : On exit, if INFO=0, WORK(1) returns the optimal (LWORK).
  LWORK = 40*nTerms   ! 13.(in) : The dimension of the array WORK. If JOBVX='V', LWORK>=4*N.
  INFO = 0   ! 14.(out) : Successful exit =0; illegal argument <0; unconverged eigenvalues >0.
  call DGEEV(JOBVL,JOBVR, N, recurMat,LDA, WR,WI, VL,LDVL, VR,LDVR, WORK,LWORK, INFO)

  if (myid==0 .and. INFO==0) then
!     print*, 'DGEEV exited successfully. Input & optimal LWORK :', LWORK, nint(WORK(1))
  else if (INFO /= 0) then
     print*, 'INFO =', INFO
     stop '! ERROR in spheroidal_waves.f90/spheroidal_angular() : DGEEV returned nonzero INFO !'
  endif

  ! Extract the requested A and DVec and make the latter orthonormal.
  A = WR(indl)
  DVec(:) = VR(:,indl)
  VR(:,indl) = VR(:,indl) / sign( dsqrt(sum(DVec(:)*DVec(:))), sum(DVec(:)) )
  DVec(:) = VR(:,indl)


end subroutine spheroidal_angular
!
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
!
subroutine spheroidal_radial(mode, lambda,m,A,R,k,Zasym, nTerms,lVec,DVec, jstart,jstop,fVec, phase)
  use distpot
  use grid_radial
  use MPI_module
  implicit none


  integer, intent(in) :: mode, lambda,m, nTerms
  real*8, intent(in) :: A,R,k,Zasym
  integer, dimension(nTerms), intent(in) :: lVec
  real*8, dimension(nTerms), intent(in) :: DVec
  integer, intent(out) :: jstart,jstop
  real*8, dimension(1:grid%nr), intent(out) :: fVec
  real*8, intent(out) :: phase

  integer :: j,jasym,jmatch, IFAIL
  real*8 :: regcut, c, rho,xi, radFn,deriv, f,fd, g,gd, norma,errorMax
  complex*16 :: h,hd, tmp
  real*8, dimension(1:grid%nr) :: rhoVec,xiVec, fdVec, ucentr,cntfug,PVec

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
  ! Definitions.

  if (k < 0d0) then   ! Forgeddaboudit !
     jstop = 1; return
  endif

  regcut = grid%regcut
  jstart = 0; jstop = grid%nr   ! Indices of the first and last nonzero values.
  rhoVec(:) = grid%gridr(:)   ! Radial coordinate equivalent to R/2*(xi-1).
  fVec(:) = 0d0; fdVec(:) = 0d0   ! Solution to the original D.E.
  phase = 0d0

  if (R == 0) then   ! Combined nuclear case: c*xi --> k*rho.
     c = k
     xiVec(:) = rhoVec(:)
  else   ! Molecular case: c is the spheroidal pseudo-momentum.
     c = R/2d0 * k
     xiVec(:) = 2d0/R * rhoVec(:) + 1d0
  end if

  ! Define the potential terms with respect to X(rho) later used in AB method.
  ucentr(:) = 0d0
  if (Zasym /= 0d0) ucentr(:) = -2d0*Zasym *(rhoVec(:)+R/2d0) /rhoVec(:) /(rhoVec(:)+R)
  cntfug(:) = A /rhoVec(:) /(rhoVec(:)+R)
  PVec(:) = (2d0*rhoVec(:) + R*dble(m+1)) /rhoVec(:) /(rhoVec(:)+R)


  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
  ! Generate the wave with an initial point/derivative, then numerically.

  f = regcut
  do while (abs(f)<=regcut .and. jstart<jstop)
     jstart = jstart + 1

     ! The first point & its derivative is calculated as Xi(xi).
     if (Zasym == 0) then   ! Netural target.
        xi = xiVec(jstart)
        call spheroidal_plane(lambda,m,R,c,xi, nTerms,lVec,DVec, f,g,fd,gd)
!print*, 'PLANE:', f,fd
!        call spheroidal_coulomb_start(m,A,R,k,Zasym, xi,regcut, f,fd)
!print*, 'COULOMB:', f,fd
!stop!!
     else   ! Charged target.
!        call spheroidal_coulomb_waves( m,A,R,k, Zasym, rho,regcut, f,fd, g,gd )
        xi = xiVec(jstart)
        call spheroidal_coulomb_start(m,A,R,k,Zasym, xi,regcut, f,fd)
     end if
  end do

  if (jstart >= jstop) then   ! Return a null wave.
     jstart = 1; jstop = 1
     return
  end if

  ! Propagate the wave in (rho) space using the AB method from this start point.
  fVec(jstart) = f; fdVec(jstart) = fd
  call AB_integ(m,R,k, ucentr,cntfug,PVec, jstart,fVec,fdVec)

!!$  ! Determine the starting point, i.e. first insignificant point.
!!$  do while ( abs(fVec(jstart))<regcut .and. jstart<jstop )
!!$     jstart = jstart + 1
!!$  end do
!!$  fVec(1:jstart-1) = 0d0; fdVec(1:jstart-1) = 0d0

!!$  ! Regenerate the wave if it starts out zero for too long a time.
!!$  if ( jstart > grid%jdouble(grid%ndouble+1) ) then   
!!$     if (Zasym == 0) then   ! What about when there is a distorting potential.
!!$        xi = xiVec(jstart)
!!$        call spheroidal_plane(lambda,m,R,c,xi, nTerms,lVec,DVec, f,g,fd,gd)
!!$     else
!!$!        call spheroidal_coulomb_waves( m,A,R,k, Zasym, rho,regcut, f,fd, g,gd )
!!$     end if
!!$
!!$     fVec(jstart) = f; fdVec(jstart) = fd
!!$     call AB_integ(m,R,k, ucentr,cntfug,PVec, jstart,fVec,fdVec)
!!$  end if

  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
  ! Determine the normalisation and phase, working in the coordinate xi.

  if (Zasym==0 .and. nstdw==0 .and. mode>0) then   ! Purely a plane wave.
     jmatch = jstop+1   ! Skip the matching loop.
     norma = 1d0   ! Normalisation should already be correct.
     phase = 0d0   ! Phase should also already be correct.
  else   ! Distorted and/or Coulomb wave--match when central potential vanishes.
     jmatch = jstart
     do while ( abs(ucentr(jmatch))>regcut .and. jmatch<jstop )
        jmatch = jmatch+1
     enddo
  end if

  ! Do matching: nonce for plane waves, once for other waves, & many for debug.
  j = 0; errorMax = 0d0
  do while (jmatch <= jstop)
     xi = xiVec(jmatch)
     radFn = fVec(jmatch)
     deriv = fdVec(jmatch)

     ! Plane wave matching to exact solutions for distorted wave phase shift.
     if (Zasym == 0) then
        call spheroidal_plane(lambda,m,R,c,xi, nTerms,lVec,DVec, f,g,fd,gd)

        ! Dmitry's magic code for determining the distorted wave phase shift.
        h = cmplx(g,f)
        hd = cmplx(gd,fd)
        tmp = (f*hd - fd*h) / (radFn*hd - deriv*h)
        norma = abs(tmp)
        tmp = tmp/norma
        phase = -norma * (radFn*fd - deriv*f ) / (g*fd - gd*f)
        tmp = cmplx( real(tmp), phase )
        phase = asin( aimag(tmp)/abs(tmp) )

     else   ! Coulomb wave.

     end if
     if (mode > 0) exit   ! Only need to match once unless debugging.

     ! Update the largest error.
     tmp = (f - radFn) * xi
     if ( abs(tmp) > abs(errorMax) ) then
        j = jmatch
        errorMax = tmp
     end if
     write(11,'(11F20.12)') xi, radFn,deriv, f,fd, sin(c*xi-lambda*asin(1d0))/c/xi, norma,phase, ucentr(jmatch),cntfug(jmatch),PVec(jmatch)

     jmatch = jmatch + 1

  end do

  if (mode < 0) then   ! Debugging.
     xi = xiVec(j)
     radFn = fVec(j); f = errorMax*xi + radFn
     errorMax = f - radFn
     write(*,'(A,I6)') 'Maximum radial error at grid point', j
     write(*,'(A10,F12.8)') 'where xi =', xi
     write(*,'(A10,F12.8)') 'analytic =', f
     write(*,'(A10,F12.8)') 'numeric =', radFn
     write(*,'(A10,F12.8)') 'diff. =', errorMax
     write(*,'(A10,F12.8,A/)') 'error =', abs(errorMax/f)*100, '%'
  end if

  ! Asymptotic form is (R/2)/xi sin(c*xi ...
  !                    = 1/(rho+R/2) sin(k*rho + c ...
  norma = norma * k   ! Factor 1/k out of each partial wave.
  fVec(:) = fVec(:) * norma   ! Normalise.
  fdVec(:) = fdVec(:) * norma


end subroutine spheroidal_radial

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
!
subroutine spheroidal_plane(lamIn,mIn,R,c,xi, nTerms,lVec,DVec, f,g,fd,gd)
  !
  use grid_radial
  implicit none

  integer, intent(in) :: lamIn,mIn, nTerms
  real*8, intent(in) :: R,c,xi
  integer, dimension(nTerms), intent(in) :: lVec
  real*8, dimension(nTerms), intent(in) :: DVec
  real*8, intent(out) :: f,g,fd,gd

  integer :: i,j, l, m, ifail
  real*8 :: pi, lam, D,tmp
  real*8, dimension(nTerms) :: coeffVec

  complex*16 :: xx,eta,zlmin, sig
  real*8 :: x
  integer :: nl, mode1,kfn,info
  complex*16, dimension(2*nTerms-1) :: fc,gc,fcp,gcp


  f = 0d0; g = 0d0; fd = 0d0; gd = 0d0
  coeffVec(:) = 0d0
  pi = acos(-1d0)
  lam = dble(lamIn); m = abs(mIn) !JSS not sure?!
        
  xx = cmplx(c*xi, 0d0)
  x = c*xi
  eta = (0d0,0d0)
  zlmin = cmplx(lVec(1), 0d0)
  nl = 2*nTerms-1
!!$x = c*xi
!!$eta = 0d0
!!$xlmin = lVec(1)
!!$lrange = 2*nTerms-1
  mode1 = 1   ! Calculate all of f, g, fd, and gd.
  kfn = 1   ! Calculate spherical Bessel functions.

!!$print*, 'xx =', xx
!!$print*, 'eta =', eta
!!$print*, 'zlmin =', zlmin
!!$print*, 'nl =', nl
!!$print*, 'mode1 =', mode1
!!$print*, 'kfn =', kfn

  ! This routine calculates spherical Bessel functions for various values of l.
  ! It can have issues when cc.f is not compiled in optimisation mode.
  call coulcc(xx,eta, zlmin,nl, fc,gc,fcp,gcp, sig, mode1,kfn,info)
  !call coul90(x,eta, xlmin,nl, rfc,rgc,rfcp,rgcp, kfn,ifail)

! need to change some variable types.

!!$print*, 'fc =', fc
!!$print*, 'gc =', gc
!!$print*, 'fcp =', fcp
!!$print*, 'gcp =', gcp
!!$print*, 'sig =', sig
!!$print*, 'info =', info

  do i = 1, nTerms
     l = lVec(i); D = DVec(i)
     if (D == 0d0) cycle

     tmp = (2d0*l + 1d0)/4d0/pi
     do j = l-m+1,l+m
        tmp = tmp * dble(j)
     end do
     coeffVec(i) = dsqrt(tmp) * D
     tmp = (-1)**((dble(l)-lam)/2d0) * coeffVec(i)

     f = f + tmp*fc(2*i-1)
     g = g + tmp*gc(2*i-1)
     fd = fd + c*tmp*fcp(2*i-1)
     gd = gd + c*tmp*gcp(2*i-1)
     if (mIn/=0 .and. R/=0d0) then
        fd = fd + dble(m)/(xi*xi*xi-xi)*tmp*fc(2*i-1)
        gd = gd + dble(m)/(xi*xi*xi-xi)*tmp*gc(2*i-1)
     end if

  end do

  tmp = 1d0 / sum(coeffVec(:))
  if (mIn/=0 .and. R/=0d0) tmp = (1d0-1d0/xi/xi)**(dble(m)/2d0) * tmp
  f = tmp * f
  g = tmp * g
  fd = tmp * fd
  gd = tmp * gd


end subroutine spheroidal_plane
!
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine AB_integ(m,R,k, ucentr,cntfug,PVec, j1,reg,regprime)
  !
  ! Solves the function over the whole grid from one point and its derivative.  
  !
  ! Coded mid-June 2012 by J.S.
  ! Incorporated grid spacing doublings and optimised early July 2012 by J.S.
  !
  use grid_radial
  implicit none

  integer, parameter :: mr_max=4
  
  integer, intent(in) :: m, j1
  real*8, intent(in) :: R, k
  real*8, dimension(grid%nr), intent(in) :: ucentr, cntfug, PVec
  real*8, dimension(grid%nr), intent(inout) :: reg,regprime

  integer :: startBlock,gridBlock, j, nr, mr
  real*8 :: h, rho, w1,w2
  real*8, dimension(-mr_max:0,2) :: f,g
  real*8, dimension(:,:), allocatable :: tmpf,tmpg
  real*8, dimension(2) :: newf,newg
  real*8, dimension(grid%nr) :: rhoVec, transVec


  nr = grid%nr
  if (j1 == nr) return
  rhoVec(:) = grid%gridr(:)
  rho = rhoVec(j1)
  
  ! Set up the Xi <--> X transformation vector ((xi+1)/(xi-1))^(m/2).
  if (m==0 .or. R==0d0) then
     transVec(:) = 1d0
  else
     transVec(:) = (1d0 + R/rhoVec(:))**(dble(m)/2d0)
  end if

  ! Transform input Xi(xi) --> X(rho) output.
  if (R > 0d0) then
     reg(j1) = transVec(j1) * reg(j1)
     regprime(j1) = 2d0/R*transVec(j1)*regprime(j1) - m*R/2d0/rho/(rho+R)*reg(j1)
  end if

  ! Find which doubling "block" in the grid we start in.
  startBlock = 1   ! Assume the first and increment until we find it.
  do while (j1 > grid%jdouble(startBlock+1))
     startBlock = startBlock + 1
  enddo
  
  ! Loop over separate grid spacings to deal with the grid doubling points.
  j = j1   ! Control the grid index via this counter rather than for loops.
  do gridBlock = startBlock, grid%ndouble+1

     h = rhoVec(j+1) - rhoVec(j)   ! Grid spacing for this block.
     f(:,:) = 0d0
     g(:,:) = 0d0
     mr = -1   ! The next loop will run at least once and increment this to 0.

     ! If we're in the first block the loop will run once to initialise vectors.
     ! Otherwise, it fills up the vectors with every second element backwards
     ! until mr_max or j1 (every other one because the spacing has doubled).
     do while ( j-2*(mr+1)>=j1 .and. mr<mr_max )
        mr = mr+1   ! Found an additional previous point to use.

        w1 = ucentr(j-2*mr) + cntfug(j-2*mr) - k*k
        w2 = -PVec(j-2*mr)

        f(-mr,1) = reg(j-2*mr)
        f(-mr,2) = regprime(j-2*mr)
        g(-mr,1) = regprime(j-2*mr)
        g(-mr,2) = w1*reg(j-2*mr) + w2*regprime(j-2*mr)
     enddo

     ! Loop for each grid index j in this block.
     do while (j < grid%jdouble(gridBlock+1))
        w1 = ucentr(j+1) + cntfug(j+1) - k*k
        w2 = -PVec(j+1)

        if (mr == mr_max) then   ! No need for temporary arrays
           call ABmethod( mr, w1,w2, h, f(:,:),g(:,:), newf(:),newg(:) )
        else   ! Make some temporary arrays.
           allocate( tmpf(-mr:0,2), tmpg(-mr:0,2) )
           tmpf(:,:) = f(-mr:0,:)
           tmpg(:,:) = g(-mr:0,:)
           call ABmethod( mr, w1,w2, h, tmpf(:,:),tmpg(:,:), newf(:),newg(:) )
           deallocate( tmpf, tmpg )
        endif

        reg(j+1) = newf(1)
        regprime(j+1) = newg(1)

        j = j+1
        mr = min( mr+1, mr_max )

        f(-mr:-1,:) = f(1-mr:0,:)   ! Shift the previous
        g(-mr:-1,:) = g(1-mr:0,:)   !  values back by one,
        f(0,:) = newf(:)            !  and append the newly
        g(0,:) = newg(:)            !  calculated values.
     enddo

  enddo

  ! Transform back from input X(rho) --> Xi(xi).
  if (R > 0d0) then
     reg(j1:nr) = reg(j1:nr) / transVec(j1:nr)
     regprime(j1:nr) = R/2d0/transVec(j1:nr)*regprime(j1:nr) + m*R*R/4d0/rhoVec(j1:nr)/(rhoVec(j1:nr)+R)*reg(j1:nr)
  end if


end subroutine AB_integ

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine AB_extend(m,A,R,k, Z, phasePts,scaling,gridLen, XPrev,XPrimePrev, xi1,xi2, X1,Xd1, X2,Xd2)
  !
  use MPI_module
  implicit none

  integer, parameter :: mr_max=4

  integer :: m, phasePts,scaling
  real*8, intent(in) :: A,R,k, Z, gridLen
  real*8, dimension(-mr_max:0), intent(inout) :: XPrev,XPrimePrev
  real*8, intent(out) :: X1,Xd1
  real*8, intent(inout) :: xi1,xi2, X2,Xd2

  integer :: nPts, ind
  real*8 :: c, rhoStart,rhoStop, h, w1,w2
  real*8, dimension(2) :: fNew,gNew
  real*8, dimension(-mr_max:0,2) :: f,g
  real*8, dimension(:), allocatable :: rhoVec, ucentr,cntfug,PVec
  

  c = k * R / 2d0
  rhoStart = R/2d0 * (xi2-1d0)
!  if(debug) print*, 'Extending radial grid from rho =', rhoStart

  h = dble(scaling) * R/2d0*(xi2-xi1) / dble(phasePts)
!  if(debug) print*, ' with a grid spacing of h =', h

  nPts = floor( dble(scaling) * R/2d0*(gridLen-1d0) / h )
!  if(debug) print*, ' over the number of grid points =', nPts

  allocate( rhoVec(-mr_max:nPts), ucentr(-mr_max:nPts), cntfug(-mr_max:nPts), PVec(-mr_max:nPts) )
  rhoVec(:) = (/( rhoStart + ind*h, ind = -mr_max,nPts )/)
  rhoStop = rhoVec(nPts)
!  if(debug) print*, ' until the value of rho =', rhoStop

!  call calc_potentials( m,A,R,c, Z, mr_max+1+nPts,rhoVec, ucentr,cntfug,PVec )

  f(:,1) = XPrev(:)
  f(:,2) = XPrimePrev(:)
  g(:,1) = XPrimePrev(:)

  do ind = -mr_max, 0
     w1 = ucentr(ind) + cntfug(ind) - k*k
     w2 = -PVec(ind)
     g(ind,2) = w1*XPrev(ind) + w2*XPrimePrev(ind)
  end do  

!  if(debug) open(9313)
  do ind = 1, nPts
     w1 = ucentr(ind) + cntfug(ind) - k*k
     w2 = -PVec(ind)
     call ABmethod( mr_max, w1,w2, h, f,g, fNew,gNew )
     
     f(-mr_max:-1,:) = f(1-mr_max:0,:)
     g(-mr_max:-1,:) = g(1-mr_max:0,:)
     f(0,:) = fNew(:)
     g(0,:) = gNew(:)

     if ( ind == nPts-phasePts ) then
        xi1 = 2d0/R * rhoVec(ind) + 1d0
        X1 = fNew(1)
        Xd1 = gNew(1)
     end if

!     if(debug) write(9313,'(6Es20.12)') 2d0/R*rhoVec(ind)+1d0, fNew(1),gNew(1), ucentr(ind),cntfug(ind),PVec(ind)
  end do
!  if(debug) close(9313)

  xi2 = 2d0/R * rhoStop + 1d0
  X2 = fNew(1)
  Xd2 = gNew(1)

  XPrev(:) = f(:,1)
  XPrimePrev(:) = g(:,1)

  deallocate ( rhoVec, ucentr, cntfug, PVec )
!  if(debug) print*


end subroutine AB_extend

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine ABmethod(mr,W21,W22,h,y,yp,res,resp)
  !
  !   Integrates for the next point (res) and its derivative (resp)
  !   given the previous mr points (y) and derivatives (yp),
  !   the grid spacing (h) and the potentials (W21 and W22) such that
  !      y" = W22 y' + W21 y
  !
  implicit none

  integer, intent(in):: mr
  real*8, intent(in):: h, W21, W22
  real*8, dimension(-mr:0,2), intent(in):: y, yp
  real*8, dimension(2), intent(out):: res, resp
  
  real*8:: W11, W12
  real*8:: Bm1
  real*8, dimension(0:4):: Ar, Br
  real*8:: ae11, ae12, ae21, ae22, Det, Re1,Re2
  integer:: m
  
  
  if(mr .gt. 4) then
     print*, 'mr > 4, diffeq.f: ABmethod(): AB method are not coded for this case.'
     stop
  endif
  
  Ar(:) = 0.0
  Br(:) = 0.0
  Bm1 = 0.0
  
  call ABmethod_coef(100+mr,Ar,Br,Bm1)
  
  Re1 = 0d0
  Re2 = 0d0
  do m=0,mr
     Re1 = Re1 + Ar(m)*y(-m,1) + h*Br(m)*yp(-m,1)
     Re2 = Re2 + Ar(m)*y(-m,2) + h*Br(m)*yp(-m,2)
  enddo
  
  W11 = 0d0   ! These are hard-coded here for
  W12 = 1d0   ! the degenerate case where g = f'.
  
  ae11 = 1d0 - Bm1 * (h*W11)
  ae12 =  -h * Bm1 * W12
  ae21 =  -h * Bm1 * W21
  ae22 = 1d0 - Bm1 * (h*W22)
  Det = ae11*ae22 - ae12*ae21
  
  res(1) = (Re1*ae22 - Re2*ae12) / Det
  res(2) = (Re2*ae11 - Re1*ae21) / Det
  
  resp(1) = W11*res(1) + W12*res(2)
  resp(2) = W21*res(1) + W22*res(2)
  
  
end subroutine ABmethod

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-!

subroutine ABmethod_coef(mr,Ar,Br,Bm1)
  !
  !   Provides the coefficients for the AB integration method.
  !
  implicit none
  integer, intent(in):: mr
  real*8, dimension(0:4), intent(out):: Ar, Br
  real*8, intent(out):: Bm1
  
  if(mr .eq. 104) then     !  for Adams-Moulton corrector-predictor method start count mr from 100.
     Ar(0) = -0.5d0
     Ar(1) =  0.5d0
     Ar(2) =  0.75d0
     Ar(3) =  0.25d0
     Ar(4) =  0d0   
     Bm1 =    1771D0/5760D0
     Br(0) =  9235D0/5760D0
     Br(1) =  5890D0/5760D0
     Br(2) =  4610D0/5760D0
     Br(3) =  35D0/5760D0
     Br(4) =  59D0/5760D0 
  elseif(mr .eq. 103) then  
     Ar(0) = 1d0
     Bm1 =   251D0/720D0
     Br(0) = 646D0/720D0
     Br(1) = -264D0/720D0
     Br(2) = 106D0/720D0
     Br(3) = -19D0/720D0
  elseif(mr .eq. 102) then  
     Ar(0) = 1d0
     Bm1 =   9D0/24D0
     Br(0) = 19D0/24D0
     Br(1) = -5D0/24D0
     Br(2) = 1D0/24D0
  elseif(mr .eq. 101) then  
     Ar(0) = 1d0
     Bm1 =   5D0/12D0
     Br(0) = 8D0/12D0
     Br(1) = -1D0/12D0
  elseif(mr .eq. 100) then  
     Ar(0) = 1d0
     Bm1 =   1D0/2D0
     Br(0) = 1D0/2D0
     
  elseif(mr .eq. 4) then
     Ar(0) = 1.0
     Br(0) = 1901.0/720.0
     Br(1) = -2774.0/720.0
     Br(2) = 2616.0/720.0     
     Br(3) = -1274.0/720.0     
     Br(4) = 251.0/720.0     
  elseif(mr .eq. 3) then
     Ar(0) = 1.0
     Br(0) = 55.0/24.0
     Br(1) = -59.0/24.0
     Br(2) = 37.0/24.0     
     Br(3) = -9.0/24.0     
  elseif(mr .eq. 2) then
     Ar(0) = 1.0
     Br(0) = 23.0/12.0
     Br(1) = -16.0/12.0
     Br(2) = 5.0/12.0
  elseif(mr .eq. 1) then
     Ar(0) = 1.0 
     Br(0) = 1.5
     Br(1) = -0.5
  endif
  
end subroutine ABmethod_coef

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine spheroidal_coulomb_start(m,A,R,k,Z, xi,regcut, f,fd)
  !
  implicit none

  integer, intent(in) :: m
  real*8, intent(in) :: A,R,k,Z, xi,regcut
  real*8, intent(out) :: f,fd

  integer :: s
  real*8 :: c2,rm,t,nu, fCoeff,fdCoeff, term,rs, pbar,qbar,rbar,tbar
  real*8, dimension(-2:10) :: gVec


  c2 = k*k*R*R/4d0
  rm = dble(m)
  t = xi - 1d0
  fCoeff = (xi*xi-1d0)**(rm/2d0)
  fdCoeff = rm*xi / (xi*xi-1d0)
  nu = 0d0
  gVec(-2) = 0d0; gVec(-1) = 0d0; gVec(0) = 1d0

  s = 0; rs = 0d0
  f = t**nu
  fd = nu * t**(nu-1d0)

  term = regcut
  do while ( abs(term) >= regcut )
     pbar = 2d0*(nu+rs+1d0)*(nu+rs+rm+1d0)
     qbar = (nu+rs)*(nu+rs+2d0*rm+1d0) + A + R*Z + rm*(rm+1d0)
     rbar = 2d0*c2 + R*Z
     tbar = c2
     gVec(s+1) = (qbar*gVec(s) + rbar*gVec(s-1) * tbar*gVec(s-2)) / (-pbar)

     s = s+1; rs = dble(s)
     term = gVec(s) * t**(rs+nu)
     f = f + term
     fd = fd + (rs+nu)*term/t
print*
print*, 's =', s
print*, 'pbar,qbar,rbar,tbar:', pbar,qbar,rbar,tbar
print*, 'g_s =', gVec(s)
print*, 'term =', term
print*, 'f =', f
print*, 'fd =', fd
  end do


end subroutine spheroidal_coulomb_start

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine spheroidal_coulomb_waves( m,A,R,k, Z, rho,regcut, f,fd, g,gd )
  !
  ! Generate a spheroidal Coulomb wave using the power series expansion of
  !   Ponomarev and Sonov (1976) J. Comp. Phys. 20: 183.
  !
  ! Calculated in coordinate t = c(xi-1) = k*rho but converted back to rho.
  !
  use MPI_module
  implicit none

  integer, parameter :: nMax = 1000   ! Maximum number of terms in the power series expansion.

  integer, intent(in) :: m
  real*8, intent(in) :: A,R,k, Z, rho,regcut
  real*8, intent(out) :: f,fd, g,gd

  integer :: n
  real*8 :: c, t, mu, Df,Dfd, Dg,Dgd
  real*8, dimension(-2:1) :: sigVec,sVec, sigdVec,sdVec


  c = k*R/2d0
  t = k*rho

  n = 0
  mu = 0

  sigVec(-2) = 1d0
  sigVec(-1) = k*R + 2d0*Z/k
  sigVec(0) = (n+mu)*(n+mu+1) + A + Z*R + c*c
  sigVec(1) = k*R * (n+mu+1) * (n+mu+1+m)

  sigdVec(-2:-1) = 0d0
  sigdVec(0) = 2*n+2*mu+1
  sigdVec(1) = k*R * (2*n+2*mu+2+m)

  sVec(-2:-1) = 0d0
  sVec(0) = 1d0
  sdVec(-2:0) = 0d0

  Df = t**dble(mu)
  Dfd = t**dble(mu-1) * mu
  Dg = 0d0
  Dgd = 0d0

  f = Df
  fd = Dfd
  g = 0d0
  gd = Dgd
!!$  if(debug) then
!!$     write(*,'(/A,I2,4(4X,A,F20.12))') 'm = ', m, 'A = ', A, 'k = ', k, 'R = ', R, 'Z = ', Z
!!$     write(*,'(4(A,Es20.12,4X))') 'rho = ', rho, 'regcut = ', regcut, 'f = ', f, 'fd = ', fd
!!$  end if

  do while ( max( max(abs(Df),abs(Dfd)), max(abs(Dg),abs(Dgd)) )>regcut .and. n<nMax )
     n = n+1
     sVec(1) = -sum( sigVec(-2:0)*sVec(-2:0) ) / sigVec(1)
     sdVec(1) = -( sum(sigdVec(-2:1)*sVec(-2:1)) + sum(sigVec(-2:0)*sdVec(-2:0)) ) / sigVec(1)
!!$     if(debug) then
!!$        write(*,'(A,4Es20.12)') 'sigVec =', sigVec
!!$        write(*,'(A,4Es20.12)') 'sVec =', sVec
!!$        write(*,'(A,4Es20.12)') 'sigdVec =', sigdVec
!!$        write(*,'(A,4Es20.12)') 'sdVec =', sdVec
!!$        write(*,'(A,2Es20.12)') 'Recurrent relations =', &
!!$             sum( sigVec(:)*sVec(:) ), sum( sigdVec(:)*sVec(:) + sigVec(:)*sdVec(:) )
!!$     end if

     sigVec(0) = sigVec(0) + 2d0*(n+mu)
     sigVec(1) = sigVec(1) + k*R*(2*n+2*mu+1+m)
     sigdVec(0) = sigdVec(0) + 2d0
     sigdVec(1) = sigdVec(1) + 2d0*k*R

     sVec(-2:0) = sVec(-1:1)
     sdVec(-2:0) = sdVec(-1:1)

     Df = sVec(0) * t**dble(n+mu)
     Dfd = sVec(0) * t**dble(n+mu-1) * (n+mu)
     Dg = sdVec(0) * t**dble(n+mu)
     Dgd = sdVec(0) * t**dble(n+mu-1) * (n+mu)

     f = f + Df
     fd = fd + Dfd
     g = g + Dg
     gd = gd + Dgd
!     if(debug) write(*,'(A,I4,6(4X,A,Es20.12))') 'n = ', n, 'Df =',Df, 'f =',f, 'Dfd =',Dfd, 'fd =',fd, 'Dg =',Dg, 'g =',g, 'Dgd =',Dgd, 'gd =',gd
  enddo
  if ( n == nMax ) print*, 'WARNING --- spheroidal_waves.f90/spheroidal_coulomb_waves() --- power series expansion did not converge for k, rho', k, rho
  
  if ( m == 0 ) then
     g = g + 2d0*log(t)*f
     gd = gd + 2d0/t*f + 2d0*log(t)*fd
  else
     g = 0d0
     gd = 0d0
  end if

  fd = k * fd   ! Derivative w.r.t. rho instead of t=k*rho
  gd = k * gd   !   so X'(rho) = d(t)/d(rho) X'(t) = k*X'(t).


end subroutine spheroidal_coulomb_waves

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine HadingerAH( l, m,A,R,k, Z, xi,regcut, XIn,dXIn, norma,phase )
  !
  ! Hadinger, G., M. Aubert-Frecon, G. Hadinger (1996) J. Phys. B: At. Mol. Opt. Phys. 29: 2951.
  ! Uses Ponomarev and Sonov convention (converges to spherical Coulomb phase in atomic limit).
  !
  implicit none

  integer, parameter :: nMax=1000

  integer, intent(in) :: l,m
  real*8, intent(in) :: A,R,k, Z, xi,regcut, XIn,dXIn
  real*8, intent(out) :: norma,phase

  integer :: n
  real*8 :: gamma,c, x, pi=acos(-1d0)
  complex*16 :: rho1, rho2, rho3, Cxi, tmp,alpha, F,dF, i=(0d0,1d0)
  complex*16, dimension(-2:1) :: Fn, dFn


  phase = 0d0
  Fn(:) = 0d0
  dFn(:) = 0d0

  gamma = Z/k
  c = k*R/2d0
  x = c * (xi+1d0) / 2d0
  Cxi = exp(cmplx( 0d0, c*xi + gamma*log(xi+1d0) )) / (xi+1d0)

  n = 0

  rho1 = cmplx( 0d0 , 4d0 )
  rho2 = cmplx( gamma*gamma - 2d0*gamma*c - c*c - A , gamma - 2d0*c*dble(m+1) )
  rho3 = cmplx( -gamma*gamma*c , -gamma*c*dble(m) )

  Fn(0) = 1d0
  dFn(0) = cmplx( -1d0/(xi+1d0) , c+gamma/(xi+1d0) )

  do while ( abs(Fn(0)-Fn(-1))>regcut .and. abs(dFn(0)-dFn(-1))>regcut .and. n<nMax )
     Fn(1) = -( (rho2/x-rho1)*Fn(0) + (rho3/x/x-rho2/x)*Fn(-1) - rho3/x/x*Fn(-2) ) / rho1
     dFn(1) = -( (rho2/x-rho1)*dFn(0) + (rho3/x/x-rho2/x)*dFn(-1) - rho3/x/x*dFn(-2) ) / rho1 + ( c/2d0/x/x*rho2*(Fn(0)-Fn(-1)) + c/x/x/x*rho3*(Fn(-1)-Fn(-2)) )/rho1

     n = n+1

     Fn(-2:0) = Fn(-1:1)
     dFn(-2:0) = dFn(-1:1)

     rho1 = rho1 + cmplx( 0d0, 4d0 )
     rho2 = rho2 + cmplx( -2d0*n, 2d0*gamma - 4d0*c )
     rho3 = rho3 + cmplx( c*dble(2*n+m-1) , -2d0*gamma*c )

  enddo
  if ( n == nMax ) print*, 'WARNING --- spheroidal_waves.f90/HadingerAH() --- power series expansion did not converge for k, xi :', k, xi

  F = Cxi * Fn(0)
  dF = Cxi * dFn(0)

  tmp = ( XIn * dF - dXIn * F )
  alpha = i/2d0 * log(-tmp/conjg(tmp))
  norma = XIn / ( 2d0 * real(exp(i*alpha)*F) )

  if ( norma < 0d0 ) alpha = alpha + pi
  norma = abs(norma)

!  phase = real(alpha) - gamma*log(c) + pi/2d0   ! Rankin and Thorson convention.
  phase = real(alpha) - gamma*log(k*R) + dble(l+1)*pi/2d0   ! Ponomarev and Somov convention.
  phase = modulo( phase+pi, 2d0*pi ) - pi


end subroutine HadingerAH

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine Coulomb_xi( l, m,A,R,k, Z, xi,regcut, XIn,dXIn, norma,phase )
  !
  ! Rankin, J. and W.R. Thorson (1978) J. Comp. Phys. 32: 437
  !
  implicit none

  integer, intent(in) :: l, m
  real*8, intent(in) :: A,R,k,Z, xi,regcut, XIn,dXIn
  real*8, intent(out) :: norma, phase

  integer, parameter :: nMax=1000
  real, parameter :: pi=acos(-1d0)

!  character*52 :: fmt='("n =",I3,4(" ;   ",A," = (",Es12.4,",",Es12.4,")"))'
  character*45 :: fmt='(A,I4,4X,A,4(" (",Es12.4,",",Es12.4,")"))'
  integer :: n, IFAIL
  real*8 :: c, f,g, fd,gd, Fn,Gn, Fnd,Gnd, beta
  real*8 :: deltac
  complex*16 :: h,hd, HPlus,HPrime, tmp
  complex*16, dimension(-2:1) :: aVec, bVec
  

  c = k * R / 2d0

  call coul90( c*xi, -Z/k, l,0, f,g, fd,gd, 0, IFAIL )
  h = cmplx( f, g )
  hd = cmplx( fd, gd )

!!$  if (debug) then
!!$     write(*,'(/3(A,F20.12,9X))') 'xi =', xi, 'XIn =', XIn,'dXIn =', dXIn
!!$     write(*,fmt) 'l =',l, 'h+, h-, h+`, h-`', h,conjg(h), hd,conjg(hd)
!!$     print*
!!$  end if


  ! Initialise the values in the recurrence relation.
  n = 0

  aVec(-2) = c*c/4d0 * l*(l+1) * h
  aVec(-1) = -c*c*Z/k * h
  aVec(0) = -(l*(l+1) + A + c*c) * h - 2*c*m * hd
  aVec(1) = 4d0*hd
!  if (debug) write(*,fmt) 'n =', n, 'aVec :', aVec

  bVec(-2:-1) = 0d0
  bVec(0) = 1d0

  HPrime = 0d0
  HPlus = 1d0
  
  ! Loop through values of n to achieve convergence.
  tmp = regcut
  do while ( abs(tmp)>=regcut .and. n<nMax )
     bVec(1) = sum( aVec(-2:0) * bVec(-2:0) ) / aVec(1)

     n = n+1
     aVec(-2) = aVec(-2) - c*c/2d0*(n-1) * h
     aVec(-1) = aVec(-1) + c*m * h - 2d0*c*c * hd
     aVec(0) = aVec(0) - 2d0*n * h
     aVec(1) = aVec(1) + 4d0 * hd
!     if (debug) write(*,fmt) 'n =', n, 'aVec :', aVec

     bVec(-2:0) = bVec(-1:1)

     tmp = bVec(0) * (2d0/c/xi)**n
     HPlus = HPlus + tmp
     HPrime = HPrime - n*tmp
!     if (debug) write(*,fmt) '    b_', n, ', tmp, HPlus, HPrime :', bVec(0), tmp, HPlus,HPrime

  end do

!  if ( n == nMax ) stop ' ! WARNING !   Coulomb (spheroidal_wave.f90) did not converge.'

  HPrime = (-h/c/xi/xi + hd/xi) * HPlus   +   h/c/xi/xi * HPrime
  HPlus = h/c/xi * HPlus
!  if (debug) write(*,fmt) 'nStop =',n, 'H+, H-, H+`, H-`', HPlus,conjg(HPlus), HPrime,conjg(HPrime)

  Fn = real(HPlus)
  Gn = imag(HPlus)
  Fnd = real(HPrime)
  Gnd = imag(HPrime)
  beta = -atan2( XIn*Fnd - dXIn*Fn , XIn*Gnd - dXIn*Gn )

  ! Determine the normalisation constant.
  norma = XIn / ( cos(beta)*Fn + sin(beta)*Gn )
  if ( norma < 0d0 ) then   ! Restrict to positive to prevent arbitrary inclusion of e^(i*pi)=-1.
     beta = beta + pi
     norma = -norma
  endif

  ! Determine the phase shift.
  phase = beta + deltac( -Z/k, l )
  phase = modulo( phase+pi, 2d0*pi ) - pi   ! Range (-pi,pi].
!  if (debug) write(*,'(/3(A,F12.8,9X))') 'norma =', norma, 'phase =', phase, 'beta =', beta

!stop!
end subroutine Coulomb_xi

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine Coulomb_x( l, m,A,R,k, Z, xi,regcut, XIn,dXIn, norma,phase )
  !
  ! Rankin, J. and W.R. Thorson (1978) J. Comp. Phys. 32: 437
  !
  implicit none

  integer, intent(in) :: l, m
  real*8, intent(in) :: A,R,k,Z, xi,regcut, XIn,dXIn
  real*8, intent(out) :: norma, phase

  integer, parameter :: nMax=1000
  real, parameter :: pi=acos(-1d0)

  character*52 :: fmt='("n =",I3,4(" ;   ",A," = (",Es12.4,",",Es12.4,")"))'
  integer :: n, IFAIL
  real*8 :: c,x, f,g, fd,gd, Fn,Gn, Fnd,Gnd, beta
  real*8 :: deltac
  complex*16 :: h,hd, aPrev,aCurr,aNext, bPrev,bCurr,bNext, HPlus,HPrime, tmp
  

  c = k * R / 2d0
  x = c * (xi+1d0)

  call coul90( x, -Z/k, l,0, f,g, fd,gd, 0, IFAIL )
  h = cmplx( g, f )
  hd = cmplx( gd, fd )

!!$  if (debug) then
!!$     write(*,'(/3(A,F20.12,9X))') 'x =', x, 'XIn =', XIn, 'dXIn =', dXIn
!!$     write(*,fmt) 0, 'h+(t)', h, 'h+`(t)', hd, 'h-(t)', cmplx(g,-f), 'h-`(t)', cmplx(gd,-fd)
!!$     print*
!!$  end if


  ! Initialise the values in the recurrence relation.
  n = 0

  aPrev = c * l*(l+1) * h
  aCurr = -( l*(l+1) + A + c*c - Z*R ) * h - 2d0*c * (m+1) * hd
  aNext = 4d0 * hd

  bPrev = 0d0
  bCurr = 1d0
!  if (debug) write(*,fmt) n, 'aPrev: ',aPrev, 'aCurr: ',aCurr, 'aNext: ',aNext, 'bCurr: ',bCurr

  HPrime = 0d0
  HPlus = 1d0
  
  ! Loop through values of n to achieve convergence.
  tmp = regcut
  do while ( abs(tmp)>=regcut .and. n<nMax )
     bNext = -( aPrev*bPrev + aCurr*bCurr ) / aNext

     n = n+1
     aPrev = aPrev + c*(2*n+m-1)*h       ! Add the difference
     aCurr = aCurr - 2*n*h - 4*c*hd   !  instead of re-calculating
     aNext = aNext + 4*hd                !  to hopefully save time.
     bPrev = bCurr
     bCurr = bNext
!     if (debug) write(*,fmt) n, 'aPrev :',aPrev, 'aCurr :',aCurr, 'aNext :',aNext, 'bCurr :',bCurr

     tmp = bCurr * (2d0/x)**n
     HPrime = HPrime - n*tmp
     HPlus = HPlus + tmp
!     if (debug) write(*,fmt) n, 'H+',HPlus, 'H+`',HPrime, 'H-',conjg(HPlus), 'H-`',conjg(HPrime)

  end do

!  if ( n == nMax ) stop ' ! WARNING !   Coulomb (spheroidal_wave.f90) did not converge.'

  HPrime = c*( -h/x + hd )/x*HPlus + c*h/x/x*HPrime
  HPlus = c*h/x*HPlus
!  if (debug) write(*,fmt) n, 'H+', HPlus, 'H+`', HPrime, 'H-', conjg(HPlus), 'H-`', conjg(HPrime)

  HPrime = HPrime * c
!  tmp = ( XIn*HPrime - dXIn*HPlus )   ! Equivalent to the Wronskian of X and F+.
!  tmp = ( 0d0 , 0.5d0 ) * log(-tmp/conjg(tmp))
!  beta = real(tmp)   ! beta = delta + c - gamma*ln(2) - pi/2 = Delta + c - (l+1)pi/2
  Gn = real(HPlus)
  Fn = imag(HPlus)
  Gnd = real(HPrime)
  Fnd = imag(HPrime)
  beta = atan2( -XIn*Gnd + dXIn*Gn , XIn*Fnd - dXIn*Fn )

  ! Determine the normalisation constant.
!  tmp = exp(cmplx( 0d0 , beta )) * HPlus
!  norma = XIn / ( 2d0*real(tmp) )
  norma = XIn / ( cos(beta)*Gn + sin(beta)*Fn )
  if ( norma < 0d0 ) then   ! Restrict to positive to prevent arbitrary inclusion of e^(i*pi)=-1.
     beta = beta + pi
     norma = -norma
  endif

  ! Determine the phase shift.
!  phase = beta + c + pi/2d0 + deltac(-Z/k,0)   ! Ponomarev and Somov convention.
  phase = -beta + c + pi/2d0 + deltac( -Z/k, l )
  phase = modulo( phase+pi, 2d0*pi ) - pi   ! Range (-pi,pi].
!  if (debug) write(*,'(/3(A,F12.8,9X))') 'norma =', norma, 'phase =', phase, 'beta =', beta

!stop!
end subroutine Coulomb_x

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine Rankin_Thorson( l, m,A,R,k, Z, xi,regcut, radFn,deriv, D,phase )
  !
  ! Rankin, J. and W.R. Thorson (1978) J. Comp. Phys. 32: 437
  !
  implicit none

  integer, intent(in) :: l, m
  real*8, intent(in) :: A,R,k,Z, xi,regcut
  real*8, intent(inout) :: radFn, deriv
  real*8, intent(out) :: D, phase

  integer, parameter :: nMax = 1000
  real, parameter :: pi = acos(-1d0)

  character*52 :: fmt='("n =",I2,4(" ;   ",A," = (",Es12.4,",",Es12.4,")"))'
  integer :: n
  real*8 :: alpha
  complex*16 :: aPrev, aCurr, aNext, bPrev, bCurr, bNext, FPlus, FPrime, tmp
  

!  if (debug) write(*,'(/3(A,F20.12,9X)/)') 'xi =', xi, 'radFn =', radFn, 'deriv =', deriv

  ! Initialise the values in the recurrence relation.
  n = 0
  aPrev = cmplx( -Z/k/2d0 , -dble(m)*Z*R/2d0 )
  aCurr = cmplx( (Z/k)**2 - Z*R - (k*R/2d0)**2 - A , Z/k - dble(m+1)*k*R )
  aNext = ( 0d0 , 4d0 )
  bPrev = 0d0
  bCurr = 1d0
  FPlus = 1d0
  FPrime = cmplx( 1d0 , -Z/k )
!  if (debug) write(*,fmt) n, 'aCurr', aCurr, 'bCurr', bCurr, 'u+', FPlus, 'u+`', FPrime
  
  ! Loop through values of n to achieve convergence.
  tmp = 1d0
  do while ( abs(tmp) > regcut .and. n < nMax )
     bNext = -( aPrev*bPrev + aCurr*bCurr ) / aNext

     n = n+1
     aPrev = aPrev + cmplx( dble(2*n+m-1)*k*R/2d0 , -Z*R )   ! Add the difference
     aCurr = aCurr + cmplx( -2d0*n , 2d0*Z/k - 2d0*k*R )     !  instead of re-calculating
     aNext = aNext + ( 0d0 , 4d0 )                           !  to hopefully save time.
     bPrev = bCurr
     bCurr = bNext

     tmp = bCurr * (4d0/k/R/(xi+1))**n
     FPlus = FPlus + tmp
     FPrime = FPrime + cmplx( n+1d0 , -Z/k )*tmp
!     if (debug) write(*,fmt) n, 'aCurr', aCurr, 'bCurr', bCurr, 'u+', FPlus, 'u+`', FPrime
  end do

  if ( n == nMax ) stop 'WARNING --- spheroidal_waves.f90/Rankin_Thorson() did not converge.'

  tmp = exp( cmplx( 0d0 , k*R*xi/2d0 + Z/k*log(xi+1d0) ) )
  FPlus = FPlus * tmp / (xi+1d0)
  FPrime = -FPrime * tmp / (xi+1d0) / (xi+1d0) + FPlus * cmplx( 0d0 , k*R/2d0 )
!  if (debug) write(*,fmt) n, 'F+', FPlus, 'F+`', FPrime, 'F-', conjg(FPlus), 'F-`', conjg(FPrime)

  tmp = ( radFn*FPrime - deriv*FPlus )   ! Equivalent to the Wronskian of X and F+.
  tmp = ( 0d0 , 0.5d0 ) * log(-tmp/conjg(tmp))
  alpha = real(tmp)   ! alpha according to Rankin and Thorson 1979.

  ! Determine the normalisation constant.
  tmp = exp(cmplx( 0d0 , alpha )) * FPlus
  D = radFn / ( 2d0*real(tmp) )
  if ( D < 0d0 ) then   ! Restrict to positive D to prevent arbitrary inclusion of e^(i*pi)=-1.
     alpha = alpha + pi
     tmp = exp(cmplx( 0d0 , alpha )) * FPlus
     D = radFn / ( 2d0*real(tmp) )
  endif

  radFn = 2d0*real(tmp)
  tmp = exp(cmplx( 0d0 , alpha )) * FPrime
  deriv = 2d0*real(tmp)

  ! Determine the phase shift.
!  phase = alpha - Z/k*log(k*R/2d0) + pi/2d0   ! Rankin and Thorson convention.
  phase = alpha - Z/k*log(k*R) + dble(l+1)*pi/2d0   ! Ponomarev and Somov convention.
  phase = modulo( phase+pi, 2d0*pi ) - pi
!  if (debug) write(*,'(/3(A,F12.8,9X))') 'D =', D, 'phase =', phase, 'alpha =', alpha


end subroutine Rankin_Thorson

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
!
!subroutine Coulomb_t( l, m,A,R,k, Z, xi,regcut, XIn,dXIn, norma,phase )
!  !
!  ! Rankin, J. and W.R. Thorson (1978) J. Comp. Phys. 32: 437
!  !
!  implicit none
!
!  integer, intent(in) :: l, m
!  real*8, intent(in) :: A,R,k,Z, xi,regcut, XIn,dXIn
!  real*8, intent(out) :: norma, phase
!
!  integer, parameter :: nMax=1000
!  real, parameter :: pi=acos(-1d0)
!
!  character*52 :: fmt='("n =",I3,4(" ;   ",A," = (",Es12.4,",",Es12.4,")"))'
!  integer :: n, IFAIL
!  real*8 :: c,t, f,g, fd,gd, beta
!  real*8 :: deltac
!  complex*16 :: h,hd, aPrev,aCurr,aNext, bPrev,bCurr,bNext, HPlus,HPrime, tmp
!!  logical :: debug=0
!  
!
!  c = k * R / 2d0
!  t = c * (xi-1d0)
!  if ( c/t > 2d0 ) then
!     norma = 1d0/0d0
!     phase = 1d0/0d0
!     return
!  end if
!
!  call coul90( t, -Z/k, l,0, f,g, fd,gd, 0, IFAIL )
!  h = cmplx( g, f )
!  hd = cmplx( gd, fd )
!
!!!$  if (debug) then
!!!$     write(*,'(/3(A,F20.12,9X))') 't =', t, 'XIn =', XIn, 'dXIn =', dXIn
!!!$     write(*,fmt) l, 'h(t)',h, 'h`(t)',hd, 'g(t)-i*f(t)',cmplx(g,-f), 'g`(t)-i*f`(t)',cmplx(gd,-fd)
!!!$     print*
!!!$  end if
!
!
!  ! Initialise the values in the recurrence relation.
!  n = 0
!  aPrev = c * l*(l+1) * h
!  aCurr = (l*(l+1) + A + c*c - Z*R ) * h  +  2d0*c * (m-1) * hd
!  aNext = -4d0 * hd
!  bPrev = 0d0
!  bCurr = 1d0
!!  if (debug) write(*,fmt) n, 'aPrev', aPrev, 'aCurr', aCurr, 'aNext', aNext, 'bCurr', bCurr
!
!  HPrime = 0d0
!  HPlus = 1d0
!  
!  ! Loop through values of n to achieve convergence.
!  tmp = regcut
!  do while ( abs(tmp) >= regcut .and. n < nMax )
!     bNext = -( aPrev*bPrev + aCurr*bCurr ) / aNext
!
!     n = n+1
!     aPrev = aPrev + c*(2*n-m-1)*h    ! Add the difference
!     aCurr = aCurr + 2*n*h - 4*c*hd   !  instead of re-calculating
!     aNext = aNext - 4*hd             !  to hopefully save time.
!     bPrev = bCurr
!     bCurr = bNext
!!     if (debug) write(*,fmt) n, 'aPrev', aPrev, 'aCurr', aCurr, 'aNext', aNext, 'bCurr', bCurr
!
!     tmp = bCurr * (2d0/t)**n
!     HPrime = HPrime - n*tmp   ! Really t*u`(t) and
!     HPlus = HPlus + tmp       !  and u(t) at this point.
!!     if (debug) write(*,fmt) n, 'H+',HPlus, 'H+`',HPrime, 'H-',conjg(HPlus), 'H-`',conjg(HPrime)
!
!  end do
!
!!  if ( n == nMax ) stop ' ! WARNING !   Coulomb (spheroidal_wave.f90) did not converge.'
!
!  HPrime = ( -h/t + hd )/t*HPlus + h/t/t*HPrime
!  HPlus = h/t*HPlus
!!  if (debug) write(*,fmt) n, 'H+', HPlus, 'H+`', HPrime, 'H-', conjg(HPlus), 'H-`', conjg(HPrime)
!
!  HPrime = HPrime * c   ! Switch to coordinate xi.
!!  tmp = ( XIn*HPrime - dXIn*HPlus )   ! Equivalent to the Wronskian of X and F+.
!!  tmp = log(-tmp/conjg(tmp)) / (0d0,2d0)
!!  beta = real(tmp)   ! beta = delta + c - gamma*ln(2) - pi/2 = Delta + c - (l+1)pi/2
!  beta = atan2( -XIn*real(HPrime) + dXIn*real(HPlus) , XIn*imag(HPrime) - dXIn*imag(HPlus) )
!
!  ! Determine the normalisation constant.
!!  tmp = exp(cmplx( 0d0 , beta )) * HPlus
!!  norma = XIn / ( 2d0*real(tmp) )
!  norma = XIn / 2d0 / ( cos(beta)*real(HPlus) + sin(beta)*imag(HPlus) )
!  if ( norma < 0d0 ) then   ! Restrict to positive to prevent arbitrary inclusion of e^(i*pi)=-1.
!     beta = beta + pi
!     norma = -norma
!  endif
!
!  ! Determine the phase shift.
!!  phase = beta - Z/k*log(k*R/2d0) + pi/2d0   ! Rankin and Thorson convention.
!!  phase = beta - Z/k*log(k*R) + dble(l+1)*pi/2d0   ! Ponomarev and Somov convention.
!  phase = -beta - c + pi/2d0 + deltac(-Z/k,l)   ! Ponomarev and Somov convention.
!  phase = modulo( phase+pi, 2d0*pi ) - pi   ! Range (-pi,pi].
!!  if (debug) write(*,'(/3(A,F12.8,9X))') 'norma =', norma, 'phase =', phase, 'beta =', beta
!
!!stop!
!end subroutine Coulomb_t

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

real*8 function appf2(ln,w,x,acc)
  !
  !  Riccati-Bessel function for small rho. Same as spherical Bessel * rho
  !
  !  Original code from numerov.f updated to F90 by J.S.
  !
  implicit none
  
  integer, intent(in) :: ln
  real*8, intent(in) :: w, x
  real*8, intent(inout) :: acc

  integer :: kmax, i, k, iprod, iterm
  real*8 :: rho, ff, sum, summ, sump, zo2k, sumold

  if (w.ne.0.0) then
     !  We use the expansion for the Riccati-Bessel function j(l,rho)
     !  j(l,rho) = rho^(l+1) sum(k=0,oo) (-1)^k 2^k (rho/2)^(2k)/k!/(2(k+l)+1)!!
     kmax = 100
     rho=w*x
     if (rho.gt.1e2) then
        appf2 = 0.0
        return
     endif
     ff = 1d0
     do i=1,ln
        ff = ff * rho / float(2*i+1)
     end do
     sum = ff
     summ = 0d0
     sump = ff
     zo2k = 1d0
     k = 1
     sumold = 0d0
     do while (k.lt.kmax.and.abs(sumold/sum-1d0).gt.1d-6)
        zo2k = zo2k  * rho * rho / 2d0 / float(k) / float(2*(ln+k)+1)
        sumold = sum
        if (mod(k,2).eq.0) then
           sump = sump + zo2k * ff
        else
           summ = summ + zo2k * ff
        endif
        sum = sump - summ
        k = k + 1
        if (abs(summ/sump-1d0).lt.1d-12) then
           k = kmax
           appf2 = 0.0
           print'("Precision loss in APPF2, returning 0.0",1pe14.4)',sum
           return
        endif
     enddo
     acc = abs(summ/sump-1d0)
     if (k.eq.kmax) print'("Possible precision loss in APPF2;result and error",1p,2e14.4)',sum,acc
     appf2 = rho * sum
  else
     !  The following is the limiting case of the Riccati-Bessel/w**(ln+1) for w=0.0
     iprod = 2 * ln + 1
     iterm = 2 * ln + 1
     do while (iterm.gt.3)
        iterm = iterm - 2
        iprod = iprod * iterm
     enddo
     appf2=x**(ln+1)/float(iprod)
  end if
  return
end function appf2

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
!!$
!!$subroutine spheroidal_plane_waves( R,cIn, mIn,lIn, xi, f,fd, g,gd )
!!$  !
!!$  ! Method of van Buren, A.L. & J.E. Boisvert (2002) Quart. App. Math. 60: 589.
!!$  !
!!$  ! Generated in the form Xi(xi) but transformed to X(rho).
!!$  !
!!$  implicit none
!!$
!!$  integer, intent(in) :: mIn,lIn
!!$  real*8, intent(in) :: R,cIn, xi
!!$  real*8, intent(out) :: f,fd, g,gd
!!$
!!$  integer :: narg, lnum, m, ioprad,iopang,iopnorm
!!$  real*16 :: c, x1, tmp
!!$  real*16, dimension(:), allocatable :: r1,r1d, r2,r2d, arg
!!$  real*16, dimension(:,:), allocatable :: s1,s1d
!!$
!!$
!!$  ! c : desired value of the size parameter (= kd/2, where k = wavenumber and d = interfocal length) (real*16)
!!$  c = cIn
!!$
!!$  ! m : desired value for m (integer)
!!$  m = mIn
!!$
!!$  ! lnum : number of values desired for l (integer)
!!$  lnum = lIn - mIn + 1
!!$
!!$  ! ioprad : (integer)
!!$  !        : =0 if radial functions are not computed
!!$  !        : =1 if radial functions of the first kind and their first derivatives are computed
!!$  !        : =2 if radial functions of both kinds and their first derivatives are computed
!!$  ioprad = 2
!!$
!!$  ! x1 : value of the radial coordinate x minus one (real*16) (a nominal value of 10.q0 can be entered for x1 if ioprad = 0) (real*16)
!!$  x1 = xi
!!$
!!$  ! r1,r1d : vectors of lnum calculated radial functions of the first kind and their first derivatives with respect to x (real*16)
!!$  allocate( r1(lnum), r1d(lnum) )
!!$
!!$  ! r2,r2d : vectors of lnum calculated radial functions of the second kind and their first derivatives with respect to x (real*16)
!!$  allocate( r2(lnum), r2d(lnum) )
!!$
!!$  ! iopang : (integer)
!!$  !        : =0 if angular functions are not computed
!!$  !        : =1 if angular functions of the first kind are computed
!!$  !        : =2 if angular functions of the first kind and their first derivatives are computed
!!$  iopang = 0
!!$
!!$  ! iopnorm : (integer)
!!$  !         : =0 if not scaled. The angular functions have the same norm as the corresponding associated legendre function [i.e., we use the Meixner-Schafke normalization scheme.] This norm becomes very large as m becomes very large.
!!$  !         : =1 if angular functions of the first kind (and their first derivatives if computed) are scaled by the square root of the normalization of the corresponding associated Legendre function. The resulting scaled angular functions have unity norm.
!!$  iopnorm = 1
!!$
!!$  ! narg : number of values of the angular coordinate eta for which angular functions are calculated (integer)
!!$  narg = 0
!!$
!!$  ! arg : vector containing the values of eta for which angular functions are desired (real*16)
!!$  allocate( arg(narg) )
!!$
!!$  ! s1,s1d : two-dimensional arrays s1(lnum,narg) and s1d(lnum,narg) that contain narg calculated angular functions and their first derivatives for each of the lnum values of l (real*16). For example, s1(10,1) is the angular function for l = m +10 -1 and the first value of eta given by arg(1)
!!$  allocate( s1(lnum,narg), s1d(lnum,narg) )
!!$
!!$
!!$  call profcn( c, m,lnum, ioprad,x1-1q0, r1,r1d, r2,r2d, iopang,iopnorm, narg,arg, s1,s1d )
!!$
!!$  f = r1(lnum)   ! Function Xi(rho) = Xi(xi).
!!$  g = r2(lnum)
!!$  fd = r1d(lnum) * 2d0/R   ! Derivative w.r.t. rho instead of xi = 2/R*rho + 1
!!$  gd = r2d(lnum) * 2d0/R   !   so Xi'(rho) = d(xi)/d(rho) Xi'(xi) = 2/R*Xi'(xi).
!!$
!!$  if ( m /= 0 ) then
!!$     tmp = ( (xi-1d0)/(xi+1d0) )**(dble(m)/2d0)
!!$
!!$     f = f / tmp   ! Changing to function X(rho)
!!$     g = g / tmp   !   from function Xi(rho) = tmp * X(rho).
!!$     fd = fd / tmp - dble(m) / (xi*xi-1d0) * f   ! Derivative d(X)/d(rho)
!!$     gd = gd / tmp - dble(m) / (xi*xi-1d0) * g   !   from d(Xi)/d(rho).
!!$  end if
!!$
!!$print*, 'f =', f
!!$print*, 'fd =', fd
!!$print*, 'g =', g
!!$print*, 'gd =', gd
!!$
!!$
!!$end subroutine spheroidal_plane_waves
!!$
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!


end module spheroidal_waves_module

