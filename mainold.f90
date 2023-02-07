!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program: H3Plus
!Purpose: performs a structure calculations for one and two electron
!         H3 molecule.
!Author: Reese Horton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!GLOBAL VARIABLES.... GLOBAL VARIABLES.... GLOBAL VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------->>>>>GLOBAL VARIABLES<<<<<-----------------
!Due to the MCCC code inheriting features and/or programming
!constructs from legacy fortran77 code, there are several 
!instances of globally defined data structures being used.
!Most often, these are declared in the preamble of a particular
!module, which is then USED in the subroutine that makes use of the
!data structure. The following is a list of such variables directly used 
!in the present version of the MCCC code. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!dpf:: format specifier for declaring double precision floats in modern fortran
!   - located in numbers.f90
!data_in:: contains program input data read in from data.in file.
!   - located in input_data.f90
!   - initialisation subroutine readin uses global variables defined in
!     target_data.f90 and input_parser.f90.
!grid:: instance of rgrid type, stores grid of radial points and integration weights, etc
!   - located in gridr.f90
program H3Plus
   use grid_radial
   use sturmian_class
   use state_class
   use input_data !Defines global data_in variable
   use basismodule
   use numbers
   use ieee_arithmetic
   implicit none
   
   type(smallinput):: indata
   integer:: num_lambda
   integer:: lblocksize
   integer:: nstates, ier
   real(dpf):: alphal, mu
   real(dpf):: Rij, cosij
   logical:: sorted, found
   real(dpf):: E1, E2, temp
   !Variables used to track program runtime
   character(len=8):: date1, date2
   character(len=10):: time1, time2
   character(len=5):: zone1, zone2
   integer, dimension(8):: values1, values2
   integer:: hr1, hr2, sec1, sec2, mins1, mins2
   integer:: hrs, mins, secs
   !Sturmian data types: original basis, auxialliary basis summed over k
   type(basis_sturmian_nr)::basis
   integer:: i1, i2
   real(dpf):: largest, smallest, largestZ
   integer:: si, sj

		!Nuclear coordinates and additional options
    type(smallinput):: indata
		!Given sturmian basis
    type(basis_sturmian_nr)::basis
    !State basis type for storing one electron states
    type(basis_state):: oneestatebasis
    integer:: nr  
		!Matrix elements, energies and CI coeficients to be calculated.
    integer:: ii, jj, kk
    real(dpf), dimension(:), allocatable:: w
    real(dpf), dimension(:,:), allocatable:: z
    complex(dpf), dimension(:,:), allocatable:: V, H, KMat, B, VPot, VPotTemp
    logical, dimension(:), allocatable:: use_list
    complex(dpf), dimension(:,:,:), allocatable:: VRadMatEl
	  !Basis function index arrays and matrix elements
    real(dpf), dimension(:,:), allocatable:: realH, realB, realK, BMat
    real(dpf), dimension(:,:,:), allocatable:: angular
    integer, dimension(:), allocatable:: k_list, l_list, m_list, sturm_ind_list 
    real(dpf), dimension(:), allocatable::  energies
    integer:: num_func, l, m, k, rad_func
    !Arrays for testing case
    integer, dimension(20):: testm, testl, testk
    integer:: numfound, u1, u2
    complex(dpf), dimension(:,:), allocatable:: tempK, tempB
	  !For renormalising CI coeffs
    integer:: n
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

   !State basis type for storing one electron states
   type(basis_state):: oneestatebasis
   integer:: inum, ncm, counter
   integer, dimension(:), allocatable:: no1, no2, mo1, mo2, phase

   !Intrinsic date and time subroutine
   call date_and_time(date1,time1,zone1,values1)

   !Legacy .f77 subroutines. Sets up common (global :( ) variables used by function YLM in plql.f
   call FAKRED
   call DFSET
 
   !Read in MCCC input data file to set up necessary data_in structure
   call readin( data_in, 10, 1)
   !Initialise data type containing rgrid and integration weights, global variable ( ): ) in grid_radial module
   call setgrids(grid)

   call readInput(indata)
   uselapack = .true.

   !Initialise sturmian data types, requires input data type defined in the MCCC code.
   if (indata%isoscop .eq. 0) then
      !call construct(basis, data_in)
      call construct_all_nr_m(basis,data_in)
   else if (indata%isoscop .eq. 1) then
      !if (indata%conroybasis .eq. 1) then
      !	 !Just for testing, remove
      !	 indata%R(1) = 1.0
      !	 indata%R(2) = 1.0
      !	 indata%R(3) = 0.0
      !	 indata%charge(1) = 1
      !	 indata%charge(2) = 1
      !	 indata%charge(3) = 0
      !
      !   call construct_all_conroy_basis(basis, grid, indata, data_in, KMatConroy)
      !   !Write some basis functions to file for testing
      !   open(70,file='conroyfunc.txt')
      !   write(70,*) "r(a_0)  psi_(k,0,0)(k=1->5)"
      !   do ii = 1, grid%nr
      !   	 write(70,*) grid%gridr(ii), (basis%b(jj)%f(ii), jj=1, 5)
      !   end do
      !   close(70)
      !else
         call construct(basis, data_in)
         !print*, indata%R(1), indata%R(2), indata%R(3)
         !print*, indata%theta(1), indata%theta(2), indata%theta(3)
         !print*, indata%phi(1), indata%phi(2), indata%phi(3)
         !print*, indata%charge(1), indata%charge(2), indata%charge(3)
      !end if
   end if


   nr = grid%nr

   !Set up expansion of potential, use formula: SUM_l=0^l=L (2l+1) = (L+1)^2
   num_lambda = (indata%lambdamax+1)**2

   allocate(VPot(nr,num_lambda))
   VPot(:,:) = 0.0_dpf
   allocate(VPotTemp(nr, num_lambda))
   do ii = 1, 3
      call getVPotNuc(grid, VPotTemp, indata%R(ii), indata%theta(ii), &
	              indata%phi(ii), indata%charge(ii), indata)
      VPot(:,:) = VPot(:,:) + VPotTemp(:,:)
   end do
   deallocate(VPotTemp)
   
   !Triatomic molecules do not have good quantum numbers such as angular momentum and parity (like H2+).
   !Instead, the clamped nuclei electronic V-matrix is diagonal in different irreducible representations of a 
   !given point group. Since we are not currently using a symmetry adapted basis, we can only loop over lm
   num_func=0
   do l = data_in%labot, data_in%latop
      !Now need to loop over l and m, m is no longer conserved.
      do m = -l, l
         num_func = num_func + data_in%nps(l) !Array nps indexed from zero
      end do
   end do
   if (num_func .eq. 0) then
      print*, "ERROR: num_func=0: no basis functions. Stopping"
      error stop
   end if

   !Defined arrays for test case
   testm(:) = indata%testm(:) 
   testl(:) = indata%testl(:)
   testk(:) = indata%testk(:) 

   !Define arrays to keep track of k,l,m values for each index
   !Indexing scheme: specify (l,m), then k goes from 1 to N_l for each such pair 
   !i.e l -> m -> k_lm
   allocate(k_list(num_func), l_list(num_func), m_list(num_func), sturm_ind_list(num_func))
   ii = 0
   rad_func = 0
   kk = 0
   do l = data_in%labot, data_in%latop
      do m = -l, l
         jj  = kk 
         do k = 1, data_in%nps(l)
            ii = ii + 1
            k_list(ii) = k
            l_list(ii) = l
	          m_list(ii) = m
	          !print*, "L,M,k,ii: ", l, m, k, ii

	          !Indexing changed for non-diatomics, indices for basis now include l, k AND m. For a given l, the same radial basis
	          !functions phi_{kl} are used for each value of m from -l to l. Introduce an array to track which radial basis 
	          !function belongs to which index.
	          jj = jj + 1

	          found = .false.
	          n = 1
	          do while (.not. found)
               !print*, basis%b(n)%l, basis%b(n)%m, basis%b(n)%k, "L,K: ", l, k
	             if ((basis%b(n)%l .eq. l) .and. (basis%b(n)%k .eq. k)) then
                        found = .true.
	             else
                        n = n +1
	             end if
	          end do

	          sturm_ind_list(ii) = n
	          !print*, jj
         end do
      end do
      !The radial part of the basis is identical for any given (l,k) pair, regardless of m
      do k = 1, data_in%nps(l)
         rad_func = rad_func + 1
      end do

      kk = kk + data_in%nps(l)
   end do
 
   allocate(H(num_func,num_func),V(num_func,num_func))
   H(:,:) = 0.0_dpf
   V(:,:) = 0.0_dpf

   allocate(realB(basis%n,basis%n))
   realB(:,:) = basis%ortint(:,:)
   call getKMatM(realK,realB,basis)
   allocate(B(num_func,num_func), KMat(num_func,num_func))
   B(:,:) = realB(:,:)
   KMat(:,:) = realK(:,:)
   deallocate(realB, realK)

   allocate(use_list(num_func))
   use_list(:) = .true.
   if ((indata%isobasisop .eq. 1) .and. (indata%isoscop .eq. 1)) then
      use_list(:) = .false.
      numfound = 0
      do ii = 1, 20
         print*, testk(ii), testl(ii), testm(ii)
         do jj = 1, num_func
      	    if (((k_list(jj) .eq. testk(ii)) .and. &
      	    (l_list(jj) .eq. testl(ii)))  .and. &
            (m_list(jj) .eq. testm(ii)))  then
      
      	       use_list(jj) = .true.
      	       numfound = numfound + 1
            end if
      	 end do
      end do
      if (numfound .ne. 20) then
      	 print*, "ERROR: basis functions with required symmetries for &
      	 &test case not found, require at least lmax=3 and N=5. &
      	 & Stopping."
      	 error stop
      end if

      !Produce K and B matrix for restricted basis, just copy over the
      !right elements from the full K and B matrices
      allocate(tempK(num_func,num_func), tempB(num_func, num_func))
      tempK(:,:) = KMat(:,:)
      tempB(:,:) = B(:,:)
      deallocate(KMat, B)
      allocate(KMat(20,20), B(20,20))
      u1 = 1
      do ii = 1, num_func
         u2 = 1
         do jj= 1, num_func
            if (use_list(ii) .and. use_list(jj)) then
               KMat(u1,u2) = tempK(ii,jj)
               B(u1,u2) = tempB(ii,jj)
               u2 = u2 + 1
            end if
         end do
         if (use_list(ii)) then
            u1 = u1 + 1
         end if
      end do
      deallocate(tempK, tempB, H, V)	
				    
      num_func = 20
      allocate(H(num_func,num_func), V(num_func,num_func))
      H(:,:) = 0.0_dpf
      V(:,:) = 0.0_dpf
   end if
 
   print*, "GET ANGULAR MATRIX ELEMENTS"
   !!Precalculate angular integrals appearing in V-matrix elements
   allocate(angular(num_lambda,num_func,num_func))
   call getAngular(num_func, angular, l_list, m_list, indata)

   !Precalculate radial matrix elements
   allocate(VRadMatEl(num_lambda,basis%n,basis%n))
   VRadMatEl(:,:,:) = 0.0_dpf
   print*, "GET RADIAL MATRIX ELEMENTS"
   call getRadMatEl(basis, VPot, grid, VRadMatEl, indata)

   !If real spherical harmonics are used, V matrix elements will have complex part zero.
   call getVMatEl(sturm_ind_list, V, VRadMatEl, num_func, angular, indata, use_list)
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
	       cosij = cos(indata%theta(ii))*cos(indata%theta(jj)) + &
	               sin(indata%theta(ii))*sin(indata%theta(jj))*cos(indata%phi(ii)-indata%phi(jj))
         Rij = sqrt(indata%R(ii)**2 + indata%R(jj)**2 - &
	       2*indata%R(ii)*indata%R(jj)*cosij)
	       !Account for degenerate case where nuclei coincide
         if (Rij .gt. 0.0_dpf) then
            realH(:,:) = realH(:,:) + realB(:,:)*dble(indata%charge(ii)*indata%charge(jj))/Rij
         end if
      end do
   end do
   print*, "NUCLEAR INTERACTION CALCULATED"

   !print*, "HMatrix"
   !do ii = 1, num_func
   !   print*, basis%b(ii)%k, basis%b(ii)%l, basis%b(ii)%m
   !end do
   !do ii = 1, num_func
   !   write(*,'(10000(ES12.5,X))'), realH(ii,:)
   !end do

   !do ii = 1, num_func
   !   do jj = ii, num_func
   !      if (abs(realH(ii,jj)) .gt. 0.0_dpf) then
   !         print*, realH(ii,jj), "l1,m1 ", basis%b(ii)%l, basis%b(ii)%m, "l2,m2 ", basis%b(jj)%l, basis%b(jj)%m
   !      end if
   !   end do
   !end do

   !do ii = 1, num_func
   !   do jj = ii, num_func
   !      if (abs(realH(ii,jj)) .gt. 0.0_dpf) then
   !         if (basis%b(ii)%m .ne. basis%b(jj)%m) then
   !            print*, "NON-ZERO matrix element"
   !            print*, realH(ii,jj), "l1,m1 ", basis%b(ii)%l, basis%b(ii)%m, "l2,m2 ", basis%b(jj)%l, basis%b(jj)%m
   !         end if
   !      end if
   !   end do
   !end do

   if (indata%harmop .eq. 0) then
      !Molecule has C_s symmetry if:
      !(1)  R_2 != R_3
      !(2)  nuclei 2 and 3 are not symmetric about z-axis
      !(3)  nucleus 1 does not lie on z-axis (should never happen)
      !Matrix elements have non-zero complex part in C_s symmetry
      if (((abs(indata%R(2) - indata%R(3)) .gt. 0.0_dpf) .or. &
         (abs(abs(indata%theta(2)) - abs(indata%theta(3))) .gt. 0.0_dpf)) .or. &
         (abs(indata%theta(1)) + abs(indata%phi(1))) .gt. 0.0_dpf) then
         print*, "WARNING: complex spherical harmonics not yet implemented for C_s geometries. Stopping."	
         error stop
      end if
   end if

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

   !do ii = 1, num_func
   !   print*, "State number: ", ii
   !   do jj = 1, num_func
   !      print*, basis%b(jj)%k, basis%b(jj)%l, basis%b(jj)%m, z(jj,ii)
   !   end do
   !  print*, "   "
   !end do

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
   call new_basis_st(oneestatebasis,num_func, .false. , 0)
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
      call construct_st(oneestatebasis%b(ii),.false.,0.0_dpf,0,0.5_dpf,w(ii),inum,num_func,z(:,ii),no1,mo1,no2,mo2,phase)
      deallocate(no1,no2,mo1,mo2,phase)
   end do

   if((.not. (sum(indata%R(:)) .gt. 0.0_dpf)) .and. (writeWaveFunc)) then
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
         write(80,*), ii, z(ii,1)
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

   deallocate(k_list, l_list, m_list)
   deallocate(H,B,V,KMat)
   deallocate(use_list)
   deallocate(w,z)    !,wf)
   if (indata%isoscop .eq. 1) then
       !if (indata%conroybasis .eq. 1) then
       !   deallocate(KMatConroy)
       !end if
   end if

   print*, "WRITE ENERGIES TO FILE"
   !Write energes of states to file
   open(80,file="1eenergies.txt")
   write(80,*) "State Energies (Ha)"
   write(80,*) "R1=", indata%R(1), " R2=", indata%R(2), "R3=", indata%R(3)
   write(80,*) "N (index), E(Ha)"
   do ii = 1, nstates
      write(80,*) ii, energies(ii)
   end do
   close(80)
















   !---------------------Perform 2e Structure --------------------!
   !one_electron_func.f90: fills oneestates type before structure12
   !is called, specifically does so in construct_1el_basis_nr
   !and Hybrid_MSCbasis functions, which also set e1me and ovlpst
   !arrays
   !Hybrid_MSCbasis actually already loops over m when copying
   !laguerre functions to TargetStates
   !rearrange(): called in construct_1el_basis_nr, I've already
   !modified rearrange, might need to check its usage still works
   !when called in construct_1el_basis_nr
  
   !Call rearrange to represent 1e states in a one-electron basis with well defined (lm)
   call rearrange(basis,data_in%latop,oneestatebasis,.false.)
  
   !Modified version of subroutine in one_electron_func.f90, constructs basis for use in 2e configs
   !call construct_1el_basis_nr_group(  )

   !Custom version of the structure12 subroutine tailored to non-linear molecules
   call structure12group(basis,oneestatebasis,basis%n,indata)

   !--------------------End Two Electron Structure------------------!

   deallocate(energies)
   call destruct_gridr(grid)

   !Compute program runtime
   print*, "COMPUTE TIME ELAPSED"
   call date_and_time(date2,time2,zone2,values2)
   read(time1(1:2), '(I20)') hr1
   read(time2(1:2), '(I20)') hr2
   read(time1(3:4), '(I20)') mins1
   read(time2(3:4), '(I20)') mins2
   read(time1(5:6), '(I20)') sec1
   read(time2(5:6), '(I20)') sec2
   hrs = abs(hr2 - hr1)
   if (hrs .gt. 0) then
      hrs = hrs -1
   end if
   mins = abs(mins2 - mins1)
   if (mins .gt. 0) then
      mins = mins -1
      if (60+sec2-sec1 .ge. 60) then
         secs = sec2 - sec1
	 mins = mins + 1
      else if (60+sec2-sec1 .lt. 60) then
	 secs = 60+sec2-sec1
      end if
   else
      secs = sec2 - sec1
   end if
   print*, "ELAPSED TIME (hrs, mins, secs): ", hrs, mins, secs

end program H3Plus
