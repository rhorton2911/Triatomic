!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program: H3Plus
!Purpose: performs a structure calculation for the 
!         H3++ molecule. First step towards H3+ structure.
!Author: Reese Horton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program H3Plus
   use grid_radial
   use sturmian_class
   use input_data !Defines global data_in variable
   use basismodule
   use numbers
	 use ieee_arithmetic
   implicit none
   
   type(smallinput):: indata
   integer:: nr  !number of radial grid points
   integer:: ii, jj, kk
   real(dp), dimension(:), allocatable:: w
   real(dp), dimension(:,:), allocatable:: z, wf !, B, H, KMat, 
   complex(dp), dimension(:,:), allocatable:: V, H, KMat, B, VPot, VPotTemp
	 logical, dimension(:), allocatable:: use_list
   complex(dp), dimension(:,:,:), allocatable:: VRadMatEl
   !real(dp), dimension(:,:), allocatable:: VReal
   !complex(dp), dimension(:,:), allocatable:: VTemp
   real(dp), dimension(:,:), allocatable:: realH, realB
   real(dp), dimension(:,:,:), allocatable:: angular
   integer, dimension(:), allocatable:: k_list, l_list, m_list, sturm_ind_list 
	 real(dp), dimension(:), allocatable::  energies
   integer:: num_func, l, m, k, rad_func
   integer:: num_lambda
   integer:: lblocksize
   integer:: N, m_max !,lmax
   integer:: nstates, ier
   real(dp):: alphal, mu
   real(dp):: Rij, cosij
   logical:: sorted, found
   real(dp):: E1, E2, temp
   !Variables used to track program runtime
   character(len=8):: date1, date2
   character(len=10):: time1, time2
   character(len=5):: zone1, zone2
   integer, dimension(8):: values1, values2
   integer:: hr1, hr2, sec1, sec2, mins1, mins2
   integer:: hrs, mins, secs
   !Sturmian data types
   type(basis_sturmian_nr)::basis
   integer:: i1, i2
	 real(dp):: largest, smallest, largestZ
	 integer:: si, sj
	 !Data required for calling lapack routine dsygv
	 real(dp), dimension(:), allocatable:: work
	 integer:: lda, ldb
	 integer:: lwork, info
	 logical:: uselapack
	 !Arrays for testing case
	 integer, dimension(20):: testm, testl, testk
	 integer:: numfound, u1, u2
	 complex(dp), dimension(:,:), allocatable:: tempK, tempB, KMatConroy

   !Intrinsic date and time subroutine
   call date_and_time(date1,time1,zone1,values1)

   !Legacy .f77 subroutine. Sets up common (global :( ) variables used by function YLM in plql.f
   call FAKRED
   call DFSET
 
   !Read in MCCC input data file to set up necessary data_in structure
   call readin( data_in, 10, 1)
   !Initialise data type containing rgrid and integration weights, global variable ( ): ) in grid_radial module
   call setgrids(grid)

   call readInput(indata)
   N = indata%N
   m_max = indata%mmax
	 uselapack = .true.

   !Initialise sturmian data types, requires input data type defined in the MCCC code.
	 if (indata%isoscop .eq. 0) then
      call construct(basis, data_in)
	 else if (indata%isoscop .eq. 1) then
			if (indata%conroybasis .eq. 1) then
				 !Just for testing, remove
				 indata%R(1) = 1.0
				 indata%R(2) = 1.0
				 indata%R(3) = 0.0
				 indata%charge(1) = 1
				 indata%charge(2) = 1
				 indata%charge(3) = 0

			   call construct_all_conroy_basis(basis, grid, indata, data_in, KMatConroy)
			   !Write some basis functions to file for testing
			   open(70,file='conroyfunc.txt')
			   write(70,*) "r(a_0)  psi_(k,0,0)(k=1->5)"
			   do ii = 1, grid%nr
			   	 write(70,*) grid%gridr(ii), (basis%b(jj)%f(ii), jj=1, 5)
			   end do
			   close(70)
				 stop
			else
         call construct(basis, data_in)
				 !print*, indata%R(1), indata%R(2), indata%R(3)
				 !print*, indata%theta(1), indata%theta(2), indata%theta(3)
				 !print*, indata%phi(1), indata%phi(2), indata%phi(3)
				 !print*, indata%charge(1), indata%charge(2), indata%charge(3)
	    end if
	 end if

   nr = grid%nr

   !Set up expansion of potential, use formula: SUM_l=0^l=L (2l+1) = (L+1)^2
   num_lambda = (indata%lambdamax+1)**2

   allocate(VPot(nr,num_lambda))
   VPot(:,:) = 0.0_dp
   allocate(VPotTemp(nr, num_lambda))
   do ii = 1, 3
      call getVPotNuc(grid, VPotTemp, indata%R(ii), indata%theta(ii), &
			indata%phi(ii), indata%charge(ii), indata)
      VPot(:,:) = VPot(:,:) + VPotTemp(:,:)
   end do
   deallocate(VPotTemp)
   
   !Triatomic molecules do not have good quantum numbers such as angular momentum and parity (like H2+).
   !Instead, the clamped nuclei electronic V-matrix is diagonal in different irreducible representations of a 
   !given point group. Since we are not currently using a symmetry adapted basis, we can only loop over l
   num_func=0
   do l = data_in%labot, data_in%latop
      !Now need to loop over l and m, m is no longer conserved.
      do m = -l, l
         num_func = num_func + data_in%nps(l) !Array nps indexed from zero
      end do
   end do
   if (num_func .eq. 0) then
      print*, "ERROR: num_func=0: no basis functions. Stopping"
      stop
   end if


	 !Defined arrays for test case
	 testm = (/0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, 2, 2, 2, -3, -3/)
	 testl = (/0, 0, 0, 0, 0, 2, 2, 2, 4, 1, 1, 1, 1, 3, 3, 2, 2, 2, 3, 3/)
	 testk = (/1, 2, 3, 4, 5, 3, 4, 5, 5, 2, 3, 4, 5, 4, 5, 3, 4, 5, 4, 5/)	

   !Define arrays to keep track of k,l,m values for each index
   !Indexing scheme: specify (l,m), then k goes from 1 to N for each such pair 
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
 
   !Calculate the matrix elements, iterating through l values
   allocate(H(num_func,num_func),B(num_func,num_func),V(num_func,num_func))
   allocate(KMat(num_func, num_func))
   KMat(:,:) = 0.0_dp
   H(:,:) = 0.0_dp
   V(:,:) = 0.0_dp
   B(:,:) = 0.0_dp

   !Define B as before, but now need to ensure it is block diagonal in m also
   B(:,:) = 0.0_dp
   do ii = 1, num_func-1
         B(ii,ii) = 1.0_dp
         l=l_list(ii)
	       m=m_list(ii)
         k=k_list(ii)
         if ((l_list(ii+1) .ne. l) .or. (m_list(ii+1) .ne. m)) then
            cycle
         end if
         !Formula for upper diagonal B matrix elements, B is symmetric
         B(ii,ii+1) = (-0.5_dp)*sqrt(1.0_dp-dble(l*(l+1))/dble((k+1+l)*(k+l)))
         B(ii+1,ii) = B(ii,ii+1)
   end do
   B(num_func,num_func)= 1.0_dp

   !K-matrix elements: analytical formula from slides
   !K matrix is block diagonal in m and in l
   !l specifies a block and m specifies a particular sub-block
   KMat(:,:) = 0.0_dp
   do ii = 1, num_func
      l = l_list(ii)
      alphal=data_in%alpha(l)
      KMat(ii,ii) = (alphal**2)
   end do

   l = -1
   jj=0
   kk=0
   lblocksize=0
   do ii = 1, num_func
      if (l .eq. l_list(ii)) then
      	 cycle
      end if
      l = l_list(ii)
      alphal = data_in%alpha(l)
      !Assuming we choose the same range of k for each m at a given l
      lblocksize = data_in%nps(l)*(2*l+1)
      jj=kk+1
      kk = kk + lblocksize
      KMat(jj:kk,jj:kk) = KMat(jj:kk,jj:kk) - 0.5_dp*(alphal**2)*B(jj:kk,jj:kk)
   end do
   mu = 1.0_dp
   KMat(:,:) = KMat(:,:)/mu

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
				 stop
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
			H(:,:) = 0.0_dp
			V(:,:) = 0.0_dp
	 end if
 
   print*, "GET ANGULAR MATRIX ELEMENTS"
   !!Precalculate angular integrals appearing in V-matrix elements
   allocate(angular(num_lambda,num_func,num_func))
   call getAngular(num_func, angular, l_list, m_list, indata)

   !Precalculate radial matrix elements
   allocate(VRadMatEl(num_lambda,rad_func,rad_func))
   VRadMatEl(:,:,:) = 0.0_dp
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

	 largest = 0.0_dp
	 do ii = 1, num_func
	    do jj = 1, num_func
	       if (isnan(realH(ii,jj))) then
						print*, "realH: ", ii,jj
						print*, realH(ii,jj)
						stop
				 end if
	       if (isnan(realB(ii,jj))) then
						print*, "realB: ", ii,jj
						stop
				 end if
	       if (.not. ieee_is_finite(realH(ii,jj))) then
						print*, "realH infinite: ", ii,jj
						print*, H(ii,jj)
						print*, realH(ii,jj)
						stop
				 end if
	       if (.not. ieee_is_finite(realB(ii,jj))) then
						print*, "realB infinite: ", ii,jj
						stop
				 end if

				 if ((abs(realH(ii,jj)) .gt. largest) .and. (abs(realH(ii,jj)) .lt. 0.0_dp)) then
				    largest = realH(ii,jj)
				 end if
			end do
	 end do
	 do ii = 1, num_func
	    do jj = 1, num_func
	       if (abs(realH(ii,jj)) .lt. 1E-8*abs(largest)) then
						realH(ii,jj) = 0.0_dp
				 end if
			end do
	 end do
	 print*, "Removed small matrix elements"

   !Call the rsg subroutine to solve the eigenvalue problem for the energies
   allocate(w(num_func),z(num_func,num_func))

   !ALTERNATIVE: USE LAPACK ZHEGVX FUNCTION FOR HERMITIAN MATRICES RATHER THAN REAL SYMMETRIC
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
			deallocate(work)
	 end if

	 !Find largest expansion coefficient, ignore those below a certain
	 !magnitude for stability/to get rid of underflow 
	 largestZ = 0.0_dp
	 do ii = 1, num_func
	    do jj = 1, num_func
	       if (abs(z(ii,jj)) .gt. abs(largestZ)) then
						largestZ = z(ii,jj)
				 end if
			end do
	 end do
	 do ii = 1, num_func
	    do jj = 1, num_func
	       if ((abs(z(ii,jj)) .lt. 1E-10*abs(largestZ)) .and. (abs(z(ii,jj)) .gt. 0.0_dp)) then
						z(ii,jj) = 0.0_dp
				 end if
			end do
	 end do

	 !smallest = 1000.0_dp 
	 !do ii = 1, num_func
	 !   do jj = 1, num_func
	 !      if ((abs(z(ii,jj)) .lt. abs(smallest)) .and. (abs(z(ii,jj)) .gt. 0.0_dp)) then 
	 ! 		  	smallest = z(ii,jj)
	 ! 		  	si = ii
	 ! 		  	sj = jj
	 ! 		 end if
	 ! 		 !if ((l_list(ii) .eq. 2) .and. ((m_list(ii) .eq. -1) .and. (k_list(ii) .eq. 1))) then
	 ! 		 !      if (jj .eq. 159) then
	 ! 		 ! 		    print*, "N=10, smallest: ", z(ii,159)
	 ! 		 ! 		    print*, l_list(ii), m_list(ii), k_list(ii)
	 ! 		 ! 		 end if
	 ! 		 !end if
	 ! 	end do
   !end do
	 !print*, smallest, si, sj
	 !!!!!print*, z(51,159)
	 !print*, l_list(si), m_list(si), k_list(si)
	 !print*, l_list(sj), m_list(sj), k_list(sj)			 
	 
!	 info = 0
!	 lwork = 8*num_func + 1
!	 allocate(alphar(num_func), alphai(num_func), beta(num_func))
!	 allocate(vl(num_func,num_func), vr(num_func,num_func), work(lwork))
!	 alphar(:) = 0.0_dp
!	 alphai(:) = 0.0_dp
!	 beta(:) = 0.0_dp
!	 vl(:,:) = 0.0_dp
!	 vr(:,:) = 0.0_dp
!	 jobvl = 'V'
!	 jobvr = 'V'
!	 work(:) = 0.0_dp
!	 !First call is a workspace query to determine the optimal size of
!	 !'work'.
!	 print*, "BEFORE FIRST"
!	 call sggev(jobvl, jobvr, num_func, realH, num_func, realB, num_func, alphar, &
!	            alphai, beta, vl, num_func, vr, num_func, work, -1, info)
!	 print*, "AFTER FIRST"
!   lwork = int(work(1))
!	 print*, info, lwork
!   !deallocate(work)
!	 !allocate(work(lwork))
!	 !!Actual calculation
!	 call sggev(jobvl, jobvr, num_func, realH, num_func, realB, num_func, alphar, &
!	            alphai, beta, vl, num_func, vr, num_func, work, lwork, info)
!
!	 print*, "COMPUTED"
!   do ii = 1, 1
!	  	print*, alphar(ii), beta(ii)
!	  	!if (beta(ii) .gt. 0.0_dp) then
!	  	!   print*, alphar(ii)/beta(ii)
!	  	!else
!	  	!	 print*, alphar(ii)
!	  	!end if
!	 end do
!	 stop
!
!	 if (info .eq. 0) then
!	  	do ii = 1, num_func
!	  		 if (abs(beta(ii)) .gt. 0.0_dp) then
!	          w(ii) = alphar(ii)/beta(ii)
!	  				z(:,ii) = vr(:,ii) !Copy eigenvectors to existing array
!	  		 else 
!	  				print*, "ERROR: beta(ii) equal to zero. Stopping."
!	  				stop
!	  		 end if
!	  	end do
!	 else
!	    print*, "ERROR: solving for eigenvalues unsuccessful. Stopping."
!	  	stop
!   end if
!	 deallocate(alphar, alphai, beta, vl, vr, work)


   deallocate(realH,realB)

   !Include effect of nucleus-nucleus interaction term by simply adding zi*zj/R_ij to energies at the end.
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
         if (Rij .gt. 0.0_dp) then
            w(:) = w(:) + dble(indata%charge(ii)*indata%charge(jj))/Rij
         end if
      end do
   end do
   print*, "NUCLEAR INTERACTION CALCULATED"

   open(80,file="energies.txt") 
   do ii = 1, num_func
      !print*, ii, w(ii)
      write(80,*) ii, w(ii)
   end do
   close(80)

   !Need to divide by r and multiply by Ylm or Xlm to get full 3D wave function
   print*, "CALCULATED WAVE FUNCTIONS"
   allocate(wf(nr,num_func))
   wf(:,:) = 0.0_dp
   do ii = 1, num_func
    	do jj=1, num_func
	       if (abs(z(ii,jj)) .lt. 1E-20 ) then
	          z(ii,jj) = 0.0_dp
	       end if
         kk = sturm_ind_list(jj)
         !MODIFY TO BE A FUNCTION OF (r,theta,phi)

	       !Basis function set to zero outside of range of r values indexed by i1, i2
	       i1 = basis%b(kk)%minf
	       i2 = basis%b(kk)%maxf
               
    	   wf(i1:i2,ii) = wf(i1:i2,ii) +  z(jj,ii)*basis%b(kk)%f(i1:i2)!*YLM
    	end do
   end do
   print*, "WAVE FUNCTIONS CALCULATED"


   !Number of energies calculated by rsg will be equal to the size of the basis used.
   !Allocate array to store energy of each state
   nstates = num_func
   allocate(energies(nstates))
   energies(:) = 0.0_dp
   nstates=0
   nstates=num_func
   energies(:) = w(:)

   deallocate(k_list, l_list, m_list)
   deallocate(H,B,V,KMat)
	 deallocate(use_list)
   deallocate(w,z,wf)
	 if (indata%isoscop .eq. 1) then
				 if (indata%conroybasis .eq. 1) then
				    deallocate(KMatConroy)
				 end if
	 end if

   !Print out the first energy in each symmetry
   !ii= 1
   !do m=0, m_max
   !   do par=-1,1,2 
   !   	 print*, "M= ", m, "par= ", par, ": ", energies(ii), " Ha"
   !   	 ii = ii + 1 
   !   end do
   !end do

   print*, "SORT ENERGIES"
   !Sort array of lowest energies for each symmetry, use bubble sort algorithm 
   !for simplicity
   sorted = .false.
   do while( .not. sorted)
      sorted = .true.
      do ii = 1, nstates-1
         E1 = energies(ii)
         E2 = energies(ii+1)
         
         if (E2 < E1) then
            sorted = .false.
            temp = energies(ii)
            energies(ii) = energies(ii+1)
            energies(ii+1) = temp
         end if
      end do
   end do

   print*, "WRITE ENERGIES TO FILE"
   !Write energes of states to file
   open(80,file="energies.txt")
   write(80,*) "State Energies (Ha)"
   write(80,*) "R1=", indata%R(1), " R2=", indata%R(2), "R3=", indata%R(3)
   write(80,*) "N (index), E(Ha)"
   do ii = 1, nstates
      write(80,*) ii, energies(ii)
   end do
   close(80)
   
   call deconstructInput(indata)

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
