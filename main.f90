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
   implicit none
   
   type(smallinput):: indata
   integer:: nr  !number of radial grid points
   integer:: ii, jj, kk
   real(dp), dimension(:), allocatable:: w
   real(dp), dimension(:,:), allocatable:: z, wf !, B, H, KMat, 
   complex(dp), dimension(:,:), allocatable:: V, H, KMat, B, VPot, VPotTemp
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
   real(dp), dimension(3):: R, theta, phi !Coordinates and charge of each nucleus
   integer, dimension(3):: charge
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

   !Intrinsic date and time subroutine
   call date_and_time(date1,time1,zone1,values1)

   !Legacy .f77 subroutine. Sets up common (global :( ) variables used by function YLM in plql.f
   call FAKRED
   call DFSET

   !Read in MCCC input data file to set up necessary data_in structure
   call readin( data_in, 10, 1)
   !Initialise data type containing rgrid and integration weights, global variable ( ): ) in grid_radial module
   call setgrids(grid)
   !Initialise sturmian data types, requires input data type defined in the MCCC code.
   call construct(basis, data_in)
 
   call readInput(indata)
   N = indata%N
   m_max = indata%mmax
   !lmax = indata%l
   !Specify charge of each nucleus
   charge(1) = indata%charge(1)
   charge(2) = indata%charge(2)
   charge(3) = indata%charge(3)
    
   !Specify distances of nuclei from the origin
   R(1) = indata%R(1)
   R(2) = indata%R(2)
   R(3) = indata%R(3)

   !Angular coordinates of each nucleus
   theta(1) = indata%theta(1)
   theta(2) = indata%theta(2)
   theta(3) = indata%theta(3)
   phi(1) = indata%phi(1)
   phi(2) = indata%phi(2)
   phi(3) = indata%phi(3)

   nr = grid%nr

   !Set up expansion of potential, use formula: SUM_l=0^l=L (2l+1) = (L+1)^2
   num_lambda = (indata%lambdamax+1)**2

   allocate(VPot(nr,num_lambda))
   VPot(:,:) = 0.0_dp
   allocate(VPotTemp(nr, num_lambda))
   do ii = 1, 3
      call getVPotNuc(grid, VPotTemp, R(ii), theta(ii), phi(ii), charge(ii), indata)
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

   !Number of energies calculated by rsg will be equal to the size of the basis used.
   !Allocate array to store energy of each state
   nstates = num_func
   allocate(energies(nstates))
   energies(:) = 0.0_dp
   nstates=0

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
   
   !!Construct a basis with the given symmetry properties using the l_list and 
   !!k_list variables
   !allocate(radbasis(nr,rad_func))
   !l=-1
   !kk=0
   !jj=0
   !do ii=1,  num_func
   !   if (l .eq. l_list(ii)) then
   !   	 cycle
   !   end if
   !   l = l_list(ii)
   !   !In principle aliglows to specify a different alpha for each l
   !   alphal = indata%alpha(l+1)
   !   jj=kk+1
   !   kk = kk+N
   !
   !   allocate(basisatl(nr, indata%N))
   !   call createBasis(basisatl, l, alphal, indata%N, grid)
   !   radbasis(:,jj:kk) = basisatl(:,1:N)
   !   deallocate(basisatl)
   !end do

   !!Write some of the basis functions to file for plotting
   !open(70,file="basis.txt")
   !write(70,*) "l, lmax,N: ", 2, lmax, N
   !if (lmax .ge. 2) then
   !   do jj=1, nr
   !      write(70,*) grid%gridr(jj), (radbasis(jj,ii), ii=2*N+1, 3*N)
   !   end do
   !end if
   !close(70)
   
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
      !Assuming we choose the same range of k for each l
      lblocksize = data_in%nps(l)*(2*l+1)
      jj=kk+1
      kk = kk + lblocksize
      KMat(jj:kk,jj:kk) = KMat(jj:kk,jj:kk) - 0.5_dp*(alphal**2)*B(jj:kk,jj:kk)
   end do
   mu = 1.0_dp
   KMat(:,:) = KMat(:,:)/mu

   print*, "GET ANGULAR MATRIX ELEMENTS"
   !!Precalculate angular integrals appearing in V-matrix elements
   allocate(angular(num_lambda,num_func,num_func))
   call getAngular(num_func, angular, l_list, m_list, indata)
   !!Follow same indexing scheme as basis, specify lambda, then v from -lambda to lambda
   !do ii = 1, num_func
   !   do jj = 1, num_func
   !      li = l_list(ii)
   !      mi = m_list(ii)
   !      lj = l_list(jj)
   !      mj = m_list(jj)
   !      lambdaind=1
   !      do lambda = 0, indata%lambdamax
   !         do q = -lambda, lambda
   !            if (indata%harmop .eq. 0) then
   !               angular(lambdaind,ii,jj) = sqrt(dble(2*lambda+1)/(4.0_dp*pi)) &
   !     	     *Yint(dble(li),dble(mi),dble(lambda),dble(q),dble(lj),dble(mj))
   !            else if (indata%harmop .eq. 1) then
   !     	  !Xint as written calculates overlap without sqrt(2lambda+1) factor, unlike Yint
   !               angular(lambdaind,ii,jj) = Xint(dble(li),dble(mi),dble(lambda),dble(q),dble(lj),dble(mj))
   !            end if

   !            !print*, li, mi, lambda, q, lj, mj
   !            !print*, angular(lambdaind,ii,jj)
   !            lambdaind = lambdaind + 1 
   !         end do
   !      end do
   !   end do
   !end do

   !Precalculate radial matrix elements
   allocate(VRadMatEl(num_lambda,rad_func,rad_func))
   VRadMatEl(:,:,:) = 0.0_dp
   print*, "GET RADIAL MATRIX ELEMENTS"
   call getRadMatEl(basis, VPot, grid, VRadMatEl, indata)

   !If real spherical harmonics are used, V matrix elements will have complex part zero.
   call getVMatEl(sturm_ind_list, V, VRadMatEl, num_func, angular, indata)
   print*, "V Matrix Elements Computed"
   !V(:,:) = (-1.0_dp)*V(:,:)

   deallocate(angular)
   deallocate(VRadMatEl)
   H(:,:) = KMat(:,:) + V(:,:)

   !Call the rsg subroutine to solve the eigenvalue problem for the energies
   allocate(w(num_func),z(num_func,num_func))

   allocate(realH(num_func,num_func), realB(num_func,num_func))
   realH(:,:) = real(H(:,:))
   realB(:,:) = real(B(:,:))

   !MODIFY: USE LAPACK ZHEGVX FUNCTION INSTEAD FOR HERMITIAN MATRICES RATHER THAN REAL SYMMETRIC
   !ALTERNATIVELY, pass only real part to rsg
   print*, "CALL RSG"
   call rsg(num_func,num_func,realH,realB,w,1,z,ier)
   deallocate(realH,realB)

   !Include effect of nucleus-nucleus interaction term by simply adding zi*zj/R_ij to energies at the end.
   !Nuclear-nuclear matrix elements, being multiples of the overlap matrix, do not affect electronic wave functions.
   print*, "ADD NUCLEAR INTERACTION ENERGY"
   do ii=1, 3
      do jj = ii+1, 3
         !Use law of cosines to compute distance between nuclei
	 cosij = cos(theta(ii))*cos(theta(jj)) + sin(theta(ii))*sin(theta(jj))*cos(phi(ii)-phi(jj))
         Rij = sqrt(R(ii)**2 + R(jj)**2 - 2*R(ii)*R(jj)*cosij)
	 !Account for degenerate case where nuclei coincide
         if (Rij .gt. 0.0_dp) then
            w(:) = w(:) + dble(charge(ii)*charge(jj))/Rij
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
	   !if (indata%harmop .eq. 0) then
	   !else if (indata%harmop .eq. 1) then
	   !end if
	   !print*, ii, jj, kk, l_list(jj), k_list(jj), m_list(jj)
           !if (ii .eq. 430) then
           !   print*, "BASIS: ", radbasis(:,kk)
           !   print*, z(jj,ii)
           !   print*, tiny(w(1))
           !   print*, "TEST"
           !   print*, z(jj,ii)*radbasis(:,kk)
           !   print*, "AFTER TEST"
           !end if

	   !Basis function set to zero outside of range of r values indexed by i1, i2
	   i1 = basis%b(kk)%minf
	   i2 = basis%b(kk)%maxf
           
    	   wf(i1:i2,ii) = wf(i1:i2,ii) +  z(jj,ii)*basis%b(kk)%f(i1:i2)!*YLM
    	end do
   end do
   print*, "WAVE FUNCTIONS CALCULATED"

   nstates=num_func
   energies(:) = w(:)

   deallocate(k_list, l_list, m_list)
   deallocate(H,B,V,KMat)
   deallocate(w,z,wf)

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
   write(80,*) "R1=", R(1), " R2=", R(2), "R3=", R(3)
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
