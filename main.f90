!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program: H3Plus
!Purpose: performs a structure calculation for the 
!         H3++ molecule. First step towards H3+ structure.
!Author: Reese Horton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program H3Plus
   use basismodule
   use numbers
   implicit none
   
   type(input):: indata
   integer:: nr  !number of radial grid points
   integer:: ii, jj, kk
   real(dp), dimension(:), allocatable:: rgrid, w
   real(dp), dimension(:,:), allocatable:: z, wf, basisatl !, B, H, KMat, 
   complex(dp), dimension(:,:), allocatable:: V, VTemp, H, KMat, B
   real(dp), dimension(:,:), allocatable:: realH, realB
   real(dp), dimension(:,:), allocatable:: radbasis
   real(dp), dimension(:,:,:), allocatable:: angular
   integer, dimension(:), allocatable:: k_list, l_list, m_list, rad_ind_list
   real(dp), dimension(:), allocatable:: weights, energies
   integer:: num_func, l, m, k, rad_func
   integer:: li, mi, lambda, q, lj, mj, num_lambda, lambdaind
   integer:: lblocksize
   integer:: lmax, N, m_max
   integer:: nstates, ier
   real(dp):: alphal, mu
   real(dp), dimension(3):: R, theta, phi !Coordinates and charge of each nucleus
   integer, dimension(3):: charge
   real(dp):: Rij, cosij
   logical:: sorted
   real(dp):: E1, E2, temp
   real(dp):: Yint !declare a type for Yint function output

   !Legacy .f77 subroutine. Sets up common (global :( ) variables used by function YLM in plql.f
   call FAKRED
    
   call readInput(indata)
   N = indata%N
   m_max = indata%mmax
   lmax = indata%l
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
  
   !Set up the radial grid over which to define the basis functions
   nr = int(indata%rmax/indata%dr)
   if (mod(nr,2) .eq. 0) then
      nr = nr + 1
   end if
   
   allocate(rgrid(nr))
   rgrid(1) = indata%dr
   do ii = 2, nr
      rgrid(ii) = rgrid(ii-1) + indata%dr
   end do
   
   !Set up integration weights for computing V-matrix elements
   allocate(weights(nr))
   weights(1) = 1.0_dp
   do ii = 2, nr-1
      weights(ii) = 2.0_dp + 2.0*mod(ii+1,2)
   end do
   weights(nr) = 1.0_dp
   weights(:) = weights(:)*indata%dr / 3.0_dp 
   
   !Triatomic molecules do not have good quantum numbers such as angular momentum and parity (like H2+).
   !Instead, the clamped nuclei electronic V-matrix is diagonal in different irreducible representations of a 
   !given point group. Since we are not currently using a symmetry adapted basis, we can only loop over l
   num_func=0
   do l = 0, lmax
      do m = -l, l
	 !Now need to loop over l and m, m is no longer conserved.
         num_func = num_func + indata%N
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
   allocate(k_list(num_func), l_list(num_func), m_list(num_func), rad_ind_list(num_func))
   ii = 0
   rad_func = 0
   kk = 0
   do l = 0, lmax
      do m = -l, l
	 jj  = kk 
         do k = 1, N
            ii = ii + 1
            k_list(ii) = k
            l_list(ii) = l
	    m_list(ii) = m
	    !print*, "L,M,k,ii: ", l, m, k, ii

	    !Indexing changed for non-diatomics, indices for basis now include l, k AND m. For a given l, the same radial basis
	    !functions phi_{kl} are used for each value of m from -l to l. Introduce an array to track which radial basis 
	    !function belongs to which index.
	    jj = jj + 1
	    rad_ind_list(ii) = jj
	    !print*, jj
         end do
      end do
      !The radial part of the basis is identical for any given (l,k) pair, regardless of m
      do k = 1, N
         rad_func = rad_func + 1
      end do
      kk = kk + N
   end do
   
   !Construct a basis with the given symmetry properties using the l_list and 
   !k_list variables
   allocate(radbasis(nr,rad_func))
   l=-1
   kk=0
   jj=0
   do ii=1,  num_func
      if (l .eq. l_list(ii)) then
      	 cycle
      end if
      l = l_list(ii)
      !In principle allows to specify a different alpha for each l
      alphal = indata%alpha(l+1)
      jj=kk+1
      kk = kk+N
   
      allocate(basisatl(nr, indata%N))
      call createBasis(basisatl, l, alphal, indata%N, rgrid)
      radbasis(:,jj:kk) = basisatl(:,1:N)
      deallocate(basisatl)
   end do

   !Write some of the basis functions to file for plotting
   open(70,file="basis.txt")
   write(70,*) "l, lmax,N: ", 0, lmax, N
   do jj=1, nr
      write(70,*) rgrid(jj), (radbasis(jj,ii), ii=1, indata%N)
   end do
   close(70)
   
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
      alphal=indata%alpha(l+1)
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
      alphal = indata%alpha(l+1)
      !Assuming we choose the same range of k for each l
      lblocksize = N*(2*l+1)
      jj=kk+1
      kk = kk + lblocksize
      KMat(jj:kk,jj:kk) = KMat(jj:kk,jj:kk) - 0.5_dp*(alphal**2)*B(jj:kk,jj:kk)
   end do
   mu = 1.0_dp
   KMat(:,:) = KMat(:,:)/mu

   !Precalculate angular integrals appearing in V-matrix elements
   num_lambda = 0
   do lambda = 0, indata%lambdamax
      num_lambda = num_lambda + 2*lambda + 1
   end do
   allocate(angular(num_lambda,num_func,num_func))
   !Follow same indexing scheme as basis, specify lambda, then v from -lambda to lambda
   do ii = 1, num_func
      do jj = 1, num_func
	 li = l_list(ii)
	 mi = m_list(ii)
	 lj = l_list(jj)
	 mj = m_list(jj)
	 lambdaind=1
         do lambda = 0, indata%lambdamax
	    do q = -lambda, lambda
	       angular(lambdaind,ii,jj) = sqrt(dble(2*lambda+1)/(4.0_dp*pi))*Yint(dble(li),dble(mi),dble(lambda),dble(q),dble(lj),dble(mj))
	       !print*, li, mi, lambda, q, lj, mj
	       !print*, angular(lambdaind,ii,jj)
	       lambdaind = lambdaind + 1 
	    end do
	 end do
      end do
   end do
	
   !Calculate V matrix elements, loop over nuclei
   allocate(VTemp(num_func,num_func))
   do ii = 1, 3
      VTemp(:,:) = 0.0_dp
      call getVMat(radbasis, rgrid, nr, rad_ind_list, VTemp, num_func, &
	           R(ii), charge(ii), theta(ii), phi(ii), angular, weights, indata)
      V(:,:) = V(:,:) + VTemp(:,:)  
   end do
   deallocate(VTemp)
   V(:,:) = (-1.0_dp)*V(:,:)

   deallocate(angular)
   H(:,:) = KMat(:,:) + V(:,:)

   !Call the rsg subroutine to solve the eigenvalue problem for the energies
   allocate(w(num_func),z(num_func,num_func))

   allocate(realH(num_func,num_func), realB(num_func,num_func))
   realH(:,:) = real(H(:,:))
   realB(:,:) = real(B(:,:))

   !MODIFY: USE LAPACK ZHEGVX FUNCTION INSTEAD FOR HERMITIAN MATRICES RATHER THAN REAL SYMMETRIC
   !ALTERNATIVELY, pass only real part to rsg
   call rsg(num_func,num_func,realH,realB,w,1,z,ier)
   deallocate(realH,realB)

   !Include effect of nucleus-nucleus interaction term by simply adding zi*zj/R_ij to energies at the end.
   !Nuclear-nuclear matrix elements, being multiples of the overlap matrix, do not affect electronic wave functions.
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

   open(80,file="energies.txt") 
   do ii = 1, num_func
      write(80,*) ii, w(ii)
   end do
   close(80)

   !Need to divide by r and multiply by Ylm to get full 3D wave function
   allocate(wf(nr,num_func))
   wf(:,:) = 0.0_dp
   do ii = 1, num_func
    	do jj=1, num_func
           kk = rad_ind_list(jj)
           !MODIFY TO BE A FUNCTION OF (r,theta,phi)
    	   wf(:,ii) = wf(:,ii) +  z(jj,ii)*radbasis(:,kk)!*YLM
    	end do
   end do

   nstates=num_func
   energies(:) = w(:)

   deallocate(k_list, l_list, m_list, rad_ind_list)
   deallocate(H,B,V,KMat)
   deallocate(w,z,wf)
   deallocate(radbasis)

   !Print out the first energy in each symmetry
   !ii= 1
   !do m=0, m_max
   !   do par=-1,1,2 
   !   	 print*, "M= ", m, "par= ", par, ": ", energies(ii), " Ha"
   !   	 ii = ii + 1 
   !   end do
   !end do

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
   deallocate(rgrid)
   deallocate(weights)

end program H3Plus
