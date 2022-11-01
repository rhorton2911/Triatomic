!Submodule: basismodule
!Purpose: contains the subroutines to construct a basis set for the structure calculation 
!         Written a part of a review of the CCC method as applied to single electron 
!         molecules. Required as a prerequisite to using configuration interaction codes.
!Author: Reese Horton


module basismodule
	 use numbers

   !Type: input
	 !Purpose: strores program input parameters read in from input file
	 type input
            integer:: N   !Basis size
	    integer:: Nalpha !nunber of exponential falloff parameters in input file
	    real(dp), dimension(:), allocatable:: alpha    !Exponential falloff parameters
	    integer:: l         !Angular momentum parameter, max l of the basis
	    real(dp):: dr       !radial grid spacing
	    real(dp):: rmax     !number of radial grid points
	    integer:: lambdamax  !Maximum lambda in expansion of potential
	    integer:: mmax       !Maximum angular momentum projection onto Rhat
	    !Properties and coordinates of each nucleus
	    integer, dimension(3):: charge  
	    real(dp), dimension(3):: R
	    real(dp), dimension(3):: theta
	    real(dp), dimension(3):: phi

	    !integer:: z1         !Charge of nucleus 1
	    !integer:: z2         !Charge of nucleus 2
	    !integer:: z3         !Charge of nucleus 3
	 end type input


	 contains



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: readInput
	 !Purpose: reads input parameters from input file
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine readInput(indata)
			implicit none
			type(input)::indata
			integer:: ii
			real(dp), dimension(:), allocatable:: temp
			integer:: counter
			real(dp):: angletemp !Temporary storage for angle in degrees

			open(20,file="input")

			read(20,*) indata%N
			read(20,*) indata%Nalpha

			allocate(indata%alpha(indata%Nalpha))
			read(20,*) (indata%alpha(ii), ii=1, indata%Nalpha)

			read(20,*) indata%l
			read(20,*) indata%dr
			read(20,*) indata%rmax
			read(20,*) indata%lambdamax
			read(20,*) indata%mmax
			read(20,*)
			read(20,*) indata%charge(1)
			read(20,*) indata%charge(2)
			read(20,*) indata%charge(3)
			read(20,*) indata%R(1)
			read(20,*) indata%R(2)
			read(20,*) indata%R(3)
			read(20,*) angletemp
			indata%theta(1) = pi*(angletemp/180.0_dp)
			read(20,*) angletemp
			indata%theta(2) = pi*(angletemp/180.0_dp)
			read(20,*) angletemp
			indata%theta(3) = pi*(angletemp/180.0_dp)
			read(20,*) angletemp
			indata%phi(1) = pi*(angletemp/180.0_dp)
			read(20,*) angletemp
			indata%phi(2) = pi*(angletemp/180.0_dp)
			read(20,*) angletemp
			indata%phi(3) = pi*(angletemp/180.0_dp)

			allocate(temp(indata%Nalpha))
			temp(:) = 0.0_dp
			temp(:) = indata%alpha(:)
			deallocate(indata%alpha)
			allocate(indata%alpha(indata%l+1))
			
			!Method chosen is to change fallof parameter
			!every three values of l
			counter = 1 
			do ii = 1, indata%l+1
			   if (mod(ii,3) .eq. 0) then
			      counter = counter + 1
			   end if
			   
			   if (counter .le. indata%Nalpha) then
			      indata%alpha(ii) = temp(counter)
			   else 
			      !Use the first alpha as default
			      indata%alpha(ii) = temp(1)
			   end if
			   end do
			deallocate(temp)
			
			close(20)
	 end subroutine readInput


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: deconstructInput
	 !Purpose: clears all memory in the input data type
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine deconstructInput(indata)
			implicit none
			type(input):: indata
	 

			deallocate(indata%alpha)
	 end subroutine deconstructInput





   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: createBasis
	 !Purpose: uses the recursion relations for the basis
	 !         functions to construct the basis set
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine createBasis(basis, lin, alphain, N, rgrid)
			use numbers
			implicit none
      real(dp), dimension(:,:):: basis
			real(dp), dimension(:):: rgrid
			integer::l, lin
			integer:: ii, N
			real(dp):: alpha, alphain
			real(dp):: coeff
			integer(kind=largeInt):: denom

			l = lin
			alpha = alphain

			!First construct without normalisation constant
			basis(:,1) = ((2.0_dp)*alpha*rgrid(:))**(l+1)*exp(-alpha*rgrid(:))
			if (N .gt. 1) then
			   basis(:,2) = (2.0_dp)*(real(l)+1.0_dp-alpha*rgrid(:))*basis(:,1)
			end if

			!Next, use the recursion relations to fill out the rest of the basis
			do ii = 3, N
				 basis(:,ii) = (2.0_dp*(dble(ii-1+l)-alpha*rgrid(:))*basis(:,ii-1) - &
				 dble(ii+2*l-1)*basis(:,ii-2))/dble(ii-1)
			end do 

			!Finally, normalise the basis functions
			coeff = 0.0_dp
			denom = 0.0_dp
			do ii =1, N
				 denom = getDenom(ii,l)
         !coeff = sqrt((alpha*dble(fact(ii-1)))/((dble(ii+l)*dble(fact(ii+2*l)))))
				 coeff = sqrt(alpha/dble(denom*(ii+l)))
				 basis(:,ii) = coeff*basis(:,ii)
			end do

	 end subroutine createBasis





	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: getDenom
	 !Purpose: calculates the denominator in the normalisation
	 !         coefficients for the basis functions to avoid
	 !         taking the ratio of two factorials.
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 integer(kind=largeInt) function getDenom(ii,l) result(res)
			use numbers
			implicit none
			integer:: ii, l
			integer:: kk

			res = 1
			do kk = ii, ii+2*l
				 res = res*kk
			end do
	 end function






	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !function: fact
	 !Purpose: calculates the factorial of the argument,
	 !         uses integers with 18 places due to large size of 
	 !         factorials
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 integer(kind=largeInt) recursive function fact(n) result(res)
			use numbers
			implicit none
			integer:: n

			if (n .lt. 0) then
				 print*, "ERROR: argument to factorial function cannot be a negative integer, stopping"
				 stop
			end if

			if (n .eq. 0) then
				 res = 1
	    end if
	 
			if (n .gt. 0) then
	       res = n*fact(n-1) 
			end if
	 end function fact




	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: getVMat
	 !Purpose: computes the VMatrix elements using the given basis.
	 !         Arrays basis and V assumed to be allocated before call
	 !         to this subroutine with dimensions V(num_func, num_func)
	 !         and basis(nr, num_func)
	 !Date last modified: 
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !subroutine getVMat(basis, rgrid,  nr, num_func, V, l_list, R1, R2, m, weights, indata)
         !    implicit none
         !    type(input):: indata
         !    integer:: nr, num_func, m
         !    real(dp), dimension(:):: rgrid, weights
         !    real(dp):: R1, R2  !Nuclear distances from origin
         !    real(dp), dimension(:,:):: basis
         !    real(dp), dimension(:), allocatable:: f, ratio1, ratio2
         !    real(dp), dimension(:,:):: V
         !    integer, dimension(:):: l_list
         !    integer:: ii, jj, li, lj, lambda
         !    integer:: z1, z2
         !    real(dp):: integral
         !    real(dp):: Yint !declare a type for Yint function output
         !    
         !    V(:,:) = 0.0_dp
         !    
         !    z1 = indata%z1
         !    z2 = indata%z2
         !    
         !    allocate(f(nr))
         !    allocate(ratio1(nr), ratio2(nr))	
         !    f(:) = 0.0_dp
         !    ratio1(:) = 0.0_dp
         !    ratio2(:) = 0.0_dp
         !    
	 !    !V-matrix elements: calculate using numerical integration
	 !    !$OMP PARALLEL DO
	 !    do ii = 1, num_func
	 !     	do jj = 1, num_func 
	 !          li = l_list(ii)
	 !     	   lj = l_list(jj)

	 !          !Calculate V-matrix element for each lambda in the expansion
	 !          lambda = 0
	 !          do while (lambda .le. indata%lambdamax)
	 !             ratio1(:) = (min(rgrid(:),R1)**lambda)/(max(rgrid(:),R1)**(lambda+1))
	 !             ratio2(:) = (min(rgrid(:),R2)**lambda)/(max(rgrid(:),R2)**(lambda+1))
	 !     	      f(:) = basis(:,jj) * (z1*ratio1(:) + z2*((-1.0_dp)**lambda)*ratio2(:)) * basis(:,ii)
	 !     	      integral = sum( f(:)*weights(:) )
	 !     	      V(ii,jj) = V(ii,jj) + integral*Yint(dble(li),dble(m),dble(lambda),0.0_dp,dble(lj),dble(m))

	 !     	      lambda = lambda + 1
	 !     	   end do
	 !     	end do
	 !    end do
	 !    !$OMP END PARALLEL DO
	 !    V(:,:) = (-1.0_dp)*V(:,:)

	 !    deallocate(ratio1, ratio2)
         !    deallocate(f)

	 !end subroutine getVMat


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: getVMat
	 !Purpose: computes VMatrix elements for a given nucleus using the given basis.
	 !         Arrays basis and V assumed to be allocated before call
	 !         to this subroutine with dimensions V(num_func,num_func)
	 !         and basis(nr,rad_func).
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine getVMat(basis, rgrid,  nr, rad_ind_list, V, num_func, R, z, theta, phi,  angular, weights, indata)
             implicit none
             type(input):: indata
             integer:: nr, num_func
	     integer, dimension(:):: rad_ind_list
             real(dp), dimension(:):: rgrid, weights
             real(dp):: R  !Nuclear distance from origin
	     real(dp):: theta, phi  !Nuclear coordinates
             real(dp), dimension(:,:):: basis  !Radial part of the basis
             real(dp), dimension(:), allocatable:: f, ratio
             complex(dp), dimension(:,:):: V
	     real(dp), dimension(:,:,:):: angular  !Precomputed angular integrals
             integer:: ii, jj, lambda, q, lambdaind
	     integer:: rii, rjj
             integer:: z !Charge of the nucleus
             real(dp):: integral
	     complex(dp):: ang_part  !Stores result of angular integration in matrix elements
	     complex(dp):: YLM  !Define output of fortran 77 function YLM
             
             V(:,:) = 0.0_dp
             
             allocate(f(nr))
             allocate(ratio(nr))
             f(:) = 0.0_dp
             ratio(:) = 0.0_dp

	     !V-matrix elements: calculate using numerical integration
	     !$OMP PARALLEL DO
	     do ii = 1, num_func
	      	do jj = 1, num_func 
		   !Get the indices of the radial basis functions
		   rii = rad_ind_list(ii)
		   rjj = rad_ind_list(jj)

	           !Calculate V-matrix element for each lambda in the expansion
		   lambdaind = 1
	           do lambda=0, indata%lambdamax
	              ratio(:) = (min(rgrid(:),R)**lambda)/(max(rgrid(:),R)**(lambda+1))
	      	      f(:) = basis(:,rjj) * dble(z) * ratio(:) * basis(:,rii)
	      	      integral = sum( f(:)*weights(:) )

		      !For non-diatomic molecules, the angular part includes weighting by 
		      !spherical harmonics of the nuclear coordinates.
		      ang_part = 0.0_dp
		      do q = -lambda, lambda
			 !Evaluate the angular integral, including spherical harmonics of the nuclear positions
			 !print*, YLM(lambda,q,theta,phi)
			 !print*, angular(lambdaind,ii,jj)
			 ang_part = ang_part + angular(lambdaind,ii,jj)*YLM(lambda,q,theta,phi)
			 lambdaind = lambdaind + 1
		      end do
	      	      V(ii,jj) = V(ii,jj) + integral*ang_part
	      	   end do
	      	end do
	     end do
	     !$OMP END PARALLEL DO

	     deallocate(ratio)
             deallocate(f)
	 end subroutine getVMat










!Deprecated functions, kept for possible future generalisation to use 
!real spherical harmonics.
!###############################################################################
!###############################################################################
	 

	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: getPLMNorm
	 !Purpose: evaluates the associated legendre polynomial multiplied by
	 !         normalisation coefficient for spherical harmonics. 
	 !         Accounts for requirement that m be positive enforced by
	 !         standard legendre polynomial subroutine.
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 real(dp) function getPLMNorm(l,m,u) result(res)
	    implicit none
	    real(dp):: u   !Argument of the polynomial
	    integer:: l, m
	    real(dp):: PLM  !Define data type for output of PLM f77 function

	    res = getCoeff(l,abs(m))*PLM(l,abs(m),u)
	    if (m .lt. 0) then
	       !Use formula for product of c_lm and P_lm with sign of m reversed
	       res = ((-1.0_dp)**m)*res
	    end if
	 end function getPLMNorm



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Function: getCoeff
      !Purpose: computes the normalisation coefficient of the spherical
      !         harmonics. The denominator in the normalisation coefficients
      !         of the spherical harmonics, needs to be stored as a
      !         very large integer due to the use of factorials in its
      !         definition.
      !Input: requires l .ge. 0 and m .ge. 0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(dp) function getCoeff(l,m) result(coeff)
         use numbers
         integer:: l,m
         integer:: ii
         integer(kind=largeInt):: denom

         if (m .lt. 0) then
            print*, "ERROR: normalisation coefficients for spherical harmonics defined only for positive m. Stopping"
            stop
         end if

         denom = 1
         if (m .gt. 0) then
            do ii = l-m+1, l+m
               denom = denom*ii
            end do
         end if

         coeff = ((-1.0_dp)**m)*sqrt(dble(2*l+1)/4*pi)*sqrt(1.0_dp/dble(denom))

      end function getCoeff



!###############################################################################
!###############################################################################

end module basismodule

