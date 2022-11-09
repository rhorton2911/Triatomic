!Submodule: basismodule
!Purpose: contains the subroutines to construct a basis set for the structure calculation 
!         Written a part of a review of the CCC method as applied to single electron 
!         molecules. Required as a prerequisite to using configuration interaction codes.
!Author: Reese Horton


module basismodule
	 use numbers

   !Type: input
	 !Purpose: strores program input parameters read in from input file
	 type smallinput
            integer:: N   !Basis size
	    integer:: Nalpha !nunber of exponential falloff parameters in input file
	    real(dp), dimension(:), allocatable:: alpha    !Exponential falloff parameters
	    integer:: l         !Angular momentum parameter, max l of the basis
	    real(dp):: dr       !radial grid spacing
	    real(dp):: rmax     !number of radial grid points
	    integer:: lambdamax  !Maximum lambda in expansion of potential
	    integer:: mmax       !Maximum angular momentum projection onto Rhat
	    integer:: harmop     !Option to use real or complex spherical harmonics (1=real, 0=complex)
	    !Properties and coordinates of each nucleus
	    integer, dimension(3):: charge  
	    real(dp), dimension(3):: R
	    real(dp), dimension(3):: theta
	    real(dp), dimension(3):: phi

	    !integer:: z1         !Charge of nucleus 1
	    !integer:: z2         !Charge of nucleus 2
	    !integer:: z3         !Charge of nucleus 3
	 end type smallinput


	 contains



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: readInput
	 !Purpose: reads input parameters from input file
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine readInput(indata)
			implicit none
			type(smallinput)::indata
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
			read(20,*) indata%harmop
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
			type(smallinput):: indata
	 

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



	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: getVPotNuc
	 !Purpose: evaluates the components in the expansion of the nuclear potential in 
	 !         real or complex spherical harmonics for use in the structure calculation.
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine getVPotNuc(rgrid, nr, VPotTemp, R, theta, phi, z, indata)
	    !use grid_radial
	    implicit none
	    real(dp), allocatable, dimension(:):: rgrid
	    integer:: nr, z
	    complex(dp), allocatable, dimension(:,:):: VPotTemp
	    real(dp):: R, theta, phi
	    type(smallinput):: indata
	    integer:: lambda, lambdaind, q
	    real(dp), allocatable, dimension(:):: f
	    complex(dp):: Ylm

	    allocate(f(nr))

	    lambdaind = 1
	    do lambda= 0, indata%lambdamax
	       f(:) = dble(z)*(min(rgrid(:),R)**lambda)/(max(rgrid(:),R)**(lambda+1))

	       do q = -lambda, lambda 
		  if (indata%harmop .eq. 1) then
		     VPotTemp(:,lambdaind) = (4.0_dp*pi/(2.0_dp*dble(lambda) + 1.0_dp))*f(:)*Xlm(lambda,q,theta,phi)
		  else if (indata%harmop .eq. 0) then
		     VPotTemp(:,lambdaind) = (4.0_dp*pi/(2.0_dp*dble(lambda) + 1.0_dp))*f(:)*YLM(lambda,q,theta,phi)
	          end if
		  lambdaind = lambdaind + 1
	       end do
	    end do
	    deallocate(f)

	 end subroutine getVPotNuc



	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: getRadMatEl
	 !Purpose: computes the radial matrix elements of the potential VPot, stores 
	 !         inside array VRadMatEl, assumed to be allocated outside of this function
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine getRadMatEl(radbasis,VPot, nr, weights, VRadMatEl, indata)
	    implicit none
	    real(dp), dimension(:,:)::radbasis
	    complex(dp), dimension(:,:):: VPot
	    integer:: nr
	    real(dp), dimension(:):: weights
	    complex(dp), dimension(:,:,:):: VRadMatEl !Radial matrix elements to compute
	    type(smallinput):: indata
	    integer:: n, m !Loop indices
	    integer:: lambda, q, lambdaind
	    integer:: numrfuncs
	    complex(dp), dimension(:), allocatable:: f

	    numrfuncs = size(radbasis(1,:))


	    print*, "COMPUTE RADIAL MATRIX ELEMENTS"
	    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(m, lambda, q, lambdaind, f) &
	    !$OMP& SHARED(VPot, weights, radbasis, indata, nr, numrfuncs, VRadMatEl)
	    do n = 1, numrfuncs
	       !print*, n
	       allocate(f(nr))
	       do m = 1, numrfuncs
		  lambdaind = 1
		  do lambda = 0, indata%lambdamax
		     do q = -lambda, lambda
	      	        f(:) = radbasis(:,n) * VPot(:,lambdaind) * radbasis(:,m)
			VRadMatEl(lambdaind,n,m) = sum(f(:)*weights(:))
			lambdaind = lambdaind + 1
		     end do
		  end do
	       end do
	       deallocate(f)
	    end do
	    !$OMP END PARALLEL DO

	 end subroutine getRadMatEl



	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: getVMatEl
	 !Purpose: calculates the potential matrix elements using the precalculated 
	 !         expansion of the potential in spherical harmonics.
	 !Note: handles both real and complex spherical harmonics. Formulas are the 
	 !      same for appropriately calculated VPot and angular integrals
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine getVMatEl(rad_ind_list, V, VRadMatEl, num_func, angular, indata)
	    implicit none
	    integer, dimension(:):: rad_ind_list
	    integer:: num_func
	    complex(dp), dimension(:,:,:):: VRadMatEl
	    complex(dp), dimension(:,:):: V
	    type(smallinput):: indata
	    integer:: lambda, lambdaind, q
	    integer:: ii,jj, n, m
	    real(dp), dimension(:,:,:):: angular



	    !!!$OMP PARALLEL DO DEFAULT(none) PRIVATE(jj, n, m, lambdaind, lambda, q, integral) &
	    !!!$OMP& SHARED(V, indata, num_func, rad_ind_list, nr, angular, pi)
	    do ii=1, num_func
	       !print*, ii
	       do jj = 1, num_func
		  n = rad_ind_list(ii)
		  m = rad_ind_list(jj)

	          !Calculate V-matrix element for each lambda in the expansion
		  lambdaind = 1
		  do lambda = 0, indata%lambdamax
		     do q = -lambda, lambda
			V(ii,jj) = V(ii,jj) + VRadMatEl(lambdaind,n,m)*angular(lambdaind,ii,jj)
			lambdaind = lambdaind + 1
		     end do
		  end do
		  if (.not. (int(real(V(ii,jj))) .eq. int(real(V(ii,jj))))) then
		     print*, "INVALID MATRIX ELEMENT (ii,jj): ", ii, jj, V(ii,jj)
		     stop
		  end if
	       end do
	    end do
	    !!!$OMP END PARALLEL DO

	 end subroutine getVMatEl


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: getVMat
	 !Purpose: computes VMatrix elements for a given nucleus using the given basis.
	 !         Arrays basis and V assumed to be allocated before call
	 !         to this subroutine with dimensions V(num_func,num_func)
	 !         and basis(nr,rad_func).
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine getVMat(basis, rgrid,  nr, rad_ind_list, V, num_func, R, z, theta, phi,  angular, weights, indata)
             implicit none
             type(smallinput):: indata
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

	     !V-matrix elements: calculate using numerical integration
	     !!!!$OMP PARALLEL DO DEFAULT(none) PRIVATE(jj, rii, rjj, lambdaind, lambda, q, ratio, f, integral, ang_part) &
	     !!!!$OMP& SHARED(V, rgrid, weights, basis, z, R, theta, phi, indata, num_func, rad_ind_list, nr, angular, pi)
	     do ii = 1, num_func
		!Allocate arrays within loop to ensure thread safety
                allocate(f(nr))
                allocate(ratio(nr))
                f(:) = 0.0_dp
                ratio(:) = 0.0_dp
		
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
			 ang_part = ang_part + angular(lambdaind,ii,jj)*YLM(lambda,q,theta,phi)

			 !thread = OMP_GET_THREAD_NUM()
			 !print*, lambdaind, thread
			 lambdaind = lambdaind + 1
		      end do
	      	      V(ii,jj) = V(ii,jj) + (4.0_dp*pi/dble(2*lambda+1))*integral*ang_part
	      	   end do
	      	end do

	        deallocate(ratio)
                deallocate(f)
	     end do
	     !!!!$OMP END PARALLEL DO


	 end subroutine getVMat



	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: XLM
	 !Purpose: evaluates the real spherical harmonics for theta and phi in radians
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         real(dp) function XLM(l, m, theta, phi) result(res)
	    use numbers
	    implicit none
	    integer:: l, m
	    real(dp):: theta, phi
	    real(dp):: RYLM

	    res = 0.0_dp
	    if (m .gt. 0) then
	       res = RYLM(l,m,theta)*sqrt(2.0_dp)*cos(dble(m)*phi)
	    else if (m .lt. 0) then
	       res = RYLM(l,abs(m),theta)*sqrt(2.0_dp)*sin(dble(abs(m))*phi)
	    else if (m .eq. 0) then
	       !Special case, m=0 equivalent to theta=0, Ylm and Xlm coincide in this case.
	       res = RYLM(l,m,theta)
	    end if
	 end function



	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: Xint
	 !Purpose: calculates an integral of the product of three real spherical harmonics
	 !         Xlm(theta,phi)
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 real(dp) function Xint(l1, mu1, l2, mu2, l3, mu3) result(res)
	    use numbers
	    implicit none
	    real(dp):: l1, mu1, l2, mu2, l3, mu3 
	    real(dp)::Yint

	    !if (mod(int(l1+l2+l3),2) .ne. 0) then
	    !   res = 0.0_dp
	    !else
	    !   if ((abs(mu2) .gt. 0.0_dp) .and. (int(mu3) .eq. 0)) then
	    !      res = 2.0_dp*sqrt((2.0_dp*l2+1.0_dp)/(4.0_dp*pi))*Yint(l1,mu2,l2,mu2,l3,0.0_dp) &
	    !            *REAL(conjg(VCoeff(mu2,mu1))*VCoeff(mu2,mu2))
	    !   else if ((int(mu2) .eq. 0) .and. (int(mu3) .eq. 0)) then
	    !      if (int(mu1) .eq. 0) then
	    !         !Xl0 and Yl0 coincide, special case. Allows use of existing function for Ylm.
	    !         res = sqrt((2.0_dp*l2+1.0_dp)/(4.0_dp*pi))*Yint(l1,0.0_dp,l2,0.0_dp,l3,0.0_dp)
	    !      else 
	    !         res = 0.0_dp
	    !      end if
	    !   else 
		  !Formula for overlap of Xlm in terms of overlap of Ylm
		  res = 2.0_dp*sqrt((2.0_dp*l2+1.0_dp)/(4.0_dp*pi))*Yint(l1,mu2+mu3,l2,mu2,l3,mu3) &
		        *REAL(conjg(VCoeff(mu2+mu3,mu1))*VCoeff(mu2,mu2) &
		        *VCoeff(mu3,mu3)) + 2.0_dp*sqrt((2.0_dp*l2+1.0_dp)/(4.0_dp*pi)) &
			*Yint(l1,mu2-mu3,l2,mu2,l3,-mu3)*REAL(conjg(VCoeff(mu2-mu3,mu1)) &
			*VCoeff(mu2,mu2)*VCoeff(-mu3,mu3))
	    !   end if
	    !end if

	 end function Xint



	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: VCoeff
	 !Purpose: evaluates the coefficients in the unitary transformation from 
	 !         complex to real spherical harmonics.
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 complex(dp) function VCoeff(m,mu) result(res)
	    implicit none
	    real(dp):: m,mu
	    complex(dp):: i

	    !if (abs(m) .gt. l) then
	    !   print*, l, m
	    !   print*, "ERROR: angular momentum projection greater than magnitude, stopping."
	    !   stop
	    !end if

	    i = (0.0_dp,1.0_dp)

	    res = kronecker(int(m),0)*kronecker(int(mu),0) + (1.0_dp/sqrt(2.0_dp))*(heavyside(mu)*kronecker(int(m),int(mu)) &
	          + heavyside(-mu)*(i)*((-1.0_dp)**m)*kronecker(int(m),int(mu)) + heavyside(-mu)*(-i)*kronecker(int(m),int(-mu)) &
		  + heavyside(mu)*((-1.0_dp)**m)*kronecker(int(m),int(-mu)))
	 end function VCoeff


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: heavyside
	 !Purpose: evaluates the heavyside function of the argument
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 real(dp) function heavyside(m) result(res)
	    implicit none
	    real(dp):: m

	    if (m .gt. 0.0_dp) then
	       res = 1.0_dp
	    else 
	       res = 0.0_dp
	    end if
	 end function heavyside


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: kronecker
	 !Purpose: evaluates the kronecker delta as a real number
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 real(dp) function kronecker(n,m) result(res)
	    implicit none
	    integer:: n,m

	    res = 0.0_dp
	    if (n .eq. m) then
	       res = 1.0_dp
	    end if
	 end function kronecker








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

