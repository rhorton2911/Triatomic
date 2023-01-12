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
      !integer:: N   !Basis size
	    !integer:: Nalpha !nunber of exponential falloff parameters in input file
	    !real(dpf), dimension(:), allocatable:: alpha    !Exponential falloff parameters
	    !integer:: l         !Angular momentum parameter, max l of the basis
	    !real(dpf):: dr       !radial grid spacing
	    !real(dpf):: rmax     !number of radial grid points
	    integer:: lambdamax  !Maximum lambda in expansion of potential
	    integer:: harmop     !Option to use real or complex spherical harmonics (1=real, 0=complex)
	    !Properties and coordinates of each nucleus
	    integer, dimension(3):: charge  
	    real(dpf), dimension(3):: R
	    real(dpf), dimension(3):: theta
	    real(dpf), dimension(3):: phi
			!Parameters for isosceles triangle testing case
			integer:: isoscop   !(1=use, 0=regular calculation)
			integer:: hernandeztest
			integer:: isobasisop
			integer:: conroybasis
			real(dpf):: R1
			real(dpf):: R2
			integer, dimension(20):: testk, testl, testm
			integer, dimension(7):: testlvec, kminvec, kmaxvec
			integer:: numtestfuncs
			real(dpf):: enparam

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
			real(dpf), dimension(:), allocatable:: temp
			integer:: counter
			real(dpf):: angletemp !Temporary storage for angle in degrees
	    !Used to set parameters in isosceles testing case
	    real(dpf):: R1, R2, L, RN2, thetaN2

			open(20,file="input")

			!read(20,*) indata%N
			!read(20,*) indata%Nalpha

			!allocate(indata%alpha(indata%Nalpha))
			!read(20,*) (indata%alpha(ii), ii=1, indata%Nalpha)

			!read(20,*) indata%l
			!read(20,*) indata%dr
			!read(20,*) indata%rmax
			read(20,*) indata%lambdamax
			read(20,*) indata%harmop
			read(20,*)
			read(20,*) indata%charge(1)
			read(20,*) indata%charge(2)
			read(20,*) indata%charge(3)
			read(20,*) indata%R(1)
			read(20,*) indata%R(2)
			read(20,*) indata%R(3)
			read(20,*) angletemp
			indata%theta(1) = pinum*(angletemp/180.0_dpf)
			read(20,*) angletemp
			indata%theta(2) = pinum*(angletemp/180.0_dpf)
			read(20,*) angletemp
			indata%theta(3) = pinum*(angletemp/180.0_dpf)
			read(20,*) angletemp
			indata%phi(1) = pinum*(angletemp/180.0_dpf)
			read(20,*) angletemp
			indata%phi(2) = pinum*(angletemp/180.0_dpf)
			read(20,*) angletemp
			indata%phi(3) = pinum*(angletemp/180.0_dpf)
			read(20,*) 
			read(20,*) indata%isoscop
	    read(20,*) indata%hernandeztest
	    read(20,*) indata%isobasisop
	    read(20,*) indata%conroybasis
			read(20,*) indata%R1
			read(20,*) indata%R2
	    
      if (indata%R1 .gt. 2.0_dpf*indata%R2) then
				 print*, "ERROR: isoscelese triangles with R1 > 2*R2 do &not exist, stopping"
				 stop
	    end if

!			allocate(temp(indata%Nalpha))
!			temp(:) = 0.0_dpf
!			temp(:) = indata%alpha(:)
!			deallocate(indata%alpha)
!			allocate(indata%alpha(indata%l+1))
!			
!			!Method chosen is to change fallof parameter
!			!every three values of l
!			counter = 1 
!			do ii = 1, indata%l+1
!			   if (mod(ii,3) .eq. 0) then
!			      counter = counter + 1
!			   end if
!			   
!			   if (counter .le. indata%Nalpha) then
!			      indata%alpha(ii) = temp(counter)
!			   else 
!			      !Use the first alpha as default
!			      indata%alpha(ii) = temp(1)
!			   end if
!			   end do
!			deallocate(temp)
!

			if (indata%isoscop .eq. 1) then
			   !Set nuclear coordinates in case where iscosceles triangle testing mode is chosen
			   R1 = indata%R1
				 R2 = indata%R2
				 if (indata%hernandeztest .eq. 1) then
				    R1 = 3.50_dpf
						R2 = 2*R1*sin(pinum/6.0_dpf)
				 end if

				 L = sqrt(R2**2 - (0.5_dpf*indata%R1)**2)

				 indata%R(1) = (2.0_dpf/3.0_dpf)*L
				 RN2 = sqrt((0.5_dpf*R1)**2 + (L/3.0_dpf)**2)
				 indata%R(2) = RN2
				 indata%R(3) = RN2    !Symmetric case

         thetaN2 = acos((RN2**2+(2.0_dpf*L/3.0_dpf)**2-R2**2)/(2.0_dpf*RN2*(2.0_dpf*L/3.0_dpf)))
				 indata%theta(1) = 0.0_dpf
				 indata%theta(2) = thetaN2  !acos gives angle in radians
				 indata%theta(3) = thetaN2
				 indata%phi(1) = 0.0_dpf
				 indata%phi(2) = pinum*90.0_dpf/180_dpf 
				 indata%phi(3) = -indata%phi(2)  !Symmetric, x axis orthogonal to plane of molecule

				 !Set up arrays storing indices of basis functions, hardcode
				 !values from Conroy's paper.
	       indata%testm = (/0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, -1, 2, 2, 2, -3, -3/)
	       indata%testl = (/0, 0, 0, 0, 0, 2, 2, 2, 4, 1, 1, 1, 1, 3, 3, 2, 2, 2, 3, 3/)
	       indata%testk = (/1, 2, 3, 4, 5, 3, 4, 5, 5, 2, 3, 4, 5, 4, 5, 3, 4, 5, 4, 5/)	
				 indata%testlvec = (/0, 2, 3, 1, 3, 2, 3/)
				 indata%kminvec = (/1, 3, 5, 2, 4, 3, 4/)
				 indata%kmaxvec = (/5, 5, 5, 5, 5, 5, 5/)
				 indata%enparam = -3.788_dpf  !energy (a.u)
				 indata%numtestfuncs = 20
			end if

			
			close(20)
	 end subroutine readInput


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: deconstructInput
	 !Purpose: clears all memory in the input data type
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !subroutine deconstructInput(indata)
	 ! 	implicit none
	 ! 	type(smallinput):: indata
	 !

	 ! 	deallocate(indata%alpha)
	 !end subroutine deconstructInput





   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: createBasis
	 !Purpose: uses the recursion relations for the basis
	 !         functions to construct the basis set
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine createBasis(basis, lin, alphain, N, ingrid)
			use numbers
			use grid_radial
			implicit none
                        real(dpf), dimension(:,:):: basis
			type(rgrid):: ingrid
			integer::l, lin
			integer:: ii, N
			real(dpf):: alpha, alphain
			real(dpf):: coeff
			integer(kind=largeInt):: denom

			l = lin
			alpha = alphain

			!First construct without normalisation constant
			basis(:,1) = ((2.0_dpf)*alpha*ingrid%gridr(:))**(l+1)*exp(-alpha*ingrid%gridr(:))
			if (N .gt. 1) then
			   basis(:,2) = (2.0_dpf)*(real(l)+1.0_dpf-alpha*ingrid%gridr(:))*basis(:,1)
			end if

			!Next, use the recursion relations to fill out the rest of the basis
			do ii = 3, N
				 basis(:,ii) = (2.0_dpf*(dble(ii-1+l)-alpha*ingrid%gridr(:))*basis(:,ii-1) - &
				 dble(ii+2*l-1)*basis(:,ii-2))/dble(ii-1)
			end do 

			!Finally, normalise the basis functions
			coeff = 0.0_dpf
			denom = 0.0_dpf
			do ii =1, N
				 denom = getDenom(ii,l)
         !coeff = sqrt((alpha*dble(fact(ii-1)))/((dble(ii+l)*dble(fact(ii+2*l)))))
				 coeff = sqrt(alpha/dble(denom*(ii+l)))
				 basis(:,ii) = coeff*basis(:,ii)
			end do

	 end subroutine createBasis






	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: construct_all_conroy_basis
	 !Purpose: constructs the basis used by Conroy in in a calculation
	 !         of H3++ structure, adapted to reflect certain properties
	 !         of the H3++ molecule. Calculates B and K matrix elements as 
	 !         part of this.
	 !         Conroy, H. 1969. J. Chem. Phys. 51, 3979: Molecular Schrodinger Equation X
	 !Note: basis indices used here and in conroy are: (k,l,m) <--> (i,j,k)
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine construct_all_conroy_basis(self, ingrid, indata, dataARG, KMat)
			use grid_radial
      use sturmian_class
			use input_data
			implicit none
			type(basis_sturmian_nr), intent(inout):: self
			type(smallinput), intent(in):: indata
			type(input), intent(in):: dataARG
			type(rgrid):: ingrid
			complex(dpf), dimension(:,:), allocatable:: KMat
			real(dpf), dimension(:,:), allocatable:: f8
			real(dpf), dimension(:), allocatable:: gridr, func, temp
      real(dpf), dimension(:), allocatable:: J, s1, s2 !Function q(r), sigma_1(r), sigma_2(r)
      real(dpf), dimension(:), allocatable:: q !Function q(r), argument of laguerre polynomial
      real(dpf), dimension(:), allocatable:: rho 
			real(dpf):: a, e, gam, r0 !Parameters alpha, epsilon, gamma and rho
			integer:: n   !Maximum order of the laguerre polynomials
			integer:: labot, latop, k, l, kmin, kmax
			real(dpf):: lagarg, lambda
			integer:: ii, jj, kn, i1, i2

      call new_basis(self,indata%numtestfuncs)
			allocate(rho(ingrid%nr), gridr(ingrid%nr))
			gridr(:) = ingrid%gridr(:)

			e = sqrt(-2*indata%enparam)
			gam = dble(sum(indata%charge))
			a = gam/e - 1.0_dpf
			rho = sum(dble(indata%charge(:))*indata%R(:))/gam
			r0 = sum(dble(indata%charge(:))*indata%R(:))/gam

			labot = dataARG%labot
			latop  = dataARG%latop
			print*, e, gam, a, rho, r0

			allocate(s1(ingrid%nr), s2(ingrid%nr), J(ingrid%nr))
			allocate(func(ingrid%nr), temp(ingrid%nr), q(ingrid%nr))
			s1(:) = 1.0_dpf/sqrt(gridr(:)**2 + r0**2)
			s2(:) = 1.0_dpf/sqrt(gridr(:)**2 + r0**2 + e**(-2))
			J(:) = (s2(:)**(-a)) * exp(-gam*rho + (gam-e)/s1(:))

			open(70,file='Jfunc.txt')
			do ii=1, ingrid%nr
				 write(70,*) gridr(ii), (s2(ii)**(-a))
			end do
			close(70)
			stop

			!Argument of laguerre polynomials
			q(:) = 2.0_dpf*e*(1.0_dpf/s1(:) - r0)

			kn = 0 !Superindex for functions in basis

			do ii = 1, size(indata%testlvec)
				 l = indata%testlvec(ii)
				 !Get min and max k for this l
				 kmin = indata%kminvec(ii)
				 kmax = indata%kmaxvec(ii) 

				 n = kmax+l+1

				 allocate(f8(grid%nr,n))
				 f8(:,:) = 0.0_dpf
				 !Calculate generalised laguerre polynomials using lagpol8 from sturmian
				 lagarg = dble(2*l+1)
				 lambda = 1.0_dpf

				 !Call with q(r) as argument as per formula
				 call lagpol8(lagarg, lambda, f8, ingrid%nr, n, q, ingrid%nr)

				 !Lagpol calculates all polynomials from 0 to n-1, choose subset used in paper
				 do k = kmin, kmax
				    !Laguerre polynomials are left unnormalised in Conroy's
						!work, using polynomials of order (k+l)
						temp(:) = f8(:,k+l+1) 

						kn = kn + 1

						i1 = 1
						i2 = grid%nr

						!Calculated radial part of basis function
						func(:) = J(:)*temp(:)*(gridr(:)**l)*(s2(:)**(k-1))
						call init_function(self%b(kn), l, k, i1, i2, func, ingrid%nr,1.0_dpf)
				 end do
				 deallocate(f8)
			end do

			!Calculate the overlap matrix for the basis via integration
			do ii = 1, indata%numtestfuncs
				 do jj = 1, indata%numtestfuncs
				    if ((indata%testl(ii) .eq. indata%testl(jj)) .and. &
						    (indata%testm(ii) .eq. indata%testm(jj))) then
				       self%ortint(ii,jj) = sum(self%b(ii)%f(:)*self%b(jj)%f(:)*ingrid%weight(:))
				    else
							 self%ortint(ii,jj) = 0.0_dpf
						end if
				 end do
			end do

	    !Calculate K matrix elements
			!call getConroyKMat(self, KMat, grid, indata, s1, s2, J, q, e, gam, a, rho, r0)

			deallocate(gridr, s1, s2, J, temp, func, q, rho)
	 end subroutine construct_all_conroy_basis


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: getConroyKMat
	 !Purpose: calculates the K matrix elements for the basis used in
	 !Conroy's paper.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine getConroyKMat(self, KMat, ingrid, indata, s1, s2, J, q, e, gam, a, rho, r0)
			use grid_radial
      use sturmian_class
			use input_data
	    implicit none
			type(basis_sturmian_nr), intent(in):: self
			complex(dpf), dimension(:,:), allocatable:: KMat
			type(smallinput):: indata
	    type(rgrid):: ingrid
	    real(dpf), dimension(:):: s1, s2, J, q
	    real(dpf), dimension(:):: rho
	    real(dpf):: e, gam, a, r0  !Parameters used in basis def
	    !Derivatives of functions s1, s2, q and J
	    real(dpf), dimension(:), allocatable:: s1p, s1pp, s2p, s2pp
	    real(dpf), dimension(:), allocatable:: qp, qpp, Jp, Jpp
	    !Derivatives of the basis functions
	    real(dpf), dimension(:), allocatable:: psi, psip, psipp
	    real(dpf), dimension(:), allocatable:: t1, t2, t3
	    real(dpf), dimension(:,:), allocatable:: f1, f2, f3
	    integer:: ii, jj, k, l, kmax, n
	    real(dpf):: lagarg, lambda
	    !K operator acting on basis function
	    real(dpf), dimension(:), allocatable:: Kb
	    real(dpf), dimension(:), allocatable:: gridr

      allocate(KMat(indata%numtestfuncs, indata%numtestfuncs))
			KMat(:,:) = 0.0_dpf

	    allocate(s1p(ingrid%nr), s1pp(ingrid%nr), s2p(ingrid%nr), s2pp(ingrid%nr))
	    allocate(qp(ingrid%nr), qpp(ingrid%nr), Jp(ingrid%nr), Jpp(ingrid%nr))
	    allocate(Kb(ingrid%nr), psi(ingrid%nr), psip(ingrid%nr), psipp(ingrid%nr))
			allocate(t1(ingrid%nr), t2(ingrid%nr), t3(ingrid%nr))
			allocate(gridr(ingrid%nr))

			gridr(:) = ingrid%gridr(:)

			do ii = 1, indata%numtestfuncs
				 do jj = 1, indata%numtestfuncs 
				    if ((indata%testl(ii) .eq. indata%testl(jj)) .and. &
						    (indata%testm(ii) .eq. indata%testm(jj))) then
							 k = indata%testk(jj)
							 l = indata%testl(jj)

							 !Calculate using formulas for K[basis func]
							 s1p(:) = -gridr(:)*(s1(:)**3) 
							 s1pp(:) = -(s1(:)**2)*(s1(:) + 3.0_dpf*gridr(:)*s1p(:))
							 s2p(:) = -gridr(:)*(s2(:)**3)
							 s2pp(:) = -(s2(:)**2)*(s2(:) + 3.0_dpf*gridr(:)*s2p(:)) 
							 qp(:) = -2.0_dpf*e*(s1(:)**(-2))*s1p(:)
							 qpp(:) = -2.0_dpf*e*(s1(:)**(-2))*(-2.0_dpf*(s1(:)**(-1))*(s1p(:)**2) &
							          + s1pp(:))

							 !Calculate laguerre polynomials used in formulas
				       !Call with q(r) as argument as per formula
							 kmax = indata%kmaxvec(indata%numtestfuncs)
							 n = kmax + l + 1
							 allocate(f1(ingrid%nr,n), f2(ingrid%nr,n), f3(ingrid%nr,n))
							 lagarg = 2*l+1
							 lambda = 1.0_dpf
				       call lagpol8(lagarg, lambda, f1, grid%nr, n, q, ingrid%nr)
							 lagarg = 2*l+2
				       call lagpol8(lagarg, lambda, f2, grid%nr, n, q, ingrid%nr)
							 lagarg = 2*l+3
				       call lagpol8(lagarg, lambda, f3, grid%nr, n, q, ingrid%nr)
							 t1(:) = f1(:,k+l+1)
							 t2(:) = f2(:,k+l-1+1)
							 if ((k .eq. 1) .and. (l .eq. 0)) then
							    t3(:) = 0.0_dpf
						   else
							    t3(:) = f3(:,k+l-2+1)
						   end if

							 psi(:) = self%b(jj)%f(:)/J(:)

							 psip(:) = -t2(:)*qp(:)*(gridr(:)**l)*(s2(:)**(k-1)) &
							           + t1(:)*(dble(l)*(gridr(:)**(l-1))*(s2(:)**(k-1)) + &
							                   ((gridr(:)**l)*dble((k-1))*(s2(:)**(k-2))*s2p(:)))

							 psipp(:) = t3(:)*(qp(:)**2)*(gridr(:)**l)*(s2(:)**(k-1)) &
							            -t2(:)*(qpp(:)*(gridr(:)**l)*(s2(:)**(k-1))  &
							            + 2.0_dpf*qp(:)*dble(l)*(gridr(:)**(l-1))*(s2(:)**(k-1)) &
							            + 2.0_dpf*qp(:)*(gridr(:)**l)*dble(k-1)*(s2(:)**(k-2))*s2p(:)) &
							            + t1(:)*(2.0_dpf*dble(l)*(gridr(:)**(l-1))*dble(k-1)*(s2(:)**(k-2))*s2p(:) &
							            + (gridr(:)**l)*dble((k-1)*(k-2))*(s2(:)**(k-3))*s2p(:) &
							            + (gridr(:)**l)*dble(k-1)*(s2(:)**(k-2))*s2pp(:))

							 Jp(:) = J(:)*((-a/s2(:))*s2p(:) - (gam-e)*(s1(:)**(-2))*s1p(:))
							 Jpp(:) = Jp(:)*(Jp(:)/J(:)) + J(:)*((a/(s2(:)**2))*(s2p(:)**2) &
							 - (a/s2(:))*s2pp(:) + 2.0_dpf*(gam-e)*(s1(:)**(-3))*(s1p(:)**2) &
							 - (s1(:)**(-2))*s1pp(:)*(gam-e)) 

							 Kb(:) = psi(:)*(Jpp(:) + (2.0_dpf/gridr(:))*Jp(:) - &
							 (dble(l*(l+1))/(gridr(:)**2))*J(:)) + &
							 psip(:)*2.0_dpf*(J(:)/gridr(:) + Jp(:)) + &
							 psipp(:)*J(:) 
							 Kb(:) = (-0.5_dpf)*Kb(:)

							 !Explicitly compute integral
							 KMat(ii,jj) = sum(self%b(ii)%f(:)*Kb(:)*grid%weight(:))
							 deallocate(f1, f2, f3)
				    else
							 KMat(ii,jj) = 0.0_dpf
						end if
				 end do
			end do

	    deallocate(s1p, s1pp, s2p, s2pp, qp, qpp, Jp, Jpp)
			deallocate(Kb, psip, psipp, t1, t2, t3, gridr, f1, f2, f3)
	 end subroutine getConroyKMat

	 



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
	 subroutine getVPotNuc(ingrid, VPotTemp, R, theta, phi, z, indata)
	    use grid_radial
	    implicit none
	    type(rgrid):: ingrid
	    !real(dpf), allocatable, dimension(:):: rgrid
	    integer:: z
	    complex(dpf), allocatable, dimension(:,:):: VPotTemp
	    real(dpf):: R, theta, phi
	    type(smallinput):: indata
	    integer:: lambda, lambdaind, q
	    real(dpf), allocatable, dimension(:):: f
	    complex(dpf):: Ylm


	    if (abs(z) .gt. 0.0_dpf) then
	       allocate(f(ingrid%nr))

	       lambdaind = 1
	       do lambda= 0, indata%lambdamax
	          f(:) = dble(z)*(min(ingrid%gridr(:),R)**lambda)/(max(ingrid%gridr(:),R)**(lambda+1))

	          do q = -lambda, lambda 
	             if (indata%harmop .eq. 1) then
	                VPotTemp(:,lambdaind) = (4.0_dpf*pinum/(2.0_dpf*dble(lambda) + 1.0_dpf))*f(:)*XLM(lambda,q,theta,phi)
	             else if (indata%harmop .eq. 0) then
	                VPotTemp(:,lambdaind) = (4.0_dpf*pinum/(2.0_dpf*dble(lambda) + 1.0_dpf))*f(:)*YLM(lambda,q,theta,phi)
	             end if
	             lambdaind = lambdaind + 1
	          end do
	       end do
	       deallocate(f)
	    else
	       VPotTemp(:,:) = 0.0_dpf
	    end if
	    !Account for difference in charge of electron and proton
	    VPotTemp(:,:) = (-1.0_dpf)*VPotTemp(:,:) 

	 end subroutine getVPotNuc


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: getAngular
	 !Purpose: computes overlap of angular part of basis functions with nuclear 
	 !         potential. Called to precalculate these matrix elements before their 
	 !         use in computing potential matrix elements
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine getAngular(num_func, angular, l_list, m_List, indata)
	    implicit none
	    type(smallinput):: indata
	    real(dpf), dimension(:,:,:):: angular
	    integer:: num_func
	    integer, dimension(:):: l_list, m_list
	    integer:: lambdaind, lambda, q
	    integer:: ii, jj
	    real(dpf):: Yint
	    integer:: li, mi, lj, mj

            !Follow same indexing scheme as basis, specify lambda, then v from -lambda to lambda
						!$OMP PARALLEL DO DEFAULT(none) PRIVATE(jj, li, mi, lj, mj, lambdaind, lambda, q) &
						!$OMP& SHARED(angular, num_func, indata, l_list, m_list,pinum)
            do ii = 1, num_func
	             !print*, ii
               do jj = 1, num_func
                  li = l_list(ii)
                  mi = m_list(ii)
                  lj = l_list(jj)
                  mj = m_list(jj)
                  lambdaind=1
                  do lambda = 0, indata%lambdamax
                     do q = -lambda, lambda
                        if (indata%harmop .eq. 0) then
                           angular(lambdaind,ii,jj) = sqrt(dble(2*lambda+1)/(4.0_dpf*pinum)) &
                 	         *Yint(dble(li),dble(mi),dble(lambda),dble(q),dble(lj),dble(mj))
                        else if (indata%harmop .eq. 1) then
                 	         !Xint as written calculates overlap without sqrt(2lambda+1) factor, unlike Yint
                           angular(lambdaind,ii,jj) = Xint(dble(li),dble(mi),dble(lambda),dble(q),dble(lj),dble(mj))
                        end if

                        !print*, li, mi, lambda, q, lj, mj
                        !print*, angular(lambdaind,ii,jj)
                        lambdaind = lambdaind + 1 
                     end do
                  end do
               end do
            end do
						!$OMP END PARALLEL DO 

	 end subroutine getAngular



	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: getRadMatEl
	 !Purpose: computes the radial matrix elements of the potential VPot, stores 
	 !         inside array VRadMatEl, assumed to be allocated outside of this function
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine getRadMatEl(basis, VPot,ingrid, VRadMatEl, indata)
	    use grid_radial
	    use sturmian_class
	    implicit none
	    type(rgrid)::ingrid
	    complex(dpf), dimension(:,:):: VPot
	    complex(dpf), dimension(:,:,:):: VRadMatEl !Radial matrix elements to compute
	    type(smallinput):: indata
	    integer:: n, m !Loop indices
	    integer:: lambda, q, lambdaind
	    integer:: numrfuncs
	    complex(dpf), dimension(:), allocatable:: f
	    type(basis_sturmian_nr):: basis
	    integer:: nr
	    integer:: i1, i2

	    numrfuncs = basis%n
	    nr = ingrid%nr


	    !$OMP PARALLEL DO DEFAULT(none) PRIVATE(m, lambda, q, lambdaind, f, i1, i2) &
	    !$OMP& SHARED(VPot, ingrid, basis, indata, nr, numrfuncs, VRadMatEl)
	    do n = 1, numrfuncs
	       !print*, n
	       allocate(f(nr))
	       f(:) = 0.0_dpf
	       do m = 1, numrfuncs
		        i1 = max(basis%b(n)%minf, basis%b(m)%minf)
		        i2 = min(basis%b(n)%maxf, basis%b(m)%maxf)
		        lambdaind = 1
		        do lambda = 0, indata%lambdamax
		           do q = -lambda, lambda
	                f(i1:i2) = basis%b(n)%f(i1:i2) * VPot(i1:i2,lambdaind) * basis%b(m)%f(i1:i2)
			            VRadMatEl(lambdaind,n,m) = sum(f(i1:i2)*ingrid%weight(i1:i2))
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
	 subroutine getVMatEl(sturm_ind_list, V, VRadMatEl, num_func, angular, indata, use_list)
	    implicit none
	    integer, dimension(:):: sturm_ind_list
	    integer:: num_func
	    complex(dpf), dimension(:,:,:):: VRadMatEl
	    complex(dpf), dimension(:,:):: V
	    type(smallinput):: indata
	    integer:: lambda, lambdaind, q
	    integer:: ii,jj, n, m
	    real(dpf), dimension(:,:,:):: angular
			logical:: use1, use2
			logical, dimension(:):: use_list


	    !!!$OMP PARALLEL DO DEFAULT(none) PRIVATE(jj, n, m, lambdaind, lambda, q, integral) &
	    !!!$OMP& SHARED(V, indata, num_func, rad_ind_list, nr, angular, pi, use_list)
	    do ii=1, num_func
	       !print*, ii
	       do jj = 1, num_func
		        !n = sturm_ind_list(ii)
		        !m = sturm_ind_list(jj)
						n = ii
						m = jj

						use1 = use_list(ii)
						use2 = use_list(jj)

						if (use1 .and. use2) then
	             !Calculate V-matrix element for each lambda in the expansion
		           lambdaind = 1
		           do lambda = 0, indata%lambdamax
		              do q = -lambda, lambda
			               V(ii,jj) = V(ii,jj) + VRadMatEl(lambdaind,n,m)*angular(lambdaind,ii,jj)
			               lambdaind = lambdaind + 1
		              end do
		           end do
						end if

		        if (.not. (int(real(V(ii,jj))) .eq. int(real(V(ii,jj))))) then
		           print*, "INVALID MATRIX ELEMENT (ii,jj): ", ii, jj, V(ii,jj)
		           stop
		        end if
	       end do
	    end do
	    !!!$OMP END PARALLEL DO

	 end subroutine getVMatEl










	 !Deprecated subroutine
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: getVMat
	 !Purpose: computes VMatrix elements for a given nucleus using the given basis.
	 !         Arrays basis and V assumed to be allocated before call
	 !         to this subroutine with dimensions V(num_func,num_func)
	 !         and basis(nr,rad_func).
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !subroutine getVMat(basis, rgrid,  nr, rad_ind_list, V, num_func, R, z, theta, phi,  angular, weights, indata)
   !          implicit none
   !          type(smallinput):: indata
   !          integer:: nr, num_func
	 !    integer, dimension(:):: rad_ind_list
   !          real(dpf), dimension(:):: rgrid, weights
   !          real(dpf):: R  !Nuclear distance from origin
	 !    real(dpf):: theta, phi  !Nuclear coordinates
   !          real(dpf), dimension(:,:):: basis  !Radial part of the basis
   !          real(dpf), dimension(:), allocatable:: f, ratio
   !          complex(dpf), dimension(:,:):: V
	 !    real(dpf), dimension(:,:,:):: angular  !Precomputed angular integrals
   !          integer:: ii, jj, lambda, q, lambdaind
	 !    integer:: rii, rjj
   !          integer:: z !Charge of the nucleus
   !          real(dpf):: integral
	 !    complex(dpf):: ang_part  !Stores result of angular integration in matrix elements
	 !    complex(dpf):: YLM  !Define output of fortran 77 function YLM
   !          
   !          V(:,:) = 0.0_dpf

	 !    !V-matrix elements: calculate using numerical integration
	 !    !!!!$OMP PARALLEL DO DEFAULT(none) PRIVATE(jj, rii, rjj, lambdaind, lambda, q, ratio, f, integral, ang_part) &
	 !    !!!!$OMP& SHARED(V, rgrid, weights, basis, z, R, theta, phi, indata, num_func, rad_ind_list, nr, angular, pi)
	 !    do ii = 1, num_func
	 ! !Allocate arrays within loop to ensure thread safety
   !             allocate(f(nr))
   !             allocate(ratio(nr))
   !             f(:) = 0.0_dpf
   !             ratio(:) = 0.0_dpf
	 ! 
	 !     	do jj = 1, num_func 
	 !    !Get the indices of the radial basis functions
	 !    rii = rad_ind_list(ii)
	 !    rjj = rad_ind_list(jj)

	 !          !Calculate V-matrix element for each lambda in the expansion
	 !    lambdaind = 1
	 !          do lambda=0, indata%lambdamax
	 !             ratio(:) = (min(rgrid(:),R)**lambda)/(max(rgrid(:),R)**(lambda+1))
	 !     	      f(:) = basis(:,rjj) * dble(z) * ratio(:) * basis(:,rii)
	 !     	      integral = sum( f(:)*weights(:) )

	 !       !For non-diatomic molecules, the angular part includes weighting by 
	 !       !spherical harmonics of the nuclear coordinates.
	 !       ang_part = 0.0_dpf
	 !       do q = -lambda, lambda
	 ! 	 !Evaluate the angular integral, including spherical harmonics of the nuclear positions
	 ! 	 ang_part = ang_part + angular(lambdaind,ii,jj)*YLM(lambda,q,theta,phi)

	 ! 	 !thread = OMP_GET_THREAD_NUM()
	 ! 	 !print*, lambdaind, thread
	 ! 	 lambdaind = lambdaind + 1
	 !       end do
	 !     	      V(ii,jj) = V(ii,jj) + (4.0_dpf*pinum/dble(2*lambda+1))*integral*ang_part
	 !     	   end do
	 !     	end do

	 !       deallocate(ratio)
   !             deallocate(f)
	 !    end do
	 !    !!!!$OMP END PARALLEL DO

	 !end subroutine getVMat



	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: XLM
	 !Purpose: evaluates the real spherical harmonics for theta and phi in radians
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(dpf) function XLM(l, m, theta, phi) result(res)
	    use numbers
	    implicit none
	    integer:: l, m
	    real(dpf):: theta, phi
	    real(dpf):: RYLM

	    res = 0.0_dpf
	    if (m .gt. 0) then
	       res = RYLM(l,m,theta)*sqrt(2.0_dpf)*cos(dble(m)*phi)
	    else if (m .lt. 0) then
	       res = RYLM(l,abs(m),theta)*sqrt(2.0_dpf)*sin(dble(abs(m))*phi)
	    else if (m .eq. 0) then
	       !Special case, m=0 equivalent to theta=0, Ylm and Xlm coincide in this case.
	       res = RYLM(l,m,theta)
	    end if
	 end function XLM



	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: Xint
	 !Purpose: calculates an integral of the product of three real spherical harmonics
	 !         Xlm(theta,phi). Also known as real gaunt coefficients.
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 real(dpf) function Xint(l1, mu1, l2, mu2, l3, mu3) result(res)
	    use numbers
	    implicit none
	    real(dpf):: l1, mu1, l2, mu2, l3, mu3 
	    real(dpf)::Yint
	    real(dpf):: val
	    real(dpf):: t1, t2
	    real(dpf), dimension(3):: muvals, lVals
	    integer:: ii, numZero
	    integer, dimension(3):: isZero
	    real(dpf):: templ, tempmu
	    logical:: found
	    !complex(dpf):: v1, v2

			!val = sqrt(3.0_dpf/(4.0_dpf*pi))*Yint(1.0_dpf,1.0_dpf,1.0_dpf,1.0_dpf,0.0_dpf,0.0_dpf)
			!print*, val
	    !stop

	    !Use complete permutation symmetry of real gaunt coefficients to
			!reduce calculation to one of several cases. 
			!Case 1: mu1, mu2, mu3 all != 0
	    !Case 2: mu3 = 0, mu2, mu1 != 0
	    !Case 3: mu2 = mu3 = 0, mu1 != 0
	    !Simpler formulas apply in each of these cases
	    if (mod(int(l1+l2+l3),2) .ne. 0) then
	       res = 0.0_dpf
	    else
				 muVals(1) = mu1
				 muVals(2) = mu2
				 muVals(3) = mu3
				 lVals(1) = l1
				 lVals(2) = l2
				 lVals(3) = l3

				 numZero = 0
				 isZero(:) = 0
				 do ii = 1, 3
				    if (int(abs(muVals(ii))) .eq. 0) then
				       isZero(ii) = 1
						end if
				 end do
			   numZero = sum(isZero)
 
	       if (numZero .eq. 1) then
						found = .false.
						ii = 1 
						do while (isZero(ii) .ne. 1)
						   ii = ii + 1
				    end do
						if (ii .ne. 3) then
						   templ = l3
						   tempmu = mu3
				       l3 = lVals(ii)
						   mu3 = muVals(ii)
							 lVals(ii) = templ
							 muVals(ii) = tempmu
							 l1 = lVals(1)
							 mu1 = muVals(1)
							 l2 = lVals(2)
							 mu2 = muVals(2)
						end if

	          res = 2.0_dpf*sqrt((2.0_dpf*l2+1.0_dpf)/(4.0_dpf*pinum))*Yint(l1,mu2,l2,mu2,l3,0.0_dpf) &
	                *REAL(conjg(VCoeff(mu2,mu1))*VCoeff(mu2,mu2))
	       else if (numZero .eq. 2) then
						res = 0.0_dpf
				 else if (numZero .eq. 3) then
            !found = .false.
            !ii = 0
            !do while (.not. found)
            !   ii = ii + 1
            !   if (int(abs(muVals(ii))) .ne. 0) then
            !      found = .true.
            !   end if
            !end do
            !if (ii .ne. 1) then
            !   templ = lVals(ii)
            !   tempmu = muVals(ii)
            !   lVals(ii) = l1
            !   muVals(ii) = mu1
            !   l1 = templ
            !   mu1 = tempmu
            !   l2 = lVals(2)
            !   mu2 = muVals(2)
            !   l3 = lVals(3)
            !   mu3 = muVals(3)
            !end if
            
	          !Xl0 and Yl0 coincide, special case. Allows use of existing function for Ylm.
	          res =  sqrt((2.0_dpf*l2+1.0_dpf)/(4.0_dpf*pinum))*Yint(l1,0.0_dpf,l2,0.0_dpf,l3,0.0_dpf)
				 else if (numZero .eq. 0) then
		        !Formula for overlap of Xlm in terms of overlap of Ylm
						t1 = 2.0_dpf*sqrt((2.0_dpf*l2+1.0_dpf)/(4.0_dpf*pinum))*Yint(l1,mu2+mu3,l2,mu2,l3,mu3) &
		              *REAL(conjg(VCoeff(mu2+mu3,mu1))*VCoeff(mu2,mu2)*VCoeff(mu3,mu3))

						t2 =2.0_dpf*sqrt((2.0_dpf*l2+1.0_dpf)/(4.0_dpf*pinum))*Yint(l1,mu2-mu3,l2,mu2,l3,-mu3) &
							 *REAL(conjg(VCoeff(mu2-mu3,mu1))*VCoeff(mu2,mu2)*VCoeff(-mu3,mu3))
				    !if ((abs(t1) .gt. 0.0_dpf) .and. (abs(t2) .gt. 0.0_dpf)) then
						!	 print*, "ERROR", t1, t2
						!	 v1 = VCoeff(mu2+mu3,mu1)
						!	 v2 = VCoeff(mu2-mu3,mu1)
            !  print*, v1, abs(mu1), abs(mu2+mu3)
						!	 print*, v2, abs(mu1), abs(mu2-mu3)
						!	 print*, "DONE"


						!	 !print*, "ERROR: ", l1, mu1, l2, mu2, l3, mu3
						!	 !print*, "Vals: ", mu1,  mu2+mu3, -mu2-mu3, mu2-mu3,mu3-mu2
						!end if

						res = t1 + t2
		        !res = 2.0_dpf*sqrt((2.0_dpf*l2+1.0_dpf)/(4.0_dpf*pi))*Yint(l1,mu2+mu3,l2,mu2,l3,mu3) &
		        !      *REAL(conjg(VCoeff(mu2+mu3,mu1))*VCoeff(mu2,mu2) &
		        !      *VCoeff(mu3,mu3)) + 2.0_dpf*sqrt((2.0_dpf*l2+1.0_dpf)/(4.0_dpf*pi)) &
			      !*Yint(l1,mu2-mu3,l2,mu2,l3,-mu3)*REAL(conjg(VCoeff(mu2-mu3,mu1)) &
			      !*VCoeff(mu2,mu2)*VCoeff(-mu3,mu3))
				 else
				    print*, "TREAT ALL CASES!!!"
						stop
	       end if
	    end if

	 end function Xint


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: VCoeff
	 !Purpose: evaluates the coefficients in the unitary transformation from 
	 !         complex to real spherical harmonics.
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 complex(dpf) function VCoeff(m,mu) result(res)
	    implicit none
	    real(dpf):: m,mu
	    complex(dpf):: i

	    !if (abs(m) .gt. l) then
	    !   print*, l, m
	    !   print*, "ERROR: angular momentum projection greater than magnitude, stopping."
	    !   stop
	    !end if

	    i = (0.0_dpf,1.0_dpf)

	    res = kronecker(int(m),0)*kronecker(int(mu),0) &
			      + (1.0_dpf/sqrt(2.0_dpf))*(heavyside(mu)*kronecker(int(m),int(mu)) &
	          + heavyside(-mu)*(i)*((-1.0_dpf)**m)*kronecker(int(m),int(mu)) &
				    + heavyside(-mu)*(-i)*kronecker(int(m),int(-mu)) &
		        + heavyside(mu)*((-1.0_dpf)**m)*kronecker(int(m),int(-mu)))
	 end function VCoeff


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: heavyside
	 !Purpose: evaluates the heavyside function of the argument
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 real(dpf) function heavyside(m) result(res)
	    implicit none
	    real(dpf):: m

	    if (m .gt. 0.0_dpf) then
	       res = 1.0_dpf
	    else 
	       res = 0.0_dpf
	    end if
	 end function heavyside


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: kronecker
	 !Purpose: evaluates the kronecker delta as a real number
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 real(dpf) function kronecker(n,m) result(res)
	    implicit none
	    integer:: n,m

	    res = 0.0_dpf
	    if (n .eq. m) then
	       res = 1.0_dpf
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
	 real(dpf) function getPLMNorm(l,m,u) result(res)
	    implicit none
	    real(dpf):: u   !Argument of the polynomial
	    integer:: l, m
	    real(dpf):: PLM  !Define data type for output of PLM f77 function

	    res = getCoeff(l,abs(m))*PLM(l,abs(m),u)
	    if (m .lt. 0) then
	       !Use formula for product of c_lm and P_lm with sign of m reversed
	       res = ((-1.0_dpf)**m)*res
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
      real(dpf) function getCoeff(l,m) result(coeff)
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

         coeff = ((-1.0_dpf)**m)*sqrt(dble(2*l+1)/4*pinum)*sqrt(1.0_dpf/dble(denom))

      end function getCoeff



!###############################################################################
!###############################################################################

end module basismodule

