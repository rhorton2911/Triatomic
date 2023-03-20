!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Module: symmetry.f90
!Purpose: contains subroutines for performing the symmetry 
!         adaptation of a basis for a point group using the
!         eigenfunction method. 
!Note: makes use of the fact that irreps are in 1-1 correspondence with
!      eigenvalues of a complete set of commuting operators (CSCO) that
!      can be constructed from the group's character table.
!Ref: Lemus, R. Symmetry (2012). 4, 667-685; doi:10.3390/sym4040667
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module symmetry

			!Stores character tables and irreducible representation labels for a given point group
			type, public:: group
			   character(len=6):: groupname  !name of the group
				 integer:: numirreps
			   character(len=6), dimension(:), allocatable:: irrlabel !Labels for the irreps
				 integer, dimension(:,:), allocatable:: chartable  !Character table for the group
				 integer, dimension(:), allocatable:: dims !dimension of each representation
				 integer, dimension(:), allocatable:: conorder !order of each conjugacy class
				 integer, dimension(:,:), allocatable:: lamtable  !Eigenvalue table for the group
				 real*8, dimension(:,:), allocatable:: euler !euler angles of the elements of classes defining irreps (alpha, beta, gamma)
				 integer, dimension(:), allocatable:: parities !equals 1 or 0 depending on whether or not, elements of class contain inversion
				 integer:: numops !number of group operations stored in euler

				 !Group operator C
				 real*8, dimension(:,:,:), allocatable:: coeffs


			end type group

			contains


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !Subroutine: init_group
	  !Purpose: initialises a group type based on input name
	  !         Character tables taken from: Bishop, M. Group Theory and
		!         Chemistry, 1993 dover edition.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		subroutine init_group(self, groupname)
		   implicit none
		   type(group):: self
       character(len=3):: groupname
			 integer:: numirrs, jj
			 real*8:: pi=4.0d0*atan(1.0d0)
			 real*8, dimension(3):: temp 

			 temp(:) = 0.0d0
			 if (trim(groupname) .eq. "D3h") then
			    numirrs = 6
			    self%numirreps = numirrs
			    allocate(self%irrlabel(numirrs), self%chartable(numirrs,numirrs), self%dims(numirrs), self%lamtable(numirrs,numirrs))
					allocate(self%conorder(numirrs))
					self%irrlabel(1) = 'A1p'
					self%irrlabel(2) = 'A2p'
					self%irrlabel(3) = 'Ep' 
					self%irrlabel(4) = 'A1pp' 
					self%irrlabel(5) = 'A2pp' 
					self%irrlabel(6) = 'Epp' 

					self%chartable(1,1:6) = (/1,  1,  1,  1,  1,  1/)
					self%chartable(2,1:6) = (/1,  1, -1,  1,  1, -1/)
					self%chartable(3,1:6) = (/2, -1,  0,  2, -1,  0/)	
					self%chartable(4,1:6) = (/1,  1,  1, -1, -1, -1/)
					self%chartable(5,1:6) = (/1,  1, -1, -1, -1,  1/)
					self%chartable(6,1:6) = (/2, -1,  0, -2,  1,  0/)

					self%dims(1:6) = (/1, 1, 2, 1, 1, 2/)
					self%conorder(1:6) = (/1, 2, 3, 1, 2, 3/)

					do jj = 1, numirrs  !classes
					   self%lamtable(:,jj) = self%chartable(:,jj)*self%conorder(jj)/self%dims(:)
				  end do

					!As per Lemus' paper, the eigenvalues of classes sigma_h and 3sigma_v can be used 
					!to uniquely label representations. Group operator C_{D_3h} = K_3 + K_4 has 
					!unique eigenvalues for all irreps.
					!Euler angles of class 3C_2 and sigma_h
					allocate(self%euler(4,3))
					self%numops = 4
					self%euler(:,:) = 0.0d0

					print*, "CHECK getEulerAngles WORKS!!!"
					!call getEulerAngles(0.0d0,0.0d0,1.0d0,pi, temp)             !sigma_h = I * R
					call axistoeuler(0.0d0, 0.0d0, pi,self%euler(1,:))
					call axistoeuler(pi/2.0d0, pi, pi, self%euler(2,:))                  !C_2,1
					call axistoeuler(pi/2.0d0, 2.0d0*pi/3.0d0, pi, self%euler(3,:))   !C_2,2
					call axistoeuler(pi/2.0d0, -2.0d0*pi/3.0d0, pi, self%euler(4,:))  !C_2,3

					allocate(self%parities(self%numops))
					self%parities(1) = 1
					self%parities(2:4) = 0
			 else if (trim(groupname) .eq. "C3v") then
					print*, "C3v"
					print*, "C3v not yet working, symmetry.f90 line 96"
					error stop
					numirrs = 3
			    self%numirreps = numirrs
			    allocate(self%irrlabel(numirrs), self%chartable(numirrs,numirrs), self%dims(numirrs), self%lamtable(numirrs,numirrs))
					allocate(self%conorder(numirrs))
					self%irrlabel(1) = 'A1'
					self%irrlabel(2) = 'A2'
					self%irrlabel(3) = 'E'

					self%chartable(1,1:3) = (/1,  1,  1/)
					self%chartable(2,1:3) = (/1,  1, -1/)
					self%chartable(3,1:3) = (/2, -1,  0/)

					self%dims(1:3) = (/1, 1, 2/)
					self%conorder(1:3) = (/1, 2, 3/)

					do jj = 1, numirrs  !classes
					   self%lamtable(:,jj) = self%chartable(:,jj)*self%conorder(jj)/self%dims(:)
				  end do

					!As per Lemus' paper, the eigenvalues of classes sigma_h and 3sigma_v can be used 
					!to uniquely label representations. Group operator C_{D_3h} = K_3 + K_4 has 
					!unique eigenvalues for all irreps.
					!Euler angles of class 3C_2 and sigma_h
					allocate(self%euler(3,3))
					self%numops = 3
					self%euler(:,:) = 0.0d0
					call axistoeuler(pi/2.0d0, pi/2.0d0, pi, self%euler(1,:))         !sigma_v,1
					call axistoeuler(pi/2.0d0, 5.0d0*pi/6.0d0, pi, self%euler(2,:))   !sigma_v,2
					call axistoeuler(pi/2.0d0, -5.0d0*pi/6.0d0, pi, self%euler(3,:))  !sigma_v,3
					print*, self%euler(1,:)
					print*, self%euler(2,:)
					print*, self%euler(3,:)

					allocate(self%parities(self%numops))
					self%parities(:) = 1
			 else
					print*, "ERROR: groups other than D3h not yet implemented, stopping."
					stop
			 end if
   end subroutine init_group


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: construct_rep
	 !Purpose: calculates the matrix representation of the operator C_G using the
	 !         wigner-d matrices, then diagonalises it to get expansion coefficients 
	 !         for symmetry adapted functions.
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine construct_rep(self, l)
			 use input_data
			 use realharms
			 use acc_wig
			 implicit none
			 !Data for the group and subgroup used to distinguish components of degenerate irreps
			 type(group):: self  !, subgroup  
			 integer:: l, m1, m2, par
			 complex*8, dimension(:,:), allocatable:: C, V, opMat
			 real*8, dimension(:,:), allocatable:: CReal, coeffs
			 real*8, dimension(:), allocatable:: eigs
			 complex*16:: i = (0.0d0,1.0d0)
			 real*8:: littled, t
			 integer:: ii, m, mu
			 integer:: numfuncs
			 integer:: m1ind, m2ind, muind, mind

			 !Construct C for the given l, 2l+1 values of m
			 numfuncs = (2*l+1)
			 if (data_in%harmop .eq. 1) then
          allocate(C(numfuncs,numfuncs), CReal(numfuncs,numfuncs), eigs(numfuncs), coeffs(numfuncs,numfuncs))
					allocate(opMat(numfuncs,numfuncs))
					C(:,:) = 0.0d0
					opMat(:,:) = 0.0d0

					!<Y_l mf | C | Y_l mi>
					!m1ind = 0
					!do m1 = -l, l
					!	 m1ind = m1ind + 1
					!   m2ind = 0
					!	 do m2 = -l, l
					!			m2ind = m2ind + 1
					!			do ii = 1, self%numops 
					!				 par = 1
					!				 par = (-1)**(l*self%parities(ii))  !If ii contains an inversion, multiply by (-1)^l
					!				 t = self%euler(ii,2) 
					!				 littled = wigd(l ,m1, m2, t)
					!				 !print*, "Vals: : ", self%euler(ii,3), self%euler(ii,2), self%euler(ii,1), littled
					!				 !D matrix formula taken from A. R. Edmonds: Angular momentum in quantum mechanics
					!         C(m1ind,m2ind) = C(m1ind,m2ind) + &
					!						dble(par) * CDEXP(-i*self%euler(ii,3)*cmplx(m2)) * littled * CDEXP(-i*self%euler(ii,1)*cmplx(m1))
					!		  end do
					!	 end do
				  !end do

          m1ind = 0
					do ii = 1, self%numops
             par = 1
             par = (-1)**(l*self%parities(ii))  !If ii contains an inversion, multiply by (-1)^l
             t = self%euler(ii,2) 
						 opMat(:,:) = 0.0d0

						 m1ind = 0
             do m1 = -l, l
             	  m1ind = m1ind + 1
                m2ind = 0
             	  do m2 = -l, l
             			 m2ind = m2ind + 1
             			 littled = wigd(l ,m1, m2, t)
             			 !print*, "Vals: : ", self%euler(ii,3), self%euler(ii,2), self%euler(ii,1), littled

             			 !D matrix formula taken from Varsholovich: Quantum Theory of Angular Momentum
									 opMat(m1ind,m2ind) = dble(par) * CDEXP(-i*self%euler(ii,3)*cmplx(m2)) * littled * CDEXP(-i*self%euler(ii,1)*cmplx(m1))
                   !C(m1ind,m2ind) = C(m1ind,m2ind) + &
             	  end do
             end do

   					!!Transform into real spherical harmonic basis: C' = V_(l) C V_(l)^(transpose)
   					!allocate(V(numfuncs,numfuncs))
   					!muind = 0
   					!do mu = -l, l
   					!	 muind = muind + 1
   					!	 mind = 0
   					!	 do m = -l, l
   					!			mind = mind + 1
   					!			V(muind,mind) = VCoeff(dble(m),dble(mu))
   					!	 end do
   				  !end do
   					!!Get C in real harmonic representation by changing basis
   					!!C = matmul(CONJG(TRANSPOSE(V)), matmul(C, V))
   					!opMat = matmul(V,matmul(opMat, CONJG(TRANSPOSE(V))))
   					!!C = CONJG(TRANSPOSE(V))* C * V
   					!deallocate(V)

						!print*, "OP: "
						!do m1 = 1, numfuncs
						! 	print*, opMat(ii,:)
						!end do
						C(:,:) = C(:,:) + opMat(:,:)
				  end do
         
					!Transform into real spherical harmonic basis: C' = V_(l) C V_(l)^(transpose)
					allocate(V(numfuncs,numfuncs))
					muind = 0
					do mu = -l, l
						 muind = muind + 1
						 mind = 0
						 do m = -l, l
								mind = mind + 1
								V(muind,mind) = VCoeff(dble(m),dble(mu))
						 end do
				  end do
					!Get C in real harmonic representation by changing basis
					!C = matmul(CONJG(TRANSPOSE(V)), matmul(C, V))
					C = matmul(V,matmul(C, CONJG(TRANSPOSE(V))))
					!C = CONJG(TRANSPOSE(V))* C * V
					deallocate(V)

				  !If done correctly, matrix will be real valued in real harmonic representation
					if (any(abs(AIMAG(C)) .gt. 10e-5)) then
						 print*, "ERROR: real harmonic representation has non-zero imaginary part, STOPPING"
						 error stop
				  end if
					CReal(:,:) = REAL(C)

					!Save eigenvector coefficients and eigenvalues
					coeffs(:,:) = 0.0d0
					eigs(:) = 0.0d0
					call project(CReal, numfuncs, coeffs, eigs) 
				  do m = 1, numfuncs
						 print*, coeffs(m,:)
				  end do
				  print*, "L: ", l
				  print*, eigs
				  stop
				 
					print*, "FINISH SYMMETRY ADAPTATION CODE, line 193 symmetry.f90"
					stop

					deallocate(C, CReal)
			else
				 print*, "ERROR: symmetry projection not yet implemented for complex spherical harmonics, STOPPING"
				 error stop
			end if

   end subroutine construct_rep



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: project
	 !Purpose: computes diagonalises the matrix representation of the operator C_II, used to
	 !         project a space of functions onto group irreps.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine project(CMat, spacedim, coeffs, eigs)
	     implicit none
			 integer:: spacedim   !Size of the basis set to be projected
			 real*8, dimension(spacedim,spacedim):: CMat
			 real*8, dimension(spacedim,spacedim):: coeffs 
			 real*8, dimension(spacedim):: eigs
       !Required inputs for lapack routine dsygvx
       real*8:: abstol
       real*8, dimension(:), allocatable:: w, work
       integer, dimension(:), allocatable:: iwork, ifail
       integer:: info, lwork, nfound
       real*8, external :: DLAMCH
			 !real*8, dimension(:,:), allocatable:: id
			 integer:: ii

			 !allocate(id(spacedim, spacedim))
			 !id(:,:) = 0.0d0
			 !do ii = 1, spacedim
			 ! 	id(ii,ii) = 1.0d0
			 !end do
			 coeffs(:,:) = 0.0d0

			 !Diagonalise the representation matrix
       allocate(work(1))
       allocate(w(spacedim))    !Stores eigenvalues
       allocate(iwork(5*spacedim))
       allocate(ifail(spacedim))
			 !allocate(alphar(spacedim), alphai(spacedim))
       lwork = -1
       abstol = 2.0d0*DLAMCH('S') !twice underflow threshold for doubles
       info = 0
       nfound = 0
        		 
       !Workspace query with lwork=-1
       !call DSYGVX(1, 'V', 'I', 'U', spacedim, CMat, spacedim, id, spacedim,    &
       !0.0d0, 0.0d0, 1, spacedim, abstol, nfound, w, coeffs, spacedim, work, lwork, &
       !iwork, ifail, info)
       call DSYEV('V', 'U', spacedim, CMat, spacedim, w, work, lwork, info)

       lwork = int(work(1))
       deallocate(work)
       allocate(work(lwork))
       !
       !Perform diagonialisation
       !call DSYGVX(1, 'V', 'I', 'U', spacedim, CMat, spacedim, id, spacedim,    &
       !0.0d0, 0.0d0, 1, spacedim, abstol, nfound, w, coeffs, spacedim, work, lwork, &
       !iwork, ifail, info) 
       call DSYEV('V', 'U', spacedim, CMat, spacedim, w, work, lwork, info)

		   eigs(:) = w(:)
			 !do ii =1, spacedim
			 ! 	print*, CMat(ii,:)
			 ! 	print*, coeffs(ii,:)
			 !end do
			 coeffs(:,:) = CMat(:,:)

			 !deallocate(id)
       deallocate(work, iwork, ifail, w)

   end subroutine project


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: getEulerAngles
	 !Purpose: using the axis-angle representation of a rotation for a given choice of basis,
	 !         computes the matrix of the rotation, followed by the euler angles of the rotation
	 !         from the matrix elements.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !subroutine getEulerAngles(nx, ny, nz, angle,euler)
	 ! 	implicit none
	 ! 	real*8, dimension(3):: euler
	 ! 	real*8:: angle   !Rotation angle
	 ! 	real*8, dimension(3,3):: R
	 ! 	real*8:: c, s
	 ! 	real*8:: nx, ny, nz  !Components of axis unit vector
	 ! 	integer:: ii
	 ! 	
	 ! 	print*, "HERE"
	 ! 	c = cos(angle)
	 ! 	s = sin(angle)

	 ! 	!!Compute rotation matrix elements
	 ! 	!R(:,:) = 0.0d0
	 ! 	!R(1,1) = c + (nx**2)*(1.0d0-c)
	 ! 	!R(1,2) = nx*ny*(1.0d0-c) - nz*s
	 ! 	!R(1,3) = nx*nz*(1.0d0-c) + ny*s
	 ! 	!R(2,1) = ny*nx*(1.0d0-c) + nz*s
	 ! 	!R(2,2) = c + (ny**2)*(1.0d0-c)
	 ! 	!R(2,3) = ny*nz*(1.0d0-c) - nx*s
	 ! 	!R(3,1) = nz*nx*(1.0d0-c) - ny*s
	 ! 	!R(3,2) = nz*ny*(1.0d0-c) + nx*s
	 ! 	!R(3,3) = c + (nz**2)*(1.0d0-c)

	 ! 	!!Compute euler angles
	 ! 	!euler(1) = atan(-R(2,3)/R(3,3))    !alpha
	 ! 	!euler(2) = atan(R(1,3)/sqrt(R(2,3)**2 + R(3,3)**2)) !beta
	 ! 	!euler(3) = atan(-R(1,2)/R(1,1))   !gamma
	 ! 	print*, "IN: ", nx, ny, nz, angle
	 ! 	!plus = 2.0d0*atan(

	 ! 	!do ii = 1, 3
	 ! 	!	 print*, R(ii,:)
	 ! 	!end do
	 ! 	!print*, euler(1), euler(2), euler(3) 

	 !end subroutine getEulerAngles



	 subroutine axistoeuler(theta, phi, angle, euler)
			implicit none
			real*8, dimension(3):: euler
			real*8:: angle   !Rotation angle
			real*8:: theta, phi
			real*8:: plus, minus
			real*8:: pi = 4.0d0*atan(1.0d0)


			plus = atan(cos(theta)*tan(angle/2.0d0))
			minus = phi - pi/2.0d0

			euler(2) = 2.0d0*asin(sin(theta)*sin(angle/2.0d0))
			euler(1) = plus + minus
			euler(3) = plus - minus

   end subroutine axistoeuler



end module symmetry
