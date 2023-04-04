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
			   character(len=6), dimension(:), allocatable:: irrlabel !Irrep labels
				 integer, dimension(:,:), allocatable:: chartable       !Character table
				 integer, dimension(:), allocatable:: dims              !dimension of each representation
				 integer, dimension(:), allocatable:: conorder          !order of each conjugacy class
				 integer, dimension(:,:), allocatable:: lamtable        !Eigenvalue table 
				 real*8, dimension(:,:), allocatable:: euler            !euler angles of the elements of classes defining irreps (alpha, beta, gamma)
				 integer, dimension(:), allocatable:: parities          !equals 1 or 0 depending on whether or not, elements of class contain inversion
				 integer, dimension(:), allocatable:: primea            !prime number coefficients of class operator
				 integer:: numops !number of group operations stored in euler

				 !Expansion coefficients of symmetry adapted harmonics
				 real*8, dimension(:,:,:), allocatable:: coeffs


			end type group


			!Canonical chain of groups: G > H1 > H2 > ...
			type, public:: groupchain
				 integer:: numgroups
				 type(group), dimension(:), allocatable:: chain
         !Frequency of occurrence of subgroup irreps in group irreps (subduced representation)
				 real*8, dimension(:,:,:), allocatable:: amat   !index=(group, group irreps, subgroup irreps)
				 integer, dimension(:,:), allocatable:: subconjarray

				 !Expansion coefficients of symmetry adapted harmonics
				 real*8, dimension(:,:,:), allocatable:: coeffs


	    end type




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

			 print*, "Constructing Symmetry Group: ", groupname

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

					call axistoeuler(0.0d0, 0.0d0, pi,self%euler(1,:))
					call axistoeuler(pi/2.0d0, pi, pi, self%euler(2,:))                  !C_2,1
					call axistoeuler(pi/2.0d0, 2.0d0*pi/3.0d0, pi, self%euler(3,:))   !C_2,2
					call axistoeuler(pi/2.0d0, -2.0d0*pi/3.0d0, pi, self%euler(4,:))  !C_2,3

					allocate(self%parities(self%numops))
					self%parities(1) = 1
					self%parities(2:4) = 0

					allocate(self%primea(self%numops))
					self%primea(:) = 1 

					!allocate(self%subnames(1))
					!self%subnames(1) = "C2v"
					!
					!allocate(self%amat(1,self%numirreps, 4))

			 else if (trim(groupname) .eq. "C2v") then
					numirrs = 4
					self%numirreps = numirrs
			    allocate(self%irrlabel(numirrs), self%chartable(numirrs,numirrs), self%dims(numirrs), self%lamtable(numirrs,numirrs))
					allocate(self%conorder(numirrs))
					self%irrlabel(1) = 'A1'
					self%irrlabel(2) = 'A2'
					self%irrlabel(3) = 'B1' 
					self%irrlabel(4) = 'B2' 

					self%chartable(1,1:4) = (/1,  1,  1,  1/)
					self%chartable(2,1:4) = (/1,  1, -1, -1/)
					self%chartable(3,1:4) = (/1, -1,  1, -1/)
					self%chartable(4,1:4) = (/1, -1, -1,  1/)

					self%dims(1:4) = (/1, 1, 1, 1/)
					self%conorder(1:4) = (/1, 1, 1, 1/)

					do jj = 1, numirrs  !classes
					   self%lamtable(:,jj) = self%chartable(:,jj)*self%conorder(jj)/self%dims(:)
				  end do

					!Choose the classes K_3 and K_4, both plane reflections. In this case, choose C_{C_2v} = K_3 + 3 K_4
					allocate(self%euler(2,3))
					self%numops = 2
					self%euler(:,:) = 0.0d0

					call axistoeuler(pi/2.0d0, pi/2.0d0, pi, self%euler(1,:))
					call axistoeuler(pi/2.0d0, 0.0d0, pi, self%euler(2,:))

					allocate(self%parities(self%numops))
					self%parities(1) = 1
					self%parities(2) = 1

					allocate(self%primea(self%numops))
					self%primea(1) = 1 
					self%primea(2) = 3

			! else if (trim(groupname) .eq. "C3v") then
			!		print*, "C3v"
			!		print*, "C3v not yet working, symmetry.f90 line 96"
			!		error stop
			!		numirrs = 3
			!    self%numirreps = numirrs
			!    allocate(self%irrlabel(numirrs), self%chartable(numirrs,numirrs), self%dims(numirrs), self%lamtable(numirrs,numirrs))
			!		allocate(self%conorder(numirrs))
			!		self%irrlabel(1) = 'A1'
			!		self%irrlabel(2) = 'A2'
			!		self%irrlabel(3) = 'E'

			!		self%chartable(1,1:3) = (/1,  1,  1/)
			!		self%chartable(2,1:3) = (/1,  1, -1/)
			!		self%chartable(3,1:3) = (/2, -1,  0/)

			!		self%dims(1:3) = (/1, 1, 2/)
			!		self%conorder(1:3) = (/1, 2, 3/)

			!		do jj = 1, numirrs  !classes
			!		   self%lamtable(:,jj) = self%chartable(:,jj)*self%conorder(jj)/self%dims(:)
			!	  end do

			!		!As per Lemus' paper, the eigenvalues of classes sigma_h and 3sigma_v can be used 
			!		!to uniquely label representations. Group operator C_{D_3h} = K_3 + K_4 has 
			!		!unique eigenvalues for all irreps.
			!		!Euler angles of class 3C_2 and sigma_h
			!		allocate(self%euler(3,3))
			!		self%numops = 3
			!		self%euler(:,:) = 0.0d0
			!		call axistoeuler(pi/2.0d0, pi/2.0d0, pi, self%euler(1,:))         !sigma_v,1
			!		call axistoeuler(pi/2.0d0, 5.0d0*pi/6.0d0, pi, self%euler(2,:))   !sigma_v,2
			!		call axistoeuler(pi/2.0d0, -5.0d0*pi/6.0d0, pi, self%euler(3,:))  !sigma_v,3
			!		print*, self%euler(1,:)
			!		print*, self%euler(2,:)
			!		print*, self%euler(3,:)

			!		allocate(self%parities(self%numops))
			!		self%parities(:) = 1


			!		allocate(self%primea(self%numops))
			!		self%primea(:) = 1
			 else
					print*, "ERROR: groups other than D3h and C2v not yet implemented, stopping."
					stop
			 end if
   end subroutine init_group


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: destruct_group
	 !Purpose: deallocates all allocated memory used by the given group object
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine destruct_group(self)
			implicit none
			type(group):: self

			if (allocated(self%irrlabel)) then
				 deallocate(self%irrlabel)
			end if
			if (allocated(self%chartable)) then
				 deallocate(self%chartable)
			end if
			if (allocated(self%dims)) then
				 deallocate(self%dims)
			end if
			if (allocated(self%conorder)) then
				 deallocate(self%conorder)
			end if
			if (allocated(self%lamtable)) then
				 deallocate(self%lamtable)
			end if
			if (allocated(self%euler)) then
				 deallocate(self%euler)
			end if
			if (allocated(self%parities)) then
				 deallocate(self%parities)
			end if
			if (allocated(self%primea)) then
				 deallocate(self%primea)
			end if
			if (allocated(self%coeffs)) then
				 deallocate(self%coeffs)
			end if

   end subroutine destruct_group


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: construct_chain
	 !Purpose: constructs a canonical chain of groups given the input group name, which
	 !         stays at the top of the chain
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine construct_chain(self, groupname)
			implicit none
			character(len=3):: groupname
			type(groupchain):: self
			integer:: ii, numirrs1, numirrs2, numirrs
			integer, dimension(:,:), allocatable:: temp

			if (TRIM(groupname) .eq. "D3h") then
				 self%numgroups = 2
				 allocate(self%chain(2))
				 call init_group(self%chain(1), "D3h")
				 call init_group(self%chain(2), "C2v")
				 allocate(self%subconjarray(1,self%chain(2)%numirreps))
				 self%subconjarray(1,1) = 1   !E
				 self%subconjarray(1,2) = 3   !C_2
				 self%subconjarray(1,3) = 6   !sigma_v(xz)
				 self%subconjarray(1,4) = 4   !sigma_v(yz) = sigma_h 

			else
				 print*, "ERROR: canonical chain of groups implemented only for D3h, STOPPING"
				 error stop
			end if


			!Calculate clebsch-gordan decomposition of each group irrep in the chain
			numirrs = self%chain(1)%numirreps
			allocate(self%amat(self%numgroups-1, numirrs, numirrs))
			self%amat(:,:,:) = 0

			do ii = 1, self%numgroups-1
				 numirrs1=self%chain(ii)%numirreps
				 numirrs2=self%chain(ii+1)%numirreps
				 allocate(temp(numirrs1, numirrs2))

				 call subdecomp(self%chain(ii), self%chain(ii+1), self%subconjarray(ii,:), temp)
				 self%amat(ii,1:numirrs1,1:numirrs2) = temp(1:numirrs1,1:numirrs2)

				 deallocate(temp)
			end do



	 end subroutine construct_chain



	 !TODO: change this so that it does this for a canonical CHAIN of groups

	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: construct_harmonics
	 !Purpose: constructs the set of expansion coefficients a^{\mu}_{lm} for the
	 !         symmetry adapted spherical harmonics for a given point group G up to 
	 !         the specified lmax.
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine construct_harmonics(self)
			use input_data
			implicit none
			type(group):: self
			integer:: maxl, maxnum, l

		  maxl = max(data_in%latop, data_in%Lpmax, data_in%l12max)
			maxnum = (maxl+1)**2
			allocate(self%coeffs(1:maxnum,1:maxnum,0:maxl))
			self%coeffs(:,:,:) = 0.0d0

			do l = 0, maxl
				 call construct_rep_group(self,l)
			end do

	 end subroutine construct_harmonics




	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: construct_rep
	 !Purpose: constructs the representation matrix of the 'group' operator 
	 !         C = C_G + C_H + .., whose diagonalisation yields expansion coefficients
	 !         for symmetry adapted functions.
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine construct_rep(self,l,maxl)
			implicit none
			type(groupchain):: self
			integer:: ii
			integer:: maxl, maxnum, l, numfuncs
			integer:: m
			real*8, dimension(:,:), allocatable:: C, CTemp, coeffs
			real*8, dimension(:), allocatable:: eigs

			maxnum = (maxl+1)**2

			allocate(self%coeffs(1:maxnum,1:maxnum,0:maxl))
			self%coeffs(:,:,:) = 0.0d0

			do l = 0, maxl
				 numfuncs = (2*l+1)**2
				 allocate(C(numfuncs,numfuncs), CTemp(numfuncs,numfuncs))
				 allocate(eigs(numfuncs), coeffs(numfuncs,numfuncs))
			   C(:,:) = 0.0d0

			   do ii = 1, self%numgroups
						CTemp(:,:) = 0.0d0
						!call construct_rep_group(self%chain(ii), l, CTemp)
						C(:,:) = C(:,:) + CTemp(:,:) 
			   end do

         !Save eigenvector coefficients and eigenvalues
         coeffs(:,:) = 0.0d0
         eigs(:) = 0.0d0
         call project(C, numfuncs, coeffs, eigs) 
         
         !Correct sign problem 
         do m = 1, numfuncs
            if (sum(coeffs(:,m)) .lt. 0.0d0) then
               coeffs(:,m) = (-1.0d0)*coeffs(:,m)
            end if
         end do
				 self%coeffs(1:numfuncs,1:numfuncs,l) = coeffs(1:numfuncs,1:numfuncs)

				 deallocate(C, CTemp)
			end do
         
	 end subroutine construct_rep



	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: construct_rep_group
	 !Purpose: calculates the matrix representation of the operator C_G using the
	 !         wigner-d matrices, then diagonalises it to get expansion coefficients 
	 !         for symmetry adapted functions.
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine construct_rep_group(self, l)
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
       allocate(C(numfuncs,numfuncs), CReal(numfuncs,numfuncs), eigs(numfuncs), coeffs(numfuncs,numfuncs))
			 allocate(opMat(numfuncs,numfuncs))
			 C(:,:) = 0.0d0
			 opMat(:,:) = 0.0d0


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

				  C(:,:) = C(:,:) + dble(self%primea(ii)) * opMat(:,:)
			 end do
         
			 if (data_in%harmop .eq. 1) then
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
			    !Get C in real harmonic representation by changing basis: C' = VCV^T
			    C = matmul(V,matmul(C, CONJG(TRANSPOSE(V))))
			    deallocate(V)

			    !If done correctly, matrix will be real valued in real harmonic representation
			    if (any(abs(AIMAG(C)) .gt. 10e-5)) then
			    	 print*, "ERROR: real harmonic representation has non-zero imaginary part, STOPPING"
			    	 error stop
			    end if
			 end if
			 CReal(:,:) = REAL(C)

	     !Save eigenvector coefficients and eigenvalues
	     coeffs(:,:) = 0.0d0
	     eigs(:) = 0.0d0
	     call project(CReal, numfuncs, coeffs, eigs) 

			 !Correct sign problem 
			 do m = 1, numfuncs
			    if (sum(coeffs(:,m)) .lt. 0.0d0) then
				     coeffs(:,m) = (-1.0d0)*coeffs(:,m)
				  end if
			 end do

			 !Set very small coefficients to zero to correct for precision loss
			 do m = 1, numfuncs
					do mu = 1, numfuncs
						 if (abs(coeffs(m,mu)) .lt. 10e-5) then
								coeffs(m,mu) = 0.0d0
						 end if
				  end do
			 end do

			 !Save expansion coefficients
			 self%coeffs(1:numfuncs,1:numfuncs,l) = coeffs(1:numfuncs,1:numfuncs)
				 
			 deallocate(C, CReal)

   end subroutine construct_rep_group



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: subdecomp
	 !Purpose: decomposes the subduced representation of a point group into irreps of
	 !         of a subgroup
	 !Def: subduced representation: restriction of a representation of a group G to 
	 !    a subgroup H < G
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine subdecomp(self, subgroup, subconjarray, amat)
			implicit none
			type(group):: self, subgroup
			integer, dimension(subgroup%numirreps):: subconjarray
			integer:: ii,jj, p
			integer:: a, order
			integer:: m
			integer, dimension(self%numirreps,subgroup%numirreps):: amat

			amat(:,:) = 0
			!Loop over point group G irreps and H irreps
			do ii = 1, self%numirreps
				 do jj = 1, subgroup%numirreps
						!#conjugacy classes = #irreps
						a = 0
						do p = 1, subgroup%numirreps 
							 !Location of the p'th conjugacy class of H within those of G
							 m = subconjarray(p)

						   a = a + subgroup%conorder(p)*subgroup%chartable(jj,p)*self%chartable(ii,m)
						end do
						order = sum(subgroup%conorder(:))
						a = a/order
						amat(ii,jj) = a
						!print*, "(mu,v) = ", self%irrlabel(ii), subgroup%irrlabel(jj),  a
	       end do
	    end do

	 end subroutine subdecomp



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: project
	 !Purpose: diagonalises the matrix representation of the operator C_II, used to
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
			 coeffs(:,:) = CMat(:,:)

			 !deallocate(id)
       deallocate(work, iwork, ifail, w)

   end subroutine project

	 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: axistoeuler
	 !Purpose: computes the euler angles of a rotation given its axis-angle representation
	 !         nhat=(theta, phi), angle
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
