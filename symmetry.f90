!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Module: symmetry.f90
!Purpose: contains subroutines for performing the symmetry 
!         adaptation of a basis for a point group using the
!         eigenfunction method. 
!Note: makes use of the fact that irreps are in 1-1 correspondence with
!      eigenvalues of a complete set of commuting operators (CSCO), that
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
       character(len=6):: groupname
			 integer:: numirrs

			 if (trim(groupname) .eq. "D3h") then
			    numirrs = 6
			    self%numirreps = numirrs
			    allocate(self%irrlabel(numirrs), chartable(numirrs,numirrs), dims(numirrs), lamtable(numirrs,numirrs))
					allocate(self%condims(numirrs))
					self%irrlabel(1) = 'A1p'
					self%irrlabel(2) = 'A2p'
					self%irrlabel(3) = 'Ep' 
					self%irrlabel(4) = 'A1pp' 
					self%irrlabel(5) = 'A2pp' 
					self%irrlabel(6) = 'Epp' 
					self%chartable(1,:) = /(1, 1, 1, 1, 1, 1)/
					self%chartable(2,:) = /(1, 1, -1, 1, 1, -1)/
					self%chartable(3,:) = /(2, -1, 0, 2, -1, 0)/
					self%chartable(4,:) = /(1, 1, 1, -1, -1, -1)/
					self%chartable(5,:) = /(1, 1, -1, -1, -1, 1)/
					self%chartable(6,:) = /(2, -1, 0, -2, 1, 0)/
					self%dims(:) = /(1, 1, 2, 1, 1, 2)/
					self%conorder(:) = /(1, 2, 3, 1, 2, 3)/

					do jj = 1, numirrs  !classes
					   self%lamtable(:,jj) = self%chartable(:,jj)*self%conorder(jj)/self%dims(:)
				  end do
			 else
					print*, "ERROR: groups other than D3h not yet implemented, stopping."
					stop
			 end if
   end subroutine init_group


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: construct_rep
	 !Purpose: calculates the matrix representation of the operator C_II using the
	 !         wigner-d matrices and characters for each group element.
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine construct_rep(self, subgroup)
			 implicit none
			 !Data for the group and subgroup used to distinguish components of degenerate irreps
			 type(group):: self, subgroup  








   end subroutine construct_rep



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Subroutine: project
	 !Purpose: computes diagonalises the matrix representation of the operator C_II, used to
	 !         project a space of functions onto group irreps.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 subroutine project(CMat, spacedim, self, coeffs)
	     implicit none
			 integer:: spacedim   !Size of the basis set to be projected
			 real(dpf), dimension(spacedim,spacedim):: CMat
			 real(dpf), dimension(spacedim,spacedim):: coeffs 
			 type(group):: self
       !Required inputs for lapack routine dsygvx
       real(dpf):: abstol
       real(dpf), dimension(:), allocatable:: w, work
       integer, dimension(:), allocatable:: iwork, ifail
       integer:: info, lwork, nfound
       real(dpf), external :: DLAMCH

       allocate(work(1))
       allocate(w(spacedim))    !Stores eigenvalues
       allocate(iwork(5*spacedim))
       allocate(ifail(spacedim))
			 if (allocated(coeffs)) then
					deallocate(coeffs)
			 end if
			 allocate(coeffs(spacedim,spacedim))
       lwork = -1
       abstol = 2.0_dpf*DLAMCH('S') !twice underflow threshold for doubles
       info = 0
       nfound = 0
        		 
       !Workspace query with lwork=-1
       call DSYGVX(1, 'V', 'I', 'U', spacedim, CMat, spacedim, b, spacedim,    &
       0.0_dpf, 0.0_dpf, 1, spacedim, abstol, nfound, w, coeffs, spacedim, work, lwork, &
       iwork, ifail, info)
       
       lwork = int(work(1))
       deallocate(work)
       allocate(work(lwork))
       
       !Perform diagonialisation
       call DSYGVX(1, 'V', 'I', 'U', spacedim, CMat, spacedim, b, spacedim,    &
       0.0_dpf, 0.0_dpf, 1, spacedim, abstol, nfound, w, coeffs, spacedim, work, lwork, &
       iwork, ifail, info) 
       deallocate(work, iwork, ifail, w)
 

   end subroutine project

	 







end module symmetry
