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
				 integer, dimension(:), allocatable:: irrdims !dimension of each representation
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



			 end if 

end module symmetry
