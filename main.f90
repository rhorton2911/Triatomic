!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program: H3Plus
!Purpose: performs a structure calculations for one and two electron
!         H3 molecule.
!Author: Reese Horton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!GLOBAL VARIABLES.... GLOBAL VARIABLES.... GLOBAL VARIABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------->>>>>GLOBAL VARIABLES<<<<<-----------------
!Due to the MCCC code inheriting features and/or programming
!constructs from legacy fortran77 code, there are several 
!instances of globally defined data structures being used.
!Most often, these are declared in the preamble of a particular
!module, which is then USED in the subroutine that makes use of the
!data structure. The following is a list of such variables directly used 
!in the present version of the MCCC code. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!dpf:: format specifier for declaring double precision floats in modern fortran
!   - located in numbers.f90
!data_in:: contains program input data read in from data.in file.
!   - located in input_data.f90
!   - initialisation subroutine readin uses global variables defined in
!     target_data.f90 and input_parser.f90.
!grid:: instance of rgrid type, stores grid of radial points and integration weights, etc
!   - located in gridr.f90
program H3Plus
   use grid_radial
   use sturmian_class
   use state_class
   use input_data     !Defines global data_in variable
	 use target_states  !Defines TargetStates global variable, TargetStates1el, etc.
	 use one_electron_func_group
   use basismodule
   use numbers
	 use vnc_module     !Stores nuclear potential vnc
   use ieee_arithmetic
   implicit none

	 type(basis_sturmian_nr):: basis 
	 !type(smallinput)::indata
   
   integer:: nr  !number of radial grid points
   integer:: ii, jj, kk
   integer:: n
   integer:: num_lambda
   real(dpf):: alphal, mu
   logical:: sorted, found
   real(dpf), dimension(:), allocatable:: func
   !Variables used to track program runtime
   character(len=8):: date1, date2
   character(len=10):: time1, time2
   character(len=5):: zone1, zone2
   integer, dimension(8):: values1, values2
   integer:: hr1, hr2, sec1, sec2, mins1, mins2
   integer:: hrs, mins, secs
   !Variables to control file io
   integer:: si, sj
	 integer:: harm


   !Intrinsic date and time subroutine
   call date_and_time(date1,time1,zone1,values1)

   !Legacy .f77 subroutines. Sets up common (global :( ) variables used by function YLM in plql.f
   call FAKRED
   call DFSET
 
   !Read in MCCC input data file to set up necessary data_in structure
   call readin( data_in, 10, 1)
   !Initialise data type containing rgrid and integration weights, global variable ( ): ) in grid_radial module
   call setgrids(grid)

	 !All data now in data.in MCCC input file
   !call readInput(indata)

   if (data_in%harmop .eq. 0) then
      !Molecule has C_s symmetry if:
      !(1)  R_2 != R_3
      !(2)  nuclei 2 and 3 are not symmetric about z-axis
      !(3)  nucleus 1 does not lie on z-axis (should never happen)
      !Matrix elements have non-zero complex part in C_s symmetry
      if (((abs(data_in%Rvec(2) - data_in%Rvec(3)) .gt. 0.0_dpf) .or. &
         (abs(abs(data_in%thetavec(2)) - abs(data_in%thetavec(3))) .gt. 0.0_dpf)) .or. &
         (abs(data_in%thetavec(1)) + abs(data_in%phivec(1))) .gt. 0.0_dpf) then
         print*, "WARNING: complex spherical harmonics not yet implemented for C_s geometries. Stopping."	
         error stop
      end if
   end if

	 !------------------------------ Set up nuclear potential ------------------------------!
	 nr = size(grid%gridr)
   call construct_vnc_group(nr,grid%gridr)

   !Set up nuclear potential stored in vnc_module, use formula: SUM_l=0^l=L (2l+1) = (L+1)^2
!	 nr = grid%nr
!   num_lambda = (indata%lambdamax+1)**2   
!   allocate(VPot(nr,num_lambda))
!   VPot(:,:) = 0.0_dpf
!   allocate(VPotTemp(nr, num_lambda))
!   do ii = 1, 3
!      call getVPotNuc(grid, VPotTemp, indata%R(ii), indata%theta(ii), &
!               indata%phi(ii), indata%charge(ii), indata)
!      VPot(:,:) = VPot(:,:) + VPotTemp(:,:)
!   end do
!   deallocate(VPotTemp)
! 
!   if (allocated(vnc)) then
!			deallocate(vnc)
!	 end if
!	 allocate(vnc(nr, num_lambda))
!	 vnc(:,:) = VPot(:,:)
!	 !lamtop_vc = data_in%ltmax
!	 lamtop_vc = num_lambda
!	 harm = indata%harmop
!	 data_in%harmop = harm
!	 deallocate(VPot)
!
   !---------------------Perform 2e Structure --------------------!
   !one_electron_func.f90: fills oneestates type before structure12
   !is called, specifically does so in construct_1el_basis_nr
   !and Hybrid_MSCbasis functions, which also set e1me and ovlpst
   !arrays
   !Hybrid_MSCbasis actually already loops over m when copying
   !laguerre functions to TargetStates
   !rearrange(): called in construct_1el_basis_nr, I've already
   !modified rearrange, might need to check its usage still works
   !when called in construct_1el_basis_nr

  
   !Modified version of subroutine in one_electron_func.f90, constructs basis for use in 2e configs
   call construct_1el_basis_nr_group(-1,basis)

   !Custom version of the structure12 subroutine tailored to non-linear molecules
   call structure12group(basis,TargetStates,basis%n)

   !--------------------End Two Electron Structure------------------!

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
