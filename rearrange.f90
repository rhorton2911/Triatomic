!--------------------------------------------------------------------------------------------------------
subroutine rearrange(bst_nr,LMAX,TargetStates,use_lag_ham1el,basissize_nr,hamlimmin,hamlimmax,lag_ham1el_m,totfuncs,m_list_nr, sturm_ind_list_nr)
!!$ OneState is a pointer to one target state, and bst_nr is its single particle basis for  target states
!!$ will return back bst_nr which contains the new single particle basis (an array of sturmian_nr structure)
  
  use input_data
  use grid_radial
  use sturmian_class
  use state_class
  use MPI_module

  implicit none
  
  type(basis_sturmian_nr), intent(inout) :: bst_nr 
  integer, intent(in):: Lmax   ! max L valus in one-electron basis bst
  type(basis_state), intent(inout):: TargetStates
  logical, intent(in):: use_lag_ham1el ! switch use lag_ham1el_m - Laguerre basis <Lp|H|L> 
	integer, intent(in):: hamlimmin, hamlimmax, basissize_nr
  real*8, dimension(basissize_nr,basissize_nr,hamlimmin:hamlimmax), intent(in), optional :: lag_ham1el_m  ! 

  integer, dimension(TargetStates%Nmax,(LMAX+1)**2):: av_l_newf ! an array that will store  l values 
  integer, dimension(TargetStates%Nmax,(LMAX+1)**2):: av_n_newf ! an array that will store  n values 
  integer, dimension(TargetStates%Nmax,(LMAX+1)**2):: av_m_newf ! an array that will store  m values 
  integer, dimension(TargetStates%Nmax):: num_n_newf
  type(basis_sturmian_nr):: bst_rearr ! these are sturmian basis (depend on k and l and m)
  type(state), pointer:: OneState !( a pointer to) one target state only
  integer:: n
  integer:: l, lfn
  real*8, pointer, dimension(:):: fn, gn
  real*8:: rCI, alpha
  type(sturmian_nr), pointer:: pn, pnnew  !  one-electron orbitals
  integer:: nc, ncm
  integer:: maxl
  integer::minf,maxf, minf_all, maxf_all
  integer:: nst, nstmax
  integer:: tot_numf ! total number of new basis (function) to allocate (for all target states)
  integer:: icntr ! counter variable for indexing bst_rearr%b(bst_rearr)
  integer:: ictr ! counter variable for counting the number of basis to allocate
  integer:: i,j, i1, j1
  integer, dimension(:), allocatable:: ma_in, na_in
  real*8, dimension(:), allocatable:: CI_in
  integer:: nam_in, num ! num = num_n_newf(nst), updated in the loop over each target state
                        ! it is used to allocate no, na_in, CI_in with the size num
  real*8, dimension(grid%nr):: arrfnew, arrgnew
  logical :: spheroidal
  character(20):: filename
  integer:: icl, icheck, licl, mst, nspm_loc, nicl
!!$------------------------------------------------
  integer:: ii, jj, nsti, nstf, ncm_i, ncm_f, msti, mstf
  integer:: icl_i, icl_f, licl_i, licl_f, nicl_i, nicl_f, nci, ncf
  integer:: n_i, n_f, lfn_i, lfn_f
  real*8::  rCI_i, rCI_f
  real*8:: sum_new,sum_two
  type(state), pointer:: OneState_i, OneState_f
  type(sturmian_nr), pointer:: pn_i, pn_f
  integer:: i_print_states
!!$-----------------------------------------------
  integer::m, mfn, mfn_i, mfn_f, micl, micl_i, micl_f

	integer:: totfuncs
	integer, dimension(totfuncs):: m_list_nr, sturm_ind_list_nr
	integer:: n_i_ind, n_f_ind
	integer:: sind

  tot_numf = 0 ! intialization
  spheroidal = data_in%calculation_type==2 .or. data_in%calculation_type==3
  nstmax = TargetStates%Nmax ! number of target states
!!$  if (myid == 0) print*,'nstmax = ', nstmax
  
  num_n_newf(:) = 0 ! initialiize to 0
  
!!$ Initialize av_l_newf array elements with -1000, can never be reached
  av_l_newf(:,:) = -1000
  av_n_newf(:,:) = -1000
  av_m_newf(:,:) = -1000
 
  do nst = 1, nstmax 

     ictr = 0 ! a counter when a new l and m is encountered in the loop over all the Laguerre basis of a particular target state, initialize ictr to zero in every target state

     OneState => TargetStates%b(nst)

     ncm = get_nam(OneState) ! ncm = nam which is the size of the array na(:), which is also the number of CI coeff
     do nc = 1,ncm ! go through every Lagueere basis (psi_kl(r)X_lm(theta,phi)) of that target state
        !!$---------------------------------------
        rCI = get_CI(OneState, nc)
				!print*, nst, rCI, get_ang_mom(bst_nr%b(sturm_ind_list_nr(OneState%na(nc)))), m_list_nr(OneState%na(nc))
        if(rCI .eq. 0.0) cycle
        !!$----------------------------------------
        n = get_na(OneState,nc,1) ! which is na(nc), where nc = 1,2,...ncm
				sind = sturm_ind_list_nr(n)
        
        pn => bst_nr%b(sind) 
        l = get_ang_mom(pn) 
	      !m = get_ang_mom_proj(pn)
				m = m_list_nr(n)

        icheck = 0  ! check if this value of  lm  is already counted
        do icl=1,ictr
           if(l .eq. av_l_newf(nst,icl)) then
	            if (m .eq. av_m_newf(nst,icl)) then
                 icheck = 1   ! found it counted
                 exit
	            end if
           endif
        enddo
        
        if(icheck .eq. 0) then !  found new lm value
           ictr = ictr + 1
           tot_numf = tot_numf + 1
           av_l_newf(nst,ictr) = l                         
           av_n_newf(nst,ictr) = tot_numf     
	         av_m_newf(nst,ictr) = m
        endif

     end do ! end nc loop
     num_n_newf(nst) = ictr ! record the total number of different lm basis for that target state

  end do ! end nst loop
  
  if (myid == 0) write(*,'("rearrange.f90:  number of sp orbitals: tot_numf = ",I5)') tot_numf

  !!$---------------------------------------------
  call new_basis_nr(bst_rearr,tot_numf) ! allocate array of structure sturmian_nr (for the new basis)
  ! for all the target states, with tot_numf as the array (basis) size
  nspm_loc = 0
  do nst = 1, nstmax ! go through every target state
     OneState => TargetStates%b(nst)
     !Angular momentum projection is not a good quantum number for H3+, loop over m
     !mst = get_ang_mom_proj(OneState)
     ncm = get_nam(OneState)

     !alpha = get_alpha( bst_nr%b(get_na(OneState,1)) )   ! Get alpha from the first sturmian.     
     
     do icl=1, num_n_newf(nst) ! go through every lm of the target state
        licl = av_l_newf(nst, icl)
        nicl = av_n_newf(nst, icl)
        micl = av_m_newf(nst, icl)
				
        nspm_loc = nspm_loc + 1

        if(nicl .ne. nspm_loc) then
           if (myid == 0) print*, '*** rearrange.f90'
           if (myid == 0) print*, '*** nspm_loc .ne. nicl:', nspm_loc, nicl
           stop
        endif

        arrfnew(:) = 0d0 ! all elements are initialized to zero
        arrgnew(:) = 0d0   ! Analytical derivative for the spheroidal case.
        
        maxf_all = 1
        minf_all = grid%nr
        
        do nc = 1, ncm 
           n = get_na(OneState,nc,1)
					 !if (data_in%good_m) then
					 ! 	sind = n
					 !else
					 !   sind = sturm_ind_list_nr(n)
				   !end if
					 sind = sturm_ind_list_nr(n)

           pn => bst_nr%b(sind) 
           lfn = get_ang_mom(pn)
           !mfn = get_ang_mom_proj(pn)
					 mfn = m_list_nr(n)

           if(lfn .eq. licl) then
	            if (mfn .eq. micl) then
                 fn => fpointer(pn) ! means fn points to bst_nr%b(sind)%f
                 alpha = get_alpha(pn)   ! Get alpha, will be the same for each (lm) pair
                 minf = get_minf(pn)
                 maxf = get_maxf(pn)
                 rCI = get_CI(OneState,nc)

                 arrfnew(minf:maxf) =  arrfnew(minf:maxf) + rCI * fn(minf:maxf)
                 if (spheroidal) then
                    gn => gpointer(pn)
                    arrgnew(minf:maxf) = arrgnew(minf:maxf) + rCI * gn(minf:maxf)
                 end if
                 
                 maxf_all = max(maxf_all, maxf)
                 minf_all = min(minf_all, minf)
	            end if
           end if
        end do ! nc

        !!$ copy arrfnew array for that particular lfn into bst_rearr%b(icntr)
        if (.not.spheroidal) then
           call init_function(bst_rearr%b(nspm_loc),licl,micl,nspm_loc,alpha,minf_all,maxf_all,arrfnew,grid%nr)
        else
        !   call init_function(bst_rearr%b(nspm_loc), icl,licl,mst,alpha, minf_all,maxf_all,arrfnew,arrgnew)
           print*, "ERROR: spheroidal coordinates not yet implemented for H3+ mode. Stopping."
           error stop
        end if

     end do   
  end do ! nst
  
  if(nspm_loc .ne. tot_numf) then
     if (myid == 0) print*, '*** rearrange.f90'
     if (myid == 0) print*, '*** nspm_loc .ne. tot_numf:', nspm_loc, tot_numf 
     stop
  endif
!
!!$------------------------------------------------------------------------------------
!!$ code for ortint(:,:)
 
  !Liam added OMP to speed this up
!  !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(TargetStates, nstmax, num_n_newf, av_l_newf) &
!  !$OMP SHARED(av_n_newf, ave_m_newf, bst_rearr, bst_nr, lag_ham1el_m) 
  do nsti = 1, nstmax
     OneState_i => TargetStates%b(nsti)
     ncm_i = get_nam(OneState_i)
     !msti = get_ang_mom_proj(OneState_i)
 
     do icl_i  = 1, num_n_newf(nsti) ! walk through all lm for that target state
        licl_i = av_l_newf(nsti, icl_i)
        nicl_i = av_n_newf(nsti, icl_i)
	      micl_i = av_m_newf(nsti, icl_i)

           do nstf = 1, nstmax
              OneState_f => TargetStates%b(nstf)
              ncm_f = get_nam(OneState_f)
              !mstf = get_ang_mom_proj(OneState_f)
 
	       !M is no longer a good quantum number
!              !!$ check the same M-state values
!              if(msti .ne. mstf) then
!                 do  icl_f  = 1, num_n_newf(nstf) ! walk through all l for that target state
!                    licl_f = av_l_newf(nstf, icl_f)
!                    nicl_f = av_n_newf(nstf, icl_f)
!
!                    bst_rearr%ortint(nicl_f,nicl_i) = 0d0
!                    bst_rearr%ham1el(nicl_f,nicl_i) = 0d0
!
!                 end do
!                 cycle  ! over this value of nstf
!              end if
!

              do icl_f  = 1, num_n_newf(nstf) ! walk through all lm for that target state
                 licl_f = av_l_newf(nstf, icl_f) !!
                 nicl_f = av_n_newf(nstf, icl_f) !!
		             micl_f = av_m_newf(nstf, icl_f)

                 if(nicl_f .gt. nicl_i) cycle

                 sum_new = 0d0
                 sum_two = 0d0

                 do nci = 1, ncm_i
                    n_i = get_na(OneState_i,nci, 1)
										n_i_ind = sturm_ind_list_nr(n_i)
										!if (data_in%good_m) then
									  !   n_i_ind = n_i
									  !else
										!   n_i_ind = sturm_ind_list_nr(n_i)
								    !end if
                    pn_i => bst_nr%b(n_i_ind)
                    lfn_i = get_ang_mom(pn_i)
		                !mfn_i = get_ang_mom_proj(pn_i)
										mfn_i = m_list_nr(n_i)

                    if(lfn_i .ne. licl_i) cycle
	                  if (mfn_i .ne. micl_i) cycle 
                    rCI_i = get_CI(OneState_i, nci)

                    do ncf = 1, ncm_f
                       n_f = get_na(OneState_f,ncf, 1)
											 n_f_ind = sturm_ind_list_nr(n_f)
										   !if (data_in%good_m) then
									     !   n_f_ind = n_f
									     !else
										   !   n_f_ind = sturm_ind_list_nr(n_f)
								       !end if
                       pn_f => bst_nr%b(n_f_ind)
                       lfn_f = get_ang_mom(pn_f)
		                   !mfn_f = get_ang_mom_proj(pn_f)
											 mfn_f = m_list_nr(n_f)

                       if(lfn_f .ne. licl_f) cycle
		                   if (mfn_f .ne. micl_f) cycle
                       rCI_f = get_CI(OneState_f, ncf)

                       !sum_new = sum_new + rCI_i * rCI_f * bst_nr%ortint(n_i,n_f)
                       sum_new = sum_new + rCI_i * rCI_f * bst_nr%ortint(n_i_ind,n_f_ind)
                       if ( use_lag_ham1el ) then
                          sum_two = sum_two + rCI_i * rCI_f * lag_ham1el_m(n_i,n_f,0) 
                          !print*, "ERROR: ham1el_m option not yet implemented for H3+ calculation. Stopping"
                          !error stop
                       else
                          sum_two = sum_two + rCI_i * rCI_f * bst_nr%ham1el(n_i,n_f)
                       end if

                    end do ! ncf

                 end do ! nci

                 bst_rearr%ortint(nicl_f,nicl_i) = sum_new
                 bst_rearr%ortint(nicl_i,nicl_f) = sum_new

                 bst_rearr%ham1el(nicl_f,nicl_i) = sum_two
                 bst_rearr%ham1el(nicl_i,nicl_f) = sum_two
              end do ! icl_f

           end do ! nstf

        end do ! icl_i

     end do ! nsti
!     !$OMP END PARALLEL DO
!!$--------------------------------------------------------------------------------------
!!$ This is the code to overwrite the bst_nr (copying bst_rearr to bst_nr)
  call destruct(bst_nr) ! deallocate old bst_nr (old basis) first
  call new_basis(bst_nr, nspm_loc) ! allocate nspm_loc new space for a new created bst_nr
  
  call copy_basis(bst_nr,bst_rearr)
  !nspm =  basis_size(bst_nr)

  call destruct(bst_rearr)
!!$-------------------------------------------------------------------------------------
!!$ Check if it's correct
!   do ii = 1, nspm
!     do jj = 1, nspm
!        if (myid == 0) print*,'bst_nr%ortint(',ii,',',jj,') = ', bst_nr%ortint(ii,jj)
!     end do
!  end do
!!$ -------------------------------------------------------------------------------
!!$ This is code to update the TargetStates links to the new basis
  
  icntr = 1
  do nst = 1, nstmax ! go through every target state

     OneState => TargetStates%b(nst)

     mst = get_ang_mom_proj(OneState)
     num = num_n_newf(nst)
     
     allocate(na_in(num), CI_in(num), ma_in(num)) 
     
     nam_in = num
     na_in(1:num) = av_n_newf(nst,1:num)
     CI_in(1:num) = 1d0   !Ci coeffs absorbed into function f

		 !ma_in(1:num) = mst
     ma_in(1:num) = av_m_newf(nst,1:num)   !Only those (lm) with non-zero CI's are recorded, in order
     
     call modify_CI_state(OneState,nam_in,na_in,CI_in, ma_in)
     
     deallocate(ma_in)
     deallocate(na_in)
     deallocate(CI_in) 
     
  end do ! nst
!!$ ------------------------------------
!!$ print to file the pseudo state

  if  (myid == 0) then
     i_print_states = 0
     if(i_print_states .eq. 1) then
        filename = 'pseudostate_1s'
        nst = 1
        call Print_pseudo_state(nst, filename)
        
        if(nstmax .ge. 2) then
           filename = 'pseudostate_2p'
           nst = 2
           call Print_pseudo_state(nst, filename)
        endif
     endif
  end if

end subroutine rearrange

!!$ -----------------------------------------------------------------------------------

subroutine Print_pseudo_state(nst, filename)

  use sturmian_class
  use target_states
  use one_electron_func
  use input_data
  use grid_radial

  implicit none

  type(state), pointer:: OneState !( a pointer to) one target state only
  type(sturmian_nr), pointer:: pn
  integer:: n, l, m
!!$    integer:: nst, nstmax
  integer:: nst
  real*8, pointer, dimension(:):: f
  real*8::rCI
  real*8:: rval, tmp1, tmp2
  integer:: nc, ncm
  integer::minf,maxf
  integer:: i
  real*8, dimension(grid%nr):: arrf
  real*8:: u1s, u2s, u2p
  character(20):: filename
  
!!$  nstmax = TargetStates%Nmax ! number of target states
!!$  do nst = 1, nstmax ! go through every target state

  OneState => TargetStates%b(nst)
  ncm = get_nam(OneState)
  arrf(:) = 0d0 ! all elements are initialized to zero
  minf = 1; maxf = grid%nr

  do nc = 1, ncm ! number of new basis

     n = get_na(OneState,nc,1)
     pn => bst_nr%b(n)
     l = get_ang_mom(pn)
	   m = get_ang_mom_proj(pn)
     f => fpointer(pn) ! means f points to bst_nr%b(n)%f
     minf = get_minf(pn)
     maxf = get_maxf(pn)
!!$ call minmaxi(arrfnew,grid%nr,i1,i2) ! no need if minf starts from 1
     rCI = get_CI(OneState,nc)

     arrf(minf:maxf) =  arrf(minf:maxf) + rCI*f(minf:maxf)

  end do ! nc
!!$  end do ! nst

!!$ Now, print to file

  open(1289, file=filename)
  
  do i=minf, maxf
     rval =  grid%gridr(i)
     u1s = 4d0*Sqrt(2d0)*exp(-2d0*rval)*rval
     u2s = 2d0*rval*(1d0-rval)*exp(-rval)
     u2p = rval*(2d0**(5d0/2d0))*(24**(-1d0/2d0))*exp(-rval)*rval
     if(nst .eq. 1) then
        write(1289,*) rval, arrf(i), u1s
     elseif(nst .eq. 2) then
        write(1289,*) rval, arrf(i), u2p
     end if

  enddo

  close(1289)

end subroutine  Print_pseudo_state
