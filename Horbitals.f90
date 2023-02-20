function H1el_Lag(ni,nj,mi,mj)
!!$ Calculates the one-electron molecular Hamiltonian matrix elements using a Laguerre basis.
!!$ Used especially for Laguerre functions with varying exponential fall offs.
!!$ Matrix of the form <i|H|j>. Arrays are of form (i,j).
!!$ This routine needs to be called before the rearrangement of H2+ ground state.
  
  use sturmian_class
  use one_electron_func
  use grid_radial
  use MPI_module
  
  implicit none
  
  real*8:: H1el_Lag
  integer, intent(in):: ni, nj     ! Indexes to Laguerre functions 
  integer, intent(in):: mi, mj         ! Angular projection of Laguerre functions
	integer:: ma
  
  real*8, dimension(:,:), pointer:: ortint
  real*8:: VLambdaR_ME_group
  type(sturmian_nr), pointer:: pi, pj           ! Laguerre basis pointers
  real*8, pointer, dimension(:):: fi, fj        ! Laguerre Basis functions
  real*8:: KE, temp_int, tmpovlp
  real*8:: alpha_i, alpha_j
  integer:: ki, kj, li, lj                      ! Laguerre funtion quatum numbers
  integer:: r1, r2, max_ri, max_rj, min_ri, min_rj
 
  ortint => bst_nr%ortint
  H1el_Lag = 0d0
  
  pi => bst_nr%b(ni)
  fi => fpointer(pi)
  ki = get_k(pi)
  li = get_ang_mom(pi)   
  alpha_i = get_alpha(pi)
  
  pj => bst_nr%b(nj)
  fj => fpointer(pj)
  kj = get_k(pj)
  lj = get_ang_mom(pj)
  alpha_j = get_alpha(pj)
  if (abs(mi) > li  .OR.  abs(mj) > lj)  return
  
!!$ One-electron Kinetic Energy Operator, diagonal in l,m
  KE = 0d0
  if (( li == lj ) .and. (mi .eq. mj)) then
     tmpovlp = bst_nr%ortint(ni,nj)
     
!!$ Kinetic energy of Laguerre functions with same alpha is done analytically.
!!$ Kinetic energy of Laguerre functions with different alpha is done 
!!$ numerically using analytic properties. 
     if ( alpha_i == alpha_j ) then
        
        KE = - alpha_j * alpha_j * tmpovlp/2d0
        if ( ki == kj ) then
           KE = KE + alpha_j * alpha_j
        end if
        
     else
        temp_int = 0d0
        KE = - alpha_j * alpha_j * tmpovlp/2d0                
        
        min_ri = get_minf(pi)
        min_rj = get_minf(pj)
        max_ri = get_maxf(pi)
        max_rj = get_maxf(pj)
        r1 = max(min_ri,min_rj)
        r2 = min(max_ri,max_rj)
        
        !$$ Analytic 1/x term now needs to be taken numerically
        temp_int = SUM( fi(r1:r2) * fj(r1:r2) * grid%weight(r1:r2) / grid%gridr(r1:r2) )
        
        KE = KE +  alpha_j * dble(kj + lj) * temp_int
        
     end if ! alpha_i == aplha_j             
  end if ! li == lj and mi==mj
  

  ! One-electron Hamiltonian for each Laguerre function 
  H1el_Lag =  KE + VLambdaR_ME_group(pi,pj,mi,mj) 

  
end function H1el_Lag

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!Function: VLambdaR_ME_group
    !Prupose: Evaluates V matrix elements between two basis functions with sturmian 
		!         parts pi, pj and angular momentum projections mi, mj. Uses predefined 
		!         nuclear potential V_{lm}(r) stored in vnc(:,:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function VLambdaR_ME_group(pi,pj,mi,mj)
      use input_data    !pinum 
      use grid_radial
      use sturmian_class
      use vnc_module
      use MPI_module
      use state_class
			use basismodule
      implicit none

      real*8:: VLambdaR_ME_group
      type(sturmian_nr), intent(in):: pi, pj
      integer, intent(in):: mi, mj
			integer:: ma

      real*8, pointer, dimension(:):: fi, fj, weight, gridr
      type(state), pointer :: coreorb
      type(sturmian_nr), pointer :: sturmj, sturmjp
      
      integer:: minfi, maxfi, minfj, maxfj, li,lj, minf, maxf, lam_min, lam_max, lcore, mcore, jj, na
      integer :: lcorep, mcorep, jjp, nap
      real*8:: tmpres, tmp, tmp1,res, ang1, ang2, CI, CIp
      real*8:: Yint
      integer:: lam, lamtmp, i1, i2, jcore, i1out, i2out
      real*8, dimension(grid%nr) :: invec, outvec
      integer, external :: OMP_GET_THREAD_NUM
      logical :: core

			integer:: lamind, q, qmin, qmax

      ma = mi

      VLambdaR_ME_group = 0d0

      weight => grid%weight

      tmpres = 0d0

      fi => fpointer(pi)
      minfi = get_minf(pi)
      maxfi = get_maxf(pi)
      li = get_ang_mom(pi)
             
      fj => fpointer(pj)
      minfj = get_minf(pj)
      maxfj = get_maxf(pj)
      lj = get_ang_mom(pj)

      minf = max(minfi,minfj)
      maxf = min(maxfi,maxfj)
      
      
!      if(data_in%iSlater .eq. 0) then
!         call fcexch_nr(pi,pj,res)
!         tmpres = tmpres + res
!      endif
      
      lam_min = abs(li-lj)
      lam_max = min(abs(li+lj),lamtop_vc)
      
      ! Can only equal lambda = 0,2,4,6,8,....
      if(data_in%good_parity .and. (-1)**lam_min .eq. -1) lam_min = lam_min + 1
      
      tmp = 0d0
			lamind = (lam_min)**2 !In case of loop over mu
      do lam=lam_min,lam_max   
         
         if(data_in%good_parity .and. (-1)**lam /= 1) cycle
         if(lam < abs(li-lj) .or. lam > li+lj) cycle
       
	 if (data_in%good_m) then
	    qmin = 0
            qmax = qmin
	 else
	    qmin = -lam
	    qmax = lam
         end if

	 !For non-diatomics, loop over lam, q
         do q = qmin, qmax
	    if (data_in%good_m) then
               tmp1 = Yint(dble(li),dble(mi),dble(lam),dble(0),dble(lj),dble(mj))
	       lamind = lam
	    else if (.not. data_in%good_m) then
	       !Running in H3+ mode, factor of 4pi/(2l+1) in partial wave expansion kept in Vlm in this case, have to compensate for
	       !Yint includeing a factor of sqrt(4pi/(2l+1))
	       if (data_in%harmop .eq. 0) then
                  tmp1 =  sqrt(dble(2*lam+1)/(4.0_dpf*pinum))*Yint(dble(li),dble(mi),dble(lam),dble(q),dble(lj),dble(mj))
               else if (data_in%harmop .eq. 1) then
                  tmp1 = Xint(dble(li),dble(mi),dble(lam),dble(q),dble(lj),dble(mj))
	       end if
	       lamind = lamind + 1
	    end if

            !  if (myid == 0) print*, i,j, lam, tmp1
            i1 = max(minvnc(lamind),minf)
            i2 = min(maxvnc(lamind),maxf)
            tmp = tmp + tmp1 * SUM(weight(i1:i2)*fi(i1:i2)*fj(i1:i2)*vnc(i1:i2,lamind))
	 end do

         if ((data_in%N_core_el > 0 .and. .not.data_in%pseudo_pot) .and. data_in%good_m) then !add core exchange ME
          do jcore=1, CoreOrbitals%Nmax
            
            coreorb => CoreOrbitals%b(jcore)
    
            do jj=1, coreorb%nam

              na = coreorb%na(jj)
              CI = coreorb%CI(jj)

              sturmj => wfcore_nr%b(na)

              lcore = sturmj%l 
              mcore = sturmj%m
               
              if(lam < abs(lcore-lj) .or. lam > lcore+lj) cycle
              ang1 = Yint(dble(lcore),dble(mcore),dble(lam),dble(mcore-ma),dble(lj),dble(ma))
              if(ang1 == 0.0d0) cycle
                
              i1 = max(minfj,sturmj%minf)
              i2 = min(maxfj,sturmj%maxf)

              invec(i1:i2) = fj(i1:i2) * sturmj%f(i1:i2)
              call form_accurate(lam,invec,i1,i2,grid%nr,outvec,i1out,i2out)

              do jjp=1, coreorb%nam
              
                nap = coreorb%na(jjp)
                CIp = coreorb%CI(jjp)

                sturmjp => wfcore_nr%b(nap)

                lcorep = sturmjp%l 
                mcorep = sturmjp%m

                if (mcorep /= mcore) error stop '*** ERROR in core exchange: mcorep /= mcore'

                if(lam < abs(lcorep-li) .or. lam > lcorep+li) cycle
                ang2 = Yint(dble(li),dble(ma),dble(lam),dble(mcore-ma),dble(lcorep),dble(mcorep))
                if(ang2 == 0.0d0) cycle
                
                i1 = max(i1out,sturmjp%minf)
                i1 = max(i1,minfi)
                i2 = min(i2out,sturmjp%maxf)
                i2 = min(i2,maxfi)

                tmp = tmp - ang1*ang2*CI*CIp*sum(fi(i1:i2)*sturmjp%f(i1:i2)*outvec(i1:i2)*weight(i1:i2))

              enddo !jjp
 
            enddo !jj

          enddo !jcore

        endif !core

      enddo !lam

      tmpres = tmpres + tmp 

      VLambdaR_ME_group = tmpres
      
    end function VLambdaR_ME_group




!!$  It is assumed that basis  bst  is made from Laguerre function of 
!!$  order  (2l+1)  that are orthogonal with weight 1/r
!!$  these functions have definite value of orbital angular momentum l
!!$ Subroutine requires Laguerre functions with the same value alpha
 subroutine Hlagorb(bst,ist,jst,ma,result)
    
    use sturmian_class
    use MPI_module
    use vnc_module, only: core_energy
    use input_data
    
    implicit none
    
    type(basis_sturmian_nr), intent(in) :: bst   ! this is Sturmian basis  
    integer, intent(in):: ist,jst
    integer, intent(in):: ma
    real*8, intent(out):: result
    
    type(sturmian_nr), pointer:: pi, pj
    integer:: li, lj
    real*8:: al, tmpovlp
    real*8::  VLambdaR_ME
    integer, external :: OMP_GET_THREAD_NUM
    
    result = 0d0
   
    pi => bst%b(ist)
    pj => bst%b(jst)
             
    li = get_ang_mom(pi)             
    lj = get_ang_mom(pj)

    if(abs(ma) .gt. li  .OR.  abs(ma) .gt. lj)  return
    if(lj .eq. li) then
       
       al = get_alpha(bst%b(ist))    ! it is assumed that it is Laguerre functions, and made for the same alpha forgiven l
       ! term with  alpha  appears only in matrix elements of the kinetic energy operator where li=lj


!!$   get kinetic energy: (using special properties of Lag func)
       tmpovlp = bst%ortint(ist,jst)

       result = -al*al*tmpovlp/2d0
       if(ist == jst) then
          result = result + al*al
       end if
    end if


!    if (myid == 0) print*, result, VLambdaR_ME(pi,pj,ma)
    result = result + VLambdaR_ME(pi,pj,ma)

!    if(myid==0) print*, 'VLAMBDA:', get_k_nr(pi),get_ang_mom(pi),get_k_nr(pj),get_ang_mom(pj),ma,VLAMBDAR_ME(pi,pj,ma)
           

  end subroutine Hlagorb
  
