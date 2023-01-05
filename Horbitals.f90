function H1el_Lag(ni,nj,ma)
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
  integer, intent(in):: ma         ! Angular porjection of matrix element   
  
  real*8, dimension(:,:), pointer:: ortint
  real*8:: VLambdaR_ME
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
  if (abs(ma) > li  .OR.  abs(ma) > lj)  return
  
!!$ One-electron Kinetic Energy Operator
  KE = 0d0
  if ( li == lj ) then
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
  end if ! li == lj
  
!!$ One-electron Hamiltonian for each Laguerre function 
  H1el_Lag =  KE + VLambdaR_ME(pi,pj,ma) 
  
end function H1el_Lag

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
  
