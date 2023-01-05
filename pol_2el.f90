module pol_2el
  
  public

  real*8:: gamma
  real*8:: rho
  real*8, dimension(:), allocatable:: pol
  integer:: i1_pol, i2_pol
  
contains
!
!
  subroutine pol_2el_initialize
    use grid_radial
    use input_data
    implicit none
  
    real*8 r
    integer:: i
    real*8:: expcut
  
    expcut = grid%expcut
    gamma = data_in%gamma
    rho = data_in%rho
    print*, gamma,rho
    
    allocate(pol(grid%nr))
          
    do i=1,grid%nr
       r = grid%gridr(i)
       pol(i) = polpot()
    end do

    call minmaxi(pol,grid%nr,i1_pol,i2_pol)
    
    pol(1:i1_pol-1) = 0d0

    write(*,'("Di-electron pol.pot.: gamma =",F9.5,", rho =",F9.5)') gamma, rho
    write(*,'("i1 =",I5,", i2 =",I5)') i1_pol, i2_pol
    !      write(*,'("Di-electron pol.pot.: gamma =",F9.5,", rho =",F9.5)') 
    !     >     gamma, rho
    !      write(*,'("i1 =",I5,", i2 =",I5)') i1, i2

contains
 function polpot()
    implicit none
    real*8:: polpot
    !
    real*8:: rat, rr0,vv, ww
    !
    rat = r/rho
    rr0 = rat*rat*rat
    if(rr0 .lt. -log(expcut)) then     
       vv =  exp(-rr0*rr0)
    else
       vv = 0.0
    end if
    ww = 1.0 - vv
    polpot = sqrt(ww)/(r*r)
    
  end function polpot


  end subroutine pol_2el_initialize
!
  subroutine pol_2el_destruct

    if(allocated(pol)) then       
       deallocate(pol)
    endif
  end subroutine pol_2el_destruct
!  
 


end module pol_2el



