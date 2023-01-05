subroutine oscstr_2e()

  use input_data 
  use target_data
  use one_electron_func
  use sturmian_class
  use target_states

  implicit none  

  integer:: nstate, npstate, Nmax, ninitial_max         ! Loops and input option controls
  ! Results
  real*8::  result_l, result_v,  result_p, result_pc, ioniz_en_au, rl, rv
  real*8:: result_p_par, result_p_perp, pol_b, pol_c
  real*8:: TranEnergy,  temp_p, Energy, Osc, result_aee_l, result_aee_v, summedosc_l, summedosc_v, result_aee_ex_l, result_aee_ex_v, tmp_l, tmp_v
  integer:: map, ma, dM, Spin, Spinp, parity, parityp   ! Quantum Numbers Selection Rules
  type(state), pointer :: statef,statei
  integer :: mstf, msti, i
  character(len=:), allocatable :: flabel, ilabel
  character(len=3) :: state_f, state_i
  
  interface
    subroutine oscstrength_2e_MOrep(nstate,npstate,TranEnergy,result_l,result_v, rl, rv)
      use grid_radial
      use input_data
      use sturmian_class
      use target_states
      use one_electron_func
      use ovlpste1me
      use state_class
      integer, intent(in):: nstate, npstate
      real*8, intent(in):: TranEnergy
      real*8, intent(out):: result_l, result_v
      real*8, optional, intent(out):: rl, rv
    end subroutine
    subroutine  oscstrength_2e_config(nstate,npstate,TranEnergy,result_l,result_v,rl,rv)
      use grid_radial 
      use sturmian_class
      use target_states
      use one_electron_func
      use state_class
      integer, intent(in):: nstate, npstate               ! Molecular state index
      real*8, intent(in):: TranEnergy
      real*8, intent(out):: result_l, result_v
      real*8, optional, intent(out) :: rl, rv
    end subroutine
  end interface
  
  Nmax = basis_size_st(TargetStates2el)                 ! Number of Molecular States Generated

  ninitial_max = min(data_in%iosc,Nmax)

  print*,"*****************************************"
  print*," OSCILLATOR STRENGTHS AND POLARISABILITY"
  print*,"*****************************************"

  open(301,file='osc.strength')
  write(301,'(A)') 'TWO-ELECTRON OSCILLATOR STRENGTHS'
  write(301,'(A)') 'Literature solutions are from S.E.Branchett and J.Tennyson'
  write(301,'(A)') 'J.Phys.B: At.Mol.Opt.Phys. 25 (1992) 2017-2026'

  if(data_in%print_dipole) open(unit=12345,file='dipole.out',action='write')
  
  ! INITIAL Molecular State Loop
  do nstate = 1, ninitial_max                           
     ! Initial State number n
     ma = NINT(TargetStates2el%b(nstate)%M)       ! 2e Molecular State Angular Projection
     Spin = NINT(TargetStates2el%b(nstate)%spin)  ! 2e Molecular State Spin
     parity = TargetStates2el%b(nstate)%parity    ! 2e Molecular State Parity

     write(*,'(A,I3,A,3I3)') 'Initial state #',nstate, ' with Parity, M, Spin :', parity,ma,Spin
     write(301,'(A,I3,A,3I3)') 'Initial state #',nstate, ' with Parity, M, Spin :', parity,ma,Spin

     result_p = 0d0
     result_pc = 0d0
     result_p_par = 0d0
     result_p_perp = 0d0
     pol_b = 0d0
     pol_c = 0d0

     summedosc_l = 0d0
     summedosc_v = 0d0

     result_aee_l= 0d0
     result_aee_v = 0d0
     result_aee_ex_l= 0d0
     result_aee_ex_v = 0d0

     ! FINAL Molecular State Loop     
     do npstate = 1, Nmax   
        
        if(npstate == nstate) cycle 
        
        map = NINT(TargetStates2el%b(npstate)%M)          
        Spinp = NINT(TargetStates2el%b(npstate)%spin)
        parityp = TargetStates2el%b(npstate)%parity
        
        dM = map - ma        
        
        TranEnergy = get_energy(TargetStates2el%b(npstate)) -  get_energy(TargetStates2el%b(nstate))

        ioniz_en_au = real( get_energy(TargetStates2el%b(npstate)) -  TargetStates2el%en_ion) 
          
        result_l= 0d0
        result_v = 0d0
        
        ! MAKE SURE BOTH MOLECULAR STATES HAVE THE SAME SPIN and difference Angular Projections <= 1
        ! Parity must change
        
        if ( Spin /= Spinp ) cycle
        if (ABS(dM) > 1 ) cycle
        if ( parity == parityp .and. data_in%good_parity) cycle 

        call oscstrength_2e_MOrep(nstate,npstate,TranEnergy,result_l,result_v,rl,rv)

        statef => TargetStates2el%b(Npstate)
        statei => TargetStates2el%b(Nstate)
        
        flabel=trim(adjustl(get_label_st(statef)))
        ilabel=trim(adjustl(get_label_st(statei)))
        mstf=get_ang_mom_proj(statef)
        msti=get_ang_mom_proj(statei)
  
        write(state_i,'(I3.3)') Nstate
        write(state_f,'(I3.3)') Npstate
        
        if(flabel == '?') then
          flabel = state_f
        elseif(mstf < 0) then
          flabel = flabel//'-'
        endif
        if(ilabel == '?') then
          ilabel = state_i
        elseif(msti < 0) then
          ilabel = ilabel//'-'
        endif

        do 
          i = index(flabel,"'")
          if(i == 0) exit
          flabel = flabel(:i-1)//'p'//flabel(i+1:)
        enddo
        
        do 
          i = index(ilabel,"'")
          if(i == 0) exit
          ilabel = ilabel(:i-1)//'p'//ilabel(i+1:)
        enddo

        if(data_in%print_dipole) write(12345,*) flabel//'.'//ilabel, rl, rv

!!$        print*,"ni",nstate,"nf",npstate,"dM",dM
!!$        print*,"L",result_l,"V",result_v

        ! Write to output.
!!$        write(*,'(4A,4I4)') get_label(TargetStates2el%b(nstate)),' --->',get_label(TargetStates2el%b(npstate)), '    i,f, dM,Spin:', nstate,npstate,dM,Spin
!!$        write(*,'(2(A,Es20.12))') 'L:', result_l, '     V:', result_v
        write(*,'(4A,2(2I4,4X),2(A,Es20.12))') get_label(TargetStates2el%b(nstate)),' --->',get_label(TargetStates2el%b(npstate)), '    i,f, dM,Spin:', nstate,npstate,dM,Spin, '    L:',result_l, '    V:',result_v

        ! Also write to osc.strengths file.
!!$        write(301,'(4A,4I4)') get_label(TargetStates2el%b(nstate)),' --->',get_label(TargetStates2el%b(npstate)), '    i,f, dM,Spin:', nstate,npstate,dM,Spin
!!$        write(301,'(2(A,Es20.12))') 'L:', result_l, '     V:', result_v
        write(301,'(4A,2(2I4,4X),2(A,Es20.12))') get_label(TargetStates2el%b(nstate)),' --> ',get_label(TargetStates2el%b(npstate)), '    i,f, dM,Spin:', nstate,npstate,dM,Spin, '    L:',result_l, '    V:',result_v

        if ( Spin /= Spinp ) cycle
        if (ABS(dM) > 1 ) cycle
        
        ! Polarizabilty calulated from Oscillator strengths and also average excitation energy (see Inokuti RMP 43(1971)297 Eq. 4.62)
        temp_p = 0d0
        tmp_l = 0d0
        tmp_v = 0d0
        
        if(parity == -parityp) then

          if ( ABS(dM) == 1 ) then
             ! Summing over f_n dm = +/-1, therefore divide by 2
  !!$ Orbital degeneracy g=2 for m> 0 states. Included in osc routine.
  !!$ therefore divide by 2 
             temp_p =  result_l / ( 2.0 * TranEnergy * TranEnergy )
             result_p_perp =  result_p_perp + 1.5 * temp_p
             tmp_l = result_l / 2d0
             tmp_v = result_v / 2d0
          else if ( dM == 0 ) then
             temp_p =  result_l / ( TranEnergy * TranEnergy )
             result_p_par =  result_p_par + 3.0 * temp_p
             tmp_l = result_l 
             tmp_v = result_v 
          end if
          ! For continuum contribution
          if (  get_energy_st(TargetStates2el%b(npstate)) > 0 ) then
             result_pc = result_pc + temp_p
          end if
          result_p = result_p + temp_p            
          if(ioniz_en_au .gt.0) then
             pol_c = pol_c + temp_p
             if(data_in%print_pol_contribution) then
               print*, '*** CONTINUUM POL '//get_label(TargetStates2el%b(nstate))//' --->'//get_label(TargetStates2el%b(npstate)), &
                ' temp_p, pol_c:', temp_p, pol_c
              endif
          else
             pol_b = pol_b + temp_p
             if(data_in%print_pol_contribution) then
               print*, '*** BOUND POL '//get_label(TargetStates2el%b(nstate))//' --->'//get_label(TargetStates2el%b(npstate)), &
                 ' temp_p, pol_b:', temp_p, pol_b
             endif
          endif

        else
          if ( ABS(dM) == 1 ) then
             temp_p =  result_l / ( 2.0 * TranEnergy * TranEnergy )
          else if ( dM == 0 ) then
             temp_p =  result_l / ( TranEnergy * TranEnergy )
          end if
        endif

        summedosc_l = summedosc_l + tmp_l
        summedosc_v = summedosc_v + tmp_v
        
        result_aee_l = result_aee_l + tmp_l*log(2.0*TranEnergy)
        result_aee_v = result_aee_v + tmp_v*log(2.0*TranEnergy)
        result_aee_ex_l = result_aee_ex_l + tmp_l*TranEnergy
        result_aee_ex_v = result_aee_ex_v + tmp_v*TranEnergy

     end do
     if (nstate == 1) then  !TODO: properly format polarizability for excited states
        ! Oscillator Strengths: R-matrix S. Branchett 92
        ! Polarizability: W. Kolos 67  TABLE II R = 1.4
        write(*,*) "Polarizability: Ref[1]:  W. Kolos 67 R = 1.4"
        write(*,*) "Parallel:",result_p_par,"Ref[1]: 6.38049"
        write(*,*) "Perpendicular:",result_p_perp,"Ref[1]: 4.57769"
        ! a_T = 1/3 * a_par + 2/3 * a_perp
        write(*,*) "Total:",result_p,"Ref[1]: 5.17862"
        write(*,*) "Total bound: ",pol_b,", cont.: ", pol_c
        if(result_p>0) write(*,*) "Total bound(%): ", pol_b/result_p*100, ", cont.(%):", pol_c/result_p*100
        write(*,*)  "Summed osc.str. (l,v)", summedosc_l, summedosc_v

        if(summedosc_l > 0 .and. summedosc_v > 0) then
          result_aee_l = exp(result_aee_l /summedosc_l)/2d0
          result_aee_v = exp(result_aee_v /summedosc_v)/2d0
          result_aee_ex_l = result_aee_ex_l /summedosc_l
          result_aee_ex_v = result_aee_ex_v /summedosc_v
          write(*,*) "Average excitation energy(l,v) (equation in RMP and just ex.energy):",result_aee_l, result_aee_v, result_aee_ex_l, result_aee_ex_v
        endif

        ! Also write to osc.strengths file.
        write(301,*) "Polarizability: Ref[1]:  W. Kolos 67 R = 1.4"
        write(301,*) "Parallel:",result_p_par,"Ref[1]: 6.38049"
        write(301,*) "Perpendicular:",result_p_perp,"Ref[1]: 4.57769"
        write(301,*) "Total:",result_p,"Ref[1]: 5.17862"
        write(301,*)  "Summed osc.str. (l,v)", summedosc_l, summedosc_v
        write(301,*) "Average excitation energy(l,v) (equation in RMP and just ex.energy):",result_aee_l, result_aee_v, result_aee_ex_l, result_aee_ex_v
     end if
  end do
  if(data_in%print_dipole) close(12345)

  write(*,'(A//)') "*****************************************"
  write(301,*) '# # #'
  close(301)  

end subroutine oscstr_2e


subroutine  oscstrength_2e(nstate,npstate,TranEnergy,result_l,result_v)
! Calculates oscillator strengths in the form of Length <np1,np2|L|n1,n2>
! f_L =  sum C_(np1,np2) C(n1,n2) <np1,np2|L|n1,n2> 

  use grid_radial 
  use sturmian_class
  use target_states
  use one_electron_func
  use state_class


  implicit none  
  integer, intent(in):: nstate, npstate               ! Molecular state index
  real*8, intent(in):: TranEnergy
  real*8, intent(out):: result_l, result_v
  ! Quantum numbers 
  integer::  la1, la2, lap1, lap2, ma, map, mp1, mp2, m1, m2,  dL, dM, Spin, Spinp
  integer:: parity, parityp
  integer::  ne_con,  nep_con, Ncon, Npcon                               ! Configuration Loops
  type(sturmian_nr), pointer:: tn1, tn2, tnp1, tnp2                      ! One-electron orbitals
  integer:: ind1, ind2, indp1, indp2
  integer::  minf, maxf, minfp, maxfp,  ir1, ir2                         ! Radial integrations
  real*8:: Coeff, temp, CI, CIp, temp_l
  real*8:: temp_v, ideriv1,  vint1, r_overlap1, overlap2                 ! Length velocity Calcs
  real*8, pointer, dimension(:)::  f1, f2, fp1, fp2                      ! f One-electron functions
  real*8, pointer, dimension(:):: weight, gridr

 
  weight => grid%weight
  gridr => grid%gridr
  result_l = 0d0
  temp_l = 0d0
  result_v = 0d0
  temp_v = 0d0


! Final State number npstate
  Npcon =  TargetStates2el%b(npstate)%nam        ! Number of A.O. and Molecular ion configurations.        
  map = NINT(TargetStates2el%b(npstate)%M )      ! Molecular state Angular Projection
  Spinp = NINT(TargetStates2el%b(npstate)%spin)  ! 2e Molecular State Spin
  parityp = TargetStates2el%b(npstate)%parity    ! 2e Molecular State Parity

! Initial State number nstate
  Ncon = TargetStates2el%b(nstate)%nam         
  ma = NINT(TargetStates2el%b(nstate)%M)       
  Spin = NINT(TargetStates2el%b(nstate)%spin)
  parity = TargetStates2el%b(nstate)%parity     
  
  dM = map - ma 

! MAKE SURE BOTH MOLECULAR STATES HAVE THE SAME SPIN and difference Angular Projections <= 1
  ! Parity must change
  if ( Spin /= Spinp ) return
  if (ABS(dM) > 1 ) return
  if ( parity == parityp) return


! Looping over FINAL Molecular State orbitals.
  do nep_con =  1, Npcon

! Quantum numbers and functions for ALPHA
     indp1 = TargetStates2el%b(npstate)%na(nep_con)       ! Final state number np. nep A.O.
     tnp1 => bst_nr%b(indp1)                              !           
     fp1 => fpointer(tnp1)                                ! One electron functions
     lap1 = get_ang_mom(tnp1)                             ! Gets angular momentum from the one-electron orbitals
     mp1 =  TargetStates2el%b(npstate)%ma(nep_con)        ! Get angular projection of A.O.   

! Quantum numbers and functions for BETA
     indp2 = TargetStates2el%b(npstate)%nb(nep_con)       
     tnp2 => bst_nr%b(indp2)                                         
     fp2 => fpointer(tnp2)                                
     lap2 = get_ang_mom(tnp2)                             
     mp2 =  TargetStates2el%b(npstate)%mb(nep_con)           
     CIp = get_CI(TargetStates2el%b(npstate),nep_con)  

! looping over INITIAL Molecular State orbitals.
     do ne_con = 1, Ncon

! Quantum numbers and functions for GAMMA
        ind1 = TargetStates2el%b(nstate)%na(ne_con)       
        tn1 => bst_nr%b(ind1)                                          
        f1 => fpointer(tn1)                               
        la1 = get_ang_mom(tn1)                            
        m1 =  TargetStates2el%b(nstate)%ma(ne_con)           


! Quantum numbers and functions for DELTA
        ind2 = TargetStates2el%b(nstate)%nb(ne_con)       
        tn2 => bst_nr%b(ind2)                                         
        f2 => fpointer(tn2)                               
        la2 = get_ang_mom(tn2)                            
        m2 =  TargetStates2el%b(nstate)%mb(ne_con)           
        CI = get_CI(TargetStates2el%b(nstate),ne_con)  
        
        dL = lap1 - la1
        dM = mp1 - m1
        
        ! Selections Rules 
        if ( ABS(dL) > 1 .OR. ABS(dM) > 1 ) cycle
       
        if ( lap2 /= la2 .OR. mp2 /= m2  ) cycle
        ! DO OVERLAPOF COORDINATE SPACE 2  
        minf = get_minf(tn2)
        maxf = get_maxf(tn2)
        minfp = get_minf(tnp2)
        maxfp = get_maxf(tnp2)      
        ir1 = max(minf,minfp)
        ir2 = min(maxf,maxfp)  
        overlap2 = SUM( fp2(ir1:ir2) *  f2(ir1:ir2) * weight(ir1:ir2) )

        ! COORDINATE 1 RADIAL INTEGRALS
        minf = get_minf(tn1)
        maxf = get_maxf(tn1)
        minfp = get_minf(tnp1)
        maxfp = get_maxf(tnp1)
        ir1 = max(minf,minfp)
        ir2 = min(maxf,maxfp) 
        ! Radial Integration <fp1|r|f1>
         r_overlap1 = SUM( fp1(ir1:ir2) *  f1(ir1:ir2)  * gridr(ir1:ir2) * weight(ir1:ir2) ) 
        
        !Velocity Integrals
        ! Derivatiie Integral Note: d/dr(phi(r)/r) = phi'(r)/r -  phi(r)/(r^2)
        ideriv1 =  SUM((( f1(ir1 + 1 : ir2 ) - f1(ir1: ir2 - 1)) /(gridr(ir1 + 1 :ir2) - gridr(ir1: ir2 - 1))) * fp1(ir1: ir2 - 1) * weight(ir1:ir2 - 1))
        vint1 = SUM( fp1(ir1:ir2) * f1(ir1:ir2) * weight(ir1:ir2) / gridr(ir1:ir2) )


        if ( la1 < lap1 ) then
           vint1 = ideriv1 - (dble(la1) + 1.0) * vint1
        else if ( la1 > lap1) then
           vint1 = ideriv1 + dble(la1) * vint1
        end if
                    

        temp = (2.0 * la1 + 2.0 + dL) * (2.0 * la1 + dL)
        Coeff = 0d0 ! Length Guage Block Varsholovich pg 145 
        if ( dM == 0) then ! Parallel Transitions <z> 
           Coeff = sqrt( ( dble(la1) + 0.5 + dble(dL) * 0.5 ) ** 2.0 - dble(m1 * m1)) 
        end if
        if ( ABS(dM) == 1) then
           if ( la1 < lap1   ) then  ! Perpendicular Transitions <x +/- iy> 
              Coeff = dble(dM)*sqrt(dble(la1+dM*m1+2)*dble(la1+dM*m1+1))
           else if ( la1 > lap1) then
              Coeff = -dble(dM)*sqrt(dble(la1-dM*m1)*dble(la1-dM*m1-1))
           end if
           Coeff = Coeff /  sqrt(2.0)                 
        end if
              
        Coeff = Coeff / sqrt(temp) 
              
        ! Multiply by 2 for 2e target. Operator L: z_1 + z_2 = 2*z_1
        result_l = temp_l + 2.0 * CIp * CI * Coeff * r_overlap1 * overlap2
        result_v = temp_v + 2.0 * CIp * CI * Coeff * vint1 * overlap2
        temp_l = result_l
        temp_v = result_v


     end do ! Initial States Orbitals
     
  end do    ! Final States Orbitals

  dM = map - ma 
! TranEnergy is in a.u. converted to Rydbergs here.
  if ( ABS(dM) == 1) then
     result_l = 4.0 * ABS(TranEnergy) * result_l * result_l / 3.0 ! Added orbital degeneracy g = 2, should test for ma=+/-2 etc 
     result_v = 4.0 * result_v * result_v / (3.0 * ABS(TranEnergy))
  else 
     result_l = 2.0 * ABS(TranEnergy) * result_l * result_l / 3.0 
     result_v = 2.0 * result_v * result_v / ( 3.0  * ABS(TranEnergy) )
  end if

 
end subroutine oscstrength_2e

subroutine  oscstrength_2e_MOrep(nstate,npstate,TranEnergy,result_l,result_v, rl, rv)
! Calculates oscillator strengths in the form of Length <np1,np2|L|n1,n2>
! f_L =  sum C_(np1,np2) C(n1,n2) <np1,np2|L|n1,n2> 

  use grid_radial
  use input_data
  use sturmian_class
  use target_states
  use one_electron_func
  use ovlpste1me
  use state_class
  implicit none

  integer, intent(in):: nstate, npstate               ! Molecular state index
  real*8, intent(in):: TranEnergy
  real*8, intent(out):: result_l, result_v
  real*8, optional, intent(out):: rl, rv
  ! Quantum numbers 
  integer:: la1, la2, lap1, lap2, ma, map, mp1, mp2, m1, m2,  dL, dM, Spin, Spinp
  integer:: parity, parityp
  integer:: n_orb, np_orb, n_orb_max, np_orb_max        ! loops nuse 
  integer:: use_norb, use_nporb                         ! index to 1e molecular orbital
  integer:: nst_con,  nstp_con, Ncon, Npcon             ! 2e configuration loops
  integer:: MO1, MO2, MOp1, MOp2                        ! index to 1e molecular orbital
  integer:: i1, i2, ip1, ip2, nam1, nam2, namp1, namp2  ! Loop over 1e molecular orbitals atomic orbitals

  logical :: spheroidal
  type(sturmian_nr), pointer:: tn1, tn2, tnp1, tnp2                      ! One-electron orbitals
  integer:: ind1, ind2, indp1, indp2
  integer:: minf, maxf, minfp, maxfp,  ir1, ir2                         ! Radial integrations
  real*8:: CI_stp, CI_st, CI_1, CI_2, CI_1p, CI_2p, temp_CI ! 2e configuration and 1e molecular orbial coefficients
  real*8:: Coeff,  temp, CI, CIp, temp_l
  real*8:: temp_v, ideriv1,  vint1, r_overlap1, overlap2                 ! Length velocity Calcs
  real*8, pointer, dimension(:)::  f1, fp1                      ! f One-electron functions
  real*8, pointer, dimension(:):: weight, gridr
  real*8, dimension(1:grid%nr) :: lVec,lVecTmp, vVec,vVecTmp

  weight => grid%weight
  gridr => grid%gridr
  result_l = 0d0
  temp_l = 0d0
  result_v = 0d0
  temp_v = 0d0

  lVec(:) = 0d0; vVec(:) = 0d0
  spheroidal = data_in%calculation_type==2 .or. data_in%calculation_type==3

! Final State number npstate
  Npcon =  TargetStates2el%b(npstate)%nam        ! Number of A.O. and Molecular ion configurations.        
  map = NINT(TargetStates2el%b(npstate)%M )      ! Molecular state Angular Projection
  Spinp = NINT(TargetStates2el%b(npstate)%spin)  ! 2e Molecular State Spin
  parityp = TargetStates2el%b(npstate)%parity    ! 2e Molecular State Parity
  np_orb_max = TargetStates2el%b(npstate)%nusemax  ! Number of orbitals used to describe this state
! Initial State number nstate
  Ncon = TargetStates2el%b(nstate)%nam
  ma = NINT(TargetStates2el%b(nstate)%M)
  Spin = NINT(TargetStates2el%b(nstate)%spin)
  parity = TargetStates2el%b(nstate)%parity
  n_orb_max = TargetStates2el%b(nstate)%nusemax

  dM = map - ma

! MAKE SURE BOTH MOLECULAR STATES HAVE THE SAME SPIN
  if ( Spin /= Spinp ) return
  
  if (ABS(dM) > 1 ) return

  if(data_in%good_parity) then
    if(parity /= parityp .and. ABS(dM) > 1) return !Dipole selection rule
  endif

  do np_orb = 1, np_orb_max
!!$ nuse defined in state_class will refer to one-electron target states
!!$ so can be used! 
     use_nporb = TargetStates2el%b(npstate)%nuse(np_orb)
     MOp1 = use_nporb ! Index 1e molecular orbital
     namp1 = get_nam(TargetStates%b(MOp1))   ! Number of atomic orbitals that make up molecular orbital 
     mp1 = get_ang_mom_proj(TargetStates%b(MOp1))
!!$ INITIAL State MO COORDINATE 1
     do n_orb = 1, n_orb_max
        use_norb =  TargetStates2el%b(nstate)%nuse(n_orb)
        MO1 = use_norb
        nam1 = get_nam(TargetStates%b(MO1))
        m1 = get_ang_mom_proj(TargetStates%b(MO1))
        dM = mp1 - m1
        ! Selections Rules 
        if (  ABS(dM) > 1 ) cycle
        overlap2 = 0d0
!!$ Looping over FINAL (p) Molecular State configorations 
        do nstp_con =  1, Npcon
           MOp1 = get_na(TargetStates2el%b(npstate),nstp_con,1) ! Index 1e molecular orbital
           if ( MOp1 /= use_nporb ) cycle   ! 
           CI_stp = get_CI(TargetStates2el%b(npstate),nstp_con) ! 2e config CI
           MOp2 = get_na(TargetStates2el%b(npstate),nstp_con,2) ! Index to 1e molecular orbital
           mp2 = get_ang_mom_proj(TargetStates%b(MOp2))      ! Angular momentum projection 1e molecular orbital
           namp2 = get_nam(TargetStates%b(MOp2))    ! Number of atomic orbitals that make up molecular orbital 

!!$ Looping over INITIAL Molecular State configurations 
           do nst_con =  1, Ncon
              MO1 = get_na(TargetStates2el%b(nstate),nst_con,1)
              if ( MO1 /= use_norb ) cycle
              CI_st = get_CI(TargetStates2el%b(nstate),nst_con)
              MO2 = get_na(TargetStates2el%b(nstate),nst_con,2)
              m2 = get_ang_mom_proj(TargetStates%b(MO2))
              nam2 = get_nam(TargetStates%b(MO2))
              if ( mp2 /= m2 ) cycle

              ! overlap 2 matrix element <varphi'_2|varphi_2> 
              if (spheroidal) then   ! Spherical needs Mark to set up ovlpst.
                 overlap2 = overlap2 + CI_stp*CI_st*ovlpst(MOp2,MO2)
                 cycle
              end if

              do ip2 = 1, namp2
!!$ FINAL STATE COORDINATE 2 
!!$ Quantum numbers and functions for BETA
                 indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
                 tnp2 => bst_nr%b(indp2)
                 lap2 = get_ang_mom(tnp2)
                 CI_2p = get_CI(TargetStates%b(MOp2),ip2)

                 do i2 = 1, nam2
                    ind2 = get_na(TargetStates%b(MO2),i2)
                    tn2 => bst_nr%b(ind2)
                    la2 = get_ang_mom(tn2)
                    CI_2 = get_CI(TargetStates%b(MO2),i2)
                    if ( lap2 /= la2 ) cycle
                    temp_CI = CI_stp * CI_2p * CI_st * CI_2
                    if ( temp_CI == 0d0 ) cycle
                    overlap2 = overlap2 + temp_CI * bst_nr%ortint(indp2,ind2)
                 end do ! i2        
              end do ! ip2  

           end do ! ne_con
        end do ! nep_con
        if ( overlap2 == 0d0 ) cycle

!!$ COORDINATE 1 INTEGRALS 
        MOp1 = use_nporb ! Index 1e molecular orbital
        MO1 = use_norb
        do ip1 = 1, namp1
!!$ Quantum numbers and functions for ALPHA
           indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
           tnp1 => bst_nr%b(indp1)   !           
           fp1 => fpointer(tnp1)     ! One electron functions
           lap1 = get_ang_mom(tnp1)  ! Gets Angular momentum A.O.
           minfp = get_minf(tnp1)
           maxfp = get_maxf(tnp1)
           CI_1p = get_CI(TargetStates%b(MOp1),ip1)

           do i1 = 1, nam1
              ! Quantum numbers and functions for GAMMA
              ind1 = get_na(TargetStates%b(MO1),i1) ! Index to underlying atomic orbitals  
              tn1 => bst_nr%b(ind1)
              f1 => fpointer(tn1)
              la1 = get_ang_mom(tn1)
              minf = get_minf(tn1)
              maxf = get_maxf(tn1)
              CI_1 = get_CI(TargetStates%b(MO1),i1)
              dL = lap1 - la1  
              dM = mp1 - m1

              if (spheroidal) then   ! Use one-electron osc.str. subroutine.
                 !call oscstr_sp_oid_Liam(tnp1,tn1, lVecTmp,vVecTmp)
                 call oscstr_sp_oid(tnp1,tn1, lVecTmp,vVecTmp)
                 Coeff = overlap2 * CI_1p * CI_1
                 lVec(:) = lVec + Coeff*lVecTmp(:)
                 vVec(:) = vVec + Coeff*vVecTmp(:)
                 cycle
              end if

              if ( ABS(dL) /= 1 ) cycle ! Selection rules
              temp_CI = CI_1p * CI_1
              if ( temp_CI == 0d0 ) cycle
              ir1 = max(minf,minfp)
              ir2 = min(maxf,maxfp)

              ! Radial Integration <fp1|r|f1>
              r_overlap1 = SUM( fp1(ir1:ir2) *  f1(ir1:ir2)  * gridr(ir1:ir2) * weight(ir1:ir2) )
              !Velocity Integrals
              ! Derivative Integral Note: d/dr(phi(r)/r) = phi'(r)/r -  phi(r)/(r^2)
              ideriv1 =  SUM((( f1(ir1 + 1 : ir2 ) - f1(ir1: ir2 - 1)) /(gridr(ir1 + 1:ir2) - gridr(ir1: ir2 - 1))) * fp1(ir1: ir2 - 1) * weight(ir1:ir2 - 1))
              vint1 = SUM( fp1(ir1:ir2) * f1(ir1:ir2) * weight(ir1:ir2) / gridr(ir1:ir2) )

              if ( lap1 == la1 + 1 ) then 
                 vint1 = ideriv1 - (dble(la1)+1d0) * vint1
              else if ( lap1 == la1 - 1 ) then
                 vint1 = ideriv1 + dble(la1) * vint1
              end if

              temp = (2d0 * dble(la1) + 2d0 + dble(dL)) * (2d0 * dble(la1) + dble(dL))
              Coeff = 0d0 ! Length Guage Block Varsholovich pg 145 
              if ( dM == 0) then ! Parallel Transitions <z> 
                 Coeff = sqrt( ( dble(la1) + dble(0.5) + dble(dL) * dble(0.5) ) ** 2d0 - dble(m1 * m1))
              end if
              if ( ABS(dM) == 1) then
                 if ( lap1 == la1 + 1 ) then ! Perpendicular Transitions <x +/- iy>
                    Coeff = sqrt(dble(la1+dM*m1+2)*dble(la1+dM*m1+1))
                 else if ( lap1 == la1 - 1) then
                    Coeff = -sqrt(dble(la1-dM*m1)*dble(la1-dM*m1-1))
                 end if
                 Coeff = Coeff /  sqrt(2d0)
              end if

              Coeff = Coeff / sqrt(temp)
              
!!$ Multiply by 2 for 2e target. Operator L: z_1 + z_2 = 2*z_1
              result_l = temp_l + 2d0 * temp_CI * Coeff * r_overlap1 * overlap2
              result_v = temp_v + 2d0 * temp_CI * Coeff * vint1 * overlap2
              temp_l = result_l
              temp_v = result_v

           end do ! i1
        end do ! ip1

     end do ! n_orb use_n
  end do ! np_orb use_np

  if (spheroidal) then
     result_l = 2d0 * sum(lVec(:) * weight(:))
     result_v = 2d0 * sum(vVec(:) * weight(:))
  end if

  if(present(rl) .and. present(rv)) then
    rl = result_l
    rv = - result_v / TranEnergy
  endif

  dM = map - ma
! TranEnergy is in a.u. converted to Rydbergs here.
  if ( ABS(dM) == 1) then
     result_l = 4d0 * ABS(TranEnergy) * result_l * result_l / 3d0 ! Added orbital degeneracy g = 2 
     result_v = 4d0 * result_v * result_v / (3d0 * ABS(TranEnergy))
  else
     result_l = 2d0 * ABS(TranEnergy) * result_l * result_l / 3d0
     result_v = 2d0 * result_v * result_v / ( 3d0  * ABS(TranEnergy) )
  end if

  return
end subroutine oscstrength_2e_MOrep

subroutine  oscstrength_2e_config(nstate,npstate,TranEnergy,result_l,result_v,rl,rv)
! Calculates oscillator strengths in the form of Length <np1,np2|L|n1,n2>
! f_L =  sum C_(np1,np2) C(n1,n2) <np1,np2|L|n1,n2> 

  use grid_radial 
  use sturmian_class
  use target_states
  use one_electron_func
  use state_class

  implicit none  
  integer, intent(in):: nstate, npstate               ! Molecular state index
  real*8, intent(in):: TranEnergy
  real*8, intent(out):: result_l, result_v
  real*8, optional, intent(out) :: rl, rv
  ! Quantum numbers 
  integer::  la1, la2, lap1, lap2, ma, map, mp1, mp2, m1, m2,  dL, dM, Spin, Spinp
  integer::  parity, parityp
  ! Configuration Loops
  integer::  ne_con,  nep_con, Ncon, Npcon,  n_orb, np_orb, n_orb_max, np_orb_max 
  integer:: nuse_norb, nuse_nporb
  type(spheroidal_fn) :: oidp,oid
  type(sturmian_nr), pointer:: tn1, tn2, tnp1, tnp2                      ! One-electron orbitals
  integer:: ind1, ind2, indp1, indp2, indp1_2, ind1_2
  integer::  minf, maxf, minfp, maxfp,  ir1, ir2                         ! Radial integrations
  real*8:: Coeff, temp, CI, CIp, temp_l
  real*8:: temp_v, ideriv1,  vint1, r_overlap1, temp_overlap2, overlap2  ! Length velocity Calcs
  real*8, pointer, dimension(:)::  f1, f2, fp1, fp2                      ! f One-electron functions
  real*8, pointer, dimension(:):: weight, gridr

  weight => grid%weight
  gridr => grid%gridr
  result_l = 0d0
  temp_l = 0d0
  result_v = 0d0
  temp_v = 0d0

  ! Final State number npstate
  Npcon =  TargetStates2el%b(npstate)%nam          ! Number of A.O. and Molecular ion configurations.        
  map = NINT(TargetStates2el%b(npstate)%M )        ! Molecular state Angular Projection
  Spinp = NINT(TargetStates2el%b(npstate)%spin)    ! 2e Molecular State Spin
  parityp = TargetStates2el%b(npstate)%parity      ! 2e Molecular State Parity
  np_orb_max = TargetStates2el%b(npstate)%nusemax  ! Number of orbitals used to describe this state

  ! Initial State number nstate
  Ncon = TargetStates2el%b(nstate)%nam         
  ma = NINT(TargetStates2el%b(nstate)%M)       
  Spin = NINT(TargetStates2el%b(nstate)%spin)
  parity = TargetStates2el%b(nstate)%parity
  n_orb_max = TargetStates2el%b(nstate)%nusemax     
  
  dM = map - ma 

  ! MAKE SURE BOTH MOLECULAR STATES HAVE THE SAME SPIN and difference Angular Projections <= 1
  ! Parity must change
  if ( Spin /= Spinp ) return
  if (ABS(dM) > 1 ) return
  if ( parity == parityp) return 


! Below sums all the overlaps for COORDINATE 2  same configuratins(1s,2s,..) in COORDINATE 1
  ! FINAL State COORDINATE 1
  do np_orb = 1, np_orb_max
   
     nuse_nporb= TargetStates2el%b(npstate)%nuse(np_orb)

     ! INITIAL State COORDINATE 1
     do n_orb = 1, n_orb_max
      
        nuse_norb =  TargetStates2el%b(nstate)%nuse(n_orb)   

        overlap2 = 0d0

        ! Looping over FINAL Molecular State orbitals. COORDINATE 2
        do nep_con =  1, Npcon          
           
           indp1_2 = TargetStates2el%b(npstate)%na(nep_con)    ! Final state number np. nep A.O.

           if ( nuse_nporb /= indp1_2 ) cycle
         
           ! Quantum numbers for ALPHA
           tnp1 => bst_nr%b(indp1_2)                           !        
           lap1 = get_ang_mom(tnp1)         ! Gets Angular momentum A.O.
           mp1 = get_ang_mom_proj(tnp1)     ! Get angular projection of A.O.   
           
           ! Quantum numbers and functions for BETA
           indp2 = TargetStates2el%b(npstate)%nb(nep_con)                
           tnp2 => bst_nr%b(indp2)                                         
           fp2 => fpointer(tnp2)                                
           lap2 = get_ang_mom(tnp2)                             
           mp2 = get_ang_mom_proj(tnp2)          
           CIp = get_CI(TargetStates2el%b(npstate),nep_con)  
           
           ! Looping over INITIAL Molecular State orbitals.  COORDINATE 2
           do ne_con = 1, Ncon       
              
              ind1_2 = TargetStates2el%b(nstate)%na(ne_con)  

              if ( nuse_norb /= ind1_2 ) cycle

              ! Quantum numbers for GAMMA
              tn1 => bst_nr%b(ind1_2)
              la1 = get_ang_mom(tn1)                            
              m1 = get_ang_mom_proj(tn1)   

              ! Quantum numbers and functions for DELTA
              ind2 = TargetStates2el%b(nstate)%nb(ne_con) 
              tn2 => bst_nr%b(ind2)                                         
              f2 => fpointer(tn2)                               
              la2 = get_ang_mom(tn2)                            
              m2 = get_ang_mom_proj(tn2)            
              CI = get_CI(TargetStates2el%b(nstate),ne_con)  

              if (basis_type==2 .or. basis_type==3) then   ! Spheroidal.
                 call convert_from_sturm_to_oid(tnp2,oidp)
                 call convert_from_sturm_to_oid(tn2,oid)
                 overlap2 = overlap2 + CIp*CI * oid_overlap(oidp,oid)
                 cycle   ! Skip the spherical part.
              endif

              dL = lap1 - la1
              dM = mp1 - m1
              
              ! Selections Rules 
              if ( ABS(dL) > 1 .OR. ABS(dM) > 1 ) cycle    
              if ( lap2 /= la2 .OR. mp2 /= m2  ) cycle
                               
              temp_overlap2 = 0d0
              ! DO OVERLAP OF COORDINATE SPACE 2  
              minf = get_minf(tn2)
              maxf = get_maxf(tn2)
              minfp = get_minf(tnp2)
              maxfp = get_maxf(tnp2)      
              ir1 = max(minf,minfp)
              ir2 = min(maxf,maxfp)  
              temp_overlap2 = SUM(fp2(ir1:ir2) * f2(ir1:ir2) * weight(ir1:ir2))
              
              overlap2 = overlap2 +  CIp * CI * temp_overlap2
              
           end do  ! INITIAL STATE COORDINATE 2
           
        end do    ! FINAL STATE COORDINATE 2

        if (overlap2 == 0d0) cycle


        ! COORDINATE 1 RADIAL INTEGRALS
        ! Quantum numbers and functions for ALPHA
        indp1 = nuse_nporb                     ! Final state number np. nep A.O.
        tnp1 => bst_nr%b(indp1)                !           
        fp1 => fpointer(tnp1)                  ! One electron functions
        lap1 = get_ang_mom(tnp1)               ! Gets Angular momentum A.O.
        mp1 = get_ang_mom_proj(tnp1)           ! Get angular projection of A.O. 
        

        ! Quantum numbers and functions for GAMMA
        ind1 = nuse_norb               
        tn1 => bst_nr%b(ind1)                                          
        f1 => fpointer(tn1)                               
        la1 = get_ang_mom(tn1)                            
        m1 = get_ang_mom_proj(tn1)

        minf = get_minf(tn1)
        maxf = get_maxf(tn1)
        minfp = get_minf(tnp1)
        maxfp = get_maxf(tnp1)
        ir1 = max(minf,minfp)
        ir2 = min(maxf,maxfp) 

        if (basis_type==2 .or. basis_type==3) then   ! Spheroidal.
           call oscstr_oid(fp1,f1, ir1,ir2, dble(lap1),dble(la1), dble(mp1),dble(m1), temp_l,temp_v)
           result_l = result_l + 2d0 * overlap2 * temp_l
           result_v = result_v + 2d0 * overlap2 * temp_v
           cycle   ! Skip the spherical part.
        endif
        
        dL = lap1 - la1
        dM = mp1 - m1
        ! Selections Rules 
        if ( ABS(dL) > 1 .OR. ABS(dM) > 1 ) cycle    

        ! Radial Integration <fp1|r|f1>
        r_overlap1 = SUM( fp1(ir1:ir2) *  f1(ir1:ir2)  * gridr(ir1:ir2) * weight(ir1:ir2) ) 
        
        !Velocity Integrals
        ! Derivatiie Integral Note: d/dr(phi(r)/r) = phi'(r)/r -  phi(r)/(r^2)
        ideriv1 =  SUM((( f1(ir1 + 1 : ir2 ) - f1(ir1: ir2 - 1)) /(gridr(ir1 + 1 :ir2) - gridr(ir1: ir2 - 1))) * fp1(ir1: ir2 - 1) * weight(ir1:ir2 - 1))
        vint1 = SUM( fp1(ir1:ir2) * f1(ir1:ir2) * weight(ir1:ir2) / gridr(ir1:ir2) )
        
        if ( la1 < lap1 ) then
           vint1 = ideriv1 - (dble(la1) + 1.0) * vint1
        else if ( la1 > lap1) then
           vint1 = ideriv1 + dble(la1) * vint1
        end if
        
        temp = (2.0 * la1 + 2.0 + dL) * (2.0 * la1 + dL)
        Coeff = 0d0 ! Length Guage Block Varsholovich pg 145 
        if ( dM == 0) then ! Parallel Transitions <z> 
           Coeff = sqrt( ( dble(la1) + 0.5 + dble(dL) * 0.5 ) ** 2.0 - dble(m1 * m1)) 
        end if
        if ( ABS(dM) == 1) then
           if ( la1 < lap1   ) then  ! Perpendicular Transitions <x +/- iy> 
              Coeff = dble(dM)*sqrt(dble(la1+dM*m1+2)*dble(la1+dM*m1+1))
           else if ( la1 > lap1) then
              Coeff = -dble(dM)*sqrt(dble(la1-dM*m1)*dble(la1-dM*m1-1))
           end if
           Coeff = Coeff /  sqrt(2.0)                 
        end if
        
        Coeff = Coeff / sqrt(temp) 
        
        ! Multiply by 2 for 2e target. Operator L: z_1 + z_2 = 2*z_1
        result_l = temp_l + 2.0 * Coeff * r_overlap1 * overlap2
        result_v = temp_v + 2.0 * Coeff * vint1 * overlap2
        temp_l = result_l
        temp_v = result_v

     end do ! INITIAL STATE COORDINATE 1
     
  end do   ! FINAL STATE COORDINATE 1
  
  if(present(rl) .and. present(rv)) then
    rl = result_l
    rv = result_v
  endif

  dM = map - ma 
  ! TranEnergy is in a.u. converted to Rydbergs here.
  if ( ABS(dM) == 1) then
     result_l = 4.0 * ABS(TranEnergy) * result_l * result_l / 3.0 ! Added orbital degeneracy g = 2, should test for ma=+/-2 etc 
     result_v = 4.0 * result_v * result_v / (3.0 * ABS(TranEnergy))
  else 
     result_l = 2.0 * ABS(TranEnergy) * result_l * result_l / 3.0 
     result_v = 2.0 * result_v * result_v / ( 3.0  * ABS(TranEnergy) )
  end if


 
end subroutine oscstrength_2e_config
