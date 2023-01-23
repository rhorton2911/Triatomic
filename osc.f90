!!$subroutine  set_up_br_ratio  
!!$  use input_data   
!!$  use state_class
!!$  implicit none  
!$
!!$  integer:: Nmax_bound, n, np
!!$  real*8:: sum1, coef, tmp, rj, rjp
!!$  
!!$!  Nmax_bound = TargetStates%Nmax_bound
!!$!  print*, 'Nmax_bound=',Nmax_bound
!!$!
!!$!  do n=1,Nmax_bound
!!$!     coef = 1.0
!!$!     rj = TargetStates%b(n)%J
!!$!     if(data_in%iosc .eq. 1) then ! convert osc. str. to emmision 
!!$!        do np=1,n
!!$!           rjp = TargetStates%b(np)%J
!!$!           coef = (2*rj + 1) / (2*rjp + 1)
!!$!           TargetStates%br_ratio(n,np) = TargetStates%br_ratio(n,np) /coef
!!$!           if(TargetStates%br_ratio(n,np) .ne. 0.0) then
!!$!         !     print*, n, np,TargetStates%br_ratio(n,np) 
!!$!           endif
!!$!        enddo
!!$!     endif
!!$!     
!!$!     sum1 = SUM(TargetStates%br_ratio(n,1:Nmax_bound))
!!$!     if(sum1 .ne. 0.0) then
!!$!        TargetStates%br_ratio(n,1:Nmax_bound) = TargetStates%br_ratio(n,1:Nmax_bound) / sum1
!!$!     endif
!!$!  enddo
!!$!
!!$!
!!$!  print*, '*** Branching ratios ***'
!!$!  do n=1,Nmax_bound
!!$!     rj = TargetStates%b(n)%J
!!$!     write(*,'(I5.3,F5.1,I5)') n,  rj , TargetStates%b(n)%parity
!!$!     do np=1,n
!!$!        if(TargetStates%br_ratio(n,np) .ne. 0.0) then
!!$!           rjp = TargetStates%b(np)%J
!!$!           write(*,'( "      " ,I5.3,F5.1,I5,F15.6)')  np, TargetStates%b(np)%J, TargetStates%b(np)%parity, TargetStates%br_ratio(n,np) 
!!$!        endif
!!$!     enddo
!!$!  enddo
!!$
!!$
!!$end subroutine set_up_br_ratio
!$$-----------------------------------------------------------------------------
subroutine oscstr1()
  use input_data 
  use target_data
  use one_electron_func
  use sturmian_class
  use target_states

  implicit none  
  type(state), pointer :: stateObj, stateObjP
  integer:: n, np, Nst, ido_osc ! Loops and input option controls
  ! Results
  real*8::  result_l, result_v,  result_p, result_pc
  real*8:: result_p_par, result_p_perp
  real*8:: TranEnergy, ioniz_en, temp_p, Energy, Osc 
  integer:: map, ma, dM, npmax, Ref_E, Ref_O, calc_type
  real*8 :: rl, rv

  INTERFACE
    subroutine oscstr_st_oid(stp,st, result_l,result_v, rl, rv)
      use grid_radial
      use one_electron_func
      use target_states
      implicit none
      type(state), intent(in) :: stp,st
      real*8, intent(out) :: result_l,result_v 
      real*8, intent(out), optional :: rl,rv 
    end subroutine
    subroutine oscstr1_st(iosc,n,np,pn,pnp,TranEnergy,result_l,result_v,rl,rv)
      use grid_radial 
      use sturmian_class
      use target_states
      use one_electron_func
      implicit none  
      integer, intent(in):: iosc
      integer, intent(in):: n, np     ! Molecular states
      type(state), intent(in):: pn, pnp   ! One electron states
      real*8, intent(in):: TranEnergy
      real*8, intent(out):: result_l, result_v
      real*8, intent(out), optional :: rl, rv
    end subroutine
  END INTERFACE


  calc_type = data_in%calculation_type
  Nst = basis_size_st(TargetStates)

  npmax = min(data_in%iosc,Nst)
  
  open(301,file='osc.strength')
  open(302,file='dip.mom')
  
  write(301,'("ONE-ELECTRON OSCILLATOR STRENGTHS  " )')
  write(301,'("Exact solutions for the Oscillator Strengths are for the R = 2.0 case" )')
  write(301,'("Ref[0]: None, Ref[1]: D.R. Bates 1954, Ref[2]: D.R. Bates 1953, Ref[3]: D.R. Bates 1951 " )')
  write(301,'("Transition (np - n)   Length      Velocity       Exact       Ref" )')
     
  write(302,'("STATIC DIPOLE POLARIZABILITY") ')
  write(302,'("Exact solution for the Polarizability are for the R =2.0 case." )')
  write(302,'("Ref[1]: D. M. Bishop 1978 , Ref[2]: , Ref[3]: " )')
  write(302,'("Polarizability calculated using oscillator strength length guages." )')

  do np = 1, npmax 
     stateObjP => TargetStates%b(np)
     map = get_ang_mom_proj(stateObjP)

     result_p = 0d0
     result_pc = 0d0
     result_p_par = 0d0
     result_p_perp = 0d0

     do n = np + 1, Nst
        stateObj => TargetStates%b(n)
        ! Loops (initial condition) will need changing if want to work out static-dipole polarizabiliy for excited states (n = 1, Nst and may need to change polarizability conditions)
        
        ma = get_ang_mom_proj(stateObj) ! Get magnetic sublevels? M = m for one-electron
        dM = map - ma 
        
        if ( get_par(stateObjP) == get_par( stateObj) .and. data_in%good_parity) cycle  ! Selection rules
        if ( ABS(dM) >= 2 )  cycle  ! Selections rules dm = +/- 1, 0

        ! Write to oscillator strength file osc.strength and polarizability to dip.mom
        if ( n == np +1) then
           write(302,'("Transition (np - n)   Excitation Energy(eV)    Polarizability(a_0^3)" )')          
        end if

        TranEnergy =  get_energy_st(stateObj) - get_energy_st(stateObjP)  ! In units of a.u., hence multiply in osc by 2 for Ryd.
        
        if (calc_type.eq.0 .or. calc_type.eq.1) then 
           call oscstr1_st(data_in%iosc,n,np,stateObj,stateObjP,TranEnergy,result_l,result_v,rl,rv)                    
        else if (calc_type==2 .or. calc_type==3) then
           call oscstr_st_oid(stateObj,stateObjP, result_l,result_v )
        end if
           
        ! Polarizabilty calulated from Oscillator strengths
        temp_p = 0d0
        ! a_T = 1/3 * a_par + 2/3 * a_perp
        if ( ABS(dM) == 1) then
           temp_p =  result_l / ( 2.0 * TranEnergy * TranEnergy)  ! Summing over f_n dm = +/-1, therefore divide by 2
           result_p_perp =  result_p_perp + 1.5 * temp_p
        else if ( dM == 0 ) then
           temp_p =  result_l / ( TranEnergy * TranEnergy )
           result_p_par =  result_p_par + 3.0 * temp_p
        end if
!        if (  get_energy_st(stateObj) > 0 ) then ! For continuum contribution
        if (  get_energy_st(stateObj) - TargetStates%en_ion > 0 ) then ! For continuum contribution
           result_pc = result_pc + temp_p
        end if
        result_p = result_p + temp_p               
        
        ! Write to oscillator strength file osc.strength and polarizability to dip.mom
        call get_target_oscstr( data_target, get_n_majconf(stateObjP), get_l_majconf(stateObjP), map, get_n_majconf(stateObj), get_l_majconf(stateObj), ma, Osc, Ref_O ) 
        write(301,'(2X,A6," -",A5,2X,3(2X,1P,Es11.4),2X,I5)')  stateObjP%label, stateObj%label, result_l, result_v, Osc, Ref_O
        ! file: dip.mom
        write(302,'(2X,A6," -",A5,11X,F8.3,15X,1P,E11.4)')   stateObjP%label, stateObj%label, data_in%eV * TranEnergy, temp_p         

        
     end do ! Loop over n-stantes 
     

     if ( np == 1) then
        write(302,'("State:",A6,", Full Spectrum Polarizability:",1P,E11.4,", Ref[1]: 2.8643 " )')  stateObjP%label, result_p
        write(302,'("State:",A6," Continuum Polarizability:",1P,E11.4 )')  stateObjP%label, result_pc
        write(302,'("State:",A6,", Parallel:",1P,E11.4,", Ref[1]: 5.0776489")')  stateObjP%label, result_p_par
        write(302,'("State:",A6,", Perpendicular:",1P,E11.4,", Ref[1]: 1.757648" )')  stateObjP%label, result_p_perp
     else
        write(302,'("State:",A6,", Full Spectrum Polarizability:",1P,E11.4,", Continuum Polarizability:",E11.4 )')  stateObjP%label, result_p, result_pc
        write(302,'("State:",A6,", Parallel:",1P,E11.4,", Perpendicular:",E11.4 )')  stateObjP%label, result_p_par, result_p_perp
     end if

     write(302,*)
  end do ! Loop over np-stantes
  
  close(301)
  close(302)
  
  
end subroutine oscstr1
!
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
!
subroutine oscstr_st_oid(stp,st, result_l,result_v, rl, rv)
  !
  use grid_radial
  use one_electron_func
  use target_states
  implicit none

  type(state), intent(in) :: stp,st
  real*8, intent(out) :: result_l,result_v   ! Length and velocity osc.str.
  real*8, intent(out), optional :: rl,rv   ! Length and velocity dipole matrix elements

  type(sturmian_nr), pointer :: sturmp,sturm
  integer :: indp,ind, dm, i1,i2
  real*8 :: lp,l, mp,m, CIp,CIpxCI, tmp_l,tmp_v, dE
  real*8, dimension(:), pointer :: weightVec, ptrp,ptr
  real*8, dimension(:), allocatable :: vecp,vec
  real*8, dimension(1:grid%nr) :: lVecTmp,vVecTmp, lVec,vVec


  result_l = 0d0; result_v = 0d0
  dm = get_ang_mom_proj(stp) - get_ang_mom_proj(st)
  if (abs(dm) > 1) return   ! We don't want your kind around here.
  if (get_par(stp) == get_par(st)) return   ! Likewise.

  lVec(:) = 0d0; vVec(:) = 0d0
  do indp = 1, get_nam(stp)
     sturmp => bst_nr%b( get_na(stp,indp) )
     CIp = get_CI(stp,indp)

     do ind = 1, get_nam(st)
        sturm => bst_nr%b( get_na(st,ind) )
        CIpxCI = CIp * get_CI(st,ind)

        call oscstr_sp_oid(sturmp,sturm, lVecTmp,vVecTmp)
        lVec(:) = lVec(:) + CIpxCI*lVecTmp(:)
        vVec(:) = vVec(:) + CIpxCI*vVecTmp(:)
     enddo
  enddo

  weightVec => grid%weight
  result_l = sum( lVec(:) * weightVec(:) )
  result_v = sum( vVec(:) * weightVec(:) )

  if(present(rl).and.present(rv)) then
    rl = result_l
    rv = result_v
  endif

  dE = abs( get_energy(stp) - get_energy(st) )
  result_l = 2d0/3d0 * dE * result_l**2
  result_v = 2d0/3d0 / dE * result_v**2
           
  if (abs(dm) == 1) then   ! Degeneracy.
     result_l = 2d0*result_l
     result_v = 2d0*result_v
  endif
  

end subroutine oscstr_st_oid
!
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
!
subroutine oscstr_sp_oid(sturmp,sturm, lVec,vVec)
  !
  use grid_radial
  use input_data
  use sturmian_class
  implicit none

  type(sturmian_nr), intent(in) :: sturmp,sturm
  real*8, dimension(1:grid%nr), intent(out) :: lVec,vVec

  integer :: i1,i2
  real*8 :: R, lp,mp, l,m, dm,C1,C3,y, Yint
  real*8, dimension(:), pointer :: rhoVec, fpVec,fVec,gVec
  real*8, dimension(:), allocatable :: ppRVec, pR2Vec

  lVec(:) = 0d0; vVec(:) = 0d0
  R = data_in%Rd
  rhoVec => grid%gridr

  lp = dble(get_ang_mom(sturmp))     ; l = dble(get_ang_mom(sturm))
  mp = dble(get_ang_mom_proj(sturmp)); m = dble(get_ang_mom_proj(sturm))
  dm = mp - m
  C1 = Yint(lp,mp, 1d0,dm, l,m); C3 = Yint(lp,mp, 3d0,dm, l,m)
  if (C1==0d0 .and. C3==0d0) return

  if (nint(lp-l) == -1) then
     y = l+1d0
  elseif (nint(lp-l) == 1) then
     y = -l
  else
     y = 0d0
  end if

  fpVec => fpointer(sturmp)
  fVec => fpointer(sturm); gVec => gpointer(sturm)
  i1 = max( get_minf(sturmp), get_minf(sturm) )
  i2 = min( get_maxf(sturmp), get_maxf(sturm) )

  allocate( ppRVec(i1:i2), pR2Vec(i1:i2) )
  ppRVec(:) = rhoVec(i1:i2) * (rhoVec(i1:i2)+R)
  pR2Vec(:) = rhoVec(i1:i2) + R/2d0

  if (nint(dm) == 0) then   ! Parallel transition.
     if (C1 /= 0d0) lVec(i1:i2) = C1 * (ppRVec(:) + R*R/10d0)
     if (C3 /= 0d0) lVec(i1:i2) = lVec(i1:i2) - C3*R*R/10d0
     !if (C3 /= 0d0) lVec(i1:i2) = lVec(i1:i2) - C3*R*R * sqrt(21.0d0)/70.0d0   ???
     lVec(i1:i2) = pR2Vec(:) * lVec(i1:i2)
  
     if (C1 /= 0d0) vVec(i1:i2) = C1 * (ppRVec(:)*gVec(i1:i2) + y*pR2Vec(:)*fVec(i1:i2))
     
  elseif ( nint(abs(dm)) == 1 ) then   ! Perpendicular transition.
     if (C1 /= 0d0) lVec(i1:i2) = C1 * (ppRVec(:) + R*R/5d0)
     if (C3 /= 0d0) lVec(i1:i2) = lVec(i1:i2) - C3*R*R/dsqrt(150d0)
     !if (C3 /= 0d0) lVec(i1:i2) = lVec(i1:i2) - C3*R*R * sqrt(14.0d0)/70.0d0  ???
     lVec(i1:i2) = -dm * dsqrt(ppRVec(:)) * lVec(i1:i2)
  
     if (C1 /= 0d0) vVec(i1:i2) = -C1*dsqrt(ppRVec(:)) * (pR2Vec(:)*gVec(i1:i2) + (abs(m)*R*R/4d0/ppRVec(:) + y)*fVec(i1:i2))
  
  end if

  deallocate(ppRVec, pR2Vec)
  lVec(i1:i2) = fpVec(i1:i2) * lVec(i1:i2) * fVec(i1:i2)
  vVec(i1:i2) = fpVec(i1:i2) * vVec(i1:i2)


end subroutine oscstr_sp_oid



subroutine oscstr_sp_oid_Liam(sturmp,sturm, lVec,vVec)
  !Liam modified from above routine to use angular_quadratic subroutine rather than explicitly coding special cases of it
  use grid_radial
  use input_data
  use sturmian_class
  implicit none

  type(sturmian_nr), intent(in) :: sturmp,sturm
  real*8, dimension(1:grid%nr), intent(out) :: lVec,vVec

  integer :: i1,i2
  real*8 :: R, lp,mp, l,m, dm,y, C1, C3, Yint
  real*8, dimension(:), pointer :: rhoVec, fpVec,fVec,gVec
  real*8, dimension(:), allocatable :: ppRVec, pR2Vec, Jvec

  lVec(:) = 0d0; vVec(:) = 0d0
  R = data_in%Rd
  rhoVec => grid%gridr

  lp = dble(get_ang_mom(sturmp))     ; l = dble(get_ang_mom(sturm))
  mp = dble(get_ang_mom_proj(sturmp)); m = dble(get_ang_mom_proj(sturm))
  dm = mp - m
  C1 = Yint(lp,mp, 1d0,dm, l,m); C3 = Yint(lp,mp, 3d0,dm, l,m)

  if (nint(lp-l) == -1) then
     y = l+1d0
  elseif (nint(lp-l) == 1) then
     y = -l
  else
     y = 0d0
  end if

  fpVec => fpointer(sturmp)
  fVec => fpointer(sturm); gVec => gpointer(sturm)
  i1 = max( get_minf(sturmp), get_minf(sturm) )
  i2 = min( get_maxf(sturmp), get_maxf(sturm) )

  allocate( ppRVec(i1:i2), pR2Vec(i1:i2), Jvec(i1:i2))
  ppRVec(:) = rhoVec(i1:i2) * (rhoVec(i1:i2)+R)
  pR2Vec(:) = rhoVec(i1:i2) + R/2d0

  if((-1)**lp /= (-1)**l) then !Electric dipole transition
    !if (C1==0d0 .and. C3==0d0) return
       
    call angular_quadratic(nint(lp),nint(mp),1,nint(dm),nint(l),nint(m),Jvec,i1,i2)
    if(all(Jvec == 0.0d0)) return

    if (nint(dm) == 0) then   ! Parallel transition.

       !Length gauge updated to use angular_quadratic routine
       lVec(i1:i2) = pR2Vec(i1:i2) * Jvec(i1:i2)

       !Have not updated velocity gauge
       if (C1 /= 0d0) vVec(i1:i2) = C1 * (ppRVec(:)*gVec(i1:i2) + y*pR2Vec(:)*fVec(i1:i2))
       
    elseif ( nint(abs(dm)) == 1 ) then   ! Perpendicular transition.
       
      !Length gauge updated to use angular_quadratic routine
      lVec(i1:i2) = -dm * dsqrt(ppRVec(:)) * JVec(i1:i2)
  
      !Have not updated velocity gauge
      if (C1 /= 0d0) vVec(i1:i2) = -C1*dsqrt(ppRVec(:)) * (pR2Vec(:)*gVec(i1:i2) + (abs(m)*R*R/4d0/ppRVec(:) + y)*fVec(i1:i2))
  
    end if

  else !if(.not.(lp==0.and.l==0)) then !Electric quadrupole transition
    !Liam added this for quadrupole transitions

    call angular_quadratic(nint(lp),nint(mp),2,nint(dm),nint(l),nint(m),Jvec,i1,i2)
    
    if(all(Jvec == 0.0d0)) return

    if(lp==0 .and. l==0) return

    vVec = 0.0d0 !Have not coded velocity gauge for quadrupole transitions...

    if (nint(dm) == 0) then
       
      lVec(i1:i2) = - 2.0d0/3.0d0*(ppRVec(i1:i2) + 2.0d0*pR2Vec(i1:i2)**2) * Jvec(i1:i2)

    elseif (nint(abs(dm)) == 1) then

      lVec(i1:i2) = dm * 2.0d0/sqrt(6.0d0) * pR2Vec(i1:i2)*dsqrt(ppRVec(i1:i2)) * Jvec(i1:i2)

    elseif (nint(abs(dm)) == 2) then

      lVec(i1:i2) = - 4.0d0/sqrt(6.0d0) * ppRVec(i1:i2) * Jvec(i1:i2)

    endif !dm

    !Add R^2 term if state_f = state_i ???

  endif !lp /= l

  deallocate(ppRVec, pR2Vec, Jvec)
  lVec(i1:i2) = fpVec(i1:i2) * lVec(i1:i2) * fVec(i1:i2)
  vVec(i1:i2) = fpVec(i1:i2) * vVec(i1:i2)

end subroutine oscstr_sp_oid_Liam

subroutine angular_quadratic(lamg,mg,l,m,lama,ma,Jvec, i1, i2)
  !Liam wrote following Eq. (E.5) in Jeremy's thesis
  use input_data
  use grid_radial
  implicit none
  integer, intent(in) :: lamg, mg, l, m, lama, ma, i1, i2
  real*8, dimension(i1:i2), intent(out) :: Jvec
  real*8, external :: Yint
  real*8 :: j_lminus2, j_l, j_lplus2, R

  R = data_in%Rd

  j_l = Yint(dble(lamg),dble(mg),dble(l),dble(m),dble(lama),dble(ma))
  j_lplus2 = Yint(dble(lamg),dble(mg),dble(l+2),dble(m),dble(lama),dble(ma))
  if(l >= 2) j_lminus2 = Yint(dble(lamg),dble(mg),dble(l-2),dble(m),dble(lama),dble(ma))

  Jvec = j_l * (grid%gridr(i1:i2)*(grid%gridr(i1:i2)+R) + R**2/2.0d0 * dble((l*l+1)+m**2-1)/dble((2*l-1)*(2*l+3))) & 
    &  - (-R)**2/4.0d0 * j_lplus2/dble(2*l+3) * sqrt(dble((l+1-m)*(l+1+m)*(l+2-m)*(l+2+m))/dble((2*l+1)*(2*l+5)))

  if(l >= 2) Jvec = Jvec - R**2/4.0d0 * j_lminus2/dble(2*l-1) * sqrt(dble((l-1-m)*(l-1+m)*(l-m)*(l+m))/dble((2*l-3)*(2*l+1)))

end subroutine angular_quadratic


subroutine oscstr_oid(vecp,vec, i1,i2, lp,l, mp,m, tmp_l,tmp_v)
  !
  use grid_radial
  use input_data
  implicit none

  real*8, dimension(i1:i2), intent(in) :: vecp,vec
  integer, intent(in) :: i1,i2
  real*8, intent(in) :: lp,l, mp,m
  real*8, intent(out) :: tmp_l,tmp_v

  integer :: calc_type, i3,i4
  real*8 :: R, dl,dm,y, C1,C3, Yint
  real*8, dimension(i1:i2) :: angVec,tmpVec, dVec
  real*8, dimension(:), pointer :: rhoVec,weightVec


  tmp_l = 0d0; tmp_v = 0d0

  calc_type = data_in%calculation_type
  R = data_in%Rd
  rhoVec => grid%gridr; weightVec => grid%weight  
  angVec(:) = 0d0; tmpVec(:) = rhoVec(i1:i2)*(rhoVec(i1:i2)+R)

  dl = lp - l; dm = mp - m
  if (dl == -1) then
     y = l+1d0
  elseif (dl == 1) then
     y = -l
  else
     y = 0d0
  end if
  C1 = Yint(lp,mp, 1d0,dm, l,m); C3 = Yint(lp,mp, 3d0,dm, l,m)
  if (C1==0d0 .and. C3==0d0) return

  i3 = i1+1; i4 = i2-1
  dVec(:) = 0d0
  dVec(i1:i4) = (vec(i3:i2) - vec(i1:i4))/(rhoVec(i3:i2) - rhoVec(i1:i4))

  if (nint(abs(dm)) == 0) then   ! Parallel transition.
     if (C1 /= 0d0) angVec(:) = C1 *(tmpVec(:)+R*R/10d0)
     if (C3 /= 0d0) angVec(:) = angVec(:) - C3 *R*R/10d0
     angVec(:) = angVec(:) * (rhoVec(:)+R/2d0)

     dVec(:) = tmpVec(:)*dVec(:) + y*(rhoVec(i1:i2)+R/2d0)*vec(i1:i2)
     if (calc_type == 3) dVec(:) = (dVec(:) - rhoVec(i1:i2)*vec(i1:i2)) /(rhoVec(i1:i2)+R)/(rhoVec(i1:i2)+R)
     dVec(:) = C1 * dVec(:)
     
  elseif ( nint(abs(dm)) == 1 ) then   ! Perpendicular transition.
     if (C1 /= 0d0) angVec(:) = -C1 *(tmpVec(:)+R*R/5d0)
     if (C3 /= 0d0) angVec(:) = angVec(:) + C3 *R*R/5d0 /dsqrt(6d0)
     angVec(:) = angVec(:) * dsqrt(tmpVec(:))

     dVec(:) = -dm*(rhoVec(i1:i2)+R/2d0)*dVec(:) + (m*R*R/4d0/tmpVec(i1:i2) - dm*y)*vec(i1:i2)
     if (calc_type == 3) dVec(:) = (dVec(:) + dm*(rhoVec(i1:i2)+R/2d0)/(rhoVec(i1:i2)+R)*vec(i1:i2))/(rhoVec(i1:i2)+R)/(rhoVec(i1:i2)+R)
     dVec(:) = C1 * dsqrt(tmpVec(i1:i2)) * dVec(:)

  end if

  if (calc_type == 3) angVec(:) = angVec(:) /(rhoVec(i1:i2)+R)/(rhoVec(i1:i2)+R)
  tmp_l = sum( vecp(i1:i2)*vec(i1:i2) * angVec(:) * weightVec(i1:i2) )
  tmp_v = sum( vecp(i1:i2)*dVec(:) * weightVec(i1:i2) )


end subroutine oscstr_oid
!
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
!
subroutine oscstr1_st(iosc,n,np,pn,pnp,TranEnergy,result_l,result_v,rl, rv) !,rlpar,rvpar,rlperp,rvperp)
  use grid_radial 
  use sturmian_class
  use target_states
  use one_electron_func


  implicit none  
  integer, intent(in):: iosc
  integer, intent(in):: n, np     ! Molecular states
  type(state), intent(in):: pn, pnp   ! One electron states
  real*8, intent(in):: TranEnergy
  real*8, intent(out):: result_l, result_v
  real*8, intent(out), optional:: rl, rv
  !real*8, intent(out):: kkkkresult_lperp, result_vperp
  !real*8, intent(out):: result_lpar, result_vpar
  !real*8, intent(out), optional :: rlpar, rvpar, rlperp, rvperp !par, perp dipole matrix elements

  integer::  Nst, Nstp, l, lp, ma, map, ne, nep,  dL, dM  !quantum numbers and loops
  type(sturmian_nr), pointer:: nn, nnp  !  one-electron orbitals
  integer:: i, ipp
  integer::  minf, maxf, minfp, maxfp,  i1, i2 ! Radial integrations
  real*8:: Coeffperp, Coeffpar, nnpRnn, temp, CI, CIp, temp_l
  real*8:: temp_v, ideriv,  vint
  real*8, pointer, dimension(:)::  f, fp ! f One-electron functions
  real*8, pointer, dimension(:):: weight, gridr
  real*8:: g !Relative degeneracy of states

  !New variables:
  real*8, dimension(3):: result_larr
  real*8, dimension(3):: result_varr

  print*, "FIX OSCSTR_1ST"
  stop
 
 ! weight => grid%weight
 ! gridr => grid%gridr
 ! result_l = 0d0
 ! temp_l = 0d0
 ! result_v = 0d0
 ! temp_v = 0d0

 ! Nstp = get_nam(pnp) 
 ! map = TargetStates%b(np)%M

 ! Nst = get_nam(pn)  ! ncm = nam which is the size of the array na(:), which is also the number of CI coeff
 ! ma = TargetStates%b(n)%M ! Get magnetic sublevels? M = m for one-electron
 ! 
 ! g =
 ! !dM = map - ma 

 ! do nep = 1, Nstp  ! Looping over the number of one-electron orbitals which make up each molecule state
 !    
 !    CIp = get_CI(pnp, nep)
 !    ipp = get_na(pnp,nep,1)
 !    nnp => bst_nr%b(ipp)
 !    fp => fpointer(nnp)
 !    lp = get_ang_mom(nnp)
 !    
 !    do ne =  1, Nst
 !       
 !       CI = get_CI(pn, ne) ! Get CI coefficient for that one-electron function
 !       i = get_na(pn,ne,1) ! which is na(ne), where nc = 1,2,...Nst
 !       nn => bst_nr%b(i)   ! Makes nn the one-electron states?
 !       f => fpointer(nn)   ! One electron functions?
 !       l = get_ang_mom(nn) ! Gets angular momentum from the one-electron orbitals
 !       
 !       dL = lp - l
 !       
 !       dM = 
 !     
 !       if ( dL == 0 ) cycle 
 !       if ( ABS(dL) > 1 ) cycle  ! Selection rules dl = +/- 1 
 ! 
 !       temp = (2 * l + 2 + dL) * (2 * l + dL)
 !       
 !       ! Extent of radial grid for integrations
 !       minf = get_minf(nn)
 !       maxf = get_maxf(nn)
 !       minfp = get_minf(nnp)
 !       maxfp = get_maxf(nnp)
 !       
 !       i1 = max(minf,minfp)
 !       i2 = min(maxf,maxfp)
 !       
 !       ! Radial Integration <fp|r|f>
 !       nnpRnn = SUM( fp(i1:i2) *  f(i1:i2)  * gridr(i1:i2) * weight(i1:i2) ) ! One-electron function normalized so that r^2 term is taken out of integrations.
 !       
 !       
 !       !Velocity Integrals
 !       ! Derivatiie Integral Note: d/dr(phi(r)/r) = phi'(r)/r -  phi(r)/(r^2)
 !       ideriv =  SUM((( f(i1+1:i2) - f(i1:i2-1)) / (gridr(i1+1:i2) - gridr(i1:i2-1))) * fp(i1:i2-1) * weight(i1:i2-1))
 !       
 !       ! Varsholovich pg 147 and look at derivative operator second term.
 !       vint = SUM( fp(i1:i2) * f(i1:i2) * weight(i1:i2) / gridr(i1:i2) )
 !       if ( l < lp ) then
 !          vint = ideriv - (l + 1.0) * vint
 !       else if ( l > lp) then
 !          vint = ideriv + dble(l) * vint
 !       end if
 !       
 !       Coeff = 0d0 ! Length Guage Block Varsholovich pg 145 
 !       if ( dM == 0) then ! Parallel Transitions <z> 
 !          Coeffpar = sqrt( ( dble(l) + 0.5 + dble(dL) * 0.5 ) ** 2.0 - dble(ma * ma)) 
 !       end if
 !       if ( ABS(dM) == 1) then
 !          if (ma < map  .AND. l < lp   ) then  ! Perpendicular Transitions <x +/- iy> 
 !             Coeffperp = sqrt(dble((l + ma + 2 ) * (l + ma + 1)))
 !          else if (  ma > map  .AND. l < lp ) then
 !             Coeffperp = - sqrt(dble((l - ma + 2) * (l - ma + 1)))  
 !          else if (  ma < map  .AND. l > lp  ) then         
 !             Coeffperp = - sqrt(dble((l - ma ) * (l - ma - 1)))  
 !          else if (  ma > map  .AND. l > lp ) then
 !             Coeffperp = sqrt(dble((l + ma ) * (l + ma - 1)))               
 !          end if
 !          Coeffperp = Coeffperp /  sqrt(2.0)                 
 !       end if
 !       
 !       Coeffpar = Coeffpar / sqrt(temp) 
 !       Coeffperp = Coeffperp / sqrt(temp) 
 !       
 !       result_lpar = temp_lpar +  CIp * CI * Coeffpar * nnpRnn 
 !       result_lperp = temp_lperp +  CIp * CI * Coeffperp * nnpRnn 
 !       result_vpar = temp_vpar +  CIp * CI * Coeffpar * vint
 !       result_vperp = temp_vperp +  CIp * CI * Coeffperp * vint
 !       temp_lpar = result_lpar
 !       temp_lperp = result_lperp
 !       temp_vpar = result_vpar
 !       temp_vperp = result_vperp
 !          
 !    end do ! ne RHS one-electron functions
 !    
 ! end do    ! nep LHS one-electron functions
 ! 
 ! if((present(rlpar).and.present(rvpar)) .and. (present(rlperp) .and. present(rvperp))) then
 !   rlpar = result_lpar
 !   rlperp = result_lperp
 !   rvpar = result_vpar
 !   rvperp = result_vperp
 ! endif

 ! !if ( ABS(dM) == 1) then
 ! !   result_l = 4.0 * ABS(TranEnergy) * result_l * result_l / 3.0 ! Added orbital degeneracy g = 2, should test for ma=+/-2 etc 
 ! !   result_v = 4.0 * result_v * result_v / (3.0 * ABS(TranEnergy))
 ! !else 
 ! !   result_l = 2.0 * ABS(TranEnergy) * result_l * result_l / 3.0 ! Need to add orbital degeneracy, prob best to do it in the coefficients
 ! !   result_v = 2.0 * result_v * result_v / ( 3.0  * ABS(TranEnergy) )
 ! !end if

 ! result_l = g * 2.0 * ABS(TranEnergy) * (result_lpar * result_lpar + result_lperp * resultl_perp) / 3.0
 ! result_v = g * 2.0 * (result_vpar * result_vpar + result_vperp * result_vperp) / ( 3.0  * ABS(TranEnergy) )
   
end subroutine oscstr1_st


