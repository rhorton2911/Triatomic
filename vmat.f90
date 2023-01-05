subroutine vdirect(Mtot,nchf,nchi,nstf,nsti,nqmf,nqmi,npk_nchf,npk_nchi,chil,formcut,nr,weight,gridr,Ldw,iBorn,vmatt)
! j-j coupling
! This is adaptation of routine for matrix element of V(1,2) electron-electron potential for 4 functions:
! <n1 n2: Mtot | V(1,2) | n1p n2p : Mtot>, where V(1,2) = sum_{lam} V_{lam}(1,2), and Mtot is total anguar moment of the two electrons.

  use sturmian_class
  use target_states
  use one_electron_func
  use channels
  use input_data
  use distpot
  use vnc_module
  use pol_2el
  use MPI_module

  implicit none

  integer, intent(in):: Mtot    ! total ang.momentum projection of two electrons - integer!
  integer, intent(in):: nchf,nchi
  integer, intent(in):: nstf,nsti
  integer, intent(in):: nqmf, nqmi
  integer, intent(in):: npk_nchf,npk_nchi
  type(basis_sturmian_nr):: chil
  real*8,intent(in):: formcut
  integer, intent(in):: nr
  real*8, dimension(nr), intent(in):: weight, gridr
!!$  real*8, dimension(nr), intent(in):: vdist
!!$  integer, intent(in):: maxvdist
!!$  real*8, dimension(nr,0:lamtop_vc), intent(in):: vdcore
!!$  integer, intent(in):: maxvdc, lamtop_vc
  integer, intent(in):: Ldw
  integer, intent(in):: iBorn
  real*8, dimension(nqmf,nqmi), intent(inout):: vmatt

  real*8:: Yint, CGC

  type(state), pointer:: state_f,state_i
  integer:: n1,n2,n1p,n2p
  integer:: l1, l2, l1p, l2p
  integer:: maxfm, minfm, i1, i2, i2form,minfun, maxfun, minfun1, maxfun1, lam, lammin, lammax
  real*8:: rlam, reslam
  real*8:: tmp, tmp1,  sum2, sum_rad, tmp_core_exch, ang1, ang2
  integer:: minf1,maxf1, minf1p,maxf1p,minf2,maxf2, minf2p,maxf2p
  real*8, pointer, dimension(:):: f1, f1p, f2, f2p
  real*8, dimension(nr):: temp, fun, fun1, fun11, ttt

  type(sturmian_nr), pointer:: pn1,pn2,pn1p,pn2p  !  one-electron orbitals

  integer:: kf, ki, kqf, kqi
  real*8:: pi, tmpsumint
  integer:: i2all, i2all_1,  la_vdcore
  real*8:: sum1, tmp_pp, rval, Vval11, Vval12, u1s,u2s, rCIf, rCIi, Zasym, rmu
  integer:: nci, ncf, ncim, ncfm, m2, m2p, m1, m1p, Mfst, Mist, mu
  integer:: i, imax

  pi = acos(-1d0)

  Zasym = data_in%Zasym

  lammin = -1; lammax = -1
  n2p = 0; n2 = 0
  maxfm = nr

  l1 = Lp_ch(nchf)   ! projectile
  m1 = Mp_ch(nchf)
  l1p = Lp_ch(nchi)   ! projectile
  m1p = Mp_ch(nchi)

  state_f => TargetStates%b(nstf)
  state_i => TargetStates%b(nsti)

  ncfm = get_nam(state_f)
  ncim = get_nam(state_i)

  Mfst = get_ang_mom_proj(state_f)
  Mist = get_ang_mom_proj(state_i)

  if(m1+Mfst  .ne. m1p+Mist) return

!!$ n2, n2p are indexes to one-electron target orbitals
!!$ while n1, n1p are indexes to continuum waves

  rmu = Mfst - Mist
  mu = Mfst - Mist


  i2all = 0
  ttt(:) = 0d0
  do ncf=1,ncfm
     m2 = get_ma(state_f,ncf,1)
     n2 = get_na(state_f,ncf,1)
     pn2 => bst_nr%b(n2)
     l2 = get_ang_mom(pn2)

     if(m2 .ne. Mfst) then
        if (myid == 0 ) print*, '>>>>>>>>>> m2 != Mstf: ', nstf, Mfst, m2, l2, n2
     endif


     f2 => fpointer(pn2)
     minf2 = get_minf(pn2)
     maxf2 = get_maxf(pn2)

     rCIf = get_CI(state_f,ncf)

     do nci=1,ncim
        m2p = get_ma(state_i,nci,1)
        n2p = get_na(state_i,nci,1)
        pn2p => bst_nr%b(n2p)
        l2p = get_ang_mom(pn2p)

        if(m2p .ne. Mist) then
           if (myid == 0 ) print*, '>>>>>>>>>>> m2p != Msti: ',  nsti, Mist, m2p, l2p, n2p
        endif

        f2p => fpointer(pn2p)
        minf2p = get_minf(pn2p)
        maxf2p = get_maxf(pn2p)

        rCIi = get_CI(state_i,nci)

        minfun = max(minf2,minf2p)
        maxfun = min(maxf2,maxf2p)

        lammin=max(abs(l1-l1p),abs(l2-l2p))
        lammax=min(l1+l1p,l2+l2p)
        if(lammin .gt. lammax) cycle

        fun(:) = 0d0
        fun(minfun:maxfun) = rCIf*rCIi*(f2(minfun:maxfun)*f2p(minfun:maxfun)) * weight(minfun:maxfun)

        do lam=lammin,lammax
           rlam = lam
           tmp1 = 0d0
           ang1 = dble((-1)**(mu))*Yint(dble(l1),dble(m1),dble(lam),-rmu,dble(l1p),dble(m1p))
           ang2 = Yint(dble(l2),dble(m2),dble(lam),rmu,dble(l2p),dble(m2p))


           tmp1 = ang1 * ang2
           if( tmp1 .eq. 0d0) then
              cycle
           endif

           maxfm = nr
           call form(lam,fun,minfun,maxfun,maxfm,temp,i1,i2form)
           i2 = i2form

           ttt(i1:i2) = ttt(i1:i2) + temp(i1:i2)*tmp1
           i2all = max(i2all,i2)

        enddo  ! end lam loop

     enddo  ! end nci loop
  enddo  ! end ncf loop

!!$ -------------------------------------------------------------------------------
!!$  Test ttt values here
!!$   if (myid == 0 ) print*,'I am in Vee'
!!$   if (myid == 0 ) print*,'ttt(200) = ',ttt(200)
!!$ -------------------------------------------------------------------------------

!!$  Subtract nuclear potential VLambdaR, but only for same atom-atom channels
  i2all_1 = 0

  if (nstf .eq. nsti) then 
     lammin = abs(l1-l1p)
     lammax = l1+l1p
     
     do lam=lammin,lammax
        if( lam.gt.lamtop_vc) exit 
        if(Mod(lam,2) .eq. 1) cycle
        
        i1 = minvnc(lam)
        i2 = maxvnc(lam)
        i2all_1 = max(i2all_1,i2)

!!$ this is nuclear() call
        if( lam.eq.0 ) then    
           ttt(1:i2) = ttt(1:i2) + ( vnc(1:i2,lam) + Zasym/gridr(1:i2))
           if(l1.le.Ldw .and. iBorn.ne.1)  then
              ! Mark: Accounting for positron scattering
              !                 ttt(1:maxvdist) = ttt(1:maxvdist) - vdist(1:maxvdist)
              ttt(1:maxvdist) = ttt(1:maxvdist) + data_in%Zproj * vdist(1:maxvdist)
           else
!!$  This section to be modified for molecular targets that have frozen core 
              !                 la_vdcore = 0   ! to be modified for general case  !!!
              !                 ttt(1:maxvdc) = ttt(1:maxvdc) + ang1*vdcore(1:maxvdc,la_vdcore) !  check if lam should be chaged to la_vdcore
              !                
           endif
        else
           ang1 = Yint(dble(l1),dble(m1),dble(lam),dble(m1-m1p),dble(l1p),dble(m1p))
           if( ang1 .eq. 0) cycle
           ttt(i1:i2) = ttt(i1:i2) + ang1 * vnc(i1:i2,lam)
        endif
           
     enddo ! end lam loop
     
  endif
  
  
  if( i2all .eq. 0 .and. i2all_1 .eq. 0) then
     return    !  this ME is zero just due to angular momentum, it never got inside lam loop
  endif
  
  call minmaxi(ttt,nr,i1,i2)
  i2all = i2
  if(i2 .lt. nr) then
     ttt(i2+1:nr) = 0d0
  endif


!!$ -----------------------------------------------------------------------------------------
!!$  print*,'nsti = ',nsti
!!$  print*, 'nstf = ',nstf
!!$ -------------------------------------------------------------------------------
!!$ Test ttt values here
!!$  print*,'I am in Ven'
!!$  print*,'ttt(200) = ',ttt(200)
!!$ -------------------------------------------------------------------------------
!!$     i = 200
!!$     rval = gridr(i)
!!$    Vval11 = -exp(-2d0*2.0*rval)*(2d0 + 1d0/rval)
!!$    Vval12 =exp(-3d0*rval)*((16d0+48d0*rval)/(27d0*sqrt(2d0)))
!!$     print*

!!$     if( nsti .eq. 1 .and. nstf .eq. 1) then
!!$        print*, '>>>>', 'i = ',i,'rval = ',rval,'ttt(i) = ', ttt(i), 'Vval11 = ', Vval11
!!$     elseif(nsti .eq. 1 .and. nstf .eq. 2) then
!!$          print*, '>>>>', 'i = ',i,'rval = ',rval,'ttt(i) = ', ttt(i), 'Vval12 = ', Vval12
!!$     end if

!!$     print*
!!$ ------------------------------------------------------------------------------------------
  do ki=1,nqmi
     kqi = npk_nchi + ki - 1
     n1p = kqi
     pn1p => chil%b(n1p)
     f1p => fpointer(pn1p)
     minf1p = get_minf(pn1p)
     maxf1p = get_maxf(pn1p)
     
     do kf=1,nqmf
        kqf = npk_nchf + kf - 1
        n1 = kqf
        pn1 => chil%b(n1)
        f1 => fpointer(pn1)
        minf1 = get_minf(pn1)
        maxf1 = get_maxf(pn1)
        
        maxfm = min(maxf1,maxf1p) ! nr as for continuum functions, but for closed channels maxf = 1
        minfm = max(minf1,minf1p)
        
        maxfm = min(maxfm,i2all)
        
        fun1(minfm:maxfm) = (f1(minfm:maxfm)*f1p(minfm:maxfm))*weight(minfm:maxfm) * ttt(minfm:maxfm)
        
        sum_rad = SUM(fun1(minfm:maxfm))
                
        do lam=lammin,lammax
           
!!$ Do not calculate core-exchange for Born ME
           tmp_core_exch = 0d0
           if( lam .eq. 0 .and. n2.eq.n2p .and. iBorn .eq. 0) then
              !              call fcexch_rel_st(pn1,pn1p,tmp_core_exch)
              sum_rad = sum_rad + tmp_core_exch
           endif
           
        enddo   ! end lam loop

! Mark: Accounting for positron scattering.         
!        vmatt(kf,ki) =  vmatt(kf,ki) + sum_rad * (2d0/pi) !!! / (rkf*rki) dealt with in scat.f
        vmatt(kf,ki) = vmatt(kf,ki)  - data_in%Zproj *  sum_rad * (2d0/pi)


     enddo    ! end kf loop
  enddo    ! end ki loop


end subroutine vdirect


!!$------------------------------------------------------------------------------------------
!!$ Make distorting potential


subroutine make_vdist(nst,vdist,maxvdist)
! This is adaptation of routine for matrix element of V(1,2) electron-electron potential for 4 functions:
! <n1 n2 : M | V(1,2) | n1p n2p : M>, where V(1,2) = sum_{lam} V_{lam}(1,2), and J is total anguar moment of the two electrons.
! Note: this is not reduced ME, < || V || > = hat(J) < | V | >
! For lam = 0

  use input_data
  use one_electron_func
  use sturmian_class
  use target_states
  use channels
  use grid_radial
  use vnc_module
  use MPI_module

  implicit none

  integer, intent(in):: nst
  real*8, dimension(grid%nr), intent(out):: vdist
  integer, intent(out):: maxvdist

  type(state), pointer:: state_f,state_i
  integer:: nstf, nsti, Mfst, Mist, ncf,nci, ncfm, ncim
  integer::  mst, m2p, l2, l2p
  integer:: n1,n2,n1p,n2p
  real*8:: rCIf, rCIi
  real*8, pointer, dimension(:):: f2, f2p
  real*8, dimension(grid%nr):: gridr, weight, temp, fun, ttt
  integer:: minvdist, lam
  type(sturmian_nr), pointer:: pn2, pn2p  !  one-electron states
  integer:: nr

  integer:: minfun, maxfun  
  integer:: minf2,maxf2, minf2p,maxf2p
  integer:: minf, maxf, minfp, maxfp,  ir1, ir2

  ! H2-like 
  real*8, pointer, dimension(:):: f1, fp1
  real*8:: CI_stp, CI_st, CI_1, CI_2, CI_1p, CI_2p, temp_CI ! 2e configuration and 1e molecular orbial coefficien
  real*8:: temp_overlap2, overlap2
  ! Configuration Loops
  integer::  nst_con,  nstp_con, Ncon, Npcon        ! 2e configuration loops
  integer::  n_orb, np_orb, n_orb_max, np_orb_max   ! loops nuse 
  integer:: use_norb, use_nporb                     ! index to 1e molecular orbital
  integer:: MO1, MO2, MOp1, MOp2                ! index to 1e molecular orbital
  integer:: i1, i2, ip1, ip2, nam1, nam2, namp1, namp2  ! Loop over 1e molecular orbitals atomic orbitals
  integer:: ind1, ind2, indp1, indp2                          ! Index to underlying atomic orbitals
  type(sturmian_nr), pointer:: tn1, tn2, tnp1, tnp2  ! One-electron orbitals
  integer:: Lap1,Mp1,Lap2,Mp2,La1,M1,La2,M2    ! A.O. QN and Molecular State M 


  nr = grid%nr
  weight(1:nr) = grid%weight(1:nr)
  gridr(1:nr) = grid%gridr(1:nr)

  ttt(:) = 0d0
  temp(:) = 0d0

  lam = 0
  vdist(:) = 0d0

  nstf = nst
  nsti = nst
  state_f => Target_Basis%b(nstf)
  state_i => Target_Basis%b(nsti)
  Mfst = get_ang_mom_proj(state_f)
  Mist = get_ang_mom_proj(state_i)

  if ( data_in%hlike ) then
  ncfm = get_nam(state_f)
  ncim = get_nam(state_i)


  do ncf=1,ncfm
     M2 = get_ma(state_f,ncf,1)
     n2 = get_na(state_f,ncf,1)
     pn2 => bst_nr%b(n2)
     l2 = get_ang_mom(pn2)

     if(M2 .ne. Mfst) then
        if (myid == 0 ) print*, 'vdist >>>>>>>>>> M2 /= Mstf: ', nstf, Mfst, M2, l2, n2
        stop
     endif

     f2 => fpointer(pn2)
     minf2 = get_minf(pn2)
     maxf2 = get_maxf(pn2)

     rCIf = get_CI(state_f,ncf)

     do nci=1,ncim
        m2p = get_ma(state_i,nci,1)
        n2p = get_na(state_i,nci,1)
        pn2p => bst_nr%b(n2p)
        l2p = get_ang_mom(pn2p)

        if(m2p .ne. Mist) then
           if (myid == 0 ) print*, 'vdist >>>>>>>>>>> m2p != Msti: ',  nsti, Mist, m2p, l2p, n2p
           stop
        endif
        
        if(l2 .ne. l2p) cycle ! this  is due to ang2: for lam=0 only l2=l2p terms are not zero

        f2p => fpointer(pn2p)
        minf2p = get_minf(pn2p)
        maxf2p = get_maxf(pn2p)

        rCIi = get_CI(state_i,nci)

        ir1 = max(minf2,minf2p)
        ir2 = min(maxf2,maxf2p)

        fun(:) = 0d0
        fun(ir1:ir2) = rCIf*rCIi*(f2(ir1:ir2)*f2p(ir1:ir2)) * weight(ir1:ir2)
        
        call form(lam,fun,ir1,ir2,nr,temp,minfun,maxfun)
                
        ttt(minfun:maxfun) = ttt(minfun:maxfun) + temp(minfun:maxfun)

     enddo

  enddo

  else ! H2-like

  ! Final State number nstf
  Npcon =  TargetStates2el%b(nstf)%nam          ! Number of 2e configurations.        
  np_orb_max = TargetStates2el%b(nstf)%nusemax  ! Number of orbitals used to describe this state
  ! Initial State number nsti
  Ncon = TargetStates2el%b(nsti)%nam
  n_orb_max = TargetStates2el%b(nsti)%nusemax

  do np_orb = 1, np_orb_max
     use_nporb = TargetStates2el%b(nstf)%nuse(np_orb)
!!$ INITIAL State MO COORDINATE 1
     do n_orb = 1, n_orb_max
        use_norb =  TargetStates2el%b(nsti)%nuse(n_orb)   
        overlap2 = 0d0
!!$ Looping over FINAL (p) Molecular State configorations 
        do nstp_con =  1, Npcon          
           MOp1 = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
           if ( MOp1 /= use_nporb ) cycle   ! 
           CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con) ! 2e config CI
           MOp2 = get_na(TargetStates2el%b(nstf),nstp_con,2) ! Index to 1e molecular orbital
           Mp2 = get_ang_mom_proj(TargetStates%b(MOp2))      ! Angular momentum projection 1e molecular orbital
           namp2 = get_nam(TargetStates%b(MOp2))    ! Number of atomic orbitals that make up molecular orbital 
!!$ Looping over INITIAL Molecular State configurations 
           do nst_con =  1, Ncon
              MO1 = get_na(TargetStates2el%b(nsti),nst_con,1)
              if ( MO1 /= use_norb ) cycle
              CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
              MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)
              M2 = get_ang_mom_proj(TargetStates%b(MO2))
              nam2 = get_nam(TargetStates%b(MO2))
              ! overlap 2 matrix element <varphi'_2|varphi_2> 
              if ( Mp2 /= M2 ) cycle  
!!$ FINAL STATE COORDINATE 2 BETA
              do ip2 = 1, namp2
                 indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
                 tnp2 => bst_nr%b(indp2)
                 Lap2 = get_ang_mom(tnp2)
                 CI_2p = get_CI(TargetStates%b(MOp2),ip2)
!!$ INITIAL STATE COORDINATE 2 DELTA
                 do i2 = 1, nam2 
                    ind2 = get_na(TargetStates%b(MO2),i2)  
                    tn2 => bst_nr%b(ind2)                                         
                    La2 = get_ang_mom(tn2)                            
                    CI_2 = get_CI(TargetStates%b(MO2),i2)        
                    ! overlap 2 matrix element <varphi'_2|varphi_2> 
                    if ( Lap2 /= La2) cycle
                    temp_CI = CI_stp * CI_2p * CI_st * CI_2 
                    if ( temp_CI == 0d0 ) cycle 
                    temp_overlap2 = bst_nr%ortint(indp2,ind2)
                    overlap2 = overlap2 + temp_CI * temp_overlap2
                 end do ! i2        
              end do ! ip2   
           end do ! nst_con
        end do ! nstp_con   
        if ( overlap2 == 0d0 ) cycle
        MOp1 = use_nporb ! Index 1e molecular orbital
        namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
        Mp1 = get_ang_mom_proj(TargetStates%b(MOp1))
        MO1 = use_norb
        nam1 = get_nam(TargetStates%b(MO1))
        M1 = get_ang_mom_proj(TargetStates%b(MO1))

!!$ lambda = 0 therefore mu=0 and Mp1 == M1
	if ( Mp1 /= M1 ) cycle
!!$ COORDINATE 1 RADIAL INTEGRALS ttt(r_0) = 2.< |V01| >.< | >
!!$ Quantum numbers and functions for ALPHA
        do ip1 = 1, namp1
           indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
           tnp1 => bst_nr%b(indp1)   !           
           fp1 => fpointer(tnp1)     ! One electron functions
           Lap1 = get_ang_mom(tnp1)  ! Gets Angular momentum A.O.
           minfp = get_minf(tnp1)
           maxfp = get_maxf(tnp1)
           CI_1p = get_CI(TargetStates%b(MOp1),ip1)
!!$ Quantum numbers and functions for GAMMA
           do i1 = 1, nam1
              ind1 = get_na(TargetStates%b(MO1),i1) ! Index to underlying atomic orbitals  
              tn1 => bst_nr%b(ind1)
              f1 => fpointer(tn1)
              La1 = get_ang_mom(tn1)
              minf = get_minf(tn1)
              maxf = get_maxf(tn1)
              CI_1 = get_CI(TargetStates%b(MO1),i1)
!!$ Electron-electron potential 2*V_{01}
              temp_CI = CI_1p * CI_1
              if ( temp_CI == 0d0 ) cycle    
!!$ lambda = 0 therefore Lap1 == La1
	      if ( Lap1 /= La1 ) cycle
              ir1 = max(minf,minfp)
              ir2 = min(maxf,maxfp)
              fun(:) = 0d0
!!$ Multiply the electron-electron potential by 2. Refer D. Fursa e-He 95 Paper
              !fun(ir1:ir2) = 2d0 * temp_CI * (f1(ir1:ir2) * fp1(ir1:ir2) * weight(ir1:ir2)) * overlap2
              fun(ir1:ir2) = 2d0 * temp_CI * (f1(ir1:ir2) * fp1(ir1:ir2)) * overlap2

              !call form(lam, fun, ir1, ir2, nr, temp, minfun, maxfun)   
              call form_accurate(lam, fun, ir1, ir2, nr, temp, minfun, maxfun)   
              ttt(minfun:maxfun) = ttt(minfun:maxfun) + temp(minfun:maxfun) 
                 
           end do ! i1
        end do ! ip1
     end do ! n_orb use_n
  end do ! np_orb use_np

  end if ! Target


!!$  Subtract nuclear potential VLambdaR (only lam=0 part), but only for same atom-atom channels
!!$  Later on make sur ethat vdcore()  is included 

!!$ Note: vnc = -2Z/r
  temp(1:nr) = vnc(1:nr,lam) + data_in%Zasym/gridr(1:nr) + ttt(1:nr) 

!  ttt(1:nr) = temp(1:nr) + ttt(1:nr) ! it has 1/r form at large r values

 ! call minmaxi(ttt,nr,i1,i2)
  call minmaxi(temp,nr,i1,i2)

!  vdist(i1:i2) = - data_in%Zproj * ttt(i1:i2)
  vdist(i1:i2) = - data_in%Zproj * temp(i1:i2)

  minvdist = 1
  maxvdist = i2

  open(unit=1000,file='vdist.out',action='write',status='replace')
  do i1=minvdist, maxvdist
    !write(1000,'(1000(ES17.10,2X))') gridr(i1), vnc(i1,lam), data_in%Zasym/gridr(i1), ttt(i1), vdist(i1)
    write(1000,*) gridr(i1), vnc(i1,lam), data_in%Zasym/gridr(i1), ttt(i1), vdist(i1)
  enddo
  close(1000)

end subroutine make_vdist

!!$------------------------------------------------------------------------------------------
!!$------------------------------------------------------------------------------------------

subroutine vexch(Mtot,nchf,nchi,nstf,nsti,nqmf,nqmi,npk_nchf,npk_nchi,chil,bst_nr,nr,weight,gridr,ncwaves,Nmax,ortchil,flchil,Etot,theta,vmatt,Rd,rspin_tot)

! This is adaptation of routine for matrix element of V(1,2) electron-electron potential for 4 functions: 
! <k_f Phi_f: Mtot | V_{0,1})P_(r0,r1) | k_i Phi_i : Mtot>, where V(1,2) = sum_{lam} V_{lam}(1,2), and Mtot is total anguar moment projection of the two electrons.

  use sturmian_class
  use target_states
  use channels
  use input_data
  use distpot
  use vnc_module
  use pol_2el
  use MPI_module

  implicit none
  
  integer, intent(in):: Mtot    ! total ang.momentum projection of two electrons - integer!
  integer, intent(in):: nchf,nchi
  integer, intent(in):: nstf,nsti
  integer, intent(in):: nqmf, nqmi
  integer, intent(in):: npk_nchf,npk_nchi
  type(basis_sturmian_nr):: chil, bst_nr
  integer, intent(in):: nr
  real*8, dimension(nr), intent(in):: weight, gridr
  integer, intent(in)::  ncwaves
  real*8, intent(in):: Etot, theta
  real*8, dimension(nqmf,nqmi), intent(inout):: vmatt
  real*8, dimension(Nmax,ncwaves), intent(in):: flchil, ortchil
  real*8, intent(in):: rspin_tot  ! For projection matrix. No dependence on spin, hence must cancel.
  integer, intent(in)::  Nmax ! Nmax Number of Molecular States
  real*8, intent(in):: Rd

  real*8:: pi, t1el, sum_rad
  real*8:: Yint,  ang0, ang1, tmp0 
  integer:: lam, lammin, lammax, mu, ir1, ir2, fr1, fr2, jr1, jr2, itest
  real*8, dimension(nr):: temp, fun, fun1,  ttt

!!$ Molecular States and Atomic Orbitals
  type(state), pointer:: state_f,state_i
  integer:: f,i, nfao, niao, fao, iao
  type(sturmian_nr), pointer:: faop, iaop
  real*8, pointer, dimension(:):: fao_f, iao_f
  integer:: Mfst, Mist, fao_l, fao_m, iao_l, iao_m
  real*8::  CIf, CIi
  integer:: minf_fao, maxf_fao, minf_iao, maxf_iao

!!$ Projectile
  integer:: kf, ki, kqf, kqi
  real*8:: kf_en, ki_en 
  type(sturmian_nr), pointer:: kfp, kip
  real*8, pointer, dimension(:):: kf_f, ki_f
  integer:: kfL, kiL, kfM, kiM
  integer:: minf_kf, maxf_kf,  minf_ki, maxf_ki

!!$ Identity Operator
  type(state), pointer:: nstatep
  integer::  nstate, nstate_M
  real*8:: Ident

  real*8 :: Z1, Z2

  Z1 = data_in%Z1
  Z2 = data_in%Z2

  pi = acos(-1d0)
  
  kfL = Lp_ch(nchf)   ! projectile
  kfM = Mp_ch(nchf)
  kiL = Lp_ch(nchi)   ! projectile
  kiM = Mp_ch(nchi)
  
  state_f => TargetStates%b(nstf)
  state_i => TargetStates%b(nsti)
  
  nfao = get_nam(state_f)
  niao = get_nam(state_i)
  
  Mfst = get_ang_mom_proj(state_f)
  Mist = get_ang_mom_proj(state_i)
  
  
  if( kfM + Mfst  .ne. kiM + Mist ) return
  
!!$ n2, n2p are indexes to one-electron target orbitals
!!$ while n1, n1p are indexes to continuum waves
  
  do kf = 1, nqmf
     kqf = npk_nchf + kf - 1
     kfp => chil%b(kqf)
     kf_en = get_energy(kfp)
     
     do ki = 1, nqmi
        kqi = npk_nchi + ki - 1
        kip => chil%b(kqi)
        ki_en = get_energy(kip)    
        
!!$ Additions for projection operator E.I_{0}.theta and [E(1-theta)-H]Pr. This needs lots of checking 
!!$ This is projection operator term  E.I_{0}.theta.           
        if ( nsti == nstf ) then
           Ident = 0d0
           do nstate = 1, Nmax  ! Loop over all molecular states      
              nstatep => TargetStates%b(nstate)
              nstate_M = get_ang_mom_proj(nstatep)
              
              ! <phi|kLM> should have same M ortherwise zero. 
              if ( nstate_M /= kiM .OR. nstate_M /= kfM ) cycle             
              
              Ident = Ident + ortchil(nstate,kqi) * ortchil(nstate,kqf)
              
           end do
           Ident = dble((-1)**(NINT(rspin_tot))) *  Etot * theta * Ident 
           
           vmatt(kf,ki) =  vmatt(kf,ki) +  Ident * (2d0/pi)
        end if

        t1el = 0d0
!!$ One-electron exchange terms V_{0}-U_{0}
        t1el = t1el - flchil(nsti,kqf) * ortchil(nstf,kqi)
!!$ One-electron exchange terms V_{1}-U_{1}
        t1el = t1el - flchil(nstf,kqi) * ortchil(nsti,kqf)
!!$ (E(1-theta) + en_ki + en_kf) + Z^2/R
        t1el = t1el - ( -Etot*(1d0 - theta) + kf_en + ki_en) * ortchil(nsti,kqf) * ortchil(nstf,kqi)
        if ( Rd /= 0d0 ) then
           t1el = t1el - ortchil(nsti,kqf) * ortchil(nstf,kqi) * Z1 * Z2 / Rd
        end if
        vmatt(kf,ki) =  vmatt(kf,ki) + t1el * (2d0 /pi)           
               
     end do ! ki loop
  end do ! kf loop
  
!!$ Electron-Electron V_{01} term 
  mu =  Mfst - kiM 
  
  ! COORDINATE 1 FINAL STATE ALPHA
  do f = 1, nfao
     fao = get_na(state_f,f,1)
     CIf = get_CI(state_f,f)
     faop => bst_nr%b(fao)
     fao_l = get_ang_mom(faop)
     fao_m = get_ma(state_f,f,1)
     
     fao_f => fpointer(faop)
     minf_fao = get_minf(faop)
     maxf_fao = get_maxf(faop)
     
     do ki = 1, nqmi
        kqi = npk_nchi + ki - 1
        kip => chil%b(kqi)
        ki_f => fpointer(kip)
        minf_ki = get_minf(kip)
        maxf_ki = get_maxf(kip)
        
        fr1 = max(minf_fao,minf_ki)
        fr2 = min(maxf_fao,maxf_ki)
        
        lammin = abs(fao_l - kiL)
        lammax = fao_l + kiL
        if ( lammin > lammax ) cycle
        
        fun(:) = 0d0
        fun(fr1:fr2) =  fao_f(fr1:fr2)  * ki_f(fr1:fr2)  * weight(fr1:fr2)       
        
        do lam = lammin, lammax
           
           ang1 = Yint(dble(fao_l),dble(fao_M),dble(lam),dble(mu),dble(kiL),dble(kiM))
           if( ang1 == 0d0) cycle
           
!!$  Test if there is a nonzero angular coef. for given value of lam.
           itest = 0
           ! COORDINATE 0 INITIAL STATE BETA
           do i = 1, niao
              iao = get_na(state_i,i,1)
              iaop => bst_nr%b(iao)
              iao_l = get_ang_mom(iaop)
              iao_m = get_ma(state_i,i,1)             
              ang0 = dble((-1)**(mu))*Yint(dble(kfL),dble(kfM),dble(lam),-dble(mu),dble(iao_l),dble(iao_m))
              
              if( ang0 /= 0d0) then
                 itest = 1
                 exit
              end if
           end do
           if (itest == 0) cycle
           
           call form(lam,fun,fr1,fr2,nr,temp,jr1,jr2)
           
           ttt(jr1:jr2) =  temp(jr1:jr2)
           
           do i = 1, niao
              iao = get_na(state_i,i,1)
              CIi = get_CI(state_i,i)
              iaop => bst_nr%b(iao)
              iao_l = get_ang_mom(iaop)
              iao_m = get_ma(state_i,i,1)
              
              ang0 = dble((-1)**(mu))*Yint(dble(kfL),dble(kfM),dble(lam),-dble(mu),dble(iao_l),dble(iao_m))
              
              if( ang0 == 0d0) cycle   
              
              iao_f => fpointer(iaop)
              minf_iao = get_minf(iaop)
              maxf_iao = get_maxf(iaop)
              
              tmp0 = CIf * CIi * ang0 * ang1
              if( tmp0 == 0d0) cycle 
              
              do kf = 1, nqmf
                 kqf = npk_nchf + kf - 1
                 kfp => chil%b(kqf)
                 kf_f => fpointer(kfp)
                 minf_kf = get_minf(kfp)
                 maxf_kf = get_maxf(kfp)
                 
                 ir1 = max(minf_iao,minf_kf)
                 ir2 = min(maxf_iao,maxf_kf)
                 
                 ir1 = max(ir1,jr1)
                 ir2 = min(ir2,jr2)
                 
                 fun1(ir1:ir2) = (kf_f(ir1:ir2) * iao_f(ir1:ir2)) * weight(ir1:ir2) * ttt(ir1:ir2) 
                 
                 sum_rad = -SUM(fun1(ir1:ir2))
                 
                 vmatt(kf,ki) =  vmatt(kf,ki) + sum_rad * tmp0 * (2d0/pi) !!! / (rkf*rki) dealt with in scat.f
                 
              end do    ! end kf loop

           end do  ! end Initial State AO loop
                      
        end do  ! end lam loop

     end do    ! end ki loop
  end do  ! end Final State AO loop

  vmatt(:,:) = -dble((-1)**(NINT(rspin_tot))) * vmatt(:,:)
  
  
end subroutine vexch

!Liam - this is a copy of the make_vdist subroutine with the old switch between MO/SP representation
subroutine make_vdist_old(nst,vdist,maxvdist)
! This is adaptation of routine for matrix element of V(1,2) electron-electron potential for 4 functions:
! <n1 n2 : M | V(1,2) | n1p n2p : M>, where V(1,2) = sum_{lam} V_{lam}(1,2), and J is total anguar moment of the two electrons.
! Note: this is not reduced ME, < || V || > = hat(J) < | V | >
! For lam = 0

  use input_data
  use one_electron_func
  use sturmian_class
  use target_states
  use channels
  use grid_radial
  use vnc_module
  use MPI_module

  implicit none

  integer, intent(in):: nst
  real*8, dimension(grid%nr), intent(out):: vdist
  integer, intent(out):: maxvdist

  type(state), pointer:: state_f,state_i
  integer:: nstf, nsti, Mfst, Mist, ncf,nci, ncfm, ncim
  integer::  mst, m2p, l2, l2p
  integer:: n1,n2,n1p,n2p
  real*8:: rCIf, rCIi
  real*8, pointer, dimension(:):: f2, f2p
  real*8, dimension(grid%nr):: gridr, weight, temp, fun, ttt
  integer:: minvdist, lam
  type(sturmian_nr), pointer:: pn2, pn2p  !  one-electron states
  integer:: nr

  integer:: minfun, maxfun  
  integer:: minf2,maxf2, minf2p,maxf2p
  integer:: minf, maxf, minfp, maxfp,  ir1, ir2

  ! H2-like 
  real*8, pointer, dimension(:):: f1, fp1
  real*8:: CI_stp, CI_st, CI_1, CI_2, CI_1p, CI_2p, temp_CI ! 2e configuration and 1e molecular orbial coefficien
  real*8:: temp_overlap2, overlap2
  ! Configuration Loops
  integer::  nst_con,  nstp_con, Ncon, Npcon        ! 2e configuration loops
  integer::  n_orb, np_orb, n_orb_max, np_orb_max   ! loops nuse 
  integer:: use_norb, use_nporb                     ! index to 1e molecular orbital
  integer:: MO1, MO2, MOp1, MOp2                ! index to 1e molecular orbital
  integer:: i1, i2, ip1, ip2, nam1, nam2, namp1, namp2  ! Loop over 1e molecular orbitals atomic orbitals
  integer:: ind1, ind2, indp1, indp2                          ! Index to underlying atomic orbitals
  type(sturmian_nr), pointer:: tn1, tn2, tnp1, tnp2  ! One-electron orbitals
  integer:: Lap1,Mp1,Lap2,Mp2,La1,M1,La2,M2    ! A.O. QN and Molecular State M 


  nr = grid%nr
  weight(1:nr) = grid%weight(1:nr)
  gridr(1:nr) = grid%gridr(1:nr)

  ttt(:) = 0d0
  temp(:) = 0d0

  lam = 0
  vdist(:) = 0d0

  nstf = nst
  nsti = nst
  state_f => Target_Basis%b(nstf)
  state_i => Target_Basis%b(nsti)
  Mfst = get_ang_mom_proj(state_f)
  Mist = get_ang_mom_proj(state_i)

  if ( data_in%hlike ) then
  ncfm = get_nam(state_f)
  ncim = get_nam(state_i)


  do ncf=1,ncfm
     M2 = get_ma(state_f,ncf,1)
     n2 = get_na(state_f,ncf,1)
     pn2 => bst_nr%b(n2)
     l2 = get_ang_mom(pn2)

     if(M2 .ne. Mfst) then
        if (myid == 0 ) print*, 'vdist >>>>>>>>>> M2 /= Mstf: ', nstf, Mfst, M2, l2, n2
        stop
     endif

     f2 => fpointer(pn2)
     minf2 = get_minf(pn2)
     maxf2 = get_maxf(pn2)

     rCIf = get_CI(state_f,ncf)

     do nci=1,ncim
        m2p = get_ma(state_i,nci,1)
        n2p = get_na(state_i,nci,1)
        pn2p => bst_nr%b(n2p)
        l2p = get_ang_mom(pn2p)

        if(m2p .ne. Mist) then
           if (myid == 0 ) print*, 'vdist >>>>>>>>>>> m2p != Msti: ',  nsti, Mist, m2p, l2p, n2p
           stop
        endif
        
        if(l2 .ne. l2p) cycle ! this  is due to ang2: for lam=0 only l2=l2p terms are not zero

        f2p => fpointer(pn2p)
        minf2p = get_minf(pn2p)
        maxf2p = get_maxf(pn2p)

        rCIi = get_CI(state_i,nci)

        ir1 = max(minf2,minf2p)
        ir2 = min(maxf2,maxf2p)

        fun(:) = 0d0
        fun(ir1:ir2) = rCIf*rCIi*(f2(ir1:ir2)*f2p(ir1:ir2)) * weight(ir1:ir2)
        
        call form(lam,fun,ir1,ir2,nr,temp,minfun,maxfun)
                
        ttt(minfun:maxfun) = ttt(minfun:maxfun) + temp(minfun:maxfun)

     enddo

  enddo

  else ! H2-like

  ! Final State number nstf
  Npcon =  TargetStates2el%b(nstf)%nam          ! Number of 2e configurations.        
  np_orb_max = TargetStates2el%b(nstf)%nusemax  ! Number of orbitals used to describe this state
  ! Initial State number nsti
  Ncon = TargetStates2el%b(nsti)%nam
  n_orb_max = TargetStates2el%b(nsti)%nusemax

  if (data_in%non_uniq ) then

  do np_orb = 1, np_orb_max
     use_nporb = TargetStates2el%b(nstf)%nuse(np_orb)
!!$ INITIAL State MO COORDINATE 1
     do n_orb = 1, n_orb_max
        use_norb =  TargetStates2el%b(nsti)%nuse(n_orb)   
        overlap2 = 0d0
!!$ Looping over FINAL (p) Molecular State configorations 
        do nstp_con =  1, Npcon          
           MOp1 = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
           if ( MOp1 /= use_nporb ) cycle   ! 
           CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con) ! 2e config CI
           MOp2 = get_na(TargetStates2el%b(nstf),nstp_con,2) ! Index to 1e molecular orbital
           Mp2 = get_ang_mom_proj(TargetStates%b(MOp2))      ! Angular momentum projection 1e molecular orbital
           namp2 = get_nam(TargetStates%b(MOp2))    ! Number of atomic orbitals that make up molecular orbital 
!!$ Looping over INITIAL Molecular State configurations 
           do nst_con =  1, Ncon
              MO1 = get_na(TargetStates2el%b(nsti),nst_con,1)
              if ( MO1 /= use_norb ) cycle
              CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
              MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)
              M2 = get_ang_mom_proj(TargetStates%b(MO2))
              nam2 = get_nam(TargetStates%b(MO2))
              ! overlap 2 matrix element <varphi'_2|varphi_2> 
              if ( Mp2 /= M2 ) cycle  
!!$ FINAL STATE COORDINATE 2 BETA
              do ip2 = 1, namp2
                 indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
                 tnp2 => bst_nr%b(indp2)
                 Lap2 = get_ang_mom(tnp2)
                 CI_2p = get_CI(TargetStates%b(MOp2),ip2)
!!$ INITIAL STATE COORDINATE 2 DELTA
                 do i2 = 1, nam2 
                    ind2 = get_na(TargetStates%b(MO2),i2)  
                    tn2 => bst_nr%b(ind2)                                         
                    La2 = get_ang_mom(tn2)                            
                    CI_2 = get_CI(TargetStates%b(MO2),i2)        
                    ! overlap 2 matrix element <varphi'_2|varphi_2> 
                    if ( Lap2 /= La2) cycle
                    temp_CI = CI_stp * CI_2p * CI_st * CI_2 
                    if ( temp_CI == 0d0 ) cycle 
                    temp_overlap2 = bst_nr%ortint(indp2,ind2)
                    overlap2 = overlap2 + temp_CI * temp_overlap2
                 end do ! i2        
              end do ! ip2   
           end do ! nst_con
        end do ! nstp_con   
        if ( overlap2 == 0d0 ) cycle
        MOp1 = use_nporb ! Index 1e molecular orbital
        namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
        Mp1 = get_ang_mom_proj(TargetStates%b(MOp1))
        MO1 = use_norb
        nam1 = get_nam(TargetStates%b(MO1))
        M1 = get_ang_mom_proj(TargetStates%b(MO1))

!!$ lambda = 0 therefore mu=0 and Mp1 == M1
	if ( Mp1 /= M1 ) cycle
!!$ COORDINATE 1 RADIAL INTEGRALS ttt(r_0) = 2.< |V01| >.< | >
!!$ Quantum numbers and functions for ALPHA
        do ip1 = 1, namp1
           indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
           tnp1 => bst_nr%b(indp1)   !           
           fp1 => fpointer(tnp1)     ! One electron functions
           Lap1 = get_ang_mom(tnp1)  ! Gets Angular momentum A.O.
           minfp = get_minf(tnp1)
           maxfp = get_maxf(tnp1)
           CI_1p = get_CI(TargetStates%b(MOp1),ip1)
!!$ Quantum numbers and functions for GAMMA
           do i1 = 1, nam1
              ind1 = get_na(TargetStates%b(MO1),i1) ! Index to underlying atomic orbitals  
              tn1 => bst_nr%b(ind1)
              f1 => fpointer(tn1)
              La1 = get_ang_mom(tn1)
              minf = get_minf(tn1)
              maxf = get_maxf(tn1)
              CI_1 = get_CI(TargetStates%b(MO1),i1)
!!$ Electron-electron potential 2*V_{01}
              temp_CI = CI_1p * CI_1
              if ( temp_CI == 0d0 ) cycle    
!!$ lambda = 0 therefore Lap1 == La1
	      if ( Lap1 /= La1 ) cycle
              ir1 = max(minf,minfp)
              ir2 = min(maxf,maxfp)
              fun(:) = 0d0
!!$ Multiply the electron-electron potential by 2. Refer D. Fursa e-He 95 Paper
              fun(ir1:ir2) = 2d0 * temp_CI * (f1(ir1:ir2) * fp1(ir1:ir2) * weight(ir1:ir2)) * overlap2

              call form(lam, fun, ir1, ir2, nr, temp, minfun, maxfun)   
              ttt(minfun:maxfun) = ttt(minfun:maxfun) + temp(minfun:maxfun) 
                 
           end do ! i1
        end do ! ip1
     end do ! n_orb use_n
  end do ! np_orb use_np

  else ! Not solving non-uniqueness: AO representation 
 
  ! FINAL State COORDINATE 1
  do np_orb = 1, np_orb_max
     use_nporb = TargetStates2el%b(nstf)%nuse(np_orb)
     ! INITIAL State COORDINATE 1
     do n_orb = 1, n_orb_max
        use_norb =  TargetStates2el%b(nsti)%nuse(n_orb)
        overlap2 = 0d0
        ! Looping over FINAL Molecular State orbitals. COORDINATE 2
        do nstp_con =  1, Npcon
           if ( use_nporb /= TargetStates2el%b(nstf)%na(nstp_con) ) cycle
           ! Quantum numbers and functions for BETA
           indp2 = TargetStates2el%b(nstf)%nb(nstp_con)
           tnp2 => bst_nr%b(indp2)
           Lap2 = get_ang_mom(tnp2)
           Mp2 = get_ang_mom_proj(tnp2)
           CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con)
           ! Looping over INITIAL Molecular State orbitals.  COORDINATE 2
           do nst_con = 1, Ncon
              if ( use_norb /= TargetStates2el%b(nsti)%na(nst_con)  ) cycle
              ! Quantum numbers and functions for DELTA
              ind2 = TargetStates2el%b(nsti)%nb(nst_con)
              tn2 => bst_nr%b(ind2)
              La2 = get_ang_mom(tn2)
              M2 = get_ang_mom_proj(tn2)
              CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
              temp_CI = CI_stp * CI_st 
              if ( temp_CI == 0d0 ) cycle
              ! Selections Rules  
              if ( Lap2 /= La2 .OR. Mp2 /= M2  ) cycle
              temp_overlap2 = bst_nr%ortint(indp2,ind2)
              overlap2 = overlap2 +  temp_CI * temp_overlap2
           end do  ! INITIAL STATE COORDINATE 2
        end do    ! FINAL STATE COORDINATE 2       
        if ( overlap2 == 0d0 ) cycle
        ! COORDINATE 1 RADIAL INTEGRALS
        ! Quantum numbers and functions for ALPHA
        indp1 = use_nporb                                  ! Final state number np. nep A.O.
        tnp1 => bst_nr%b(indp1)                             !           
        fp1 => fpointer(tnp1)                               ! One electron functions
        Lap1 = get_ang_mom(tnp1)                            ! Gets Angular momentum A.O.
        Mp1 = get_ang_mom_proj(tnp1)                        ! Get angular projection of A.O. 
        minfp = get_minf(tnp1)
        maxfp = get_maxf(tnp1)
        ! Quantum numbers and functions for GAMMA
        ind1 = use_norb
        tn1 => bst_nr%b(ind1)
        f1 => fpointer(tn1)
        La1 = get_ang_mom(tn1)
        M1 = get_ang_mom_proj(tn1)
        minf = get_minf(tn1)
        maxf = get_maxf(tn1)

!!$ lambda = 0 therefore Lap1 == La1
        if ( Lap1 /= La1 .OR. Mp1 /= M1 ) cycle
        ir1 = max(minf,minfp)
        ir2 = min(maxf,maxfp)
        fun(:) = 0d0
!!$ Multiply the electron-electron potential by 2. Refer D. Fursa e-He 95 Paper
        fun(ir1:ir2) = 2d0 * (f1(ir1:ir2) * fp1(ir1:ir2) * weight(ir1:ir2)) * overlap2
        call form(lam, fun, ir1, ir2, nr, temp, minfun, maxfun)
        ttt(minfun:maxfun) = ttt(minfun:maxfun) + temp(minfun:maxfun)
     end do ! INITIAL STATE COORDINATE 1    
  end do  ! FINAL STATE COORDINATE 1



  end if ! non-uniquesness MO representation?


  end if ! Target


!!$  Subtract nuclear potential VLambdaR (only lam=0 part), but only for same atom-atom channels
!!$  Later on make sur ethat vdcore()  is included 

!!$ Note: vnc = -2Z/r
  temp(1:nr) = vnc(1:nr,lam) + data_in%Zasym/gridr(1:nr)   

  ttt(1:nr) = temp(1:nr) + ttt(1:nr) ! it has 1/r form at large r values

  call minmaxi(ttt,nr,i1,i2)

  vdist(i1:i2) = - data_in%Zproj * ttt(i1:i2)

  minvdist = 1
  maxvdist = i2

end subroutine make_vdist_old

!!$------------------------------------------------------------------------------------------
