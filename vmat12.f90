module vmat_subroutines

  contains

!Liam renamed this from vdirect_2e_MOrep to vdirect_2e since we always assume MO representation now
subroutine vdirect_2e(Mtot,rspin_tot, nchf,nchi,nstf,nsti,nqmf,nqmi,npk_nchf,npk_nchi,chil,nr,weight,gridr,Zproj,Zasym,iBorn,&
    &vmatt,numint)
  
  ! sum C^{np}_(np1,np2)  C^{n}_(n1,n2) < np0, np1, np2 : Mtot | V(0,1) | n0, n1, n2 : Mtot >
  ! Note: All integrals in target states coordinate space and projectile angular
  ! terms 
  ! multiply radial potentials are taken first. 
  ! This is because we need to eliminate the lambda = 0 term for the
  ! electron-electron 
  ! potential of the projectile and target states.
  ! V_{lambda=0}(0,1) = int (sin(kr))^2.phi_i.phi_f/r = oo. 
  ! Hence we eliminate the problem by using the projectile-nuclear term
  
  use sturmian_class
  use target_states
  use one_electron_func
  use ovlpste1me
  use channels
  use input_data
  use distpot
  use vnc_module
  use pol_2el
  implicit none
  integer, intent(in):: Mtot
  real*8, intent(in) :: rspin_tot
  integer, intent(in):: nchf, nchi        ! Channel number
  integer, intent(in):: nstf, nsti        ! State for particular channel
  integer, intent(in):: nqmf, nqmi        ! Number of k-grid points for particular channel
  integer, intent(in):: npk_nchf,npk_nchi ! On-shell k-grid point for particular channel
  integer, intent(in):: nr                ! Max r points
  real*8, dimension(nr), intent(in):: weight, gridr
  real*8, intent(in):: Zproj, Zasym      ! Z projectile charge, Z Aymptotic charge of target
  integer, intent(in):: iBorn             ! Born = 1
  real*8, dimension(nqmf,nqmi), intent(out):: vmatt
  real*8, optional, intent(inout) :: numint
  
  type(basis_sturmian_nr):: chil
  
  integer:: ki, kf, kqi, kqf
  real*8, dimension(nr):: ttt, temp, fun0, fun1
  
  ! Configuration Loops
  integer::  nst_con,  nstp_con, Ncon, Npcon        ! 2e configuration loops
  integer::  n_orb, np_orb, n_orb_max, np_orb_max   ! loops nuse 
  integer:: use_norb, use_nporb                     ! index to 1e molecular orbital
  integer:: MO1, MO2, MOp1, MOp2                ! index to 1e molecular orbital
  integer:: i, i1,i2, ip1, nam1,namp1  ! Loop over 1e molecular orbitals atomic orbitals
  
  type(sturmian_nr), pointer:: tn0,tnp0, tn1,tnp1  ! One-electron orbitals
  integer:: ind1, indp1                          ! Index to underlying atomic orbitals
  real*8, pointer, dimension(:)::  fp0, f0, fp1, f1  ! underlying atomic orbital radial functions
  ! Radial integrations
  integer::  minf, maxf, minfp, maxfp,  ir1, ir2, irt1,irt2
  integer::  ir2_nuc, ir2_e_e, ir_all
  real*8:: CI_stp, CI_st, CI_1, CI_2, CI_1p, CI_2p, temp_CI ! 2e configuration and 1e molecular orbial coefficients
  real*8:: sum_rad, pi, overlap2, temp_overlap2
  
  real*8:: Yint
  real*8:: ang0, ang1                                   ! Angular integrals for projectile and A.O.(1) 
  integer:: lambda, lamMin0,lamMax0, lamMin1,lamMax1, lambda_min, lambda_max, lambda_step, mu
  integer:: Lap0, La0, Mp0, M0                          ! Projectile Quantum Numbers (QN)  
  integer:: Lap1,Mp1,La1,M1,Mp12,M12    ! A.O. QN and Molecular State M 
  integer:: Spin, Spinp

  logical :: load

  integer :: icount
  icount = 0
  
  if(present(numint)) then
    load = .true.
    numint = 0
  else
    load = .false.
  endif

  vmatt(:,:) = 0d0
  pi = acos(-1d0)
  ir2_e_e = 0
  ir2_nuc = 0

  ! Projectile QN
  Lap0 = Lp_ch(nchf)
  Mp0 = Mp_ch(nchf)
  La0 = Lp_ch(nchi)
  M0 = Mp_ch(nchi)

  lamMin0 = abs(Lap0 - La0)
  lamMax0 = Lap0 + La0
  mu = M0 - Mp0

  
  ! Final State number nstf
  Npcon =  TargetStates2el%b(nstf)%nam          ! Number of 2e configurations.
  Mp12 = NINT(TargetStates2el%b(nstf)%M )       ! Molecular state Angular Projection
  Spinp = NINT(TargetStates2el%b(nstf)%spin)    !  Molecular State Spin
  np_orb_max = TargetStates2el%b(nstf)%nusemax  ! Number of orbitals used to describe this state
  
  ! Initial State number nsti
  Ncon = TargetStates2el%b(nsti)%nam         
  M12 = NINT(TargetStates2el%b(nsti)%M)       
  Spin = NINT(TargetStates2el%b(nsti)%spin)
  n_orb_max = TargetStates2el%b(nsti)%nusemax     
  
  ! Conserved QN. Selection Rules. Parity should be taken care of in channels.f90
  if (Spin /= Spinp) return   ! No direct interaction between states of different spins.
  if (rspin_tot==1.5d0 .and. Spin==0) return   ! Both states must be triplet for direct interaction in the S=3/2 partial wave.
  if( Mp12 + Mp0  /=  M12 + M0 ) return
  
  ttt(:) = 0d0


  ! FINAL (p) State Molecular Orbital (MO) COORDINATE 1
!!$ Allows evaluation of overlap for coordinate 2 first, then evaluation of complex
!!$ matrix element <varphi'|V01|varphi> is done 
  do np_orb = 1, np_orb_max
!!$ nuse defined in state_class will refer to one-electron target states
!!$ so can be used! 
     MOp1 = TargetStates2el%b(nstf)%nuse(np_orb)
     namp1 = get_nam(TargetStates%b(MOp1))
     Mp1 = get_ang_mom_proj(TargetStates%b(MOp1))
     ! INITIAL State MO COORDINATE 1
     do n_orb = 1, n_orb_max
        MO1 =  TargetStates2el%b(nsti)%nuse(n_orb)   
        nam1 = get_nam(TargetStates%b(MO1))
        M1 = get_ang_mom_proj(TargetStates%b(MO1))
        if (Mp1 - M1 /= mu) cycle

        overlap2 = 0d0
        
!!$ Looping over FINAL (p) Molecular State configorations 
        do nstp_con =  1, Npcon          
           use_nporb = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
           if ( MOp1 /= use_nporb ) cycle   ! 
           CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con) ! 2e config CI
           MOp2 = get_na(TargetStates2el%b(nstf),nstp_con,2) ! Index to 1e molecular orbital
           
!!$ Looping over INITIAL Molecular State configurations 
           do nst_con =  1, Ncon
              use_norb = get_na(TargetStates2el%b(nsti),nst_con,1)
              if ( MO1 /= use_norb ) cycle
              CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
              MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)

              ! overlap 2 matrix element <varphi'_2|varphi_2> 
              overlap2 = overlap2 + CI_stp*CI_st*ovlpst(MOp2,MO2)
              !if((nchf==3.or.nchf==4).and.nchi==1) write(100+nchf,*) 'MOp2, MO2, mp2, m2, ovlpst:',  MOp2, MO2, TargetStates%b(MOp2)%m, TargetStates%b(MO2)%m,  ovlpst(MOp2, MO2)
           end do ! nst_con
        end do ! nstp_con   
        !if((nchf==3.or.nchf==4).and.nchi==1) write(100+nchf,*) 'MOp1, MO1, Mp, M, overlap:', MOp1, MO1, Mp1, M1, overlap2
        !if((nchf==3.or.nchf==4).and.nchi==1) write(100+nchf,*) 
        if ( overlap2 == 0d0 ) cycle

        
!!$ COORDINATE 1 RADIAL INTEGRALS ttt(r_0) = 2.< |V01| >.< | >
        do ip1 = 1, namp1
!!$ Quantum numbers and functions for ALPHA
           indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
           tnp1 => bst_nr%b(indp1)   !           
           fp1 => fpointer(tnp1)     ! One electron functions
           Lap1 = get_ang_mom(tnp1)  ! Gets Angular momentum A.O.
           minfp = get_minf(tnp1)
           maxfp = get_maxf(tnp1)
           CI_1p = get_CI(TargetStates%b(MOp1),ip1)
           
           do i1 = 1, nam1
              ! Quantum numbers and functions for GAMMA
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

              lamMin1 = abs(Lap1 - La1)
              lamMax1 = Lap1 + La1
              lambda_min = max(lamMin0, lamMin1)
              lambda_max = min(lamMax0, lamMax1)
              lambda_max = min(data_in%ltmax,lambda_max)
              if ( lambda_min > lambda_max ) cycle

              ir1 = max(minf,minfp)
              ir2 = min(maxf,maxfp)
              fun1(:) = 0d0
!!$ Multiply the electron-electron potential by 2. Refer D. Fursa e-He 95 Paper


              !LHS: removed weight so we can use form_accurate below
!              fun1(ir1:ir2) = 2.0 * temp_CI * (f1(ir1:ir2) * fp1(ir1:ir2) * weight(ir1:ir2)) * overlap2
              fun1(ir1:ir2) = 2.0 * temp_CI * (f1(ir1:ir2) * fp1(ir1:ir2)) * overlap2
             
              if(.not.load) then
                do lambda = lambda_min, lambda_max

                   ang0 = dble((-1)**(mu))*Yint(dble(Lap0),dble(Mp0),dble(lambda),dble(-mu),dble(La0),dble(M0))
                   ang1 = Yint(dble(Lap1),dble(Mp1),dble(lambda),dble(mu),dble(La1),dble(M1))
                   if (ang0 == 0d0 .OR. ang1 == 0d0) cycle
                 
                   call form_accurate(lambda, fun1, ir1, ir2, nr, temp, irt1, irt2)
                   temp(irt1:irt2) = ang0*ang1 * temp(irt1:irt2)

                   ttt(irt1:irt2) = ttt(irt1:irt2) + temp(irt1:irt2)
                   icount = icount + 1
                    !if((nchf==3.or.nchf==4).and.nchi==1) write(1000+nchf,*) 'lambda, icount, temp, i1, i2:', lambda, icount,sum(temp(irt1:irt2)), irt1, irt2
                   ! Checking of angular terms non-zero
                   ir2_e_e = max(ir2_e_e,irt2)                 

                end do ! lambda
              else
                numint = numint + (ir2-ir1+1)
                if(ir2-ir1+1 < 0) error stop '1. IR2-IR1+1:'
                ir2_e_e = ir2
              endif !load

           end do ! i1
        end do ! ip1

     end do ! n_orb use_n
  end do ! np_orb use_np


  ! Projectile-Nuclear Potential
  if ( nsti == nstf) then

     !Liam added: lambda bounds here should only account for projectile l 
     lambda_min = LamMin0
     lambda_max = LamMax0

     if(data_in%good_parity) then
       lambda_step = 2
     else
       lambda_step = 1
     endif

     do lambda = lambda_min, lambda_max, lambda_step

        if ( lambda >  lamtop_vc) exit

        ir1 = minvnc(lambda)
        ir2 = maxvnc(lambda)      
        
        ! Checking of angular terms non-zero
        ir2_nuc = max(ir2_nuc,ir2)
    
        if ( lambda == 0) then
           ttt(1:ir2) = ttt(1:ir2) + vnc(1:ir2,lambda) + Zasym/gridr(1:ir2)

           ! Distorting potential -U_{0}
           if ( La0 <= data_in%Ldw ) then
              if(iBorn .eq. 0) then   
                 ttt(1:maxvdist) =  ttt(1:maxvdist) + Zproj * vdist(1:maxvdist)
              endif
           end if
        else
           ang0 = 0d0
           ang0 = Yint(dble(Lap0),dble(Mp0),dble(lambda),dble(mu),dble(La0),dble(M0))         
           if ( ang0 == 0d0 ) cycle
           ttt(ir1:ir2) = ttt(ir1:ir2) + ang0 * vnc(ir1:ir2,lambda)
        end if       

     end do
  end if


  if ( ir2_nuc == 0 .AND. ir2_e_e == 0 ) then
     return !  this ME is zero just due to angular momentum, it never got inside lambda loop
  end if

  call minmaxi(ttt,nr,ir1,ir2)
  ir_all = ir2
  if ( ir2 < nr ) then
     ttt(ir2+1:nr) = 0d0
  end if


  do ki=1,nqmi
     kqi = npk_nchi + ki - 1
     tn0 => chil%b(kqi)
     f0 => fpointer(tn0)
     minf = get_minf(tn0)
     maxf = get_maxf(tn0)
     
     do kf=1,nqmf
        kqf = npk_nchf + kf - 1
        tnp0 => chil%b(kqf)
        fp0 => fpointer(tnp0)
        minfp = get_minf(tnp0)
        maxfp = get_maxf(tnp0)
        
        ir2 = min(maxf,maxfp) ! nr as for continuum functions, but for closed channels maxf = 1
        ir1 = max(minf,minfp)

        if(iR2 < iR1) cycle

        if(.not.load) then
          ir2 = min(ir2,ir_all)
          fun0(ir1:ir2) = (f0(ir1:ir2) * fp0(ir1:ir2)) * weight(ir1:ir2) * ttt(ir1:ir2)  
          sum_rad = SUM( fun0(ir1:ir2) ) 
          vmatt(kf,ki) = vmatt(kf,ki) - Zproj * sum_rad * (2d0/pi)
        else
          numint = numint + (ir2-ir1+1)
        endif

     end do ! end kf loop
  end do ! end ki loop


end subroutine vdirect_2e

!Liam renamed this from vexch_2e_MOrep to vexch_2e since we always assume MO representation now
subroutine vexch_2e(Mtot,nchf,nchi,nstf,nsti,nqmf,nqmi,npk_nchf,npk_nchi,chil,Nmax_1e_orb,ncwaves,ortchil,flchil,nr,weight,gridr,&
    &Zproj,Etot,theta,Rd,rspin_tot,vmatt,numint)
  
! sum C^{np}_(np1,np2)  C^{n}_(n1,n2) < np0, np1, np2 : Mtot S |  | n0, n1, n2 : Mtot S >
! nstf = sum C^{np}_(np1,np2) | np1 np2 : m_f pi_f sp
!!$ Potential terms:  -2(E(1-theta)-H)P_{01}-E.I_0.theta, H=H0+H1+H2+V01+V02+V12+Z^2/R
  
  use sturmian_class
  use target_states
  use one_electron_func
  use channels
  use input_data
  use distpot
  use vnc_module
  use pol_2el
  use vmat_exch_module

  implicit none
  integer, intent(in):: Mtot
  integer, intent(in):: nchf, nchi        ! Channel number
  integer, intent(in):: nstf, nsti        ! State for particular channel
  integer, intent(in):: nqmf, nqmi        ! Number of k-grid points for particular channel
  integer, intent(in):: npk_nchf,npk_nchi ! On-shell k-grid point for particular channel
  integer, intent(in):: Nmax_1e_orb       ! Total number of target states one-electron orbitals      
!  real*8, dimension(Nmax_1e_orb,Nmax_1e_orb), intent(in):: ham1el_1elorb ! <np|H|n>  for the atomic orbitals
  integer, intent(in):: ncwaves           ! Number of on-, off- and bound projectile function ncwaves=nkgmax    
  real*8, dimension(Nmax_1e_orb,ncwaves), intent(in):: ortchil, flchil  ! ortchil=  <np|kLM>, flchil= <np|V-z.Zasy/r-U_dis|kLM>
  integer, intent(in):: nr                ! Max r points
  real*8, dimension(nr), intent(in):: weight, gridr
  real*8, intent(in):: Zproj        ! Z projectile charge
  real*8, intent(in):: Etot, theta        ! For non-uniqueness 
  real*8, intent(in):: Rd                 ! Internuclear distance      
  real*8, intent(in):: rspin_tot  ! Total Spin of the system. For H2 s_i=0 <=> S=1/2, s_i=1 <=> S=1/2,3/2
  real*8, dimension(nqmf,nqmi), intent(inout):: vmatt
  real*8, optional, intent(inout) :: numint

  type(basis_sturmian_nr):: chil
  
  integer:: ki, kf, kqi, kqf   ! Scroll and Index of projectile-waves
  real*8:: ki_en, kf_en ! Projectile-wave energy
  real*8, dimension(nr):: epot_i, fun1, fun2
! Configuration Loops
  integer::  nst_con,  nstp_con, Ncon, Npcon        ! 2e configuration loops
  integer::  n_orb, np_orb, n_orb_max, np_orb_max   ! loops nuse 
  integer:: use_norb, use_nporb                     ! index to 1e molecular orbital
  integer:: MO_0, MO1, MO2, MOp1, MOp2                ! index to 1e molecular orbital
  integer:: i0, i1, i2, ip1, ip2, nam0, nam1, nam2, namp1, namp2  ! Loop over 1e molecular orbitals atomic orbitals
  type(sturmian_nr), pointer:: tn0,tnp0,tn1, tn2, tnp1, tnp2  ! pointers to one-electron atomic orbitals
  integer:: ind0, ind1, ind2, indp2, indp1                      ! index to underlying one-electron atomic orbital
  real*8:: CI_stp, CI_st, CI_0, CI_1, CI_2, CI_1p, CI_2p, temp_CI ! 2e configuration and 1e molecular orbial coefficients
  real*8, pointer, dimension(:)::  fp0, f0, fp1, f1, fp2, f2  ! underlying atomic orbital radial functions
! Radial integrations
  integer:: minf0, maxf0, minfp0, maxfp0, minf1, maxf1, minfp1, maxfp1
  integer:: minf2, maxf2, minfp2, maxfp2
  integer:: ir1, ir2, jr1, jr2, or1, or2
  real*8:: pi, overlap2, temp_overlap2, temp, overlap1, theta_unique
  real*8, allocatable, dimension(:,:):: epot_lambda
  integer, allocatable, dimension(:):: epot_lambda_iminr, epot_lambda_imaxr 
  real*8, dimension(nr) :: temp_fun, temp2_fun, f0weight
  
  logical:: Logic_V12,  Logic_V02
  real*8:: Yint, COF6J
  real*8:: ang0, ang1, ang2                        ! Angular integrals for projectile and A.O.(1) 
  integer:: lambda, lambda_min, lambda_max, mu
  integer:: lambda_temp_min, lambda_temp_max, lambda_count
  integer:: kLap0, kLa1, kMp0, kM1                      ! Projectile Quantum Numbers (QN)  
  integer:: Lap1,Mp1,Lap2,Mp2,La0,M0,La2,M2,Mp_st,M_st  ! A.O. QN and Molecular State M 
  integer:: Spin, Spinp, tot_Spin
  integer:: La1, M1             ! QN and index for  projection operator term E.I.theta
  real*8 :: Z1, Z2
  logical :: load

  if(present(numint)) then
    load = .true.
    numint = 0
  else
    load = .false.
  endif

  Z1 = data_in%Z1
  Z2 = data_in%Z2

  pi = acos(-1d0)

  ! Projectile QN
  kLap0 = Lp_ch(nchf)  ! <kpLM(r_0)| 
  kMp0 = Mp_ch(nchf)   ! <kpLM(r_0)| 
  kLa1 = Lp_ch(nchi)   ! |kLM(r_1)>
  kM1 = Mp_ch(nchi)    ! |kLM(r_1)>
  
  ! Final State number nstf
  Npcon =  TargetStates2el%b(nstf)%nam          ! Number of 2e configurations.        
  Mp_st = NINT(TargetStates2el%b(nstf)%M )      ! Molecular state Angular Projection
  Spinp = NINT(TargetStates2el%b(nstf)%spin) ! 2e Molecular State Spin
  np_orb_max = TargetStates2el%b(nstf)%nusemax  ! Number of orbitals used to describe this state
  
  ! Initial State number nsti
  Ncon = TargetStates2el%b(nsti)%nam         
  M_st = NINT(TargetStates2el%b(nsti)%M)       
  Spin = NINT(TargetStates2el%b(nsti)%spin)
  n_orb_max = TargetStates2el%b(nsti)%nusemax    
  
  ! Conserved QN. Selection Rules. Parity should be taken care of in
  ! channels.f90
  if( Mp_st + kMp0  /=  M_st + kM1 ) return
!!$ Cp*C*<np2|n2>(- <kp0|V_0-z.Zasy/r-U_0|n0><np1|k1>
!!$ - <kp0|n0><np1|V_1-z.Zasy/r-U_1|k1> - <kp0,np1|V_01|n0,k1>)
!!$ FINAL (p) State Molecular Orbital (MO) COORDINATE 1
  do np_orb = 1, np_orb_max
     use_nporb = TargetStates2el%b(nstf)%nuse(np_orb)
!!$ INITIAL State MO COORDINATE 0 
     do n_orb = 1, n_orb_max
        use_norb =  TargetStates2el%b(nsti)%nuse(n_orb)   
        overlap2 = 0d0
!!$ Looping over FINAL (p) Molecular State configorations 
        do nstp_con =  1, Npcon          
           !MOp1 = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
           MOp1 = TargetStates2el%b(nstf)%na(nstp_con) ! Index 1e molecular orbital
           if ( MOp1 /= use_nporb ) cycle   ! 
           !CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con) ! 2e config CI
           CI_stp = TargetStates2el%b(nstf)%CI(nstp_con) ! 2e config CI
           !MOp2 = get_na(TargetStates2el%b(nstf),nstp_con,2) ! Index to 1e molecular orbital
           MOp2 = TargetStates2el%b(nstf)%nb(nstp_con) ! Index to 1e molecular orbital
           !Mp2 = get_ang_mom_proj(TargetStates%b(MOp2))      ! Angular momentum projection 1e molecular orbital
           Mp2 = TargetStates%b(MOp2)%m      ! Angular momentum projection 1e molecular orbital
           !namp2 = get_nam(TargetStates%b(MOp2))    
           namp2 = TargetStates%b(MOp2)%nam  
!!$ Looping over INITIAL Molecular State configurations 
           do nst_con =  1, Ncon
              !MO_0 = get_na(TargetStates2el%b(nsti),nst_con,1)
              MO_0 = TargetStates2el%b(nsti)%na(nst_con)
              if ( MO_0 /= use_norb ) cycle
              !CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
              CI_st = TargetStates2el%b(nsti)%CI(nst_con)
              !MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)
              MO2 = TargetStates2el%b(nsti)%nb(nst_con)
              !M2 = get_ang_mom_proj(TargetStates%b(MO2))
              M2 = TargetStates%b(MO2)%m
              !nam2 = get_nam(TargetStates%b(MO2))
              nam2 = TargetStates%b(MO2)%nam
              ! overlap 2 matrix element <varphi'_2|varphi_2> 
              if ( Mp2 /= M2 ) cycle  
!!$ FINAL STATE COORDINATE 2 for BETA
              do ip2 = 1, namp2
                 !indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
                 indp2 = TargetStates%b(MOp2)%na(ip2) ! Index to underlying atomic orbitals 
                 tnp2 => bst_nr%b(indp2)
                 !Lap2 = get_ang_mom(tnp2)
                 Lap2 = tnp2%l
                 !CI_2p = get_CI(TargetStates%b(MOp2),ip2)
                 CI_2p = TargetStates%b(MOp2)%CI(ip2)
!!$ INITIAL STATE COORDINATE 2 for DELTA                
                 do i2 = 1, nam2 
                    !ind2 = get_na(TargetStates%b(MO2),i2)  
                    ind2 = TargetStates%b(MO2)%na(i2) 
                    tn2 => bst_nr%b(ind2)                                         
                    !La2 = get_ang_mom(tn2)                            
                    La2 = tn2%l
                    !CI_2 = get_CI(TargetStates%b(MO2),i2)        
                    CI_2 = TargetStates%b(MO2)%CI(i2)
                    ! overlap 2 matrix element <varphi'_2|varphi_2> 
                    if ( Lap2 /= La2) cycle
                    temp_CI = CI_stp * CI_2p * CI_st * CI_2 
                    if ( temp_CI == 0d0 ) cycle 
                    temp_overlap2 = bst_nr%ortint(indp2,ind2)
                    overlap2 = overlap2 + temp_CI * temp_overlap2
                 end do ! i2        
              end do ! ip2   
           end do ! nst_con - 2e config
        end do ! nstp_con
!!$     Here do one electron matrix elements and V_01 term:
!!$     Cp*C*<np2|n2>(-<kp0|V_0-z.Zasy/r-U_0|n0><np1|k1>
!!$     - <kp0|n0><np1|V_1-z.Zasy/r-U_1|k1> - <kp0,np1|V_01|n0,k1>)
        if ( overlap2 == 0d0 ) cycle 

        MOp1 = use_nporb ! get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
        !Mp1 = get_ang_mom_proj(TargetStates%b(MOp1))
        Mp1 = TargetStates%b(MOp1)%m
        !namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
        namp1 = TargetStates%b(MOp1)%nam             ! Number of atomic orbitals that make up molecular orbital 
        MO_0 = use_norb ! get_na(TargetStates2el%b(nsti),nst_con,1)
        !M0 = get_ang_mom_proj(TargetStates%b(MO_0))
        M0 = TargetStates%b(MO_0)%m
        !nam0 = get_nam(TargetStates%b(MO_0))
        nam0 = TargetStates%b(MO_0)%nam
!!$ Quantum numbers and functions for ALPHA             
        do ip1 = 1, namp1
           !indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
           indp1 = TargetStates%b(MOp1)%na(ip1)  ! Index to underlying atomic orbitals 
           tnp1 => bst_nr%b(indp1)   !           
           fp1 => fpointer(tnp1)     ! One electron functions
           !Lap1 = get_ang_mom(tnp1)  ! Gets Angular momentum A.O.
           Lap1 = tnp1%l  ! Gets Angular momentum A.O.
           !minfp1 = get_minf(tnp1)
           minfp1 = tnp1%minf
           !maxfp1 = get_maxf(tnp1)
           maxfp1 = tnp1%maxf
           !CI_1p = get_CI(TargetStates%b(MOp1),ip1)
           CI_1p = TargetStates%b(MOp1)%CI(ip1)
!!$ Calculating -Cp*C*<np2|n2><kp0,np1|V_01|n0,k1>
!!!$ Inital COORDINATE 1 |kLM(r_1)>
           do ki = 1, nqmi
              kqi = npk_nchi + ki - 1
              tn1 => chil%b(kqi)
              !ki_en = get_energy(tn1)
              ki_en = tn1%en
              f1 => fpointer(tn1)
              !minf1 = get_minf(tn1)
              minf1 = tn1%minf
              !maxf1 = get_maxf(tn1)
              maxf1 = tn1%maxf
              or1 = max(minf1,minfp1)
              or2 = min(maxf1,maxfp1)
              
              !LHS: removed weight so we can use form_accurate below
!              fun1(or1:or2) = ( fp1(or1:or2) * f1(or1:or2) * weight(or1:or2) )
              fun1(or1:or2) = ( fp1(or1:or2) * f1(or1:or2) )
!!$ Quantum numbers and functions for GAMMA                  
              do i0 = 1, nam0
                 !ind0 = get_na(TargetStates%b(MO_0),i0) ! Index to underlying atomic orbitals  
                 ind0 = TargetStates%b(MO_0)%na(i0) ! Index to underlying atomic orbitals  
                 tn0 => bst_nr%b(ind0)
                 f0 => fpointer(tn0)
                 !La0 = get_ang_mom(tn0)
                 La0 = tn0%l
                 !minf0 = get_minf(tn0)
                 minf0 = tn0%minf
                 !maxf0 = get_maxf(tn0)
                 maxf0 = tn0%maxf
                 !CI_0 = get_CI(TargetStates%b(MO_0),i0)
                 CI_0 = TargetStates%b(MO_0)%CI(i0)
                 temp_CI = CI_1p * CI_0
                 if ( temp_CI == 0d0 ) cycle
                 lambda_min = max(ABS(kLap0 - La0),ABS(Lap1 - kLa1))
                 lambda_max = min(kLap0 + La0, Lap1 + kLa1)
                 lambda_max = min(data_in%ltmax,lambda_max)

                 f0weight(minf0:maxf0) = f0(minf0:maxf0) * weight(minf0:maxf0)

                 ! mu should be the same for both integrals in coordinate 0 and 1
                 mu = Mp1 - kM1
                 do lambda = lambda_min, lambda_max
                    ang0 = 0d0
                    ang1 = 0d0
                    ang0 = dble((-1)**(mu))*Yint(dble(kLap0),dble(kMp0),dble(lambda),dble(-mu),dble(La0),dble(M0))
                    ang1 = Yint(dble(Lap1),dble(Mp1),dble(lambda),dble(mu),dble(kLa1),dble(kM1))
                    if ( ang0 * ang1 == 0d0 ) cycle

                    
                    if(.not.load) then
                      ! Integrate over coordinate 1 
                      call form_accurate(lambda,fun1,or1,or2,nr,epot_i,jr1,jr2)
                      epot_i(jr1:jr2) = epot_i(jr1:jr2) * f0weight(jr1:jr2) !multiply by f0 and weight here rather than in the sum below
                    else
                      numint = numint + or2-or1+1
                      jr1 = or1
                      jr2 = or2
                    endif

!!$ COORDINATE 0 <kpLM(r_0)|
                    do kf = 1, nqmf
                       kqf = npk_nchf + kf - 1
                       tnp0 => chil%b(kqf)
                       fp0 => fpointer(tnp0)
                       !minfp0 = get_minf(tnp0)
                       minfp0 = tnp0%minf
                       !maxfp0 = get_maxf(tnp0)
                       maxfp0 = tnp0%maxf
                       ir1 = max(minf0,minfp0)
                       ir2 = min(maxf0,maxfp0)
                       ir1 = max(jr1,ir1)
                       ir2 = min(jr2,ir2)
                       if(ir2<ir1) cycle
                       ! Integrate over coordinate 0 
                       if(.not.load) then
                         temp = SUM(fp0(ir1:ir2) * epot_i(ir1:ir2)) 
                         ! overlap2 has CI coefficients
                         temp = temp * ang0 * ang1 * overlap2
                         vmatt(kf,ki) = vmatt(kf,ki) - temp_CI * temp
                       else
                         numint = numint + ir2-ir1+1
                       endif
                    end do ! end kf loop
                 end do ! lambda
!!$ Calculating Cp*C*<np2|n2>(-<kp0|V_0-z.Zasy/r-U_0|n0><np1|k1> - <kp0|n0><np1|V_1-z.Zasy/r-U_1|k1> ) 
                 do kf = 1, nqmf
                    kqf = npk_nchf + kf - 1
                    tnp0 => chil%b(kqf)
                    vmatt(kf,ki) = vmatt(kf,ki) - temp_CI * flchil(ind0,kqf) * ortchil(indp1,kqi) * overlap2  
                    vmatt(kf,ki) = vmatt(kf,ki) - temp_CI * ortchil(ind0,kqf) * flchil(indp1,kqi) * overlap2
                 end do  ! end kf loop

              end do ! i0
           end do ! ki loop
        end do ! ip1
     end do ! n_orb
  end do ! np_orb


!!$ Taking care of non-uniqueness
!!$ Cp*C*<np2|n2><kp0|n0><np1|k1>(E(1-theta)-e_k'-e_k-Z^2/R) 
!!$ Looping over FINAL (p) Molecular State configorations 
  do nstp_con =  1, Npcon
     !CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con) ! 2e config CI
     CI_stp = TargetStates2el%b(nstf)%CI(nstp_con) ! 2e config CI
     !MOp2 = get_na(TargetStates2el%b(nstf),nstp_con,2) ! Index to 1e molecular orbital
     MOp2 = TargetStates2el%b(nstf)%nb(nstp_con) ! Index to 1e molecular orbital
     !Mp2 = get_ang_mom_proj(TargetStates%b(MOp2))      ! Angular momentum projection 1e molecular orbital
     Mp2 = TargetStates%b(MOp2)%m      ! Angular momentum projection 1e molecular orbital
     !namp2 = get_nam(TargetStates%b(MOp2))
     namp2 = TargetStates%b(MOp2)%nam
     !MOp1 = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
     MOp1 = TargetStates2el%b(nstf)%na(nstp_con) ! Index 1e molecular orbital
     !Mp1 = get_ang_mom_proj(TargetStates%b(MOp1))
     Mp1 = TargetStates%b(MOp1)%m
     !namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
     namp1 = TargetStates%b(MOp1)%nam             ! Number of atomic orbitals that make up molecular orbital 
!!$ Looping over INITIAL Molecular State configurations 
     do nst_con =  1, Ncon
        !CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
        CI_st = TargetStates2el%b(nsti)%CI(nst_con)
        !MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)
        MO2 = TargetStates2el%b(nsti)%nb(nst_con)
        !M2 = get_ang_mom_proj(TargetStates%b(MO2))
        M2 = TargetStates%b(MO2)%m
        !nam2 = get_nam(TargetStates%b(MO2))
        nam2 = TargetStates%b(MO2)%nam
        !MO_0 = get_na(TargetStates2el%b(nsti),nst_con,1)
        MO_0 = TargetStates2el%b(nsti)%na(nst_con)
        !M0 = get_ang_mom_proj(TargetStates%b(MO_0))
        M0 = TargetStates%b(MO_0)%m
        !nam0 = get_nam(TargetStates%b(MO_0))
        nam0 = TargetStates%b(MO_0)%nam
        overlap2 = 0d0
!!$ Solve for non-uniqueness if core orbital
        theta_unique = 1d0
        if(data_in%non_uniq) then
          if (is_core_MO(MO2) == 2) theta_unique = (1d0 - theta)
          if (is_core_MO(MO2) == 1 .AND. is_core_MO(MOp1) /= 0 ) theta_unique = (1d0 - theta)
        endif
!!$      write(*,'(4I4,F10.5)') MOp1,MOp2,MO_0,MO2,theta_unique
!!$ overlap 2 matrix element <varphi'_2|varphi_2> 
        if ( Mp2 /= M2 ) cycle
!!$ FINAL STATE COORDINATE 2 for BETA
        do ip2 = 1, namp2
           !indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
           indp2 = TargetStates%b(MOp2)%na(ip2) ! Index to underlying atomic orbitals 
           tnp2 => bst_nr%b(indp2)
           !Lap2 = get_ang_mom(tnp2)
           Lap2 = tnp2%l
           !CI_2p = get_CI(TargetStates%b(MOp2),ip2)
           CI_2p = TargetStates%b(MOp2)%CI(ip2)
!!$ INITIAL STATE COORDINATE 2 for DELTA                
           do i2 = 1, nam2
              !ind2 = get_na(TargetStates%b(MO2),i2)
              ind2 = TargetStates%b(MO2)%na(i2)
              tn2 => bst_nr%b(ind2)
              !La2 = get_ang_mom(tn2)
              La2 = tn2%l
              !CI_2 = get_CI(TargetStates%b(MO2),i2)
              CI_2 = TargetStates%b(MO2)%CI(i2)
              ! overlap 2 matrix element <varphi'_2|varphi_2> 
              if ( Lap2 /= La2) cycle
              temp_CI = CI_stp * CI_2p * CI_st * CI_2
              if ( temp_CI == 0d0 ) cycle
              temp_overlap2 = bst_nr%ortint(indp2,ind2)
              overlap2 = overlap2 + temp_CI * temp_overlap2
           end do ! i2        
        end do ! ip2   
        if ( overlap2 == 0d0 ) cycle
        !MOp1 = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
        MOp1 = TargetStates2el%b(nstf)%na(nstp_con) ! Index 1e molecular orbital
        !Mp1 = get_ang_mom_proj(TargetStates%b(MOp1))
        Mp1 = TargetStates%b(MOp1)%m
        !namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
        namp1 = TargetStates%b(MOp1)%nam             ! Number of atomic orbitals that make up molecular orbital 

        !MO_0 = get_na(TargetStates2el%b(nsti),nst_con,1)
        MO_0 = TargetStates2el%b(nsti)%na(nst_con)
        !M0 = get_ang_mom_proj(TargetStates%b(MO_0))
        M0 = TargetStates%b(MO_0)%m
        !nam0 = get_nam(TargetStates%b(MO_0))
        nam0 = TargetStates%b(MO_0)%nam
!!$ Calculating Cp*C*<np2|n2>(<kp0|n0><np1|k1>(E(1-theta)-e_k'-e_k-Z^2/R) 
        if ( kM1 /= Mp1 .OR.  kMp0 /= M0 ) cycle
!!$ Quantum numbers and functions for ALPHA             
        do ip1 = 1, namp1
           !indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
           indp1 = TargetStates%b(MOp1)%na(ip1)  ! Index to underlying atomic orbitals 
           tnp1 => bst_nr%b(indp1)   !           
           !Lap1 = get_ang_mom(tnp1)  ! Gets Angular momentum A.O.
           Lap1 = tnp1%l  ! Gets Angular momentum A.O.
           !CI_1p = get_CI(TargetStates%b(MOp1),ip1)
           CI_1p = TargetStates%b(MOp1)%CI(ip1)
           if (  kLa1 /= Lap1 ) cycle
!!$ Quantum numbers and functions for GAMMA                  
           do i0 = 1, nam0
              !ind0 = get_na(TargetStates%b(MO_0),i0) ! Index to underlying atomic orbitals  
              ind0 = TargetStates%b(MO_0)%na(i0) ! Index to underlying atomic orbitals  
              tn0 => bst_nr%b(ind0)                 
              !La0 = get_ang_mom(tn0)
              La0 = tn0%l
              !CI_0 = get_CI(TargetStates%b(MO_0),i0)
              CI_0 = TargetStates%b(MO_0)%CI(i0)
              if ( kLap0 /=  La0 ) cycle
              temp_CI = CI_1p * CI_0
              if ( temp_CI == 0d0 ) cycle
!!!$ Initial COORDINATE 1 |kLM(r_1)>
              do ki = 1, nqmi
                 kqi = npk_nchi + ki - 1
                 tn1 => chil%b(kqi)
                 !ki_en = get_energy(tn1)
                 ki_en = tn1%en
!!$ Calculating Cp*C*<np2|n2>(<kp0|n0><np1|k1>(E(1-theta)-e_k'-e_k-Z^2/R) 
                 do kf = 1, nqmf
                    kqf = npk_nchf + kf - 1
                    tnp0 => chil%b(kqf)
                    !kf_en = get_energy(tnp0)
                    kf_en = tnp0%en
                    temp = ortchil(ind0,kqf) * ortchil(indp1,kqi) * (Etot*(theta_unique)-kf_en-ki_en)

                    vmatt(kf,ki) = vmatt(kf,ki) + temp_CI * overlap2 * temp
                    if ( Rd /= 0d0 ) then
                       vmatt(kf,ki) = vmatt(kf,ki) - temp_CI * Z1*Z2/Rd * ortchil(ind0,kqf) * ortchil(indp1,kqi) * overlap2

                    end if
                 end do  ! end kf loop
              end do ! ki loop
           end do ! i0
        end do ! ip1
     end do ! nst_con - 2e config
  end do ! nstp_con


!!$ Calculating -Cp*C*<np1|k1><kp0,np2|V_02|n0,n2>
!!$ Looping over FINAL Molecular State orbitals
!!$ Looping over FINAL (p) Molecular State configorations 
  do nstp_con =  1, Npcon
     !CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con) ! 2e config CI
     CI_stp = TargetStates2el%b(nstf)%CI(nstp_con) ! 2e config CI
     !MOp1 = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
     MOp1 = TargetStates2el%b(nstf)%na(nstp_con) ! Index 1e molecular orbital
     !Mp1 = get_ang_mom_proj(TargetStates%b(MOp1))      ! Angular momentum projection 1e molecular orbital
     Mp1 = TargetStates%b(MOp1)%m      ! Angular momentum projection 1e molecular orbital
     !MOp2 = get_na(TargetStates2el%b(nstf),nstp_con,2) ! Index to 1e molecular orbital
     MOp2 = TargetStates2el%b(nstf)%nb(nstp_con) ! Index to 1e molecular orbital
     !Mp2 = get_ang_mom_proj(TargetStates%b(MOp2))      ! Angular momentum projection 1e molecular orbital
     Mp2 = TargetStates%b(MOp2)%m      ! Angular momentum projection 1e molecular orbital
     !namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
     namp1 = TargetStates%b(MOp1)%nam             ! Number of atomic orbitals that make up molecular orbital 
     namp2 = TargetStates%b(MOp2)%nam
     ! Check that coordinate 1 overlaps are non-zero
     Logic_V02 = .FALSE.
     do ki = 1, nqmi
        kqi = npk_nchi + ki - 1
        do ip1 = 1, namp1
           !indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
           indp1 = TargetStates%b(MOp1)%na(ip1)  ! Index to underlying atomic orbitals 
           if ( ortchil(indp1,kqi) /= 0d0 ) then
              Logic_V02=.TRUE.
              exit
           end if
        end do ! ip1
        if (Logic_V02) exit
     end do ! ki
     if ( .NOT. Logic_V02 ) cycle
!!$ Looping over INITIAL Molecular State configurations 
     do nst_con =  1, Ncon
        !CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
        CI_st = TargetStates2el%b(nsti)%CI(nst_con)
        !MO_0 = get_na(TargetStates2el%b(nsti),nst_con,1)
        MO_0 = TargetStates2el%b(nsti)%na(nst_con)
        !M0 = get_ang_mom_proj(TargetStates%b(MO_0))
        M0 = TargetStates%b(MO_0)%m
        !MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)
        MO2 = TargetStates2el%b(nsti)%nb(nst_con)
        !M2 = get_ang_mom_proj(TargetStates%b(MO2))
        M2 = TargetStates%b(MO2)%m
        !nam0 = get_nam(TargetStates%b(MO_0))
        nam0 = TargetStates%b(MO_0)%nam
        !nam2 = get_nam(TargetStates%b(MO2))
        nam2 = TargetStates%b(MO2)%nam
        mu = Mp2 - M2
!!$ FINAL STATE COORDINATE 2 for BETA      
        do ip2 = 1, namp2
           !indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
           indp2 = TargetStates%b(MOp2)%na(ip2) ! Index to underlying atomic orbitals 
           tnp2 => bst_nr%b(indp2)
           fp2 => fpointer(tnp2)
           !Lap2 = get_ang_mom(tnp2)
           Lap2 = tnp2%l
           !CI_2p = get_CI(TargetStates%b(MOp2),ip2)
           CI_2p = TargetStates%b(MOp2)%CI(ip2)
           !minfp2 = get_minf(tnp2)
           minfp2 = tnp2%minf
           !maxfp2 = get_maxf(tnp2) 
           maxfp2 = tnp2%maxf
!!$ INITIAL STATE COORDINATE 2 for DELTA          
           do i2 = 1, nam2
              !ind2 = get_na(TargetStates%b(MO2),i2)
              ind2 = TargetStates%b(MO2)%na(i2)
              tn2 => bst_nr%b(ind2)
              f2 => fpointer(tn2)
              !La2 = get_ang_mom(tn2)
              La2 = tn2%l
              !CI_2 = get_CI(TargetStates%b(MO2),i2)
              CI_2 = TargetStates%b(MO2)%CI(i2)
              temp_CI = CI_stp * CI_2p * CI_st * CI_2
              if ( temp_CI == 0d0 ) cycle
              !minf2 = get_minf(tn2)
              minf2 = tn2%minf
              !maxf2 = get_maxf(tn2)
              maxf2 = tn2%maxf
              
              or1 = max(minfp2,minf2)     
              or2 = min(maxfp2,maxf2)
              
              !LHS: removed weight so we can use form_accurate below
              !fun2(or1:or2) = fp2(or1:or2) * f2(or1:or2) * weight(or1:or2)
              fun2(or1:or2) = fp2(or1:or2) * f2(or1:or2)
!!$ INITIAL STATE COORDINATE 0 for GAMMA
!!$ Determine bounds of lambda and calculate form once.
              lambda_count = 0
              do i0 = 1, nam0
                 !ind0 = get_na(TargetStates%b(MO_0),i0) ! Index to underlying atomic orbitals  
                 ind0 = TargetStates%b(MO_0)%na(i0) ! Index to underlying atomic orbitals  
                 tn0 => bst_nr%b(ind0)
                 !La0 = get_ang_mom(tn0)
                 La0 = tn0%l
                 !CI_0 = get_CI(TargetStates%b(MO_0),i0)
                 CI_0 = TargetStates%b(MO_0)%CI(i0)
                 temp_CI = CI_stp * CI_2p * CI_st * CI_0 * CI_2
                 if ( temp_CI == 0d0 ) cycle
                 lambda_temp_min = max(ABS(kLap0 - La0), ABS(Lap2 - La2))
                 lambda_temp_max = min(kLap0 + La0, Lap2 + La2)
                 if ( lambda_count == 0 ) then
                    lambda_max = lambda_temp_max
                    lambda_min = lambda_temp_min
                    lambda_count = 1
                 end if
                 if ( lambda_temp_max > lambda_max ) lambda_max = lambda_temp_max
                 if ( lambda_temp_min < lambda_min ) lambda_min = lambda_temp_min
              end do ! i0
              if ( allocated(epot_lambda) ) deallocate(epot_lambda)
              if ( allocated(epot_lambda_imaxr)) deallocate(epot_lambda_imaxr,epot_lambda_iminr)
              allocate(epot_lambda(lambda_min:lambda_max,nr),epot_lambda_imaxr(lambda_min:lambda_max),epot_lambda_iminr(lambda_min:lambda_max) )
              lambda_max = min(data_in%ltmax,lambda_max)
              epot_lambda_iminr = nr
              epot_lambda_imaxr = 1
              do lambda = lambda_min, lambda_max
!                 ang0 = 0d0
                 ang2 = 0d0
!!$ No defintion of La0 here 
!                 ang0 = dble((-1)**(mu))*Yint(dble(kLap0),dble(kMp0),dble(lambda),dble(-mu),dble(La0),dble(M0))
                 ang2 = Yint(dble(Lap2),dble(Mp2),dble(lambda),dble(mu),dble(La2),dble(M2))
                 if (  ang2 == 0d0 ) cycle
                 call form_accurate(lambda,fun2,or1,or2,nr,epot_i,jr1,jr2)
                 epot_lambda(lambda,jr1:jr2) = epot_i(jr1:jr2) 
                 epot_lambda_iminr(lambda) = jr1
                 epot_lambda_imaxr(lambda) = jr2

              end do ! lambda 

!!$  INITIAL STATE COORDINATE 0 for GAMMA
              do i0 = 1, nam0
                 !ind0 = get_na(TargetStates%b(MO_0),i0) ! Index to underlying atomic orbitals  
                 ind0 = TargetStates%b(MO_0)%na(i0) ! Index to underlying atomic orbitals  
                 tn0 => bst_nr%b(ind0)
                 f0 => fpointer(tn0)
                 !La0 = get_ang_mom(tn0)
                 La0 = tn0%l
                 !minf0 = get_minf(tn0)
                 minf0 = tn0%minf
                 !maxf0 = get_maxf(tn0)
                 maxf0 = tn0%maxf
                 !CI_0 = get_CI(TargetStates%b(MO_0),i0)
                 CI_0 = TargetStates%b(MO_0)%CI(i0)
                 temp_CI = CI_stp * CI_2p * CI_st * CI_0 * CI_2
                 if ( temp_CI == 0d0 ) cycle
                 lambda_min = max(ABS(kLap0 - La0),ABS(Lap2 - La2))
                 lambda_max = min(kLap0 + La0, Lap2 + La2)
!!$ Could include and if statement here to only calculate this part if it is required for this i0
!!$ May not be required if ang0 * ang2 == 0d0 for all lambda...
                 temp_fun(minf0:maxf0) = f0(minf0:maxf0) * weight(minf0:maxf0)
                 do lambda = lambda_min, lambda_max
                    ang0 = 0d0
                    ang2 = 0d0
                    ang0 = dble((-1)**(mu))*Yint(dble(kLap0),dble(kMp0),dble(lambda),dble(-mu),dble(La0),dble(M0))
                    ang2 = Yint(dble(Lap2),dble(Mp2),dble(lambda),dble(mu),dble(La2),dble(M2))
                    if ( ang0 * ang2 == 0d0 ) cycle
                    jr1 = max(minf0,epot_lambda_iminr(lambda))
                    jr2 = min(maxf0,epot_lambda_imaxr(lambda))
                    if(jr2<0) then
                      print*, 'imaxr, maxf0, jr2:', epot_lambda_imaxr(lambda), maxf0, jr2
                      error stop
                    endif
                    
                    temp2_fun(jr1:jr2) = temp_fun(jr1:jr2) * epot_lambda(lambda,jr1:jr2) 
!!$ COORDINATE 0 <kpLM(r_0)|
                    do kf = 1, nqmf
                       kqf = npk_nchf + kf - 1
                       tnp0 => chil%b(kqf)
                       fp0 => fpointer(tnp0)
                       !minfp0 = get_minf(tnp0)
                       minfp0 = tnp0%minf
                       !maxfp0 = get_maxf(tnp0)
                       maxfp0 = tnp0%maxf
                       ir1 = max(jr1,minfp0)
                       ir2 = min(jr2,maxfp0)
                       if(ir2<ir1) cycle
!!$ Integrate over coordinate 0, temp = <kf|(CI_stp.CI_2p.CI_st.CI_0.CI_2.<varphi'_2|V_02|varphi_2>)|varphi_0> 
!!$ temp2_fun0 = epot_lambda(lambda,jr1:jr2) * f0(minf0:maxf0) * weight(minf0:maxf0)
!!$                   temp = SUM(fp0(ir1:ir2) * epot_lambda(lambda,ir1:ir2) * f0(ir1:ir2) * weight(ir1:ir2))
                       temp = SUM(fp0(ir1:ir2) * temp2_fun(ir1:ir2))
                       temp = temp_CI * temp * ang0 * ang2
!!$ Overlap coordinate 1
                       do ki = 1, nqmi
                          kqi = npk_nchi + ki - 1
                          do ip1 = 1, namp1
                             !indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
                             indp1 = TargetStates%b(Mop1)%na(ip1)
                             !CI_1p = get_CI(TargetStates%b(MOp1),ip1)
                             CI_1p = TargetStates%b(Mop1)%CI(ip1)
!!$ CI_p.<varphi'_1|ki>.(<kf|(CI_stp.CI_2p.CI_st.CI_0.CI_2.<varphi'_2|V_02|varphi_2>)|varphi_0>)                        
                             vmatt(kf,ki)= vmatt(kf,ki) - CI_1p * temp * ortchil(indp1,kqi)
                          end do ! ip1
                       end do ! ki cord 1
                    end do ! kf cord 0
                 end do ! lambda 
              end do ! i0 
           end do ! i2 
        end do !ip2
     end do ! nst_con
  end do ! nstp_con


!!$ Calculating -Cp*C*<kp0|n0><np1,np2|V_12|k1,n2>
!!$ Looping over FINAL Molecular State orbitals
  do nstp_con =  1, Npcon
     !CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con) ! 2e config CI
     CI_stp = TargetStates2el%b(nstf)%CI(nstp_con) ! 2e config CI
     !MOp1 = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
     MOp1 = TargetStates2el%b(nstf)%na(nstp_con) ! Index 1e molecular orbital
     !Mp1 = get_ang_mom_proj(TargetStates%b(MOp1))      ! Angular momentum projection 1e molecular orbital
     Mp1 = TargetStates%b(MOp1)%m      ! Angular momentum projection 1e molecular orbital
     !MOp2 = get_na(TargetStates2el%b(nstf),nstp_con,2) ! Index to 1e molecular orbital
     MOp2 = TargetStates2el%b(nstf)%nb(nstp_con) ! Index to 1e molecular orbital
     !Mp2 = get_ang_mom_proj(TargetStates%b(MOp2))      ! Angular momentum projection 1e molecular orbital
     Mp2 = TargetStates%b(MOp2)%m      ! Angular momentum projection 1e molecular orbital
     !namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
     namp1 = TargetStates%b(MOp1)%nam             ! Number of atomic orbitals that make up molecular orbital 
     !namp2 = get_nam(TargetStates%b(MOp2))
     namp2 = TargetStates%b(MOp2)%nam
!!$ Looping over INITIAL Molecular State configurations 
     do nst_con =  1, Ncon
        !CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
        CI_st = TargetStates2el%b(nsti)%CI(nst_con)
        !MO_0 = get_na(TargetStates2el%b(nsti),nst_con,1)
        MO_0 = TargetStates2el%b(nsti)%na(nst_con)
        !M0 = get_ang_mom_proj(TargetStates%b(MO_0))
        M0 = TargetStates%b(MO_0)%m
        !MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)
        MO2 = TargetStates2el%b(nsti)%nb(nst_con)
        !M2 = get_ang_mom_proj(TargetStates%b(MO2))
        M2 = TargetStates%b(MO2)%m
        !nam0 = get_nam(TargetStates%b(MO_0))
        nam0 = TargetStates%b(MO_0)%nam
        !nam2 = get_nam(TargetStates%b(MO2))
        nam2 = TargetStates%b(MO2)%nam
        mu = Mp2 - M2
!!$ Check that coordinate 0 overlaps are non-zero
        Logic_V12 = .FALSE.
        do kf = 1, nqmf
           kqf = npk_nchf + kf - 1
           do i0 = 1, nam0
              !ind0 = get_na(TargetStates%b(MO_0),i0) ! Index to underlying atomic orbitals 
              ind0 = TargetStates%b(MO_0)%na(i0) ! Index to underlying atomic orbitals 
              if ( ortchil(ind0,kqf) /= 0d0 ) then
                 Logic_V12=.TRUE.
                 exit
              end if
           end do ! i0
           if ( Logic_V12 ) exit
        end do ! kf
        if ( .NOT. Logic_V12 ) cycle
!!$ FINAL STATE COORDINATE 2 BETA
        do ip2 = 1, namp2
           !indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
           indp2 = TargetStates%b(MOp2)%na(ip2) ! Index to underlying atomic orbitals 
           tnp2 => bst_nr%b(indp2)
           fp2 => fpointer(tnp2)
           !Lap2 = get_ang_mom(tnp2)
           Lap2 = tnp2%l
           !CI_2p = get_CI(TargetStates%b(MOp2),ip2)
           CI_2p = TargetStates%b(MOp2)%CI(ip2)
           !minfp2 = get_minf(tnp2)
           minfp2 = tnp2%minf
           !maxfp2 = get_maxf(tnp2)
           maxfp2 = tnp2%maxf
!!$ INITIAL STATE COORDINATE 2 DELTA
           do i2 = 1, nam2
              !ind2 = get_na(TargetStates%b(MO2),i2)
              ind2 = TargetStates%b(MO2)%na(i2)
              tn2 => bst_nr%b(ind2)
              f2 => fpointer(tn2)
              !La2 = get_ang_mom(tn2)
              La2 = tn2%l
              !minf2 = get_minf(tn2)
              minf2 = tn2%minf
              !maxf2 = get_maxf(tn2)
              maxf2 = tn2%maxf
              !CI_2 = get_CI(TargetStates%b(MO2),i2)
              CI_2 = TargetStates%b(MO2)%CI(i2)
              temp_CI = CI_stp * CI_2p * CI_st * CI_2
              if ( temp_CI == 0d0 ) cycle
              or1 = max(minfp2,minf2)
              or2 = min(maxfp2,maxf2)
             
              !LHS: removed weight so we can use form_accurate below
!              fun2(or1:or2) = fp2(or1:or2) * f2(or1:or2) * weight(or1:or2)
              fun2(or1:or2) = fp2(or1:or2) * f2(or1:or2)
!!$ FINAL STATE COORDINATE 1 ALPHA
!!$ Determine bounds of lambda and calculate form once.
              lambda_count = 0
              do ip1 = 1, namp1
                 !indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
                 indp1 = TargetStates%b(MOp1)%na(ip1)  ! Index to underlying atomic orbitals 
                 tnp1 => bst_nr%b(indp1)   !           
                 !Lap1 = get_ang_mom(tnp1)  ! Gets Angular momentum A.O.
                 Lap1 = tnp1%l  ! Gets Angular momentum A.O.
                 !CI_1p = get_CI(TargetStates%b(MOp1),ip1)
                 CI_1p = TargetStates%b(MOp1)%CI(ip1)
                 temp_CI = CI_stp * CI_1p * CI_2p * CI_st * CI_2
                 if ( temp_CI == 0d0 ) cycle
                 lambda_temp_min = max(ABS(Lap1 - kLa1), ABS(Lap2 - La2))
                 lambda_temp_max = min(Lap1 + kLa1, Lap2 + La2)
                 if ( lambda_count == 0 ) then
                    lambda_max = lambda_temp_max
                    lambda_min = lambda_temp_min
                    lambda_count = 1
                 end if
                 if ( lambda_temp_max > lambda_max ) lambda_max = lambda_temp_max
                 if ( lambda_temp_min < lambda_min ) lambda_min = lambda_temp_min
              end do ! ip1
              if ( allocated(epot_lambda) ) deallocate(epot_lambda)
              if ( allocated(epot_lambda_imaxr)) deallocate(epot_lambda_imaxr,epot_lambda_iminr)
              allocate(epot_lambda(lambda_min:lambda_max,nr),epot_lambda_imaxr(lambda_min:lambda_max),epot_lambda_iminr(lambda_min:lambda_max) )
              lambda_max = min(data_in%ltmax,lambda_max)
              do lambda = lambda_min, lambda_max
                 ang2 = 0d0
                 ang2 = Yint(dble(Lap2),dble(Mp2),dble(lambda),dble(mu),dble(La2),dble(M2))
                 if ( ang2 == 0d0 ) cycle
                 call form_accurate(lambda,fun2,or1,or2,nr,epot_i,jr1,jr2)
                 epot_lambda(lambda,jr1:jr2) = epot_i(jr1:jr2)
                 epot_lambda_iminr(lambda) = jr1
                 epot_lambda_imaxr(lambda) = jr2
              end do ! lambda 
!!$ FINAL STATE COORDINATE 1 ALPHA
              do ip1 = 1, namp1
                 !indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
                 indp1 = TargetStates%b(MOp1)%na(ip1)  ! Index to underlying atomic orbitals 
                 tnp1 => bst_nr%b(indp1)   !           
                 fp1 => fpointer(tnp1)     ! One electron functions
                 !Lap1 = get_ang_mom(tnp1)  ! Gets Angular momentum A.O.
                 Lap1 = tnp1%l  ! Gets Angular momentum A.O.
                 !minfp1 = get_minf(tnp1)
                 minfp1 = tnp1%minf
                 !maxfp1 = get_maxf(tnp1)
                 maxfp1 = tnp1%maxf
                 !CI_1p = get_CI(TargetStates%b(MOp1),ip1)
                 CI_1p = TargetStates%b(MOp1)%CI(ip1)
                 temp_CI = CI_stp * CI_1p * CI_2p * CI_st * CI_2
                 if ( temp_CI == 0d0 ) cycle
                 lambda_min = max(ABS(Lap1 - kLa1), ABS(Lap2 - La2))
                 lambda_max = min(Lap1 + kLa1, Lap2 + La2)
!!$ Could include and if statement here to only calculate this part if it is required for this i0
!!$ May not be required if ang1 * ang2 == 0d0 for all lambda...
                 temp_fun(minfp1:maxfp1) = fp1(minfp1:maxfp1) * weight(minfp1:maxfp1)
                 do lambda = lambda_min, lambda_max
                    ang1 = 0d0
                    ang2 = 0d0
                    ang1 = dble((-1)**(mu))*Yint(dble(Lap1),dble(Mp1),dble(lambda),dble(-mu),dble(kLa1),dble(kM1))
                    ang2 = Yint(dble(Lap2),dble(Mp2),dble(lambda),dble(mu),dble(La2),dble(M2))
                    if ( ang1 * ang2 == 0d0 ) cycle
!!$ Integrate over coordinate 2 
                    jr1 = max(minfp1,epot_lambda_iminr(lambda))
                    jr2 = min(maxfp1,epot_lambda_imaxr(lambda))
                    temp2_fun(jr1:jr2) = temp_fun(jr1:jr2) * epot_lambda(lambda,jr1:jr2)
!!$ Initial COORDINATE 1 |kLM(r_1)>           
                    do ki = 1, nqmi
                       kqi = npk_nchi + ki - 1
                       tn1 => chil%b(kqi)
                       f1 => fpointer(tn1)
                       !minf1 = get_minf(tn1)
                       minf1 = tn1%minf
                       !maxf1 = get_maxf(tn1)
                       maxf1 = tn1%maxf
                       ir1 = max(jr1,minf1)
                       ir2 = min(jr2,maxf1)
                       if(ir2< ir1) cycle
!!$ Integrate over coordinate 1 
!                       temp = SUM(fp1(ir1:ir2) * epot_i(ir1:ir2) * f1(ir1:ir2) * weight(ir1:ir2))
		       temp = SUM( f1(ir1:ir2) * temp2_fun(ir1:ir2) )
                       temp = temp_CI * temp * ang1 * ang2
!!$ Overlap coordinate 0 
                       do kf = 1, nqmf
                          kqf = npk_nchf + kf - 1
                          do i0 = 1, nam0
                             !ind0 = get_na(TargetStates%b(MO_0),i0) ! Index to underlying atomic orbitals  
                             ind0 = TargetStates%b(MO_0)%na(i0) ! Index to underlying atomic orbitals  
                             !CI_0 = get_CI(TargetStates%b(MO_0),i0)
                             CI_0 = TargetStates%b(MO_0)%CI(i0)
                             vmatt(kf,ki)= vmatt(kf,ki) - CI_0 * temp * ortchil(ind0,kqf)
                          end do ! i0
                       end do ! kf
                    end do ! ki cord 1        
                 end do ! lambda 
              end do ! ip1
           end do ! i2 
        end do ! ip2
     end do ! nst_con
  end do ! nstp_con


!!$ Calculating -Cp*C*<kp0|n0><np1|k1><np2|H_2|n2>
!!$ Looping over FINAL Molecular State orbitals
  do nstp_con =  1, Npcon
     !CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con) ! 2e config CI
     CI_stp = TargetStates2el%b(nstf)%CI(nstp_con) ! 2e config CI
     !MOp1 = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
     MOp1 = TargetStates2el%b(nstf)%na(nstp_con) ! Index 1e molecular orbital
     !MOp2 = get_na(TargetStates2el%b(nstf),nstp_con,2) ! Index to 1e molecular orbital
     MOp2 = TargetStates2el%b(nstf)%nb(nstp_con) ! Index to 1e molecular orbital
     !namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
     namp1 = TargetStates%b(MOp1)%nam             ! Number of atomic orbitals that make up molecular orbital 
     !namp2 = get_nam(TargetStates%b(MOp2))
     namp2 = TargetStates%b(MOp2)%nam
!!$ Looping over INITIAL Molecular State configurations 
     do nst_con =  1, Ncon
        !CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
        CI_st = TargetStates2el%b(nsti)%CI(nst_con)
        !MO_0 = get_na(TargetStates2el%b(nsti),nst_con,1)
        MO_0 = TargetStates2el%b(nsti)%na(nst_con)
        !MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)
        MO2 = TargetStates2el%b(nsti)%nb(nst_con)
        !nam0 = get_nam(TargetStates%b(MO_0))
        nam0 = TargetStates%b(MO_0)%nam
        !nam2 = get_nam(TargetStates%b(MO2))
        nam2 = TargetStates%b(MO2)%nam
!!$ FINAL STATE COORDINATE 2 for BETA      
        do ip2 = 1, namp2
           !indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
           indp2 = TargetStates%b(MOp2)%na(ip2) ! Index to underlying atomic orbitals 
           !CI_2p = get_CI(TargetStates%b(MOp2),ip2)
           CI_2p = TargetStates%b(MOp2)%CI(ip2)
!!$ INITIAL STATE COORDINATE 2 for DELTA          
           do i2 = 1, nam2
              !ind2 = get_na(TargetStates%b(MO2),i2)
              ind2 = TargetStates%b(MO2)%na(i2)
              !CI_2 = get_CI(TargetStates%b(MO2),i2)
              CI_2 = TargetStates%b(MO2)%CI(i2)
              if ( bst_nr%ham1el(indp2,ind2) == 0d0 ) cycle
!!$ Overlap coordinate 1 <np1|k1>
              do ki = 1, nqmi
                 kqi = npk_nchi + ki - 1
                 do ip1 = 1, namp1
                    !indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
                    indp1 = TargetStates%b(MOp1)%na(ip1)  ! Index to underlying atomic orbitals 
                    !CI_1p = get_CI(TargetStates%b(MOp1),ip1)
                    CI_1p = TargetStates%b(MOp1)%CI(ip1)
                    if (ortchil(indp1,kqi) == 0d0) cycle
!!$ Overlap coordinate 0 <kp0|n0>
                    do kf = 1, nqmf
                       kqf = npk_nchf + kf - 1
                       do i0 = 1, nam0
                          !ind0 = get_na(TargetStates%b(MO_0),i0) ! Index to underlying atomic orbitals  
                          ind0 = TargetStates%b(MO_0)%na(i0) ! Index to underlying atomic orbitals  
                          !CI_0 = get_CI(TargetStates%b(MO_0),i0)
                          CI_0 = TargetStates%b(MO_0)%CI(i0)
                          if (ortchil(ind0,kqf) == 0d0 ) cycle
                          temp_CI = CI_stp * CI_1p * CI_2p * CI_st * CI_0 * CI_2
                          if ( temp_CI == 0d0 ) cycle
                          temp = temp_CI * ortchil(ind0,kqf) * ortchil(indp1,kqi) 
                          vmatt(kf,ki)= vmatt(kf,ki) - temp * bst_nr%ham1el(indp2,ind2) 
                       end do ! i0
                    end do ! kf coord 0
                 end do ! ip1
              end do ! ki coord 1
           end do ! i2 
        end do ! ip2
     end do ! nst_con
  end do ! nstp_con


!!$ Spin and other coefficients
!!$ Need to be done before projection operator overlaps
!!$ + 2. Spin_coeff.(E(1-theta)-H)P_01, vmatt=(E(1-theta)-H)P_01 
  temp = sqrt(dble(2*Spinp+1)) * sqrt(dble(2*Spin+1))
  temp = COF6J(0.5d0,0.5d0,dble(Spinp),0.5,dble(rspin_tot),dble(Spin)) * temp
  temp = 2d0 * dble((-1)**(Spinp+Spin+1)) * temp

  vmatt(:,:) = temp * vmatt(:,:)

!!$ Projection operator part of non-uniqueness
!!$ 2(E-theta+theta)P_01
!!$ -2 E * theta * <Phi_f|Phi_i> * <kp| I_0 |k>
!  if (nsti == nstf .AND. theta /= 0d0 ) then
  if ( theta /= 0d0 .AND.  Spinp == Spin ) then
!!$ Looping over FINAL (p) Molecular State configorations 
     do nstp_con =  1, Npcon
        !CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con) ! 2e config CI
        CI_stp = TargetStates2el%b(nstf)%CI(nstp_con) ! 2e config CI
        !MOp1 = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
        MOp1 = TargetStates2el%b(nstf)%na(nstp_con) ! Index 1e molecular orbital
        !Mp1 = get_ang_mom_proj(TargetStates%b(MOp1))      ! Angular momentum projection 1e molecular orbital
        Mp1 = TargetStates%b(MOp1)%m      ! Angular momentum projection 1e molecular orbital
        !MOp2 = get_na(TargetStates2el%b(nstf),nstp_con,2) ! Index to 1e molecular orbital
        MOp2 = TargetStates2el%b(nstf)%nb(nstp_con) ! Index to 1e molecular orbital
        !Mp2 = get_ang_mom_proj(TargetStates%b(MOp2))      ! Angular momentum projection 1e molecular orbital
        Mp2 = TargetStates%b(MOp2)%m      ! Angular momentum projection 1e molecular orbital
        !namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
        namp1 = TargetStates%b(MOp1)%nam             ! Number of atomic orbitals that make up molecular orbital 
        !namp2 = get_nam(TargetStates%b(MOp2))
        namp2 = TargetStates%b(MOp2)%nam

!!$ Looping over INITIAL Molecular State configurations 
        do nst_con =  1, Ncon
           !CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
           CI_st = TargetStates2el%b(nsti)%CI(nst_con)
           !MO1 = get_na(TargetStates2el%b(nsti),nst_con,1)
           MO1 = TargetStates2el%b(nsti)%na(nst_con)
           !M1 = get_ang_mom_proj(TargetStates%b(MO1))
           M1 = TargetStates%b(MO1)%m
           !MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)
           MO2 = TargetStates2el%b(nsti)%nb(nst_con)
           !M2 = get_ang_mom_proj(TargetStates%b(MO2))
           M2 = TargetStates%b(MO2)%m
           !nam1 = get_nam(TargetStates%b(MO1))
           nam1 = TargetStates%b(MO1)%nam
           !nam2 = get_nam(TargetStates%b(MO2))
           nam2 = TargetStates%b(MO2)%nam
           if ( Mp1 /= M1 .OR. Mp2 /= M2 ) cycle
           if ( is_core_MO(MO2) == 0 ) then ! cycle
!                print*,"Proj M02 delta",MO2
                cycle
           end if
           if ( is_core_MO(MO2) == 1 .AND. is_core_MO(MOp1) == 0 ) then ! cycle
!                print*,"Proj M02 delta, MOp1",MO2, MOp1
                cycle
           end if
!!$ overlap 2 matrix element CI_stp.CI_2p.CI_st.CI_2.<np2|n2> 
           do ip2 = 1, namp2 
              !indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
              indp2 = TargetStates%b(MOp2)%na(ip2) ! Index to underlying atomic orbitals 
              tnp2 => bst_nr%b(indp2)
              fp2 => fpointer(tnp2)
              !Lap2 = get_ang_mom(tnp2)
              Lap2 = tnp2%l
              !CI_2p = get_CI(TargetStates%b(MOp2),ip2)
              CI_2p = TargetStates%b(MOp2)%CI(ip2)
              !minfp2 = get_minf(tnp2)
              minfp2 = tnp2%minf
              !maxfp2 = get_maxf(tnp2)
              maxfp2 = tnp2%maxf
              do i2 = 1, nam2
                 !ind2 = get_na(TargetStates%b(MO2),i2)
                 ind2 = TargetStates%b(MO2)%na(i2)
                 tn2 => bst_nr%b(ind2)
                 f2 => fpointer(tn2)
                 !La2 = get_ang_mom(tn2)
                 La2 = tn2%l
                 !minf2 = get_minf(tn2)
                 minf2 = tn2%minf
                 !maxf2 = get_maxf(tn2)
                 maxf2 = tn2%maxf
                 !CI_2 = get_CI(TargetStates%b(MO2),i2)
                 CI_2 = TargetStates%b(MO2)%CI(i2)
                 if (Lap2 /= La2 ) cycle
                 temp_CI = CI_stp * CI_2p * CI_st * CI_2
                 if ( temp_CI == 0d0 ) cycle
                 ir1 = max(minf2,minfp2)
                 ir2 = min(maxf2,maxfp2)
                 if(ir2<ir1) cycle
!                 overlap2 = temp_CI * SUM( fp2(ir1:ir2) *  f2(ir1:ir2) * weight(ir1:ir2))
!!$ MARK: Use the below but needs testing and checking MARK
                 overlap2 = temp_CI * bst_nr%ortint(indp2,ind2)

                 if ( overlap2 == 0d0 ) cycle
!!$ overlap1=CI_1p.CI_1.<np1|n1>
                 do ip1 = 1, namp1
                    !indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
                    indp1 = TargetStates%b(MOp1)%na(ip1)  ! Index to underlying atomic orbitals 
                    tnp1 => bst_nr%b(indp1)   !           
                    fp1 => fpointer(tnp1)     ! One electron functions
                    !Lap1 = get_ang_mom(tnp1)  ! Gets Angular momentum A.O.
                    Lap1 = tnp1%l  ! Gets Angular momentum A.O.
                    !minfp1 = get_minf(tnp1)
                    minfp1 = tnp1%minf
                    !maxfp1 = get_maxf(tnp1)
                    maxfp1 = tnp1%maxf
                    !CI_1p = get_CI(TargetStates%b(MOp1),ip1)
                    CI_1p = TargetStates%b(MOp1)%CI(ip1)
                    do i1 = 1, nam1
                       !ind1 = get_na(TargetStates%b(MO1),i1)  ! Index to underlying atomic orbitals 
                       ind1 = TargetStates%b(MO1)%na(i1)  ! Index to underlying atomic orbitals 
                       tn1 => bst_nr%b(ind1)   !           
                       f1 => fpointer(tn1)     ! One electron functions
                       !La1 = get_ang_mom(tn1)  ! Gets Angular momentum A.O.
                       La1 = tn1%l  ! Gets Angular momentum A.O.
                       if ( Lap1 /= La1 ) cycle
                       !minf1 = get_minf(tn1)
                       minf1 = tn1%minf
                       !maxf1 = get_maxf(tn1)
                       maxf1 = tn1%maxf
                       !CI_1 = get_CI(TargetStates%b(MO1),i1)
                       CI_1 = TargetStates%b(MO1)%CI(i1)
                       temp_CI = CI_1p * CI_1
                       if ( temp_CI == 0d0 ) cycle
                       ir1 = max(minf1,minfp1)
                       ir2 = min(maxf1,maxfp1)
                       if(ir2<ir1) cycle
!                       overlap1 = temp_CI * SUM( fp1(ir1:ir2) *  f1(ir1:ir2) * weight(ir1:ir2))
                       overlap1 = temp_CI * bst_nr%ortint(indp1,ind1)
                       if ( overlap1 == 0d0 ) cycle
                       do ind0 = 1, TargState_Orthonorm%Nmax  ! Loop over all one-electron target states 
!!$ Note the difference between is_core_MO and is_core_POMO
!!$ CORE ORBITALS CHECKED HERE FOR MULTI CORE CASE 
                          if (is_core_MO(MO2) == 1 .AND. is_core_POMO(ind0) == 0) cycle

                          do ki = 1, nqmi
                             kqi = npk_nchi + ki - 1
                             do kf = 1, nqmf
                                kqf = npk_nchf + kf - 1
                                temp = ovlp_OrtnrmMO_k(ind0,kqi) * ovlp_OrtnrmMO_k(ind0,kqf)
                                temp = temp * overlap1 * overlap2
                                vmatt(kf,ki) = vmatt(kf,ki) - 2d0 * temp * Etot * theta 
                             end do ! kf
                          end do ! ki
!                          print*,"n, <n|k_1>",ind0,ovlp_OrtnrmMO_k(ind0,1)
                       end do ! ind0 I_0 
                    end do ! i1
                 end do ! ip1
              end do ! i2
           end do !ip2
        end do ! nst_con
     end do ! nstp_con
  end if ! non-uniqueness? and spin

  vmatt(:,:) = vmatt(:,:) * (2d0/pi) !!! / (rkf*rki) dealt with in scat.f 

end subroutine vexch_2e


subroutine vdirect_2e_oid(rspin_tot, nchf,nchi, nqmf,nqmi,npk_nchf,npk_nchi,cont, nr,rhoVec,weightVec, vmatt, numint)
  !  
  !   Spheroidal Direct V-matrix Elements.
  !
  use channels
  use input_data
  use one_electron_func
  use ovlpste1me
  use spheroidal_class
  use sturmian_class
  use target_states
  implicit none

  real*8, intent(in) :: rspin_tot   ! Total spin for the current partial wave.
  integer, intent(in) :: nchf, nchi   ! Final and initial channel numbers.
  integer, intent(in) :: nqmf, nqmi   ! Number of k-grid points per channel.
  integer, intent(in) :: npk_nchf, npk_nchi   ! On-shell points per channel.
  type(spheroidal_basis), intent(in) :: cont   ! Spheroidal continuum waves.
  integer, intent(in) :: nr   ! Number of points on the radial grid.
  real*8, dimension(1:nr), intent(in) :: rhoVec, weightVec   ! Radial vectors.
  real*8, dimension(nqmf,nqmi), intent(out) :: vmatt   ! Direct V-matrix.
  real*8, optional, intent(inout) :: numint

  type(spheroidal_fn), pointer :: wavef,wavei
  type(state), pointer :: statef,statei, orbC,orbA
  type(sturmian_nr), pointer :: sturmC,sturmA
  logical :: elastic
  integer :: nf,ni, lammin,lammax,jlam,nlam, l,m, nstf,nsti, usef,usei, jA,jB,jC,jD, cfgf,cfgi, spC,spA, indC,indA, i1,i2, j1,j2, jf,ji, termf,termi, i
  integer, dimension(:), pointer :: lfVec,liVec
  real*8 :: pi,R,Zproj,Zplus, Spin, lf,li, mf,mi, lam,mu, overlap2, CIf,CIi, lC,mC, lA,mA, rhoCoeff,etaCoeff, CIC,CIA, Df,Di, coeff,factor
  real*8, dimension(1:nr) :: ppRVec,ZVec, inVec,outVec
  real*8, dimension(:), pointer :: vecC,vecA, vecf,veci
  real*8, dimension(:,:), allocatable :: formMat, rhoMat,etaMat
  logical :: load

logical :: faster=.true.
integer, dimension(:), allocatable :: i1Vec,i2Vec
real*8, dimension(:,:,:), allocatable :: rhoCoeffMat,etaCoeffMat

  if(present(numint)) then
    load = .true.
    numint = 0
  else
    load = .false.
  endif

  vmatt(:,:) = 0d0
  pi = acos(-1d0)
  R = data_in%Rd
  Zproj = data_in%Zproj
  Zplus = data_in%Z1 + data_in%Z2
  ppRVec(:) = rhoVec(:) * (rhoVec(:)+R)

  ! Get the largest l-vectors from the last (highest momentum) k-grid point.
  wavef => get_spheroidal_fn(cont, npk_nchf+nqmf-1)
  wavei => get_spheroidal_fn(cont, npk_nchi+nqmi-1)
  nf = get_num_terms(wavef); ni = get_num_terms(wavei)
  lfVec => get_l_vec(wavef); liVec => get_l_vec(wavei)
  lf = dble(lfVec(nf)); li = dble(liVec(ni))
  mf = dble(get_m(wavef)); mi = dble(get_m(wavei))

  nstf = st_ch(nchf); statef => TargetStates2el%b(nstf)
  nsti = st_ch(nchi); statei => TargetStates2el%b(nsti)

  Spin = get_spin(statei)
  if (get_spin(statef) /= Spin) return   ! Conservation of spin.
  if (rspin_tot==1.5d0 .and. Spin==0d0) return   ! Must be triplet for S=3/2.
  if (mf + get_ang_mom_proj(statef) /= mi + get_ang_mom_proj(statei)) return

  ! Set up the limits of lambda and mu.
  mu = mi - mf; m = nint(mu)
  lammin = abs(m) + modulo(nint(lf+li+mu), 2)
  lammax = nint(lf+li) + 2   ! Technically correct, but you don't need...
  lammax = min(lammax, 2*get_max_L(bst_nr)+2)    ! ...the higher lambda's.
  !lammax = lammin + data_in%Lpmax   ! New rule to speed things along.
  if (lammin > lammax) return
  if (modulo(lammin+lammax,2) == 1) lammax=lammax-1   ! Match the parity. 
  nlam = (lammax-lammin)/2 + 1   ! Number of lambdas to calculate for.
  allocate( formMat(nr,nlam), rhoMat(nf,ni),etaMat(nf,ni) )
  formMat(:,:) = 0d0   ! This will contain the form results for each lambda.

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!


  ! First step -- integrate in space 1 & 2 to obtain form vectors in space 0.
  do usef = 1, get_nusemax(statef)
     jC = get_nuse(statef,usef)
     orbC => TargetStates%b(jC)
     mC = dble(get_ang_mom_proj(orbC))
     do usei = 1, get_nusemax(statei)
        jA = get_nuse(statei,usei)
        orbA => TargetStates%b(jA)
        mA = dble(get_ang_mom_proj(orbA))
        if (mC - mA /= mu) cycle   ! Conservation of orbital angular momentum projection.

        ! Evaluate space 2 overlaps.
        overlap2 = 0d0
        do cfgf = 1, get_nam(statef)
           if (get_na(statef,cfgf,1) /= jC) cycle
           jD = get_na(statef,cfgf,2)
           CIf = get_CI(statef,cfgf)
           do cfgi = 1, get_nam(statei)
              if (get_na(statei,cfgi,1) /= jA) cycle
              jB = get_na(statei,cfgi,2)
              CIi = get_CI(statei,cfgi)

              overlap2 = overlap2 + CIf*CIi*ovlpst(jD,jB)
           end do
        end do
        if(abs(overlap2) < data_in%expcut) overlap2 = 0.0d0
        if (overlap2 == 0d0) cycle


        ! Evaluate space 1 overlaps.
        do spC = 1, get_nam(orbC)
           indC = get_na(orbC,spC); sturmC => bst_nr%b(indC)
           CIC = get_CI(orbC,spC); vecC => fpointer(sturmC)
           lC = dble(get_ang_mom(sturmC))
           do spA = 1, get_nam(orbA)
              indA = get_na(orbA,spA); sturmA => bst_nr%b(indA)
              CIA = get_CI(orbA,spA); vecA => fpointer(sturmA)
              lA = dble(get_ang_mom(sturmA))

              i1 = max( get_minf(sturmC), get_minf(sturmA) )
              i2 = min( get_maxf(sturmC), get_maxf(sturmA) )
              ! 2x for two-electron calculation (V01=V02).
              coeff = 2d0 * overlap2 * CIC*CIA

              do jlam = 1, nlam
                 l = lammin + (jlam-1)*2; lam = dble(l) 

                 call angular_coeffs_3(lC,mC, lam,mu, lA,mA, rhoCoeff,etaCoeff)
                 if (rhoCoeff==0d0 .and. etaCoeff==0d0) cycle
!                 angVec(i1:i2) = rhoCoeff*ppRVec(i1:i2) + etaCoeff*R*R

                 ! NOTE: Spheroidal target states are calculated as-is.
                 !       In spherical, 1/r has been factored out.
                 inVec(:) = 0d0
!                 inVec(i1:i2) = vecC(i1:i2)*vecA(i1:i2) * angVec(i1:i2)*weightVec(i1:i2)


                 if(.not.load) then

                   !LHS: replaced old spheroidal form subroutine with more accurate form routine.
                   !     - the more accurate form does not need integration weights in inVec 
                   inVec(i1:i2) = vecC(i1:i2)*vecA(i1:i2) * (rhoCoeff*ppRVec(i1:i2)+etaCoeff*R*R)
                   call form_spheroidal_accurate(l,abs(m), inVec,i1,i2,nr, outVec,j1,j2)
                   !inVec(i1:i2) = vecC(i1:i2)*vecA(i1:i2) * (rhoCoeff*ppRVec(i1:i2)+etaCoeff*R*R) * weightVec(i1:i2)
                   !call form_spheroidal(l,abs(m), inVec,i1,i2,nr, outVec,j1,j2)
                 
                   formMat(j1:j2,jlam) = formMat(j1:j2,jlam) + (2d0*lam+1d0)*coeff*outVec(j1:j2)
                 else
                   numint = numint + (i2-i1+1)
                   j1 = i1
                   j2 = i2
                 endif

              end do ! jlam

           end do ! spA
        end do ! spC

     end do ! usei
  end do ! usef

  coeff = -Zproj * 2d0/pi
  ZVec(:) = Zplus * (rhoVec(:) + R/2d0)


if (faster) then !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
! THIS IS THE FASTER METHOD.

  do jlam = 1, nlam
     lam = dble( lammin + (jlam-1)*2 )
     elastic = nstf==nsti .and. lam==0d0

     inVec(:) = formMat(:,jlam)
     call minmaxi(inVec,nr, i1,i2)
     ! Check that the form vector is nonzero. Otherwise, we can skip...
     ! ...that is, unless we are about to incorporate V_0 into it (elastic).
     if (i1>=i2 .and. .not.elastic) cycle

     do jf = 1, nf
        lf = dble(lfVec(jf))
        do ji = 1, ni
           li = dble(liVec(ji))

           call angular_coeffs_3(lf,mf, lam,-mu, li,mi, rhoCoeff,etaCoeff)
           rhoMat(jf,ji) = rhoCoeff
           etaMat(jf,ji) = etaCoeff
        end do
     end do

     inVec = inVec * weightVec

     do jf = 1, nqmf
        wavef => get_spheroidal_fn(cont, npk_nchf+jf-1)
        vecf => get_rad_vec(wavef)
        do ji = 1, nqmi
           wavei => get_spheroidal_fn(cont, npk_nchi+ji-1)
           veci => get_rad_vec(wavei)

           j1 = max( get_rad_min(wavef), get_rad_min(wavei) )
           j2 = min( get_rad_max(wavef), get_rad_max(wavei) )
           if (j1 >= j2) cycle

           factor = 0d0; rhoCoeff = 0d0; etaCoeff = 0d0
           do termf = 1, get_num_terms(wavef)
              l = lfVec(termf)
              if (l /= get_term_l(wavef,termf)) stop 'lf error'
              Df = get_term_D(wavef,termf)

              do termi = 1, get_num_terms(wavei)
                 if (liVec(termi) /= get_term_l(wavei,termi)) stop 'li error'
                 Di = get_term_D(wavei,termi)
                 rhoCoeff = rhoCoeff + Df*Di*rhoMat(termf,termi)
                 etaCoeff = etaCoeff + Df*Di*etaMat(termf,termi)
                 if (elastic .and. liVec(termi)==l) factor = factor + Df*Di
              end do
           end do

           if (rhoCoeff==0d0 .and. etaCoeff==0d0 .and. factor==0d0) cycle
         
           if(.not.load) then
             !projectile-nuclear term:
             if (factor /= 0d0) vmatt(jf,ji) = vmatt(jf,ji) - coeff*factor*sum( vecf(j1:j2)*veci(j1:j2) * ZVec(j1:j2) * weightVec(j1:j2) )
             if(data_in%Z1/=data_in%Z2) error stop '*** Check Eq. (4.100) of Jeremy''s thesis, think we''re missing something here (vdirect_2e_oid)'
          
             !Now update integration bounds to account for target wave function bounds too
             j1 = max(i1,j1); j2 = min(i2,j2)
             if (j1 >= j2) cycle
             !vmatt(jf,ji) = vmatt(jf,ji) + coeff*sum( vecf(j1:j2)*veci(j1:j2) * (rhoCoeff*ppRVec(j1:j2)+etaCoeff*R*R) * inVec(j1:j2) * weightVec(j1:j2) )
             vmatt(jf,ji) = vmatt(jf,ji) + coeff*sum( vecf(j1:j2)*veci(j1:j2) * (rhoCoeff*ppRVec(j1:j2)+etaCoeff*R*R) * inVec(j1:j2) )
             
           else
             numint = numint + (j2-j1+1)
             if(factor /= 0.0d0) numint = numint + (j2-j1+1)
           endif

        end do ! ji
     end do ! jf

  end do ! jlam
  deallocate(formMat, rhoMat,etaMat)

else !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!
! THIS IS THE SLOWER METHOD.

  allocate( i1Vec(nlam),i2Vec(nlam), rhoCoeffMat(nf,nlam,ni),etaCoeffMat(nf,nlam,ni) )

  do jlam = 1, nlam
     inVec(:) = formMat(:,jlam)
     call minmaxi(inVec,nr, i1,i2)
     i1Vec(jlam) = i1
     i2Vec(jlam) = i2
  end do

  do jf = 1, nf
     lf = dble( lfVec(jf) )
     do jlam = 1, nlam
        lam = dble( lammin + 2*(jlam-1) )
        do ji = 1, ni
           li = dble( liVec(ji) )
           
           call angular_coeffs_3(lf,mf, lam,-mu, li,mi, rhoCoeff,etaCoeff)
           rhoCoeffMat(jf,jlam,ji) = rhoCoeff
           etaCoeffMat(jf,jlam,ji) = etaCoeff
        end do
     end do
  end do

  do jf = 1, nqmf
     wavef => get_spheroidal_fn(cont, npk_nchf+jf-1)
     vecf => get_rad_vec(wavef)
     do ji = 1, nqmi
        wavei => get_spheroidal_fn(cont, npk_nchi+ji-1)
        veci => get_rad_vec(wavei)

        i1 = max( get_rad_min(wavef), get_rad_min(wavei) )
        i2 = min( get_rad_max(wavef), get_rad_max(wavei) )
        if (i1 >= i2) cycle

        outVec(:) = 0d0
        factor = 0d0
        do jlam = 1, nlam
           lam = dble( lammin + 2*(jlam-1) )
           elastic = nstf==nsti .and. lam==0d0

           j1 = max(i1, i1Vec(jlam)); j2 = min(i2, i2Vec(jlam))
           if (j1 >= j2) cycle

           rhoCoeff = 0d0; etaCoeff = 0d0
           do termf = 1, get_num_terms(wavef)
              lf = dble( get_term_l(wavef,termf) )
              Df = get_term_D(wavef,termf)
              do termi = 1, get_num_terms(wavei)
                 li = dble( get_term_l(wavei,termi) )
                 Di = get_term_D(wavei,termi)
                 
                 if (elastic .and. lf==li) factor = factor + Df*Di
                 rhoCoeff = rhoCoeff + Df*Di*rhoCoeffMat(termf,jlam,termi)
                 etaCoeff = etaCoeff + Df*Di*etaCoeffMat(termf,jlam,termi)
              end do
           end do

           outVec(j1:j2) = outVec(j1:j2) + formMat(j1:j2,jlam) * (rhoCoeff*ppRVec(j1:j2)+etaCoeff*R*R)

        end do

        vmatt(jf,ji) = coeff*sum( vecf(i1:i2)*veci(i1:i2) * outVec(i1:i2) * weightVec(i1:i2) )
        if (factor /= 0d0) vmatt(jf,ji) = vmatt(jf,ji) - coeff*factor*sum( vecf(i1:i2)*veci(i1:i2) * ZVec(i1:i2) * weightVec(i1:i2) )

     end do ! ji
  end do ! jf
  deallocate(formMat, i1Vec,i2Vec, rhoCoeffMat,etaCoeffMat)

end if !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!


end subroutine vdirect_2e_oid

subroutine vexch_2e_oid(nchf,nchi, nqmf,nqmi,npk_nchf,npk_nchi,cont_oid, nOrb,nWaves,ortchil,flchil, nr,rhoVec,weightVec, Etot,R,Zproj,S,theta, vmatt, P01mat, numint)
  !
  ! We start with the exchange potential:
  !   V^E = 2 (E - H_0 - H_1 - H_2 - V_01 - V_02 - V_12 - Z^2/R) P_r0r1
  !
  ! and divide it into several stages while enforcing antisymmetry:
  !   V^E = -2 [V_0 + V_1 + V_01] P_r0r1                 : two-two
  !         +2 [E(1-theta) - K_0 - K_1 - Z^2/R] P_r0r1   : scalars
  !         -2 [V_02 + V_12] P_r0r1                      : one-three
  !         -2 [H_2] P_r0r1                              : analytical
  !         -2 E theta I_0                               : projection
  !
  ! JSS
  !
  use channels
  use input_data
  use one_electron_func
  use ovlpste1me
  use spheroidal_class
  use state_class
  use sturmian_class
  use target_states
  use vmat_exch_module
  use grid_radial
  implicit none

  integer, intent(in) :: nchf,nchi, nqmf,nqmi, npk_nchf,npk_nchi, nOrb,nWaves,nr
  type(spheroidal_basis), intent(in) :: cont_oid
  real*8 :: Etot, R, Zproj, S, theta
  real*8, dimension(nOrb,nWaves), intent(in) :: ortchil,flchil
  real*8, dimension(nr), intent(in) :: rhoVec,weightVec
  real*8, dimension(nqmf,nqmi), intent(out) :: vmatt
  real*8, dimension(nqmf,nqmi), optional, intent(out) :: P01mat
  real*8, optional, intent(inout) :: numint

  type(state), pointer :: stf,sti, stA,stB,stC,stD
  type(spheroidal_fn), pointer :: wavef,wavei
  type(sturmian_nr), pointer :: spA,spB,spC,spD
  logical :: non_uniq
  integer :: nstf,nsti, nusef,nusei, jf,ji, nconfigf,nconfigi, configf,configi
  integer :: lfMax,liMax, lAMax,lBMax,lCMax,lDMax, lamMin,lamMax, lamMin_temp, lamMax_temp, jlf,jli, jlam,nf,ni,nlam, l,m
  integer :: ind, indA,indB,indC,indD, jkf,jki, indf,indi, jA,jB,jC,jD, nA,nB,nC,nD, orbA,orbB,orbC,orbD
  integer :: parf,pari,par, a1,a2, b1,b2, c1,c2, d1,d2, e1,e2, i1,i2
  integer, dimension(:), pointer :: lfVec,liVec
  real*8 :: pi, mf,mi, mA,mB,mC,mD, lam,mu, lf,li, lA,lB,lC,lD, rhoCoeff,etaCoeff, ef,ei,E, Spinf,Spini
  real*8 :: overlap0,overlap1,overlap2, CIf,CIi, CIA,CIB,CIC,CID, Df,Di, coeff, tmp, theta_unique, COF6J
  real*8, dimension(nr) :: ppRVec, inVec,outVec, formVec, invecP, invecQ
  real*8, dimension(:), pointer :: vecf,veci, vecA,vecB,vecC,vecD
  real*8, dimension(:), allocatable :: rhoVecf,etaVecf, rhoVeci,etaVeci
  real*8, dimension(:,:,:), allocatable :: rhoVecf2,etaVecf2, rhoVeci2,etaVeci2
  real*8, dimension(nr) :: vecAweight, vecCP, vecCQ
  real*8 :: Rsquared, Z1, Z2
  logical :: save_P01

  integer :: nn, ii
  real*8, dimension(:), pointer :: Pvec, Qvec
  integer :: imin, imax
  real*8 :: tmp1, tmp2, tmp3, tmp4, tmp5
  real*8, dimension(:,:,:), allocatable :: ppRvecf

  logical :: load

  save_P01 = .false. !TODO: remove this and temp_logical in argument list when issue is fixed
                     !      - for some reason PRESENT(P01mat) is always true even when P01mat is not supplied...

  vmatt(:,:) = 0d0
  if(save_P01) then
    P01mat = 0.0d0
  endif

  if(present(numint)) then
    load = .true.
    numint = 0
  else
    load = .false.
  endif

  non_uniq = data_in%non_uniq
  Z1 = data_in%Z1
  Z2 = data_in%Z2
  pi = acos(-1d0)
  ppRVec(:) = rhoVec(:) * (rhoVec(:)+R)   ! Used for the angular quadratics.
  Rsquared = R*R
  
  ! Parameters of the target states.
  nstf = st_ch(nchf); stf => TargetStates2el%b(nstf)
  nusef = get_nusemax(stf); nconfigf = get_nam(stf)
  Spinf = dble( get_spin(stf) )
  nsti = st_ch(nchi); sti => TargetStates2el%b(nsti)
  nusei = get_nusemax(sti); nconfigi = get_nam(sti)
  Spini = dble( get_spin(sti) )

  ! Parameters of the spheroidal continuum waves.
  wavef => get_spheroidal_fn(cont_oid, npk_nchf+nqmf-1)   ! The last waves
  wavei => get_spheroidal_fn(cont_oid, npk_nchi+nqmi-1)   ! (highest momentum).
  nf = get_num_terms(wavef); ni = get_num_terms(wavei)    ! Hence most terms.
  lfVec => get_l_vec(wavef); mf = dble( get_m(wavef) )
  liVec => get_l_vec(wavei); mi = dble( get_m(wavei) )
  lfMax = lfVec(nf); liMax = liVec(ni)
  parf = (-1)**lfMax; pari = (-1)**liMax
  
  ! Conservation of total parity Pi and total angular momentum projection M.
  if ( data_in%good_parity .and. parf*get_par(stf) /= pari*get_par(sti) ) return
  if ( mf+get_ang_mom_proj(stf) /= mi+get_ang_mom_proj(sti) ) return

  allocate(rhoVeci(ni),etaVeci(ni))
  allocate(rhoVecf(nf),etaVecf(nf))
  
  lamMax = 0
  
  !Determine lamMax so we can allocate rhovecf2 and etavecf2
  do jf = 1, nusef
     orbC = get_nuse(stf,jf); stC => TargetStates%b(orbC)
     mC = dble( get_ang_mom_proj(stC) )
     
     mu = mC - mi 
     m = nint(mu)
        
     lamMin = abs(m) 
     lamMax = max(lamMax, lamMin + data_in%Lpmax + 1)
  enddo
  !TODO
  lammax = 2*get_max_L(bst_nr)+2    ! ...the higher lambda's.
  
  allocate(rhovecf2(bst_nr%n,nf,0:lamMax), etavecf2(bst_nr%n,nf,0:lamMax)) !, ppRVecf(grid%nr,nqmf,nA))
  
  rhovecf2 = 0.0d0; etavecf2 = 0.0d0
  do ji = 1, nusei
    orbA = get_nuse(sti,ji); stA => TargetStates%b(orbA)
    nA = get_nam(stA)
    mA = dble( get_ang_mom_proj(stA) )
    do jA = 1, nA
      indA = get_na(stA,jA); spA => bst_nr%b(indA)
      lA = dble( get_ang_mom(spA) )
      do jlf = 1, nf
        lf = dble( lfVec(jlf) )
        do l=nint(abs(mA-mf)), lamMax
          lam = dble(l)
          call angular_coeffs_3(lf,mf, lam,-(mA-mf), lA,mA, rhoCoeff,etaCoeff)
          rhoVecf2(indA,jlf,l) = rhoCoeff 
          etaVecf2(indA,jlf,l) = etaCoeff 
        enddo
      enddo !jlf
    end do !jA
  enddo !ji


  ! Do the "two-two" potentials--named for the electron-electron interaction
  ! between two projectile waves and two "atomic" single particle functions.
  !   < kf lf mf, stf | -2 [V_0 + V_1 + V_01] P_r0r1 | ki li mi, sti >
  !
  ! Once we expand the states, we need the following terms:
  !   < kf lf mf | alpha lA mA >       = ortchil(indA,indf)
  !   < kf lf mf | V_0 | alpha lA mA > = flchil(indA,indf)
  !   < gamma lC mC | ki li mi >       = ortchil(indC,indi)
  !   < gamma lC mC | V_1 | ki li mi > = flchil(indC,indi)
  !   < delta | beta >                 = ovlpst(orbD,orbB)
  !   < kf lf mf, gamma lC mC | V_01 | alpha lA mA, ki li mi >
  !
  ! V_0 and V_1 should contain the asymptotic charge and distorting potential.

  do jf = 1, nusef
     orbC = get_nuse(stf,jf); stC => TargetStates%b(orbC)
     nC = get_nam(stC)
     mC = dble( get_ang_mom_proj(stC) )
     par = get_par(stC) * pari
     
     mu = mC - mi 
     m = nint(mu)
        
     lamMin = abs(m) 
     lamMax = lamMin + data_in%Lpmax + 1
     !TODO
     lammax = 2*get_max_L(bst_nr)+2    ! ...the higher lambda's.
     if (lamMin > lamMax) cycle
     nlam = lamMax - lamMin + 1
    
     ! Loop through the lambdas that will give nonzero contributions.
     do jlam = 1, nlam
        l = lammin + jlam-1; lam = dble(l)
   
        !Liam added - the form subroutine will do nothing if l > ltmax anyway so lets just skip 
        !everything below from the start if this condition is met
        if(l > grid%ltmax) cycle 

        nn = abs(m)
        do ii = 1, l
          nn = nn + ii
        end do
           
        PVec => grid%LegendreP(:,nn)
        QVec => grid%LegendreQ(:,nn)

        ! Single particle functions for electron 1 final orbital.
        do jC = 1, nC
           indC = get_na(stC,jC); spC => bst_nr%b(indC)
           CIC  = get_CI(stC,jC); lC = dble( get_ang_mom(spC) )
           vecC => fpointer(spC)
           a1 = get_minf(spC); a2 = get_maxf(spC)
           
           ! Get all possible l combinations for the angular quadratic.
           do jli = 1, ni
              li = dble( liVec(jli) )
              call angular_coeffs_3(lC,mC, lam,mu, li,mi, rhoCoeff,etaCoeff)
              rhoVeci(jli) = rhoCoeff
              etaVeci(jli) = etaCoeff
           end do

           vecCP(a1:a2) = vecC(a1:a2)*Pvec(a1:a2)
           vecCQ(a1:a2) = vecC(a1:a2)*Qvec(a1:a2)

           ! Spheroidal continuum waves for electron 1 initial channel.
           do jki = 1, nqmi
              indi = npk_nchi + jki - 1; wavei => cont_oid%b(indi)
              veci => get_rad_vec(wavei)
              b1 = max(a1, get_rad_min(wavei))
              b2 = min(a2, get_rad_max(wavei))
              if (b1 >= b2) cycle
              
              ! Build the angular quadratic from component l's.
              rhoCoeff = 0d0; etaCoeff = 0d0
              do jli = 1, get_num_terms(wavei)
                 Di = get_term_D(wavei,jli)
                 rhoCoeff = rhoCoeff + Di*rhoVeci(jli)
                 etaCoeff = etaCoeff + Di*etaVeci(jli)
              end do
              if (rhoCoeff==0d0 .and. (etaCoeff==0d0.or.R==0d0)) cycle
!              angVec(b1:b2) = rhoCoeff*ppRVec(b1:b2) + etaCoeff*R*R

              
              if(.not.load) then
                ! Integrate in space 1 for a form vector in space 0.
                
                !LHS: replaced old spheroidal form subroutine with more accurate form routine.
                !     - the more accurate form does not need integration weights in inVec 
                !inVec(b1:b2) = vecC(b1:b2)*veci(b1:b2) * (rhoCoeff*ppRVec(b1:b2)+etaCoeff*R*R)
                inVec(b1:b2) = veci(b1:b2) * (rhoCoeff*ppRVec(b1:b2)+etaCoeff*R*R)
                inVecP(b1:b2) = vecCP(b1:b2) * inVec(b1:b2)
                inVecQ(b1:b2) = vecCQ(b1:b2) * inVec(b1:b2)

                !call form_spheroidal_accurate(l,abs(m), inVec,b1,b2,nr, outVec,c1,c2)
                call form_spheroidal_accurate2(l,abs(m), inVecP,inVecQ,b1,b2,nr, outVec,c1,c2,load)
                !inVec(b1:b2) = vecC(b1:b2)*veci(b1:b2) * (rhoCoeff*ppRVec(b1:b2)+etaCoeff*R*R) * weightVec(b1:b2)
                !call form_spheroidal(l,abs(m), inVec,b1,b2,nr, outVec,c1,c2)
                
                if (c1 >= c2) cycle

                outvec(c1:c2) = outvec(c1:c2) * weightVec(c1:c2)
              else
                numint = numint + (b2-b1+1)
                c1 = b1
                c2 = b2
              endif !load
     
              do ji = 1, nusei

                orbA = get_nuse(sti,ji); stA => TargetStates%b(orbA)
                nA = get_nam(stA)
                mA = dble( get_ang_mom_proj(stA) )
                if (data_in%good_parity .and. get_par(stA)*parf /= par) cycle   ! Enforcing parity conservation.


                ! Determine the limits on lambda in wavef/orbA and orbC/wavei integrals.
                lAMax = get_ang_mom( bst_nr%b(get_na(stA,nA)) )
                lamMin_temp = abs(m) + modulo(lfMax+lAMax+m,2)
                lamMax_temp = lamMin_temp + data_in%Lpmax
                lammax_temp = 2*get_max_L(bst_nr)+2    ! ...the higher lambda's.
                !TODO
                if(mod(lamMin_temp + lamMax_temp,2)==1) lamMax_temp = lamMax_temp - 1
                !lambda sum originally started at lamMin = abs(m) + modulo(lfMax+lAMax+m,2) and went in steps of two
                !now it starts at abs(m) and goes in steps of one so we need to cycle any lambdas that aren't needed for this nA:
                if(mod(l,2) /= mod(lfMax+lAMax,2)) cycle
                if(l < lamMin_temp .or. l > lamMax_temp) cycle
                if (lamMin_temp > lamMax_temp) cycle
                ! Enforce (mf+mC) = (mi+mA) ===> (M'-mB) = (M-mD).
                !mu = mA - mf   ! Coordinate 0: << lf mf | lam -mu | lA mA >>
                !if (mC - mi /= mu) cycle   ! Coordinate 1: << lC mC | lam mu | li mi >>
                if(mA - mf /= mu) cycle

                ! Overlap of the second electron functions through orbital overlaps.
                overlap2 = 0d0
                do configf = 1, nconfigf
                   if (get_na(stf,configf,1) /= orbC) cycle
                   orbD = get_na(stf,configf,2)
                   CIf = get_CI(stf,configf)
                   do configi = 1, nconfigi
                      if (get_na(sti,configi,1) /= orbA) cycle
                      orbB = get_na(sti,configi,2)
                      CIi = get_CI(sti,configi)

                      overlap2 = overlap2 + CIf*CIi*ovlpst(orbB,orbD)
                   end do
                end do
                if(abs(overlap2)<data_in%expcut) overlap2 = 0.0d0
                if (overlap2 == 0d0) cycle   ! Enforces mB = mD ===> M' = M.

                 ! Single particle functions for electron 0 initial orbital.
                 do jA = 1, nA
                    indA = get_na(stA,jA); spA => bst_nr%b(indA)
                    CIA  = get_CI(stA,jA); lA = dble( get_ang_mom(spA) )
                    coeff = -(2d0*lam+1d0) * overlap2 * CIC*CIA
                    vecA => fpointer(spA)
                    d1 = max(c1, get_minf(spA)); d2 = min(c2, get_maxf(spA))
                    if (d1 >= d2) cycle

                    if(.not.load) vecAweight(d1:d2) = vecA(d1:d2)*outVec(d1:d2)

                    ! Spheroidal continuum waves for electron 0 final channel.
                    do jkf = 1, nqmf
                       indf = npk_nchf + jkf - 1; wavef => cont_oid%b(indf)
                       vecf => get_rad_vec(wavef)
                       e1 = max(d1, get_rad_min(wavef))
                       e2 = min(d2, get_rad_max(wavef))
                       if (e1 >= e2) cycle
                       
                       ! Build the angular quadratic from component l's.
                       rhoCoeff = 0d0; etaCoeff = 0d0
                       do jlf = 1, get_num_terms(wavef)
                          Df = get_term_D(wavef,jlf)
                          rhoCoeff = rhoCoeff + Df*rhoVecf2(indA,jlf,l)
                          etaCoeff = etaCoeff + Df*etaVecf2(indA,jlf,l)
                       end do
                       if (rhoCoeff==0d0 .and. (etaCoeff==0d0.or.R==0)) cycle

                       if(.not.load) then
                         !tmp = coeff * sum(vecf(e1:e2)*vecA(e1:e2)*outVec(e1:e2) * (rhoCoeff*ppRVec(e1:e2)+etaCoeff*R*R) * weightVec(e1:e2))
                         !tmp = coeff * dot_product(vecf_local_times_angular(e1:e2,jkf,jA), vecAweight(e1:e2))
                         tmp = coeff * sum(vecf(e1:e2) * vecAweight(e1:e2) * (rhoCoeff*ppRVec(e1:e2)+etaCoeff*Rsquared))
                       else
                         numint = numint + (e2-e1+1)
                       endif !load

                       if (lam == 0 .and. .not.load) then
                          tmp = tmp + coeff*flchil(indA,indf)*ortchil(indC,indi)
                          tmp = tmp + coeff*ortchil(indA,indf)*flchil(indC,indi)
                       end if

                       if(.not.load) vmatt(jkf,jki) = vmatt(jkf,jki) + tmp
                    end do ! jkf
                 end do !jA
              end do ! jki
           end do ! jC
        end do ! jlam
     end do ! ji
  end do ! jf

  ! Evaluate the "scalar" terms such as kinetic energy and nuclear-nuclear
  ! interaction, as well as the first part of non-uniqueness.
  !  < kf lf mf, stf | [E(1-theta) - K_0 - K_1 - Z^2/R] P_r0r1 | ki li mi, sti >
  !
  ! Using the operators transforms the terms in square brackets to the scalar
  !  [E(1-theta) - k_0^2 - k_1^2 - Z^2/R
  ! and leaves the overlaps
  !  < kf lf mf | alpha lA mA >   = ortchil(indA,indf)
  !  < gamma lC mC | ki li mi >   = ortchil(indC,indi)
  !  < delta lB mB | beta lB mB > = ovlpst(orbD,orbB)

  do configf = 1, nconfigf
     orbC = get_na(stf,configf,1); stC => TargetStates%b(orbC)
     mC = dble( get_ang_mom_proj(stC) ); nC = get_nam(stC)
     orbD = get_na(stf,configf,2); CIf = get_CI(stf,configf)
     if (data_in%good_parity .and. get_par(stC) /= pari) cycle
     if (mC /= mi) cycle

     do configi = 1, nconfigi
        orbA = get_na(sti,configi,1); stA => TargetStates%b(orbA)
        mA = dble( get_ang_mom_proj(stA) ); nA = get_nam(stA)
        orbB = get_na(sti,configi,2); CIi = get_CI(sti,configi)
        if (data_in%good_parity .and. get_par(stA) /= parf) cycle
        if (mA /= mf) cycle

        overlap2 = CIf*CIi * ovlpst(orbD,orbB)
        if (overlap2 == 0d0) cycle

        ! NEED TO CONFIRM THIS LOGIC.
!JSS        theta_unique = 1d0
        theta_unique = 1d0-theta
        if (non_uniq) then
           theta_unique = 1d0
           if (is_core_MO(orbB) == 2) theta_unique = 1d0 - theta
           if (is_core_MO(orbB)==1 .and. is_core_MO(orbC)/=0) theta_unique = 1d0 - theta
        end if

        do jC = 1, nC
           indC = get_na(stC,jC)
           CIC = get_CI(stC,jC)
           do jki = 1, nqmi
              indi = npk_nchi + jki - 1
              ei = get_energy( cont_oid%b(indi) )
              overlap1 = CIC * ortchil(indC,indi)
              if (overlap1 == 0d0) cycle

              do jA = 1, nA
                 indA = get_na(stA,jA)
                 CIA = get_CI(stA,jA)
                 do jkf = 1, nqmf
                    indf = npk_nchf + jkf - 1
                    ef = get_energy( cont_oid%b(indf) )
                    overlap0 = CIA * ortchil(indA,indf)
                    if (overlap0 == 0d0) cycle

                    E = Etot*theta_unique - ei - ef
                    if (R > 0d0) E = E - Z1*Z2/R
                    tmp = overlap0*overlap1*overlap2
                    vmatt(jkf,jki) = vmatt(jkf,jki) + E*tmp
                    if(save_P01) P01mat(jkf,jki) = P01mat(jkf,jki) + tmp
                 end do ! jkf
              end do ! jA
           end do ! jki
        end do ! jC
     end do ! configi
  end do ! configf

  ! Do the first of the "one-three" potentials--named for the electron-electron
  ! interaction between one projectile wave and two single particle functions.
  !   < kf lf mf, stf | -2 V_02 P_r0r1 | ki li mi, sti >
  !
  ! Once we expand the states, we need the following terms:
  !   < gamma lC mC | ki li mi >       = ortchil(indC,indi)
  !   < kf lf mf, delta lD mD | V_02 | alpha lA mA, beta lB mB >

  do configf = 1, nconfigf
     orbC = get_na(stf,configf,1); stC => TargetStates%b(orbC)
     orbD = get_na(stf,configf,2); stD => TargetStates%b(orbD)
     mC = dble( get_ang_mom_proj(stC) ); nC = get_nam(stC)
     mD = dble( get_ang_mom_proj(stD) ); nD = get_nam(stD)
     CIf = get_CI(stf,configf)

     overlap1 = 0d0
     jC = 0
     do while (overlap1==0d0 .and. jC<nC)
        jC = jC + 1; indC = get_na(stC,jC)
        jki = 0
        do while (overlap1==0d0 .and. jki<nqmi)
           jki = jki + 1; indi = npk_nchi + jki - 1
           overlap1 = ortchil(indC,indi)
        end do
     end do
     if (overlap1==0d0) cycle        


     do configi = 1, nconfigi
        orbA = get_na(sti,configi,1); stA => TargetStates%b(orbA)
        orbB = get_na(sti,configi,2); stB => TargetStates%b(orbB) 
        mA = dble( get_ang_mom_proj(stA) ); nA = get_nam(stA)
        mB = dble( get_ang_mom_proj(stB) ); nB = get_nam(stB)
        CIi = get_CI(sti,configi)

        par = get_par(stA) * parf
        mu = mA - mf; m = nint(mu)
        if (data_in%good_parity .and. get_par(stD)*get_par(stB) /= par) cycle   ! Conservation of parity.
        if (mD - mB /= mu) cycle   ! Conservation of orb. ang. mom. proj.

        lDMax = get_ang_mom( bst_nr%b(get_na(stD,nD)) )
        lBMax = get_ang_mom( bst_nr%b(get_na(stB,nB)) )
        lamMax = lDMax + lBMax + 2
        lamMin = abs(m) + modulo(m+lamMax,2)
        nlam = (lamMax - lamMin)/2 + 1

        do jlam = 1, nlam
           l = lamMin + 2*(jlam-1); lam = dble(l)
           formVec(:) = 0d0

           d1 = nr; d2 = 1
           do jD = 1, nD
              indD = get_na(stD,jD); spD => bst_nr%b(indD)
              CID = get_CI(stD,jD); lD = dble( get_ang_mom(spD) )
              vecD => fpointer(spD)
              a1 = get_minf(spD); a2 = get_maxf(spD)
   
              do jB = 1, nB
                 indB = get_na(stB,jB); spB => bst_nr%b(indB)
                 CIB = get_CI(stB,jB); lB = dble( get_ang_mom(spB) )
                 vecB => fpointer(spB)
                 b1 = max(a1, get_minf(spB)); b2 = min(a2, get_maxf(spB))
                 if (b1 >= b2) cycle

                 call angular_coeffs_3(lD,mD, lam,mu, lB,mB, rhoCoeff,etaCoeff)
                 if (rhoCoeff==0d0 .and. (etaCoeff==0d0.or.R==0d0)) cycle

!                 inVec(b1:b2) = vecD(b1:b2) * vecB(b1:b2) * angVec(b1:b2) * weightVec(b1:b2)
                 
                 !LHS: replaced old spheroidal form subroutine with more accurate form routine.
                 !     - the more accurate form does not need integration weights in inVec 
                 inVec(b1:b2) = vecD(b1:b2) * vecB(b1:b2) * (rhoCoeff*ppRVec(b1:b2) + etaCoeff*R*R)
                 call form_spheroidal_accurate(l,abs(m), inVec,b1,b2,nr, outVec,c1,c2)
                 !inVec(b1:b2) = vecD(b1:b2) * vecB(b1:b2) * (rhoCoeff*ppRVec(b1:b2) + etaCoeff*R*R) * weightVec(b1:b2)
                 !call form_spheroidal(l,abs(m), inVec,b1,b2,nr, outVec,c1,c2)
                 if (c1 >= c2) cycle

                 formVec(c1:c2) = formVec(c1:c2) + outVec(c1:c2)
                 d1 = min(c1,d1); d2 = max(c2,d2)

              end do ! jB
           end do ! jD
           if (d1 >= d2) cycle

           do jA = 1, nA
              indA = get_na(stA,jA); spA => bst_nr%b(indA)
              CIA = get_CI(stA,jA); lA = dble( get_ang_mom(spA) )
              vecA => fpointer(spA)
              a1 = max(d1, get_minf(spA)); a2 = min(d2, get_maxf(spA))
              if (a1 >= a2) cycle

              coeff = -(2d0*lam+1d0) * CIf*CIi * CIA * CID*CIB
              rhoVecf(:) = 0d0; etaVecf(:) = 0d0
              do jlf = 1, nf
                 lf = dble( lfVec(jlf) )
                 call angular_coeffs_3(lf,mf, lam,-mu, lA,mA, rhoCoeff,etaCoeff)
                 rhoVecf(jlf) = rhoCoeff
                 etaVecf(jlf) = etaCoeff
              end do

              do jkf = 1, nqmf
                 indf = npk_nchf + jkf - 1; wavef => cont_oid%b(indf)
                 vecf => get_rad_vec(wavef)
                 b1 = max(a1, get_rad_min(wavef)); b2 = min(a2, get_rad_max(wavef))
                 if (b1 >= b2) cycle

                 rhoCoeff = 0d0; etaCoeff = 0d0
                 do jlf = 1, get_num_terms(wavef)
                    Df = get_term_D(wavef,jlf)
                    rhoCoeff = rhoCoeff + Df*rhoVecf(jlf)
                    etaCoeff = etaCoeff + Df*etaVecf(jlf)
                 end do
                 if (rhoCoeff==0d0 .and. (etaCoeff==0d0.or.R==0d0)) cycle
!                 angVec(b1:b2) = rhoCoeff*ppRVec(b1:b2) + etaCoeff*R*R
!                 tmp = coeff * sum(vecf(b1:b2)*vecA(b1:b2) * formVec(b1:b2)*angVec(b1:b2)*weightVec(b1:b2))
                 tmp = coeff * sum(vecf(b1:b2)*vecA(b1:b2) * formVec(b1:b2)*(rhoCoeff*ppRVec(b1:b2)+etaCoeff*R*R)*weightVec(b1:b2))

                 do jC = 1, nC
                    indC = get_na(stC,jC); CIC = get_CI(stC,jC)
                    do jki = 1, nqmi
                       indi = npk_nchi + jki - 1
                       overlap1 = CIC * ortchil(indC,indi)
                       vmatt(jkf,jki) = vmatt(jkf,jki) + overlap1*tmp
                    end do ! jki
                 end do ! jC

              end do ! jkf
           end do ! jA
        end do ! jlam
     end do ! configi
  end do ! configf

  ! Do the second of the "one-three" potentials--named for the electron-electron
  ! interaction between one projectile wave and two single particle functions.
  !   < kf lf mf, stf | -2 V_12 P_r0r1 | ki li mi, sti >
  !
  ! Once we expand the states, we need the following terms:
  !   < kf lf mf, alpha lA mA >       = ortchil(indA,indf)
  !   < gamma lC mC, delta lD mD | V_12 | ki li mi, beta lB mB >

  do configi = 1, nconfigi
     orbA = get_na(sti,configi,1); stA => TargetStates%b(orbA)
     orbB = get_na(sti,configi,2); stB => TargetStates%b(orbB)
     mA = dble( get_ang_mom_proj(stA) ); nA = get_nam(stA)
     mB = dble( get_ang_mom_proj(stB) ); nB = get_nam(stB)
     CIi = get_CI(sti,configi)
     
     overlap0 = 0d0
     jA = 0
     do while (overlap0==0d0 .and. jA<nA)
        jA = jA + 1; indA = get_na(stA,jA)
        jkf = 0
        do while (overlap0==0d0 .and. jkf<nqmf)
           jkf = jkf + 1; indf = npk_nchf + jkf - 1
           overlap0 = ortchil(indA,indf)
        end do
     end do
     if (overlap0==0d0) cycle        

     do configf = 1, nconfigf
        orbC = get_na(stf,configf,1); stC => TargetStates%b(orbC)
        orbD = get_na(stf,configf,2); stD => TargetStates%b(orbD)
        mC = dble( get_ang_mom_proj(stC) ); nC = get_nam(stC)
        mD = dble( get_ang_mom_proj(stD) ); nD = get_nam(stD)
        CIf = get_CI(stf,configf)

        par = pari * get_par(stC)
        mu = mi - mC; m = nint(mu)
        if (data_in%good_parity .and. get_par(stD)*get_par(stB) /= par) cycle   ! Conservation of parity.
        if (mD - mB /= mu) cycle   ! Conservation of orb. ang. mom. proj.

        lDMax = get_ang_mom( bst_nr%b(get_na(stD,nD)) )
        lBMax = get_ang_mom( bst_nr%b(get_na(stB,nB)) )
        lamMax = lDMax + lBMax + 2
        lamMin = abs(m) + modulo(m+lamMax,2)
        nlam = (lamMax - lamMin)/2 + 1

        do jlam = 1, nlam
           l = lamMin + 2*(jlam-1); lam = dble(l)
           formVec(:) = 0d0

           d1 = nr; d2 = 1
           do jD = 1, nD
              indD = get_na(stD,jD); spD => bst_nr%b(indD)
              CID = get_CI(stD,jD); lD = dble( get_ang_mom(spD) )
              vecD => fpointer(spD)
              a1 = get_minf(spD); a2 = get_maxf(spD)
              do jB = 1, nB
                 indB = get_na(stB,jB); spB => bst_nr%b(indB)
                 CIB = get_CI(stB,jB); lB = dble( get_ang_mom(spB) )
                 vecB => fpointer(spB)
                 b1 = max(a1, get_minf(spB)); b2 = min(a2, get_maxf(spB))
                 if (b1 >= b2) cycle

                 call angular_coeffs_3(lD,mD, lam,mu, lB,mB, rhoCoeff,etaCoeff)
                 if (rhoCoeff==0d0 .and. etaCoeff==0d0) cycle
!                 angVec(b1:b2) = rhoCoeff*ppRVec(b1:b2) + etaCoeff*R*R

                 !LHS: replaced old spheroidal form subroutine with more accurate form routine.
                 !     - the more accurate form does not need integration weights in inVec 
                 inVec(b1:b2) = vecD(b1:b2) * vecB(b1:b2) * (rhoCoeff*ppRVec(b1:b2)+etaCoeff*R*R)
                 call form_spheroidal_accurate(l,abs(m), inVec,b1,b2,nr, outVec,c1,c2)
                 !inVec(b1:b2) = vecD(b1:b2) * vecB(b1:b2) * (rhoCoeff*ppRVec(b1:b2)+etaCoeff*R*R) * weightVec(b1:b2)
                 !call form_spheroidal(l,abs(m), inVec,b1,b2,nr, outVec,c1,c2)
                 if (c1 >= c2) cycle

                 formVec(c1:c2) = formVec(c1:c2) + outVec(c1:c2)
                 d1 = min(c1,d1); d2 = max(c2,d2)

              end do ! jB
           end do ! jD
           if (d1 >= d2) cycle

           do jC = 1, nC
              indC = get_na(stC,jC); spC => bst_nr%b(indC)
              CIC = get_CI(stC,jC); lC = dble( get_ang_mom(spC) )
              vecC => fpointer(spC)
              a1 = max(d1, get_minf(spC)); a2 = min(d2, get_maxf(spC))
              if (a1 >= a2) cycle

              coeff = -(2d0*lam+1d0) * CIf*CIi * CIC * CID*CIB

              rhoVeci(:) = 0d0; etaVeci(:) = 0d0
              do jli = 1, ni
                 li = dble( liVec(jli) )
                 call angular_coeffs_3(lC,mC, lam,-mu, li,mi, rhoCoeff,etaCoeff)
                 rhoVeci(jli) = rhoCoeff
                 etaVeci(jli) = etaCoeff
              end do

              do jki = 1, nqmi
                 indi = npk_nchi + jki - 1; wavei => cont_oid%b(indi)
                 veci => get_rad_vec(wavei)
                 b1 = max(a1, get_rad_min(wavei)); b2 = min(a2, get_rad_max(wavei))
                 if (b1 >= b2) cycle

                 rhoCoeff = 0d0; etaCoeff = 0d0
                 do jli = 1, get_num_terms(wavei)
                    Di = get_term_D(wavei,jli)
                    rhoCoeff = rhoCoeff + Di*rhoVeci(jli)
                    etaCoeff = etaCoeff + Di*etaVeci(jli)
                 end do
                 if (rhoCoeff==0d0 .and. etaCoeff==0d0) cycle
!                 angVec(b1:b2) = rhoCoeff*ppRVec(b1:b2) + etaCoeff*R*R

                 tmp = coeff * sum(vecC(b1:b2)*veci(b1:b2) * formVec(b1:b2) * (rhoCoeff*ppRVec(b1:b2)+etaCoeff*R*R) * weightVec(b1:b2))

                 do jA = 1, nA
                    indA = get_na(stA,jA); CIA = get_CI(stA,jA)
                    do jkf = 1, nqmf
                       indf = npk_nchf + jkf - 1
                       overlap0 = CIA * ortchil(indA,indf)
                       vmatt(jkf,jki) = vmatt(jkf,jki) + overlap0*tmp 
                    end do ! jkf
                 end do ! jA
              end do ! jki
           end do ! jC
        end do ! jlam
     end do ! configi
  end do ! configf


  ! Do the electron 2 Hamiltonian term.
  !   < kf lf mf, stf | -H_2 P_r0r1 | ki li mi, sti >
  ! This has the relatively simple terms
  !   < kf lf mf | alpha lA mA >         = ortchil(indA,indf)
  !   < gamma lC mC | ki li mi >         = ortchil(indC,indi)
  !   < delta lD mD | H_2 | beta lB mB > = e1me(orbD,orbB)

  do configf = 1, nconfigf
     orbC = get_na(stf,configf,1); stC => TargetStates%b(orbC)
     orbD = get_na(stf,configf,2)
     CIf = get_CI(stf,configf)
     do configi = 1, nconfigi
        orbA = get_na(sti,configi,1); stA => TargetStates%b(orbA)
        orbB = get_na(sti,configi,2)
        CIi = get_CI(sti,configi)

        tmp = -CIf*CIi * e1me(orbD,orbB)
        if (tmp == 0d0) cycle

        do jC = 1, get_nam(stC)
           indC = get_na(stC,jC); CIC = get_CI(stC,jC)
           do jki = 1, nqmi
              indi = npk_nchi + jki - 1
              overlap1 = ortchil(indC,indi)

              do jA = 1, get_nam(stA)
                 indA = get_na(stA,jA); CIA = get_CI(stA,jA)
                 do jkf = 1, nqmf
                    indf = npk_nchf + jkf - 1
                    overlap0 = ortchil(indA,indf)

                    vmatt(jkf,jki) = vmatt(jkf,jki) + overlap0*overlap1*tmp
                 end do ! jkf
              end do ! jA
           end do ! jki
        end do ! jC
     end do ! configi
  end do ! configf


  coeff = -Zproj * 2d0/pi
! NEED TO CHECK WHERE THIS FACTOR OF 2 COMES FROM???
  tmp = 2d0*coeff * (-1)**(nint(Spinf+Spini+1)) * dsqrt((2d0*Spinf+1d0)*(2d0*Spini+1d0)) * COF6J(0.5d0,0.5d0,Spinf,0.5d0,S,Spini)
  vmatt(:,:) = tmp * vmatt(:,:)

  if(save_P01) P01mat = tmp*P01mat

  ! Do the "projection" term for when theta is non-zero.
  !   < kf lf mf, stf | -2 E theta I_0 | ki li mi, sti >

  if (non_uniq .and. theta/=0d0 .and. Spinf==Spini) then
     do configf = 1, nconfigf
        indC = get_na(stf,configf,1); indD = get_na(stf,configf,2)
        CIf = get_CI(stf,configf)
        do configi = 1, nconfigi
           indA = get_na(sti,configi,1); indB = get_na(sti,configi,2)
           CIi = get_CI(sti,configi)

           if (is_core_MO(indB) == 0) cycle
           if ( is_core_MO(indB)==1 .and. is_core_MO(indC)==0 ) cycle

           tmp = -2d0 * coeff * Etot*theta * CIf*CIi * ovlpst(indC,indA)*ovlpst(indD,indB)
           if (tmp == 0d0) cycle

           do ind = 1, TargState_Orthonorm%Nmax
              if (is_core_MO(indB)==1 .and. is_core_POMO(ind)==0) cycle
             
              do jkf = 1, nqmf
                 indf = npk_nchf + jkf - 1
                 do jki = 1, nqmi
                    indi = npk_nchi + jki - 1

                    overlap0 = ovlp_OrtnrmMO_k(ind,indf) * ovlp_OrtnrmMO_k(ind,indi)
                    vmatt(jkf,jki) = vmatt(jkf,jki) + overlap0 * tmp

                  end do ! jki
              end do ! jkf
           end do ! ind
        end do ! configi
     end do ! configf
  end if


end subroutine vexch_2e_oid

         






subroutine vexch_2e_MOrep_old(Mtot,nchf,nchi,nstf,nsti,nqmf,nqmi,npk_nchf,npk_nchi,chil,Nmax_1e_orb,ncwaves,ortchil,flchil,nr,weight,gridr,Zproj,Etot,theta,Rd,rspin_tot,vmatt)

! sum C^{np}_(np1,np2)  C^{n}_(n1,n2) < np0, np1, np2 : Mtot S |  | n0, n1, n2 : Mtot S >
! nstf = sum C^{np}_(np1,np2) | np1 np2 : m_f pi_f sp
!!$ Potential terms:  -2(E(1-theta)-H)P_{01}-E.I_0.theta, H=H0+H1+H2+V01+V02+V12+Z^2/R


  use sturmian_class
  use target_states
  use one_electron_func
  use channels
  use input_data
  use distpot
  use vnc_module
  use pol_2el
  use vmat_exch_module

  implicit none
  integer, intent(in):: Mtot
  integer, intent(in):: nchf, nchi        ! Channel number
  integer, intent(in):: nstf, nsti        ! State for particular channel
  integer, intent(in):: nqmf, nqmi        ! Number of k-grid points for particular channel
  integer, intent(in):: npk_nchf,npk_nchi ! On-shell k-grid point for particular channel
  integer, intent(in):: Nmax_1e_orb       ! Total number of target states one-electron orbitals      
!  real*8, dimension(Nmax_1e_orb,Nmax_1e_orb), intent(in):: ham1el_1elorb ! <np|H|n>  for the atomic orbitals
  integer, intent(in):: ncwaves           ! Number of on-, off- and bound projectile function ncwaves=nkgmax    
  real*8, dimension(Nmax_1e_orb,ncwaves), intent(in):: ortchil, flchil  ! ortchil=  <np|kLM>, flchil= <np|V-z.Zasy/r-U_dis|kLM>
  integer, intent(in):: nr                ! Max r points
  real*8, dimension(nr), intent(in):: weight, gridr
  real*8, intent(in):: Zproj        ! Z projectile charge
  real*8, intent(in):: Etot, theta        ! For non-uniqueness 
  real*8, intent(in):: Rd                 ! Internuclear distance      
  real*8, intent(in):: rspin_tot  ! Total Spin of the system. For H2 s_i=0 <=> S=1/2, s_i=1 <=> S=1/2,3/2
  real*8, dimension(nqmf,nqmi), intent(inout):: vmatt

  type(basis_sturmian_nr):: chil

  integer:: ki, kf, kqi, kqf   ! Scroll and Index of projectile-waves
  real*8:: ki_en, kf_en ! Projectile-wave energy
  real*8, dimension(nr):: epot_i, fun1, fun2
! Configuration Loops
  integer::  nst_con,  nstp_con, Ncon, Npcon        ! 2e configuration loops
  integer::  n_orb, np_orb, n_orb_max, np_orb_max   ! loops nuse 
  integer:: use_norb, use_nporb                     ! index to 1e molecular orbital
  integer:: MO_0, MO1, MO2, MOp1, MOp2                ! index to 1e molecular orbital
  integer:: i0, i1, i2, ip1, ip2, nam0, nam1, nam2, namp1, namp2  ! Loop over 1e molecular orbitals atomic orbitals
  type(sturmian_nr), pointer:: tn0,tnp0,tn1, tn2, tnp1, tnp2  ! pointers to one-electron atomic orbitals
  integer:: ind0, ind1, ind2, indp2, indp1                      ! index to underlying one-electron atomic orbital
  real*8:: CI_stp, CI_st, CI_0, CI_1, CI_2, CI_1p, CI_2p, temp_CI ! 2e configuration and 1e molecular orbial coefficients
  real*8, pointer, dimension(:)::  fp0, f0, fp1, f1, fp2, f2  ! underlying atomic orbital radial functions
! Radial integrations
  integer:: minf0, maxf0, minfp0, maxfp0, minf1, maxf1, minfp1, maxfp1
  integer:: minf2, maxf2, minfp2, maxfp2
  integer:: ir1, ir2, jr1, jr2, or1, or2
  real*8:: pi, overlap2, temp_overlap2, temp, overlap1, theta_unique

  logical:: Logic_V12,  Logic_V02
  real*8:: Yint, COF6J
  real*8:: ang0, ang1, ang2                        ! Angular integrals for projectile and A.O.(1) 
  integer:: lambda, lambda_min, lambda_max, mu
  integer:: kLap0, kLa1, kMp0, kM1                      ! Projectile Quantum Numbers (QN)  
  integer:: Lap1,Mp1,Lap2,Mp2,La0,M0,La2,M2,Mp_st,M_st  ! A.O. QN and Molecular State M 
  integer:: Spin, Spinp, tot_Spin
  integer:: La1, M1             ! QN and index for  projection operator term E.I.theta
  real*8 :: Z1, Z2

  Z1 = data_in%Z1
  Z2 = data_in%Z2

  pi = acos(-1d0)

  ! Projectile QN
  kLap0 = Lp_ch(nchf)  ! <kpLM(r_0)| 
  kMp0 = Mp_ch(nchf)   ! <kpLM(r_0)| 
  kLa1 = Lp_ch(nchi)   ! |kLM(r_1)>
  kM1 = Mp_ch(nchi)    ! |kLM(r_1)>

  ! Final State number nstf
  Npcon =  TargetStates2el%b(nstf)%nam          ! Number of 2e configurations.        
  Mp_st = NINT(TargetStates2el%b(nstf)%M )      ! Molecular state Angular Projection
  Spinp = NINT(TargetStates2el%b(nstf)%spin) ! 2e Molecular State Spin
  np_orb_max = TargetStates2el%b(nstf)%nusemax  ! Number of orbitals used to describe this state

  ! Initial State number nsti
  Ncon = TargetStates2el%b(nsti)%nam
  M_st = NINT(TargetStates2el%b(nsti)%M)
  Spin = NINT(TargetStates2el%b(nsti)%spin)
  n_orb_max = TargetStates2el%b(nsti)%nusemax

  ! Conserved QN. Selection Rules. Parity should be taken care of in
  ! channels.f90
  if( Mp_st + kMp0  /=  M_st + kM1 ) return

!!$     Here do one electron matrix elements and V_01 term:
!!$     Cp*C*<np2|n2>(<kp0|n0><np1|k1>(E-e_k'-e_k-Z^2/R) 
!!$     - <kp0|V_0-z.Zasy/r-U_0|n0><np1|k1>
!!$     - <kp0|n0><np1|V_1-z.Zasy/r-U_1|k1> - <kp0,np1|V_01|n0,k1>)
!!$ FINAL (p) State Molecular Orbital (MO) COORDINATE 1
!  print*,nstf,nsti
!  print*,"<MOp1, MOp2|...|MO_0, MO2>, theta_unique"
  do np_orb = 1, np_orb_max
!!$ nuse defined in state_class will refer to one-electron target states
!!$ so can be used! 
     use_nporb = TargetStates2el%b(nstf)%nuse(np_orb)

!!$ INITIAL State MO COORDINATE 0 
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
           namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
           namp2 = get_nam(TargetStates%b(MOp2))
!!$ Looping over INITIAL Molecular State configurations 
           do nst_con =  1, Ncon
              MO_0 = get_na(TargetStates2el%b(nsti),nst_con,1)
              if ( MO_0 /= use_norb ) cycle
              CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
              MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)
              M2 = get_ang_mom_proj(TargetStates%b(MO2))
              nam0 = get_nam(TargetStates%b(MO_0))
              nam2 = get_nam(TargetStates%b(MO2))

              overlap2 = 0d0

!!$ Solve for non-uniqueness if core orbital
              theta_unique = 1d0
              if (is_core_MO(MO2) == 2) theta_unique = (1d0 - theta)
              if (is_core_MO(MO2) == 1 .AND. is_core_MO(MOp1) /= 0 ) theta_unique = (1d0 - theta)
!              write(*,'(4I4,F10.5)') MOp1,MOp2,MO_0,MO2,theta_unique

              ! overlap 2 matrix element <varphi'_2|varphi_2> 
              if ( Mp2 /= M2 ) cycle
!!$ FINAL STATE COORDINATE 2 for BETA
              do ip2 = 1, namp2
                 indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
                 tnp2 => bst_nr%b(indp2)
                 fp2 => fpointer(tnp2)
                 Lap2 = get_ang_mom(tnp2)
                 CI_2p = get_CI(TargetStates%b(MOp2),ip2)
                 minfp2 = get_minf(tnp2)
                 maxfp2 = get_maxf(tnp2)
!!$ INITIAL STATE COORDINATE 2 for DELTA                
                 do i2 = 1, nam2
                    ind2 = get_na(TargetStates%b(MO2),i2)
                    tn2 => bst_nr%b(ind2)
                    f2 => fpointer(tn2)
                    La2 = get_ang_mom(tn2)
                    minf2 = get_minf(tn2)
                    maxf2 = get_maxf(tn2)
                    CI_2 = get_CI(TargetStates%b(MO2),i2)
                    ir1 = max(minf2,minfp2)
                    ir2 = min(maxf2,maxfp2)
                    if(ir2 < ir1) cycle
                    ! overlap 2 matrix element <varphi'_2|varphi_2> 
                    if ( Lap2 /= La2) cycle
                    temp_CI = CI_stp * CI_2p * CI_st * CI_2
                    if ( temp_CI == 0d0 ) cycle
!                    temp_overlap2 = SUM( fp2(ir1:ir2) *  f2(ir1:ir2) * weight(ir1:ir2))
!!$ MARK: Use the below but needs testing and checking MARK
                    temp_overlap2 = bst_nr%ortint(indp2,ind2)
                    overlap2 = overlap2 + temp_CI * temp_overlap2
                 end do ! i2        
              end do ! ip2   
!!$ Mark: Could code ovlpst for 1e Target States/Molecular orbitals
!!$              overlap2 = TargetStates%ortint ovlpst
!!$ This has no been done however,

!!$     Here do one electron matrix elements and V_01 term:
!!$     Cp*C*<np2|n2>(<kp0|n0><np1|k1>(E-e_k'-e_k-Z^2/R) 
!!$     - <kp0|V_0-z.Zasy/r-U_0|n0><np1|k1>
!!$     - <kp0|n0><np1|V_1-z.Zasy/r-U_1|k1> - <kp0,np1|V_01|n0,k1>)
              if ( overlap2 == 0d0 ) cycle
!!$ MO_0 = get_na(TargetStates2el%b(nsit),nst_con,1) = use_norb 
!!$ index 1e molecular orbital 
              M0 = get_ang_mom_proj(TargetStates%b(MO_0))
              Mp1 = get_ang_mom_proj(TargetStates%b(MOp1))
!!$ Quantum numbers and functions for ALPHA             
              do ip1 = 1, namp1
                 indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
                 tnp1 => bst_nr%b(indp1)   !           
                 fp1 => fpointer(tnp1)     ! One electron functions
                 Lap1 = get_ang_mom(tnp1)  ! Gets Angular momentum A.O.
                 minfp1 = get_minf(tnp1)
                 maxfp1 = get_maxf(tnp1)
                 CI_1p = get_CI(TargetStates%b(MOp1),ip1)
!!$ Calculating -Cp*C*<np2|n2><kp0,np1|V_01|n0,k1>
!!!$ Inital COORDINATE 1 |kLM(r_1)>
                 do ki = 1, nqmi
                    kqi = npk_nchi + ki - 1
                    tn1 => chil%b(kqi)
                    ki_en = get_energy(tn1)
                    f1 => fpointer(tn1)
                    minf1 = get_minf(tn1)
                    maxf1 = get_maxf(tn1)
                    or1 = max(minf1,minfp1)
                    or2 = min(maxf1,maxfp1)
                    
                    !LHS: removed weight so we can use form_accurate below
!                    fun1(or1:or2) = ( fp1(or1:or2) * f1(or1:or2) * weight(or1:or2) )
                    fun1(or1:or2) = ( fp1(or1:or2) * f1(or1:or2))
!!$ Quantum numbers and functions for GAMMA                  
                    do i0 = 1, nam0
                       ind0 = get_na(TargetStates%b(MO_0),i0) ! Index to underlying atomic orbitals  
                       tn0 => bst_nr%b(ind0)
                       f0 => fpointer(tn0)
                       La0 = get_ang_mom(tn0)
                       minf0 = get_minf(tn0)
                       maxf0 = get_maxf(tn0)
                       CI_0 = get_CI(TargetStates%b(MO_0),i0)
                       temp_CI = CI_1p * CI_0
                       if ( temp_CI == 0d0 ) cycle
                       lambda_min = max(ABS(kLap0 - La0),ABS(Lap1 - kLa1))
                       lambda_max = min(kLap0 + La0, Lap1 + kLa1)
                       ! mu should be the same for both integrals in coordinate 0 and 1
                       mu = Mp1 - kM1
                       lambda_max = min(data_in%ltmax,lambda_max)
                       do lambda = lambda_min, lambda_max
                          ang0 = 0d0
                          ang1 = 0d0
                          ang0 = dble((-1)**(mu))*Yint(dble(kLap0),dble(kMp0),dble(lambda),dble(-mu),dble(La0),dble(M0))
                          ang1 = Yint(dble(Lap1),dble(Mp1),dble(lambda),dble(mu),dble(kLa1),dble(kM1))
                          if ( ang0 * ang1 == 0d0 ) cycle
                          ! Integrate over coordinate 1 
                          call form_accurate(lambda,fun1,or1,or2,nr,epot_i,jr1,jr2)
!!$ COORDINATE 0 <kpLM(r_0)|
                          do kf = 1, nqmf
                             kqf = npk_nchf + kf - 1
                             tnp0 => chil%b(kqf)
                             fp0 => fpointer(tnp0)
                             minfp0 = get_minf(tnp0)
                             maxfp0 = get_maxf(tnp0)
                             ir1 = max(minf0,minfp0)
                             ir2 = min(maxf0,maxfp0)
                             ir1 = max(jr1,ir1)
                             ir2 = min(jr2,ir2)
                             if(ir2<ir1) cycle
                             ! Integrate over coordinate 0 
                             temp = SUM(fp0(ir1:ir2) * epot_i(ir1:ir2) * f0(ir1:ir2) * weight(ir1:ir2))
                             ! overlap2 has CI coefficients
                             temp = temp * ang0 * ang1 * overlap2
                             vmatt(kf,ki) = vmatt(kf,ki) - temp_CI * temp
                          end do ! end kf loop
                       end do ! lambda

!!$ Calculating Cp*C*<np2|n2>(<kp0|n0><np1|k1>(E(1-theta)-e_k'-e_k-Z^2/R) 
!!$ -<kp0|V_0-z.Zasy/r-U_0|n0><np1|k1> - <kp0|n0><np1|V_1-z.Zasy/r-U_1|k1> ) 
                       do kf = 1, nqmf
                          kqf = npk_nchf + kf - 1
                          tnp0 => chil%b(kqf)
                          kf_en = get_energy(tnp0)

                          temp = ortchil(ind0,kqf) * ortchil(indp1,kqi) * (Etot*(theta_unique)-kf_en-ki_en)
                          vmatt(kf,ki) = vmatt(kf,ki) + temp_CI * overlap2 * temp
                          if ( Rd /= 0d0 ) then
                             vmatt(kf,ki) = vmatt(kf,ki) - temp_CI * Z1*Z2/Rd * ortchil(ind0,kqf) * ortchil(indp1,kqi) * overlap2
                          end if
                          vmatt(kf,ki) = vmatt(kf,ki) - temp_CI * flchil(ind0,kqf) * ortchil(indp1,kqi) * overlap2
                          vmatt(kf,ki) = vmatt(kf,ki) - temp_CI * ortchil(ind0,kqf) * flchil(indp1,kqi) * overlap2
                       end do  ! end kf loop
                    end do ! i0
                 end do ! ki loop
              end do ! ip1
           end do ! nst_con - 2e config
        end do ! nstp_con
     end do ! n_orb
  end do ! np_orb

!!$ Calculating -Cp*C*<np1|k1><kp0,np2|V_02|n0,n2>
!!$ Looping over FINAL Molecular State orbitals
!!$ Looping over FINAL (p) Molecular State configorations 
  do nstp_con =  1, Npcon
     CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con) ! 2e config CI
     MOp1 = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
     Mp1 = get_ang_mom_proj(TargetStates%b(MOp1))      ! Angular momentum projection 1e molecular orbital
     MOp2 = get_na(TargetStates2el%b(nstf),nstp_con,2) ! Index to 1e molecular orbital
     Mp2 = get_ang_mom_proj(TargetStates%b(MOp2))      ! Angular momentum projection 1e molecular orbital
     namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
     namp2 = get_nam(TargetStates%b(MOp2))

     ! Check that coordinate 1 overlaps are non-zero
     Logic_V02 = .FALSE.
     do ki = 1, nqmi
        kqi = npk_nchi + ki - 1
        do ip1 = 1, namp1
           indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
           if ( ortchil(indp1,kqi) /= 0d0 ) then
              Logic_V02=.TRUE.
              exit
           end if
        end do ! ip1
        if (Logic_V02) exit
     end do ! ki
     if ( .NOT. Logic_V02 ) cycle
!!$ Looping over INITIAL Molecular State configurations 
     do nst_con =  1, Ncon
        CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
        MO_0 = get_na(TargetStates2el%b(nsti),nst_con,1)
        M0 = get_ang_mom_proj(TargetStates%b(MO_0))
        MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)
        M2 = get_ang_mom_proj(TargetStates%b(MO2))
        nam0 = get_nam(TargetStates%b(MO_0))
        nam2 = get_nam(TargetStates%b(MO2))
!!$ FINAL STATE COORDINATE 2 for BETA      
        do ip2 = 1, namp2
           indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
           tnp2 => bst_nr%b(indp2)
           fp2 => fpointer(tnp2)
           Lap2 = get_ang_mom(tnp2)
           CI_2p = get_CI(TargetStates%b(MOp2),ip2)
           minfp2 = get_minf(tnp2)
           maxfp2 = get_maxf(tnp2)
!!$ INITIAL STATE COORDINATE 2 for DELTA          
           do i2 = 1, nam2
              ind2 = get_na(TargetStates%b(MO2),i2)
              tn2 => bst_nr%b(ind2)
              f2 => fpointer(tn2)
              La2 = get_ang_mom(tn2)
              CI_2 = get_CI(TargetStates%b(MO2),i2)
              minf2 = get_minf(tn2)
              maxf2 = get_maxf(tn2)
              or1 = max(minfp2,minf2)
              or2 = min(maxfp2,maxf2)
              
              !LHS: removed weight so we can use form_accurate below             
!              fun2(or1:or2) = fp2(or1:or2) * f2(or1:or2) * weight(or1:or2)
              fun2(or1:or2) = fp2(or1:or2) * f2(or1:or2)
!!$ Mark: This loop can be moved inside the kf routine. min and max lambda
!!$ can be determined through an io loop here. Then the integration
!!$ can be done inside the kf loop after form. Lambda can be checked there. 
!!$  INITIAL STATE COORDINATE 0 for GAMMA
              do i0 = 1, nam0
                 ind0 = get_na(TargetStates%b(MO_0),i0) ! Index to underlying atomic orbitals  
                 tn0 => bst_nr%b(ind0)
                 f0 => fpointer(tn0)
                 La0 = get_ang_mom(tn0)
                 minf0 = get_minf(tn0)
                 maxf0 = get_maxf(tn0)
                 CI_0 = get_CI(TargetStates%b(MO_0),i0)
                 temp_CI = CI_stp * CI_2p * CI_st * CI_0 * CI_2
                 if ( temp_CI == 0d0 ) cycle

                 lambda_min = max(ABS(kLap0 - La0),ABS(Lap2 - La2))
                 lambda_max = min(kLap0 + La0, Lap2 + La2)
                 mu = Mp2 - M2
                 lambda_max = min(data_in%ltmax,lambda_max)
                 do lambda = lambda_min, lambda_max
                    ang0 = 0d0
                    ang2 = 0d0
                    ang0 = dble((-1)**(mu))*Yint(dble(kLap0),dble(kMp0),dble(lambda),dble(-mu),dble(La0),dble(M0))
                    ang2 = Yint(dble(Lap2),dble(Mp2),dble(lambda),dble(mu),dble(La2),dble(M2))
                    if ( ang0 * ang2 == 0d0 ) cycle
                    ! Integrate over coordinate 2 
                    call form_accurate(lambda,fun2,or1,or2,nr,epot_i,jr1,jr2)
!!$ COORDINATE 0 <kpLM(r_0)|
                    do kf = 1, nqmf
                       kqf = npk_nchf + kf - 1
                       tnp0 => chil%b(kqf)
                       fp0 => fpointer(tnp0)
                       minfp0 = get_minf(tnp0)
                       maxfp0 = get_maxf(tnp0)
                       ir1 = max(minfp0,minf0)
                       ir2 = min(maxfp0,maxf0)
                       ir1 = max(jr1,ir1)
                       ir2 = min(jr2,ir2)
                       if(ir2<ir1) cycle
!!$ Integrate over coordinate 0, temp = <kf|(CI_stp.CI_2p.CI_st.CI_0.CI_2.<varphi'_2|V_02|varphi_2>)|varphi_0> 
                       temp = SUM(fp0(ir1:ir2) * epot_i(ir1:ir2) * f0(ir1:ir2) * weight(ir1:ir2))
                       temp = temp_CI * temp * ang0 * ang2
!!$ Overlap coordinate 1
                       do ki = 1, nqmi
                          kqi = npk_nchi + ki - 1
                          do ip1 = 1, namp1
                             indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
                             CI_1p = get_CI(TargetStates%b(MOp1),ip1)
!!$ CI_p.<varphi'_1|ki>.(<kf|(CI_stp.CI_2p.CI_st.CI_0.CI_2.<varphi'_2|V_02|varphi_2>)|varphi_0>)                        
                             vmatt(kf,ki)= vmatt(kf,ki) - CI_1p * temp * ortchil(indp1,kqi)
                          end do ! ip1
                       end do ! ki cord 1
                    end do ! kf cord 0
                 end do ! lambda 
              end do ! i0 
           end do ! i2 
        end do !ip2
     end do ! nst_con
  end do ! nstp_con

!!$ Calculating -Cp*C*<kp0|n0><np1,np2|V_12|k1,n2>
!!$ Looping over FINAL Molecular State orbitals
  do nstp_con =  1, Npcon
     CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con) ! 2e config CI
     MOp1 = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
     Mp1 = get_ang_mom_proj(TargetStates%b(MOp1))      ! Angular momentum projection 1e molecular orbital
     MOp2 = get_na(TargetStates2el%b(nstf),nstp_con,2) ! Index to 1e molecular orbital
     Mp2 = get_ang_mom_proj(TargetStates%b(MOp2))      ! Angular momentum projection 1e molecular orbital
     namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
     namp2 = get_nam(TargetStates%b(MOp2))
!!$ Looping over INITIAL Molecular State configurations 
     do nst_con =  1, Ncon
        CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
        MO_0 = get_na(TargetStates2el%b(nsti),nst_con,1)
        M0 = get_ang_mom_proj(TargetStates%b(MO_0))
        MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)
        M2 = get_ang_mom_proj(TargetStates%b(MO2))
        nam0 = get_nam(TargetStates%b(MO_0))
        nam2 = get_nam(TargetStates%b(MO2))
!!$ Check that coordinate 0 overlaps are non-zero
        Logic_V12 = .FALSE.
        do kf = 1, nqmf
           kqf = npk_nchf + kf - 1
           do i0 = 1, nam0
              ind0 = get_na(TargetStates%b(MO_0),i0) ! Index to underlying atomic orbitals 
              if ( ortchil(ind0,kqf) /= 0d0 ) then
                 Logic_V12=.TRUE.
                 exit
              end if
           end do ! i0
           if ( Logic_V12 ) exit
        end do ! kf
        if ( .NOT. Logic_V12 ) cycle
!!$ FINAL STATE COORDINATE 2 BETA
        do ip2 = 1, namp2
           indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
           tnp2 => bst_nr%b(indp2)
           fp2 => fpointer(tnp2)
           Lap2 = get_ang_mom(tnp2)
           CI_2p = get_CI(TargetStates%b(MOp2),ip2)
           minfp2 = get_minf(tnp2)
           maxfp2 = get_maxf(tnp2)
!!$ INITIAL STATE COORDINATE 2 DELTA
           do i2 = 1, nam2
              ind2 = get_na(TargetStates%b(MO2),i2)
              tn2 => bst_nr%b(ind2)
              f2 => fpointer(tn2)
              La2 = get_ang_mom(tn2)
              minf2 = get_minf(tn2)
              maxf2 = get_maxf(tn2)
              CI_2 = get_CI(TargetStates%b(MO2),i2)
              or1 = max(minfp2,minf2)
              or2 = min(maxfp2,maxf2)
              
              !LHS: removed weight so we can use form_accurate below
              !fun2(or1:or2) = fp2(or1:or2) * f2(or1:or2) * weight(or1:or2)
              fun2(or1:or2) = fp2(or1:or2) * f2(or1:or2)
!!$ FINAL STATE COORDINATE 1 ALPHA
              do ip1 = 1, namp1
                 indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
                 tnp1 => bst_nr%b(indp1)   !           
                 fp1 => fpointer(tnp1)     ! One electron functions
                 Lap1 = get_ang_mom(tnp1)  ! Gets Angular momentum A.O.
                 minfp1 = get_minf(tnp1)
                 maxfp1 = get_maxf(tnp1)
                 CI_1p = get_CI(TargetStates%b(MOp1),ip1)
                 temp_CI = CI_stp * CI_1p * CI_2p * CI_st * CI_2
                 if ( temp_CI == 0d0 ) cycle

                 lambda_min = max(ABS(Lap1 - kLa1), ABS(Lap2 - La2))
                 lambda_max = min(Lap1 + kLa1, Lap2 + La2)
                 mu = Mp2 - M2
                 lambda_max = min(data_in%ltmax,lambda_max)
                 do lambda = lambda_min, lambda_max
                    ang1 = 0d0
                    ang2 = 0d0
                    ang1 = dble((-1)**(mu))*Yint(dble(Lap1),dble(Mp1),dble(lambda),dble(-mu),dble(kLa1),dble(kM1))
                    ang2 = Yint(dble(Lap2),dble(Mp2),dble(lambda),dble(mu),dble(La2),dble(M2))
                    if ( ang1 * ang2 == 0d0 ) cycle
!!$ Integrate over coordinate 2 
                    call form_accurate(lambda,fun2,or1,or2,nr,epot_i,jr1,jr2)
!!$ Initial COORDINATE 1 |kLM(r_1)>           
                    do ki = 1, nqmi
                       kqi = npk_nchi + ki - 1
                       tn1 => chil%b(kqi)
                       f1 => fpointer(tn1)
                       minf1 = get_minf(tn1)
                       maxf1 = get_maxf(tn1)
                       ir1 = max(minfp1,minf1)
                       ir2 = min(maxfp1,maxf1)
                       ir1 = max(jr1,ir1)
                       ir2 = min(jr2,ir2)
                       if(ir2<ir1) cycle
!!$ Integrate over coordinate 1 
                       temp = SUM(fp1(ir1:ir2) * epot_i(ir1:ir2) * f1(ir1:ir2) * weight(ir1:ir2))
                       temp = temp_CI * temp * ang1 * ang2
!!$ Overlap coordinate 0 
                       do kf = 1, nqmf
                          kqf = npk_nchf + kf - 1
                          do i0 = 1, nam0
                             ind0 = get_na(TargetStates%b(MO_0),i0) ! Index to underlying atomic orbitals  
                             CI_0 = get_CI(TargetStates%b(MO_0),i0)
                             vmatt(kf,ki)= vmatt(kf,ki) - CI_0 * temp * ortchil(ind0,kqf)
                          end do ! i0
                       end do ! kf
                    end do ! ki cord 1        
                 end do ! lambda 
              end do ! ip1
           end do ! i2 
        end do ! ip2
     end do ! nst_con
  end do ! nstp_con

!!$ Calculating -Cp*C*<kp0|n0><np1|k1><np2|H_2|n2>
!!$ Looping over FINAL Molecular State orbitals
  do nstp_con =  1, Npcon
     CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con) ! 2e config CI
     MOp1 = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
     MOp2 = get_na(TargetStates2el%b(nstf),nstp_con,2) ! Index to 1e molecular orbital
     namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
     namp2 = get_nam(TargetStates%b(MOp2))
!!$ Looping over INITIAL Molecular State configurations 
     do nst_con =  1, Ncon
        CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
        MO_0 = get_na(TargetStates2el%b(nsti),nst_con,1)
        MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)
        nam0 = get_nam(TargetStates%b(MO_0))
        nam2 = get_nam(TargetStates%b(MO2))
!!$ FINAL STATE COORDINATE 2 for BETA      
        do ip2 = 1, namp2
           indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
           CI_2p = get_CI(TargetStates%b(MOp2),ip2)
!!$ INITIAL STATE COORDINATE 2 for DELTA          
           do i2 = 1, nam2
              ind2 = get_na(TargetStates%b(MO2),i2)
              CI_2 = get_CI(TargetStates%b(MO2),i2)
              if ( bst_nr%ham1el(indp2,ind2) == 0d0 ) cycle
!!$ Overlap coordinate 1 <np1|k1>
              do ki = 1, nqmi
                 kqi = npk_nchi + ki - 1
                 do ip1 = 1, namp1
                    indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
                    CI_1p = get_CI(TargetStates%b(MOp1),ip1)
                    if (ortchil(indp1,kqi) == 0d0) cycle
!!$ Overlap coordinate 0 <kp0|n0>
                    do kf = 1, nqmf
                       kqf = npk_nchf + kf - 1
                       do i0 = 1, nam0
                          ind0 = get_na(TargetStates%b(MO_0),i0) ! Index to underlying atomic orbitals  
                          CI_0 = get_CI(TargetStates%b(MO_0),i0)
                          if (ortchil(ind0,kqf) == 0d0 ) cycle
                          temp_CI = CI_stp * CI_1p * CI_2p * CI_st * CI_0 * CI_2
                          if ( temp_CI == 0d0 ) cycle
                          temp = temp_CI * ortchil(ind0,kqf) * ortchil(indp1,kqi)
                          vmatt(kf,ki)= vmatt(kf,ki) - temp * bst_nr%ham1el(indp2,ind2)
                       end do ! i0
                    end do ! kf coord 0
                 end do ! ip1
              end do ! ki coord 1
           end do ! i2 
        end do ! ip2
     end do ! nst_con
  end do ! nstp_con


!!$ Spin and other coefficients
!!$ Need to be done before projection operator overlaps
!!$ + 2. Spin_coeff.(E(1-theta)-H)P_01, vmatt=(E(1-theta)-H)P_01 
  temp = sqrt(dble(2*Spinp+1)) * sqrt(dble(2*Spin+1))
  temp = COF6J(0.5d0,0.5d0,dble(Spinp),0.5d0,dble(rspin_tot),dble(Spin)) * temp
  temp = 2d0 * dble((-1)**(Spinp+Spin+1)) * temp

  vmatt(:,:) = temp * vmatt(:,:)

!!$ Projection operator part of non-uniqueness
!!$ 2(E-theta+theta)P_01
!!$ -2 E * theta * <Phi_f|Phi_i> * <kp| I_0 |k>
!  if (nsti == nstf .AND. theta /= 0d0 ) then
  if ( theta /= 0d0 .AND.  Spinp == Spin ) then
!!$ Looping over FINAL (p) Molecular State configorations 
     do nstp_con =  1, Npcon
        CI_stp = get_CI(TargetStates2el%b(nstf),nstp_con) ! 2e config CI
        MOp1 = get_na(TargetStates2el%b(nstf),nstp_con,1) ! Index 1e molecular orbital
        Mp1 = get_ang_mom_proj(TargetStates%b(MOp1))      ! Angular momentum projection 1e molecular orbital
        MOp2 = get_na(TargetStates2el%b(nstf),nstp_con,2) ! Index to 1e molecular orbital
        Mp2 = get_ang_mom_proj(TargetStates%b(MOp2))      ! Angular momentum projection 1e molecular orbital
        namp1 = get_nam(TargetStates%b(MOp1))             ! Number of atomic orbitals that make up molecular orbital 
        namp2 = get_nam(TargetStates%b(MOp2))

!!$ Looping over INITIAL Molecular State configurations 
        do nst_con =  1, Ncon
           CI_st = get_CI(TargetStates2el%b(nsti),nst_con)
           MO1 = get_na(TargetStates2el%b(nsti),nst_con,1)
           M1 = get_ang_mom_proj(TargetStates%b(MO1))
           MO2 = get_na(TargetStates2el%b(nsti),nst_con,2)
           M2 = get_ang_mom_proj(TargetStates%b(MO2))
           nam1 = get_nam(TargetStates%b(MO1))
           nam2 = get_nam(TargetStates%b(MO2))
           if ( Mp1 /= M1 .OR. Mp2 /= M2 ) cycle
           if ( is_core_MO(MO2) == 0 ) then ! cycle
!                print*,"Proj M02 delta",MO2
                cycle
           end if
           if ( is_core_MO(MO2) == 1 .AND. is_core_MO(MOp1) == 0 ) then ! cycle
!                print*,"Proj M02 delta, MOp1",MO2, MOp1
                cycle
           end if
!!$ overlap 2 matrix element CI_stp.CI_2p.CI_st.CI_2.<np2|n2> 
           do ip2 = 1, namp2
              indp2 = get_na(TargetStates%b(MOp2),ip2) ! Index to underlying atomic orbitals 
              tnp2 => bst_nr%b(indp2)
              fp2 => fpointer(tnp2)
              Lap2 = get_ang_mom(tnp2)
              CI_2p = get_CI(TargetStates%b(MOp2),ip2)
              minfp2 = get_minf(tnp2)
              maxfp2 = get_maxf(tnp2)
              do i2 = 1, nam2
                 ind2 = get_na(TargetStates%b(MO2),i2)
                 tn2 => bst_nr%b(ind2)
                 f2 => fpointer(tn2)
                 La2 = get_ang_mom(tn2)
                 minf2 = get_minf(tn2)
                 maxf2 = get_maxf(tn2)
                 CI_2 = get_CI(TargetStates%b(MO2),i2)
                 if (Lap2 /= La2 ) cycle
                 temp_CI = CI_stp * CI_2p * CI_st * CI_2
                 if ( temp_CI == 0d0 ) cycle
                 ir1 = max(minf2,minfp2)
                 ir2 = min(maxf2,maxfp2)
                 if(ir2<ir1) cycle
!                 overlap2 = temp_CI * SUM( fp2(ir1:ir2) *  f2(ir1:ir2) * weight(ir1:ir2))
!!$ MARK: Use the below but needs testing and checking MARK
                 overlap2 = temp_CI * bst_nr%ortint(indp2,ind2)

                 if ( overlap2 == 0d0 ) cycle
!!$ overlap1=CI_1p.CI_1.<np1|n1>
                 do ip1 = 1, namp1
                    indp1 = get_na(TargetStates%b(MOp1),ip1)  ! Index to underlying atomic orbitals 
                    tnp1 => bst_nr%b(indp1)   !           
                    fp1 => fpointer(tnp1)     ! One electron functions
                    Lap1 = get_ang_mom(tnp1)  ! Gets Angular momentum A.O.
                    minfp1 = get_minf(tnp1)
                    maxfp1 = get_maxf(tnp1)
                    CI_1p = get_CI(TargetStates%b(MOp1),ip1)
                    do i1 = 1, nam1
                       ind1 = get_na(TargetStates%b(MO1),i1)  ! Index to underlying atomic orbitals 
                       tn1 => bst_nr%b(ind1)   !           
                       f1 => fpointer(tn1)     ! One electron functions
                       La1 = get_ang_mom(tn1)  ! Gets Angular momentum A.O.
                       if ( Lap1 /= La1 ) cycle
                       minf1 = get_minf(tn1)
                       maxf1 = get_maxf(tn1)
                       CI_1 = get_CI(TargetStates%b(MO1),i1)
                       temp_CI = CI_1p * CI_1
                       if ( temp_CI == 0d0 ) cycle
                       ir1 = max(minf1,minfp1)
                       ir2 = min(maxf1,maxfp1)
                       if(ir2<ir1) cycle
!                       overlap1 = temp_CI * SUM( fp1(ir1:ir2) *  f1(ir1:ir2) * weight(ir1:ir2))
                       overlap1 = temp_CI * bst_nr%ortint(indp1,ind1)
                       if ( overlap1 == 0d0 ) cycle
                       do ind0 = 1, TargState_Orthonorm%Nmax  ! Loop over all one-electron target states 
!!$ Note the difference between is_core_MO and is_core_POMO
!!$ CORE ORBITALS CHECKED HERE FOR MULTI CORE CASE 
                          if (is_core_MO(MO2) == 1 .AND. is_core_POMO(ind0) == 0) cycle

                          do ki = 1, nqmi
                             kqi = npk_nchi + ki - 1
                             do kf = 1, nqmf
                                kqf = npk_nchf + kf - 1
                                temp = ovlp_OrtnrmMO_k(ind0,kqi) * ovlp_OrtnrmMO_k(ind0,kqf)
                                temp = temp * overlap1 * overlap2
                                vmatt(kf,ki) = vmatt(kf,ki) - 2d0 * temp * Etot * theta
                             end do ! kf
                          end do ! ki
!                          print*,"n, <n|k_1>",ind0,ovlp_OrtnrmMO_k(ind0,1)
                       end do ! ind0 I_0 
                    end do ! i1
                 end do ! ip1
              end do ! i2
           end do !ip2
        end do ! nst_con
     end do ! nstp_con
  end if ! non-uniqueness? and spin

  vmatt(:,:) = vmatt(:,:) * (2d0/pi) !!! / (rkf*rki) dealt with in scat.f 


end subroutine vexch_2e_MOrep_old

!Liam renamed this from vdirect_2e to vdirect_2e_sp since we always assume MO representation now
subroutine vdirect_2e_sp(Mtot,rspin_tot, nchf,nchi,nstf,nsti, nqmf,nqmi,npk_nchf,npk_nchi,chil, nr,weight,gridr,Zproj,Zasym,iBorn,vmatt)

! sum C^{np}_(np1,np2)  C^{n}_(n1,n2) < np0, np1, np2 : Mtot | V(0,1) | n0, n1, n2 : Mtot >
! Note: All integrals in target states coordinate space and projectile angular terms 
! multiply radial potentials are taken first. 
! This is because we need to eliminate the lambda = 0 term for the electron-electron 
! potential of the projectile and target states.
! V_{lambda=0}(0,1) = int (sin(kr))^2.phi_i.phi_f/r = oo. 
! Hence we eliminate the problem by using the projectile-nuclear term

  use sturmian_class
  use target_states
  use one_electron_func
  use channels
  use input_data
  use distpot
  use vnc_module
  use pol_2el

  implicit none
  integer, intent(in):: Mtot
  real*8, intent(in) :: rspin_tot
  integer, intent(in):: nchf, nchi        ! Channel number
  integer, intent(in):: nstf, nsti        ! State for particular channel
  integer, intent(in):: nqmf, nqmi        ! Number of k-grid points for particular channel
  integer, intent(in):: npk_nchf,npk_nchi ! On-shell k-grid point for particular channel
  integer, intent(in):: nr                ! Max r points
  real*8, dimension(nr), intent(in):: weight, gridr
  real*8, intent(in):: Zproj, Zasym      ! Z projectile charge, Z Aymptotic charge of target
  integer, intent(in):: iBorn             ! Born = 1
  real*8, dimension(nqmf,nqmi), intent(inout):: vmatt

  type(basis_sturmian_nr):: chil
  
  integer:: ki, kf, kqi, kqf, i
  real*8, dimension(nr):: ttt, temp, fun0, fun1
  ! Configuration Loops
  integer::  ne_con,  nep_con, Ncon, Npcon, n_orb, np_orb, n_orb_max, np_orb_max 
  integer:: nuse_norb, nuse_nporb                         
  type(sturmian_nr), pointer:: tn0,tnp0,tn1, tn2, tnp1, tnp2     ! One-electron orbitals
  integer:: ind1, ind2, indp1, indp2, indp1_2, ind1_2
  real*8, pointer, dimension(:)::  fp0, f0, fp1, f1, fp2, f2
  ! Radial integrations
  integer::  minf, maxf, minfp, maxfp,  ir1, ir2, irt1,irt2   
  integer::  ir2_nuc, ir2_e_e, ir_all
  real*8:: CIp, CI                                         ! Moleculart state Configuration coefficients
  real*8:: sum_rad, pi, overlap2,  temp_overlap2
  
  real*8:: Yint
  real*8:: ang0, ang1                                   ! Angular integrals for projectile and A.O.(1) 
  integer:: lambda, lambda_min, lambda_max, lambda_step, mu, par
  integer:: Lap0, La0, Mp0, M0                          ! Projectile Quantum Numbers (QN)  
  integer:: Lap1,Mp1,Lap2,Mp2,La1,M1,La2,M2,Mp12,M12    ! A.O. QN and Molecular State M 
  integer:: Spin, Spinp


  pi = acos(-1d0)
  ir2_e_e = 0
  ir2_nuc = 0

  ! Projectile QN
  Lap0 = Lp_ch(nchf)  
  Mp0 = Mp_ch(nchf)
  La0 = Lp_ch(nchi)  
  M0 = Mp_ch(nchi)

  ! Final State number nstf
  Npcon =  TargetStates2el%b(nstf)%nam          ! Number of A.O. and Molecular ion configurations.        
  Mp12 = NINT(TargetStates2el%b(nstf)%M )        ! Molecular state Angular Projection
  Spinp = NINT(TargetStates2el%b(nstf)%spin)    ! 2e Molecular State Spin
  np_orb_max = TargetStates2el%b(nstf)%nusemax  ! Number of orbitals used to describe this state

  ! Initial State number nsti
  Ncon = TargetStates2el%b(nsti)%nam         
  M12 = NINT(TargetStates2el%b(nsti)%M)       
  Spin = NINT(TargetStates2el%b(nsti)%spin)
  n_orb_max = TargetStates2el%b(nsti)%nusemax     


! Conserved QN. Selection Rules. Parity should be taken care of in channels.f90
  if (Spin /= Spinp) return   ! No direct interaction between states of different spins.
  if (rspin_tot==1.5d0 .and. Spin==0) return   ! Both states must be triplet for direct interaction in the S=3/2 partial wave.
  if ( Mp12 + Mp0  /=  M12 + M0 ) return

  ttt(:) = 0d0

  ! FINAL State COORDINATE 1
  do np_orb = 1, np_orb_max
   
     nuse_nporb = TargetStates2el%b(nstf)%nuse(np_orb)
     
     ! INITIAL State COORDINATE 1
     do n_orb = 1, n_orb_max
        
        nuse_norb =  TargetStates2el%b(nsti)%nuse(n_orb)   
        
        ! Below sums all the overlaps for COORDINATE 2  same configuratins(1s,2s,..) in COORDINATE 1
        overlap2 = 0d0
        
        ! Looping over FINAL Molecular State orbitals. COORDINATE 2
        do nep_con =  1, Npcon          
           
           indp1_2 = TargetStates2el%b(nstf)%na(nep_con)    ! Final state number np. nep A.O.
           if ( nuse_nporb /= indp1_2 ) cycle
           
           ! Quantum numbers and functions for BETA
           indp2 = TargetStates2el%b(nstf)%nb(nep_con)                
           tnp2 => bst_nr%b(indp2)                                         
           fp2 => fpointer(tnp2)                                
           Lap2 = get_ang_mom(tnp2)                             
           Mp2 = get_ang_mom_proj(tnp2)          
           CIp = get_CI(TargetStates2el%b(nstf),nep_con)  
           
           ! Looping over INITIAL Molecular State orbitals.  COORDINATE 2
           do ne_con = 1, Ncon       
              
              ind1_2 = TargetStates2el%b(nsti)%na(ne_con)  
              if ( nuse_norb /= ind1_2 ) cycle            

              ! Quantum numbers and functions for DELTA
              ind2 = TargetStates2el%b(nsti)%nb(ne_con) 
              tn2 => bst_nr%b(ind2)                                         
              f2 => fpointer(tn2)                               
              La2 = get_ang_mom(tn2)                            
              M2 = get_ang_mom_proj(tn2)            
              CI = get_CI(TargetStates2el%b(nsti),ne_con)  
              
              ! Selections Rules  
              if ( Lap2 /= La2 .OR. Mp2 /= M2  ) cycle
              
              temp_overlap2 = 0d0
              ! DO OVERLAP OF COORDINATE SPACE 2  
              minf = get_minf(tn2)
              maxf = get_maxf(tn2)
              minfp = get_minf(tnp2)
              maxfp = get_maxf(tnp2)      
              ir1 = max(minf,minfp)
              ir2 = min(maxf,maxfp)  
              if(ir2<ir1) cycle
!              temp_overlap2 = SUM( fp2(ir1:ir2) *  f2(ir1:ir2) * weight(ir1:ir2) )
              temp_overlap2 = bst_nr%ortint(indp2,ind2)
              
              overlap2 = overlap2 +  CIp * CI * temp_overlap2
              
           end do  ! INITIAL STATE COORDINATE 2           
        end do    ! FINAL STATE COORDINATE 2       
        if ( overlap2 == 0d0 ) cycle 


        ! COORDINATE 1 RADIAL INTEGRALS
        ! Quantum numbers and functions for ALPHA
        indp1 = nuse_nporb                                  ! Final state number np. nep A.O.
        tnp1 => bst_nr%b(indp1)                             !           
        fp1 => fpointer(tnp1)                               ! One electron functions
        Lap1 = get_ang_mom(tnp1)                            ! Gets Angular momentum A.O.
        Mp1 = get_ang_mom_proj(tnp1)                        ! Get angular projection of A.O. 
        
        
        ! Quantum numbers and functions for GAMMA
        ind1 = nuse_norb               
        tn1 => bst_nr%b(ind1)                                          
        f1 => fpointer(tn1)                               
        La1 = get_ang_mom(tn1)                            
        M1 = get_ang_mom_proj(tn1)  

        ! Electron-electron potential 2*V_{01}
        lambda_min = max(ABS(Lap0 - La0),ABS(Lap1-La1))
        lambda_max = min(Lap0 + La0, Lap1 + La1)
        ! mu should be the same for both integrals in coordinate 0 and 1, otherwise == 0.
        mu = Mp1 - M1
        if ( lambda_min > lambda_max ) cycle
              
        minf = get_minf(tn1)
        maxf = get_maxf(tn1)
        minfp = get_minf(tnp1)
        maxfp = get_maxf(tnp1)
        ir1 = max(minf,minfp)
        ir2 = min(maxf,maxfp)

        fun1(:) = 0d0
        ! Multiply the electron-electron potential by 2. Refer D. Fursa e-He 95 Paper
        
        !LHS: removed weight so we can use form_accurate below
        !fun1(ir1:ir2) = 2.0 * (f1(ir1:ir2) * fp1(ir1:ir2) * weight(ir1:ir2)) * overlap2
        fun1(ir1:ir2) = 2.0 * (f1(ir1:ir2) * fp1(ir1:ir2)) * overlap2

        lambda_max = min(data_in%ltmax,lambda_max)
        do lambda = lambda_min, lambda_max
           ang0 = 0d0
           ang1 = 0d0
           ang0 = dble((-1)**(mu))*Yint(dble(Lap0),dble(Mp0),dble(lambda),dble(-mu),dble(La0),dble(M0))
           ang1 = Yint(dble(Lap1),dble(Mp1),dble(lambda),dble(mu),dble(La1),dble(M1))

           if ( ang0 == 0d0 .OR. ang1 == 0d0) cycle
          
           call form_accurate(lambda, fun1, ir1, ir2, nr, temp, irt1, irt2)   

           ttt(irt1:irt2) = ttt(irt1:irt2) + temp(irt1:irt2) * ang0 * ang1
           ! Checking of angular terms non-zero
           ir2_e_e = max(ir2_e_e,irt2)

        end do ! lambda
     end do ! INITIAL STATE COORDINATE 1    
  end do  ! FINAL STATE COORDINATE 1


  ! Projectile-Nuclear Potential
  if (nsti == nstf) then

     lambda_min = ABS(Lap0 - La0)
     lambda_max = Lap0 + La0
     mu = Mp0 - M0

!!$ These lines are to make sure that pw Born and analytical Born produce the same values
!!$ Currently only lambda=0 (spherical) term is used in analytical Born part, do the same here
     if(iBorn .eq. 1) then   
!        lambda_max = 0
     endif
!
     if(data_in%good_parity) then 
       lambda_step = 2
     else
       lambda_step = 1
     endif

     do lambda = lambda_min, lambda_max, lambda_step

        if (lambda > lamtop_vc) exit

        ir1 = minvnc(lambda)
        ir2 = maxvnc(lambda)      
        
        ! Checking of angular terms non-zero
        ir2_nuc = max(ir2_nuc,ir2)
    
        if ( lambda == 0) then
           ttt(1:ir2) = ttt(1:ir2) + ( vnc(1:ir2,lambda) + Zasym/gridr(1:ir2) )
           ! Distorting potential -U_{0}
           if (La0<=data_in%Ldw .and. iBorn/=1) then
              if(iBorn .eq. 0) then   
                 ttt(1:maxvdist) =  ttt(1:maxvdist) + Zproj * vdist(1:maxvdist)
              endif
           end if
        else
           ang0 = 0d0
           ang0 = Yint(dble(Lap0),dble(Mp0),dble(lambda),dble(mu),dble(La0),dble(M0))         
           if ( ang0 == 0d0 ) cycle
           ttt(ir1:ir2) = ttt(ir1:ir2) + ang0 * vnc(ir1:ir2,lambda)
           
        end if
     end do
  end if

  if ( ir2_nuc == 0 .AND. ir2_e_e == 0 ) then
     return !  this ME is zero just due to angular momentum, it never got inside lambda loop
  end if

  call minmaxi(ttt,nr,ir1,ir2)
  ir_all = ir2
  if ( ir2 < nr ) then
     ttt(ir2+1:nr) = 0d0
  end if

  do ki=1,nqmi
     kqi = npk_nchi + ki - 1
     tn0 => chil%b(kqi)
     f0 => fpointer(tn0)
     minf = get_minf(tn0)
     maxf = get_maxf(tn0)
     
     do kf=1,nqmf

       !Liam added: V-matrix blocks with same final/initial channels are symmetric
        if(nchf == nchi .and. kf < ki) then
          vmatt(kf,ki) = vmatt(ki,kf)
          cycle
        endif

        kqf = npk_nchf + kf - 1
        tnp0 => chil%b(kqf)
        fp0 => fpointer(tnp0)
        minfp = get_minf(tnp0)
        maxfp = get_maxf(tnp0)
        
        ir2 = min(maxf,maxfp) ! nr as for continuum functions, but for closed channels maxf = 1
        ir1 = max(minf,minfp)
        
        ir2 = min(ir2,ir_all)

        fun0(ir1:ir2) = (f0(ir1:ir2) * fp0(ir1:ir2)) * weight(ir1:ir2) * ttt(ir1:ir2)  

        sum_rad = SUM( fun0(ir1:ir2) ) 

        vmatt(kf,ki) = vmatt(kf,ki) - Zproj * sum_rad * (2d0/pi)

     end do ! end kf loop

  end do ! end ki loop

end subroutine vdirect_2e_sp


!Liam renmed this from vexch_2e to vexch_2e_sp since we always use MO representation now
subroutine vexch_2e_sp(Mtot,nchf,nchi,nstf,nsti,nqmf,nqmi,npk_nchf,npk_nchi,chil,Nmax_1e_orb,ncwaves,ortchil,flchil,nr,weight,gridr,Zproj,Etot,theta,Rd,rspin_tot,vmatt)
  
  ! sum C^{np}_(np1,np2)  C^{n}_(n1,n2) < np0, np1, np2 : Mtot S |  | n0, n1, n2 : Mtot S >
  ! nstf = sum C^{np}_(np1,np2) | np1 np2 : m_f pi_f sp
!!$ Potential terms:  -2(E(1-theta)-H)P_{01}-E.I_0.theta, H=H0+H1+H2+V01+V02+V12+Z^2/R


!!$****** This subroutine contains some bugs and wrong logic. 
!!$ Follow the logic of the tried and tested MOrep subroutine
  
  use sturmian_class
  use target_states
  use one_electron_func
  use channels
  use input_data
  use distpot
  use vnc_module
  use pol_2el
  use vmat_exch_module
  
  implicit none
  integer, intent(in):: Mtot
  integer, intent(in):: nchf, nchi        ! Channel number
  integer, intent(in):: nstf, nsti        ! State for particular channel
  integer, intent(in):: nqmf, nqmi        ! Number of k-grid points for particular channel
  integer, intent(in):: npk_nchf,npk_nchi ! On-shell k-grid point for particular channel
  integer, intent(in):: Nmax_1e_orb       ! Total number of target states one-electron orbitals      
!  real*8, dimension(Nmax_1e_orb,Nmax_1e_orb), intent(in):: ham1el_1elorb ! <np|H|n>  for the atomic orbitals
  integer, intent(in):: ncwaves           ! Number of on-, off- and bound projectile function ncwaves=nkgmax    
  real*8, dimension(Nmax_1e_orb,ncwaves), intent(in):: ortchil, flchil  ! ortchil=  <np|kLM>, flchil= <np|V-z.Zasy/r-U_dis|kLM>
  integer, intent(in):: nr                ! Max r points
  real*8, dimension(nr), intent(in):: weight, gridr
  real*8, intent(in):: Zproj        ! Z projectile charge
  real*8, intent(in):: Etot, theta        ! For non-uniqueness 
  real*8, intent(in):: Rd                 ! Internuclear distance      
  real*8, intent(in):: rspin_tot  ! Total Spin of the system. For H2 s_i=0 <=> S=1/2, s_i=1 <=> S=1/2,3/2
  real*8, dimension(nqmf,nqmi), intent(inout):: vmatt
  
  
  type(basis_sturmian_nr):: chil
  
  integer:: ki, kf, kqi, kqf   ! Scroll and Index of projectile-waves
  real*8:: ki_en, kf_en ! Projectile-wave energy
  real*8, dimension(nr):: epot_i, fun1, fun2
  ! Configuration Loops
  integer::  ne_con,  nep_con, Ncon, Npcon, n_orb, np_orb, n_orb_max, np_orb_max 
  integer:: nuse_norb, nuse_nporb                         
  type(sturmian_nr), pointer:: tn0,tnp0,tn1, tn2, tnp1, tnp2  ! One-electron orbitals and projectile functions
  integer:: ind0, ind2, indp2, indp1                      ! index to one-electron orbital
  real*8, pointer, dimension(:)::  fp0, f0, fp1, f1, fp2, f2
  ! Radial integrations
  integer:: minf0, maxf0, minfp0, maxfp0, minf1, maxf1, minfp1, maxfp1
  integer:: minf2, maxf2, minfp2, maxfp2
  integer:: ir1, ir2, jr1, jr2, or1, or2   
  real*8:: CIp, CI                                 ! Moleculart state Configuration coefficients
  real*8:: pi, overlap2, temp_overlap2, temp
  real*8:: overlap12
  
  logical:: Logic_V12,  Logic_V02
  real*8:: Yint, COF6J!
  real*8:: ang0, ang1, ang2                        ! Angular integrals for projectile and A.O.(1) 
  integer:: lambda, lambda_min, lambda_max, lambda_step, mu
  integer:: kLap0, kLa1, kMp0, kM1                      ! Projectile Quantum Numbers (QN)  
  integer:: Lap1,Mp1,Lap2,Mp2,La0,M0,La2,M2,Mp_st,M_st  ! A.O. QN and Molecular State M 
  integer:: Spin, Spinp, tot_Spin
  integer:: ind1, La1, M1             ! QN and index for  projection operator term E.I.theta
  integer :: j  
  real*8 :: Z1, Z2

  real*8, dimension(nr) :: f0weight
  
  Z1 = data_in%Z1
  Z2 = data_in%Z2
  
  pi = acos(-1d0)
  vmatt(:,:) = 0d0
  
  ! Projectile QN
  kLap0 = Lp_ch(nchf)  ! <kpLM(r_0)| 
  kMp0 = Mp_ch(nchf)   ! <kpLM(r_0)| 
  kLa1 = Lp_ch(nchi)   ! |kLM(r_1)>
  kM1 = Mp_ch(nchi)    ! |kLM(r_1)>
  
  ! Final State number nstf
  Npcon =  TargetStates2el%b(nstf)%nam          ! Number of A.O. and Molecular ion configurations.        
  Mp_st = NINT(TargetStates2el%b(nstf)%M )      ! Molecular state Angular Projection
  Spinp = NINT(TargetStates2el%b(nstf)%spin) ! 2e Molecular State Spin
  np_orb_max = TargetStates2el%b(nstf)%nusemax  ! Number of orbitals used to describe this state
  
  ! Initial State number nsti
  Ncon = TargetStates2el%b(nsti)%nam         
  M_st = NINT(TargetStates2el%b(nsti)%M)       
  Spin = NINT(TargetStates2el%b(nsti)%spin)
  n_orb_max = TargetStates2el%b(nsti)%nusemax    
  
  ! Conserved QN. Selection Rules. Parity should be taken care of in channels.f90
  if( Mp_st + kMp0  /=  M_st + kM1 ) return
  
  ! FINAL State COORDINATE 1
  do np_orb = 1, np_orb_max
     nuse_nporb = TargetStates2el%b(nstf)%nuse(np_orb)
     
     ! INITIAL State COORDINATE 0 
     do n_orb = 1, n_orb_max

        nuse_norb =  TargetStates2el%b(nsti)%nuse(n_orb)   
        
        ! Below sums all the overlaps for COORDINATE 2  same
        ! configuratins(1s,2s,..) in COORDINATE 1
        overlap2 = 0d0
        
        ! Looping over FINAL Molecular State orbitals. COORDINATE 2
        do nep_con =  1, Npcon          
           
           indp1 = TargetStates2el%b(nstf)%na(nep_con)    ! Final state number np. nep A.O.
           
           if ( nuse_nporb /= indp1 ) cycle
           
           ! Quantum numbers and functions for BETA
           indp2 = TargetStates2el%b(nstf)%nb(nep_con)                
           tnp2 => bst_nr%b(indp2)                                         
           fp2 => fpointer(tnp2)                                
           Lap2 = get_ang_mom(tnp2)                             
           Mp2 = get_ang_mom_proj(tnp2)          
           CIp = get_CI(TargetStates2el%b(nstf),nep_con)  
           
           ! Looping over INITIAL Molecular State orbitals.  COORDINATE 2
           do ne_con = 1, Ncon       
              
              ind0 = TargetStates2el%b(nsti)%na(ne_con)  
              
              if ( nuse_norb /= ind0 ) cycle            
              
              ! Quantum numbers and functions for DELTA
              ind2 = TargetStates2el%b(nsti)%nb(ne_con) 
              tn2 => bst_nr%b(ind2)                                         
              f2 => fpointer(tn2)                               
              La2 = get_ang_mom(tn2)                            
              M2 = get_ang_mom_proj(tn2)            
              CI = get_CI(TargetStates2el%b(nsti),ne_con)  
              
              ! Selections Rules  
              if ( Lap2 /= La2 .OR. Mp2 /= M2  ) cycle
              
              temp_overlap2 = 0d0
              ! DO OVERLAP OF COORDINATE SPACE 2  
              minf2 = get_minf(tn2)
              maxf2 = get_maxf(tn2)
              minfp2 = get_minf(tnp2)
              maxfp2 = get_maxf(tnp2)      
              ir1 = max(minf2,minfp2)
              ir2 = min(maxf2,maxfp2)  
!              temp_overlap2 = SUM( fp2(ir1:ir2) *  f2(ir1:ir2) * weight(ir1:ir2) )
              temp_overlap2 = bst_nr%ortint(indp2,ind2)
              
              overlap2 = overlap2 +  CIp * CI * temp_overlap2
              
           end do  ! INITIAL STATE COORDINATE 2
        end do    ! FINAL STATE COORDINATE 2       

!!$     Here do one electron matrix elements and V_01 term:
!!$     Cp*C*<np2|n2>(<kp0|n0><np1|k1>(E-e_k'-e_k-Z^2/R) - <kp0|V_0-z.Zasy/r-U_0|n0><np1|k1>
!!$     - <kp0|n0><np1|V_1-z.Zasy/r-U_1|k1> - <kp0,np1|V_01|n0,k1>)
        if ( overlap2 == 0d0 ) cycle 
        
        ! FINAL STATE COORDINATE 1 
        ! Quantum numbers and functions for ALPHA
        indp1 = nuse_nporb            ! Final state number np. nep A.O.
        tnp1 => bst_nr%b(indp1)       !           
        fp1 => fpointer(tnp1)         ! One electron functions
        Lap1 = get_ang_mom(tnp1)      ! Gets Angular momentum A.O.
        Mp1 = get_ang_mom_proj(tnp1)  ! Get angular projection of A.O. 
        minfp1 = get_minf(tnp1)
        maxfp1 = get_maxf(tnp1)
        
        
        ! INITIAL STATE COORDINATE 0 
        ! Quantum numbers and functions for GAMMA
        ind0 = nuse_norb               
        tn0 => bst_nr%b(ind0)                                          
        f0 => fpointer(tn0)                               
        La0 = get_ang_mom(tn0)                            
        M0 = get_ang_mom_proj(tn0)  
        minf0 = get_minf(tn0)
        maxf0 = get_maxf(tn0)

        f0weight(minf0:maxf0) = f0(minf0:maxf0) * weight(minf0:maxf0)
        
!!$ Calculating -Cp*C*<np2|n2><kp0,np1|V_01|n0,k1>
!!!$ Inital COORDINATE 1 |kLM(r_1)>
        do ki=1,nqmi
           kqi = npk_nchi + ki - 1
           tn1 => chil%b(kqi)
           ki_en = get_energy(tn1)
           f1 => fpointer(tn1)
           minf1 = get_minf(tn1)
           maxf1 = get_maxf(tn1)
           
           or1 = max(minf1,minfp1)
           or2 = min(maxf1,maxfp1)
           !LHS: removed weights from fun1 so we can use form_accurate below
!           fun1(or1:or2) = ( fp1(or1:or2) * f1(or1:or2) * weight(or1:or2) )
           fun1(or1:or2) = ( fp1(or1:or2) * f1(or1:or2) )
 
           lambda_min = max(ABS(kLap0 - La0),ABS(Lap1 - kLa1))
           lambda_max = min(kLap0 + La0, Lap1 + kLa1)
           lambda_max = min(data_in%ltmax,lambda_max)
           ! mu should be the same for both integrals in coordinate 0 and 1
           mu = Mp1 - kM1
           do lambda = lambda_min, lambda_max

              ang0 = 0d0
              ang1 = 0d0
              ang0 = dble((-1)**(mu))*Yint(dble(kLap0),dble(kMp0),dble(lambda),dble(-mu),dble(La0),dble(M0))
              ang1 = Yint(dble(Lap1),dble(Mp1),dble(lambda),dble(mu),dble(kLa1),dble(kM1))             
              if ( ang0 * ang1 == 0d0 ) cycle
              
              ! Integrate over coordinate 1 
              call form_accurate(lambda,fun1,or1,or2,nr,epot_i,jr1,jr2)
              
              !LHS: added this for speedup 
              epot_i(jr1:jr2) = epot_i(jr1:jr2) * f0weight(jr1:jr2)

!!!$ COORDINATE 0 <kpLM(r_0)|
              do kf=1,nqmf
                 kqf = npk_nchf + kf - 1
                 tnp0 => chil%b(kqf)
                 fp0 => fpointer(tnp0)
                 minfp0 = get_minf(tnp0)
                 maxfp0 = get_maxf(tnp0)
                 
                 ir1 = max(minf0,minfp0)
                 ir2 = min(maxf0,maxfp0)
                 ir1 = max(jr1,ir1)
                 ir2 = min(jr2,ir2)

                 ! Integrate over coordinate 0 
                 
                 !LHS: Below line is original 
                 !temp = SUM(fp0(ir1:ir2) * epot_i(ir1:ir2) * f0(ir1:ir2) * weight(ir1:ir2)) 
                 
                 temp = SUM(fp0(ir1:ir2) * epot_i(ir1:ir2))

                 ! overlap2 has CI coefficients
                 temp = temp * ang0 * ang1 * overlap2
                 vmatt(kf,ki) = vmatt(kf,ki) - temp

              end do ! end kf loop
           end do ! lambda

           do kf=1,nqmf
              kqf = npk_nchf + kf - 1
              tnp0 => chil%b(kqf)
              kf_en = get_energy(tnp0)    
              
              temp = ortchil(ind0,kqf) * ortchil(indp1,kqi) * (Etot*(1d0-theta)-kf_en-ki_en)
              vmatt(kf,ki) = vmatt(kf,ki) + overlap2 * temp 
              if ( Rd /= 0d0 ) then
                 vmatt(kf,ki) = vmatt(kf,ki) - Z1*Z2/Rd * ortchil(ind0,kqf) * ortchil(indp1,kqi) * overlap2 
              end if
              vmatt(kf,ki) = vmatt(kf,ki) - flchil(ind0,kqf) * ortchil(indp1,kqi) * overlap2  
              vmatt(kf,ki) = vmatt(kf,ki) - ortchil(ind0,kqf) * flchil(indp1,kqi) * overlap2
              
           end do  ! end kf loop
        end do ! end ki loop
     end do ! INITIAL STATE COORDINATE 0   
  end do  ! FINAL STATE COORDINATE 1


!!$ Calculating -Cp*C*<np1|k1><kp0,np2|V_02|n0,n2>
!!$ Looping over FINAL Molecular State orbitals
  do nep_con =  1, Npcon  
     ! Quantum numbers and functions for ALPHA
     indp1 = TargetStates2el%b(nstf)%na(nep_con)
     
     ! FINAL STATE COORDINATE 2 
     ! Quantum numbers and functions for BETA
     indp2 = TargetStates2el%b(nstf)%nb(nep_con)
     tnp2 => bst_nr%b(indp2)
     fp2 => fpointer(tnp2)
     Lap2 = get_ang_mom(tnp2)
     Mp2 = get_ang_mom_proj(tnp2)
     CIp = get_CI(TargetStates2el%b(nstf),nep_con)
     minfp2 = get_minf(tnp2)
     maxfp2 = get_maxf(tnp2) 
     
     ! Check that coordinate 1 overlaps are non-zero
     Logic_V02=.FALSE.
     do ki=1,nqmi
        kqi = npk_nchi + ki - 1
        if ( ortchil(indp1,kqi) /= 0d0 ) then
           Logic_V02=.TRUE.
           exit
        end if
     end do ! ki
     if ( .NOT. Logic_V02 ) cycle
     
     ! Looping over INITIAL Molecular State orbitals
     do ne_con = 1, Ncon       
        
        ! INITIAL STATE COORDINATE 0 
        ! Quantum numbers and functions for GAMMA
        ind0 = TargetStates2el%b(nsti)%na(ne_con) 
        tn0 => bst_nr%b(ind0)
        f0 => fpointer(tn0)
        La0 = get_ang_mom(tn0)
        M0 = get_ang_mom_proj(tn0)
        minf0 = get_minf(tn0)
        maxf0 = get_maxf(tn0)
        
        ! INITIAL STATE COORDINATE 2 
        ! Quantum numbers and functions for DELTA
        ind2 = TargetStates2el%b(nsti)%nb(ne_con) 
        tn2 => bst_nr%b(ind2)                                         
        f2 => fpointer(tn2)                               
        La2 = get_ang_mom(tn2)                            
        M2 = get_ang_mom_proj(tn2)            
        CI = get_CI(TargetStates2el%b(nsti),ne_con)  
        minf2 = get_minf(tn2)
        maxf2 = get_maxf(tn2)
        
        or1 = max(minfp2,minf2)     
        or2 = min(maxfp2,maxf2)
        
        !LHS: removed weight so we can use form_accurate below
!        fun2(or1:or2) = fp2(or1:or2) * f2(or1:or2)  * weight(or1:or2)
        fun2(or1:or2) = fp2(or1:or2) * f2(or1:or2)
        
        lambda_min = max(ABS(kLap0 - La0),ABS(Lap2 - La2))
        lambda_max = min(kLap0 + La0, Lap2 + La2)
        lambda_max = min(data_in%ltmax,lambda_max)
        mu = Mp2 - M2
        do lambda = lambda_min, lambda_max
           ang0 = 0d0
           ang2 = 0d0
           ang0 = dble((-1)**(mu))*Yint(dble(kLap0),dble(kMp0),dble(lambda),dble(-mu),dble(La0),dble(M0))
           ang2 = Yint(dble(Lap2),dble(Mp2),dble(lambda),dble(mu),dble(La2),dble(M2))
           
           if ( ang0 * ang2 == 0d0 ) cycle
           
           ! Integrate over coordinate 2 
           call form_accurate(lambda,fun2,or1,or2,nr,epot_i,jr1,jr2)
           
!!!$ COORDINATE 0 <kpLM(r_0)|
           do kf=1,nqmf
              kqf = npk_nchf + kf - 1
              tnp0 => chil%b(kqf)
              fp0 => fpointer(tnp0)
              minfp0 = get_minf(tnp0)
              maxfp0 = get_maxf(tnp0)
              
              ir1 = max(minfp0,minf0)
              ir2 = min(maxfp0,maxf0)
              ir1 = max(jr1,ir1)
              ir2 = min(jr2,ir2)
              
              ! Integrate over coordinate 0 
              temp = SUM(fp0(ir1:ir2) * epot_i(ir1:ir2) * f0(ir1:ir2) * weight(ir1:ir2)) 
              temp = CIp * CI * temp * ang0 * ang2
              
! Overlap coordinate 1
              do ki=1,nqmi
                 kqi = npk_nchi + ki - 1
                 vmatt(kf,ki)= vmatt(kf,ki) - temp * ortchil(indp1,kqi) 
              end do ! ki

           end do ! kf
        end do ! lambda 
     end do ! INITIAL STATE 
  end do ! FINAL STATE


!!$ Calculating -Cp*C*<kp0|n0><np1,np2|V_12|k1,n2>
!!$ Looping over INITIAL Molecular State orbitals
  do ne_con = 1, Ncon        
     ! INITIAL STATE COORDINATE 0 
     ! Quantum numbers and functions for GAMMA
     ind0 = TargetStates2el%b(nsti)%na(ne_con) 
     
     ! INITIAL STATE COORDINATE 2 
     ! Quantum numbers and functions for DELTA
     ind2 = TargetStates2el%b(nsti)%nb(ne_con) 
     tn2 => bst_nr%b(ind2)                                         
     f2 => fpointer(tn2)                               
     La2 = get_ang_mom(tn2)                            
     M2 = get_ang_mom_proj(tn2)            
     CI = get_CI(TargetStates2el%b(nsti),ne_con)  
     minf2 = get_minf(tn2)
     maxf2 = get_maxf(tn2)
     
     ! Check that coordinate 0 overlaps are non-zero
     Logic_V12=.FALSE.
     do kf=1,nqmf
        kqf = npk_nchf + kf - 1
        if ( ortchil(ind0,kqf) /= 0d0 ) then
           Logic_V12=.TRUE.
           exit
        end if
     end do ! kf
     if ( .NOT. Logic_V12 ) cycle
     
!!$ Looping over FINAL Molecular State orbitals
     do nep_con =  1, Npcon  
        ! FINAL STATE COORDINATE 1 
        ! Quantum numbers and functions for ALPHA
        indp1 = TargetStates2el%b(nstf)%na(nep_con)
        tnp1 => bst_nr%b(indp1)       !           
        fp1 => fpointer(tnp1)         ! One electron functions
        Lap1 = get_ang_mom(tnp1)      ! Gets Angular momentum A.O.
        Mp1 = get_ang_mom_proj(tnp1)  ! Get angular projection of A.O. 
        minfp1 = get_minf(tnp1)
        maxfp1 = get_maxf(tnp1)
        
        ! FINAL STATE COORDINATE 2 
        ! Quantum numbers and functions for BETA
        indp2 = TargetStates2el%b(nstf)%nb(nep_con)
        tnp2 => bst_nr%b(indp2)
        fp2 => fpointer(tnp2)
        Lap2 = get_ang_mom(tnp2)
        Mp2 = get_ang_mom_proj(tnp2)
        CIp = get_CI(TargetStates2el%b(nstf),nep_con)
        minfp2 = get_minf(tnp2)
        maxfp2 = get_maxf(tnp2) 
        
        or1 = max(minfp2,minf2)     
        or2 = min(maxfp2,maxf2)

        !LHS: removed weight so we can use form_accurate below
!        fun2(or1:or2) = fp2(or1:or2) * f2(or1:or2) * weight(or1:or2)
        fun2(or1:or2) = fp2(or1:or2) * f2(or1:or2)
        
        lambda_min = max(ABS(Lap1 - kLa1), ABS(Lap2 - La2))
        lambda_max = min(Lap1 + kLa1, Lap2 + La2)
        lambda_max = min(data_in%ltmax,lambda_max)
        mu = Mp2 - M2
        do lambda = lambda_min, lambda_max
           ang1 = 0d0
           ang2 = 0d0
           ang1 = dble((-1)**(mu))*Yint(dble(Lap1),dble(Mp1),dble(lambda),dble(-mu),dble(kLa1),dble(kM1))
           ang2 = Yint(dble(Lap2),dble(Mp2),dble(lambda),dble(mu),dble(La2),dble(M2))
           if ( ang1 * ang2 == 0d0 ) cycle

           ! Integrate over coordinate 2 
           call form_accurate(lambda,fun2,or1,or2,nr,epot_i,jr1,jr2)

           ! Initial COORDINATE 1 |kLM(r_1)>           
           do ki = 1, nqmi
              kqi = npk_nchi + ki - 1
              tn1 => chil%b(kqi)
              f1 => fpointer(tn1)
              minf1 = get_minf(tn1)
              maxf1 = get_maxf(tn1)
              
              ir1 = max(minfp1,minf1)
              ir2 = min(maxfp1,maxf1)
              ir1 = max(jr1,ir1)
              ir2 = min(jr2,ir2)            
              ! Integrate over coordinate 1 
              temp = SUM(fp1(ir1:ir2) * epot_i(ir1:ir2) * f1(ir1:ir2) * weight(ir1:ir2)) 
              temp = CIp * CI * temp * ang1 * ang2
              
              ! Overlap coordinate 0 
              do kf=1,nqmf
                 kqf = npk_nchf + kf - 1
                 vmatt(kf,ki)= vmatt(kf,ki) - temp * ortchil(ind0,kqf) 
              end do ! kf
           end do ! ki
        end do ! lambda
     end do  ! FINAL STATE 
  end do ! INITIAL STATE    


!!$ Calculating -Cp*C*<kp0|n0><np1|k1><np2|H_2|n2>
!!$ Looping over INITIAL Molecular State orbitals
  do ne_con = 1, Ncon
     ind0 = TargetStates2el%b(nsti)%na(ne_con)
     ind2 = TargetStates2el%b(nsti)%nb(ne_con)
     CI = get_CI(TargetStates2el%b(nsti),ne_con)
!!$ Looping over FINAL Molecular State orbitals
     do nep_con =  1, Npcon
        indp1 = TargetStates2el%b(nstf)%na(nep_con)
        indp2 = TargetStates2el%b(nstf)%nb(nep_con)
        CIp = get_CI(TargetStates2el%b(nstf),nep_con)
        if ( bst_nr%ham1el(indp2,ind2) == 0d0 ) cycle

        do ki = 1, nqmi
           kqi = npk_nchi + ki - 1
           if (ortchil(indp1,kqi) == 0d0) cycle

        do kf=1,nqmf
           kqf = npk_nchf + kf - 1
           if (ortchil(ind0,kqf) == 0d0 ) cycle
                 temp = CIp * CI * ortchil(ind0,kqf) * ortchil(indp1,kqi) 
                 vmatt(kf,ki)= vmatt(kf,ki) - temp * bst_nr%ham1el(indp2,ind2) 
        end do ! kf
        end do ! ki
   end do  ! FINAL STATE 
  end do ! INITIAL STATE 


!!$ Spin and other coefficients
!!$ Need to be done before projection operator overlaps
!!$ + 2. Spin_coeff.(E(1-theta)-H)P_01, vmatt=(E(1-theta)-H)P_01 
  temp = sqrt(dble(2*Spinp+1)) * sqrt(dble(2*Spin+1))
  temp = COF6J(0.5d0,0.5d0,dble(Spinp),0.5d0,dble(rspin_tot),dble(Spin)) * temp
  temp = 2d0 * dble((-1)**(Spinp+Spin+1)) * temp

  vmatt(:,:) = temp * vmatt(:,:)


!!$ Projection operator part of non-uniqueness
!!$ 2(E-theta+theta)P_01
!!$ -2 E * theta * <Phi_f|Phi_i> * <kp| I_0 |k>
  if ( nsti == nstf .AND. theta /= 0d0) then  ! elastic channel only???
     ! Looping over INITIAL Molecular State orbitals.  COORDINATE 2
     do ne_con = 1, Ncon
        ! Quantum numbers and functions for DELTA
        ind2 = TargetStates2el%b(nsti)%nb(ne_con)
        tn2 => bst_nr%b(ind2)
        f2 => fpointer(tn2)
        La2 = get_ang_mom(tn2)
        M2 = get_ang_mom_proj(tn2)
        CI = get_CI(TargetStates2el%b(nsti),ne_con)

!        if ( is_core_orb(ind2) == 0 ) cycle

        ! INITIAL STATE COORDINATE 0 
        ! Quantum numbers and functions for GAMMA
        ind1 = TargetStates2el%b(nsti)%na(ne_con)
        tn1 => bst_nr%b(ind1)
        f1 => fpointer(tn1)
        La1 = get_ang_mom(tn1)
        M1 = get_ang_mom_proj(tn1)
        minf1 = get_minf(tn1)
        maxf1 = get_maxf(tn1)
      ! Selections Rules  
        if ( Lap1 /= La1 .OR. Mp1 /= M1  ) cycle

        ! Looping over FINAL Molecular State orbitals. COORDINATE 2
        do nep_con =  1, Npcon
           ! FINAL STATE COORDINATE 1 
           ! Quantum numbers and functions for ALPHA
           indp1 = TargetStates2el%b(nstf)%na(nep_con)    ! Final state number np. nep A.O.
           tnp1 => bst_nr%b(indp1)       !           
           fp1 => fpointer(tnp1)         ! One electron functions
           Lap1 = get_ang_mom(tnp1)      ! Gets Angular momentum A.O.
           Mp1 = get_ang_mom_proj(tnp1)  ! Get angular projection of A.O. 
           minfp1 = get_minf(tnp1)
           maxfp1 = get_maxf(tnp1)

!           if ( is_core_orb(ind2) /= 2 .AND. is_core_orb(nuse_nporb) == 0 ) cycle

          ! Quantum numbers and functions for BETA
          indp2 = TargetStates2el%b(nstf)%nb(nep_con)
          tnp2 => bst_nr%b(indp2)
          fp2 => fpointer(tnp2)
          Lap2 = get_ang_mom(tnp2)
          Mp2 = get_ang_mom_proj(tnp2)
          CIp = get_CI(TargetStates2el%b(nstf),nep_con)
  
          ! Selections Rules  
          if ( Lap2 /= La2 .OR. Mp2 /= M2  ) cycle

!          temp_overlap2 = 0d0
!          ! DO OVERLAP OF COORDINATE SPACE 2  
!          minf2 = get_minf(tn2)
!          maxf2 = get_maxf(tn2)
!          minfp2 = get_minf(tnp2)
!          maxfp2 = get_maxf(tnp2)
!          ir1 = max(minf2,minfp2)
!          ir2 = min(maxf2,maxfp2)
!         temp_overlap2 = SUM( fp2(ir1:ir2) *  f2(ir1:ir2) * weight(ir1:ir2))
          overlap12 = CIp * CI * bst_nr%ortint(indp1,ind1)*bst_nr%ortint(indp2,ind2) 

          if (overlap12 == 0d0 ) cycle

          do ind0 = 1, TargState_Orthonorm%Nmax  ! Loop over all one-electron target states 
!!$ Note the difference between is_core_orb and is_core_MO
!             if ( is_core_orb(ind2) /= 2 .AND. is_core_MO(ind0) == 0) cycle
             do ki = 1, nqmi
                kqi = npk_nchi + ki - 1
                do kf = 1, nqmf
                   kqf = npk_nchf + kf - 1
                   temp = ovlp_OrtnrmMO_k(ind0,kqi) * ovlp_OrtnrmMO_k(ind0,kqf)     
                   vmatt(kf,ki) = vmatt(kf,ki) - 2d0 * temp * Etot * theta * overlap12
                end do ! kf
             end do ! ki
          end do ! ind0 I_0 



       end do  ! INITIAL STATE CONFIGURATIONS

    end do    ! FINAL STATE CONFIGURATIONS

  end if ! theta non-uniqueness
!  do ind0 = 1, Nmax_1e_orb  ! Loop over all one-electron orbitals
!  do ki = 1, nqmi
!     kqi = npk_nchi + ki - 1
!     if (ortchil(ind0,kqi) == 0d0) cycle
!          do kf=1,nqmf
!           kqf = npk_nchf + kf - 1
!           if (ortchil(ind0,kqf) == 0d0 ) cycle
!
!          temp = Etot * theta * ortchil(ind0,kqf) * ortchil(ind0,kqi)
!          vmatt(kf,ki) = vmatt(kf,ki) - 2d0 * temp
!         end do ! kf
!  end do ! ki
!  end do ! one-electron orbital loop
!  end if ! elastic channel non-uniqueness


  vmatt(:,:) = vmatt(:,:) * (2d0/pi) !!! / (rkf*rki) dealt with in scat.f 
  
end subroutine vexch_2e_sp

end module vmat_subroutines
