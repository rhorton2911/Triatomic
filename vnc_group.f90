!!$ vnc(r) = VLambda(r)  + Vcore_electrons(r) + pol(r)
!!$  r -> 0,    vnc(r) -> -Z/r
!!$  r -> inf,  vcentre(r) -> -(Z - N_core_el)/r = Zas/r 
!!$
!!$ vdcore(r) = -(N_core_el)/r + Vcore_electrons(r) + pol(r)

  subroutine construct_vnc_group(nr, gridr, indata)
    use input_data
    use vnc_module
    use MPI_module
		use basismodule
    implicit none

    integer, intent(in):: nr
    real*8, dimension(nr), intent(in):: gridr
!  
    real*8:: polpot
    integer:: nsize, n, lam, lcore, lbar, npr_core, kappa, l
    real*8, dimension(nr):: v, temp
    integer:: i, i1, i2,  ntmp, lamtop_vc_set
    real*8:: r, j
    real*8, dimension(0:20):: corep, r0, exchlocal
    real*8:: Z1, Z2, Zasym
    logical :: core

		type(smallinput):: indata  !Contains nuclear coords for H3+
		integer:: vncsize
		integer:: vlim1, vlim2

    if(data_in%pseudo_pot) then
      data_in%Z1 = data_in%Z1 - data_in%N_core_el
      core_energy = data_in%core_energy
    endif

    Z1 = data_in%Z1
    Z2 = data_in%Z2
    Zasym = data_in%Zasym
    if (.not. data_in%hlike) then
       Zasym = data_in%Zasym + 1
    endif

    temp(:) = 0.0d0
    exchlocal(:) = 0.0d0
    corep(:) = 0.0d0
    r0(:) = 2.0d0


    lamtop_vc = data_in%ltmax  ! change in future if required...
		if (data_in%good_m) then
       allocate(vnc(nr,0:lamtop_vc),minvnc(0:lamtop_vc),maxvnc(0:lamtop_vc))
			 vlim1 = 0
			 vlim2 = lamtop_vc
	  else 
			 vncsize = (data_in%ltmax+1)**2   !now loops over (lam, mu)
			 allocate(vnc(nr,vncsize),minvnc(vncsize),maxvnc(vncsize)) 
			 vlim1 = 1
			 vlim2 = vncsize
	  end if
    vnc(:,:) = 0d0

    if(myid == 0 .and. data_in%good_parity) then
       write(*,*) "Setting up nuclear potential: only even terms (good parity)"
    elseif((myid == 0 .and. .not.data_in%good_parity) .and. data_in%good_m) then
       write(*,*) "Setting up nuclear potential: even and odd terms (no good parity)"
    elseif(myid == 0) then 
       write(*,*) "Setting up nuclear potential: even and odd terms, with a loop over projection mu (no good parity or m)"
    endif

    call  VLambdaMuR_group(nr,gridr,vlim1,vlim2,lamtop_vc, data_in%Rd, Z1, Z2, vnc,minvnc,maxvnc,lamtop_vc_set, indata)
   
    if ((data_in%N_core_el > 0 .and..not.data_in%pseudo_pot) .and. (data_in%good_m)) then
      call make_vdcore
      if(ltop_vc > lamtop_vc) then
        if(myid==0)then
          print*, '*** ERROR: max lambda for vdcore > lamtop for nuclear potential'
          print*, ltop_vc, '>', lamtop_vc
        endif
        error stop
      endif

      do lam=0, ltop_vc
        vnc(:,lam) = vnc(:,lam) + vdcore(:,lam)
      enddo

    else
      minvdc = minval(minvnc)
      maxvdc = maxval(maxvnc)
    endif

    if( lamtop_vc_set .ne. lamtop_vc ) then

       lamtop_vc = lamtop_vc_set   ! interesting :)

    endif

    !   if (myid == 0)  print*, 'module vnc: test vc: ', vnc(1,0),vnc(1,0)*gridr(1), Za

!!$ set up central joint nuclear potential to be used in analytical Born subtraction (p.w. vdirect and anaytical FBA  vdirect)     
    allocate(vcentre_jn(nr))
    vcentre_jn(:) = 0d0
    vcentre_jn(1:nr) = -(Z1+Z2) / gridr(1:nr)

!!$    if (myid == 0) then
!!$       open(124,file='tmpfile')
!!$       do n=1,nr,4
!!$          write(124,'(F15.6,2F20.8)') gridr(n), vcentre(n,0),  vdcore(n,0)
!!$       end do
!!$       close(124)
!!$    end if
    
  if(myid == 0) write(*,*)

!  open(unit=100,file='vnc.out',action='write',status='replace')
!  do i=minval(minvnc), maxval(maxvnc)
!    write(100,*) gridr(i), vnc(i,0:2), - data_in%corep/(2.0d0*gridr(i)**4)*(1.0d0-exp(-(gridr(i)/data_in%r0)**2))**2 - data_in%corep/(2.0d0*data_in%Rd**4)
!  enddo
!  close(100)

  end subroutine construct_vnc_group
 !------------------------------------------------------------------------------

    subroutine VLambdaMuR_group(maxr,gridr,vlim1,vlim2,lambda_max, Rd, Z1_in, Z2_in, vlr,minvnc,maxvnc, lamtop_vc_set, indata)
      use MPI_module
      use input_data
			use basismodule
			use grid_radial
      implicit none

			type(smallinput):: indata  !Contains nuclear coordinates
      
      integer, intent(in):: maxr
      real*8, dimension(maxr),intent(in):: gridr
      integer, intent(in):: lambda_max 
      real*8, intent(in):: Rd
      real*8, intent(in):: Z1_in, Z2_in  
			integer:: vlim1, vlim2
      real*8, dimension(maxr,vlim1:vlim2), intent(out):: vlr 
      integer, dimension(vlim1:vlim2), intent(out):: minvnc,maxvnc 
      integer, intent(out):: lamtop_vc_set
      logical :: good_parity
      real*8 :: Z1, Z2

      integer::  i, lambda, iRd, i1, i2, lambda_max_local, iRd1, iRd2
      real*8:: const1, r, Rdh, Zcoef, Rd1, Rd2, origin
      real*8, dimension(maxr):: temp

			integer:: num_lambda, harm, ii, nr
      complex*16, dimension(:,:), allocatable:: VPot, VPotTemp

      if(data_in%pseudo_pot .and. .false.) then
        Z1 = Z1_in - data_in%N_core_el
      else
        Z1 = Z1_in
      endif
      Z2 = Z2_in


      lambda_max_local = lambda_max

      origin = data_in%origin

      good_parity = data_in%good_parity

      vlr(:,:) = 0d0
      maxvnc(:) = 1

      Rdh = Rd / 2d0  !  distance from COM to a nuclear   (only when data_in%origin = 0.5

      !Liam added for origin not at geometric centre: Rd1, Rd2 distances from origin to Z1 and Z2
      Rd1 = Rd * origin
      Rd2 = Rd * (1.0d0 - origin)
      
      if ((Rd .eq. 0d0) .and. (data_in%good_m)) then ! joint nuclear, spherically symmetric case, only lam=0 term
         
         iRd = 1
         lambda_max_local = 0

			else if (data_in%good_m) then
         do i = 1, maxr
            r = gridr(i)
            if(r .gt. Rdh) exit
         enddo
         iRd = i
        
         !Liam added: get locations in radial grid where r = Rd1, Rd2
         do i = 1, maxr
            r = gridr(i)
            if(r .gt. Rd1) exit
         enddo
         iRd1 = i
         do i = 1, maxr
            r = gridr(i)
            if(r .gt. Rd2) exit
         enddo
         iRd2 = i

      endif

			if (data_in%good_m) then    !Diatomic molecule/atom, loop over lambda
         do lambda = 0, lambda_max_local

            if(good_parity .and. (-1)**lambda /= 1) cycle

            Zcoef = Z1 + Z2*(-1)**lambda

            temp(:) = 0d0

            if((origin == 0.5d0 .or. Rd == 0.0d0) .and. (data_in%good_m)) then !simpler expansion when origin is at geometric centre
            
               do i = 1, iRd-1
                  r = gridr(i)
                  temp(i) = -Zcoef * (r**lambda/Rdh**(lambda + 1))
               end do
               do i = iRd, maxr
                  r = gridr(i)
                  temp(i) = -Zcoef * (Rdh**lambda/r**(lambda + 1))           
               end do

			   	 else if (data_in%good_m) then !Liam added for arbitrary origin
             
               !Contribution from Z1:
               do i = 1, iRd1-1
                  r = gridr(i)
                  temp(i) = -Z1 * (r**lambda/Rd1**(lambda + 1))
               end do
               do i = iRd1, maxr
                  r = gridr(i)
                  temp(i) = -Z1 * (Rd1**lambda/r**(lambda + 1))           
               end do

               !Contribution from Z2:
               do i = 1, iRd2-1
                  r = gridr(i)
                  temp(i) = temp(i) - (-1)**lambda * Z2 * (r**lambda/Rd2**(lambda + 1))
               end do
               do i = iRd2, maxr
                  r = gridr(i)
                  temp(i) = temp(i) - (-1)**lambda * Z2 * (Rd2**lambda/r**(lambda + 1))           
               end do

            endif

            if(data_in%pseudo_pot .and. lambda <= 2) then
              temp = temp + data_in%pseudo_pot_B(lambda) * exp(-data_in%pseudo_pot_Beta(lambda)*gridr**2)
              if(data_in%corep > 0.0d0) then
                if(lambda == 0) temp = temp - data_in%corep/(2.0d0*gridr**4)*(1.0d0-exp(-(gridr/data_in%r0)**2))**2 &
                                       - data_in%corep/(2.0d0*data_in%Rd**4)
                if(lambda == 1) temp = temp + data_in%corep/(gridr**2*data_in%Rd**2)*(1.0d0-exp(-(gridr/data_in%r0)**2))
             endif
            endif

            call minmaxi(temp,maxr,i1,i2)
            vlr(i1:i2,lambda) = temp(i1:i2)
            minvnc(lambda) = i1
            maxvnc(lambda) = i2
         end do !lambda
	    else if (.not. data_in%good_m) then !H3+ mode
				 !Subroutine loops over lambda, q
         !Number of (lam,mu) pairs: use formula SUM_l=0^l=L (2l+1) = (L+1)^2 nr = grid%nr
				 nr = size(gridr)
         num_lambda = (lambda_max_local+1)**2   
				 !indata%lambdamax = lambda_max_local
         allocate(VPot(nr,num_lambda))
         VPot(:,:) = 0.0_dpf

         allocate(VPotTemp(nr,num_lambda))
         do ii = 1, 3
            call getVPotNuc(grid, VPotTemp, indata%R(ii), indata%theta(ii), &
                     indata%phi(ii), indata%charge(ii), indata)
            VPot(:,:) = VPot(:,:) + VPotTemp(:,:)
         end do
         deallocate(VPotTemp)
       
				 do i = 1, num_lambda
            call minmaxi(real(VPot(:,i)),maxr,i1,i2)
            vlr(i1:i2,i) = real(VPot(i1:i2,i))
            minvnc(i) = i1
            maxvnc(i) = i2
				 end do

      	 harm = indata%harmop
      	 data_in%harmop = harm
      	 deallocate(VPot)
			end if

      lamtop_vc_set = lambda_max_local

    end subroutine VLambdaMuR_group
   
!-----------------------------------------------------------------------------------------
!!$  the potential has rank lam (up to max value lamtop_vc) and its projection=0
!!$  so that the magnetic sublevel projections of pi and pj must be the same
    function VLambdaR_ME(pi,pj,ma)
      use input_data
      use grid_radial
      use sturmian_class
      use vnc_module
      use MPI_module
      use state_class
      implicit none

      real*8:: VLambdaR_ME
      type(sturmian_nr), intent(in):: pi, pj
      integer, intent(in):: ma

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

      VLambdaR_ME = 0d0

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
      do lam=lam_min,lam_max

         if(data_in%good_parity .and. (-1)**lam /= 1) cycle

         if(lam < abs(li-lj) .or. lam > li+lj) cycle
         tmp1 = Yint(dble(li),dble(ma),dble(lam),dble(0),dble(lj),dble(ma))
         !  if (myid == 0) print*, i,j, lam, tmp1
         i1 = max(minvnc(lam),minf)
         i2 = min(maxvnc(lam),maxf)
         tmp = tmp + tmp1 * SUM(weight(i1:i2)*fi(i1:i2)*fj(i1:i2)*vnc(i1:i2,lam))
        if(data_in%N_core_el > 0 .and. .not.data_in%pseudo_pot) then !add core exchange ME
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

      VLambdaR_ME = tmpres
      
    end function VLambdaR_ME

subroutine make_vdcore
  use vnc_module
  use input_data
  use grid_radial
  use state_class
  implicit none

  integer :: j, jj, jjp, lam, lj, ljp, mj, mjp, i1, i2, i1out, i2out, na, nap
  real*8 :: ang, CI, angp, CIp
  type(sturmian_nr), pointer :: sturmj, sturmjp
  type(state), pointer :: orbj
  real*8, external :: Yint
  real*8, dimension(grid%nr) :: invec, outvec
  real*8, parameter :: pi = acos(-1.0d0)

  ltop_vc = 2*get_max_l(wfcore_nr) !selection rule from <Y1|Y2|Y3>, wfcore_nr holds core orbitals (vnc_module.f90)
    
  allocate(vdcore(grid%nr,0:ltop_vc))
  vdcore = 0.0d0

  do j=1, CoreOrbitals%Nmax

    orbj => CoreOrbitals%b(j)

    do jj=1, orbj%nam

      na = orbj%na(jj)
      CI = orbj%CI(jj)
      sturmj => wfcore_nr%b(na)
      lj = sturmj%l 
      mj = sturmj%m
    
      do jjp=1, orbj%nam

        nap = orbj%na(jjp)
        CIp = orbj%CI(jjp)
        sturmjp => wfcore_nr%b(nap)
        ljp = sturmjp%l 
        mjp = sturmjp%m
       
        if(mjp /= mj) error stop '*** ERROR in make_vdcore: mjp /= mj'

        i1 = max(sturmjp%minf,sturmj%minf)
        i2 = min(sturmjp%maxf,sturmj%maxf)
      
        invec(i1:i2) = sturmjp%f(i1:i2)*sturmj%f(i1:i2)

        do lam=0, ltop_vc
          ang = Yint(dble(ljp),dble(mjp),dble(lam),0.0d0,dble(lj),dble(mj)) !* sqrt(2.0d0*pi/dble(2*lam+1))
          if(ang == 0.0d0) cycle
          call form_accurate(lam,invec,i1,i2,grid%nr,outvec,i1out,i2out)
          vdcore(i1out:i2out,lam) = vdcore(i1out:i2out,lam) + outvec(i1out:i2out)*ang*CI*CIp
        enddo !lam

      enddo !jjp
    enddo !jj
  enddo !j

  vdcore = 2.0d0 * vdcore

  if(data_in%corep > 0.0d0) then
    vdcore(:,0) = vdcore(:,0) - data_in%corep/(2.0d0*grid%gridr**4)*(1.0d0-exp(-(grid%gridr/data_in%r0)**2))**2 &
    - data_in%corep/(2.0d0*data_in%Rd**4)

    vdcore(:,1) = vdcore(:,1) + data_in%corep/(grid%gridr**2*data_in%Rd**2)*(1.0d0-exp(-(grid%gridr/data_in%r0)**2))
  endif

end subroutine make_vdcore
