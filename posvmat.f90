 
subroutine psformation(Mtot,Etotp,nchf,nchi,nstf,nsti,nqmf,nqmi,Qlp,Q0p,Alp,aq,ionsh,nqgen,nqim,xp,wp,imax,h2tops,vmatt)

  use sturmian_class
  use target_states
  use channels
  use input_data
  use kgrid 
  use Positronium
  use one_electron_func
  use  state_class 
  use vmat_module

      
  
  implicit none
  integer, parameter:: maxl=2*20+1
  integer, intent(in):: Mtot
  integer, intent(in):: nchf, nchi        ! Channel number
  integer, intent(in):: nstf, nsti        ! State for particular channel
  integer, intent(in):: nqmf, nqmi        ! Number of k-grid points for particular channel
  real*8, dimension(nqmf,nqmi), intent(inout):: vmatt
  type(basis_sturmian_nr):: chil
  integer:: ki, kf, kqi, kqf
  real*8:: pi
  integer:: Spin, Spinp
  integer:: i,nst_ps
  
  real*8:: gp(Ps%iglp), Bnlp(0:4,Ps%iglp),Dnlp(0:4,Ps%iglp)
  real*8:: Qlp(nqmax,2*nqmax*Ps%igpm),Q0p(nqmax,2*nqmax*Ps%igpm),&
     &       Alp(nqmax,2*nqmax*Ps%igpm)
  real*8 aq(0:nqmax),qq(nqmax),xp(2*nqmax*Ps%igpm+nqmax)&
     &                ,wp(2*nqmax*Ps%igpm)
  integer:: JK,ionsh,nqim,imax,nqgen,nchi1,Li,la,iglp,mnqmf,mnqmi,iq,iz,igz,nn,ip,nnmin,nnmax,npos,intp,nnmax1,llf
  real*8:: qb,qa,qb2,qa2,qbqa,pa2,pp,pp2,pb2,qbpp, qbE,Etot,Etotp,sumf0z,xx(1:5),yy(1:5),xz,ddp,pdp,nd1,temp,efactor
  real*8:: xgz(Ps%igzm),wgz(Ps%igzm),pl(0:20,Ps%igzm),f0z(Ps%igzm), &
       & f1z(5,Ps%igzm,imax),&
       & fpqb1(5,0:19,2*nqmax*Ps%igpm+nqmax) 
  real*8:: resPs0(Ps%igzm),resPs1(Ps%igzm),resBnlp(Ps%igzm),resDnlp(Ps%igzm),resPs0f,resPs1f,gpmin,gpmax 
  real*8:: res0(0:19),sumlb1,flb1,flb2,sumlb2,sumla1,fla1,clb,qal2,fla2,cla,suml1,fl1,suml2,fl2,reslam,flam,allcoefs,cgc1,cgc2,cgc3,cgc4
  real*8:: wigner,resiz,cof3j,resi,rlf,rli,restest,c0,result,resl1,fla,fma,rMtot,rMi,rMf
  integer:: lb1,lb2,  lab2,l1,l2,lam,lf,ipl,n_ps,Mi,Mf,ma
  real*8, dimension(0:maxl):: factrl,sqrfct,hat
  real*8:: fpqb(0:19,2*nqmax*Ps%igpm+nqmax),gnl0(Ps%iglp),gnl1(Ps%iglp),rJK,cof6j,cgc,cgc0,resJK
  real*8:: alfa2, alfa,res1,res1b,res1c,res2,res3,res4,res12,fp,resAlp,resal
  integer:: ising,igp,mllf,mla,ka,kam, la1,la2 
  integer:: ic, in,nc, Lia1,Lia2,n1,n2,M1,M2
  logical:: h2tops ! to mark the channel from H2 to Ps or opposite
  real*8:: res0ps,res1ps
  real*8:: cof9j,resl2,resla1,reslb1,reska,reskam, Clb1,Cla1,Cl1,Cl2,Clam,Cka,Ckam,fka,fkam,fmla,fmllf,fllf
  real*8:: CI,const,sum_Bnlp,sum_Dnlp
  type(sturmian_nr), pointer:: tn1,tn2,tn2_ion  ! One-electron orbitals
  real*8:: dwp(2*nqmax*Ps%igpm)
  integer:: iqstep,iqstep1,iqs, nu
  integer:: icalam(0:19),itest(40,0:19)
  real*8:: res1_nu_lam(40,0:19)
  integer:: illi
!!!
  real*8:: rlam2, a(ncmax1), f0(0:ncmax1-1),f1(0:ncmax1-1)
  integer:: Nlf,k,nfa
  real*8:: pc0(0:ncmax1-1,ncmax1,0:3),pc1(0:ncmax1-1,ncmax1,0:3)
  real*8::resPs0old,resPs1old

  pi = acos(-1d0)
  Li = Lp_ch(nchi); Mi=Mp_ch(nchi); Lf = Lp_ch(nchf);Mf=Mp_ch(nchf);

  if(Li.gt.8.or.Lf.gt.8) RETURN

  rMtot=float(Mtot); rLi=float(Li);rMi=float(Mi); rLf=float(Lf); rMf=float(Mf);
  
  iglp = Ps%iglp; gp(1:iglp) = Ps%gp(1:iglp);

        do illi=1,Ps%icount_l(nsti)
            la = Ps%l_nsti(nsti,illi)    
!       write(15,*)
!       write(15,'(a1,4i5)'),'#',nsti,Ps%icount_l(nsti) ,la,illi

           do ip=1,iglp
        Bnlp(la,ip) = Ps%Bnlp(nsti,la,ip)
        Dnlp(la,ip) = Ps%Dnlp(nsti,la,ip) 
 
!       write(15,'(4e20.10)') gp(ip),Bnlp(la,ip),Dnlp(la,ip),Ps%Amg(ip)

           enddo
        enddo

  ddp=1.e+08; pdp=0.0065d0;  nd1=1; intp=2*nd1+1; nnmin=INT(Log(1.E-04*ddp)/pdp); nnmax1=INT(Log(2.0e+02*ddp)/pdp)
  gpmin=gp(nnmin)
  gpmax=gp(nnmax1)
  
  nst_ps = st2_ch(nchf); !gnl0(1:iglp) = Ps%gnl0(1:iglp,nst_ps); gnl1(1:iglp) = Ps%gnl1(1:iglp,nst_ps); 
  llf = Ps%l(nst_ps);fllf=float(llf); n_ps = Ps%n(nst_ps); mllf = Mtot - Mf; fmllf = float(mllf)
  
!    IF((Mf-Mi).ne.(mla-mllf)) RETURN
!  IF((-1)**(Li+la) /= (-1)**(Lf+llf)) RETURN ! assumes FC structure of H2
 ! nst_ps is the order number of Ps-states
!       write(21,*)
!       write(21,'(a1,2i5)') '#', n_ps,llf
!      do ip = 1, iglp
!         pb2 = gp(ip)
!         call  pswfp(n_ps,llf,pb2,res0ps,res1ps) 
!         write(21,'(3e20.10,i5)') pb2,res0ps,res1ps,ip
!      enddo

  Etot = Etotp - Target_Basis%en_ion
  igz = Ps%igzm; xgz(1:igz) = Ps%xgz(1:igz); wgz(1:igz) = Ps%wgz(1:igz);pl(0:20,1:igz)=Ps%pl(0:20,1:igz);

  factrl(0:maxl)=Ps%factrl(0:maxl); sqrfct(0:maxl)=Ps%sqrfct(0:maxl); hat(0:maxl)=Ps%hat(0:maxl);

!   c0 = 1.d0/pi * (-1)**(Lf)*hat(Li)*sqrfct(2*la+1)*hat(llf)*sqrfct(2*llf+1)*(-1)**((la-llf+Li-Lf)/2)

  alfa=1.d0    ! 
  alfa2=alfa*alfa
!!!!!!!!!!!!!!!!!!!!!!!!!!
      Nlf=npsstates(2,llf)
      rlam2=rlambda(2,llf)
      nfa = n_ps - llf
      npos = nfa
!!       print*,'Nlf:',Nlf,nfa,llf,rlam2
      if(Nlf.ne.0) then
         do k=1,Nlf
            a(k)=cknd(k,nfa,llf)
!!           print*,'a(k):',a(k)
          end do

         select case(llf)
         case(0)
            call coefsS(llf,Nlf,a,f0,f1)
         case(1)
            call coefsP(llf,Nlf,a,f0,f1)
         case(2)
            call coefsD(llf,Nlf,a,f0,f1)
         case default
            stop 'coeffs are missing for lb>2'
         end select ! llf

         do k=0,Nlf-1
            pc0(k,nfa,llf)= f0(k)
            pc1(k,nfa,llf)= f1(k)
!!         print*, 'pc0&pc1:', k,nfa,llf,pc0(k,nfa,llf),pc1(k,nfa,llf)
         end do
       endif
!!!!!!!!!!!!!!!!!!!!!!!!
!       STOP'test1'

  igp = Ps%igpm
! iqstep below is to speed up the q-integral when there are no singularities
  iqstep = 1
  iqs = 0 
  dwp(:) = wp(:) !0.d0
  IF(MOD(igp,3).eq.0) THEN
     iqstep = 3
     iqs = 1
     dwp(:) = wp(:)
   if(iqs==1) then
     do i=1+iqs,2*igp*nqim,iqstep
        dwp(i) = wp(i-1) + wp(i) + wp(i+1)
     enddo
   endif
  ELSE
     STOP 'change igp to 6 or 9 or 12, it must be multiple of 3'
  ENDIF

!      PRINT*,'l of target = ',la,nsti,nstf

!      RETURN  
  
!**************************         
  

        DO kf = 1, nqmf      ! 1 - for born  ! kfstart,kfstop
!c           print*, kf, mnqmf,gridk(kf,nchf)
           qb=dble(gridk(kf,nchf))
           if(qb.lt.0.) cycle !.or.qb.gt.20.) cycle
           qb2=qb*qb
           !c          print '(i3,$)', kf
           iq=0
           sumf0z=0.
           do ki=1,nqmi ! 1 - for born !kistart, kistop
              qa=dble(gridk(ki,nchi))
              !c         print*, 'gridk', nchf,npos,qb,qa
              qa2=qa*qa
              if(qa.lt.0.) cycle !.or.qa.gt.20) cycle
              !!      IF((abs(qa-qb).lt.0.1.or.abs(qa-qb/2).lt.0.1).and.qa.gt.1) CYCLE
              iq = iq+1
              
              qa2=qa*qa
              qbqa=qb*qa
              qbE=qa2 - Etot 

              result = 0.d0
              
!  do ic = 1, TargetStates2el%b(nsti)%nusemax 
!     
!     nc =  TargetStates2el%b(nsti)%nuse(ic)   
!     
!     
!     do in = 1, TargetStates2el%b(nsti)%nam   !get_nam(TargetStates2el%b(nsti))
!        
!        n1 = TargetStates2el%b(nsti)%nb(in) ! get_na(TargetStates2el%b(nsti),in,1)
!        tn1 => bst_nr%b(n1)                                          
!        Lia1 = get_ang_mom(tn1)                            
!        M1 = get_ang_mom_proj(tn1)  
!        
!        n2 =  TargetStates2el%b(nsti)%na(in) ! nc ! get_na(TargetStates2el%b(nsti),in,2)
!        tn2 => bst_nr%b(n2)                                         
!        
!        Lia2 = get_ang_mom(tn2)                            
!        M2 = get_ang_mom_proj(tn2)            
!        CI = get_CI(TargetStates2el%b(nsti),in)
!        IF(CI==0d0) cycle
!        const = CI * sqrt(2.d0) 
!        la = max0(Lia1,Lia2); fla = float(la); mla = Mtot-Mi; fmla=float(mla)

!     print*,'in posvmat:', nqmax,Ps%igpm,Ps%iglp,ki,kf

 
        DO illi = 1,Ps%icount_l(nsti)
           la = Ps%l_nsti(nsti,illi);  fla = float(la); mla = Mtot-Mi; fmla=float(mla)
        const = 1d0
        c0 = 1.d0/pi * (-1)**(Lf)*hat(Li)*sqrfct(2*la+1)*hat(llf)*sqrfct(2*llf+1)*(-1)**((la-llf+Li-Lf)/2)
   
              f0z(:)=0.d0
              if(iq.eq.1) f1z(illi,:,:)=0.d0
              resPs0(:)=0.d0
              resPs1(:)=0.d0
                
              resDnlp(:)=0.d0
              do iz=1,igz
                 xz=xgz(iz)
                 pa2 =qb2 + qa2 - 2.d0*qbqa*xz
                  
                 nn=INT(Log(pa2*ddp)/pdp)
                 if(nn.le.2.or.nn.gt.iglp-2) cycle
                 
                 IF(pa2.gt.gpmin.and.pa2.lt.gpmax) THEN
                    yy(1:intp)=Bnlp(la,nn-INT(nd1):nn+INT(nd1))
                    xx(1:intp)=gp(nn-INT(nd1):nn+INT(nd1))
                    call intrpl(intp,xx,yy,1,pa2,ResBnlp(iz))
                    
                    yy(1:intp)=Dnlp(la,nn-INT(nd1):nn+INT(nd1))
                    call intrpl(intp,xx,yy,1,pa2,ResDnlp(iz))    
                 ENDIF
                 
                 IF(pa2.le.gpmin) THEN
                    ResBnlp(iz)=Bnlp(la,nnmin)
                    ResDnlp(iz)=Dnlp(la,nnmin)
                 ELSEIF(pa2.ge.gpmax) THEN
                    ResBnlp(iz)=0d0 !Bnlp(nnmax1,ic,in)*gp(nnmax1)/pa2
                    ResDnlp(iz)=0d0 !Dnlp(nnmax1,ic,in)*gp(nnmax1)/pa2
                 ENDIF

              enddo
              
              do iz=1,igz 
                 xz=xgz(iz)
                 pb2 = 0.25d0*qb2 + qa2 - qbqa*xz
                 
!!!        call  pswfp(n_ps,llf,pb2,resPs0old,resPs1old) ! this only works upto npos = 2, l_Ps=1 
!!!              print*,'old:', resPs0old,resPs1old

        call f0zpart(Nlf,rlam2,2.d0,nfa,llf,nfa,pb2,resPs0(iz),resPs1(iz),pc0,pc1)
!!!        IF(abs(resPs0old-resPs0(iz))/abs(resPs0old).gt.0.0000001) print*,'new.vs.old:', nfa,llf,resPs0(iz),resPs0old
                 efactor = 0.25d0*qb2 + pb2 - Etot
                 temp = (efactor*resPs1(iz)-resPs0(iz))*resBnlp(iz) + resPs1(iz)*resDnlp(iz) 
                 f0z(iz)=f0z(iz)+temp
              enddo
!       print*,'before f1z calcs.', iq,igz

              if(iq.eq.1) then
                 do ip=1+iqs,2*nqim*igp,iqstep 
                    pp=xp(ip)
                    if(pp.gt.100.) cycle
                    qbpp=qb*pp
                    pp2=pp*pp
                    ResBnlp(:) = 0.d0
                    do iz=1,igz
                       xz=xgz(iz)
                       pa2 =qb2 + pp2 - 2.d0*qbpp*xz
                       nn=INT(Log(pa2*ddp)/pdp)
                       if(nn.le.2.or.nn.gt.iglp-2)  cycle
                       IF(pa2.gt.gpmin.and.pa2.lt.gpmax) THEN
                          yy(1:intp)=Bnlp(la,nn-INT(nd1):nn+INT(nd1))
                          xx(1:intp)=gp(nn-INT(nd1):nn+INT(nd1))
                          call intrpl(intp,xx,yy,1,pa2,ResBnlp(iz))
                       ENDIF
                       IF(pa2.le.gpmin) THEN
                          ResBnlp(iz)=Bnlp(la,nnmin)
                       ELSEIF(pa2.ge.gpmax) THEN
                          ResBnlp(iz)=0.d0 
                       ENDIF
!                    enddo
!                    do iz=1,igz!!
!                       xz=xgz(iz)
                       pb2=0.25*qb2 + pp**2.d0 - qbpp*xz
!!!                       call pswfp(n_ps,llf,pb2,resPs0old,resPs1old) 
                       call f0zpart(Nlf,rlam2,2.d0,nfa,llf,nfa,pb2,resPs0f,resPs1f,pc0,pc1)
!!!        IF(abs(resPs0old-resPs0f)/abs(resPs0old).gt.0.0000001) print*,'new.vs.old:', nfa,llf,resPs0(iz),resPs0old
 
                       f1z(illi,iz,ip)=resPs1f*resBnlp(iz)
                       f1z(illi,iz,ip-1) =   f1z(illi,iz,ip)
                       f1z(illi,iz,ip+1) =   f1z(illi,iz,ip)
                    enddo
                 enddo
!!!!
                 do ip=2*nqim*igp+1,2*nqim*igp + nqim   
                    pp=xp(ip)
                    if(pp.gt.100.) cycle
                    qbpp=qb*pp
                    pp2=pp*pp
                    ResBnlp(:) = 0.d0
                    do iz=1,igz
                       xz=xgz(iz)
                       pa2 =qb2 + pp2 - 2.d0*qbpp*xz
                       nn=INT(Log(pa2*ddp)/pdp)
                       if(nn.le.2.or.nn.gt.iglp-2)  cycle
                       IF(pa2.gt.gpmin.and.pa2.lt.gpmax) THEN
                          yy(1:intp)=Bnlp(la,nn-INT(nd1):nn+INT(nd1))
                          xx(1:intp)=gp(nn-INT(nd1):nn+INT(nd1))
                          call intrpl(intp,xx,yy,1,pa2,ResBnlp(iz))
                       ENDIF
                       IF(pa2.le.gpmin) THEN
                          ResBnlp(iz)=Bnlp(la,nnmin)
                       ELSEIF(pa2.ge.gpmax) THEN
                          ResBnlp(iz)=0d0
                       ENDIF
!                    enddo
!                    do iz=1,igz!!
!                       xz=xgz(iz)
                       pb2=0.25*qb2 + pp**2.d0 - qbpp*xz
                       call pswfp(n_ps,llf,pb2,resPs0f,resPs1f)
                       f1z(illi,iz,ip)=resPs1f*resBnlp(iz)
                    enddo
                 enddo
              endif
              
              icalam(:) = 0
              res0(:) = 0.d0

              itest(:,:)=0
              
              res1_nu_lam(:,:) = 0d0

              reslb1 = 0.d0

!       print*,'Start angular mom. sums:', la,llf 

              do lb1=0, llf
                 flb1 = float(lb1)
                 lb2 = llf - lb1
                 flb2 = float(lb2)
                 Clb1 = 1.d0 / sqrfct(2*lb1) / sqrfct(2*lb2) / 2d0**lb1 /hat(lb2)
                 
                 resla1 = 0.d0
                 do la1=0, la
                    fla1 = float(la1)
                    la2 = la - la1
                    fla2 = float(la2)
                    
                    lab2 = lb2+la2 
                    
                    Cla1 = qb**(la1+lb1) / sqrfct(2*la1) / sqrfct(2*la2)/hat(la1) 
                    
                    resl1 = 0.d0
                    do l1 = iabs(la1-lb1),la1+lb1,2
                       fl1 = float(l1)
                       Cl1 = hat(l1) * cgc0(fl1,flb1,fla1)  !!!/hat(l1)
                       
                       resl2 = 0.d0
                       do l2 = iabs(lb2-la2), lab2,2  
                          fl2 = float(l2)
                          Cl2 = hat(l2) *  cgc0(fl2,fla2,flb2) 
                          
                          
                          reslam = 0.d0
                          do lam=max0(iabs(Lf-l1),iabs(Li-l2)), min0(Lf+l1,Li+l2),2
                             flam=float(lam)
                             
                             Clam = (-1)**lam * (2*lam+1) *  cgc0(rLf,flam,fl1) * cgc0(rLi,flam,fl2)
                             
                             reska = 0.d0
!                             do ka = imax0(iabs(Li-Lf),iabs(llf-la)), min0(Lf+Li,llf+la)
                             do ka = max0(iabs(Li-Lf),iabs(llf-la)), min0(Lf+Li,llf+la)
                                fka = float(ka)
                                Cka =(-1)**ka * (2*ka+1) * cof6j(fl2,flam,rLi,rLf,fka,fl1)
                                Cka = Cka * cof9j(fla, fla1,fla2,fllf,flb1,flb2,fka,fl1,fl2)
                                if(Cka==0d0) cycle
                                
                                kam = Mf - Mi
                                fkam = float(kam)
                                
                                Ckam =  cgc(rLi,rMi,fka,fkam,rLf,rMf) 
                                Ckam = Ckam * cgc(fllf,fmllf,fka,fkam,fla,fmla)
                                
                                reska = reska + Cka * Ckam 
                             enddo
                             
                             Clam = Clam  * reska
                             
                             if(Clam == 0d0) cycle
                             
                             if(icalam(lam).ne.1) then
                                
                                resiz=0.d0
                                do iz=1,igz
                                   resiz=resiz+f0z(iz)*pl(lam,iz)
                                enddo
                                res0(lam) = resiz
!!!!!
                                if(iq.eq.1) then
                                   fpqb1(illi,lam,:)=0.d0
                                   do ipl=1+iqs,2*nqim*igp,iqstep !imax
                                      pp=xp(ipl)
                                      resi=0.d0
                                      if(pp.gt.100.d0) cycle
                                      do iz=1,igz
                                         resi=resi+f1z(illi,iz,ipl)*pl(lam,iz)
                                      enddo
                                      fpqb1(illi,lam,ipl)=resi       
                                      fpqb1(illi,lam,ipl-1)=resi       
                                      fpqb1(illi,lam,ipl+1)=resi       
                                   enddo

                                   do ipl=2*nqim*igp+1,2*nqim*igp + nqim !imax
                                      pp=xp(ipl)
                                      resi=0.d0
                                      if(pp.gt.100.d0) cycle
                                      do iz=1,igz
                                         resi=resi+f1z(illi,iz,ipl)*pl(lam,iz)
                                      enddo
                                      fpqb1(illi,lam,ipl)=resi       
                                   enddo

                                endif
                                   fpqb(lam,:) =fpqb1(illi,lam,:) 
                                   icalam(lam) = 1
                                endif
!!!!!!
!      GO TO 15
! Integration over composite mesh:
!c             if(gridk(1,nchi).gt.0.d0) then! Rav: subtraction
!                           method must be used for both open and closed
!                           channels. Some changes have been made to
!                           make it work. 
                            res1 = 0.d0
!         IF(Lia2.eq.-100) THEN

      nu=lab2+1 
      IF(itest(nu,lam).eq.0) THEN       

             ising=ki    ! singular point
             if(ki.le.ionsh) ising=ki-1
             if(ki.eq.1) ising=ionsh
         if(gridk(1,nchi).lt.0.) ising=ki-1
!        print*,'Statrt Q Integr', nqmi,ising,ki,kf,ionsh,nqgen 

            if(ising.ne.nqim) then   !i.e. singularity is not at last k-mesh
            res1=0.d0
!*   integrals coming before singularity
            do i=1+iqs,2*(ising-1)*igp,iqstep
             pp=xp(i)
             fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i)) !-Alp(ki,i))
             res1=res1+fp * dwp(i) 
             enddo
!*   integral with singularities
             res2=0.d0
             do i=2*(ising-1)*igp+1,(2*ising-1)*igp
             pp=xp(i)
             fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i)) &
     &        - fpqb(lam,2*nqim*igp+ising)*qa**(lab2+1)*Q0p(ki,i)
             res2=res2+wp(i)*fp
             enddo

            res1b=2*qa*Log(2.)-qa*Log(qa-aq(ising-1))+2*qa*Log(qa)- &
     &            qa*Log(aq(ising-1)+qa)- &
     &          aq(ising-1)*Log((aq(ising-1)+qa)/(-aq(ising-1)+qa))

           res1b=res1b*fpqb(lam,2*nqim*igp+ising)*qa**(lab2+1)
           res2=res2+res1b
     
           res3=0.d0
           do i=(2*ising-1)*igp+1, 2*ising*igp
           pp=xp(i)
           fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i)) &
     &        - fpqb(lam,2*nqim*igp+ising)*qa**(lab2+1)*Q0p(ki,i) 
           res3=res3+wp(i)*fp
           enddo

           res1c=-2*qa*Log(2.)+qa*Log(aq(ising)-qa)- &
     &           2*qa*Log(qa)+qa*Log(aq(ising)+qa) &
     &          +aq(ising)*Log((aq(ising)+qa)/(aq(ising)-qa))
           res1c=res1c*fpqb(lam,2*nqim*igp+ising)*qa**(lab2+1) 
           res3=res3+res1c
!c           print*,'res3:', res3,res1c    
!* Integrals coming after singularity
          res4=0.d0
          do i=2*ising*igp+1+iqs,2*nqim*igp, iqstep
             pp=xp(i)
             fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i)) 
             res4=res4+fp * dwp(i)
          enddo
!c          print*,'RES::', res1,res2,res3,res4
          res1=res1+res4+res2+res3
          else
!*  singularity at last k-mesh point
          res1=0.d0
          do i=1+iqs,2*(nqim-1)*igp, iqstep
             pp=xp(i)
             fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i))
             res1=res1+fp * dwp(i)
             enddo

!*   integ with singularities
            res2=0.d0
          do i=2*(nqim-1)*igp+1,(2*nqim-1)*igp
            pp=xp(i)
            fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i)) &
     &        -fpqb(lam,2*nqim*igp+nqim)*qa**(lab2+1)*Q0p(ki,i)
           res2=res2+wp(i)*fp

          enddo

        res1b=2*qa*Log(2.)-qa*Log(qa-aq(nqim-1)) &
     &        +2*qa*Log(qa)-qa*Log(aq(nqim-1)+qa) &
     &        -aq(nqim-1)*Log((aq(nqim-1)+qa)/(-aq(nqim-1)+qa))
        res1b=res1b*fpqb(lam,2*nqim*igp+nqim)*qa**(lab2+1)      ! check lab2+1??
        res2=res2+res1b

            res3=0.d0
            do i=(2*nqim-1)*igp+1,2*nqim*igp
               pp=xp(i)
           fp=fpqb(lam,i)*pp**(lab2+1)*(Qlp(ki,i)) &
     &      -fpqb(lam,2*nqim*igp+nqim)*qa**(lab2+1)*(qa2+alfa2)**2.d0/ &
     &         (pp*pp+alfa2)/(pp*pp+alfa2)*Q0p(ki,i)
           res3=res3+wp(i)*fp
                      enddo

          res1c=(-2*qa*ATan(qa/alfa)+Pi*qa)/(2.*alfa*(alfa2+qa2))
          res1c=res1c*fpqb(lam,2*nqim*igp+nqim)* &
     &           qa**(lab2+1)*(qa2+alfa2)*(qa2+alfa2)
                   ! CHECK qa**lab2 or **lab2+1
          res3=res3+res1c

          res1=res1 +res2+res3
          endif
!           else
!
!*i.e. on-shell point is negative
!
!        res1=0.d0
!c        print*,'max',2*nqim*igp,imax,nqim,igp
!       do i=1,2*nqim*igp
!           pp=xp(i)
!           fp=fpqb(lam,i) *pp**(lab2+1)*(Qlp(ki,i))  !-Alp(ki,i))
!           res1=res1+wp(i)*fp
!          enddo
!        ENDIF
! 15  CONTINUE
        resAlp=0.d0

!***  this additional term with Alp doesn't contain singularity 1/q^2

          do i=1+iqs,2*nqim*igp, iqstep
            pp=xp(i)
            fp=fpqb(lam,i) *pp**(lab2+1)*Alp(ki,i)  
            resAlp=resAlp+fp * dwp(i)
          enddo

         res1 = res1 + resAlp

         res1_nu_lam(nu,lam)= res1
         itest(nu,lam) = 1
      ENDIF ! itest(nu,lam)==0       
         res1 = res1_nu_lam(nu,lam)

!!  Integration over q-mesh is done!

                    res12 = res1/pi + res0(lam) * qa**(lab2+1)
                    reslam=reslam + Clam * res12
!!!**********
                       enddo !lam

!!*************
                       resl2 = resl2 + reslam * Cl2
                    enddo !l2
                    
                    resl1 = resl1 + resl2 * Cl1
                 enddo !l1
                 
                 resla1 = resla1 + resl1 * Cla1
              enddo !la1
              
              reslb1 = reslb1 + resla1 * Clb1
           enddo !lb1
           
           result =  result + const * reslb1*c0*qb*sqrt(2.d0)
!     if(la==1) write(*,'(a7,7i4,a3,2i4,e20.10)') 'vmatt:',la,ma,Li,Mi,Lf,Mf,Mtot,' $$$',kf,ki,result 
           
!  enddo !in     !n_orb
!enddo !ic
      ENDDO !illi



           IF(h2tops) then
              vmatt(kf,ki) = vmatt(kf,ki) + result
           ELSE
              vmatt(ki,kf) = vmatt(ki,kf) + result
           ENDIF
        enddo !ki
     enddo  ! kf
     



RETURN
END Subroutine psformation


subroutine pswfp(n,l,p2,res0,res1)
integer, intent(in):: n,l
real*8, intent(in)::p2
real*8, intent(out)::res0,res1
real*8:: pp,ppll

  pp = dsqrt(p2)
  ppll = pp**l

if(n.eq.1.and.l.eq.0) then
  res0 = dsqrt(8.d0)/(4*p2+1.d0) / ppll
  res1 = dsqrt(128.d0)/(4*p2+1.d0)**2.d0 / ppll
endif

if(n.eq.2.and.l.eq.0) then
  res0 = (64*p2-4)/(16*p2+1.d0)**2/ ppll
  res1 = (64*(16*p2-1))/(16*p2+1.d0)**3.d0 / ppll
endif

if(n.eq.2.and.l.eq.1) then
  res0 = 32.d0*pp/sqrt(3.d0)/(16*p2+1.d0)**2/ ppll
  res1 = 512.d0*pp/sqrt(3.d0)/(16d0*p2+1d0)**3/ppll
endif



if(n.eq.3.and.l.eq.0) then
  res0 = 2d0*sqrt(6d0)*(1296d0*pp**5-120d0*pp**3+pp)/(36*p2+1.d0)**3 / ppll
  res1 = 72d0*sqrt(6d0)*(1296d0*pp**5-120d0*pp**3+pp)/(36*p2+1.d0)**4 / ppll
endif

if(n.gt.3.or.l.gt.1) STOP 'Ps-wavefunctions have not yet been coded for cases l>1 or n>3'

return
end subroutine


!c*****************************************
!C this is for Ps-Ps direct V m.e.
!!  uses fd0 for R < r_{max} of (wafe function)    
!c*****************************************
  subroutine pspsdme(Mtot,nsti,nstf,nchi,nchf,nqmi,nqmf,npk_nchi,npk_nchf,chil,nr,weight,gridr,vmatt)

  use sturmian_class
  use channels
  use input_data
  use kgrid 
  use Positronium
  use grid_radial
   
  
  implicit none
  integer, parameter:: maxl=2*20+1
  integer, intent(in):: Mtot,nr
  integer, intent(in):: nchf, nchi        ! Channel number
  integer, intent(in):: nstf, nsti        ! State for particular channel
  integer, intent(in):: nqmf, nqmi        ! Number of k-grid points for particular channel
  integer, intent(in):: npk_nchi, npk_nchf ! On-shell k-grid point for particular channel
  real*8, dimension(nqmf,nqmi), intent(inout):: vmatt
  type(basis_sturmian_nr):: chil
  integer:: ki, kf, kqi, kqf
  real*8:: pi
  real*8:: xgz(Ps%igzm),wgz(Ps%igzm),pl(0:20,Ps%igzm)
  real*8, dimension(0:maxl):: factrl,sqrfct,hat
  real*8,dimension(nr), intent(in):: weight,gridr
  real*8,dimension(nr):: fun, temp,fun0,ttt
  type(sturmian_nr), pointer:: tn0,tnp0
  real*8, pointer, dimension(:)::  fp0, f0
  real*8:: sum_rad,const,rnorm,cgc0,cgc
  integer:: Li,Mi,lia,mlia,Lf,Mf,lfa,mlfa,lam,mlam,minfun,maxfun,ir1,ir2,minf,maxf,i1,i2,ifdstop,i,minfp,maxfp
  real*8::  rLi,rMi,flia,fmlia,rLf,rMf,flfa,fmlfa,flam,fmlam,rMtot
  real*8, dimension(nr):: xin,yin,xout,yout

  pi = acos(-1d0)
  rnorm = 2.d0 / pi 
  
!   RETURN
  
 
!    weight(1:nr) = grid%weight(1:nr)
!    nr = grid%nr
!    rmesh(1:nr) = grid%gridr(1:nr)


  Li = Lp_ch(nchi); Mi=Mp_ch(nchi); Lf = Lp_ch(nchf);Mf=Mp_ch(nchf);
  rMtot=float(Mtot); rLi=float(Li);rMi=float(Mi); rLf=float(Lf); rMf=float(Mf);
  lia = Ps%l(nsti);flia=float(lia);  mlia = Mtot - Mi; fmlia = float(mlia)
  lfa = Ps%l(nstf);flfa=float(lfa);  mlfa = Mtot - Mf; fmlfa = float(mlfa)
  mlam = Mi - Mf ! == mlfa - mlia
  fmlam = float(mlam)


   IF(((-1)**(lia+lfa).eq.1).or.((-1)**(Li+Lf).eq.1)) then
!         write(17,*)
!         write(17,*)'# lam cycle:',lia+lfa,Li+Lf
!         write(17,'(a8,3i4,a8,3i4)') '# nsti=',nsti,Li,lia,'nstf=',nstf,lf,lfa
         return
   ENDIF

   IF((-1)**(Li+lia).ne.(-1)**(Lf+lfa))  then
!         write(17,*)
!         write(17,*)'# parity is not the same for i&f'
!         write(17,*) '# nsti=',nsti,li,lia,'=>','nstf=',nstf,lf,lfa
         return
   ENDIF

  factrl(0:maxl)=Ps%factrl(0:maxl); sqrfct(0:maxl)=Ps%sqrfct(0:maxl); hat(0:maxl)=Ps%hat(0:maxl);
  

      minfun = 1
      maxfun = min(Ps%maxf(nsti),Ps%maxf(nsti))
      ifdstop = maxfun
      do i = minfun, maxfun
         fun(i) = Ps%f(nsti,i) * Ps%f(nstf,i) * weight(i)
      end do 

      do i = 1,nr
         if(gridr(i).le.20.d0) ifdstop = i
      end do 

      ifdstop=min(ifdstop,maxfun)

      IF(ifdstop.gt.maxfun) THEN
       PRINT*, '',ifdstop,maxfun
       STOP 'check ifdstop'
      ENDIF



     ttt(1:nr) = 0.d0

      DO lam = max0(iabs(lia-lfa),iabs(Li-Lf)), min0(lia+lfa,Li+Lf), 2
         flam = float(lam) 
         temp(:)=0.0
         
         if((-1)**lam.eq.1) CYCLE
             
         const =  cgc0(rLf,flam,rLi) *  cgc0(flia,flam,flfa)
         const = const * cgc(rLf,rMf,flam,fmlam,rLi,rMi) * cgc(flia,fmlia,flam,fmlam,flfa,fmlfa)     
!         if(const.eq.0.d0) cycle
         const = const*(-1)**(Li+Lf+lia+lfa)
         const = const * hat(lf) * hat(lia) / hat(Li) / hat(lfa)
         const = 2.0 * const * (1-(-1)**lam)

         temp(1:nr) = 0.d0 

         call form(lam,fun,minfun,maxfun,nr,temp,i1,i2)

             do i = i1, i2
               xin(i-i1+1) = gridr(i)
               yin(i-i1+1) = temp(i) * xin(i-i1+1) ** (lam+1)
               xout(i-i1+1) = gridr(i) * 2d0
            enddo
            if (i2-i1+1.gt.nr) then
               print*,'I2,I1,LT:',i2,i1,lam
               stop 'problem in call to intrpl'
            endif
            call intrpl(i2-i1+1,xin,yin,i2-i1+1,xout,yout)
            do i = i1, i2
               if (xout(i-i1+1).gt.xin(i2-i1+1))  yout(i-i1+1) = yin(i2-i1+1)
               temp(i) = 2.0 * (yout(i-i1+1) / xout(i-i1+1)**(lam+1))
            enddo
!         write(17,*) 
!         write(17,'(a7,5i6)'),'# lam=',lam,nsti,nstf,i1,i2
!         do i =1,nr
!         write(17,'(i6,2e20.10)')i, gridr(i),temp(i)*(1-(-1)**lam)
!         enddo
!         STOP


        do i = i1,i2
        ttt(i) = const * temp(i) + ttt(i)   
        enddo

!         write(17,*) 
!         write(17,'(a7,5i6,e20.10)'),'# lam=',lam,nsti,nstf,i1,i2,const
!         do i =1,i2
!         write(17,'(i6,2e20.10)')i, gridr(i),ttt(i)
!         enddo

     ENDDO !lam



  do ki=1,nqmi
     kqi = npk_nchi + ki - 1
     tn0 => chil%b(kqi)
     f0 => fpointer(tn0)
     minf = get_minf(tn0)
     maxf = get_maxf(tn0)

!     write(19,*) 
!     write(19,*)'#', f0(1)     

     sum_rad=0d0

     do kf=1,nqmf
        kqf = npk_nchf + kf - 1
        tnp0 => chil%b(kqf)
        fp0 => fpointer(tnp0)
        minfp = get_minf(tnp0)
        maxfp = get_maxf(tnp0)
        
        ir2 = min(maxf,maxfp,i2) ! nr as for continuum functions, but for closed channels maxf = 1
        ir1 = max(minf,minfp,i1)
        
        fun0(ir1:ir2) = (f0(ir1:ir2) * fp0(ir1:ir2)) * weight(ir1:ir2) * ttt(ir1:ir2)  

        sum_rad = SUM( fun0(ir1:ir2) ) 

        vmatt(kf,ki) = vmatt(kf,ki) + sum_rad * rnorm

!        write(19,'(2i5,3e20.10,2i6)') kf,ki, vmatt(kf,ki),ttt(100),rnorm,ir1,ir2

     end do ! end kf loop

  end do ! end ki loop

      return 
      end subroutine
    


       
!!!! Rav: below subroutine added for Ps-channels


subroutine makevstat

  use input_data
  use channels
  use target_states
  use grid_radial
  use sturmian_class
  use one_electron_func
  use Positronium
    

real*8, dimension(grid%nr):: fun,formf,gridr,weight,temp 
! real*8, dimension(gridr%nr,)
real*8:: formf0,rr,V0,Rd, ve2s_gp(Ps%iglp), gp(Ps%iglp)
type(sturmian_nr), pointer:: tn2  ! One-electron orbitals
real*8, pointer, dimension(:):: f2, f2p
integer i,j, nr,i1,i2,ni2,li2,Ncon,ne_con,ind2,minf,maxf,mingp,maxgp,ip,iglp
real*8:: xx(5),yy(5),y, x,r1,r2,rr2
integer:: ig, igz, lam


     Rd = data_in%Rd  
     iglp = Ps%iglp; gp(1:iglp) = Ps%gp(1:iglp);      ddp=1.e+08; pdp=0.0065d0;

       fun(:) = 0.d0
       formf(:) = 0.d0  
       nr = grid%nr
          ni2 = 1
          li2 = 0
         nsti = 1 
         Ncon = TargetStates2el%b(nsti)%nam         
         ne_con = 1  !, Ncon
         ind2 = TargetStates2el%b(nsti)%nb(ne_con) 
          tn2 => bst_nr%b(ind2)                                         
           f2 => fpointer(tn2)                               
          minf = get_minf(tn2)
          maxf = get_maxf(tn2)

     fun(minf:maxf) = f2(minf:maxf)*f2(minf:maxf) * grid%weight(minf:maxf)


        call form(0,fun,minf,maxf,nr,formf,i1,i2)

            formf0 = formf(1)

!     write(18,*)'#',i1,i2

      do i = 1, i2
        rr=grid%gridr(i) 
       IF(rr.ge.Rd/2d0) THEN 
         V0 = 2d0/rr
       ELSE
         V0 = 4d0/Rd
       ENDIF
        Ps%ve2st(i) = (V0 - formf(i)) * rr
        Ps%va2st(i) = (-V0 + formf(i)) * rr
!        write(18,'(i5,4e20.10)') i, grid%gridr(i), Ps%ve2st(i), Ps%va2st(i), f2(i)    
      enddo

         maxgp=1
        do ip=1,iglp
        if(gp(ip).ge.grid%gridr(i2-3)) cycle
        if(gp(ip).le.grid%gridr(4)) mingp=ip
         maxgp=ip
        enddo

      ve2s_gp(:)=0.d0

      call intrpl(i2,grid%gridr(1),Ps%ve2st(1),iglp,gp,ve2s_gp)

!!      ve2s_gp(1:mingp) = Ps%ve2st(1)/grid%gridr(1) * gp(1:mingp)
!!      ve2s_gp(maxgp+1:iglp) = Ps%ve2st(i2)/grid%gridr(i2) * gp(maxgp+1:iglp)
!
!      write(19,*)'#', mingp,maxgp
!      do ip=1,iglp
!         
!         write(19,'(2e20.10,i6)') gp(ip),ve2s_gp(ip),ip
!
!      enddo

!      DO nsti = 1, Ps%Nmax
!         DO nstf = 1, Ps%Nmax
!            IF((-1)**(Ps%l(nsti) + Ps%l(nstf)).eq.1) CYCLE
!
!            maxfun = min(Ps%maxf(nsti),Ps%maxf(nsti))
!            do i = minfun, maxfun
!               fun(i) = Ps%f(nsti,i) * Ps%f(nstf,i) * grid%weight(i)
!            end do
!
!
!            DO i = 1, grid%nr
!               ri = gridr(i)
!               DO j = 1, grid%nr
!                  rj = gridr(j)
!
!                  res0=0.0d0                
!                  res1=0.0d0                
!                  r1=max(ri,rj/2.)
!                  r2=min(ri,rj/2.)
!
!
!                  DO ig=1,igz            
!                     z=xgz(ig)
!                     rr2=r1*r1+r2*r2-2*r1*r2*z             
!                     rr=dsqrt(rr2)
!                     
!                     nn=INT(Log(rr*ddp)/pdp)
!                     yy(1:intp)=ve2s_gp(nn-INT(nd1):nn+INT(nd1))
!                     xx(1:intp)=gp(nn-INT(nd1):nn+INT(nd1))
!                     call intrpl(intp,xx,yy,1,rr,y)
!
!                     res1=res1+y/rr*pl(lam,ig)   
!                  END DO ! ig  z-integration
!
!    
!                  res=((2.*lam+1.)/2.)*res1
!                  
!                  ResInt=ResInt+fun(j)*res
!
!               END DO  ! j   \rho - integration
!
!               fd0(i,lam,npsi,npsf) = ResInT  
!
!            END DO ! i R -
!
!         ENDDO ! nstf
!
!      ENDDO ! nsti
!
return 
end subroutine makevstat


subroutine h2ff(Nmax1)

  use input_data
  use channels
  use target_states
  use grid_radial
  use sturmian_class
  use one_electron_func
  use Positronium
  use  state_class 
  use MPI_module

integer, parameter:: maxl=2*20+1  
type(sturmian_nr), pointer:: tn1,tn2,tn2_ion  ! One-electron orbitals
real*8, pointer, dimension(:):: f1,f2
integer:: Ncon,M12,Spin,n_orb_max,nsti,Nmax1,n_orb,nuse_norb,minf,maxf&
& ,mini,maxi,ne_con,ind2,ind1_2,La2,maxf_ion,minf_ion, ne_ion,ip

real*8:: Bnl(0:4,grid%nr),Dnl(0:4,grid%nr) 

real*8:: CI,ort_n2wion,rr,const,hep,r,pp,pp2,ppll,res1,res2,res3,res4,arg,f1sp,besskr,Rd,V0
real*8, dimension(grid%nr):: formf, gridr, weight,fun,Fn2,wfion,formf1,wfion1,fun1,fun2
integer:: lli,iglp,in,n1,n2,i, n2_ion,ic,nc,l_ion, nicm
logical:: done_l(0:4)
integer:: icount1_l
real*8:: flli,fmlli,fLa1,fLa2,fM1,fM2,clli,clli0, cgc,cgc0
real*8, dimension(0:maxl):: factrl,sqrfct,hat


  minf_ion=0; maxf_ion=0
  l_ion=0; n2_ion=0
  Rd = data_in%Rd 
  iglp = Ps%iglp
     
   nicm =   get_nam(TargetStates2el%b(1)) ! 4  !TargetStates2el%nam 
   Ps%nam_max = nicm

   IF(allocated(Ps%Amg))  THEN
      deallocate(Ps%Amg, Ps%Bnlp, Ps%Dnlp)
      deallocate(Ps%l_nsti,Ps%icount_l)
      
   ELSE
      allocate(Ps%Amg(iglp), Ps%Bnlp(Nmax1,0:4,iglp), Ps%Dnlp(Nmax1,0:4,iglp))
      allocate(Ps%l_nsti(Nmax1,5),Ps%icount_l(Nmax1))
   ENDIF

 print*,'nicm, Bnlp_mem:',nicm, Nmax1,float(5*Nmax1*iglp)*8.d0*2d0/1000000.d0,'Mb' 
 Ps%Bnlp(:,:,:) = 0d0; 
 Ps%Dnlp(:,:,:) = 0d0; 
 Ps%l_nsti(:,:) = 0
 Ps%icount_l(:) = 0 

 nr = grid%nr
 gridr(1:nr)=grid%gridr(1:nr)
 weight(1:nr)=grid%weight(1:nr)

 hat(0:maxl)=Ps%hat(0:maxl);

 nsti = 1
 ne_ion = 1
 ne_con = ne_ion                   
 ind2 = TargetStates2el%b(nsti)%nb(ne_con) 

 nam = get_nam(TargetStates2el%b(nsti))
 do i=1,1 !nam
    n1 = get_na(TargetStates2el%b(nsti),i,1)  
    tn2 => bst_nr%b(n1)                                         
    f2 => fpointer(tn2)    
    minf_ion = get_minf(tn2)
    maxf_ion = get_maxf(tn2)
    wfion1(1:maxf_ion) = f2(1:maxf_ion)
    wfion1(maxf_ion+1:nr) = 0.d0
    l_ion = get_ang_mom(tn2)                            
  
    n2 = get_na(TargetStates2el%b(nsti),i,1) 
    tn2 => bst_nr%b(n2)                                         
    tn2_ion => bst_nr%b(n2)                                         
    f2 => fpointer(tn2)    
  
 enddo 

!   tn2 => bst_nr%b(ind2)                                         
!   f2 => fpointer(tn2)                               
 n2_ion = n2
 minf_ion = get_minf(tn2)
 maxf_ion = get_maxf(tn2)
 wfion(1:maxf_ion) = f2(1:maxf_ion)
 wfion(maxf_ion+1:nr) = 0.d0
    
!!!RETURN

 PRINT*,'start Bnlp calcs:'
!   !Dont use "OMP private" if there are derived or allocated data inside the OMP loop

 
!$OMP PARALLEL DO &
!$OMP SCHEDULE(dynamic)&
!$OMP DEFAULT(none)&
!$OMP private(nsti,Bnl,Dnl,Fn2,done_l,icount1_l)&
!$OMP private(ic,nc,in,n1,n2,tn1,tn2,f1,f2,mini,maxi,fun1,fun2,minf,maxf)&
!$OMP private(La1,La2,M1,M2,Ci,formf,formf1,minfun,maxfun)&
!$OMP private(fun,i,r,hep,i1,i2,ort_n2wion,rr,V0,lli,const)&
!$OMP private(illi,ip,pp2,pp,ppll,res1,res2,arg,besskr,res4,f1sp)&
!$OMP shared(Nmax1,minf_ion,maxf_ion,n2_ion)&
!$OMP shared(wfion,l_ion)&
!$OMP shared(iglp)&
!$OMP shared(nr,gridr,weight,Rd)&
!$OMP shared(TargetStates2el, bst_nr,Ps,hat)&
!$OMP private(flli,fmlli,fLa1,fLa2,fM1,fM2,clli,clli0) 
     DO  nsti = 1, Nmax1 
          
   Bnl(:,1:nr)=0.d0; Dnl(:,1:nr)=0.d0; Fn2(1:nr)=0.d0;
   done_l(:)=.false.
   icount1_l = 0
!   print*,'Nmax1,nsti:',Nmax1,nsti

     do ic = 1, TargetStates2el%b(nsti)%nusemax 
        
        nc =  TargetStates2el%b(nsti)%nuse(ic)   


        do in = 1, TargetStates2el%b(nsti)%nam   !get_nam(TargetStates2el%b(nsti))

          n1 = TargetStates2el%b(nsti)%nb(in) ! get_na(TargetStates2el%b(nsti),in,1)

           tn1 => bst_nr%b(n1)                                           
           f1 => fpointer(tn1)                               
           mini = get_minf(tn1)
           maxi = get_maxf(tn1)

             fun1(1:maxi) = f1(1:maxi)
             fun1(maxi+1:nr) = 0d0

           La1 = get_ang_mom(tn1)                            
           M1 = get_ang_mom_proj(tn1)  

             n2 =  TargetStates2el%b(nsti)%na(in) ! nc ! get_na(TargetStates2el%b(nsti),in,2)
             tn2 => bst_nr%b(n2)                                         
             f2 => fpointer(tn2)                               
             minf = get_minf(tn2)
             maxf = get_maxf(tn2)

             fun2(1:maxf) = f2(1:maxf)
             fun2(maxf+1:nr) = 0d0

             La2 = get_ang_mom(tn2)                            
             M2 = get_ang_mom_proj(tn2)            
             CI = get_CI(TargetStates2el%b(nsti),in)
             IF(CI==0d0) cycle
             const = CI * sqrt(2.d0)

       fun(:)=0.  
       formf(:)=0.

          minfun = max(minf,minf_ion)
          maxfun = min(maxf_ion,maxf)
       do i = minfun, maxfun
          r=gridr(i)
          hep=wfion(i) !!  ion w.f.   For He: 4.*sqrt(2.)*exp(-2.*r)
          fun(i)=hep*fun2(i)*weight(i)/(2*La2+1)
       enddo
      
      call form(La2,fun,minfun,maxfun,nr,formf,i1,i2)
   IF((n2 /= nc)) CYCLE
!   IF(La2*La1.gt.0d0) CYCLE !Assumes frozen core (or only S-states for the core orbital)
!!!   Ps%l_of_t(nsti)=max0(La1,La2)


!        write(17,*)
!       write(17,'(3i5)') nsti,ic,in
!       write(17,'(a2,3i5)')'a:',TargetStates2el%b(nsti)%nusemax,n1,La1
!       write(17,'(a2,3i5)') 'b:',TargetStates2el%b(nsti)%nam,n2,La2

!       ort_n2wion = 0.d0
!       if(n2.eq.n2_ion) ort_n2wion = 1.d0
        
       ort_n2wion = SUM(f2(minfun:maxfun) *  wfion(minfun:maxfun) * weight(minfun:maxfun))
       IF(La2.ne.l_ion) ort_n2wion = 0d0 ! l_ion ==0 if FC

   formf(i2+1:nr) = ort_n2wion / gridr(i2+1:nr)
   IF(n2.eq.n2_ion) formf1(:)=formf(:)
!   write(17,'(a1,3i5,2e20.10,2i6)')'#', n1,n2, TargetStates2el%b(nsti)%nusemax, CI, ort_n2wion,i1,i2

 
!      write(17,'(a1,3i5)'),'#',La1,La2,nsti
!      write(17,'(a1,3i5)'),'#',M1,M2,M1+M2

    fLa1=float(La1); fLa2=float(La2); fM1=float(M1); fM2=float(M2)
    do lli=abs(La1-La2), La1+La2
      flli=float(lli); fmlli=float(M1+M2)
      if(lli.lt.M1+M2) cycle
      clli= cgc(fLa1,fM1,fLa2,fM2,flli,fmlli)
      clli0=cgc0(fLa1,fLa2,flli)
      if(clli*clli0.eq.0.d0) cycle
      const=const*clli*clli0*hat(La1)*hat(La2)/hat(lli)

      do i=1, nr
        rr=gridr(i) 

       IF(rr.ge.Rd/2d0) THEN 
         V0 = 2d0/rr
       ELSE
         V0 = 4d0/Rd
       ENDIF

!      lli = max0(La1,La2) !Ps%l_of_t(nsti)

      Fn2(i) = -V0 * ort_n2wion + formf(i) 
      Bnl(lli,i) = Bnl(lli,i) + const * ort_n2wion * fun1(i)
      Dnl(lli,i) =  Dnl(lli,i) + const * Fn2(i) * fun1(i) 
!      write(17,'(i5,5e20.10)') i,rr,Fn2(i),Bnl(i),Dnl(i),formf(i)
      enddo
       if(.not.done_l(lli)) then
          icount1_l = icount1_l + 1
          done_l(lli)=.true.
          Ps%l_nsti(nsti,icount1_l) = lli
       endif

     enddo !lli

                        
      enddo !in     !n_orb
     enddo !ic

!            enddo !ne_con
  

!      write(17,'(a1,6i5)'),'#',n2_ion,get_nam(TargetStates2el%b(nsti)), get_na(TargetStates2el%b(nsti),1,1),&
! &   get_na(TargetStates2el%b(nsti),1,2),  get_na(TargetStates2el%b(nsti),2,1), get_na(TargetStates2el%b(nsti),2,2)
!     write(17,*)

!!       lli = max0(La1,La2) !Ps%l_of_t(nsti)
      
!      write(17,'(a1,2i5)'),'#',nsti,lli
!       write(17,*)
!      do i=1, nr 
!        rr=gridr(i) 
!      write(17,'(i5,5e20.10)') i,rr,Bnl(i),Dnl(i),wfion(i),formf1(i)
!      enddo
      Ps%icount_l(nsti) = icount1_l  
!       write(18,*)
 
      DO illi=1,icount1_l
         lli = Ps%l_nsti(nsti,illi)

!        write(18,'(3i5)') nsti, illi,lli
             

      DO ip=1,iglp  ! dense p2-grid for interpolation

         pp2=Ps%gp(ip)
         pp=sqrt(pp2)
         ppll = pp**lli

        res1=0.
        res2=0.

        do i=1,nr 
         rr=gridr(i)
         arg = rr * pp

        call sbessel(arg,lli,besskr) 

         rr = rr*weight(i)
         res1 = res1 + besskr*Bnl(lli,i)*rr
         res2 = res2 + besskr*Dnl(lli,i)*rr
       enddo

        Ps%Bnlp(nsti,lli,ip) = res1/ppll
        Ps%Dnlp(nsti,lli,ip) = res2/ppll

       ENDDO  
      ENDDO
!!**********************************************************8

      IF(nsti.eq.1) then

       Ps%Amg(:)=0.

       DO ip=1,iglp ! dense p2-grid for interpolation
 
         pp2=Ps%gp(ip)
         pp=sqrt(pp2)

        res4=0.0

        do i=1,nr
         rr=gridr(i)
         arg = rr * pp

       IF(rr.ge.Rd/2d0) THEN 
         V0 = 1d0/rr
       ELSE
         V0 = 2d0/Rd
       ENDIF


        f1sp = (V0 - formf1(i)) * rr
  
        res4=res4+ f1sp*dsin(pp*rr)/pp*weight(i)
       enddo

        Ps%Amg(ip)=res4
        ENDDO        
       ENDIF

   ENDDO
!$OMP END PARALLEL DO
 
  
RETURN
END SUBROUTINE h2ff

  
!*---------------------------------------------------

      subroutine sbessel(X, L, BES)
      implicit real*8 (a-h,o-z)
      implicit integer(i-n) 
!     modifyed spherical Bessel function 
      if (L .lt. 0) then                       
         BES = 0d0                           
         return                                
      end if
      
      if (X .eq. 0d0) then
         BES = 0d0 
         if (L .eq. 0)  BES = 1d0          
         return
      end if
      
      select case(L)
      case(0)
         if (x.lt.1.e-2) then            
            x2 = x*x
            BES = 1.0 - (x2*(1.0 - ((1.0 - x2/42.)*x2)/20.))/6.
         else
            BES = sin(x)/x
         end if
      case(1)
         if (x.lt.1.e-2) then            
            x2 = x*x
            BES = (x*(1.-(x2*(1.-(x2*(1.-x2/54.))/28.))/10.))/3.
         else         
            BES = (sin(x)/x - cos(x))/x
         end if
      case(2)
         if (x.lt.1.e-2) then            
            x2 = x*x
            BES = (x2*(1.0 - ((1.0 - x2/36.)*x2)/14.))/15.
         else
            x2 = 3d0/x/x
            BES = sin(x)/x*(x2-1d0) - cos(x)*x2
         end if
      case(3)
         x2 = x*x
         if (x.lt.1.e-2) then
            BES =(x**3*(1.-((1.-x2/44.)*x2)/18.))/105. 
         else            
            BES = sin(x)/x2*(15d0/x2-6d0) - cos(x)/x*(15d0/x2-1d0)
         end if
      case default   ! this is from bes.f to save one call
         arg =X - dacos(-1d0) * dble(L) /2d0 
         m1 = int(L/2d0)                         
         m2 = int(abs((L-1)/2d0))                
         if ((X .gt. 1d0) .and. (X .gt. dble(L))) then 
            c  = 1d0
            pp = 1d0
            X2 = 2d0 * X
            d  = dble(L*(L+1)) / X2
            qq = d
            if (m1 .ne. 0) then
               do i = 1, m1
                  i2 = 2*i
                  pp = - pp * dble(L-i2+1) * dble(L-i2+2)/ (dble(i2) * dble(i2-1)) * dble(L+i2) *  dble(L+i2-1) / X2**2   
                  c  = c + pp
               end do
            end if
            if (m2 .ne. 0) then
               do i = 1, m2
                  i2 = 2*i
                  qq = - qq * dble(L+i2) * dble(L+i2+1) / (dble(i2) * dble(i2+1)) * dble(L-i2) * dble(L-i2+1) / X2**2
                  d = d + qq
               end do
            end if   
            BES = (c*sin(arg)  + d*cos(arg))/X
         else                 ! small |X|
            ii = 0
            BES = 0
            a1 = 1/dble(2*L+1)
            if (l .ne. 0) then
               do i = 1, L
                  a1 = a1 * 2d0 / dble(i+L) * X
               end do
            end if   
            if (abs(a1) .lt. 10e-11 ) goto 20
 10         BES = BES + a1
            ii = ii + 1
            Lii = L+ii
            a1 = a1 * (-1d0) * X**2 * dble(Lii) / dble(ii) / dble(2*Lii+1) / dble(2*Lii)         
            if (abs(a1) .gt. 10e-11) goto 10
 20         continue
         end if         
      end select  
      return
      end subroutine sbessel
            
!*-------------------------------------------------

!     This subroutine calculates Q_L+A_L functions
!     used in poshevmat
       subroutine qlandal(gp,nqmi,nchi,Li,Qlp,Q0p,aq,ionsh,&
     &           nqgen ,nqim ,xp,wp,imax,Alp)
!      include 'par.f'
!      include 'par.pos'
      use kgrid 
      use Positronium 
      use input_data

      implicit real*8 (a-h,o-z)       
      implicit integer(i-n)
      integer, parameter:: maxl=20 !, kmax=281

      real*8:: xgz(Ps%igzm),wgz(Ps%igzm),pl(0:maxl,Ps%igzm)
      real*8:: xgp(Ps%igpm),wgp(Ps%igpm)
      integer:: igz,igp
      real*8:: zasym
      real*8:: Qlp(nqmax,2*nqmax*Ps%igpm),Q0p(nqmax,2*nqmax*Ps%igpm)&
     &           ,gp(Ps%iglp),Amg(Ps%iglp),Alp(nqmax,2*nqmax*Ps%igpm)
      real*8:: aq(0:nqmax),qq(nqmax)
      real*8:: xp(2*nqmax*Ps%igpm+nqmax),wp(2*nqmax*Ps%igpm)
      real*8:: xx(6),yy(6)
      real*8:: resA,resRd,resRd0, Rd, errQl

!       PAUSE'Entered qlandal'
!       if(Li.gt.lstoppos+latop)  RETURN
       if(nqmi.eq.1) return

       Rd = data_in%Rd ;  Rd2 = Rd / 2d0
       igp = Ps%igpm; xgp(1:igp) = Ps%xgp(1:igp); wgp(1:igp) = Ps%wgp(1:igp);
       igz = Ps%igzm; xgz(1:igz) = Ps%xgz(1:igz); wgz(1:igz) = Ps%wgz(1:igz);
       pl(0:20,1:igz)=Ps%pl(0:20,1:igz);
       iglp = Ps%iglp; gp(1:iglp) = Ps%gp(1:iglp); Amg(1:iglp) = Ps%Amg(1:iglp)

       knqm=nqmi
!!       nqm=-1 !  1  for on-shell!
       call qmesh(gridk(1,nchi),igp,xgp,wgp,knqm,imax,xp,wp,aq,ionsh,nqgen,nqim)
          imaxp=2*nqgen*igp 
          mnqmi=nqmi

      ddp=1.e+08
      pdp=0.0065d0
      nd1=2
      intp=2*nd1+1
       nnmin=INT(Log(1.E-04*ddp)/pdp) !1000
       gpmin=gp(nnmin)
       gpmax=400. !gp(iglp-5) 
       
       Q0p=0.d0; Qlp=0.d0;

       do iq=1,nqmi  ! 1- for on-shell and born      
          qa=dble(gridk(iq,nchi))
          if(qa.lt.0.) cycle 
           qa2=qa*qa

!       write(19,*)
!       write(19,'(a1,e20.10)')'#', Rd

        do i=1, imaxp 
          pp=dble(xp(i))
!          if(pp.gt.40.) cycle
          pp2=pp*pp

!      PRINT*, 'i=1,imaxp loop', i

           p2q2=pp2+qa2
           ppqa2=2.d0*pp*qa
           arg=p2q2/ppqa2

       call funleg(arg,Li,Q0p(iq,i),Qlp(iq,i))

            resA = 0.d0
            resRd = 0.d0
            resRd0 = 0.d0
!       IF(Z.eq.2d0) THEN
!         do iz=1,igz
!            qapp2=p2q2-ppqa2*xgz(iz)
!c             Akq1= (1-256./(16.+qapp2)**2.d0)/qapp2 ! for He only
!             Akq=(32.+qapp2)/(256.+32.*qapp2+qapp2*qapp2)
!             resA=resA+Akq*pl(Li,iz)
!          enddo
!        ENDIF
         
       IF(Rd.gt.0d0) THEN
         do iz=1, igz
            qapp2=p2q2-ppqa2*xgz(iz)
            qapp = sqrt(qapp2)
            qapp3 = qapp*qapp2
            Rd2qapp = Rd2*qapp
            Rd2qapp3 = Rd2*qapp3

        nn=INT(Log(qapp2*ddp)/pdp)
        if(nn.lt.3.or.nn.gt.iglp-2) cycle

        IF(qapp2.gt.gpmin) THEN
          yy(1:intp)=Amg(nn-nd1:nn+nd1)
          xx(1:intp)=gp(nn-nd1:nn+nd1)
          call intrpl(intp,xx,yy,1,qapp2,Akq) ! for both Mg and He
        ELSEIF(qapp2.le.gpmin) THEN 
           Akq=(32.+qapp2)/(256.+32.*qapp2+qapp2*qapp2) ! if p->0 (r->oo) it should converge to He
        ENDIF
 
           resA = resA+Akq*pl(Li,iz)
           resRd = resRd + dsin(Rd2qapp)/Rd2qapp3 * pl(Li,iz)
           resRd0 = resRd0 + 1d0/qapp2 * pl(Li,iz)
         enddo
      ELSEIF(Rd.eq.0d0) THEN
         do iz=1, igz
            qapp2=p2q2-ppqa2*xgz(iz)
            qapp = sqrt(qapp2)
            qapp3 = qapp*qapp2
            Rd2qapp = Rd2*qapp
            Rd2qapp3 = Rd2*qapp3

        nn=INT(Log(qapp2*ddp)/pdp)
        if(nn.lt.3.or.nn.gt.iglp-2) cycle

!        IF(qapp2.gt.gpmin) THEN
          yy(1:intp)=Amg(nn-nd1:nn+nd1)
          xx(1:intp)=gp(nn-nd1:nn+nd1)
          call intrpl(intp,xx,yy,1,qapp2,Akq) ! for both Mg and He
!        ELSEIF(qapp2.le.gpmin) THEN 
!           Akq=(32.+qapp2)/(256.+32.*qapp2+qapp2*qapp2) ! if p->0 (r->oo) it should converge to He
!        ENDIF
 
           resA = resA+Akq*pl(Li,iz)
           resRd = resRd + 1d0/qapp2 * pl(Li,iz)
           resRd0 = resRd
         enddo

       ENDIF
 
!        write(19,'(5e20.10)') qa,pp,Qlp(iq,i),resRd*qa*pp, qa*pp*resA 
              errQl=Qlp(iq,i)-resRd0*qa*pp
              Qlp(iq,i) = resRd*qa*pp + errQl
              Alp(iq,i) = qa*pp*resA   

       enddo 
      enddo 
!       STOP 'Ql are made'
        RETURN
    end subroutine qlandal
  

!c*****************************************************************
!c  this subroutine is to create composite mesh
!c  for q-integration in posvmat
!c*****************************************************************
      subroutine qmesh(gki,igp,xgp,wgp,nqmi,imax,xp,wp,aq,ionsh,nqgen,nqim)
      implicit real*8 (a-h,o-z)       
      implicit integer (i-n)
!      include 'par.f'
!      include 'par.pos'
!!!      integer, parameter::! kmax=281
      real*8:: xgp(igp),wgp(igp)
      real*8:: gki(nqmi)
      real*8:: aq(0:nqmi),qq(nqmi)
      real*8:: xp(2*nqmi*igp+nqmi),wp(2*nqmi*igp)

!      RETURN

      if(gki(1).gt.0.d0) then
!*     putting on-shell point to its place

         ionsh=1
         do i=2,nqmi
            if(gki(i).lt.gki(1)) then
               qq(i-1)=dble(gki(i))
               ionsh=i
            else
               qq(i)=dble(gki(i))
            endif
         end do
         qq(ionsh)=dble(gki(1))

         do ji=1,nqmi
            j0=2*nqmi*igp+ji
            xp(j0)=qq(ji)
         end do
         nqgen=nqmi
         nqim=nqmi
         imax=2*nqmi*igp+nqmi

      else
!*     on-shell point is negative (no on-shell point)
  
        ionsh=-1
        if (nqmi.eq.1) then
         nqim=1
         nqgen=nqim
         imax=2*nqim*igp+nqim
         else 
         do i=2,nqmi
           qq(i-1)=dble(gki(i))
           end do
              nqim=nqmi-1
              nqgen=nqim
              imax=2*nqim*igp+nqim
         do ji=1,nqim
            j0=2*nqim*igp+ji
            xp(j0)=qq(ji)
         end do
        endif

!*### the way as below gives wrong nqim in some cases.. why?
!         do i=2,nqmi
!            qq(i-1)=gki(i)
!         end do
!         nqim=nqmi-1
!         nqgen=nqim
!         imax=2*nqim*igp
      endif

!*     end points
      aq(0)=0.d0
      do i=1,nqgen-1
         aq(i)=(qq(i+1)+qq(i))*0.5d0
!*     forming composite mesh
         do ji=1,igp
            j1=2*(i-1)*igp+ji
            xp(j1)=aq(i-1)+ (qq(i)-aq(i-1))*xgp(ji)
            wp(j1)=(qq(i)-aq(i-1))*wgp(ji)
            j2=(2*i-1)*igp+ji
            xp(j2)=aq(i)-(aq(i)-qq(i))*xgp(igp+1-ji)
            wp(j2)=(aq(i)-qq(i))*wgp(igp+1-ji)
         end do
      end do
!*     last 2 intervals
      do ji=1,igp
         j1=2*(nqgen-1)*igp+ji
         xp(j1)=aq(nqgen-1)+(qq(nqgen)-aq(nqgen-1))*xgp(ji)
         wp(j1)=(qq(i)-aq(i-1))*wgp(ji)
!*     last interval
         j2=(2*nqgen-1)*igp+ji
         xp(j2)=qq(nqgen)/xgp(igp+1-ji)
         wp(j2)=qq(nqgen)*wgp(igp+1-ji)/xgp(igp+1-ji)/xgp(igp+1-ji)
      end do
!*     composite mesh is ready 
      return
      end subroutine qmesh

!*----------------------------------------------

      subroutine funleg(z, Lla, Q0, result)      
!C     returns the Legendre polinomial of the second kind - Ql
!C     for given Lla and z with the use of the recurrence
!C     realtion
            
!      include 'par.f'
!      include 'par.pos'
      implicit real*8 (a-h,o-z)
      implicit integer(i-n)
      integer, parameter:: ltmax=20, maxl=2*20+1
      real *8, dimension (0:Lla) :: QN

!*     arg of Legendre function is always .ge. 1
!*     min[z]=1 at pp=qa
!*     Mathematica cannot handle this region.
            
      Q0 = 0.5d0 * log((z+1d0)/(z-1d0))
      Q0sp = real(Q0)      
      
      zc = 0d0
      select case (Lla)
      case (0:4)
         zc = 4.0d0
      case (5)
         zc = 3.0d0
      case (6:7)
         zc = 2.5d0
      case (8:)
         zc = 2.0d0
      end select
      
!C     use analitical expression for small z
      if (z.lt.zc) then      
         select case(Lla)
         case(0)
            result = Q0
         case(1)
            result = z * Q0 - 1d0
         case(2)
            z3 = 3d0 * z
            result = ((z3 * z - 1d0) * Q0- z3)/2d0
         case(3)
            z2 = z*z
            result = ((5d0*z2 - 3d0)*z*Q0 - (5d0 * z2-4d0/3d0))/2d0
         case(4)
            z2 = z*z
            result = (((35d0*z2 - 30d0)*z2 + 3d0)*Q0 &
     &           - (35d0*z2 - 55d0/3d0)*z)/8d0
         case(5)
            z2 = z * z
            z4 = z2 * z2
            result = -0.5333333333333333 + 6.125*z2 - 7.875*z4 + 0.9375 &
     &           *z*(1. - 4.666666666666666 * z2 + 4.2 * z4) * 2d0 * Q0
         case(6)
            z2 = z * z
            z4 = z2 * z2
            result= (-2.8875 + 14.875*z2 - 14.4375*z4)*z - 0.15625  &
     &           *(1d0-21d0*z2 + 62.99999999999999*z4-46.2*z2*z4)*2d0 &
     &           *Q0
         case(7)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            result=0.4571428571428571-10.6125*z2+34.375*z4 - 26&
                &.81249999999999*z6 - 1.09375*z*(1.-9.*z2+19.8*z4-12&
                &.25714285714285*z6)*2.d0*Q0
         case(8)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            z8 = z4 * z4
            result=(3.383705357142857 - 32.9140625*z2 + 77.0859375*z4 - &
     &           50.2734375*z6)*z + 0.13671875* (1. - 36.*z2 + 198.*z4 - &
     &           343.1999999999999*z6 + 183.8571428571428*z8)*2d0*Q0
         case(9)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            z8 = z4 * z4
            result=-0.4063492063492063 + 15.82477678571428*z2 - 92&

                &.7265625*z4 + 169.4401041666666*z6 - 94.9609375*z8 + 1&
                &.23046875*z* (1. - 14.66666666666666*z2 + 57.2*z4 - 81&
                &.7142857142857*z6 + 38.58730158730158*z8)*2d0*Q0
         case(10)
            z2 = z * z
            z4 = z2 * z2
            z6 = z2 * z4
            z8 = z4 * z4
            result=(-3.817398313492063 + 59.68973214285714*z2 - 245&
                &.5578125*z4 + 367.1822916666666*z6 - 180.4257812499999&
     &           *z8)*z - 0.123046875*(1. - 55.*z2 + 476.6666666666667 &
     &           *z4 - 1430.*z6 + 1736.428571428571*z8 - 733.15873015873 &
     &           *z8*z2)*2.d0*Q0

         case default
               result = 0.d0  !QN(Lla)
         end select        
      else           
         result = FlegQ(Lla, z) !  series expansion for large z                        
      end if                    ! z

      return
      end subroutine funleg     


      real*8 function FlegQ(l, z)
!c------------------------------------------------------     
!c     returns the Legendre function of the second kind 
!c            
!c            f 1
!c         1  !      P(l,x) 
!c     Q = -  | dx ---------- 
!c         2  |     (z - x)
!c            J -1     
!c             
!c     where P(n,x) is the Legendre polinomial 
!c     and the argument z is large
!c-------------------------------------------------------
!      use apar                  ! maxx, nmax, lmax1, maxq 
!     implicit real*8 (a-h,o-z)     
!      include 'par.f'
!      include 'par.for'
      real*8, intent(in) :: z      
      real *8 Q, zn
      parameter (jmaxC = 5, lmaxC = 20)
      real*8, dimension (0:lmaxC, 0:jmaxC) :: c
      data c/1.,0.3333333333333333,0.1333333333333333,0&
          &.05714285714285714,0.0253968253968254,0.01154401154401154,0&
          &.005328005328005328,0.002486402486402486,0.001170071758307052&
          &,0.0005542445170928143,0.0002639259605203878,0&
          &.0001262254593793159,0.00006058822050207163,0&
          &.00002917210616766412,0.00001408308573611371,6&
          &.814396323925989d-6,3.303949732812601d-6,1.604775584508978d-6&
          &,7.807016357070702d-7,3.803418225239573d-7,1.855325963531499d&
          &-7,0.3333333333333333,0.2,0.1142857142857143,0&
          &.06349206349206349,0.03463203463203463,0.01864801864801865,0&
          &.009945609945609946,0.005265322912381736,0.002771222585464072&
          &,0.001451592782862133,0.0007573527562758953,0&
          &.0003938234332634656,0.0002042047431736488,0&
          &.0001056231430208528,0.00005451517059140791,0&
          &.00002808357272890711,0.0000144429802605808,7&
          &.416665539217167d-6,3.803418225239573d-6,1.948092261708074d-6&
          &,9.966983664552936d-7,0.2,0.1428571428571429,0&
          &.09523809523809524,0.06060606060606061,0.0372960372960373,0&
          &.02237762237762238,0.01316330728095434,0.007620862110026197,0&
          &.004354778348586398,0.00246139645789666,0.001378382016422129&
          &,0.0007657677869011831,0.0004224925720834113,0&
          &.0002316894750134836,0.000126376077280082,0&
          &.00006860415623775879,0.00003708332769608583,0&
          &.00001996794568250776,0.00001071450743939441,5&
          &.731015607117938d-6,3.056541657129567d-6,0.1428571428571429,0&
          &.1111111111111111,0.08080808080808081,0.05594405594405594,0&
          &.0372960372960373,0.02413273001508296,0.01524172422005239,0&
          &.009435353088603863,0.005743258401758873,0.003445955041055324&
          &,0.002042047431736488,0.001197062287569665,0&
          &.0006950684250404509,0.0004001909113869263,0&
          &.0002286805207925293,0.0001297916469363004,0&
          &.00007321580083586177,0.00004107227851767856,0&
          &.00002292406242847175,0.00001273559023803986,7&
          &.045220131681626d-6,0.1111111111111111,0.09090909090909091,0&
          &.06993006993006993,0.05128205128205128,0.03619909502262443,0&
          &.02476780185758514,0.01651186790505676,0.01076860950329789,0&
          &.006891910082110647,0.004339350792440037,0.002693390147031747&
          &,0.001650787509471071,0.001000477278467316,0&
          &.0006002863670803894,0.0003569270290748261,0&
          &.0002104954274031026,0.0001232168355530357,0&
          &.00007163769508897422,0.00004139066827362955,0&
          &.00002377761794442549,0.00001358721025395742,0&
          &.09090909090909091,0.07692307692307692,0.06153846153846154,0&
          &.04705882352941176,0.0346749226006192,0.02476780185758514,0&
          &.01722977520527662,0.0117162471395881,0.007810831426392067,0&
          &.00511744127936032,0.003301575018942142,0.002101002284781363&
          &,0.001320630007576857,0.0008209321668721001,0&
          &.0005051890257674462,0.0003080420888825892,0&
          &.000186258007231333,0.0001117548043387998,0&
          &.00006657733024439136,0.00003940290973647652,0&
          &.00002317818219792737/

      Q = 0d0          
            
!c     series expansion              
!c         zn = 1/z**(l+1)
         do j = jmaxC, 0, -1                                
            n = l+2*j+1
            zn = 1.0d0/z**n
            Q = Q + c(l,j) * zn
         end do         
!c$$$    end if                    ! z
      FlegQ = Q      
      return
      end function FlegQ

 
