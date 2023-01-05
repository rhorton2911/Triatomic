C  The following routine gets the regular solution of the following
C  differential equation. ECMN is the energy in Rydbergs.
C        2
C       d S   [ LN(LN+1)       V(X)          ]
C       --- - [ --------  + --------- - ECMN ] S =0
C         2   [   X*X           X            ]
C       dX      
C
C  Note that UCENTR(X) = V(X)/X and CNTFUG(X) = LN(LN+1)/X/X
      
      subroutine regular(ln,ecmn,eta,ucentr,cntfug,ldw,gridx,nx,
     >   jdouble,njdouble,regcut,expcut,reg,jstart,jstop,pshift,DELTA)

      use MPI_module
       
      implicit real*8 (a-h,o-z)
      
      real*8:: pshift,DELTA
      dimension jdouble(njdouble),ucentr(nx),reg(nx),gridx(nx),
     >     cntfug(nx),f(nx),g(nx)
      complex phase,sigc
      complex*16 xx,eta1,zlmin,cfc(9),cgc(9),cfcp(9),cgcp(9),sig(9),ci,
     >     coulphase
! Mark: Added for COUL90
      integer:: ANG_L
      real*8, dimension(0:100):: FCR, GCR, FCPR, GCPR
      data ci/(0d0,1d0)/
      
      logical match
      character file*10
      
      match = .false.
      
      pshift = 0d0
      DELTA = 0d0
      phase = (0,0)
      sigc = (0,0)
      
      do j = 1, nx
         reg(j) = 0.0d0
      enddo 

      if( ecmn .lt. 0 ) then
         jstart = 1
         jstop = 1
         return
      endif
      wnn=sqrt(abs(ecmn))          
      
      if (ecmn.le.0.0.or.ln.gt.3.or.ucentr(1).ne.0.0) then
         if ( eta <= 0d0) then
            jstart = jfirst1(ln,ecmn,gridx,nx,jdouble,njdouble,regcut)
         else 
            jstart = jfirst_eta(ln,eta,gridx,nx,njdouble,jdouble,regcut) 
         end if
c     print*,'jfirst1: jstart=', jstart , ln, 
c     >        sqrt(ecmn)*gridx(jstart)
 20      if (jstart.gt.10) then
            
! Mark: Solutions to Coulomb waves should be for eta /= 0. Mark added coulrepw for checking.
	    if (eta == 0d0) then
               s1=appf1(ln,wnn,gridx(jstart-1),acc)
            else if (eta > 0d0 ) then
               s1=coulrepw(ln,eta,wnn*gridx(jstart-1),y_sum,acc)
            else
               s1 = coulrho(ln,eta,wnn*gridx(jstart-1),acc)
            end if
            if (abs(s1).gt.regcut.or.acc.lt.1e-6) then
               jstart = jstart - 10
               go to 20
            endif 
C     make sure that we are not near a point where DX doubles
            do jd = 2, njdouble - 1
               if (abs(jstart-jdouble(jd)).lt.3) jstart=jdouble(jd)+3
            end do       
         endif 
         jstop  = jlast1(ecmn,gridx,nx,jdouble,njdouble,expcut)
c$$$  if (ecmn.gt.200.0) then
c$$$  jstop = nx
c$$$  jstart = nx + 1
c$$$  endif 
         
c$$$  print*,'jstart=', jstart

         jstartorig = jstart
         if (ucentr(1) .eq. 0.0) then
c$$$  if (ucentr(1) .eq. 0.0 .and. wnn .gt. 2.0) then
            rho = wnn * gridx(jstart)
            acc = 1.0
c$$$  C  Need to take care that jstart stays less than jstop. For this reason
c$$$  C  the wnn > 2 condition was added above.
            do while (rho.lt.ln.and.acc.gt.1e-6.and.jstart.lt.jstop)
               reg(jstart) = appf1(ln,wnn,gridx(jstart),acc)
               jstart = jstart + 1
               rho = wnn * gridx(jstart)
            enddo
C     make sure that we are not near a point where DX doubles
            do jd = 2, njdouble - 1
               if (abs(jstart-jdouble(jd)).lt.3) jstart=jdouble(jd)-3
            end do
            jstart = max(jstart,1)
         endif 
         
         if (jstart.gt.jstop) return
         
         if (jstart.eq.1) then
            s1 = 0.0d0
         else
            s1=appf1(ln,wnn,gridx(jstart-1),acc)            
            if (eta > 0d0 ) then 
               s1 = coulrepw(ln,eta,wnn*gridx(jstart-1),y_sum,acc)	
            else if (eta < 0d0 ) then
               s1 = coulrho(ln,eta,wnn*gridx(jstart-1),acc)
            end if
               
         end if 

         s2=appf1(ln,wnn,gridx(jstart),acc)
         if (eta > 0d0) then 
            s2 = coulrepw(ln,eta,wnn*gridx(jstart),y_sum,acc)
         else if ( eta < 0d0 ) then
            s2 = coulrho(ln,eta,wnn*gridx(jstart),acc)
         end if



c$$$  if (eta.eq.0.0.and.abs((s2-s2old)/(s2old+1e-10)).gt.1e-4)then
c$$$  print*,'S2,S2OLD,rho,l,jstart',s2,s2old,wnn*gridx(jstart),
c$$$  >         jstart
c$$$  endif 


      
         call numerovf(ln,ecmn,ucentr,cntfug,gridx,nx,
     >        jdouble,njdouble,s1,s2,reg,jstart,jstop)


         jstart = jstartorig
         jmatch = min(nx,jstop)
         rho = wnn * gridx(jmatch)
         tmp = (ecmn - ucentr(jmatch)) / ecmn - 1.0
         
         
         
         if (eta.eq.0.0) then
            
            do while (abs(tmp).lt.sqrt(regcut).and.rho.gt.1.0)
c$$$  do while (abs(tmp).lt.sqrt(regcut).and.jmatch.gt.jstart+40)
               jmatch = jmatch - 1
               rho = wnn * gridx(jmatch)
               tmp = (ecmn-ucentr(jmatch)) /
     >              (ecmn) - 1.0



            enddo 
c$$$            if (jmatch.eq.min(nx,jstop).and.jmatch.gt.jstart+40.and.
c$$$     >         ecmn.gt.1.0.and.abs(ucentr(jmatch)/
c$$$     >         (cntfug(jmatch)+1e-30)).gt.2e-2) then
c$$$               print*,'WARNING: abs(ucentr(jmatch)/cntfug(jmatch)) > ',
c$$$     >            '2e-2',ucentr(jmatch)/(cntfug(jmatch)+1e-30)
c$$$            endif 
         else
            asympot = 2.0d0*eta*ecmn/wnn/gridx(jmatch)
            tmp=abs(ucentr(jmatch)/asympot-1.0d0)



            do while (tmp.lt.1e-5.and.rho.gt.0.1.and.jmatch.gt.jstart)
c$$$            do while (tmp.lt.1e-5.and.jmatch.gt.jstart+40)
               jmatch = jmatch - 1
               rho = wnn * gridx(jmatch)
               asympot = 2.0d0*eta*ecmn/wnn/gridx(jmatch)
               tmp = abs(ucentr(jmatch)/asympot-1.0d0)
            enddo 

c$$   
      

c$$$ Mark: Changed for repulsive waves
            if (jmatch.eq.min(nx,jstop) .AND. eta <= 0d0) then
!            if (jmatch.eq.min(nx,jstop)) then
               if (myid == 0) then
                  print*,'WARNING: asymptotic potential must be ',
     >                 'Coulomb ECMN,jmatch,tmp:',ecmn,jmatch,tmp
               end if
            end if

         end if
         

         if (match.or.ucentr(1).ne.0.0.or.ecmn.lt.0.0) then
c$$$            print*,'JMATCH, RMATCH,U',jmatch,gridx(jmatch),
c$$$     >         ucentr(jmatch)
            if (abs(ecmn-34.527).lt.-1e-2) then 
               do jm = jmatch-50, jmatch+400
                  do j = jstart, jstop
                     f(j) = reg(j)
                  enddo 
                  jmat = jm
                  call matchcc(ln,eta,ecmn,jstart,jstop,gridx,jmat,f,
     >                 phase,sigc,j1,j2,pshift,DELTA,ldw)
                  if ( myid == 0) print*,j1,j2,jmatch,phase,ucentr(jmat)
               enddo
            else
               call matchcc(ln,eta,ecmn,jstart,jstop,gridx,jmatch,reg,
     >            phase,sigc,j1,j2,pshift,DELTA,ldw)

c$$$         print*,sqrt(ecmn),jmatch,phase,reg(jstart),reg(jstart+1)
               diff = abs(sigc-coulphase(dble(eta),ln))
               if (diff.gt.1e-5) then
                  if ( myid == 0) print*,'Coulphase difference:',sigc,
     >               coulphase(dble(eta),ln)
c$$$                  stop 'Coulphase difference > 1e-5'
               endif 
               if (ln.gt.ldw) then
                  diff = abs(real(phase)-1.0)
                  if (diff.gt.1e-5) then
                     if ( myid == 0) then
                        print*,'Expected phase to be 1 for L, ECMN:',
     >                       phase,ln,ecmn
                        write(file,'(i1,"_",1p,e8.2)') ln,ecmn
                        mode1 = 12
                        kfn = 0
                        zlmin = cmplx(ln)
                        nl = 1
                        eta1 = eta
                        open(57,file=file)
                        do j = jstart, jstop
                           xx=sqrt(cmplx(dble(ecmn)))*dble(gridx(j))

! Mark: Using subroutine COUL90
!                        call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,
!     >                     sig,mode1,kfn,ifai)
                        
                           ifai = 0
                           rho = REAL(xx)  
                           Lang = dble(ln)
                           ANG_L = INT(Lang)
                           call coul90(rho, eta, dble(Lang), 0, FCR,GCR, 
     >                          FCPR,GCPR,kfn,ifai)
                           cfc(1) = CMPLX(FCR(ANG_L))
                           cgc(1) = CMPLX(GCR(ANG_L))+
     >                          ci*CMPLX(FCR(ANG_L))
                           
                           write(57,'(5e14.6,i8)') gridx(j), reg(j),
     >                          real(cfc(1)),
     >                          coulrho(ln,eta,wnn*gridx(j),acc),acc,j
c$$$  reg(j) = real(cfc(1))
                        enddo
                        close(57)
                     end if !myid
                  endif 
               endif 
            endif                
         
! Mark: Need to check, this should not be here? Test distorted wave for He
            if (eta.eq.0.0.and.ln.le.3.and.ucentr(1).ne.0.0) then
               call ricbessel(wnn,ln,gridx,jstop,jstart,f)
               call ricbessel(wnn,-ln-1,gridx,jstop,jstart,g)
               rph = real(phase)
               aph = aimag(phase)
               do j = jmatch, jstop
                  reg(j) = rph * f(j) + aph * g(j)
               enddo 
            endif 
         else
            phase = (1.0d0,0.0d0)
            sigc = (1.0d0,0.0d0)
         endif 
            

         do while (abs(reg(jstart)).lt.regcut)
            jstart = jstart + 1
            if (jstart.gt.jstop) return
         end do
         if (match.and.ln.gt.ldw) then
c$$$         if (match.and.ln.gt.ldw.and.eta.ne.0.0) then
            if (abs(aimag(phase)).gt.1e-5) then
               if ( myid == 0) then
                  print*,'Old phase (NOT reset to 1) for L, ecmn:',
     >                 phase,ln,ecmn
               end if
            end if
c$$$            phase = (1.0,0.0)
            mode1 = 12
            kfn = 0
            zlmin = cmplx(ln)
            nl = 1
            eta1 = eta

c$$$            write(file,'(i1,"_",1p,e8.2)') ln,ecmn
c$$$            open(57,file=file)
            acc = 1.0
            nerrs = 0
            do j = jstart, jstop
               if (acc.gt.1e-6) then
                  cfc(1) = cmplx(coulrho(ln,eta,wnn*gridx(j),acc))
               else 
                  xx=sqrt(cmplx(dble(ecmn)))*dble(gridx(j))
! Mark: Using subroutine COUL90
!                        call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,
!     >                     sig,mode1,kfn,ifai)
                        
                  rho = REAL(xx)  
                  L = dble(ln)
                  ANG_L = INT(L)
                  call coul90(rho, eta, dble(L), 0, FCR,GCR,FCPR,GCPR, 
     >                 kfn,ifai)
                  cfc(1) = CMPLX(FCR(ANG_L))
                  cgc(1) = CMPLX(GCR(ANG_L))+ci*CMPLX(FCR(ANG_L))
               endif 
c$$$               if (wnn*gridx(j).lt.1.0) then
c$$$                  tmp = appf1(ln,wnn,gridx(j),acc)
c$$$               else
c$$$                  tmp = 0.0
c$$$               endif 
c$$$               write(57,'(10e14.6)') gridx(j), reg(j),
c$$$     >            real(cfc(1)), tmp
               if (abs((reg(j)+real(cfc(1)))/
     >            (reg(j)-real(cfc(1))+1e-20)).lt.1e-3) then
                  if ( myid == 0) then
                     print*,ln,ecmn,gridx(j), phase, 
     >                    reg(j), real(cfc(1))
c$$$                  cfc(1) = abs(real(cfc(1))) * reg(j)/abs(reg(j))
                     print*,'cfc of wrong sign?:',reg(j-2),reg(j-1),
     >                    reg(j),reg(j+1),reg(j+2)
                  end if !myid
               endif
               if (abs((reg(j)-real(cfc(1)))/(reg(j)+1e-20)).gt.0.1.
     >            and.abs(reg(j)).gt.0.1)
     >            nerrs = nerrs + 1
c$$$     >            ,print*,'>10% error in reg(j):',ln,ecmn,j,reg(j),
c$$$     >            real(cfc(1))
               reg(j) = real(cfc(1))
            enddo
            if (nerrs.gt.0) print*,'WARNING; L, ECMN, NERRS:',
     >         ln,ecmn,nerrs
c$$$  close(57)
         endif 

c$$$      if (abs(reg(jstart)/regcut).gt.1e+4.and.jstart.gt.1) then
c$$$         print*,'JSTART should be smaller. JSTART, REG(JSTART), LN,'
c$$$     >      //' ETA, ECMN, JSTOP',jstart,reg(jstart),ln,eta,ecmn,jstop
c$$$      end if 
      else
         phase = (1.0d0,0.0d0)
         sigc  = (1.0d0,0.0d0)
         jstart = 1
         jstop = nx
         call ricbessel(wnn,ln,gridx,jstop,jstart,reg)
      end if

 
      return
      end
!!$------------------------------------------------------------------------------      
      subroutine ricbessel(wnn,ln,gridx,jstop,jstart,reg)
      use MPI_module
      implicit real*8 (a-h,o-z)

      real*8 gridx(jstop), reg(jstop)
      real*8 rk, rhosmall
      rhosmall = 1d-3
      if (ln.eq.0) then
         jstart = 1
         do i = jstart, jstop
            rk = wnn * gridx(i)
            reg(i) = sin(rk)
         enddo
      else if (ln.eq.1) then
         jstart = 1
         do i = jstart, jstop
            rk = wnn * gridx(i)
            if (rk.lt.rhosmall) then
               reg(i) = rk * rk / 3d0 - rk * rk * rk * rk / 30d0
            else 
               reg(i) = sin(rk) / rk - cos(rk)
            endif 
         enddo
      else if (ln.eq.2) then
         jstart = 1
         do i = jstart, jstop
            rk = wnn * gridx(i)
            if (rk.lt.rhosmall) then
               reg(i) = rk * rk * rk / 15d0
            else 
               reg(i) = sin(rk) * (3d0/rk/rk - 1d0) - 3d0 * cos(rk)/rk
            endif 
         enddo
      else if (ln.eq.3) then
         jstart = 1
         do i = jstart, jstop
            rk = wnn * gridx(i)
            if (rk.lt.rhosmall) then
               reg(i) = rk * rk * rk * rk / 105d0
            else 
               reg(i) = ((15d0/rk/rk/rk-6d0/rk) * sin(rk) -
     >            (15d0/rk/rk - 1d0) * cos(rk))
            endif
c$$$            print*,appf1(ln,wnn,gridx(i),acc),reg(i)
         enddo 
      else if (ln.eq.-1) then
         do i = jstart, jstop
            rk = wnn * gridx(i)
            reg(i) = cos(rk)
         enddo
      elseif (ln.eq.-2) then
          do i = jstart, jstop
            rk = wnn * gridx(i)
            reg(i) = sin(rk) + cos(rk) / rk
         enddo
      elseif (ln.eq.-3) then
          do i = jstart, jstop
            rk = wnn * gridx(i)
            reg(i) = 3d0 * sin(rk) / rk + (3d0/rk/rk - 1d0) * cos(rk)
         enddo
      elseif (ln.eq.-4) then
          do i = jstart, jstop
            rk = wnn * gridx(i)
            reg(i) = - sin(rk) * (1d0 - 15d0/rk/rk)
     >         - cos(rk) * (6d0/rk - 15d0/rk/rk/rk)
         enddo
      else 
          if ( myid == 0) 
     >        print*,'Riccati-Bessel functions are not",
     >        " coded for |L| > 3.'
      end if
c$$$      real rho
c$$$      complex*16 xx,eta1,zlmin,cfc(1),cgc(1),cfcp(1),cgcp(1),sig(1)
c$$$
c$$$      xx = rho
c$$$      eta1 = (0.0,0.0)
c$$$      zlmin = dcmplx(l)
c$$$      nl = 1
c$$$      ifai = 0
c$$$      mode1 = 4
c$$$      kfn = 0
c$$$      call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,kfn,ifai)
c$$$      ricbessel = cfc(1)
      return
      end
!!$------------------------------------------------------------------------------      

      subroutine numerovf(ln,en,ucentr,cntfug,gridx,nx,
     >   jdouble,njdouble,rs1,rs2,reg,jstart,jstop)
      implicit real*8 (a-h,o-z)

      real*8 ucentr(nx),reg(nx),gridx(nx),cntfug(nx),en,rs1,rs2,
     >   appf1, acc
      dimension jdouble(njdouble)

      if (jstart.ge.jstop) return
      ecmn = en
      x = gridx(jstart)
      dx= gridx(min(jstart+1,nx)) - x
      f1 = 0.0d0
      f2 = ucentr(jstart)+cntfug(jstart)-ecmn
      h1=dx*dx
      h2=h1/12d0
      s1 = rs1
      s2 = rs2
      

      if (jstart.eq.1) then
C  We get here if the integration is going to start from the first X point.
C  This means that S1 is the solution of the differential equation at X=0
C  S2 is the solution at the first X point. 
         s1=0d0
         t1=0d0
C  LN = 1 is a special case
         if (ln.eq.1) t1=-h1/18d0
         if (ecmn.ne.0d0) t1=t1*ecmn
      else
         j=jstart-1
         f1=ucentr(j)+cntfug(j)-ecmn
         t1=s1*(1d0-h2*f1)
      end if

      reg(jstart) = s2
      t2=s2*(1d0-h2*f2)
      
      istart=2
      do while (jstart.gt.jdouble(istart).and.istart.lt.njdouble)
         istart=istart+1
      end do
      istart=istart-1
      istop=njdouble-1
C  JDOUBLE(ISTART) points to the first doubling of DX that happens after JSTART
C  JDOUBLE(ISTOP) points to the last doubling of DX that happens before JSTOP

C    integration loop
      f0=0d0; f3=0d0; s0=0d0; s3=0d0
      do i=istart,istop
         j1=max(jstart,jdouble(i))+1
         j2=min(jstop,jdouble(i+1))
         do j=j1,j2
            f3 = ucentr(j)+cntfug(j)-ecmn
            t3 = 2.0d0*t2-t1+h1*f2*s2
            s3 = t3/(1d0-h2*f3)
c$$$            if (j.lt.jstart+10) then
c$$$               test=appf1(ln,sqrt(abs(en)),gridx(j),acc)
c$$$               print*,test/s3
c$$$               s3 = test
c$$$            endif 
            reg(j)=s3      

            t1=t2
            t2=t3
            f0=f1
            f1=f2
            f2=f3
            s0=s1
            s1=s2
            s2=s3
         end do
         dx=2d0*dx
         h1=dx*dx
         h2=h1/12d0
         t2=s3*(1d0-h2*f3)
         t1=s0*(1d0-h2*f0)
      end do
      return
      end
!!$------------------------------------------------------------------------------      




      function jfirst_eta(l, eta, gridx, nx, njdouble,jdouble, regcut)  
      use MPI_module
      implicit none
      real*8, intent(in)::  gridx(nx)
      integer, intent(in)::  jdouble(njdouble)
      real*8,intent(in):: eta
      real*8, intent(inout):: regcut
      integer,intent(in):: l, nx, njdouble
      
      real*8::  rho, acc
      real*8:: s1, y_sum,  coulrepw
      integer:: jfirst_eta, j, jd


      s1 = 0d0
 10   j = 1    
      do while ( j < nx )
         rho = gridx(j) / eta
         s1 = coulrepw(l,eta,rho,y_sum,acc)
      
         if ( ABS(s1)  > regcut ) then
            jfirst_eta = j
            exit
         end if
         
         j = j + 1
         
      end do

      if (acc < 1e-10) then
         regcut = regcut / 10.0d0
         if ( myid == 0) then 
            print*,"regcut, acc,j:",regcut,acc,j
            print '("regcut redefined. S1, REGCUT:",1p,e13.4,e8.1)', 
     >           s1, regcut
         end if
         go to 10
      end if

! Mark: Remove after testing
!      print*, j, jdouble(2)

      do jd = 2, njdouble - 1
         if (abs(j-jdouble(jd)) < 3) j = jdouble(jd) + 3
      end do 

      jfirst_eta = j
      
      end function jfirst_eta

! Mark: This function is a combination of yours and Igors function      
      function coulrepw(l, eta, rho, y_sum, acc)

      use MPI_module
      implicit none
      real*8, intent(in) :: eta, rho
      integer, intent(in):: l
      real*8, intent(out):: y_sum
      real*8, intent(out):: acc
      
      logical:: converged
      integer:: lp, n, nmax = 1000, lmax = 40
      real*8:: rhopow, term, sum, pi, delta = 1d-10
      real*8:: coulrepw, termn, termp
      real*8, allocatable, dimension(:,:):: A
      real*8, allocatable, dimension(:):: C

      pi = ACOS(-1d0)

      allocate( A(0:nmax,0:lmax), C(0:lmax) )
      if ( l > lmax ) then 
         if ( myid == 0) 
     >        print*, 'numerov.f: coulrho(): l > lmax:',l, lmax
         stop
      end if

! If we ever need to calculate very large eta (>120) over large rmax. 
! Comment statemensts out  below
      if( rho < 1d0) then
         if( nint(dble(l)*abs(log10(rho)) ) > 30) then
            coulrepw = 0d0
            acc = 1d-30
            return
         endif
      endif
    
! Making everything in terms of LN functions. 
! When eta > 120.0 e^(2*pi*eta) can't be calculated. 
      if ( eta > 100.0 ) then !  exp(2*pi*eta) - 1 ~  exp(2*pi*eta)
         C(0) = DLOG(2d0*pi*eta) - 2d0 * pi * eta
      else
         C(0) = DLOG(2d0*pi*eta / (exp(2d0*pi*eta)-1d0))
      end if
      C(0) = C(0) / 2d0

      rhopow = rho
      do lp = 1, l
         rhopow = rhopow * rho 
         C(lp) = sqrt(lp * lp + eta * eta)/ (lp*(2d0*lp+1d0))
         C(lp) = C(lp-1) + DLOG(C(lp))
      end do
      A(l+1,l) = 1d0
      A(l+2,l) = eta / (l+1d0)

      n = l + 1
      sum = 0d0   
      termp = 0d0
      termn = 0d0
      coulrepw = 0d0
      converged = .FALSE.
      do while (.NOT. converged .AND. n < nmax )
     
         term =  A(n,l) * rhopow
         if (term > 0d0) then
            termp = termp + term
         else 
            termn = termn + term
         end if
         sum = sum + term

         n = n + 1
         rhopow = rhopow * rho
         converged = ABS(term/sum) < delta 
 
         A(n+1,l) =  2d0 * eta * A(n,l) - A(n-1,l) 
         A(n+1,l) = A(n+1,l) / ((n-l)*(n+l+1))
      end do

! Used for checking power series expansion
      y_sum = sum
      if (n == nmax .OR. termp == -termn) then
         acc = 1d-30
      else 
         acc = abs(termn/termp+1d0)
      endif 

   
! Need to square F in order to have F>0 hence to take ln(F). Then divide both sides by 2.
      coulrepw = DLOG(ABS(sum)) + C(l) 
      coulrepw = exp(coulrepw)
! Need to keep correct sign of sum so we have correct F=+/-|F|.
      coulrepw = sign(coulrepw,sum) 
  

      end function coulrepw




      function coulrho(l,seta,srho,acc)
      use MPI_module
      implicit real*8 (a-h,o-z)

      parameter (nmax=1000, lmax = 40)

      real*8 seta, srho, coulrho, acc
      dimension C(0:lmax),A(nmax,0:lmax)
      logical converged
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
      data delta /1d-10/
     
      
      if( l .gt. lmax) then
         if (myid == 0) print*,'numerov.f: coulrho(): l > lmax:',l, lmax
         stop
      endif


      coulrho = 0.0d0
      eta = seta
      rho = srho
      if (eta.eq.0d0) then
         C(0) = 1d0
      else 
         C(0) = sqrt(2d0*pi*eta/(exp(2d0*pi*eta)-1d0))
      endif
      do lp = 1, l
         C(lp) = C(lp-1)*sqrt(lp*lp+eta*eta)/lp/(2*lp+1)
      enddo 
      converged = .false.
      A(l+1,l) = 1d0
      A(l+2,l) = eta/(l+1)
      n = l+1
      sum = 0d0
      termp = 0.0d0
      termn = 0.0d0
      
      if( rho .lt. 1d0) then
         if( nint(dble(l)*abs(log10(rho)) ) .gt. 30) then
            coulrho = 0d0
            acc = 1d-30
            return
         endif
      endif

      rhop = rho**l
      do while (.not.converged.and.n.lt.nmax)
         rhop = rhop * rho
         term = A(n,l) * rhop * C(l)
         if (term.ge.0.0) then
            termp = termp + term
         else
            termn = termn + term
         endif 
         sum = sum + term
         converged = abs(term/sum).lt.delta.and.abs(A(n,l)).gt.0d0
         n = n + 1
         A(n+1,l) = (2d0*eta*A(n,l)-A(n-1,l))/(n-l)/(n+1+l)
         if (termp.gt.1d30) then
            n = nmax
c$$$            print*,'Exiting before convergence'
         endif 
      enddo
      if (n.eq.nmax.or.termp.eq.-termn) then
         acc = 1d-30
      else
         acc = abs(termn/termp+1d0)
      endif 
      coulrho = sum 
      if (abs(sum) .le. 1d-30) then
         if ( myid == 0) print*, 'sum=',sum, ', l=', l, ', rho=',srho
         stop 'COULRHO <= 1d-30'
      endif 
      return
      end
!!$------------------------------------------------------------------------------      
C  Riccati-Bessel function for small rho. Same as spherical Bessel * rho
      function appf1(ln,w,x,acc)
      use MPI_module
      implicit real*8 (a-h,o-z)

      real*8 appf1, w, x, acc
      if (w.ne.0.0) then
C  We use the expansion for the Riccati-Bessel function j(l,rho)
C  j(l,rho) = rho^(l+1) sum(k=0,oo) (-1)^k 2^k (rho/2)^(2k)/k!/(2(k+l)+1)!!
         kmax = 100
         rho=w*x
         if (rho.gt.1e2) then
            appf1 = 0.0
            return
         endif
         ff = 1d0
         do i=1,ln
            ff = ff * rho / float(2*i+1)
         end do
         sum = ff
         summ = 0d0
         sump = ff
         zo2k = 1d0
         k = 1
         sumold = 0d0
         do while (k.lt.kmax.and.abs(sumold/sum-1d0).gt.1d-6)
            zo2k = zo2k  * rho * rho / 2d0 / float(k) /
     >         float(2*(ln+k)+1)
            sumold = sum
            if (mod(k,2).eq.0) then
               sump = sump + zo2k * ff
            else
               summ = summ + zo2k * ff
            endif 
            sum = sump - summ
c$$$            print '(i4,3e20.14)', k, sum, sump, summ
            k = k + 1
            if (abs(summ/sump-1d0).lt.1d-12) then
               k = kmax
               appf1 = 0.0
               if ( myid == 0) 
     >              print'("Precision loss in APPF1, returning 0.0",
     >            1pe14.4)',sum
               return
            endif 
         enddo
         acc = abs(summ/sump-1d0)
         if ( myid == 0) then
            if (k.eq.kmax)
     >           print'("Possible precision loss in", 
     >           "APPF1;result and error",1p,2e14.4)',sum,acc
         end if
         appf1 = rho * sum
      else
C  The following is the limiting case of the Riccati-Bessel/w**(ln+1) for w=0.0
         iprod = 2 * ln + 1
         iterm = 2 * ln + 1
         do while (iterm.gt.3)
            iterm = iterm - 2
            iprod = iprod * iterm
         enddo 
         appf1=x**(ln+1)/float(iprod)
      end if 
      return
      end

C  This function returns the index of the first X for which the regular
C  solution is < REGSTART.
      function jfirst1(l,ecmn,gridx,nx,jdouble,njdouble,regcut)
      use MPI_module
      implicit real*8 (a-h,o-z)

      dimension gridx(nx),jdouble(njdouble)
      w=sqrt(abs(ecmn))
      if (w.ne.0.0) then
         tmp=0.0
         do i=1,l
            tmp=tmp+log(dble(2*i+1))
         end do
      else
         tmp=log(dble(2*l+1))
         w=1.0
      end if
 10   xstart=exp((tmp+log(regcut))/(l+1))/w
      j=max(int(xstart*nx/gridx(nx)),1)
c      print*, 'j,xstart:', j,xstart, regcut
      j = min(j,nx)
      do while (gridx(j).lt.xstart.and.j.lt.nx)
         j=j+1
      end do
      s1 = appf1(l,w,gridx(j),acc)
      if (acc.lt.1e-10) then
         regcut = regcut / 10.0
         if ( myid == 0)
     >        print '("regcut redefined. S1, REGCUT:",1p,e13.4,e8.1)',
     $      s1,regcut
         go to 10
      endif 

C  make sure that we are not near a point where DX doubles

      do jd = 2, njdouble - 1
         if (abs(j-jdouble(jd)).lt.3) j = jdouble(jd) + 3
      end do       
      
c$$$      if (j.eq.nx) j = nx + 1
      jfirst1=j
      return
      end

C  This function returns the index of the last X for
C  which EXP(-WNN*X) < EXPCUT
      function jlast1(ecmn,gridx,nx,jdouble,njdouble,expcut)
      implicit real*8 (a-h,o-z)

      dimension gridx(nx),jdouble(njdouble)
      if (ecmn.ge.-0.005) then
         jlast1 = nx
      else
         wnn = sqrt(-ecmn)
         xmax=-log(expcut)/wnn
         j=min(int(xmax*nx/gridx(nx)),nx)
         do while(gridx(j).gt.xmax.and.j.gt.1)
            j=j-1
         end do
         do while (gridx(j)**(1.0/wnn)*exp(-gridx(j)*wnn).gt.expcut
     >      .and.j.lt.nx)
            j=j+1
         end do 

C  make sure that we are not near a point where DX doubles
         do jd = 2, njdouble - 1
            if (abs(j-jdouble(jd)).lt.3) j = jdouble(jd) + 3
         end do       

         jlast1=j
      end if 
      return
      end
      
C  The following routine solves the differential equation below backwards
C  from S(JSTOP)=S3 and S(JSTOP-1)=S2
C        2
C       D S   [ LN(LN+1)       V(X)          ]
C       --- - [ --------  + --------- - ECMN ] S =0
C         2   [   X*X           X            ]
C       DX      
C
C  Note that V(X)/X = UCENTR(X)
      subroutine numerovb(ln,ecmn,ucentr,cntfug,gridx,nx,
     >   jdouble,njdouble,s2,s3,chi,jstart,jstop)
      implicit real*8 (a-h,o-z)

      dimension ucentr(nx),cntfug(nx),gridx(nx),jdouble(njdouble),
     >   chi(nx)

C  We have to be careful not to lose precission in defining DX for large X.
      dx= gridx(jstop)-gridx(jstop-1)
      n = nint(dx/gridx(1))
      dx = gridx(1) * float(n)
      
      h2= dx*dx
      h2d= h2/12.d0
      xl= gridx(jstop)
      xlm1= gridx(jstop-1)
      wnn=sqrt(abs(ecmn))

      f3= ucentr(jstop)+cntfug(jstop)-ecmn
      f2= ucentr(jstop-1)+cntfug(jstop-1)-ecmn
      t3= (1d0 -h2d*f3)*s3
      t2= (1d0 -h2d*f2)*s2

      chi(jstop)= s3
      chi(jstop-1)= s2

      istart=njdouble-1
      do while (istart.gt.1.and.jdouble(istart).gt.jstop)
         istart=istart-1
      end do
      istop=istart
      istart=istart+1
      do while (istop.gt.1.and.jdouble(istop).gt.jstart)
         istop=istop-1
      end do
      istop=istop+1
C  JDOUBLE(ISTOP) points to the first doubling of DX that happens after JSTOP
C  JDOUBLE(ISTART) points to the last doubling of DX that happens before JSTART
      do i=istart,istop,-1
         j1=min(jstop,jdouble(i))-2
         j2=max(jstart,jdouble(i-1))
         do j=j1,j2,-1
            f1= ucentr(j)+cntfug(j)-ecmn
            t1= 2.0d0*t2 +h2*f2*s2 -t3
            s1= t1/(1.0d0- h2d*f1)
            t3=t2
            t2=t1
            f3=f2
            f2=f1
            s3=s2
            s2=s1
            chi(j)= s1
         end do
         if (j2.eq.jstart) return
         j=j2-1
         dx= dx/2.0d0
         h2= dx*dx
         h2d= h2/12.0d0
         f1= ucentr(j)+cntfug(j)-ecmn
         s1= s2*(36.0d0 +33.0d0*h2*f2) +s3*(-12.0d0+5.0d0*h2*f3)
         s1= s1/(24.0d0 +2.*h2*f1)
         t2= s2*(1. -h2d*f2)
         t1= s1*(1.0d0 -h2d*f1)
         t3=t2
         t2=t1
         f3=f2
         f2=f1
         s3=s2
         s2=s1
         chi(j)= s1
      end do
      
      return
      end

C  This routine matches the input REG solution to the apropriate
C  boundary conditions. It matches at two points of REG to the coulomb
C  functions F and H = G+iF. The two points are chosen where REG is not
C  too small.
      subroutine matchcc(l,eta,ecmn,jstart,jstop,gridx,jm,reg,
     >   phase,sigc,j1,j2,pshift,DELTA,ldw)
      use MPI_module
      implicit real*8 (a-h,o-z)

      dimension reg(jstop),gridx(jstop)
      real*8 dx, deta, xlmin, fc(0:1), gc(0:1), fcp(0:1), gcp(0:1)
      complex phase,sigc
      complex*16 xx,eta1,zlmin,cfc(9),cgc(9),cfcp(9),cgcp(9),sig(9),ci,
     >   f1,f2,f3,h1,h2,h3
! Mark: Added for COUL90
      real*8, dimension(0:100):: FCR, GCR, FCPR, GCPR
      real*8:: DELTAC
      integer:: ANG_L

      data pi,ci/3.1415927,(0d0,1d0)/

      if ( l >= 100 ) then
         if ( myid == 0) 
     >        print*,"Specified array size for Coulomb waves L>100"
         stop
      end if

!Mark: Add the following condition for repulsive waves. No need for normalizing 
!      if (jstart.eq.jstop) then
      if (jstart == jstop .AND. eta <= 0d0 ) then
         phase = (1.0d0,0.0d0)
         sigc = (1.0d0,0.0d0)
         return
      else if (jstart == jstop .AND. eta > 0d0 ) then
         sig(1) =  DELTAC(eta,dble(l))
         phase = (1.0d0, 1.0d0)
         sigc = exp ( ci * sig(1) ) 
         DELTA = sig(1)
         return
      endif
c$$$      if (reg(jstart).eq.0.0.and.reg(jstop).eq.0.0.and.reg(jm).eq.0.0)
c$$$     >   then
c$$$         phase = (1.0,0.0)
c$$$         sigc = (1.0,0.0)
c$$$         return
c$$$      endif
      mode1=12
      kfn=0
      deta = dble(eta)
      xlmin = dble(0)
      lrange = l
      if (ecmn.gt.0.0) then
         eta1=cmplx(dble(eta),0d0)
      else
         eta1=cmplx(0d0,-dble(eta))
      end if 
      zlmin=cmplx(dble(l),0d0)
      nl=2
      wnn=sqrt(abs(ecmn))
      if (wnn.lt.1e-15) wnn=1.0d0
      dx=5.0d0*pi/4.0d0/wnn

c$$$      j=min(jm+40,jstop)
c$$$C  Find the smallest value of REG closest to REG(JM)
c$$$      do while (abs(reg(j)).gt.abs(reg(j-1)).and.j-jm.gt.20)
c$$$         j=j-1
c$$$      end do
c$$$
c$$$C  Find the previous largest value of REG
c$$$      do while (abs(reg(j)).lt.abs(reg(j-1)).and.jm-j.lt.40)
c$$$         j=j-1
c$$$      end do




! Mark: I have changed the matching point to the larger valus of r.
! Somewhere past the matching and turning point are well behaved.
! This was because COULCC and FCOUL break down at the smaller values 
! of rho, for larger values of eta. For higher energy functions 0 <eta < 0.1
! I find that the numerov.f routine may propogate errors at large values of r.

! Mark: There is also a solution below for the high energy functions. We test
! the distorting phase is zero, if it is not, it loops back to point 60, hence
! does matching just outside the asymptotic region of the distorting potential.


      if (eta > 0d0) then
         turnpoint = eta + SQRT(eta * eta + (dble(l)*dble(l+1))) 
         do j = jstart, jstop
            tempr = gridx(j)*SQRT(ABS(ecmn))-turnpoint
            if ( tempr > 0d0 ) exit
         end do
         if ( j > jm ) jm = j
      end if

 60   if ( jm >= jstop .AND. eta > 0d0 ) then
         jm = jstop - 41
      end if

      j = jm
      j1 = jm
      j2 = jm
      regmax = reg(jm)
      regmin = abs(reg(jm))
      do while (j.lt.min(jm+40,jstop))
         j = j + 1
         if (reg(j).gt.regmax) then
            regmax = reg(j)
            j1 = j
         endif
         if (abs(reg(j)).lt.regmin) then
            regmin = abs(reg(j))
            j2 = j
         endif
      enddo
      if (j1.eq.j2) j1 = j2 - 2
c$$$      print*,'JM, J1, J2:', jm,j1,j2
C  First matching point
      if (j1.lt.jstart) j1=min(jstart+1,jstop)
      ifai=0
 10   xx=sqrt(cmplx(dble(ecmn)))*dble(gridx(j1))
      dx = real(xx)

      if (abs(ecmn).lt.1e-30) then
C  Here to evaluate the Coulomb G.F. at zero energy
         xx=dcmplx(sqrt(8.0d0*gridx(j1)))
         kfn=2
         zlmin=2*l+1
      end if 


! Mark: If the coulomb wave is zero within the distorting potential region 
! then there is no distorted wave phase shift.  There is also no need for 
! normalization (scaling)

! Mark: Using  subroutine COUL90. Far better than COULCC      
         xx = sqrt(cmplx(dble(ecmn)))*dble(gridx(j1))

         if ( eta < 0 ) then
            call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,
     >           kfn,ifai)

         else if ( eta >= 0 ) then
            rho = REAL(xx)  
            L = dble(l)
            ANG_L = INT(L)
            call coul90(rho, eta, dble(L), 0, FCR,GCR,FCPR,GCPR, 
     >           kfn,ifai)

            sig(1) =  DELTAC(eta,dble(L))
            cfc(1) = CMPLX(FCR(ANG_L))
            cgc(1) = CMPLX(GCR(ANG_L))+ci*CMPLX(FCR(ANG_L))  
         end if

c$$$      call coul90(dx, deta, xlmin, lrange, fc, gc, fcp, gcp, kfn, ifai)
c$$$      print*,'CFC,fc,CGC,GC:',cfc(1),fc(l),cgc(1),gc(l)
c$$$      print*,'|CGC-(GC+iFC)|:',abs(cgc(1)-cmplx(gc(l),fc(l)))
         if (ifai.ne.0.and.j1.ne.jstop) then
            ifai=0
            j1=jstop
            j = j1
            goto 10
         else if (ifai.ne.0) then
             if ( myid == 0) 
     >           print*,'IFAIL, ECMN, J1 in routine MATCH:',
     >           ifai,ecmn,j1
            jstart=jstop
            jstop =jstop -1
            phase = (1.0d0,0.0d0)
            do j = 1, jstop
               reg(j) = 0.0
            enddo 
            return
         end if 
         if (abs(ecmn).gt.1e-30) then
            f1=cfc(1)
            h1=cgc(1)
c$$$  print*,j1,cgc(1) * cfcp(1) - cgcp(1) * cfc(1)
         else if (eta.ne.0.0) then
            f1=cfc(1)*sqrt(gridx(j1))*pi
            h1=cgc(1)*sqrt(gridx(j1))*pi
         else
            f1=appf1(l,0d0,gridx(j1),acc)  
            h1=cmplx(gridx(j1)**(-l),real(f1))
         end if 
         
c$$$  C  Choose the other matching point
c$$$  do while (gridx(j)+dx.gt.gridx(j1).and.j1-j.lt.20)
c$$$  j=j-1
c$$$  end do
c$$$  j2=j
         if (j2.lt.jstart.or.j1.eq.j2) j2=j1-1
 20      xx=sqrt(cmplx(dble(ecmn)))*dble(gridx(j2))
         if (abs(ecmn).lt.1e-30) then
C     Here to evaluate the Coulomb G.F. at zero energy
            xx=dcmplx(sqrt(8.0*gridx(j2)))
            kfn=2
            zlmin=2*l+1
         end if 

! Mark: Using subroutine COUL90
         if ( eta < 0 ) then
            call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,
     >           kfn,ifai)
            
         else if ( eta >= 0 ) then
            
            rho = REAL(xx)  
            L = dble(l)
            ANG_L = INT(L)
            call coul90(rho, eta, dble(L), 0, FCR,GCR,FCPR,GCPR, 
     >           kfn,ifai)
            
            cfc(1) = CMPLX(FCR(ANG_L))
            cgc(1) = CMPLX(GCR(ANG_L))+ci*CMPLX(FCR(ANG_L))
         end if


         if (ifai.ne.0) then
            if (j2.ne.j1-1) then
               ifai=0
               j2 = j1 - 1
               j = j2
               goto 20
            else if (j1.ne.jstop) then
               ifai=0
               j1 = jstop
               j = j1
               goto 10
            else
               if ( myid == 0)
     >              print*,"IFAIL, ECMN, J1, J2 in",
     >              " routine MATCH:",ifai,ecmn,j1,j2
               jstart=jstop
               jstop =jstop -1
               phase = (1.0d0,0.0d0)
               do j = 1, jstop
                  reg(j) = 0.0d0
               enddo 
               return
            endif 
         end if 
         if (abs(ecmn).gt.1e-30) then
            f2=cfc(1)
            h2=cgc(1)
c$$$  print*,j2,cgc(1) * cfcp(1) - cgcp(1) * cfc(1)
         else if (eta.ne.0.0) then
            f2=cfc(1)*sqrt(gridx(j2))*pi
            h2=cgc(1)*sqrt(gridx(j2))*pi
         else
            f2=appf1(l,0d0,gridx(j2),acc)
            h2=cmplx(gridx(j2)**(-l),real(f2))
         end if 
         
C     Find real and imaginary parts of the phase
         phase=(f1*h2-f2*h1)/(reg(j1)*h2-reg(j2)*h1)
         const=abs(phase)
         phase=phase/const


! Mark: Testing accuracy of Coulomb regular irregular 
! functions for normailization. Should choose different matching point
! Expected phase to be 1 for L > ldw or when Wave start outside distorting
! region.
! This may mean the the Coulomb waves calculated by numerov are inaccurate
         if (l.gt.ldw) then
            diff = abs(real(phase)-1.0)
            if (diff.gt.1e-5) then
! Mark: Uncomment the next few lines if run to trouble with small positive eta.
!               jm = jmtemp
               if (myid==0 ) then
                  print*,"eta",eta,"L",l,"jstart"
                  print*,"Renormalizing, may need to increase qcut. 
     >                 Or problems with numerov.f or COUL90"
               end if
! This chooses the original jm as the matching point. Will continue running if uncomment
!               go to 60
            end if
         end if


c$$$  if (abs(phase+(1.0,0.0)).lt.1e-3) then
c$$$  phase = - phase
c$$$  const = - const
c$$$  endif
         if (reg(j1).eq.0.0.or.reg(j2).eq.0.0)
     >        stop 'REG(J1) or REG(J2) = 0 in MATCH'
         
         do j=jstart,jstop
            reg(j)=reg(j)*const
         end do
         
         y =  (reg(j1)*f2-reg(j2)*f1) / (f2*h1-f1*h2)
         y =  (reg(j1)*real(f2)-reg(j2)*real(f1)) /
     >        (real(f2)*real(h1)-real(f1)*real(h2))
         phase = cmplx(real(phase),y)
         if (abs((y-aimag(phase))/(abs(y)+1e-2)).gt.1e-2) then
            if (myid==0 )
     >           print*,'Warning different imaginary parts of phase',
     >           aimag(phase),y,(reg(j1)*f2-reg(j2)*f1) / (f2*h1-f1*h2),
     >           f2,h1,f1,h2,reg(j1),reg(j2),(f2*h1-f1*h2)
         endif 
c$$$  print*,'l, ecmn, const, phase:',l,ecmn, const, phase
c$$$  if (abs(phase-(1.0,0.0)).lt.1e-5) then
c$$$  if (abs(const-1.0).gt.0.2) then
c$$$  print*,'CAUTION: CONST <> 1.0; const,phase,ecmn,l',
c$$$  >         const,phase,ecmn,l
c$$$  c$$$            stop 'CONST <> 1.0'
c$$$  endif 
c$$$  endif
         
         
      sigc = exp(ci * sig(1))
c$$$      phase = phase * sigc
c$$$      jm = j1


!!$  D.Fursa



      DELTA = sig(1)   
      r1 = real(phase)/abs(phase)
      r2 = aimag(phase)/abs(phase)

      pshift = asin(r2)
      if(abs(pshift) .lt. 1e-10) pshift = 0d0
      if( r1 * cos(pshift) .lt. 0 ) then
         pshift = pi - pshift
      endif

! Mark: DELTA == Sigma_L and pshift == delta_L (distorted part)


!      print*,"eta",eta
!      print*,phase,r1,r2,pshift
!      print*, eta, pshift, DELTA

!!$


C     We now test the matching procedure by calling COULCC again
      j3=(j1 + j2) / 2
      xx=sqrt(cmplx(dble(ecmn)))*dble(gridx(j3))
      if (abs(ecmn).lt.1e-30) then
C     Here to evaluate the Coulomb G.F. at zero energy
         xx=dcmplx(sqrt(8.0d0*gridx(j3)))
         kfn=2
         zlmin=2*l+1
      end if 
      
! Mark: Using subroutine  COUL90
      if ( eta < 0 ) then
         call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,
     >        kfn,ifai)
         
      else if ( eta >= 0 ) then
         
         call coul90(REAL(xx), eta, dble(L), 0, FCR,GCR,FCPR,GCPR, 
     >        kfn,ifai)
         ANG_L = INT(L)      
         cfc(1) = CMPLX(FCR(ANG_L))
         cgc(1) = CMPLX(GCR(ANG_L))+ci*CMPLX(FCR(ANG_L))
      end if

      if (j3.ne.j1.and.j3.ne.j2.and.ifai.eq.0) then
         if (abs(ecmn).gt.1e-30) then
            f3=cfc(1)
            h3=cgc(1)
c$$$  print*,j3,cgc(1) * cfcp(1) - cgcp(1) * cfc(1)
         else if (eta.ne.0.0) then ! Mark: May need to change this?
            f3=cfc(1)*sqrt(gridx(j3))*pi
            h3=cgc(1)*sqrt(gridx(j3))*pi
         else
            f3=appf1(l,0d0,gridx(j3),acc)
            h3=cmplx(gridx(j3)**(-l),real(f3))
         end if
         test = abs(phase-(f1*h3-f3*h1)/(reg(j1)*h3-reg(j3)*h1))
         if (test.gt.1e-2) then
            if (myid==0 ) then
               print*,""
               print*,'Matching process has problems:',
     >              test, l, ecmn,j1,j2,j3,phase,
     >              (f1*h3-f3*h1)/(reg(j1)*h3-reg(j3)*h1)
               print*, eta, "ecmn",ecmn
               print*,""
            endif ! myid 
         end if
      endif
      return
      end




C  Make the Green's function on GRIDX
      subroutine makegreen(ln,eproj,eta,ucentr,cntfug,gridx,nx,
     >   jdouble,njdouble,regcut,expcut,reg,beta,jstart,jstop,phase)
      implicit real*8 (a-h,o-z)

      dimension ucentr(nx),cntfug(nx),gridx(nx),jdouble(njdouble),
     >   reg(nx),beta(nx)
      complex*16 phase,xx,eta1,zlmin,cfc(2),cgc(2),cfcp(2),
     >   cgcp(2),sig(2)

      ldw = -1
      call regular(ln,eproj,eta,ucentr,cntfug,ldw,gridx,nx,
     >   jdouble,njdouble,regcut,expcut,reg,jstart,jstop,pshift,DELTA)

      if (jstart.ge.jstop) return
      
      xlmin=ln
      zlmin=cmplx(dble(xlmin),0d0)
      phase=cmplx(pshift,0d0)
      ifail=1
      nl=1
      
C  Now to find the iregular solution. We take care so as not to start
C  at values that are too small. Use routine COULCC to get the G function
! Mark: COULCC is not good. Could be replaced with subroutine COUL90
      j=jstop
      if (eproj.lt.0.0) then
         phase=phase/(0.0d0,1.0d0)
         mode1=12
         kfn=0
         eta1=eta/(0.0,1.0)
         do while (j.ge.jstop-1)
            xx=sqrt(cmplx(dble(eproj)))*dble(gridx(j))
            call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,
     >         kfn,ifail)
            beta(j)=abs(cgc(1))
            j=j-1
         end do 
         phase=phase*cgc(1)/abs(cgc(1))
         rph=real(phase)
         aph=aimag(phase)
         if (abs(aph).gt.1e-5) then
            print*,'Problems in making Green''s functions'
            print*,'Imaginary part of phase is too big, expect 0',aph
         end if 
      else if (eproj.eq.0.0) then
         do while (j.ge.jstop-1)
            if (eta.eq.0.0) then
               beta(j)=gridx(j)**(-ln)
            else
C  Here to evaluate the Coulomb G.F. at zero energy
               xx=dcmplx(sqrt(8.0*gridx(j)))
               mode1=2
               kfn=2
               xlmin=2*ln+1
               zlmin=cmplx(dble(xlmin),0d0)
               call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,
     >            mode1,kfn,ifail)
               beta(j)=-real(cgc(1))*sqrt(gridx(j))
            end if
            j=j-1
         end do 
      else
         jj=jstop
         do while (j.ge.jstop-1.or.abs(beta(jj)).lt.0.1)
            xx=sqrt(cmplx(dble(eproj)))*dble(gridx(j))
            eta1=cmplx(dble(eta),0d0)
            mode1=2
            kfn=0
            call coulcc(xx,eta1,zlmin,nl,cfc,cgc,cfcp,cgcp,sig,mode1,
     >         kfn,ifail)
c$$$            print*,'G.F. test',cfc(1)*cgc(2)-cfc(2)*cgc(1)
            beta(j)=real(cgc(1))
            jj=j
            j=j-1
         end do
      end if 
      
      j=j+2
      s2=beta(j-1)
      s3=beta(j)
      call numerovb(ln,eproj,ucentr,cntfug,gridx,nx,
     >   jdouble,njdouble,s2,s3,beta,jstart,j)
      t2=1.0/real(phase)
      t=t2*aimag(phase)
      if (abs(t).gt.1e-4) then
         DO J=JSTART,JSTOP
            BETA(J)=beta(j)*t2-t*reg(j)
         end do
      else
C  T2 may be -1.0         
         DO J=JSTART,JSTOP
            BETA(J)=beta(j)*t2
         end do
      end if 
      RETURN
      END

      subroutine getnewal2(small,diffmin,diffmax,alstep,position,
     >   etot,psen2,nps,diff,it,al)
      implicit real*8 (a-h,o-z)

      dimension psen2(nps)
      n = 1 ; nsmall = 0
      almin = al; almax = al   ! Probably not right, but just to initialise.
      small = 1e10
      do while (psen2(n).lt.position.and.n.lt.nps)
         n = n + 1
      enddo
      distance = psen2(n) - psen2(n-1)
      nf = n - 1
      if (psen2(n).lt.position) nf = n
      slowery = 
     >   (position + distance / 2.0 + 4.0*psen2(n))/5.0
      do n = 1, nps
         diff = (slowery-psen2(n)) / slowery
         if (abs(diff).lt.small) then
            small = abs(diff)
            nsmall = n
c$$$                     print*, n, diff, almin, al, almax, psen2(n) * ry
         endif
      enddo
               
      if (psen2(nps).lt.etot) then
c$$$         al = al + alstep*it
         al = al - alstep
         small = 0.0
      else 
         diff = (slowery-psen2(nsmall)) / slowery
         if (diff.lt.0.0) then
            almax = al
            diffmax = diff
            if (diffmin.gt.0.0) then
c$$$  al = (almax + almin) / 2.0
               al = (almax * diffmin - almin * diffmax) /
     >            (diffmin - diffmax)
            else
               al = al - alstep
               almin = al
            endif 
         else
            almin = al
            diffmin = diff
            if (diffmax.lt.0.0) then
               al = (almax * diffmin - almin * diffmax) /
     >            (diffmin - diffmax)
c$$$  al = (almax + almin) / 2.0
            else
               al = al + alstep
               almax = al
            endif 
         endif
      endif
      print*,'reset AL, slowery to:', al,slowery,it
      return
      end
      
