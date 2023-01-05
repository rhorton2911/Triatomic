module Positronium

  real*8,public:: faclog(1000) 
  integer, parameter, public:: lamax1=20, ncmax1=700
!  real*8,dimension(:,:), allocatable:: ps2
!  real*8,dimension(:), allocatable::  psen2
 real*8:: cknd(ncmax1,ncmax1,0:lamax1),rlambda(2,0:lamax1) 
 integer:: npsstates(2,0:lamax1) 
  public:: Ps_structure

  type Ps_st_data    
     integer, dimension(:,:), allocatable:: istoppsinb 
     real*8, dimension(:,:,:), allocatable:: psinb
     real*8, dimension(:,:), allocatable::   enpsinb,gnl0,gnl1,f
     real*8, dimension(:), allocatable:: en, l,n, m,maxf
     integer:: Nmax, nchm, nchm_open   
     real*8, dimension(:),allocatable:: ve2st,va2st ! static potentials from H2+ ion
     integer:: iglp
     real*8, dimension(:),allocatable:: gp,xgz,wgz,xgp,wgp
     real*8, dimension(:,:), allocatable:: pl(:,:)
     integer:: igpm,igzm,nam_max
     real*8, dimension(:),allocatable:: Amg
     real*8, dimension(:,:,:),allocatable:: Bnlp,Dnlp
     real*8, dimension(:),allocatable:: factrl,sqrfct,hat
     integer, dimension(:), allocatable:: l_of_t,icount_l
     integer, dimension(:,:),allocatable:: l_nsti

  end type Ps_st_data

  type(Ps_st_data):: Ps


contains   
      

subroutine Ps_structure

  use input_data
  use MPI_module
  use grid_radial  
 
  implicit none  
  integer, parameter:: maxl=2*20+1
  real*8:: enpsinb_temp, tsum,ry,au
  integer n, npos, lp, i, Nmax,nst , iglp,ip,mp
  real*8:: chi(grid%nr),rr,arg,besj,term,res0,res1,pp,ppll,pp2
  real*8:: uplane(grid%nr)
  real*8::  ps2(grid%nr,ncmax1),psen2(ncmax1)
  integer,dimension(ncmax1):: minps2,maxps2 
  integer k,nstart,n1
  real*8:: ckndtmp(ncmax1,ncmax1)
         uplane(:)=0.d0
  
        au=data_in%eV  !  27.2116
        ry=au/2.d0

     Ps%igzm =  data_in%igz
     Ps%igpm = data_in%igp

     allocate(Ps%xgz(Ps%igzm),Ps%wgz(Ps%igzm),Ps%pl(0:20,Ps%igzm), Ps%xgp(Ps%igpm),Ps%wgp(Ps%igpm))
     call polleg(20,Ps%igzm,Ps%igpm,Ps%xgz,Ps%wgz,Ps%pl,Ps%xgp,Ps%wgp)
     allocate(Ps%l_of_t(1000))

       iglp = 4000
       Ps%iglp = iglp
       allocate(Ps%gp(Ps%iglp))
       do i = 1, iglp
        Ps%gp(i) = 1.e-8*exp(0.0065*i)
       enddo
  
 
      allocate(Ps%psinb(grid%nr,abs(data_in%nptop(data_in%lptop)),0:data_in%lptop))
      allocate(Ps%istoppsinb(abs(data_in%nptop(data_in%lptop)),0:data_in%lptop))
      allocate(Ps%enpsinb(abs(data_in%nptop(data_in%lptop)),0:data_in%lptop))
      allocate(Ps%ve2st(grid%nr),Ps%va2st(grid%nr))


      Ps%ve2st(:) = 0.d0; Ps%va2st(:) = 0.d0
     
      Nmax=0

      do lp = data_in%lpbot, data_in%lptop
        print*,'Defining positronium eigenstates for LP:',lp,data_in%nptop(lp)
         print*,' l Np N      e(Ry)     e(eV)     ovlp    maxr' 

            nstart = lp+1 ! 
            n=0
         do npos =data_in%npbot(lp), abs(data_in%nptop(lp))!
            n = n + 1
              do mp = -lp, lp
                 Nmax = Nmax+1
              enddo
!            call rnl(0,npos,lp,Ps%psinb(1:grid%nr,n,lp),enpsinb_temp,&
!     &         Ps%istoppsinb(n,lp),grid%nr,grid%gridr,grid%expcut)
!            Ps%enpsinb(n,lp)= enpsinb_temp/2.d0 ! a.u. 
            npsstates(2,lp) = data_in%npsp(lp)-lp
            rlambda(2,lp) = abs(data_in%alphap(lp) * 2.0)
            ps2(:,:) = 0.0
            psen2(:) = 0.0
            ckndtmp(:,:)=0.d0


              call makeps(0.5,.false.,abs(data_in%alphap(lp)),lp,grid%expcut,data_in%npsp(lp)-lp,ps2,psen2, &
     &                    minps2, maxps2,grid%gridr,uplane, 1, grid%nr,ckndtmp)

          do k=1,npsstates(2,lp)
           do n1=1,ncmax1
          cknd(k,n1,lp)=ckndtmp(k,n1)
          enddo ! n1
!!          print*,'Ps::a(k):',k,npos,lp,ckndtmp(k,npos),cknd(k,npos,lp)
          end do


             Ps%psinb(1:grid%nr,n,lp)=ps2(1:grid%nr,npos-lp)
             Ps%enpsinb(n,lp) = psen2(npos-lp) /2.d0 ! in a.u.
             Ps%istoppsinb(n,lp) = maxps2(npos-lp)

            tsum = 0.0
            do i = 1, Ps%istoppsinb(n,lp)
               tsum = tsum + Ps%psinb(i,n,lp) * Ps%psinb(i,n,lp) * grid%weight(i)
            enddo
            print'(3i3,3f10.5,i8)',lp,npos,n,Ps%enpsinb(n,lp),&
     &         Ps%enpsinb(n,lp)*au,tsum,Ps%istoppsinb(n,lp)

         enddo

      enddo 
   
       allocate(Ps%en(Nmax),Ps%n(Nmax),Ps%l(Nmax),Ps%m(Nmax))
       allocate(Ps%gnl0(iglp,Nmax),Ps%gnl1(iglp,Nmax),Ps%f(Nmax,1:grid%nr),Ps%maxf(Nmax))

       nst=0
       do lp = data_in%lpbot,data_in%lptop
        n=0 
         do npos =data_in%npbot(lp), abs(data_in%nptop(lp))
          n=n+1
           do mp = -lp, lp, 1
              nst=nst+1
              Ps%en(nst)= Ps%enpsinb(n,lp)
              Ps%n(nst)=npos
              Ps%l(nst)=lp
              Ps%m(nst) = mp
              Ps%maxf(nst) = Ps%istoppsinb(n,lp)
              Ps%f(nst,1:grid%nr) = Ps%psinb(1:grid%nr,n,lp)
           enddo !mp
         do ip=1,Ps%iglp,20
            pp2=Ps%gp(ip)
            pp=dsqrt(pp2)

            res0=0.d0; res1=0.d0
         do i = 1, Ps%istoppsinb(n,lp)
            rr = grid%gridr(i)
            arg = rr * pp

            call sbessel(arg, lp, besj)
            term = Ps%psinb(i,npos,lp) * besj * grid%weight(i)
            res0 = res0 + term  ! CHECK SIGN, + or - ?
            res1 = res1 + term * rr
! if(nst==1) print*,'Ps-wf(r):', rr, Ps%psinb(i,npos,lp)/rr, 1.d0/sqrt(2.d0)*exp(-rr/2.d0),i
! if(nst==2) write(*,'(a8,3e20.10,i6)')'Pswfr:', rr, Ps%psinb(i,npos,lp)/rr, 1.d0/4*(1d0-rr/4d0)*exp(-rr/4d0),i
         enddo
!            if(nst.eq.2) stop'nst cycle'
           ppll = pp ** lp
           res0 = res0 / ppll
           res1 = res1 / ppll

           Ps%gnl0(ip,nst) = res0
           Ps%gnl1(ip,nst) = res1

!        if(nst.eq.1) then
!         res0 = dsqrt(8.d0)/(4*pp2+1.d0)/ ppll
!         res1 = dsqrt(128.d0)/(4*pp2+1.d0)**2.d0 / ppll
!        write(*,'(a6,5e20.10)')'Pswf:',pp2,Ps%gnl0(ip,nst),res0,Ps%gnl1(ip,nst),res1
!        endif
!        if(nst.eq.2) then
!         res0 = (64*pp2-4)/(16*pp2+1.d0)/ ppll
!         res1 = (64*(16*pp2-1))/(16*pp2+1.d0)**3.d0 / ppll
!        write(*,'(a6,5e20.10)')'Pswf:',pp2,Ps%gnl0(ip,nst),res0,Ps%gnl1(ip,nst),res1
!        endif

          enddo !ip
         enddo
!            npsstates(2,lp) = data_in%npsp(lp)-lp
!            rlambda(2,lp) = abs(data_in%alphap(lp) * 2.0)

       enddo

       Ps%Nmax = Nmax

      allocate(Ps%sqrfct(0:maxl),Ps%factrl(0:maxl), Ps%hat(0:maxl))
  
      call fctrls(Ps%factrl,Ps%sqrfct,Ps%hat)
        
     return
     END SUBROUTINE Ps_structure
 

!  The following routine returns the hydrogenic radial functions.
      subroutine rnl(nznuc,n,l,chi,en,jstop,maxr,rmesh,expcut)
      implicit real*8 (a-h,o-z)
      implicit integer (j-n)
!      include 'par.f'
!      COMMON/MESHRR/ MESHR,RMESH(MAXR,3)
!      common/smallr/ formcut,regcut,expcut

      real*8 rmesh(maxr),rweight(maxr),chi(maxr) 
      real*8  rlrho, rho, const2, sumpos, sumneg
!      intrinsic fact

       meshr=maxr
      if (nznuc.gt.0) then
         Z = float(nznuc)
         en = - (Z / float(n)) ** 2  ! Ry
      else
!  Here we define positronium          
         Z = 0.5
         en = - 0.5 / float(n)**2  ! Ry
      endif 
      do j=1,meshr
         chi(j)=0.0
      end do 
      const1 = sqrt(Z * fact(n-l-1) * fact(n+l)) / float(n)
      rho = 0.0
      exprho = exp(- rho / 2.0)
      j = 0
      do while (exprho.gt.expcut.and.j.lt.meshr)
         j = j + 1
         rlrho = 0.0
         rho = 2d0 * Z * rmesh(j) / float(n)
         exprho = exp(- rho / 2.0)
         do k = 0, n - l - 1
            const2=(-1)**k/fact(n-l-1-k)/fact(2*l+1+k)/fact(k)
            rlrho = rlrho + const2 * rho ** k
         enddo
         tmp = chi(j)
         chi(j) = const1 * rlrho * rho ** (l + 1) * exprho
      enddo
      jstop = j
      return
      end subroutine rnl

      real*8 function fact(k)
      integer k
      fact=1d0
      do i=2, k
         fact=fact*float(i)
      end do
      return
      end function fact


!!!!
      subroutine polleg(limit,igz,igp,xgzR,wgzR,plR,xgp,wgp)

      implicit real*8 (a-h,o-z)       
      implicit integer (i-n)
      real*8::  xgz(igz),wgz(igz),pl(0:limit,igz)
      integer:: igz,igp,igpa,igza,igzR
      real*8:: xgzR(igz),wgzR(igz),plR(0:limit,igz)
      real*8:: xgp(igp),wgp(igp)
      real*8:: endpts(2),b(igz)
      real*8:: pl8(0:20)

      real*8:: xga(igz), wga(igz)
      real*8:: xgb(igp), wgb(igp)
      real*8 :: x1, x2, z

!***** for z integral
      x1=-1.d0; x2=1.d0

!C     gauss-legendre quadrature is used
!C     gauleg utilizes either cgqf or the code from Numerical recepies

!C     print*, '   (1.1) polleg'
      call gauleg(x1,x2,xga,wga,igz)

      do ig = 1,igz
         z = xga(ig)
         xgz(ig) = z
         wgz(ig) = wga(ig)
         pl8(0)  = 1.0d0
         pl(0,ig) = 1.d0*wga(ig)
         pl(1,ig) = z*wga(ig)
         pl8(1)=z
         do i=2,limit
            pl8(i)=((2*i-1)*z*pl8(i-1)-(i-1)*pl8(i-2))/dble(i)
            pl(i,ig)=pl8(i)*wga(ig)
         end do
      enddo
!cRRRR***** for z integral
      zeps=0.03d0
      x1=-1.d0; x2=1.d0-zeps
      igz1=igz/4*3
      igz2=igz/4
      call gauleg(x1,x2,xga,wga,igz1)
      x1=x2; x2=1.d0

      do ig = 1,igz1
         z = xga(ig)
         xgzR(ig) = z
         wgzR(ig) = wga(ig)
         pl8(0)  = 1.0d0
         plR(0,ig) = 1.d0*wga(ig)
         plR(1,ig) = z*wga(ig)
         pl8(1)=z
         do i=2,limit
            pl8(i)=((2*i-1)*z*pl8(i-1)-(i-1)*pl8(i-2))/dble(i)
            plR(i,ig)=pl8(i)*wga(ig)
         end do
      enddo

      call gauleg(x1,x2,xga,wga,igz2)
         ig=0
      do ig2 = igz1+1,igz
         ig=ig+1
         z = xga(ig)
         xgzR(ig2) = z
         wgzR(ig2) = wga(ig)
         pl8(0)  = 1.0d0
         plR(0,ig2) = 1.d0*wga(ig)
         plR(1,ig2) = z*wga(ig)
         pl8(1)=z
         do i=2,limit
            pl8(i)=((2*i-1)*z*pl8(i-1)-(i-1)*pl8(i-2))/dble(i)
            plR(i,ig2)=pl8(i)*wga(ig)
         end do
      enddo

      x1=0.d0; x2=1.d0
      
      call gauleg(x1,x2,xgb,wgb,igp)
      
      do ig=1,igp
         xgp(ig)=xgb(ig)
         wgp(ig)=wgb(ig)
      enddo

      return
      end subroutine polleg


      subroutine gauleg(x1,x2,x,w,n)
      implicit real*8 (a-h,o-z)
      implicit integer(i-n)
      real*8:: x1,x2,x(n),w(n),wf(n*2),iwf(n*2)
      real*8, parameter:: eps=1.d-14

      pi = acos(-1d0)
      call cgqf(n,x,w,1,0d0,0d0,x1,x2,0,2*n,wf,2*n,iwf,ier)
      if (ier.eq.0) return
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
         z=cos(pi*(dble(i)-0.25d0)/(dble(n)+0.5d0))
 1       continue
         p1=1.d0
         p2=0.d0
         do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(dble(j)-1.d0)*p3)/dble(j)
 11      continue
         pp=n*(z*p1-p2)/(z*z-1.d0)
         z1=z
         z=z1-p1/pp
         if(abs(z-z1).gt.eps)go to 1
         x(i)=xm-xl*z
         x(n+1-i)=xm+xl*z
         w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
         w(n+1-i)=w(i)
 12   continue
      return
      end subroutine gauleg

!*------------------------------------------------------

      subroutine fctrls(factrl,sqrfct,hat) !Qltable

!      include 'par.f'
!      include 'par.pos'
      integer, parameter:: maxl=2*20+1
!      implicit real*8 (a-h,o-z)    
!c$$$      common/const/pi
      real*8, dimension(0:maxl):: factrl,sqrfct,hat
!     >   Qlfactor(0:ltmax)
!      common/fik/fik(0:19,0:19)
!      common/gausz/xgz(igzm),wgz(igzm),pl(0:maxl,igzm),igza
!      common/gauszR/xgzR(igzm),wgzR(igzm),plR(0:maxl,igzm),igzR
!       common/gausp/xgp(igpm),wgp(igpm),igpa
!      common/funQl/Qlarray(0:ltmax,ich)
!      common/cheb/x(ich),w(ich)
!      real*8 dfactrl(0:ltmax),dfactrl2(0:ltmax),arg(ich)
      
!c      open(99,file='matritsa')

      pi = acos(-1d0)
!      if(igz.gt.igzm .or. igp.gt.igpm) then
!         print*, igz,igzm,igp,igpm
!         stop'increase igzm and/or igpm in par.pos'
!      endif
!      igza=igz
!      igpa=igp
!      igzR=igza
!
!      do i=1,ich
!         t=x(i)
!         arg(i)=(t*t+1.d0)/t/2.d0
!      end do

      factrl(0)=1.d0
      hat(0)=1.d0
      sqrfct(0)=1.d0
      do i=1,maxl
         factrl(i)=factrl(i-1)*dble(i)
         hat(i)=sqrt(dble(2*i+1))
         sqrfct(i)=sqrt(factrl(i))
      end do

!      do k=0,19
!         do i=0,19
!            fik(i,k)=factrl(k+i+1)/(factrl(i+2)*factrl(i))/
!     >         (factrl(k+2)*factrl(k))
!         enddo
!      enddo    
!      print*
!      print*,' posVmat: using ',igz,' points for z and ',
!     >   igp,'*nqmi*2 for p integral'
!      print*,' Forming tables of Ql and Pl'
!      dfactrl(0)=1.d0
!      dfactrl2(0)=1.d0
!      Qlfactor(0)=1.d0
!      do i=1,ich
!         Qlarray(0,i)=1.d0
!      end do
!      do il=1,lstopm
!         dfactrl(il)=dfactrl(il-1) * dble(il)
!         dfactrl2(il)=dfactrl2(il-1) * dble(2*il+1)
!         Qlfactor(il)=dfactrl(il)/dfactrl2(il)
!*         do i=1,ich
!*            Q0= 0.5d0*log((arg(i)+1.d0)/(arg(i)-1.d0))
!*            call funleg(arg(i),il,resQl)
!*            Qlarray(il,i)=resQl/Q0
!*         end do
!      enddo
!      call polleg(2*lstopm+1,igz,igp)
!      print*
      return
      end subroutine fctrls !Qltable












end module Positronium 

!   include 'laguer.f'

!!include 'iqpackd.f' 
