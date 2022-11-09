!  This subroutine is used to calculate e-e Coulomb potential.
!  This routine returns  temp(r) = int dr' fun(r') * f(r<) * g(R>), where
!  the range of integration is zero to infinity. 
! fun(r) : is passed via argument list together with its start and last points, it already contains the Simson's integration weights.
! f(r) = r**l, and g(r) = 1/r**(l+1): are passed via module rpow.

! l : angular momentum
! nr : size of r-grid
! maxfm : index to the maximum r value for which temp(i) need to be calculated.

subroutine form(l,fun,minfun,maxfun,maxfm,temp,i1,i2)
  use grid_radial 
  implicit none
  
  integer, intent(in):: l
  real*8, dimension(grid%nr), intent(in):: fun
  integer, intent(in):: minfun, maxfun, maxfm
  real*8, dimension(grid%nr), intent(out):: temp
  integer, intent(out):: i1, i2
  
  real*8, dimension(0:grid%nr):: t1, t2
  integer:: i
  integer:: if1, if2, istop, it1max, it2min, i1a, i2a


  if(l .gt. grid%ltmax) then
     i1=2
     i2=1     
     return
  endif

! Set correct limits for integration: \int dr' fun(r') f(r')
  if1 = max(minfun,grid%i1_f(l))
  if2 = min(maxfun,grid%i2_g(l))


  if (if2 .le. if1 + 2) then
     i1=2
     i2=1
     return
  end if

! do not initialise array temp to zero for faster calculation:
!  temp = 0.0

  istop=min(if2,maxfm)

!  Find the integral for fun(r') * f(r') from 0 to istop. The function fun 
!  already contains the Simson's integration weights. 
! limits for which integral is required are:  if1:istop
         t1(if1 - 1) = 0.0
         do i=if1,istop
            t1(i) = t1(i-1) + fun(i)*grid%rpow_f(i,l)
         end do 
! account for the case when:   maxfun < it1max=min(maxfm,i2_g),
! note that in this case istop = maxfun
         it1max = min(maxfm,grid%i2_g(l))
         if(maxfun .lt. it1max) then
            t1(istop+1:it1max) = t1(istop)      
         else
            it1max = istop
         endif
! limits now became:  if1:it1max
!         write(*,'("t1=",3F15.8)') t1(maxfm),t1(maxfm-1), t1(maxfm-10)

!  Find the integral of fun(r') * g(r') from infinity to if1
! limits for which integral should be taken:  if1:if2
         t2(if2) = 0.0
         do i=if2,if1,-1
            t2(i-1) = t2(i) + fun(i)*grid%rpow_g(i,l)
         end do 
! account for the case when:   minfun > i1_f(l)
         if(minfun .gt. grid%i1_f(l)) then
            t2(grid%i1_f(l):if1-1) = t2(if1)
            it2min = grid%i1_f(l)
         else
            it2min = if1
         endif
! limits now became:  it2min:if2, but the upper limit is istop, so the final limits are:
!  it2min:istop

!  Make the form factor by summing two parts

         temp(if1:it1max) = t1(if1:it1max)*grid%rpow_g(if1:it1max,l)

! if array tmp is initialised to zero then:
!         temp(it2min:istop) = temp(it2min:istop) +  t2(it2min:istop)*rpow_f(it2min:istop,l)
! if not:
         if(if1 .le. it2min) then
            if(it1max .ge. istop) then
               temp(it2min:istop) = temp(it2min:istop) +  t2(it2min:istop)*grid%rpow_f(it2min:istop,l)
            else
               temp(it2min:it1max) = temp(it2min:it1max) +  t2(it2min:it1max)*grid%rpow_f(it2min:it1max,l)
               temp(it1max:istop) =  t2(it1max:istop)*grid%rpow_f(it1max:istop,l)
            end if
         else
            if(it1max .ge. istop) then
               temp(if1:istop) = temp(if1:istop) +  t2(if1:istop)*grid%rpow_f(if1:istop,l)
               temp(it2min:if1) =  t2(it2min:if1)*grid%rpow_f(it2min:if1,l)
            else
               temp(if1:it1max) = temp(if1:it1max) +  t2(if1:it1max)*grid%rpow_f(if1:it1max,l)
               temp(it2min:if1) =  t2(it2min:if1)*grid%rpow_f(it2min:if1,l)
               temp(it1max:istop) =  t2(it1max:istop)*grid%rpow_f(it1max:istop,l)
            end if
         endif
! limits: 
         i1 = min(if1,it2min)
         i2 = max(it1max,istop)

         i1a = i1
         i2a = i2
         call minmaxi_form(temp,grid%nr,i1a,i2a)
         
         i1 = max(i1,i1a)
         i2 = min(i2,i2a)

end subroutine form
!
!!$-----------------------------------------------------------------
!

subroutine form_accurate(l,fun,minfun,maxfun,maxfm,temp,i1,i2)
  !Liam wrote more accurate form subroutine based on Igor's routine in the atomic code
  use grid_radial 
  implicit none
  
  integer, intent(in):: l
  real*8, dimension(grid%nr), intent(in):: fun
  integer, intent(in):: minfun, maxfun, maxfm
  real*8, dimension(grid%nr), intent(out):: temp
  integer, intent(out):: i1, i2
  
  real*8, dimension(0:grid%nr):: t1, t2
  integer:: i, n
  integer:: if1, if2, istart, istop
  real*8, dimension(:), pointer :: dr, dr_on_three, dr_three_on_eight
  real*8, dimension(grid%nr) :: ffun, gfun

  !radial grid spacings
  dr => grid%dr
  dr_on_three => grid%dr_on_three
  dr_three_on_eight => grid%dr_three_on_eight

  if(l .gt. grid%ltmax) then
     i1=2
     i2=1     
     return
  endif

  if1 = grid%i1_f(l) 
  if2 = grid%i2_g(l)
  !minfun, maxfun limits of fun
  !if1, if2 will be range of integration
  if1 = max(minfun,if1)
  if2 = min(maxfun,if2)
  if2 = min(if2, maxfm)

  !ffun = f(r)*fun(r) etc
  ffun(if1:if2) = grid%rpow_f(if1:if2,l)*fun(if1:if2)
  gfun(if1:if2) = grid%rpow_g(if1:if2,l)*fun(if1:if2)
 
  !integration indices:
    !integral is required over if1 to if2 based on radial limits of ffun and gfun
    !we want to iterate only over even indices so we set
    istart=(if1+1)/2*2! increase if1 by 1 if it is odd
    istop=if2/2*2     ! decrease if2 by 1 if it is odd
  
  if (istop - istart.lt.4) then
    i1=2
    i2=1
    return
  end if

  !The iteration below starts at istart+2 since each step refers to the two
  !below it.
  !If if1 is odd then the iteration starts at istart+2 = if1+3
  !We need to set t1(if1), and also t1(if1+1) and t1(if1+2) for the iteration.
  !If if1 is odd we only need to set t1(if1) and t1(if1+1) so 
  !We assume ffun(if1-1) = 0

  if(istart==if1) ffun(istart-1) = 0.0d0

  !Integrate fun(r')*P(r') from 0 to r -> store in t1
  t1(istart)=(4.0d0*ffun(istart-1)+ffun(istart))*dr_on_three(istart) !simpson with ffun(if1-1) = 0
  
  do i=istart+2, istop-2,2
    !even points: Simpson 1/3 rule:
      t1(i)=t1(i-2)+(ffun(i-2)+4.0d0*ffun(i-1)+ffun(i))*dr_on_three(i)
    !odd points: Simpson 3/8 rule:
      t1(i+1)=t1(i-2)+(ffun(i-2)+3.0d0*ffun(i-1)+3.0d0*ffun(i)+ffun(i+1))*dr_three_on_eight(i+1)
  enddo
  
  if(istart==if1+1) then !if1 is odd 
    t1(if1)=t1(if1+1)*0.75d0 - t1(if1+3)/8.0d0
    t1(if1+2)=t1(if1+1)*0.75d0 + t1(if1+3)*0.375d0 !t1(if1-1)=0
  else
    !if1 is even, ffun(istart-1) is zero
    t1(if1+1)=t1(if1)*0.75d0 + t1(if1+2)*0.375d0 !t1(if1-1)=0
  endif
  
  !set last point, t1(istop-1) was set in last iteration above
  t1(istop)=t1(istop-2)+(ffun(istop-2)+4.0d0*ffun(istop-1)+ffun(istop))*dr_on_three(istop)
  if(if2==istop+1) then !if2 is odd so t1(if2) has not been set
    t1(if2)=t1(if2-3)+(ffun(if2-3)+3.0d0*ffun(if2-2)+3.0d0*ffun(if2-1)+ffun(if2))*dr_three_on_eight(if2) !newton 3/8
  endif

  !above does not work when grid spacing has doubled.
  !1/3 rule works at doubling point since it is always at an even index.
  do n=2, grid%ndouble+1 !cycle through each grid doubling
    i=grid%jdouble(n)+1  !index of nth doubling
    if(i<istart.or.i>istop) cycle
    if(i>if2-2) then !can't refer to t1(i+2)
      t1(i)=t1(i-1)
      cycle
    endif
    !3/8 rule breaks down at the first index after doubling
    t1(i)=t1(i+2)-(ffun(i)+4.0d0*ffun(i+1)+ffun(i+2))*dr_on_three(i+2)
  enddo

  !Iterate backwards for integral from r to infinity
  if(istop .ne. if2) t2(if2) = 0.0d0
  t2(istop) = 0.0d0
  do i=istop,istart+2,-2
    t2(i-2)=t2(i)+(gfun(i-2)+4.0d0*gfun(i-1)+gfun(i))*dr_on_three(i)
    t2(i-3)=t2(i)+(gfun(i-3)+3.0d0*gfun(i-2)+3.0d0*gfun(i-1)+gfun(i))*dr_three_on_eight(i)
  enddo
  
  !fix issues at doubling points again
  do n=2, grid%ndouble+1
    i=grid%jdouble(n)-1
    if(i<istart) cycle
    t2(i)=t2(i-2)-(gfun(i-2)+4.0d0*gfun(i-1)+gfun(i))*dr_on_three(i)
    if(i==istart+1) t2(i)=t2(istart)
  enddo
  t2(istop-1)=t2(istop-3)-(gfun(istop-3)+4.0d0*gfun(istop-2)+gfun(istop-1))*dr_on_three(istop-1)

  temp(if1:if2) = t1(if1:if2)*grid%rpow_g(if1:if2,l) + t2(if1:if2)*grid%rpow_f(if1:if2,l) 
  
  i=if2
  do while (abs(temp(i)) > grid%formcut .and. i<min(maxfm,grid%i2_g(l)))
    i=i+1
    temp(i)=t1(if2)*grid%rpow_g(i,l)
  enddo
  i2=i
  do while(abs(temp(i2)) < grid%formcut .and. i2 > if1)
    i2=i2-1
  enddo
  i=if1
  do while(abs(temp(i))>grid%formcut .and. i > grid%i1_f(l))
    i=i-1
    temp(i)=t2(if1)*grid%rpow_f(i,l)
  enddo
  i1=i
  do while(abs(temp(i1)) < grid%formcut .and. i1 <= i2)
    i1=i1+1
  enddo
  
end subroutine form_accurate


subroutine minmaxi_form(f,nr,i1,i2)
  use grid_radial
  implicit none
  integer, intent(in):: nr
  real*8, intent(in), dimension(nr)::  f
  integer, intent(inout):: i1, i2

  integer:: i

  i=i1
  do while (i.lt.nr.and.abs(f(i)).lt.grid%formcut)
     i=i+1
  end do
  i1=i
  i=i2
  do while (i.gt.1.and.abs(f(i)).lt.grid%formcut)
     i=i-1
  end do
  i2=i
!  if(i2.lt.i1)i2=i1 !cb18.10.10
  if(i2.lt.i1) then  ! function is zero on the grid
     i1 = 2
     i2 = 1 
  endif
  return
end subroutine minmaxi_form

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

subroutine form_spheroidal(l,m,fun,minfun,maxfun,maxfm,temp,i1,i2)
  !
  ! Essentially the same as above but with additional m argument and replacing
  ! i1_f -> i1_P, i2_g -> i2_Q, rpow_f -> LegendreP, and rpow_g -> LegendreQ.
  !
  ! Unfortunately it was necessary to create this function as my initial idea of
  ! m being an optional argument doesn't work without making this into a module.
  !
  ! 15/10/12 by J.S.
  !
  use grid_radial 
  implicit none
  
  integer, intent(in):: l, m
  real*8, dimension(grid%nr), intent(in):: fun
  integer, intent(in):: minfun, maxfun, maxfm
  real*8, dimension(grid%nr), intent(out):: temp
  integer, intent(out):: i1, i2
  
  real*8, dimension(0:grid%nr):: t1, t2
  real*8, dimension(:), allocatable :: PVec, QVec
  integer:: i, n
  integer:: if1, if2, istop, it1max, it2min, i1a, i2a


  if(l .gt. grid%ltmax) then
     i1=2
     i2=1     
     return
  endif

  ! Translate the 2D (l,m) pair to the 1D (n=sum{l}+m) index.
  n = abs(m)
  do i = 1, l
     n = n + i
  end do

! Set correct limits for integration: \int dr' fun(r') f(r')
  if1 = max(minfun,grid%i1_P(n))
  if2 = min(maxfun,grid%i2_Q(n))
  allocate( PVec(1:grid%nr), QVec(1:grid%nr) )
  PVec(:) = grid%LegendreP(:,n)
  QVec(:) = grid%LegendreQ(:,n)

  if (if2 .le. if1 + 2) then
     i1=2
     i2=1
     return
  end if

! do not initialise array temp to zero for faster calculation:
!  temp = 0.0

  istop=min(if2,maxfm)


!  Find the integral for fun(r') * f(r') from 0 to istop. The function fun 
!  already contains the Simson's integration weights. 
! limits for which integral is required are:  if1:istop
         t1(if1 - 1) = 0.0
         do i=if1,istop
            t1(i) = t1(i-1) + fun(i)*PVec(i)
         end do 
! account for the case when:   maxfun < it1max=min(maxfm,i2_g),
! note that in this case istop = maxfun
         it1max = min(maxfm,grid%i2_Q(n))
         if(maxfun .lt. it1max) then
            t1(istop+1:it1max) = t1(istop)      
         else
            it1max = istop
         endif
! limits now became:  if1:it1max
!         write(*,'("t1=",3F15.8)') t1(maxfm),t1(maxfm-1), t1(maxfm-10)

!  Find the integral of fun(r') * g(r') from infinity to if1
! limits for which integral should be taken:  if1:if2
         t2(if2) = 0.0
         do i=if2,if1,-1
            t2(i-1) = t2(i) + fun(i)*QVec(i)
         end do 
! account for the case when:   minfun > i1_P(n)
         if(minfun .gt. grid%i1_P(n)) then
            t2(grid%i1_P(n):if1-1) = t2(if1)
            it2min = grid%i1_P(n)
         else
            it2min = if1
         endif
! limits now became:  it2min:if2, but the upper limit is istop, so the final limits are:
!  it2min:istop


!  Make the form factor by summing two parts
         temp(if1:it1max) = t1(if1:it1max)*QVec(if1:it1max)

! if array tmp is initialised to zero then:
!         temp(it2min:istop) = temp(it2min:istop) +  t2(it2min:istop)*LegendreP(it2min:istop,l)
! if not:
         if(if1 .le. it2min) then
            if(it1max .ge. istop) then
               temp(it2min:istop) = temp(it2min:istop) + t2(it2min:istop)*PVec(it2min:istop)
            else
               temp(it2min:it1max) = temp(it2min:it1max) + t2(it2min:it1max)*PVec(it2min:it1max)
               temp(it1max:istop) =  t2(it1max:istop)*PVec(it1max:istop)
            end if
         else
            if(it1max .ge. istop) then
               temp(if1:istop) = temp(if1:istop) + t2(if1:istop)*PVec(if1:istop)
               temp(it2min:if1) =  t2(it2min:if1)*PVec(it2min:if1)
            else
               temp(if1:it1max) = temp(if1:it1max) + t2(if1:it1max)*PVec(if1:it1max)
               temp(it2min:if1) =  t2(it2min:if1)*PVec(it2min:if1)
               temp(it1max:istop) =  t2(it1max:istop)*PVec(it1max:istop)
            end if
         endif
! limits: 
         i1 = min(if1,it2min)
         i2 = max(it1max,istop)

         i1a = i1
         i2a = i2
         call minmaxi_form(temp,grid%nr,i1a,i2a)
         
         i1 = max(i1,i1a)
         i2 = min(i2,i2a)

         deallocate(PVec, QVec)

end subroutine form_spheroidal

subroutine form_spheroidal_accurate(l,m,fun,minfun,maxfun,maxfm,temp,i1,i2)
  !LHS: added accurate form subroutine - following Igor's cc/form.f
  use grid_radial 
  implicit none
  
  integer, intent(in):: l, m
  real*8, dimension(grid%nr), intent(in):: fun
  integer, intent(in):: minfun, maxfun, maxfm
  real*8, dimension(grid%nr), intent(out):: temp
  integer, intent(out):: i1, i2
  
  real*8, dimension(0:grid%nr):: t1, t2
  real*8, dimension(:), pointer :: PVec, QVec
  integer:: i, n
  integer:: if1, if2, iP1, iQ2, istart, istop
  real*8, dimension(:), pointer :: dr_on_three, dr_three_on_eight
  real*8, dimension(grid%nr) :: Pfun, Qfun 

  !radial grid spacings
  dr_on_three => grid%dr_on_three
  dr_three_on_eight => grid%dr_three_on_eight

  if(l .gt. grid%ltmax) then
     i1=2
     i2=1     
     return
  endif

  ! Translate the 2D (l,m) pair to the 1D (n=sum{l}+m) index.
  n = abs(m)
  do i = 1, l
     n = n + i
  end do

  !iP1, iQ2 limits of Legendre P,Q functions
  iP1 = grid%i1_P(n)
  iQ2 = grid%i2_Q(n)
  !minfun, maxfun limits of fun
  !if1, if2 will be range of integration
  if1 = max(minfun,iP1)
  if2 = min(maxfun,iQ2)
  if2 = min(if2, maxfm)

  PVec => grid%LegendreP(:,n)
  QVec => grid%LegendreQ(:,n)

  !temp = 0.0d0
  !t1=0.0d0
  !t2=0.0d0

  !Pfun = P(r)*fun(r) etc
  Pfun(if1:if2) = PVec(if1:if2)*fun(if1:if2)
  Qfun(if1:if2) = QVec(if1:if2)*fun(if1:if2)
 
  !integration indices:
    !integral is required over if1 to if2 based on radial limits of Pfun and Qfun
    !we want to iterate only over even indices so we set
    istart=(if1+1)/2*2! increase if1 by 1 if it is odd
    istop=if2/2*2     ! decrease if2 by 1 if it is odd
  
  if (istop - istart.lt.4) then
    i1=2
    i2=1
    return
  end if

  !The iteration below starts at istart+2 since each step refers to the two
  !below it.
  !If if1 is odd then the iteration starts at istart+2 = if1+3
  !We need to set t1(if1), and also t1(if1+1) and t1(if1+2) for the iteration.
  !If if1 is odd we only need to set t1(if1) and t1(if1+1) so 
  !We assume Pfun(if1-1) = 0

  if(istart==if1) Pfun(istart-1) = 0.0d0

  !Integrate fun(r')*P(r') from 0 to r -> store in t1
  t1(istart)=(4.0d0*Pfun(istart-1)+Pfun(istart))*dr_on_three(istart) !simpson with Pfun(if1-1) = 0
  do i=istart+2, istop-2,2
    !even points: Simpson 1/3 rule:
      t1(i)=t1(i-2)+(Pfun(i-2)+4.0d0*Pfun(i-1)+Pfun(i))*dr_on_three(i)
    !odd points: Simpson 3/8 rule:
      t1(i+1)=t1(i-2)+(Pfun(i-2)+3.0d0*(Pfun(i-1)+Pfun(i))+Pfun(i+1))*dr_three_on_eight(i+1)
  enddo
  
  if(istart==if1+1) then !if1 is odd 
    t1(if1)=t1(if1+1)*0.75d0 - t1(if1+3)/8.0d0
    t1(if1+2)=t1(if1+1)*0.75d0 + t1(if1+3)*0.375d0 !t1(if1-1)=0
  else
    !if1 is even, Pfun(istart-1) is zero
    t1(if1+1)=t1(if1)*0.75d0 + t1(if1+2)*0.375d0 !t1(if1-1)=0
  endif
  
  !set last point, t1(istop-1) was set in last iteration above
  t1(istop)=t1(istop-2)+(Pfun(istop-2)+4.0d0*Pfun(istop-1)+Pfun(istop))*dr_on_three(istop)
  if(if2==istop+1) then !if2 is odd so t1(if2) has not been set
    t1(if2)=t1(if2-3)+(Pfun(if2-3)+3.0d0*Pfun(if2-2)+3.0d0*Pfun(if2-1)+Pfun(if2))*dr_three_on_eight(if2) !newton 3/8
  endif

  !above does not work when grid spacing has doubled.
  !1/3 rule works at doubling point since it is always at an even index.
  do n=2, grid%ndouble+1 !cycle through each grid doubling
    i=grid%jdouble(n)+1  !index of nth doubling
    if(i<istart.or.i>istop) cycle
    if(i>if2-2) then !can't refer to t1(i+2)
      t1(i)=t1(i-1)
      cycle
    endif
    !3/8 rule breaks down at the first index after doubling
    t1(i)=t1(i+2)-(Pfun(i)+4.0d0*Pfun(i+1)+Pfun(i+2))*dr_on_three(i+2)
  enddo

  !Iterate backwards for integral from r to infinity
  if(istop .ne. if2) t2(if2) = 0.0d0
  t2(istop) = 0.0d0
  do i=istop,istart+2,-2
    t2(i-2)=t2(i)+(Qfun(i-2) + 4.0d0*Qfun(i-1) +Qfun(i))*dr_on_three(i)
    t2(i-3)=t2(i)+(Qfun(i-3) + 3.0d0*(Qfun(i-2)+Qfun(i-1)) + Qfun(i))*dr_three_on_eight(i)
  enddo
  
  !fix issues at doubling points again
  do n=2, grid%ndouble+1
    i=grid%jdouble(n)-1
    if(i<istart) cycle
    t2(i)=t2(i-2)-(Qfun(i-2)+4.0d0*Qfun(i-1)+Qfun(i))*dr_on_three(i)
    if(i==istart+1) t2(i)=t2(istart)
  enddo
  t2(istop-1)=t2(istop-3)-(Qfun(istop-3)+4.0d0*Qfun(istop-2)+Qfun(istop-1))*dr_on_three(istop-1)

  temp(if1:if2) = t1(if1:if2)*QVec(if1:if2) + t2(if1:if2)*PVec(if1:if2) 
  
  i=if2
  do while (abs(temp(i)) > grid%formcut .and. i<min(maxfm,iQ2))
    i=i+1
    temp(i)=t1(if2)*QVec(i)
  enddo
  i2=i
  do while(abs(temp(i2)) < grid%formcut .and. i2 > if1)
    i2=i2-1
  enddo
  i=if1
  do while(abs(temp(i))>grid%formcut .and. i > iP1)
    i=i-1
    temp(i)=t2(if1)*PVec(i)
  enddo
  i1=i
  do while(abs(temp(i1)) < grid%formcut .and. i1 <= i2)
    i1=i1+1
  enddo
  
  !deallocate(PVec, QVec)

end subroutine form_spheroidal_accurate

subroutine form_spheroidal_accurate2(l,m,Pfun,Qfun,minfun,maxfun,maxfm,temp,i1,i2,load)
  !LHS: added accurate form subroutine - following Igor's cc/form.f
  use grid_radial 
  implicit none
  
  integer, intent(in):: l, m
  real*8, dimension(grid%nr), intent(inout):: Pfun, Qfun
  integer, intent(in):: minfun, maxfun, maxfm
  real*8, dimension(grid%nr), intent(out):: temp
  integer, intent(out):: i1, i2
  logical, intent(in) :: load
  
  real*8, dimension(0:grid%nr):: t1, t2
  real*8, dimension(:), pointer :: PVec, QVec
  integer:: i, n
  integer:: if1, if2, iP1, iQ2, iP2, iQ1, istart, istop
  real*8, dimension(:), pointer :: dr_on_three, dr_three_on_eight
  
  !radial grid spacings
  dr_on_three => grid%dr_on_three
  dr_three_on_eight => grid%dr_three_on_eight

  if(l .gt. grid%ltmax) then
     i1=2
     i2=1     
     return
  endif

  ! Translate the 2D (l,m) pair to the 1D (n=sum{l}+m) index.
  n = abs(m)
  do i = 1, l
     n = n + i
  end do

  !iP1, iQ2 limits of Legendre P,Q functions
  iP1 = grid%i1_P(n)
  !iP2 = grid%i2_P(n)
  iQ2 = grid%i2_Q(n)
  !iQ1 = grid%i1_Q(n)
  !minfun, maxfun limits of fun
  !if1, if2 will be range of integration
  if1 = max(minfun,iP1)
  if2 = min(maxfun,iQ2)
  if2 = min(if2, maxfm)

  if(load) then
    i1 = if1
    i2 = if2
    return
  endif

  PVec => grid%LegendreP(:,n)
  QVec => grid%LegendreQ(:,n)

  !temp = 0.0d0
  !t1=0.0d0
  !t2=0.0d0

  !Pfun = P(r)*fun(r) etc
  !Pfun(if1:if2) = PVec(if1:if2)*fun(if1:if2)
  !Qfun(if1:if2) = QVec(if1:if2)*fun(if1:if2)
 
  !integration indices:
    !integral is required over if1 to if2 based on radial limits of Pfun and Qfun
    !we want to iterate only over even indices so we set
    istart=(if1+1)/2*2! increase if1 by 1 if it is odd
    istop=if2/2*2     ! decrease if2 by 1 if it is odd
  
  if (istop - istart.lt.4) then
    i1=2
    i2=1
    return
  end if

  !The iteration below starts at istart+2 since each step refers to the two
  !below it.
  !If if1 is odd then the iteration starts at istart+2 = if1+3
  !We need to set t1(if1), and also t1(if1+1) and t1(if1+2) for the iteration.
  !If if1 is odd we only need to set t1(if1) and t1(if1+1) so 
  !We assume Pfun(if1-1) = 0

  if(istart==if1) Pfun(istart-1) = 0.0d0

  !Integrate fun(r')*P(r') from 0 to r -> store in t1
  t1(istart)=(4.0d0*Pfun(istart-1)+Pfun(istart))*dr_on_three(istart) !simpson with Pfun(if1-1) = 0
  do i=istart+2, istop-2,2
    !even points: Simpson 1/3 rule:
      t1(i) = (Pfun(i-2)+4.0d0*Pfun(i-1)+Pfun(i))*dr_on_three(i)
    !odd points: Simpson 3/8 rule:
      t1(i+1) = (Pfun(i-2)+3.0d0*(Pfun(i-1)+Pfun(i))+Pfun(i+1))*dr_three_on_eight(i+1)
  enddo
  do i=istart+2, istop-2,2
    t1(i) = t1(i-2) + t1(i)
    t1(i+1) = t1(i-2) + t1(i+1) 
  enddo
  
  if(istart==if1+1) then !if1 is odd 
    t1(if1)=t1(if1+1)*0.75d0 - t1(if1+3)/8.0d0
    t1(if1+2)=t1(if1+1)*0.75d0 + t1(if1+3)*0.375d0 !t1(if1-1)=0
  else
    !if1 is even, Pfun(istart-1) is zero
    t1(if1+1)=t1(if1)*0.75d0 + t1(if1+2)*0.375d0 !t1(if1-1)=0
  endif
  
  !set last point, t1(istop-1) was set in last iteration above
  t1(istop)=t1(istop-2)+(Pfun(istop-2)+4.0d0*Pfun(istop-1)+Pfun(istop))*dr_on_three(istop)
  if(if2==istop+1) then !if2 is odd so t1(if2) has not been set
    t1(if2)=t1(if2-3)+(Pfun(if2-3)+3.0d0*Pfun(if2-2)+3.0d0*Pfun(if2-1)+Pfun(if2))*dr_three_on_eight(if2) !newton 3/8
  endif

  !above does not work when grid spacing has doubled.
  !1/3 rule works at doubling point since it is always at an even index.
  do n=2, grid%ndouble+1 !cycle through each grid doubling
    i=grid%jdouble(n)+1  !index of nth doubling
    if(i<istart.or.i>istop) cycle
    if(i>if2-2) then !can't refer to t1(i+2)
      t1(i)=t1(i-1)
      cycle
    endif
    !3/8 rule breaks down at the first index after doubling
    t1(i)=t1(i+2)-(Pfun(i)+4.0d0*Pfun(i+1)+Pfun(i+2))*dr_on_three(i+2)
  enddo

  !Iterate backwards for integral from r to infinity
  if(istop .ne. if2) t2(if2) = 0.0d0
  t2(istop) = 0.0d0
  do i=istop-2,istart,-2
    t2(i) = (Qfun(i) + 4.0d0*Qfun(i+1) +Qfun(i+2))*dr_on_three(i+2)
    t2(i-1) = (Qfun(i-1) + 3.0d0*(Qfun(i)+Qfun(i+1)) + Qfun(i+2))*dr_three_on_eight(i+2)
  enddo
  do i=istop-2,istart,-2
    t2(i)   = t2(i) + t2(i+2)
    t2(i-1) = t2(i-1) + t2(i+2)
  enddo
  
  !fix issues at doubling points again
  do n=2, grid%ndouble+1
    i=grid%jdouble(n)-1
    if(i<istart) cycle
    t2(i)=t2(i-2)-(Qfun(i-2)+4.0d0*Qfun(i-1)+Qfun(i))*dr_on_three(i)
    if(i==istart+1) t2(i)=t2(istart)
  enddo
  t2(istop-1)=t2(istop-3)-(Qfun(istop-3)+4.0d0*Qfun(istop-2)+Qfun(istop-1))*dr_on_three(istop-1)

  temp(if1:if2) = t1(if1:if2)*QVec(if1:if2) + t2(if1:if2)*PVec(if1:if2) 

!  temp(iQ1:iQ2) = t1(iQ1:iQ2)*Qvec(iQ1:iQ2)
!  temp(iP1:iQ1-1) = 0.0d0; temp(iQ2+1:iP2) = 0.0d0
!  temp(iP1:iP2) = temp(iP1:iP2) + t2(iP1:iP2)*Pvec(iP1:iP2)
  
  i=if2
  do while (abs(t1(if2)*Qvec(i)) > grid%formcut .and. i<min(maxfm,iQ2))
    i=i+1
    temp(i)=t1(if2)*QVec(i)
  enddo
  i2=i
  do while(abs(temp(i2)) < grid%formcut .and. i2 > if1)
    i2=i2-1
  enddo
  i=if1
  do while(abs(temp(i))>grid%formcut .and. i > iP1)
    i=i-1
    temp(i)=t2(if1)*PVec(i)
  enddo
  i1=i
  do while(abs(temp(i1)) < grid%formcut .and. i1 <= i2)
    i1=i1+1
  enddo
  
  !deallocate(PVec, QVec)

end subroutine form_spheroidal_accurate2
