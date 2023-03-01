!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Module: realHarms
!Purpose: contains functions for the evaluation of real spherical 
!         harmonics and their coupling coefficients.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module realHarms


	 contains

	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: XLM
	 !Purpose: evaluates the real spherical harmonics for theta and phi in radians
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(dpf) function XLM(l, m, theta, phi) result(res)
	    use numbers
	    implicit none
	    integer:: l, m
	    real(dpf):: theta, phi
	    real(dpf):: RYLM

	    res = 0.0_dpf
	    if (m .gt. 0) then
	       res = RYLM(l,m,theta)*sqrt(2.0_dpf)*cos(dble(m)*phi)
	    else if (m .lt. 0) then
	       res = RYLM(l,abs(m),theta)*sqrt(2.0_dpf)*sin(dble(abs(m))*phi)
	    else if (m .eq. 0) then
	       !Special case, m=0 equivalent to theta=0, Ylm and Xlm coincide in this case.
	       res = RYLM(l,m,theta)
	    end if
	 end function XLM


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: Xint
	 !Purpose: calculates an integral of the product of three real spherical harmonics
	 !         Xlm(theta,phi). Also known as real gaunt coefficients.
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 real(dpf) function Xint(l1, mu1, l2, mu2, l3, mu3) result(res)
      use ieee_arithmetic
	    use numbers
	    implicit none
	    real(dpf):: l1, mu1, l2, mu2, l3, mu3 
	    real(dpf):: Yint
	    real(dpf):: val
	    real(dpf):: t1, t2
	    real(dpf), dimension(3):: muvals, lVals
	    integer:: ii, numZero
	    integer, dimension(3):: isZero
	    real(dpf):: templ, tempmu
	    logical:: found
	    !complex(dpf):: v1, v2

	    !val = sqrt(3.0_dpf/(4.0_dpf*pi))*Yint(1.0_dpf,1.0_dpf,1.0_dpf,1.0_dpf,0.0_dpf,0.0_dpf)
	    !print*, val
	    !stop

	    !Use complete permutation symmetry of real gaunt coefficients to
	    !reduce calculation to one of several cases. 
	    !Case 1: mu1, mu2, mu3 all != 0
	    !Case 2: mu3 = 0, mu2, mu1 != 0
	    !Case 3: mu2 = mu3 = 0, mu1 != 0
	    !Simpler formulas apply in each of these cases
	    if (mod(int(l1+l2+l3),2) .ne. 0) then
	       res = 0.0_dpf
	    else
         muVals(1) = mu1
         muVals(2) = mu2
         muVals(3) = mu3
         lVals(1) = l1
         lVals(2) = l2
         lVals(3) = l3
               
         numZero = 0
         isZero(:) = 0
         do ii = 1, 3
            if (int(abs(muVals(ii))) .eq. 0) then
               isZero(ii) = 1
            end if
         end do
	       numZero = sum(isZero)

	       if (numZero .eq. 1) then
		        found = .false.
		        ii = 1 
		        do while (isZero(ii) .ne. 1)
		           ii = ii + 1
		        end do
            if (ii .ne. 3) then
               templ = l3
               tempmu = mu3
               l3 = lVals(ii)
               mu3 = muVals(ii)
               lVals(ii) = templ
               muVals(ii) = tempmu
               l1 = lVals(1)
               mu1 = muVals(1)
               l2 = lVals(2)
               mu2 = muVals(2)
		        end if

	          res = 2.0_dpf*sqrt((2.0_dpf*l2+1.0_dpf)/(4.0_dpf*pinum))*Yint(l1,mu2,l2,mu2,l3,0.0_dpf) &
	                *REAL(conjg(VCoeff(mu2,mu1))*VCoeff(mu2,mu2))
	       else if (numZero .eq. 2) then
	          res = 0.0_dpf
	       else if (numZero .eq. 3) then
                  !found = .false.
                  !ii = 0
                  !do while (.not. found)
                  !   ii = ii + 1
                  !   if (int(abs(muVals(ii))) .ne. 0) then
                  !      found = .true.
                  !   end if
                  !end do
                  !if (ii .ne. 1) then
                  !   templ = lVals(ii)
                  !   tempmu = muVals(ii)
                  !   lVals(ii) = l1
                  !   muVals(ii) = mu1
                  !   l1 = templ
                  !   mu1 = tempmu
                  !   l2 = lVals(2)
                  !   mu2 = muVals(2)
                  !   l3 = lVals(3)
                  !   mu3 = muVals(3)
                  !end if

            
	          !Xl0 and Yl0 coincide, special case. Allows use of existing function for Ylm.
	          res =  sqrt((2.0_dpf*l2+1.0_dpf)/(4.0_dpf*pinum))*Yint(l1,0.0_dpf,l2,0.0_dpf,l3,0.0_dpf)
	       else if (numZero .eq. 0) then
		        !Formula for overlap of Xlm in terms of overlap of Ylm
            t1 = 2.0_dpf*sqrt((2.0_dpf*l2+1.0_dpf)/(4.0_dpf*pinum))*Yint(l1,mu2+mu3,l2,mu2,l3,mu3) &
		             *REAL(conjg(VCoeff(mu2+mu3,mu1))*VCoeff(mu2,mu2)*VCoeff(mu3,mu3))
            t2 = 2.0_dpf*sqrt((2.0_dpf*l2+1.0_dpf)/(4.0_dpf*pinum))*Yint(l1,mu2-mu3,l2,mu2,l3,-mu3) & 
                 *REAL(conjg(VCoeff(mu2-mu3,mu1))*VCoeff(mu2,mu2)*VCoeff(-mu3,mu3))
                  !if ((abs(t1) .gt. 0.0_dpf) .and. (abs(t2) .gt. 0.0_dpf)) then
                      	!	 print*, "ERROR", t1, t2
                      	!	 v1 = VCoeff(mu2+mu3,mu1)
                      	!	 v2 = VCoeff(mu2-mu3,mu1)
                  !  print*, v1, abs(mu1), abs(mu2+mu3)
                        !	 print*, v2, abs(mu1), abs(mu2-mu3)
                        !	 print*, "DONE"
                        
                        
                        !	 !print*, "ERROR: ", l1, mu1, l2, mu2, l3, mu3
                        !	 !print*, "Vals: ", mu1,  mu2+mu3, -mu2-mu3, mu2-mu3,mu3-mu2
                        !end if
                        
                  res = t1 + t2
		        !res = 2.0_dpf*sqrt((2.0_dpf*l2+1.0_dpf)/(4.0_dpf*pi))*Yint(l1,mu2+mu3,l2,mu2,l3,mu3) &
		        !      *REAL(conjg(VCoeff(mu2+mu3,mu1))*VCoeff(mu2,mu2) &
		        !      *VCoeff(mu3,mu3)) + 2.0_dpf*sqrt((2.0_dpf*l2+1.0_dpf)/(4.0_dpf*pi)) &
			      !*Yint(l1,mu2-mu3,l2,mu2,l3,-mu3)*REAL(conjg(VCoeff(mu2-mu3,mu1)) &
			      !*VCoeff(mu2,mu2)*VCoeff(-mu3,mu3))
         else
            print*, "TREAT ALL CASES!!!"
            stop
	       end if
	    end if
	 end function Xint



	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: VCoeff
	 !Purpose: evaluates the coefficients in the unitary transformation from 
	 !         complex to real spherical harmonics.
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 complex(dpf) function VCoeff(m,mu) result(res)
	    implicit none
	    real(dpf):: m,mu
	    complex(dpf):: i

	    !if (abs(m) .gt. l) then
	    !   print*, l, m
	    !   print*, "ERROR: angular momentum projection greater than magnitude, stopping."
	    !   stop
	    !end if

	    i = (0.0_dpf,1.0_dpf)

	    res = kronecker(int(m),0)*kronecker(int(mu),0) &
			      + (1.0_dpf/sqrt(2.0_dpf))*(heavyside(mu)*kronecker(int(m),int(mu)) &
	          + heavyside(-mu)*(i)*((-1.0_dpf)**m)*kronecker(int(m),int(mu)) &
				    + heavyside(-mu)*(-i)*kronecker(int(m),int(-mu)) &
		        + heavyside(mu)*((-1.0_dpf)**m)*kronecker(int(m),int(-mu)))
	 end function VCoeff



	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: heavyside
	 !Purpose: evaluates the heavyside function of the argument
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 real(dpf) function heavyside(m) result(res)
	    implicit none
	    real(dpf):: m

	    if (m .gt. 0.0_dpf) then
	       res = 1.0_dpf
	    else 
	       res = 0.0_dpf
	    end if
	 end function heavyside


	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 !Function: kronecker
	 !Purpose: evaluates the kronecker delta as a real number
	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	 real(dpf) function kronecker(n,m) result(res)
	    implicit none
	    integer:: n,m

	    res = 0.0_dpf
	    if (n .eq. m) then
	       res = 1.0_dpf
	    end if
	 end function kronecker

end module realHarms
