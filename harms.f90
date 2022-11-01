!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Module: harms.f90
!Purpose: contains functions for evaluating the spherical harmonics
!         Bizarrely, functions to do this do not come standard with 
!         any fortran compilers (Come in people, you've had since 1977!)
!Author: Reese Horton, Curtin University
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module harms


   contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !function: Ylm
      !Purpose: returns the value of the spherical harmonic with given 
      !         arguments and l,m values as a double precision real
      !Input: theta and phi in radians, l and m as integers
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      complex(dp) function Ylm(theta, phi, l, m) result(harm)
	 use numbers
	 implicit none
	 real(dp):: ctheta !cosine of theta
         real(dp):: theta, phi
	 integer:: l, m
	 integer(kind=largeInt):: denom

	 ctheta = cos(theta)
	 denom = getDenom(l,m)
	 coeff = sqrt(      /dble(denom))
	 harm = coeff*legendre(ctheta,l,m)*exp(i*dble(m)*phi)

      end function Ylm



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !function: legendre
      !Purpose: computes the value of the associated legendre polynomial
      !         with the given argument
      !Input: value u from -1, 1 with integers (l,m)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(dp) function legendre(u,l,m) result(leg)
         use numbers
	 real(dp):: u  !Argument lying on the unit interval
	 integer:: l, m   !Angular momentum quantum numbers

	 if (abs(u) .gt. 1) then
	    print*, "ERROR: argument to legendre function has magnitude greater than 1. Stopping."
	    stop
	 end if

      end function legendre


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Function: getCoeff
      !Purpose: computes the normalisation coefficient of the spherical 
      !         harmonics. The denominator in the normalisation coefficients 
      !         of the spherical harmonics, needs to be stored as a
      !         very large integer due to the use of factorials in its
      !         definition.
      !Input: requires l .ge. 0 and m .ge. 0
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(dp) function getCoeff(l,m) result(coeff)
	 use numbers
	 integer:: l,m
	 integer:: ii
	 integer(kind=largeInt):: denom

	 if (m .lt. 0) then
	    print*, "ERROR: normalisation coefficients for spherical harmonics defined only for positive m. Stopping"
	    stop
	 end if

	 denom = 1
	 if (m .gt. 0) then
	    do ii = l-m+1, l+m
	       denom = denom*ii
	    end do
	 end if

	 coeff = ((-1.0dp)**m)*sqrt(dble(2*l+1)/4*pi)*sqrt(1.0_dp/dble(denom)) 

      end function

end module harms

