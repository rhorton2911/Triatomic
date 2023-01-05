MODULE accuracy
! Define floating-point working precision, idp
  IMPLICIT NONE
  INTEGER, PARAMETER        :: idp = SELECTED_REAL_KIND(15, 307)
END MODULE  accuracy

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

MODULE input_output
!***begin prologue     input_output
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           
!***author             schneider, b. i.(nsf)
!***source             
!***purpose            
!
!***references

!***routines called    
!***end prologue       input_output

  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                   Some Constants and Units
!
   INTEGER                                       :: inp=8
   INTEGER                                       :: iout=9
   INTEGER                                       :: rows_to_print
   INTEGER                                       :: columns_to_print
   INTEGER                                       :: eigenvectors_to_print
   LOGICAL                                       :: print_parameter
   CHARACTER(LEN=8), DIMENSION(:), ALLOCATABLE   :: rowlab, collab

!
!
END MODULE input_output

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

MODULE Data_Module
!***begin prologue     Data_Module
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)
!***keywords           time, finite element dvr, orthogonal polynomial
!***author             schneider, b. i.(nsf)
!***source             dvrlib
!***purpose            global shared variables for dvr library
!***description        this routine defines the global variables
!***                   and data needed for the dvrlib
!
!***references

!***routines called    
!***end prologue       Data_Module
  use numbers
  USE accuracy
  USE input_output
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!---------------------------------------------------------------------
!                   Some Constants and Units
!---------------------------------------------------------------------
!
   REAL(idp)                              ::          pi    = 3.141592653589793238462643383276D0 
   REAL(idp)                              ::          two_pi= 6.283185307179586476925286766552D0 
   REAL(idp)                              ::          zero    =  0.D0
   REAL(idp)                              ::          quarter = .25D0
   REAL(idp)                              ::          half    = .50D0
   REAL(idp)                              ::          third   = .33333333333333333333333333333D0
   REAL(idp)                              ::          fourth  = .25000000000000000000000000000D0
   REAL(idp)                              ::          fifth   = .20000000000000000000000000000D0
   REAL(idp)                              ::          sixth   = .16666666666666666666666666666D0
   REAL(idp)                              ::          seventh = .14285714285714285714285714285D0
   REAL(idp)                              ::          eighth  = .12500000000000000000000000000D0
   REAL(idp)                              ::          ninth   = .11111111111111111111111111111D0
   REAL(idp)                              ::          tenth   = .10000000000000000000000000000D0
   REAL(idp)                              ::          one     = 1.0D0
   REAL(idp)                              ::          two     = 2.0D0
   REAL(idp)                              ::          three   = 3.0D0
   REAL(idp)                              ::          four    = 4.0D0
   REAL(idp)                              ::          five    = 5.0D0
   REAL(idp)                              ::          six     = 6.0D0
   REAL(idp)                              ::          seven   = 7.0D0
   REAL(idp)                              ::          eight   = 8.0D0
   REAL(idp)                              ::          nine    = 9.0D0
   REAL(idp)                              ::          ten     = 10.D0
   REAL(idp)                              ::          nrzero  = 1.D-06
   INTEGER                                ::          int_zero      = 0
   INTEGER                                ::          int_one       = 1
   INTEGER                                ::          int_two       = 2
   INTEGER                                ::          int_three     = 3
   INTEGER                                ::          int_four      = 4
   INTEGER                                ::          int_five      = 5
   INTEGER                                ::          int_six       = 6
   INTEGER                                ::          int_seven     = 7
   INTEGER                                ::          int_eight     = 8
   INTEGER                                ::          int_nine      = 9
   INTEGER                                ::          int_ten       = 10
   INTEGER                                ::          int_eleven    = 11
   INTEGER                                ::          int_twelve    = 12
   INTEGER                                ::          int_thirteen  = 13
   INTEGER                                ::          int_fourteen  = 14
   INTEGER                                ::          int_fifteen   = 15
   INTEGER                                ::          int_sixteen   = 16
   INTEGER                                ::          int_seventeen = 17
   INTEGER                                ::          int_eighteen  = 18
   INTEGER                                ::          int_nineteen  = 19
   INTEGER                                ::          int_twenty    = 20
!
!                    hbar in joule-sec
!
   REAL(idp)                              ::          hbar = 1.054571596D-34
!
!                    electron mass in kg
!
   REAL(idp)                              ::          massau = 9.10938188D-31
!
!                    bohr radius in meters
!
   REAL(idp)                              ::          lenau = 5.291772083D-11
!
!                    time for an electron to make one bohr orbit in seconds
!
   REAL(idp)                              ::          timau    = 2.418884326D-17
   REAL(idp)                              ::          efieldau = 5.14220624D+11
   REAL(idp)                              ::          electric_field_to_intensity = 3.509338D+16
   REAL(idp)                              ::          peak_electric_field = .2849540283D-03
   REAL(idp)                              ::          pmass    = 1.67262158D-27
   REAL(idp)                              ::          massn2p  = 1.00137841887D0
   REAL(idp)                              ::          au_in_ev = 27.211396132D0
!
!
END MODULE Data_Module

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

!***********************************************************************
                           MODULE Matrix_Print
!
                           use numbers
													 use accuracy

                           IMPLICIT NONE
!
!
                         INTERFACE Print_Matrix
                   MODULE PROCEDURE Print_Matrix_d,                                    &
                                    Print_Matrix_z,                                    &
                                    Print_Triangle_Matrix_d,                           &
                                    Print_Triangle_Matrix_z,                           &          
                                    Print_Vector_d,                                    &
                                    Print_Vector_z                                     
                          END INTERFACE Print_Matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           Contains
!***********************************************************************
!***********************************************************************
!deck Print_Matrix_d
  subroutine Print_Matrix_d(a,n,m,iout,frmt,title,collab,rowlab)
  implicit none
  REAL(idp), DIMENSION(:,:)                   :: a
  INTEGER                                  :: n
  INTEGER                                  :: m
  INTEGER                                  :: iout
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: rowlab
  INTEGER                                  :: ibeg
  INTEGER                                  :: iend
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: k
  INTEGER                                  :: ii
  IF (PRESENT(title)) THEN
      !write(iout,1) title
  END IF
  ibeg=0
  IF (PRESENT(frmt)) THEN
      do i=1,m,5
         iend=min(ibeg+5,m)
         ibeg=ibeg+1
         IF (PRESENT(collab)) THEN
             !write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         IF ( PRESENT(rowlab) ) THEN
              IF (frmt == 'f') then
                  do j=1,n
                     write (iout,4) rowlab(j),  (a(j,k),k=ibeg,iend)
                 end do
              else
                  do j=1,n
	             write (iout,5) rowlab(j),  (a(j,k),k=ibeg,iend)
                  end do
              end if
         ELSE
              IF (frmt == 'f') then
                  do j=1,n
                     write (iout,6) (a(j,k),k=ibeg,iend)
                  end do
              else
                  do j=1,n
	             write (iout,7)  (a(j,k),k=ibeg,iend)
                  end do
              end if
         END IF
         ibeg=iend
      end do
  ELSE
      do i=1,m,5
         iend=min(ibeg+5,m)
         ibeg=ibeg+1
         IF (PRESENT(collab)) THEN
             !write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         IF ( PRESENT(rowlab) ) THEN
              do j=1,n
                 write (iout,5) rowlab(j),  (a(j,k),k=ibeg,iend)
              end do
         ELSE
              do j=1,n
                 write (iout,7) (a(j,k),k=ibeg,iend)
              end do
         END IF
         ibeg=iend
      end do
  END IF
1 format(/,a80,/)
2 format(15x,5(3x,a8,9x))
3 format(1x,'col',5(5x,i6,9x))
4 format(4x,a6,5(1pf20.12))
5 format(4x,a6,5(1pe20.12))
6 format(4x,5(1pf20.12))
7 format(4x,5(1pe20.12))
end subroutine  Print_Matrix_d
!***********************************************************************
!***********************************************************************
!deck Print_Matrix_z
  subroutine Print_Matrix_z(a,n,m,iout,frmt,title,collab,rowlab)
  implicit none
  COMPLEX(idp), DIMENSION(:,:)               :: a
  INTEGER                                  :: n
  INTEGER                                  :: m
  INTEGER                                  :: iout
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: rowlab
  INTEGER                                  :: ibeg
  INTEGER                                  :: iend
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: k
  INTEGER                                  :: ii
  IF (PRESENT(title)) THEN
      !write(iout,1) title
  END IF
  ibeg=0
  IF (PRESENT(frmt)) THEN
      do i=1,m,2
         iend=min(ibeg+2,m)
         ibeg=ibeg+1
         IF (PRESENT(collab)) THEN
             !write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         IF ( PRESENT(rowlab) ) THEN
              IF (frmt == 'f') then
                  do j=1,n
                     write (iout,4) rowlab(j),  (a(j,k),k=ibeg,iend)
                 end do
              else
                  do j=1,n
	             write (iout,5) rowlab(j),  (a(j,k),k=ibeg,iend)
                  end do
              end if
         ELSE
              IF (frmt == 'f') then
                  do j=1,n
                     write (iout,6) (a(j,k),k=ibeg,iend)
                  end do
              else
                  do j=1,n
	             write (iout,7)  (a(j,k),k=ibeg,iend)
                  end do
              end if
         END IF
         ibeg=iend
      end do
  ELSE
      do i=1,m,2
         iend=min(ibeg+2,m)
         ibeg=ibeg+1
         IF (PRESENT(collab)) THEN
             !write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         IF ( PRESENT(rowlab) ) THEN
              do j=1,n
                 write (iout,6) rowlab(j),  (a(j,k),k=ibeg,iend)
              end do
         ELSE
              do j=1,n
                 write (iout,7) (a(j,k),k=ibeg,iend)
              end do
         END IF
         ibeg=iend
      end do
  END IF
1 format(/,a80,/)
2 format(15x,2(15x,a8,18x))
3 format(1x,'col(Re,Im)',4(9x,i6,25x))
4 format(4x,a6,4(1pf20.12))
5 format(4x,a6,4(1pe20.12))
6 format(4x,4(1pf20.12))
7 format(4x,4(1pe20.12))
end subroutine  Print_Matrix_z
!***********************************************************************
!***********************************************************************
!deck Print_Triangle_Matrix_d
  subroutine Print_Triangle_Matrix_d(a,n,iout,frmt,title,collab)
  implicit none
  REAL(idp), DIMENSION (:)                    :: a
  INTEGER                                  :: n
  INTEGER                                  :: iout
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  INTEGER                                  :: ibeg
  INTEGER                                  :: iend
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: ii
  IF (PRESENT(title)) THEN
      !write(iout,1) title
  END IF
  ibeg = 0
  count = 0
  IF (PRESENT(frmt)) THEN
      DO i = 1, n, 5
         iend = min(ibeg+5,n)
         ibeg = ibeg + 1
         IF (PRESENT(collab)) THEN
             !write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         IF (frmt == 'f') THEN
             DO j = ibeg, iend
                write (iout,4) a(count + 1:count + j)
                count=count+j
             END DO
         ELSE
             DO j = ibeg, iend
                write (iout,5) a(count + 1:count + j)
                count=count+j
             END DO
         END IF
         ibeg = iend
      END DO
  ELSE
      DO i = 1, n, 5
         iend = min(ibeg+5,n)
         ibeg = ibeg + 1
         IF (PRESENT(collab)) THEN
             !write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         DO j = ibeg, iend
            write (iout,5) a(count + 1:count + j)
            count=count+j
         END DO
         ibeg = iend
      END DO
  END IF
1 format(/,a80,/)
2 format(15x,5(3x,a8,9x))
3 format(1x,'col',5(5x,i6,9x))
4 format(10x,5(1pf20.12))
5 format(10x,5(1pe20.12))
end subroutine  Print_Triangle_Matrix_d
!***********************************************************************
!***********************************************************************
!deck Print_Triangle_Matrix_z
  subroutine Print_Triangle_Matrix_z(a,n,iout,frmt,title,collab)
  implicit none
  COMPLEX(idp), DIMENSION (:)                :: a
  INTEGER                                  :: n
  INTEGER                                  :: iout
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  INTEGER                                  :: ibeg
  INTEGER                                  :: iend
  INTEGER                                  :: count
  INTEGER                                  :: i
  INTEGER                                  :: j
  INTEGER                                  :: ii
  IF (PRESENT(title)) THEN
      !write(iout,1) title
  END IF
  ibeg = 0
  count = 0
  IF (PRESENT(frmt)) THEN
      DO i = 1, n, 5
         iend = min(ibeg+5,n)
         ibeg = ibeg + 1
         IF (PRESENT(collab)) THEN
             !write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         IF (frmt == 'f') THEN
             DO j = ibeg, iend
                write (iout,4) a(count + 1:count + j)
                count=count+j
             END DO
         ELSE
             DO j = ibeg, iend
                write (iout,5) a(count + 1:count + j)
                count=count+j
             END DO
         END IF
         ibeg = iend
      END DO
  ELSE
      DO i = 1, n, 5
         iend = min(ibeg+5,n)
         ibeg = ibeg + 1
         IF (PRESENT(collab)) THEN
             !write(iout,2) collab(ibeg:iend)
         ELSE
             write (iout,3) (ii,ii=ibeg,iend)
         END IF
         DO j = ibeg, iend
            write (iout,5) a(count + 1:count + j)
            count=count+j
         END DO
         ibeg = iend
      END DO
  END IF
1 format(/,a80,/)
2 format(15x,3(3x,a8,9x))
3 format(1x,'col',3(5x,i6,9x))
4 format(10X,5(1pf20.12))
5 format(10x,5(1pe20.12))
end subroutine  Print_Triangle_Matrix_z
!***********************************************************************
!***********************************************************************
!deck Print_Vector_d
  subroutine Print_Vector_d(a,iout,frmt,title,collab)
  implicit none
  REAL(idp), DIMENSION (:)                    :: a
  INTEGER                                  :: iout
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  IF (PRESENT(title)) THEN
      !write(iout,1) title
  END IF
  IF (PRESENT(frmt)) THEN
      IF (PRESENT(collab)) THEN
          !write(iout,2) collab
      ELSE
          write (iout,3)
      END IF
      IF (frmt == 'f') THEN
          write (iout,4) a(:)
      ELSE
          write (iout,5) a(:)
      END IF
  ELSE
      IF (PRESENT(collab)) THEN
          !write(iout,2) collab
      ELSE
          write (iout,3)
      END IF
      write (iout,5) a(:)
  END IF
1 format(/,a80,/)
2 format(18x,a8)
3 format(18x,'col',5x,i6)
4 format(10x,1pf20.12)
5 format(10x,1pe20.12)
end subroutine  Print_Vector_d
!***********************************************************************
!***********************************************************************
!deck Print_Vector_z
  subroutine Print_Vector_z(a,iout,frmt,title,collab)
  implicit none
  COMPLEX(idp), DIMENSION (:)                :: a
  INTEGER                                  :: iout
  CHARACTER(LEN=*), OPTIONAL               :: frmt
  CHARACTER(LEN=*), OPTIONAL               :: title
  CHARACTER(LEN=*), DIMENSION(:), OPTIONAL :: collab
  IF (PRESENT(title)) THEN
      !write(iout,1) title
  END IF
  IF (PRESENT(frmt)) THEN
      IF (PRESENT(collab)) THEN
          !write(iout,2) collab
      ELSE
          write (iout,3)
      END IF
      IF (frmt == 'f') THEN
          write (iout,4) a(:)
      ELSE
          write (iout,5) a(:)
      END IF
  ELSE
      IF (PRESENT(collab)) THEN
          !write(iout,2) collab
      ELSE
          write (iout,3)
      END IF
      write (iout,5) a(:)
  END IF
1 format(/,a80,/)
2 format(18x,a8)
3 format(18x,'col',5x,i6)
4 format(10x,2(1pf20.12))
5 format(10x,2(1pe20.12))
end subroutine  Print_Vector_z
!***********************************************************************
!***********************************************************************
END MODULE Matrix_Print

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

MODULE Special_Functions
!***begin prologue     Special_Functions
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)                                                                                            
!***keywords           associated legendre functions                                                              
!***author             schneider, b. i.(nsf)                                                                                        
!***source                                                                                                                    
!***purpose            Compute P_lm(x) and Q_lm(x) for all x                                                                     
!***description        For P_lm and Q_lm with x=(-1,1) upward recursion is stable.
!***                   For P_lm for abs(x) > 1 upward recursion is also stable
!***                   For Q_lm for abs(x) > 1 upward recursion is unstable and we use
!***                   downward recursion starting at a high value of and then renormalizing 
!***                   using the analytically known first function.  This is called Millers
!***                   algorithm and is well know in computing Bessel functions.  
!***                                                                                                  
!                                                                                                                                   
!***references                                                                                                                      
!***routines called                                                                                                                 
!***end prologue       Special_Functions                                                                               
!
  USE accuracy
  USE Data_Module
  IMPLICIT NONE
  REAL(idp), DIMENSION(:),                          &
             ALLOCATABLE                  :: x
  REAL(idp), DIMENSION(:),                          &
             ALLOCATABLE                  :: y
  INTEGER                                 :: m_max = 1
  INTEGER                                 :: l_max = 5
  INTEGER                                 :: n_points = 1
  LOGICAL                                 :: normalize = .true.
  LOGICAL                                 :: Derivative = .false.
  LOGICAL                                 :: Print_Functions = .false.
  LOGICAL                                 :: Print_Wronskian = .false.
  LOGICAL                                 :: input_values = .true.
  LOGICAL                                 :: test_wron = .false.
  REAL(idp)                               :: norm
  REAL(idp)                               :: arg
  REAL(idp)                               :: scale_factor
  REAL(idp)                               :: log_factor
  REAL(idp)                               :: wron
  REAL(idp), DIMENSION(:), ALLOCATABLE    :: Factor
  INTEGER                                 :: l
  INTEGER                                 :: m
  INTEGER                                 :: m_sign
  INTEGER                                 :: s_fac
  REAL(idp)                               :: smallest = tiny(1.d0) * 1.d+04
!  REAL(idp)                               :: smallest = 2.5d-300
  REAL(idp)                               :: biggest  = huge(1.d0) * 1.d-04
  REAL(idp)                               :: eps = 1.d-10
  REAL(idp)                               :: upper=1.d0
  REAL(idp)                               :: lower=-1.d0
  REAL(idp)                               :: step
  CHARACTER (LEN=8), DIMENSION(:),                  &
                     ALLOCATABLE          :: row_label
  CHARACTER (LEN=8), DIMENSION(:),                  &
                     ALLOCATABLE          :: col_label
  CHARACTER(LEN=48)                       :: title='Test Calculation'
  CHARACTER(LEN=24)                       :: Control='compute_functions'
  CHARACTER(LEN=24)                       :: recur = 'Miller'
  CHARACTER(LEN=16)                       :: Directive = 'regular'
!
  TYPE Xi
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F_Small
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F_Large
  END TYPE Xi
  TYPE Eta
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F_1
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F_2
  END TYPE Eta

  TYPE Reg_L
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: DF
  END TYPE Reg_L

  TYPE Reg_M
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: DF
  END TYPE Reg_M

  TYPE Reg_LM
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: DF
       real(idp), dimension(:,:,:), allocatable :: P
       TYPE(Xi)                           :: Xi  
       TYPE(Eta)                          :: Eta  
  END TYPE Reg_LM

  TYPE Irreg_L
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: DF
  END TYPE Irreg_L

  TYPE Irreg_M
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: DF
  END TYPE Irreg_M

  TYPE Irreg_LM
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: DF
       real(idp), dimension(:,:,:), allocatable :: Q
       TYPE(Xi)                           :: Xi  
       TYPE(Eta)                          :: Eta  
  END TYPE Irreg_LM

  TYPE Up
       CHARACTER(LEN=24)                  :: Dir
  END TYPE UP

  TYPE Down_A
       CHARACTER(LEN=24)                  :: Dir
  END TYPE Down_A

  TYPE Down_B
       CHARACTER(LEN=24)                  :: Dir
  END TYPE Down_B

  TYPE Down
       TYPE(Down_A)                       :: A  
       TYPE(Down_B)                       :: B  
  END TYPE Down

  TYPE CF_Legendre
       CHARACTER(LEN=24)                  :: Dir
  END TYPE CF_Legendre
 
  TYPE Coefficients
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: a
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: b
  END TYPE Coefficients
!
  TYPE Legendre_Functions
       TYPE(Reg_L)                        :: R_L
       TYPE(Reg_M)                        :: R_M
       TYPE(Reg_LM)                       :: R_LM
       TYPE(Irreg_L)                      :: I_L
       TYPE(Irreg_M)                      :: I_M
       TYPE(Irreg_LM)                     :: I_LM
       TYPE(Up)                           :: U
       TYPE(Down)                         :: D
       TYPE(Coefficients)                 :: C_Q
  END TYPE Legendre_Functions
!
  TYPE(Legendre_Functions)                :: Leg
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Special_Functions

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

MODULE Prolate_Functions
!***begin prologue     Prolate_Functions
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)                                                                                            
!***keywords           prolate functions                                                              
!***author             schneider, b. i.(nsf)                                                                                        
!***source                                                                                                                    
!***purpose            
!***description        
!***                   
!***                   
!***                   
!***                   
!***                   
!***                                                                                                  
!                                                                                                                                   
!***references                                                                                                                      
!***routines called                                                                                                                 
!***end prologue       Prolate_Functions                                                                               
!
  USE accuracy
  USE Data_Module
  IMPLICIT NONE
  INTEGER                                 :: lorder
  INTEGER                                 :: morder
  INTEGER                                 :: mabs
  INTEGER                                 :: lsum
  INTEGER                                 :: msum
  INTEGER                                 :: meo
  REAL(idp)                               :: a
  REAL(idp)                               :: R
  REAL(idp)                               :: radius_moeq
  REAL(idp)                               :: xi1
  REAL(idp)                               :: xi2
  REAL(idp)                               :: eta1
  REAL(idp)                               :: eta2
  REAL(idp)                               :: varphi1
  REAL(idp)                               :: varphi2
  REAL(idp)                               :: x1
  REAL(idp)                               :: y1
  REAL(idp)                               :: z1
  REAL(idp)                               :: rho1
  REAL(idp)                               :: x2
  REAL(idp)                               :: y2
  REAL(idp)                               :: z2
  REAL(idp)                               :: rho2
  REAL(idp)                               :: xi_small
  REAL(idp)                               :: xi_large
  REAL(idp)                               :: facm
  REAL(idp)                               :: vardm
  REAL(idp)                               :: dl21
  REAL(idp)                               :: temp
  REAL(idp)                               :: csum_real
  REAL(idp)                               :: csum_imag
  REAL(idp)                               :: varphi_diff
  REAL(idp)                               :: ctemp_real
  REAL(idp)                               :: ctemp_imag
  REAL(idp)                               :: dx
  REAL(idp)                               :: dy
  REAL(idp)                               :: dz
  REAL(idp)                               :: rsqr
  REAL(idp)                               :: r1
  REAL(idp)                               :: r2
  REAL(idp)                               :: r_12
  REAL(idp)                               :: r_12_invs
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

END MODULE Prolate_Functions
MODULE Lentz_Thompson
!***begin prologue     Lentz_Thompson
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)                                                                                            
!***keywords           Lentz_Thompson                                                              
!***author             schneider, b. i.(nsf)                                                                                        
!***source                                                                                                                    
!***purpose            Compute continued fractions using Lentz-Thompson algoritm                                                                     
!***description        
!***                   
!***                   
!***                   
!***                   
!***                   
!***                                                                                                  
!                                                                                                                                   
!***references                                                                                                                      
!***routines called                                                                                                                 
!***end prologue       Lentz_Thompson                                                                               
!
!                          Needed Modules
!
  USE accuracy
  USE Data_Module
  USE input_output
  USE Matrix_Print
  USE Special_Functions
  IMPLICIT NONE
!
!
                           INTERFACE Continued_Fractions
             MODULE PROCEDURE Continued_Fraction_Legendre
                       END INTERFACE Continued_Fractions                                                                    
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                      
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Continued_Fraction_Legendre
!***begin prologue     Continued_Fraction_Legendre  
!***date written       100206   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           continued fractions
!***author             schneider, barry (NSF)
!***source             
!***purpose            Compute continued fraction
!***description        
!***                   
!***references         none
!***routines called
!***end prologue       Continued_Fraction_Legendre  
  Subroutine Continued_Fraction_Legendre(CFL,f,x,nu,mu) 
  IMPLICIT NONE
  TYPE(CF_Legendre)                 :: CFL
  REAL(idp)                         :: x 
  INTEGER                           :: nu 
  INTEGER                           :: mu 
  INTEGER                           :: n_0 
  REAL(idp)                         :: a 
  REAL(idp)                         :: b 
  REAL(idp)                         :: f        
  REAL(idp)                         :: C        
  REAL(idp)                         :: D        
  REAL(idp)                         :: Del        
  REAL(idp)                         :: test        
  INTEGER                           :: count         
  INTEGER                           :: iwrite         
!
  f = smallest
  C = f
  D = zero
  test = one
  a = one 
  n_0 = nu
  b = zero
!  !write(iout,1) f, C, D, test, a, b, n_0
  count = 0
!  iwrite = 0
  DO While (test > eps )              
     count = count + int_one
!     iwrite = iwrite + int_one
     b = ( ( n_0 + n_0 + one ) * x ) / ( n_0 + mu )
     D = b + a * D
     IF ( D == zero ) THEN
          D = smallest
     END IF
     C = b + a / C
     IF ( C == zero ) THEN
          C = smallest
     END IF
     D = 1.d0 / D
     Del = C * D
     f = f * Del
     test = abs ( Del - one )
     a = - ( n_0 - mu + one ) / ( n_0 + mu )
     n_0 = n_0 + 1
!     IF ( iwrite == 50 ) THEN
!          iwrite = 0
!          !write(iout,2) count, f, C, D, test, a, b, n_0
!     END IF
  END DO
  !write(iout,3) nu, mu, count, test, f 
1 Format(5x, 'Initial Values',/,10x,                                    &
         'f    = ',d20.12,1x,'C = ',d15.8,1x,'D = ',d15.8,/,10x,        &
         'test = ',d15.8,1x,'a  = ',d15.8,1x,'b = ',d15.8,1x,'n_0 = ',i10)
2 Format(5x, 'Loop Values',5x,'Count = ',i5/,10x,                       &
         'f    = ',d20.12,1x,'C = ',d15.8,1x,'D = ',d15.8,/,10x,        &
         'test = ',d15.8,1x,'a  = ',d15.8,1x,'b = ',d15.8,1x,'n_0 = ',i10)
3 Format(1x,'l = ',i3,1x,'m = ',i3,1x,'Iteration Count = ',i5,          &
         /,1x,'Convergence = ',1pe19.11,1x,'Final Value of Continued Fraction = ',1pe19.11)
END SUBROUTINE Continued_Fraction_Legendre
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Lentz_Thompson

!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=!

MODULE Associated_Legendre_Functions
!***begin prologue     Associated_Legendre_Functions
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)                                                  
!***keywords           associated legendre functions                                         
!***author             schneider, b. i.(nsf)                                                 
!***source                                                                                    
!***purpose            Compute P_lm(x) and Q_lm(x) for all x                                 
!***description        See subroutine Info in driver codes for description.
!***references
!***routines called
!***end prologue       Associated_Legendre_Functions
!
!                          Needed Modules
!
  use numbers
  USE accuracy
  USE Data_Module
  USE input_output
  USE Matrix_Print
  USE Special_Functions
  USE Prolate_Functions
  USE Lentz_Thompson
  IMPLICIT NONE
                                                                                                  
                           INTERFACE Legendre
             MODULE PROCEDURE Legendre                              
                       END INTERFACE Legendre
!                                                                                                  
                           INTERFACE Legendre_Recursion
             MODULE PROCEDURE Upward_Regular_Legendre_Recursion_L,                         &
                              Upward_Regular_Legendre_Recursion_LM,                        &
                              Upward_Irregular_Legendre_Recursion_LM,                      &
                              Downward_Irregular_Legendre_Recursion_LM_A,                  &  
                              Downward_Irregular_Legendre_Recursion_LM_B
                       END INTERFACE Legendre_Recursion
!
                           INTERFACE Initialize
             MODULE PROCEDURE Initialize_Regular_L,                                        &
                              Initialize_Regular_LM,                                       &
                              Initialize_Irregular_L,                                      &
                              Initialize_Irregular_LM                              
                       END INTERFACE Initialize
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Legendre_Functions 
!***begin prologue     Legendre_Functions  
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of P_LM(x) functions using various recursion
!***                   relations.  
!***
!***references         The needed relations can be found at the following two references;
!***                   1. http://en.wikipedia.org/wiki/Associated_Legendre_function
!***                   2. Abramowitz and Stegun, Handbook of Mathematical Functions
!***                           UNITED STATES DEPARTMENT OF COMMERCE
!***                                       With
!***                   Formulas, Graphs, and Mathematical Tables
!***                   Edited by Milton Abramowitz and Irene A. Stegun
!***                   National Bureau of Standards, Applied Mathematics Series - 55
!***                   Issued June 1964, Tenth Printing, December 1972, with corrections
!***                   For sale by the Superintendent of Documents, U.S. Government Printing Office 
!***                   Washington, D.C. 20402 - Price $11.35 domestic postpaid, or $10.50 GPO Bookstore 
!***
!***                   Some comments.  The recursion relationships when the argument is inside or 
!***                   ouside the cut [-1,1] are slightly different.  Reference 2. does not contain 
!***                   all of the recurrances needed.  This is noted in the subroutines.
!                      P_LM(z)[Q_LM(z)] are the regular[irregular] associated legendre 
!***                   functions.  z are the real values of the argument.  Between - 1 and + 1 upward
!                      recursion can be used for both the regular and irregular function.
!                      For other values of z upward recursion is fine for the regular
!                      function but downward recusion must be used for the irregular function.
!***routines called
!***end prologue       Legendre_Functions
      Subroutine Legendre ( R_LM, I_LM ) 
      IMPLICIT NONE
      TYPE(Reg_LM), OPTIONAL                         :: R_LM
      TYPE(Irreg_LM), OPTIONAL                       :: I_LM
      TYPE(Up)                                       :: U
      TYPE(Down)                                     :: D
      TYPE(Down_A)                                   :: A
      TYPE(Down_B)                                   :: B
      INTEGER                                        ::  i
      INTEGER                                        ::  j
      INTEGER                                        ::  lrow=1
      INTEGER                                        ::  lcol=1
!      CHARACTER(LEN=3)                               ::  itoc
      CHARACTER(LEN=16)                              ::  fptoc
!
!
!----------------------------------------------------------------------c
!
!
!----------------------------------------------------------------------c
!
!----------------------------------------------------------------------c

!         Set print labels
!
  ALLOCATE(col_label(int_zero:m_max))
  DO i=int_zero,m_max
!     col_label(i) = 'm = '//itoc(i)
     write(col_label(i),'(A,I3)') 'm = ', i
  END DO
  ALLOCATE(row_label(int_zero:l_max))
  DO i=int_zero,l_max
!     row_label(i) = 'l = '//itoc(i)
     write(row_label(i),'(A,I3)') 'l =', i
  END DO
!
!        Allocate some space for storage of often used variables
!
  ALLOCATE(y(int_one:int_twenty))
  DO i = int_one, n_points
!----------------------------------------------------------------------c
!
!    The definition and calculation of Legendre functions depends if the
!    argument is ( -1.,1.) or outside that range.  Using s_fac allows a uniform
!    treatment.
     arg = x(i)
     s_fac = int_one
     IF ( abs(arg) > one ) THEN
          s_fac = - int_one
     END IF
     y(int_one) = one - arg * arg
     y(int_two) =  s_fac * y(int_one)
     y(int_three) = sqrt ( y(int_two) )   
     y(int_four) = y(int_three)   
     y(int_five) = arg * arg
     y(int_six) = y(int_five) * arg
     y(int_seven) = y(int_six) * arg
     y(int_eight) = y(int_seven) * arg 
     y(int_nine) = y(int_eight) * arg 
     y(int_ten) = y(int_nine) * arg 
     y(int_eleven) = y(int_ten) * arg 
     y(int_twelve) = y(int_eleven) * arg 
     IF ( PRESENT(R_LM) ) THEN
!
!
          Call Legendre_Recursion ( Leg%R_LM )
          IF (Print_Functions) THEN
              !write(iout,1) arg
              title='Regular Associated Legendre Functions'
              !write(iout,2) title
              Call Print_Matrix(Leg%R_LM%F, l_max + int_one, m_max + int_one, iout, frmt='e',          &
                                collab=col_label, rowlab=row_label )
          END IF
          IF ( abs(arg) <= one) THEN
               IF (normalize) THEN
!
!                  Normalize
!
                   Call Renormalize( Leg%R_LM%F )
                   IF (Print_Functions) THEN
                       title='Normalized Regular Associated Legendre Functions'
                       !write(iout,2) title
                       Call Print_Matrix(Leg%R_LM%F, l_max + int_one, m_max + int_one, iout, frmt='e',  &
                                         collab=col_label, rowlab=row_label )
                   END IF
               END IF
          END IF
     END IF
!----------------------------------------------------------------------c
!
!----------------------------------------------------------------------c
     IF ( PRESENT(I_LM) ) THEN
!
          IF ( abs(arg) == one) THEN
               !write(iout,3)
               return
          END IF
!          log_factor = log ( abs ( ( arg + one )/( arg - one ) ) )
          log_factor = log ( arg + one ) - log( abs(arg - one ) )
!----------------------------------------------------------------------c
!
!         Starting values for Q_LM upward recursion.
!----------------------------------------------------------------------c
!
          IF( abs(arg) < one ) THEN
              Call Legendre_Recursion( Leg%I_LM, Leg%U )
          ELSE
!----------------------------------------------------------------------c
!           Recur downward for m = 0,1 using either Millers algorithm  c
!           or the continued fraction.
!----------------------------------------------------------------------c
!
               IF ( Leg%D%A%Dir == 'Miller' ) THEN
                    Call Legendre_Recursion(  Leg%I_LM, Leg%D, Leg%D%A )
               ELSE IF (Leg%D%B%Dir == 'Wronskian' ) THEN
                    Call Legendre_Recursion(  Leg%I_LM, Leg%D, Leg%D%B )
               END IF
          END IF
          IF (Print_Functions) THEN
              title='Irregular Associated Legendre Functions'
              !write(iout,1) arg
              !write(iout,2) title
              Call Print_Matrix(Leg%I_LM%F, l_max + int_one, m_max + int_one, iout, frmt='e',   &
                                collab=col_label, rowlab=row_label )
          END IF
     END IF 
     IF ( PRESENT (R_LM) .and. PRESENT(I_LM) .and. test_wron ) THEN
          Call Test_Wronskian ( Leg%R_LM , Leg%I_LM )
     END IF
     Leg%R_LM%P(:,:,i) = Leg%R_LM%F(:,:)
     Leg%I_LM%Q(:,:,i) = Leg%I_LM%F(:,:)
  END DO
  DEALLOCATE( y )
  DEALLOCATE(col_label)
  DEALLOCATE(row_label)
1 Format(/,25x,'Argument = ',f15.8)
2 Format(/,25x,a48)
3 Format(/,25x,'Cannot Compute Irregular Function for Argument One')
END SUBROUTINE Legendre
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Initialize_Regular_L 
!***begin prologue     Initialize_Regular_L    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation starting values of regular Legendre functions
!***                   for upward L recursion when only one M value needed.
!***references         none
!                      
!***routines called
!***end prologue       Initialize_Regular_L         
      Subroutine Initialize_Regular_L ( R_L )
      IMPLICIT NONE
      TYPE ( Reg_L )                              :: R_L
      INTEGER                                     :: n_1
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              
!         Starting values for P_LM.
!              
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!       Use: 
!              P_MM = - s_fac * ( 2*M - 1) *  sqrt ( s_fac * ( 1 - x* x ) ) * P_(M-1)(M-1)
!       To step up in M after initializing at one.  
!
  Leg%R_L%F(m) = one
!
!              Overwrite starting value until the result is obtained.
!
  n_1 = int_one  
  DO l = int_one, m
     Leg%R_L%F(m) = - s_fac * n_1 * y(int_three) * Leg%R_L%F(m)
     n_1 = n_1 + int_two
  END DO
!
!              Calculate the second P term.
!
  IF ( l_max > m ) THEN
!
!          Now calculate:
!                 P_(M+1)M
!
       Leg%R_L%F(m+1) = ( m + m + int_one) * arg * Leg%R_L%F(m)
  END IF
END SUBROUTINE Initialize_Regular_L  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Initialize_Regular_LM 
!***begin prologue     Initialize_Regular_LM    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation starting values of regular Legendre functions
!***                   for upward L recursion when multiple M values are needed.
!***references         none
!                      
!***routines called
!***end prologue       Initialize_Regular_LM         
      Subroutine Initialize_Regular_LM ( R_LM )
      IMPLICIT NONE
      TYPE ( Reg_LM )                             :: R_LM
      INTEGER                                     :: n_1
      INTEGER                                     :: n_2

!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              
!         Starting values for P_LM.
!              
!----------------------------------------------------------------------c
!
!
!----------------------------------------------------------------------c
!       Use: 
!              P_MM = - s_fac * ( 2*M - 1) *  sqrt ( s_fac * ( 1 - x* x ) ) * P_(M-1)(M-1)
!       To step up in M after initializing at one.  
!
  n_1 = int_one 
  DO m = int_zero, m_max              
!
!               Initialize first value.
!
     Leg%R_LM%F(m,m) = one
     n_2 = int_one 
     DO l = int_one, m
        Leg%R_LM%F(m,m) = - s_fac * n_2 * y(int_three) * Leg%R_LM%F(m,m)
        n_2 = n_2 + int_two
     END DO
!
!               Calculate second value.
!
     IF (l_max /= m) THEN
!
!          Now calculate:
!                 P_(M+1)M
!
         Leg%R_LM%F(m+int_one,m) = n_1 * arg * Leg%R_LM%F(m,m)
         n_1 = n_1 + int_two
     END IF
  END DO
END SUBROUTINE Initialize_Regular_LM  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Initialize_Irregular_L 
!***begin prologue     Initialize_Irregular_L    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation starting values of regular Legendre functions
!***                   for upward L recursion when only one M is needed.
!***references         none
!                      
!***routines called
!***end prologue       Initialize_Irregular_L         
      Subroutine Initialize_Irregular_L ( I_L )
      IMPLICIT NONE
      TYPE ( Irreg_L )                            :: I_L
      REAL(idp)                                   :: I_0
!
!
!         Starting values for Q_LM upward recursion.
!              Use:       
!                    Q_00 = .5 * ln ( abs( ( z + 1) /( z - 1))
!                    Q_10 = z * Q_00 - 1
!                    Q_01 = - 1 / sqrt ( s_fac * ( 1 - z * z ) )
!                    Q_11 = - sqrt ( s_fac * ( 1 - z * z ) ) * ( Q_00 + z / ( 1 - z * z ) )
!
!----------------------------------------------------------------------c
!
  IF ( m == int_zero) THEN
       Leg%I_L%F(int_zero) = half * log_factor
       Leg%I_L%F(int_one) = arg * Leg%I_L%F(int_zero) - one
  END IF
  IF ( m == int_one) THEN
       I_0 = half * log_factor
       Leg%I_L%F(int_zero) =  - one / y(int_three)
       Leg%I_L%F(int_one) =  - s_fac * y(int_three) * ( I_0 + arg / y(int_one) )
  END IF
!
END SUBROUTINE Initialize_Irregular_L  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Initialize_Irregular_LM 
!***begin prologue     Initialize_Irregular_LM    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation starting values of regular Legendre functions
!***                   for upward L recursion when multiple M values needed.
!***references         none
!                      
!***routines called
!***end prologue       Initialize_Irregular_LM         
      Subroutine Initialize_Irregular_LM ( I_LM )
      IMPLICIT NONE
      TYPE ( Irreg_LM )                           :: I_LM
!
!
!         Starting values for Q_LM upward recursion.
!              Use:       
!                    Q_00 = .5 * ln ( abs( ( z + 1) /( z - 1))
!                    Q_10 = z * Q_00 - 1
!                    Q_01 = - 1 / sqrt ( s_fac * ( 1 - z * z ) )
!                    Q_11 = - sqrt ( s_fac * ( 1 - z * z ) ) * ( Q_00 + z / ( 1 - z * z ) )
!
!----------------------------------------------------------------------c
!
  Leg%I_LM%F(int_zero,int_zero) = half * log_factor
  Leg%I_LM%F(int_one,int_zero) = arg * Leg%I_LM%F(int_zero,int_zero) - one
  IF ( m_max > int_zero) THEN
       Leg%I_LM%F(int_zero,int_one) =  - one / y(int_three)
       Leg%I_LM%F(int_one,int_one) =  - s_fac * y(int_three) * ( Leg%I_LM%F(int_zero,int_zero) + arg / y(int_one) )
  END IF
!
END SUBROUTINE Initialize_Irregular_LM  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Upward_Regular_Legendre_Recursion_L 
!***begin prologue     Upward_Regular_Legendre_Recursion_L   
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of Q_LM(z) functions using forward recursion in L.
!***                   ( L + 1 - M) Q_(L+1)M = ( 2*L + 1) z Q_LM - ( L + M ) Q_(L-1)M
!***                   Recursion started with the explicit forms of Q_MM and Q_(M+1)M
!***                   The forward recursion is stable for both the regular and irregular
!***                   functions as long as abs(z) <= one and L is not huge.  It is also stable
!***                   for the regular functions for abx(z) > one.  It is NOT stable for the
!***                   irregular functions under these conditions and backward recursion is 
!***                   required.
!***references         none
!                      
!***routines called
!***end prologue       Upward_Regular_Legendre_Recursion_L         
      Subroutine Upward_Regular_Legendre_Recursion_L ( R_L )
      IMPLICIT NONE
      TYPE ( Reg_L )                              :: R_L
      INTEGER                                     ::  l
      INTEGER                                     :: n_1
      INTEGER                                     :: n_2
      INTEGER                                     :: n_3
!
!
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              
!              Starting values for P_LM.
!              
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!
!
  Call Initialize ( Leg%R_L ) 
!
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              Get the other L values by upward recursion in l
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!
!         The upward recursion.  Stable for all values of z
!         This is just the standard l recursion in textbooks.
!
  n_1 = m + m + int_three 
  n_2 = m + m + int_one
  n_3 = int_two
  DO l = m + int_two, l_max
    Leg%R_L%F(l) = ( n_1 * arg * Leg%R_L%F(l - int_one)                &
                                         -                             &
                     n_2 * Leg%R_L%F(l - int_two) ) / n_3  
    n_1 = n_1 + int_two
    n_2 = n_2 + int_one
    n_3 = n_3 + int_one
  END DO
END SUBROUTINE Upward_Regular_Legendre_Recursion_L 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Upward_Regular_Legendre_Recursion_LM 
!***begin prologue     Upward_Regular_Legendre_Recursion_LM    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of Q_LM(z) functions using forward recursion in L.
!***                   ( L + 1 - M) Q_(L+1)M = ( 2*L + 1) z Q_LM - ( L + M ) Q_(L-1)M
!***                   Recursion started with the explicit forms of Q_MM and Q_(M+1)M
!***                   The forward recursion is stable for both the regular and irregular
!***                   functions as long as abs(z) <= one and L is not huge.  It is also stable
!***                   for the regular functions for abx(z) > one.  It is NOT stable for the
!***                   irregular functions under these conditions and backward recursion is 
!***                   required.
!***references         none
!                      
!***routines called
!***end prologue       Upward_Regular_Legendre_Recursion_LM          
      Subroutine Upward_Regular_Legendre_Recursion_LM ( R_LM )
      IMPLICIT NONE
      TYPE ( Reg_LM )                               :: R_LM
      INTEGER                                       ::  l
      INTEGER                                       :: n_1
      INTEGER                                       :: n_2
      INTEGER                                       :: n_3
!
!
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              
!              Starting values for P_LM.
!              
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!
!
  Call Initialize ( Leg%R_LM ) 
!
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              Get the other L values by upward recursion
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!
!
   DO m = int_zero, m_max
!
!         The upward recursion.  Stable for all values of z
!         This is just the standardv reecursion in textbooks.
!
      n_1 = m + m + int_three 
      n_2 = m + m + int_one
      n_3 = int_two
      DO l = m + int_two, l_max
         Leg%R_LM%F(l,m) = ( n_1 * arg * Leg%R_LM%F(l - int_one,m)             &
                                       -                                       &
                             n_2 * Leg%R_LM%F(l - int_two,m) ) / n_3
         n_1 = n_1 + int_two
         n_2 = n_2 + int_one
         n_3 = n_3 + int_one
      END DO
   END DO
END SUBROUTINE Upward_Regular_Legendre_Recursion_LM  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Upward_Irregular_Legendre_Recursion_LM 
!***begin prologue     Upward_Irregular_Legendre_Recursion_LM    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of Q_LM(z) functions using forward recursion in L.
!***                   This is used when abs(arg) is < one.
!***                   ( L + 1 - M) Q_(L+1)M = ( 2*L + 1) z Q_LM - ( L + M ) Q_(L-1)M
!***                   Recursion started with the explicit forms of Q_MM and Q_(M+1)M
!***                   The forward recursion is stable for both the regular and irregular
!***                   functions as long as abs(z) <= one and L is not huge.  It is also stable
!***                   for the regular functions for abx(z) > one.  It is NOT stable for the
!***                   irregular functions under these conditions and backward recursion is required.
!***references         none
!                      
!***routines called
!***end prologue       Upward_Irregular_Legendre_Recursion_LM          
      Subroutine Upward_Irregular_Legendre_Recursion_LM ( I_LM, U )
      IMPLICIT NONE
      TYPE ( Irreg_LM )                             :: I_LM
      TYPE ( Up )                                   :: U
      INTEGER                                       ::  l
      INTEGER                                       :: n_1
      INTEGER                                       :: n_2
      INTEGER                                       :: n_3
!----------------------------------------------------------------------c
!
!         Starting values for Q_00, Q0_1, Q_10 and Q_11 
!         for upward recursion.
!----------------------------------------------------------------------c
!
  !write(iout,*) 'upward recursion'
  Call Initialize( Leg%I_LM )
!
!----------------------------------------------------------------------c
!         Get other values by upward recursion in L starting with
!         the M=0,1 values and then to all L and M by upward
!         recursion on both variables.
!----------------------------------------------------------------------c
!
!         The upward recursion.  Stable for all values of z <= one
!
!----------------------------------------------------------------------c
!         Step up to get all Q_l0 and Q_l1
!----------------------------------------------------------------------c

  n_1 = int_three 
  n_2 = int_one
  n_3 = int_two      
  DO l = int_two, l_max
     Leg%I_LM%F(l,int_zero) = ( n_1 * arg * Leg%I_LM%F(l - int_one,int_zero)          &
                                                  -                                   &
                         n_2 * Leg%I_LM%F(l - int_two,int_zero) ) / n_3
     n_1 = n_1 + int_two
     n_2 = n_2 + int_one
     n_3 = n_3 + int_one
  END DO
  IF ( m_max > int_zero ) THEN
!
       n_1 = int_three 
       n_2 = int_two
       n_3 = int_one      
       DO l = int_two, l_max
          Leg%I_LM%F(l,int_one) = ( n_1 * arg * Leg%I_LM%F(l - int_one,int_one)      &
                                                  -                                  &
                              l * Leg%I_LM%F(l - int_two,int_one) ) / ( l - int_one )
          n_1 = n_1 + int_two
          n_2 = n_2 + int_one
          n_3 = n_3 + int_one
       END DO
  END IF
!
!----------------------------------------------------------------------c
!             Now for each L value, step up in M if needed.
!----------------------------------------------------------------------c
  IF ( m_max > int_one ) THEN
       DO l = int_zero, l_max
!
!                     The upward recursion in m
!
          n_1 = - int_two 
          n_2 = l + int_one
          n_3 = l      
          DO m = int_two, m_max
             Leg%I_LM%F(l,m) = ( - m - m + int_two ) * arg * Leg%I_LM%F(l,m-int_one) / y(3)    &
                                                     -                                   &
                               s_fac * ( l + m - int_one) * ( l - m + int_two) * Leg%I_LM%F(l,m-int_two)
             n_1 = n_1 - int_two
             n_2 = n_2 + int_one
             n_3 = n_3 - int_one
          END DO
       END DO
  END IF
END SUBROUTINE Upward_Irregular_Legendre_Recursion_LM  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Downward_Irregular_Legendre_Recursion_LM_A 
!***begin prologue     Downward_Irregular_Legendre_Recursion_LM_A    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of Q_LM(z) functions using backward recursion
!***                   and the modified Miller algorithm.
!***                   This is used when abs(arg) is > one.
!***                   ( L + 1 - M) T_LM = ( 2*L + 3) z T_(L+1)M - ( L - M + 2 ) T_(L+2)M
!***                   Starting at a large value of L set the last value to zero and the
!***                   next to last to one.  The recur downward which is the stable direction.
!***                   The T's are proportional to the desired Q functions.  
!***                   The proportionality constant is determined by the known value of Q00.
!***                   This allows us to compute the Q's for m=0. The process is repeated 
!***                   for Q_01
!***references         none
!                      
!***routines called
!***end prologue       Downward_Irregular_Legendre_Recursion_LM_A         
      Subroutine Downward_Irregular_Legendre_Recursion_LM_A ( I_LM, D, A )
      IMPLICIT NONE
      TYPE ( Irreg_LM )                             :: I_LM
      TYPE(CF_Legendre)                             :: CFL
      TYPE ( Down )                                 :: D
      TYPE ( Down_A )                               :: A
      REAL(idp)                                     :: I_2 
      REAL(idp)                                     :: I_1 
      REAL(idp)                                     :: I_0 
      INTEGER                                       :: l
      INTEGER                                       :: n_1
      INTEGER                                       :: n_2
      INTEGER                                       :: n_3
      !write(iout,*) 'downward recursion'
!
!     Compute continued fraction Q_L_max / Q_(L_max-1) for m=0
!
      m = int_zero
      Call Continued_Fraction_Legendre(CFL,I_2,arg,l_max,m)    
      Leg%I_LM%F(l_max,m) = I_2 * smallest
      Leg%I_LM%F(l_max-int_one,m) = smallest
!
!     Downward recursion for m = 0
!
      n_1 = m - l_max 
      n_2 = l_max + l_max - int_one
      n_3 = l_max + m - int_one
      DO l = l_max-2, int_zero, - int_one
         Leg%I_LM%F(l,m) = ( n_2 * arg * Leg%I_LM%F(l+int_one,m) + n_1 * Leg%I_LM%F(l+int_two,m) ) / n_3
         n_1 = n_1 + int_one
         n_2 = n_2 - int_two
         n_3 = n_3 - int_one
      END DO
!
!     Renormalize using known value of first member,
!
      scale_factor =  half * log_factor  / Leg%I_LM%F(int_zero,m)
      Leg%I_LM%F(int_zero:l_max,m) = scale_factor * Leg%I_LM%F(int_zero:l_max,m)
      IF ( m_max > int_zero) THEN
!
!         Downward recursion for m = 1
!
!
!         Compute continued fraction Q_L_max / Q_(L_max-1) for m=1
!
          m = int_one
          Call Continued_Fraction_Legendre(CFL,I_2,arg,l_max,m)
          Leg%I_LM%F(l_max,m) = I_2 * smallest
          Leg%I_LM%F(l_max-int_one,m) = smallest
!
!         Downward recursion 
!
          n_1 = m - l_max 
          n_2 = l_max + l_max - int_one
          n_3 = l_max + m - int_one
          DO l = l_max-2, int_zero, - int_one
             Leg%I_LM%F(l,m) = ( n_2 * arg * Leg%I_LM%F(l+int_one,m) + n_1 * Leg%I_LM%F(l+int_two,m) ) / n_3
             n_1 = n_1 + int_one
             n_2 = n_2 - int_two
             n_3 = n_3 - int_one
          END DO
!
!         Renormalize using known value of first member.
!
          scale_factor = ( - one / y(int_three) ) / Leg%I_LM%F(int_zero,m)
          Leg%I_LM%F(int_zero:l_max,m) = scale_factor * Leg%I_LM%F(int_zero:l_max,m)
      END IF
!
!     For each l value, step up in m
!
      DO l = int_zero, l_max
         n_1 = - int_two
         n_2 = l + int_one
         n_3 = l
!
!        The upward recursion in m
!
         DO m = int_two, m_max
            Leg%I_LM%F(l,m) = n_1 * arg * Leg%I_LM%F(l,m - int_one) / y(3)       &
                                           -                                     &
                              s_fac * n_2 * n_3 * Leg%I_LM%F(l, m - int_two)
            n_1 = n_1 - int_two
            n_2 = n_2 + int_one
            n_3 = n_3 - int_one
         END DO
      END DO
!1 Format(/,10x,'Argument = ',e15.8,1x,'Top L = ',i5)
!2 Format(/,25x,a48)
END SUBROUTINE Downward_Irregular_Legendre_Recursion_LM_A  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Downward_Irregular_Legendre_Recursion_LM_B 
!***begin prologue     Downward_Irregular_Legendre_Recursion_LM_B    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of Q_LM(z) functions using a novel idea of Segura and Gil;
!***                   This is used when abs(arg) is < one.
!***                   1. continued fraction for Q_L0/Q_(L-1)0 and Q_L1/Q_(L-1)1
!***                   2. upward recursion for P_L0 and P_L1
!***                   3. the wronskian  P_L0 * Q_(L-1)0 - P_(L-1)0 * Q_L0  = 1 / L 
!***                      wronskian  P_L1 * Q_(L-1)1 - P_(L-1)1 * Q_L1  =  - 1 / L 
!***                      to compute the two  highest values of Q_L0 and Q_L1.  
!***                   4. downward recursion for all Q_L0 and Q_L1.  
!***                   5. upward recursion in M for all QLM.
!***references         Gil and Segura CPC 1998
!                      
!***routines called
!***end prologue       Downward_Irregular_Legendre_Recursion_LM_B         
      Subroutine Downward_Irregular_Legendre_Recursion_LM_B ( I_LM, D, B ) 
      IMPLICIT NONE
      TYPE ( Irreg_LM )                             :: I_LM
      TYPE ( Reg_L )                                :: R_L
      TYPE ( Down )                                 :: D
      TYPE ( Down_B )                               :: B
      TYPE(CF_Legendre)                             :: CFL
      REAL(idp)                                     :: I_2 
      REAL(idp)                                     :: I_1 
      REAL(idp)                                     :: I_0 
      REAL(idp)                                     :: cf 
      INTEGER                                       :: l
      INTEGER                                       :: n_1
      INTEGER                                       :: n_2
      INTEGER                                       :: n_3
      !write(iout,*) 'downward recursion'
      m = int_zero
!
!                      Recur up for P_L0
!
      Call initialize(Leg%R_L)
      Call Legendre_Recursion ( Leg%R_L )
!
!                      Get continued fraction
!
      Call Continued_Fraction_Legendre(CFL,cf,arg,l_max,m)       
      Leg%I_LM%F(l_max-int_one,m) = 1.d0                            &
                                         /                          &
                      ( l_max * ( Leg%R_L%F(l_max) - Leg%R_L%F(l_max-int_one) * cf ) )      
      Leg%I_LM%F(l_max,m) = cf * Leg%I_LM%F(l_max-int_one,m)            
!
!     Downward recursion for m = 0
!
      n_1 = l_max + l_max - int_one
      n_2 = l_max 
      n_3 = l_max - int_one
      DO l = l_max, int_two, - int_one
         Leg%I_LM%F(l-int_two,m) = ( n_1 * arg * Leg%I_LM%F(l-int_one,m)         &
                                                      -                          &
                                            n_2 * Leg%I_LM%F(l,m) ) / n_3 
         n_1 = n_1 - int_two
         n_2 = n_2 - int_one
         n_3 = n_3 - int_one
      END DO
      IF ( m_max > int_zero) THEN
           m = int_one
!
!                      Recur up for P_L1
!
           Call initialize(Leg%R_L)
           Call Legendre_Recursion ( Leg%R_L )
!
!                      Get continued fraction
!
           Call Continued_Fraction_Legendre(CFL,cf,arg,l_max,m)       
           Leg%I_LM%F(l_max-int_one,m) = - l_max / ( Leg%R_L%F(l_max) - Leg%R_L%F(l_max-int_one) * cf )   
           Leg%I_LM%F(l_max,m) = cf * Leg%I_LM%F(l_max-int_one,m)            
!           !write(iout,*) Leg%I_LM%F(l_max,m), Leg%I_LM%F(l_max-int_one,m)
!
!                      Downward recursion for m = 1
!
           n_1 = l_max + l_max - int_one
           n_2 = l_max - int_one
           n_3 = l_max  
           DO l = l_max, int_two, - int_one
              Leg%I_LM%F(l-int_two,m) = ( n_1 * arg * Leg%I_LM%F(l-int_one,m)         &
                                                    -                                 &
                                          n_2 * Leg%I_LM%F(l,m) ) / n_3 
              n_1 = n_1 - int_two
              n_2 = n_2 - int_one
              n_3 = n_3 - int_one
           END DO
      END IF

!
!     For each l value, step up in m
!
      DO l = int_zero, l_max
         n_1 = - int_two
         n_2 = l + int_one
         n_3 = l
!
!        The upward recursion in m
!
         DO m = int_two, m_max
            Leg%I_LM%F(l,m) = n_1 * arg * Leg%I_LM%F(l,m - int_one) / y(3)       &
                                           -                                     &
                              s_fac * n_2 * n_3 * Leg%I_LM%F(l, m - int_two)
            n_1 = n_1 - int_two
            n_2 = n_2 + int_one
            n_3 = n_3 - int_one
         END DO
      END DO
1 Format(/,5x,'m = ',i3,1x,'Continued Fraction = ',d15.8)
2 Format(/,5x,'Starting Values of Legendre Functions')
3 Format(/,5x,5e15.8)
END SUBROUTINE Downward_Irregular_Legendre_Recursion_LM_B  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Test_Wronskian 
!***begin prologue     Test_Wronskian     
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        
!***references         Gil and Segura CPC 1998
!                      
!***routines called
!***end prologue       Test_Wronskian         
      Subroutine Test_Wronskian ( R_LM, I_LM ) 
      IMPLICIT NONE
      TYPE ( Reg_LM )                               :: R_LM
      TYPE ( Irreg_LM )                             :: I_LM
      REAL(idp)                                     :: wr_calc
      REAL(idp)                                     :: max_err
      REAL(idp)                                     :: diff
      m_sign = 1
      DO m = 0, m_max
         max_err = zero
         DO l = m + int_one, l_max
            wron = m_sign * Factor( l + m - int_one )/ Factor (l - m )
            wr_calc = Leg%R_LM%F(l,m) * Leg%I_LM%F(l-int_one,m)            &
                                      -                                    &
                      Leg%R_LM%F(l-int_one,m) * Leg%I_LM%F(l,m)            
            diff = ( wron - wr_calc )/ wron
            max_err = max(max_err,abs(diff))
            IF (Print_Wronskian) THEN
                !write(iout,1) m, l, wron, wr_calc
            END IF
         END DO
         write (iout,2) m, max_err
         IF (abs(arg) > one ) THEN
             m_sign = - m_sign
         END IF
      END DO
!
1 Format(/,5x,'m = ',i3,1x,'l = ',i3,/,10x,                                &
              'Exact Wronskian = ',1pe20.12,1x,'Computed Wronskian = ',1pe20.12)
2 Format(/,5x,'m = ',i3,1x,'Maximum Relative Error in Wronskian for all Computed l Values = ',1pe20.12)
END SUBROUTINE Test_Wronskian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Factorials 
!***begin prologue     Factorials    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            Factorials 
!***description        Factorials
!***                     
!***references         none
!                      
!***routines called
!***end prologue       Factorials         
      Subroutine Factorials 
      IMPLICIT NONE
      INTEGER                  :: i
      Factor(0) = int_one
      Factor(int_one) = int_one
      DO i = int_two, l_max + m_max
         Factor(i) = i * Factor( i - int_one )
      END DO
END SUBROUTINE Factorials 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Renormalize 
!***begin prologue     Renormalize   
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        Normalized regular Legendre functions from their un-normalized values.
!***                   This is only used on the cut (-1,+1).
!***                     
!***references         none
!                      
!***routines called
!***end prologue       Renormalizen        
      Subroutine Renormalize ( F_lm )
      IMPLICIT NONE
      REAL(idp), DIMENSION(0:l_max,int_zero:m_max)            ::  F_lm
      INTEGER                                          :: l_fac
      DO m = int_zero, m_max
         l_fac = m + m + int_one
         DO l = m , l_max
            norm = sqrt ( half * l_fac * Factor ( l - m ) / Factor ( l + m ) )
            l_fac = l_fac +int_two
            F_lm(l,m) = norm * F_lm(l,m)
         END DO
      END DO
END SUBROUTINE Renormalize 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Associated_Legendre_Functions
