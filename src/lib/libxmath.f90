MODULE LibXMath
   !
   ! Extended mathematics
   !
   ! ____________________________________________________
   !
   ! Developed For:
   !
   ! RIVM - National Institute for Public Health and the Environment
   !
   ! PO Box 1,
   ! NL - 3720 BA Bilthoven
   ! The Netherlands
   ! Website: http://www.rivm.nl/
   !
   ! ____________________________________________________
   !
   IMPLICIT NONE

   PRIVATE

   PUBLIC :: Float, Vector2, Inc,Vector,dp,&
      & OPERATOR(*), OPERATOR(/), OPERATOR(+), OPERATOR(-),&
      & SparseMatrix,SparseMatrixElement,Matrix2SparseMatrix,SparseMatrix2Matrix,&
      & SparseLogicalMatrixElement,SparseLogicalMatrix,LMatrix2SparseLMatrix,SparseLMatrix2LMatrix
   !
   ! From now on only use "REAL(float)" for reals
   !
   ! The number of bytes associated with 6 digits precision, usually stored in 4 bytes, default REAL
   !
   INTEGER,PARAMETER :: Float  = SELECTED_REAL_KIND(12, 60)
   !
   ! The number of bytes associated with 15 digits precision, usually stored in 8 bytes, default DOUBLE PRECISION
   !
   INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)
   !
   ! 2D problems can use the following type:
   !
   TYPE Vector2
      REAL(Float) :: x(2)
   END TYPE Vector2

   TYPE Vector
      REAL(Float) :: x(3)
   END TYPE Vector

   INTERFACE OPERATOR (*)
      MODULE PROCEDURE MultiplyFloatVector2
      MODULE PROCEDURE MultiplyFloatVector
      MODULE PROCEDURE MultiplyVector2Float
      MODULE PROCEDURE MultiplyVectorFloat
      MODULE PROCEDURE MultiplyVector2Vector2
      MODULE PROCEDURE MultiplyVectorVector
   END INTERFACE

   INTERFACE OPERATOR (/)
      MODULE PROCEDURE DivideVector2Float
      MODULE PROCEDURE DivideVectorFloat
      MODULE PROCEDURE DivideVector2Vector2
      MODULE PROCEDURE DivideVectorVector
   END INTERFACE

   INTERFACE OPERATOR (+)
      MODULE PROCEDURE AddVector2Vector2
      MODULE PROCEDURE AddVectorVector
   END INTERFACE

   INTERFACE OPERATOR (-)
      MODULE PROCEDURE SubtractVector2Vector2
      MODULE PROCEDURE SubtractVectorVector
   END INTERFACE

   TYPE SparseMatrixElement
      INTEGER :: i,j
      REAL(Float) :: x
   END TYPE SparseMatrixElement

   TYPE SparseMatrix
      INTEGER :: N,NMax ! N is the number of elements, NMax >=N can be used to allocate on beforehand
      TYPE(SparseMatrixElement), ALLOCATABLE, DIMENSION(:) :: Element
   END TYPE SparseMatrix

   TYPE SparseLogicalMatrixElement
      INTEGER :: i,j
   END TYPE SparseLogicalMatrixElement

   TYPE SparseLogicalMatrix
      INTEGER :: N,NMax ! N is the number of elements, NMax >=N can be used to allocate on beforehand
      TYPE(SparseLogicalMatrixElement), ALLOCATABLE, DIMENSION(:) :: Element
   END TYPE SparseLogicalMatrix
   !
   ! The following parameter can be set for debugging this MODULE
   !
   INTEGER, PARAMETER :: DebugLevel = 0
   !
   ! Below this line you'll get the real implementations
   !
CONTAINS
   FUNCTION MultiplyFloatVector2( xFloat, xVector2 )
      !
      ! Multiply a REAL(Float) with a Vector2
      !
      TYPE( Vector2 ) MultiplyFloatVector2, DumVec2
      TYPE( Vector2 ), INTENT( IN ) :: xVector2
      REAL( Float ), INTENT( IN ) :: xFloat

      DumVec2 = Vector2((/xFloat*xVector2%x(1),xFloat*xVector2%x(2)/))
      MultiplyFloatVector2 = DumVec2
   END FUNCTION MultiplyFloatVector2



   FUNCTION MultiplyFloatVector( xFloat, xVector )
      !
      ! Multiply a REAL(Float) with a Vector
      !
      TYPE( Vector ) MultiplyFloatVector, DumVec
      TYPE( Vector ), INTENT( IN ) :: xVector
      REAL( Float ), INTENT( IN ) :: xFloat

      DumVec = Vector((/xFloat*xVector%x(1),xFloat*xVector%x(2),xFloat*xVector%x(3)/))
      MultiplyFloatVector = DumVec
   END FUNCTION MultiplyFloatVector


   FUNCTION MultiplyVector2Float( xVector2,xFloat )
      !
      ! Multiply a Vector2 with a REAL(Float)
      !
      TYPE( Vector2 ) MultiplyVector2Float
      TYPE( Vector2 ), INTENT( IN ) :: xVector2
      REAL( Float ), INTENT( IN ) :: xFloat

      MultiplyVector2Float = xFloat*xVector2
   END FUNCTION MultiplyVector2Float


   FUNCTION MultiplyVectorFloat( xVector,xFloat )
      !
      ! Multiply a Vector with a REAL(Float)
      !
      TYPE( Vector ) MultiplyVectorFloat
      TYPE( Vector ), INTENT( IN ) :: xVector
      REAL( Float ), INTENT( IN ) :: xFloat

      MultiplyVectorFloat = xFloat*xVector
   END FUNCTION MultiplyVectorFloat


   FUNCTION DivideVector2Float( xVector2,xFloat )
      !
      ! Divide a Vector2 by a REAL(Float)
      !
      TYPE( Vector2 ) DivideVector2Float
      TYPE( Vector2 ), INTENT( IN ) :: xVector2
      REAL( Float ), INTENT( IN ) :: xFloat

      DivideVector2Float = xVector2*(1./xFloat)
   END FUNCTION DivideVector2Float


   FUNCTION DivideVectorFloat( xVector,xFloat )
      !
      ! Divide a Vector by a REAL(Float)
      !
      TYPE( Vector ) DivideVectorFloat
      TYPE( Vector ), INTENT( IN ) :: xVector
      REAL( Float ), INTENT( IN ) :: xFloat

      DivideVectorFloat = xVector*(1./xFloat)
   END FUNCTION DivideVectorFloat


   FUNCTION MultiplyVector2Vector2( xVector2, yVector2 )
      !
      ! Multiply two Vector2 type of variables
      !
      TYPE( Vector2 ) MultiplyVector2Vector2, DumVec2
      TYPE( Vector2 ), INTENT( IN ) :: xVector2,yVector2

      DumVec2 = Vector2((/xVector2%x(1)*yVector2%x(1),xVector2%x(2)*yVector2%x(2)/))
      MultiplyVector2Vector2 = DumVec2
   END FUNCTION MultiplyVector2Vector2


   FUNCTION DivideVector2Vector2( xVector2, yVector2 )
      !
      ! Divide two Vector2 type of variables
      !
      TYPE( Vector2 ) DivideVector2Vector2, DumVec2
      TYPE( Vector2 ), INTENT( IN ) :: xVector2,yVector2

      DumVec2 = Vector2((/xVector2%x(1)/yVector2%x(1),xVector2%x(2)/yVector2%x(2)/))
      DivideVector2Vector2 = DumVec2
   END FUNCTION DivideVector2Vector2


   FUNCTION DivideVectorVector( xVector, yVector )
      !
      ! Divide two Vector type of variables
      !
      TYPE( Vector ) DivideVectorVector, DumVec
      TYPE( Vector ), INTENT( IN ) :: xVector,yVector

      DumVec = Vector((/xVector%x(1)/yVector%x(1),xVector%x(2)/yVector%x(2),xVector%x(3)/yVector%x(3)/))
      DivideVectorVector = DumVec
   END FUNCTION DivideVectorVector


   FUNCTION MultiplyVectorVector( xVector, yVector )
      !
      ! Multiply two Vector type of variables
      !
      TYPE( Vector ) MultiplyVectorVector, DumVec
      TYPE( Vector ), INTENT( IN ) :: xVector,yVector

      DumVec = Vector((/xVector%x(1)*yVector%x(1),xVector%x(2)*yVector%x(2),xVector%x(3)*yVector%x(3)/))
      MultiplyVectorVector = DumVec
   END FUNCTION MultiplyVectorVector


   FUNCTION AddVector2Vector2( xVector2, yVector2 )
      !
      ! Add two Vector2 type of variables
      !
      TYPE( Vector2 ) AddVector2Vector2, DumVec2
      TYPE( Vector2 ), INTENT( IN ) :: xVector2,yVector2

      DumVec2 = Vector2((/xVector2%x(1)+yVector2%x(1),xVector2%x(2)+yVector2%x(2)/))
      AddVector2Vector2 = DumVec2
   END FUNCTION AddVector2Vector2


   FUNCTION AddVectorVector( xVector, yVector )
      !
      ! Add two Vector type of variables
      !
      TYPE( Vector ) AddVectorVector, DumVec
      TYPE( Vector ), INTENT( IN ) :: xVector,yVector

      DumVec = Vector((/xVector%x(1)+yVector%x(1),xVector%x(2)+yVector%x(2),xVector%x(3)+yVector%x(3)/))
      AddVectorVector = DumVec
   END FUNCTION AddVectorVector


   FUNCTION SubtractVector2Vector2( xVector2, yVector2 )
      !
      ! Subtract two Vector2 type of variables
      !
      TYPE( Vector2 ) SubtractVector2Vector2
      TYPE( Vector2 ), INTENT( IN ) :: xVector2,yVector2

      SubtractVector2Vector2 = xVector2 + (-1.0_Float)*yVector2
   END FUNCTION SubtractVector2Vector2


   FUNCTION SubtractVectorVector( xVector, yVector )
      !
      ! Subtract two Vector type of variables
      !
      TYPE( Vector ) SubtractVectorVector
      TYPE( Vector ), INTENT( IN ) :: xVector,yVector

      SubtractVectorVector = xVector + (-1.0_Float)*yVector
   END FUNCTION SubtractVectorVector



   SUBROUTINE Inc(x)
      INTEGER, INTENT(INOUT) :: x
      x = x+1
   END SUBROUTINE Inc



   FUNCTION Matrix2SparseMatrix(A)
      !
      ! Convert regular matrix to sparse matrix
      !
      REAL(Float), DIMENSION(:,:), INTENT(IN) :: A

      TYPE(SparseMatrix) :: ASparse,Matrix2SparseMatrix
      INTEGER :: N,i,j,NNonZero,iElement

      N = SIZE(A,1)
      !
      ! Count non-zero elements
      !
      NNonZero = 0
      DO i = 1,N
         DO j = 1,N
            IF (A(j,i).NE.0._Float) NNonZero = NNonZero + 1
         ENDDO
      ENDDO

      ASparse%N = NNonZero
      !
      ! Fill sparse array
      !
      IF (ALLOCATED(ASparse%Element)) DEALLOCATE(ASparse%Element)
      ASparse%NMax = NNonZero
      ALLOCATE(ASparse%Element(NNonZero))

      iElement = 0

      DO i = 1,N
         DO j = 1,N
            IF (A(j,i).NE.0._Float) THEN
            iElement = iElement + 1
            ASparse%Element(iElement)%i = i
            ASparse%Element(iElement)%j = j
            ASparse%Element(iElement)%x = A(j,i)
            ENDIF
         ENDDO
      ENDDO
      !
      ! Pass result
      !
      Matrix2SparseMatrix = ASparse
   END FUNCTION Matrix2SparseMatrix



   SUBROUTINE SparseMatrix2Matrix(ASparse,A)
      !
      ! Convert sparse matrix to regular matrix
      !
      TYPE(SparseMatrix), INTENT(IN) :: ASparse
      REAL(Float), DIMENSION(:,:), INTENT(INOUT) :: A

      INTEGER :: iElement

      INTEGER, PARAMETER :: DebugLevel = 0
      !
      ! Fill array
      !
      A = 0._Float

      IF (DebugLevel.GT.0) THEN
         WRITE(*,'(A,L1,1X,I0,1X,I0)') 'SparseMatrix2Matrix: ',ALLOCATED(ASparse%Element),ASparse%NMax,ASparse%N
      ENDIF

      DO iElement = 1,ASparse%N
         IF (DebugLevel.GT.0) THEN
            WRITE(*,'(A,I4,1X,I4,1X,EN20.10)') 'i,j,x = ',ASparse%Element(iElement)
         ENDIF
         A(ASparse%Element(iElement)%j,ASparse%Element(iElement)%i) = ASparse%Element(iElement)%x
      ENDDO
   END SUBROUTINE SparseMatrix2Matrix



   FUNCTION LMatrix2SparseLMatrix(A)
      !
      ! Convert regular matrix to sparse matrix
      !
      LOGICAL, DIMENSION(:,:), INTENT(IN) :: A

      TYPE(SparseLogicalMatrix) :: LASparse,LMatrix2SparseLMatrix
      INTEGER :: N,i,j,NNonZero,iElement

      N = SIZE(A,1)
      !
      ! Count non-zero elements
      !
      NNonZero = 0
      DO i = 1,N
         DO j = 1,N
            IF (A(j,i)) NNonZero = NNonZero + 1
         ENDDO
      ENDDO

      LASparse%N = NNonZero
      !
      ! Fill sparse array
      !
      IF (ALLOCATED(LASparse%Element)) DEALLOCATE(LASparse%Element)
      LASparse%NMax = NNonZero
      ALLOCATE(LASparse%Element(NNonZero))

      iElement = 0

      DO i = 1,N
         DO j = 1,N
            IF (A(j,i)) THEN
            iElement = iElement + 1
            LASparse%Element(iElement)%i = i
            LASparse%Element(iElement)%j = j
            ENDIF
         ENDDO
      ENDDO
      !
      ! Pass result
      !
      LMatrix2SparseLMatrix = LASparse
   END FUNCTION LMatrix2SparseLMatrix



   SUBROUTINE SparseLMatrix2LMatrix(LASparse,LA)
      !
      ! Convert sparse matrix to regular matrix
      !
      TYPE(SparseLogicalMatrix), INTENT(IN) :: LASparse
      LOGICAL, DIMENSION(:,:), INTENT(INOUT) :: LA

      INTEGER :: iElement

      INTEGER, PARAMETER :: DebugLevel = 0
      !
      ! Fill array
      !
      LA = .FALSE.

      IF (DebugLevel.GT.0) THEN
         WRITE(*,'(A,L1,1X,I0,1X,I0)') 'SparseLMatrix2LMatrix: ',ALLOCATED(LASparse%Element),LASparse%NMax,LASparse%N
      ENDIF

      DO iElement = 1,LASparse%N
         IF (DebugLevel.GT.0) THEN
            WRITE(*,'(A,I4,1X,I4)') 'i,j = ',LASparse%Element(iElement)%i,LASparse%Element(iElement)%j
         ENDIF
         LA(LASparse%Element(iElement)%j,LASparse%Element(iElement)%i) = .TRUE.
      ENDDO
   END SUBROUTINE SparseLMatrix2LMatrix
END MODULE LibXMath
