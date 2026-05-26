module LibXMath
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
   implicit none

   private

   public :: Float, Vector2, Inc,Vector,dp,&
      & OPERATOR(*), OPERATOR(/), OPERATOR(+), OPERATOR(-),&
      & SparseMatrix,SparseMatrixElement,Matrix2SparseMatrix,SparseMatrix2Matrix,&
      & SparseLogicalMatrixElement,SparseLogicalMatrix,LMatrix2SparseLMatrix,SparseLMatrix2LMatrix
   !
   ! From now on only use "REAL(float)" for reals
   !
   ! The number of bytes associated with 6 digits precision, usually stored in 4 bytes, default real
   !
   integer,parameter :: Float  = selected_real_kind(12, 60)
   !
   ! The number of bytes associated with 15 digits precision, usually stored in 8 bytes, default DOUBLE PRECISION
   !
   integer, parameter  :: dp = selected_real_kind(12, 60)
   !
   ! 2D problems can use the following type:
   !
   type Vector2
      real(Float) :: x(2)
   end type Vector2

   type Vector
      real(Float) :: x(3)
   end type Vector

   interface OPERATOR (*)
      module PROCEDURE MultiplyFloatVector2
      module PROCEDURE MultiplyFloatVector
      module PROCEDURE MultiplyVector2Float
      module PROCEDURE MultiplyVectorFloat
      module PROCEDURE MultiplyVector2Vector2
      module PROCEDURE MultiplyVectorVector
   end interface

   interface OPERATOR (/)
      module PROCEDURE DivideVector2Float
      module PROCEDURE DivideVectorFloat
      module PROCEDURE DivideVector2Vector2
      module PROCEDURE DivideVectorVector
   end interface

   interface OPERATOR (+)
      module PROCEDURE AddVector2Vector2
      module PROCEDURE AddVectorVector
   end interface

   interface OPERATOR (-)
      module PROCEDURE SubtractVector2Vector2
      module PROCEDURE SubtractVectorVector
   end interface

   type SparseMatrixElement
      integer :: i,j
      real(Float) :: x
   end type SparseMatrixElement

   type SparseMatrix
      integer :: N,NMax ! N is the number of elements, NMax >=N can be used to allocate on beforehand
      type(SparseMatrixElement), allocatable, dimension(:) :: Element
   end type SparseMatrix

   type SparseLogicalMatrixElement
      integer :: i,j
   end type SparseLogicalMatrixElement

   type SparseLogicalMatrix
      integer :: N,NMax ! N is the number of elements, NMax >=N can be used to allocate on beforehand
      type(SparseLogicalMatrixElement), allocatable, dimension(:) :: Element
   end type SparseLogicalMatrix
   !
   ! The following parameter can be set for debugging this module
   !
   integer, parameter :: DebugLevel = 0
   !
   ! Below this line you'll get the real implementations
   !
contains
   function MultiplyFloatVector2( xFloat, xVector2 )
      !
      ! Multiply a real(Float) with a Vector2
      !
      type( Vector2 ) MultiplyFloatVector2, DumVec2
      type( Vector2 ), intent( in ) :: xVector2
      real( Float ), intent( in ) :: xFloat

      DumVec2 = Vector2((/xFloat*xVector2%x(1),xFloat*xVector2%x(2)/))
      MultiplyFloatVector2 = DumVec2
   end function MultiplyFloatVector2



   function MultiplyFloatVector( xFloat, xVector )
      !
      ! Multiply a real(Float) with a Vector
      !
      type( Vector ) MultiplyFloatVector, DumVec
      type( Vector ), intent( in ) :: xVector
      real( Float ), intent( in ) :: xFloat

      DumVec = Vector((/xFloat*xVector%x(1),xFloat*xVector%x(2),xFloat*xVector%x(3)/))
      MultiplyFloatVector = DumVec
   end function MultiplyFloatVector


   function MultiplyVector2Float( xVector2,xFloat )
      !
      ! Multiply a Vector2 with a real(Float)
      !
      type( Vector2 ) MultiplyVector2Float
      type( Vector2 ), intent( in ) :: xVector2
      real( Float ), intent( in ) :: xFloat

      MultiplyVector2Float = xFloat*xVector2
   end function MultiplyVector2Float


   function MultiplyVectorFloat( xVector,xFloat )
      !
      ! Multiply a Vector with a real(Float)
      !
      type( Vector ) MultiplyVectorFloat
      type( Vector ), intent( in ) :: xVector
      real( Float ), intent( in ) :: xFloat

      MultiplyVectorFloat = xFloat*xVector
   end function MultiplyVectorFloat


   function DivideVector2Float( xVector2,xFloat )
      !
      ! Divide a Vector2 by a real(Float)
      !
      type( Vector2 ) DivideVector2Float
      type( Vector2 ), intent( in ) :: xVector2
      real( Float ), intent( in ) :: xFloat

      DivideVector2Float = xVector2*(1./xFloat)
   end function DivideVector2Float


   function DivideVectorFloat( xVector,xFloat )
      !
      ! Divide a Vector by a real(Float)
      !
      type( Vector ) DivideVectorFloat
      type( Vector ), intent( in ) :: xVector
      real( Float ), intent( in ) :: xFloat

      DivideVectorFloat = xVector*(1./xFloat)
   end function DivideVectorFloat


   function MultiplyVector2Vector2( xVector2, yVector2 )
      !
      ! Multiply two Vector2 type of variables
      !
      type( Vector2 ) MultiplyVector2Vector2, DumVec2
      type( Vector2 ), intent( in ) :: xVector2,yVector2

      DumVec2 = Vector2((/xVector2%x(1)*yVector2%x(1),xVector2%x(2)*yVector2%x(2)/))
      MultiplyVector2Vector2 = DumVec2
   end function MultiplyVector2Vector2


   function DivideVector2Vector2( xVector2, yVector2 )
      !
      ! Divide two Vector2 type of variables
      !
      type( Vector2 ) DivideVector2Vector2, DumVec2
      type( Vector2 ), intent( in ) :: xVector2,yVector2

      DumVec2 = Vector2((/xVector2%x(1)/yVector2%x(1),xVector2%x(2)/yVector2%x(2)/))
      DivideVector2Vector2 = DumVec2
   end function DivideVector2Vector2


   function DivideVectorVector( xVector, yVector )
      !
      ! Divide two Vector type of variables
      !
      type( Vector ) DivideVectorVector, DumVec
      type( Vector ), intent( in ) :: xVector,yVector

      DumVec = Vector((/xVector%x(1)/yVector%x(1),xVector%x(2)/yVector%x(2),xVector%x(3)/yVector%x(3)/))
      DivideVectorVector = DumVec
   end function DivideVectorVector


   function MultiplyVectorVector( xVector, yVector )
      !
      ! Multiply two Vector type of variables
      !
      type( Vector ) MultiplyVectorVector, DumVec
      type( Vector ), intent( in ) :: xVector,yVector

      DumVec = Vector((/xVector%x(1)*yVector%x(1),xVector%x(2)*yVector%x(2),xVector%x(3)*yVector%x(3)/))
      MultiplyVectorVector = DumVec
   end function MultiplyVectorVector


   function AddVector2Vector2( xVector2, yVector2 )
      !
      ! Add two Vector2 type of variables
      !
      type( Vector2 ) AddVector2Vector2, DumVec2
      type( Vector2 ), intent( in ) :: xVector2,yVector2

      DumVec2 = Vector2((/xVector2%x(1)+yVector2%x(1),xVector2%x(2)+yVector2%x(2)/))
      AddVector2Vector2 = DumVec2
   end function AddVector2Vector2


   function AddVectorVector( xVector, yVector )
      !
      ! Add two Vector type of variables
      !
      type( Vector ) AddVectorVector, DumVec
      type( Vector ), intent( in ) :: xVector,yVector

      DumVec = Vector((/xVector%x(1)+yVector%x(1),xVector%x(2)+yVector%x(2),xVector%x(3)+yVector%x(3)/))
      AddVectorVector = DumVec
   end function AddVectorVector


   function SubtractVector2Vector2( xVector2, yVector2 )
      !
      ! Subtract two Vector2 type of variables
      !
      type( Vector2 ) SubtractVector2Vector2
      type( Vector2 ), intent( in ) :: xVector2,yVector2

      SubtractVector2Vector2 = xVector2 + (-1.0_Float)*yVector2
   end function SubtractVector2Vector2


   function SubtractVectorVector( xVector, yVector )
      !
      ! Subtract two Vector type of variables
      !
      type( Vector ) SubtractVectorVector
      type( Vector ), intent( in ) :: xVector,yVector

      SubtractVectorVector = xVector + (-1.0_Float)*yVector
   end function SubtractVectorVector



   subroutine Inc(x)
      integer, intent(inout) :: x
      x = x+1
   end subroutine Inc



   function Matrix2SparseMatrix(A)
      !
      ! Convert regular matrix to sparse matrix
      !
      real(Float), dimension(:,:), intent(in) :: A

      type(SparseMatrix) :: ASparse,Matrix2SparseMatrix
      integer :: N,i,j,NNonZero,iElement

      N = SIZE(A,1)
      !
      ! Count non-zero elements
      !
      NNonZero = 0
      do i = 1,N
         do j = 1,N
            if (A(j,i) /= 0._Float) NNonZero = NNonZero + 1
         enddo
      enddo

      ASparse%N = NNonZero
      !
      ! Fill sparse array
      !
      if (ALLOCATED(ASparse%Element)) deallocate(ASparse%Element)
      ASparse%NMax = NNonZero
      allocate(ASparse%Element(NNonZero))

      iElement = 0

      do i = 1,N
         do j = 1,N
            if (A(j,i) /= 0._Float) then
            iElement = iElement + 1
            ASparse%Element(iElement)%i = i
            ASparse%Element(iElement)%j = j
            ASparse%Element(iElement)%x = A(j,i)
            endif
         enddo
      enddo
      !
      ! Pass result
      !
      Matrix2SparseMatrix = ASparse
   end function Matrix2SparseMatrix



   subroutine SparseMatrix2Matrix(ASparse,A)
      !
      ! Convert sparse matrix to regular matrix
      !
      type(SparseMatrix), intent(in) :: ASparse
      real(Float), dimension(:,:), intent(inout) :: A

      integer :: iElement

      integer, parameter :: DebugLevel = 0
      !
      ! Fill array
      !
      A = 0._Float

      if (DebugLevel > 0) then
         write(*,'(A,L1,1X,I0,1X,I0)') 'SparseMatrix2Matrix: ',ALLOCATED(ASparse%Element),ASparse%NMax,ASparse%N
      endif

      do iElement = 1,ASparse%N
         if (DebugLevel > 0) then
            write(*,'(A,I4,1X,I4,1X,EN20.10)') 'i,j,x = ',ASparse%Element(iElement)
         endif
         A(ASparse%Element(iElement)%j,ASparse%Element(iElement)%i) = ASparse%Element(iElement)%x
      enddo
   end subroutine SparseMatrix2Matrix



   function LMatrix2SparseLMatrix(A)
      !
      ! Convert regular matrix to sparse matrix
      !
      logical, dimension(:,:), intent(in) :: A

      type(SparseLogicalMatrix) :: LASparse,LMatrix2SparseLMatrix
      integer :: N,i,j,NNonZero,iElement

      N = SIZE(A,1)
      !
      ! Count non-zero elements
      !
      NNonZero = 0
      do i = 1,N
         do j = 1,N
            if (A(j,i)) NNonZero = NNonZero + 1
         enddo
      enddo

      LASparse%N = NNonZero
      !
      ! Fill sparse array
      !
      if (ALLOCATED(LASparse%Element)) deallocate(LASparse%Element)
      LASparse%NMax = NNonZero
      allocate(LASparse%Element(NNonZero))

      iElement = 0

      do i = 1,N
         do j = 1,N
            if (A(j,i)) then
            iElement = iElement + 1
            LASparse%Element(iElement)%i = i
            LASparse%Element(iElement)%j = j
            endif
         enddo
      enddo
      !
      ! Pass result
      !
      LMatrix2SparseLMatrix = LASparse
   end function LMatrix2SparseLMatrix



   subroutine SparseLMatrix2LMatrix(LASparse,LA)
      !
      ! Convert sparse matrix to regular matrix
      !
      type(SparseLogicalMatrix), intent(in) :: LASparse
      logical, dimension(:,:), intent(inout) :: LA

      integer :: iElement

      integer, parameter :: DebugLevel = 0
      !
      ! Fill array
      !
      LA = .FALSE.

      if (DebugLevel > 0) then
         write(*,'(A,L1,1X,I0,1X,I0)') 'SparseLMatrix2LMatrix: ',ALLOCATED(LASparse%Element),LASparse%NMax,LASparse%N
      endif

      do iElement = 1,LASparse%N
         if (DebugLevel > 0) then
            write(*,'(A,I4,1X,I4)') 'i,j = ',LASparse%Element(iElement)%i,LASparse%Element(iElement)%j
         endif
         LA(LASparse%Element(iElement)%j,LASparse%Element(iElement)%i) = .TRUE.
      enddo
   end subroutine SparseLMatrix2LMatrix
end module LibXMath
