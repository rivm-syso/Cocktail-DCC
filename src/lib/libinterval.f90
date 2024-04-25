MODULE LibInterval
   !
   ! RIVM - National Institute for Public Health and the Environment
   !
   ! PO Box 1,
   ! NL - 3720 BA Bilthoven
   ! The Netherlands
   USE libxmath
   USE libutil

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: IntervalSpecsType,&
   & ExponentialIntervalType,ExponentialIntervalInterpolate

   TYPE IntervalSpecsType
      INTEGER :: N
      REAL(Float) :: XMin,XMax
   END TYPE IntervalSpecsType

   TYPE ExponentialIntervalType
      REAL(Float) :: FirstDelay,DelayGrowthFactor
      TYPE(IntervalSpecsType) :: IntervalSpecs
      REAL(Float), ALLOCATABLE, DIMENSION(:) :: Values
   END TYPE ExponentialIntervalType
   !
   ! The following parameter can be set for debugging this MODULE
   !
   INTEGER, PARAMETER :: DebugLevel = 0

CONTAINS
   FUNCTION ExponentialIntervalInterpolate(x,Interval,InterpolationWay,Error)
      !
      ! Give the interpolated value in a point arbitrary x>=0 of an interval with values at x=0 and at exponentially
      ! increasing x-values starting from a first delay.
      ! The function supports:
      !      InterpolationWay = 0: zeroth-order,
      !      InterpolationWay = 1: linear and
      !      InterpolationWay = 2: powerlaw interpolation.
      !
      REAL(Float), INTENT(IN) :: x
      TYPE(ExponentialIntervalType), INTENT(IN) :: Interval
      INTEGER, INTENT(IN) :: InterpolationWay
      LOGICAL, INTENT(OUT) :: Error
      REAL(Float) :: ExponentialIntervalInterpolate

      REAL(Float) :: Dum,LocX,CellX,ThePower,XLeft,XRight,YLeft,YRight,dx
      INTEGER :: i,NX,iLocX
      LOGICAL :: BothAreNonZero

      INTEGER :: DebugLevel = 0

      Error = .FALSE.

      NX = Interval%IntervalSpecs%N
      IF (DebugLevel.GT.0) WRITE(*,'(A,I0)') 'ExponentialIntervalInterpolate: NX = ',NX

      IF (x.LT.0._Float) THEN
         Dum = Interval%Values(0)
         Error = .TRUE. ! You are violating the conditions!
         IF (DebugLevel.GT.0) WRITE(*,'(A,F15.5)') 'ExponentialIntervalInterpolate: 1: ',Dum
      ELSE IF (x.GT.Interval%IntervalSpecs%XMax) THEN
         Dum = Interval%Values(NX)
         Error = .TRUE. ! You are violating the conditions!
         IF (DebugLevel.GT.0) WRITE(*,'(A,F15.5)') 'ExponentialIntervalInterpolate: 2: ',Dum
      ELSE IF (x.LT.Interval%FirstDelay) THEN
         !
         ! For any x in [0,FirstDelay> use first order interpolation y = ax + b, irrespective of the given interpolation order
         !
         XLeft  = 0._Float
         XRight = Interval%FirstDelay
         YLeft  = Interval%Values(0)
         YRight = Interval%Values(1)
         CellX = (x-XLeft)/(XRight-XLeft)

         Dum = (1._Float-CellX) * YLeft + CellX *YRight
         !
         ! Find the indices of the two samples around x: iLocX and iLocX+1
         !
      ELSE
         LocX = 1._Float + LOG(x/Interval%FirstDelay)/LOG(Interval%DelayGrowthFactor)
         iLocX = MAX(MIN(INT(LocX),NX-1),1)
         IF (DebugLevel.GT.0) WRITE(*,'(A,I0)') 'ExponentialIntervalInterpolate: iLocX = ',iLocX
         CALL FLUSH(6)

         XLeft  = Interval%FirstDelay*Interval%DelayGrowthFactor**(iLocX-1)
         XRight = XLeft*Interval%DelayGrowthFactor
         YLeft  = Interval%Values(iLocX)
         YRight = Interval%Values(iLocX+1)
         IF (DebugLevel.GT.0) THEN
            WRITE(*,'(A,F15.2)') 'ExponentialIntervalInterpolate: XLeft  = ',XLeft
            WRITE(*,'(A,F15.2)') 'ExponentialIntervalInterpolate: XRight = ',XRight
            WRITE(*,'(A,F15.2)') 'ExponentialIntervalInterpolate: YLeft  = ',YLeft
            WRITE(*,'(A,F15.2)') 'ExponentialIntervalInterpolate: YRight = ',YRight
            CALL FLUSH(6)
         ENDIF
         !
         ! Powerlaw interpolation will not work if 1 of the y-values is 0. In that case use linear interpolation.
         !
         BothAreNonZero = (ABS(YLeft*YRight).GT.0._Float)
         IF (DebugLevel.GT.0) WRITE(*,'(A,L1)') 'ExponentialIntervalInterpolate: BothAreNonZero? ',BothAreNonZero
         CALL FLUSH(6)
         !
         ! Nearest neighbour interpolation
         !
         IF (InterpolationWay.EQ.0) THEN

            CellX = (x-XLeft)/(XRight-XLeft)

            IF (CellX.LE.0.5_Float) THEN
               Dum = YLeft
            ELSE
               Dum = YRight
            ENDIF
               IF (DebugLevel.GT.0) THEN
               WRITE(*,'(A,F10.2,A,EN20.10)') 'ExponentialIntervalInterpolate: Order 0, CellX = ',CellX,'    Dum = ',Dum
               CALL FLUSH(6)
            ENDIF
         ELSE IF ((InterpolationWay.EQ.1).OR.(.NOT.BothAreNonZero)) THEN
            !
            ! Linear interpolation y = ax + b
            !
            CellX = (x-XLeft)/(XRight-XLeft)

            Dum = (1._Float-CellX) * YLeft + CellX *YRight
            IF (DebugLevel.GT.0) THEN
               WRITE(*,'(A,F10.2,A,EN20.10)') 'ExponentialIntervalInterpolate: Linear: CellX = ',CellX,'    Dum = ',Dum
               CALL FLUSH(6)
            ENDIF
         ELSE IF (InterpolationWay.EQ.2) THEN
            !
            ! Powerlaw y = ax^b
            !
            IF (BothAreNonZero) THEN
               ThePower = LOG(YRight/YLeft)/LOG(XRight/XLeft)

               Dum = YLeft * (x/XLeft)**ThePower
               IF (DebugLevel.GT.0) THEN
                  WRITE(*,'(A,F10.2,A,EN20.10)') 'ExponentialIntervalInterpolate: Powerlaw: ThePower = ',&
                  & ThePower,'     Dum = ',Dum
                  CALL FLUSH(6)
               ENDIF
            ELSE
               WRITE(*,'(A)') 'ExponentialIntervalInterpolate: You are not supposed to end up here...'
            ENDIF
         ELSE
            WRITE(*,'(A,I0)') 'Invalid value for InterpolationWay : ',InterpolationWay
            CALL EXIT()
         ENDIF
      ENDIF
      !
      ! Neglect very small values to avoid the possible writing of a number with an exponent that has 3 digits:
      ! the "E" is omitted in such cases, making the output unreadable for other programs...
      !
      IF (ABS(Dum).LT.1.E-95_Float) Dum = 0._Float

      ExponentialIntervalInterpolate = Dum

   END FUNCTION ExponentialIntervalInterpolate
END MODULE LibInterval
