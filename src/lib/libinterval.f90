module LibInterval
   !
   ! RIVM - National Institute for Public Health and the Environment
   !
   ! PO Box 1,
   ! NL - 3720 BA Bilthoven
   ! The Netherlands
   use libxmath, only: Float

   implicit none

   private

   public :: IntervalSpecsType,&
   & ExponentialIntervalType,ExponentialIntervalInterpolate

   type IntervalSpecsType
      integer :: N
      real(Float) :: XMin,XMax
   end type IntervalSpecsType

   type ExponentialIntervalType
      real(Float) :: FirstDelay,DelayGrowthFactor
      type(IntervalSpecsType) :: IntervalSpecs
      real(Float), allocatable, dimension(:) :: Values
   end type ExponentialIntervalType
   !
   ! The following parameter can be set for debugging this module
   !
   integer, parameter :: DebugLevel = 0

contains
   function ExponentialIntervalInterpolate(x,Interval,InterpolationWay,Error)
      !
      ! Give the interpolated value in a point arbitrary x>=0 of an interval with values at x=0 and at exponentially
      ! increasing x-values starting from a first delay.
      ! The function supports:
      !      InterpolationWay = 0: zeroth-order,
      !      InterpolationWay = 1: linear and
      !      InterpolationWay = 2: powerlaw interpolation.
      !
      real(Float), intent(in) :: x
      type(ExponentialIntervalType), intent(in) :: Interval
      integer, intent(in) :: InterpolationWay
      logical, intent(out) :: Error
      real(Float) :: ExponentialIntervalInterpolate

      real(Float) :: Dum,LocX,CellX,ThePower,XLeft,XRight,YLeft,YRight
      integer :: NX,iLocX
      logical :: BothAreNonZero

      integer :: DebugLevel = 0

      Error = .FALSE.

      NX = Interval%IntervalSpecs%N
      if (DebugLevel > 0) write(*,'(A,I0)') 'ExponentialIntervalInterpolate: NX = ',NX

      if (x < 0._Float) then
         Dum = Interval%Values(0)
         Error = .TRUE. ! You are violating the conditions!
         if (DebugLevel > 0) write(*,'(A,F15.5)') 'ExponentialIntervalInterpolate: 1: ',Dum
      else if (x > Interval%IntervalSpecs%XMax) then
         Dum = Interval%Values(NX)
         Error = .TRUE. ! You are violating the conditions!
         if (DebugLevel > 0) write(*,'(A,F15.5)') 'ExponentialIntervalInterpolate: 2: ',Dum
      else if (x < Interval%FirstDelay) then
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
      else
         LocX = 1._Float + LOG(x/Interval%FirstDelay)/LOG(Interval%DelayGrowthFactor)
         iLocX = MAX(MIN(INT(LocX),NX-1),1)
         if (DebugLevel > 0) write(*,'(A,I0)') 'ExponentialIntervalInterpolate: iLocX = ',iLocX
         call FLUSH(6)

         XLeft  = Interval%FirstDelay*Interval%DelayGrowthFactor**(iLocX-1)
         XRight = XLeft*Interval%DelayGrowthFactor
         YLeft  = Interval%Values(iLocX)
         YRight = Interval%Values(iLocX+1)
         if (DebugLevel > 0) then
            write(*,'(A,F15.2)') 'ExponentialIntervalInterpolate: XLeft  = ',XLeft
            write(*,'(A,F15.2)') 'ExponentialIntervalInterpolate: XRight = ',XRight
            write(*,'(A,F15.2)') 'ExponentialIntervalInterpolate: YLeft  = ',YLeft
            write(*,'(A,F15.2)') 'ExponentialIntervalInterpolate: YRight = ',YRight
            call FLUSH(6)
         endif
         !
         ! Powerlaw interpolation will not work if 1 of the y-values is 0. In that case use linear interpolation.
         !
         BothAreNonZero = (abs(YLeft*YRight) > 0._Float)
         if (DebugLevel > 0) write(*,'(A,L1)') 'ExponentialIntervalInterpolate: BothAreNonZero? ',BothAreNonZero
         call FLUSH(6)
         !
         ! Nearest neighbour interpolation
         !
         if (InterpolationWay == 0) then

            CellX = (x-XLeft)/(XRight-XLeft)

            if (CellX <= 0.5_Float) then
               Dum = YLeft
            else
               Dum = YRight
            endif
               if (DebugLevel > 0) then
               write(*,'(A,F10.2,A,EN20.10)') 'ExponentialIntervalInterpolate: Order 0, CellX = ',CellX,'    Dum = ',Dum
               call FLUSH(6)
            endif
         else if ((InterpolationWay == 1).OR.(.NOT.BothAreNonZero)) then
            !
            ! Linear interpolation y = ax + b
            !
            CellX = (x-XLeft)/(XRight-XLeft)

            Dum = (1._Float-CellX) * YLeft + CellX *YRight
            if (DebugLevel > 0) then
               write(*,'(A,F10.2,A,EN20.10)') 'ExponentialIntervalInterpolate: Linear: CellX = ',CellX,'    Dum = ',Dum
               call FLUSH(6)
            endif
         else if (InterpolationWay == 2) then
            !
            ! Powerlaw y = ax^b
            !
            if (BothAreNonZero) then
               ThePower = LOG(YRight/YLeft)/LOG(XRight/XLeft)

               Dum = YLeft * (x/XLeft)**ThePower
               if (DebugLevel > 0) then
                  write(*,'(A,F10.2,A,EN20.10)') 'ExponentialIntervalInterpolate: Powerlaw: ThePower = ',&
                  & ThePower,'     Dum = ',Dum
                  call FLUSH(6)
               endif
            else
               write(*,'(A)') 'ExponentialIntervalInterpolate: You are not supposed to end up here...'
            endif
         else
            write(*,'(A,I0)') 'Invalid value for InterpolationWay : ',InterpolationWay
            call exit()
         endif
      endif
      !
      ! Neglect very small values to avoid the possible writing of a number with an exponent that has 3 digits:
      ! the "E" is omitted in such cases, making the output unreadable for other programs...
      !
      if (abs(Dum) < 1.E-95_Float) Dum = 0._Float

      ExponentialIntervalInterpolate = Dum

   end function ExponentialIntervalInterpolate
end module LibInterval
