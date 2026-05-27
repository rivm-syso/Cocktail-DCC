module LibPinPoint
   ! ____________________________________________________
   !
   ! Developed For:
   !
   ! VLH   Centre for Environmental Safety and Security
   ! RIVM - National Institute for Public Health and the Environment
   !
   ! PO Box 1,
   ! NL - 3720 BA Bilthoven
   ! The Netherlands
   !
   ! ____________________________________________________
   !
   use LibDCC, only: DCCEffectiveDose, DCCGround, DCCInhalation, DCCThyroidInhalation, &
      InitLibDCC, &
      InAir, Age_Adult, Inhalation_public_adult, &
      ReadTissueDCCs, ReadDCCEffectiveDose, ReadGroundDCCs, ReadInhalationDCCs, ReadThyroidInhalationDCCs
   use LibENDF, only: MaxNuclides, NNuclides, NuclideSpecs, GetNuclideNumber, &
      MotherDaughterMatrix, &
      MassNuc2NucMass, EnsureHyphen, &
      TransitionMatrixPath
   use LibUtil, only: DefaultLength, ScratchFile, &
      AllUpCase
   use LibXMath, only: Float, SparseMatrix

   implicit none

   private

   public :: CocktailType,InitLibPinPoint,NStartingTimes,StartingTimeName,&
   & AvailableDelay,TransitionMatrix,GetPinpointCocktails,&
   & MakePinpointDoseRates,ReadRIVMSourceTermFile,SaveCocktail2,FirstDelay,DelayGrowthFactor,MatureCocktail

   type CocktailType
      character(DefaultLength) :: MyName,MyDirectory
      real(Float), dimension(MaxNuclides) :: x
   end type CocktailType

   real(Float), parameter :: FirstDelay = 60._Float ! seconds
   real(Float), parameter :: DelayGrowthFactor = 1.15_Float
   integer, parameter :: NStartingTimes = 200 ! There is also pinpoint 0, which is 0

   character(15), dimension(0:NStartingTimes) :: StartingTimeName
   real(Float), dimension(0:NStartingTimes) :: AvailableDelay ! In seconds

   type(SparseMatrix), dimension(0:NStartingTimes) :: TransitionMatrix

contains
   subroutine InitLibPinPoint(UseICRP)
      !
      ! Prepare some settings before use of this library
      !
      logical, intent(in) :: UseICRP

      integer :: iStartingTime

      integer, parameter :: DebugLevel = 0
      !
      ! Get the DCCs
      !
      call InitLibDCC(UseICRP) ! includes reading data from ICRP or ENDF
      call ReadTissueDCCs()
      call ReadDCCEffectiveDose()
      call ReadGroundDCCs()
      call ReadInhalationDCCs()
      call ReadThyroidInhalationDCCs()

      AvailableDelay(0) = 0
      AvailableDelay(1) = 60
      do iStartingTime = 2,NStartingTimes
         AvailableDelay(iStartingTime) = DelayGrowthFactor*AvailableDelay(iStartingTime-1)
      enddo ! loop over starting times

      if (DebugLevel > 0) write(*,*)
      if (DebugLevel > 0) write(*,'(A)') 'Considering the following exponentially distributed delays:'
      do iStartingTime = 0,NStartingTimes
         write(StartingTimeName(iStartingTime),'(E11.5)') AvailableDelay(iStartingTime)
         if (DebugLevel > 0) write(*,'(A)') StartingTimeName(iStartingTime)//'second'
      enddo ! loop over starting times

      call GetTransitionMatrices(UseICRP)
   end subroutine InitLibPinPoint



   subroutine ReadRIVMSourceTermFile(FName,MyCocktail)
      !
      ! Read sourceterm with RIVM file format
      !
      type(CocktailType), intent(out) :: MyCocktail
      character(*), intent(in) :: FName

      character(10) :: MyName
      character(DefaultLength) :: ALine
      integer :: iNuclide,TheIndex,iCharacter
      logical :: Ready,IsNotABackSlash
      real(Float) :: MyActivity,ScalingFactor

      integer, parameter :: DebugLevel = 0
      !
      ! Extract name of the source term from filename
      !
      TheIndex = INDEX(FName,'.RIVMSource')
      if (TheIndex == 0) then
         write(*,'(A)') 'Sourcterm "'//trim(FName)//'" is not of type .RIVMSource! Exiting!!'
         call exit()
      endif
      !
      ! Extract directory
      !
      MyCocktail%MyDirectory = FName
      iCharacter = len_trim(MyCocktail%MyDirectory)
      IsNotABackSlash = .NOT.(MyCocktail%MyDirectory(iCharacter:iCharacter) == '/')
      do while (IsNotABackSlash.AND.(iCharacter > 0))
         MyCocktail%MyDirectory(iCharacter:iCharacter) = ' '
         iCharacter = iCharacter - 1
         if (iCharacter > 0) IsNotABackSlash = .NOT.(MyCocktail%MyDirectory(iCharacter:iCharacter) == '/')
      enddo
      !
      ! Extract name of sourceterm
      !
      MyCocktail%MyName = FName((len_trim(MyCocktail%MyDirectory)+1):(TheIndex-1))
      !
      ! Read file
      !
      MyCocktail%x = 0._Float

      open(ScratchFile,file=FName,form='FORMATTED',action='READ')
      !
      ! Read header to find meta-info, e.g. scaling factor for activity
      !
      ScalingFactor = 1._Float

      Ready = .FALSE.
      do while (.NOT.Ready)
         read(ScratchFile,'(A)') ALine
         if (DebugLevel > 1) write(*,'(A)') 'ALine = "'//trim(ALine)//'"'

         TheIndex = INDEX(ALine,'<scalingfactor>')
         if (TheIndex > 0) then
            read(ALine((TheIndex+15):len_trim(ALine)),*) ScalingFactor
            if (DebugLevel > 1) write(*,'(A,EN15.5)') 'Found scaling factor ',ScalingFactor
         endif

         Ready = (ALine(1:1) /= '!')
      enddo
      !
      ! Read the nuclide vector and take into account the optional scaling factor
      !
      write(*,'(A)') 'Reading source term at t=0 from '//trim(FName)

      Ready = .FALSE.

      do while (.NOT.Ready)
         if (DebugLevel > 1) write(*,'(A)') 'ALine = "'//trim(ALine)//'"'

         if (ALine(1:1) /= '!') then
            read(ALine,*) MyName,MyActivity

            call MassNuc2NucMass(MyName)
            call AllUpCase(MyName)
            call EnsureHyphen(MyName)

            iNuclide = GetNuclideNumber(MyName)

            if (iNuclide == 0) then
               write(*,'(A)') 'Could not recognize nuclide "'//trim(MyName)//'"! Exiting!'
               call exit()
            endif
            MyCocktail%x(iNuclide) = ScalingFactor * MyActivity ! [Bq/10kt]
            if (DebugLevel > 0) then
               write(*,'(A,A,I4,A,EN20.10,A)') MyName,' = nuclide ',iNuclide,&
               & ' and has activity ',MyCocktail%x(iNuclide),' Bq/10kt'
            endif
         endif ! not a comment line stat starts with !

         read(ScratchFile,'(A)',end=10) ALine

         Ready = (len_trim(ALine) == 0)
      enddo ! loop until ready

      10 continue
      close(ScratchFile)
   end subroutine ReadRIVMSourceTermFile



   subroutine GetTransitionMatrices(UseICRP)
      !
      ! Get transition matrices for all pinpoints
      !
      logical, intent(in) :: UseICRP

      character(DefaultLength) :: FName,ALine
      integer :: iStartingTime,iLine
      character(6) :: PinPointName

      integer, parameter :: DebugLevel = 0

      write(*,'(A)') 'Getting transition matrices...'
      do iStartingTime = 1,NStartingTimes
         !
         ! Open file
         !
         write(PinPointName,'(I3.3)') iStartingTime

         if (UseICRP) then
            FName = 'SparseMatrix_'//trim(StartingTimeName(iStartingTime))//'second.dat'
         else
            FName = 'SparseMatrix_ENDF_'//trim(PinPointName)//'.dat'
         endif

         if (DebugLevel > 0) write(*,'(A)') 'Reading '//trim(FName)
         call FLUSH(6)
         open(ScratchFile, file=TransitionMatrixPath() // '/' // trim(FName), action='READ')

         if (.NOT.UseICRP) then
            do iLine=1,7
            read(ScratchFile,*) ! Skip header
            enddo
         endif
         !
         ! Count number of sparse matrix elements to allocate space for
         !
         TransitionMatrix(iStartingTime)%N = 0
         do
            read(ScratchFile,'(A)',end=10) ALine
            TransitionMatrix(iStartingTime)%N = TransitionMatrix(iStartingTime)%N + 1
         enddo
         10  continue

         if (DebugLevel > 0) write(*,'(A,I0,A)') trim(FName)//' has ',TransitionMatrix(iStartingTime)%N,' elements'

         REWIND(ScratchFile)
         !
         ! Rewind file and read matrix elements
         !
         TransitionMatrix(iStartingTime)%NMax = TransitionMatrix(iStartingTime)%N
         allocate(TransitionMatrix(iStartingTime)%Element(TransitionMatrix(iStartingTime)%NMax))

         TransitionMatrix(iStartingTime)%N = 0

         if (.NOT.UseICRP) then
            do iLine=1,7
            read(ScratchFile,*) ! Skip header
            enddo
         endif

         do
            read(ScratchFile,'(A)',end=20) ALine
            TransitionMatrix(iStartingTime)%N = TransitionMatrix(iStartingTime)%N + 1
            read(ALine,*) TransitionMatrix(iStartingTime)%Element(TransitionMatrix(iStartingTime)%N)
         enddo
         20  continue
         close(ScratchFile)
      enddo ! loop over starting times
   end subroutine GetTransitionMatrices



   subroutine SaveCocktail2(Cocktail,Tag)
      !
      ! Save a cocktail of radionuclides to file
      !
      type(CocktailType), intent(in) :: Cocktail
      character(*), intent(in) :: Tag

      real(Float), parameter :: NegligibleActivity = 0._Float ! 10.E9_Float ! Bq/10kt, seems to characterize difference between tables 2 and 5
      character(DefaultLength) :: FName
      integer :: iNuclide

      FName = trim(Cocktail%MyName)//'_'//trim(Tag)//'.dat'
      open(ScratchFile,file=trim(Cocktail%MyDirectory)//trim(FName),action='WRITE')
      write(ScratchFile,'(A)') 'Nuclide    Activity[Bq]'
      do iNuclide = 1,NNuclides
         if (Cocktail%x(iNuclide) > NegligibleActivity) then
            write(ScratchFile,'(A,1X,EN11.2)') NuclideSpecs(iNuclide)%NuclideName,Cocktail%x(iNuclide)
         endif
      enddo
      close(ScratchFile)
   end subroutine SaveCocktail2




   subroutine MatureCocktail(InCocktail,TheMatrix,WithProgeny,OutCocktail)
      !
      ! Estimate the cocktail at a given delay from the cocktail at t=0 and the transition matrix for the given delay
      !
      type(CocktailType), intent(in) :: InCocktail
      type(SparseMatrix), intent(in) :: TheMatrix
      logical, intent(in) :: WithProgeny
      type(CocktailType), intent(out) :: OutCocktail

      integer :: iElement,iMother,iDaughter,iNuclide
      real(Float) :: x

      integer, parameter :: DebugLevel = 0

      if (DebugLevel > 0) then
         do iNuclide = 1,NNuclides
            if (InCocktail%x(iNuclide) > 0._Float) then
               write(*,'(A,A,A,I4,A,EN20.10,A)') 'MatureCocktail: ',&
            & NuclideSpecs(iNuclide)%NuclideName,' = nuclide ',iNuclide,&
               &   ' and has activity ',InCocktail%x(iNuclide),' Bq'
            endif
         enddo
      endif ! debug

      OutCocktail%x = 0._Float
      OutCocktail%MyName = trim(InCocktail%MyName)//'_aged'
      OutCocktail%MyDirectory = InCocktail%MyDirectory

      do iElement = 1,TheMatrix%N
         iMother = TheMatrix%Element(iElement)%i
         iDaughter = TheMatrix%Element(iElement)%j
         x = TheMatrix%Element(iElement)%x
         if (WithProgeny.OR.(iMother == iDaughter)) then

            if ((DebugLevel > 0).AND.(InCocktail%x(iMother) > 0._Float)) then
               write(*,'(A,A,A,A,A,EN15.5,A,EN15.5,A,EN15.5)') 'Maturization of ',&
               & NuclideSpecs(iDaughter)%NuclideName,&
               & ' because of initial presence of ',&
               & NuclideSpecs(iMother)%NuclideName ,&
               & ' : ',&
               & OutCocktail%x(iDaughter),&
               & ' + ',x*InCocktail%x(iMother),' = ',&
               & OutCocktail%x(iDaughter) + x*InCocktail%x(iMother)
            endif ! debug

            OutCocktail%x(iDaughter) = OutCocktail%x(iDaughter) + x*InCocktail%x(iMother)
            if (OutCocktail%x(iDaughter) < 1.E-95_Float) OutCocktail%x(iDaughter) = 0._Float ! to prevent exponents with 3 digits, which print ugly...
         endif
      enddo ! loop over elements
   end subroutine MatureCocktail



   subroutine MakePinpointDoseRates(MyStartCocktail,WithDaughters,MotherNature,DaughterNature)
      !
      ! Estimate dose contributions at all pinpoints
      ! Add dose contribution to the head of chain, assume ground activity to be fully on top = "planar"
      !
      type(CocktailType), intent(in) :: MyStartCocktail
      logical, intent(in) :: WithDaughters
      integer, intent(in) :: MotherNature,DaughterNature ! 0 = any, 1 = noble gas, 2 = not a noble gas

      real(Float), dimension(0:NStartingTimes,MaxNuclides) :: PinAirDoseRate,PinGroundDoseRate,&
      & PinInhalationDoseRate,PinThyroidInhalationDoseRate
      integer, dimension(MaxNuclides) :: ParticipatingNuclide
      integer :: NParticipatingNuclides
      real(Float), dimension(MaxNuclides) :: SumAirDoseRate
      integer :: iPinpoint,iNuclide,iDaughter,iMother,iElement,jNuclide,NProgeny
      real(Float) :: DaughterActivity,x
      real(Float), dimension(0:NStartingTimes) :: TotalAirDoseRate,TotalGroundDoseRate,TotalImmersionDoseRate,&
      & TotalThyroidDoseRate
      character(DefaultLength) :: FName,WithDaughterString,NaturesName
      logical :: NaturesOkay,MotherNatureOkay,DaughterNatureOkay

      integer, parameter :: AnyNature = 0
      integer, parameter :: NobleNature = 1
      integer, parameter :: NonNobleNature = 2

      integer, parameter :: DebugLevel = 1

      PinAirDoseRate = 0._Float
      PinGroundDoseRate = 0._Float
      PinInhalationDoseRate = 0._Float
      PinThyroidInhalationDoseRate = 0._Float
      SumAirDoseRate = 0._Float

      if (WithDaughters) then
         WithDaughterString = 'withprogeny'
      else
         WithDaughterString = 'noprogeny'
      endif

      do iPinpoint = 0,NStartingTimes

         if (iPinpoint == 0) then
            ! According to pinpoint 0, there is only the head of chain, which is any fission product available at t=0.
            do iNuclide = 1,NNuclides

               MotherNatureOkay =      (MotherNature == AnyNature)&
               &                  .OR.((MotherNature == NobleNature   ).AND.(NuclideSpecs(iNuclide)%NuclideGroup == 1)) &
               &                  .OR.((MotherNature == NonNobleNature).AND.(NuclideSpecs(iNuclide)%NuclideGroup /= 1))
               DaughterNatureOkay =    (DaughterNature == AnyNature)&
               &                  .OR.((DaughterNature == NobleNature   ).AND.(NuclideSpecs(iNuclide)%NuclideGroup == 1)) &
               &                  .OR.((DaughterNature == NonNobleNature).AND.(NuclideSpecs(iNuclide)%NuclideGroup /= 1))

               NaturesOkay = MotherNatureOkay .AND. DaughterNatureOkay

               if (NaturesOkay) then
                  DaughterActivity = MyStartCocktail%x(iNuclide)

                  PinAirDoseRate(iPinpoint,iNuclide) = PinAirDoseRate(iPinpoint,iNuclide) &
                  & +   DaughterActivity &
                  &   * DCCEffectiveDose(InAir,Age_Adult)%x(iNuclide) &
                  &   * 1.E-9_Float/3600._Float ! Conversion of nSv/h to Sv/s

                  PinGroundDoseRate(iPinpoint,iNuclide) = PinGroundDoseRate(iPinpoint,iNuclide) &
                  & +   DaughterActivity &
                  &   * DCCGround(1,Age_Adult)%x(iNuclide) &
                  &   * 1.E-9_Float/3600._Float ! Conversion of nSv/h to Sv/s
                  !
                  ! Units in the statement below are in Bq * (Sv/Bq) * (h/s)
                  ! Multiplication with the thinning factor in /m3 and the breathing rate in m3/h
                  ! gives the inhalation dose rate in Sv/s.
                  !
                  PinInhalationDoseRate(iPinpoint,iNuclide) = PinInhalationDoseRate(iPinpoint,iNuclide) &
                  & +   DaughterActivity &
                  &   * DCCInhalation(Inhalation_public_adult)%x(iNuclide) &
                  &   * (1._Float/3600._Float)

                  PinThyroidInhalationDoseRate(iPinpoint,iNuclide) = PinThyroidInhalationDoseRate(iPinpoint,iNuclide) &
                  & +   DaughterActivity &
                  &   * DCCThyroidInhalation%x(iNuclide) &
                  &   * (1._Float/3600._Float)

                  SumAirDoseRate(iNuclide) = SumAirDoseRate(iNuclide) + PinAirDoseRate(iPinpoint,iNuclide)
               endif ! natures okay
            enddo ! loop over nuclides
         else
            !
            ! For any pinpoint > 0, you have the full decay chain with progeny under each head of chain (formed at t=0)
            !
            do iElement = 1,TransitionMatrix(iPinpoint)%N
               iMother = TransitionMatrix(iPinpoint)%Element(iElement)%i
               iDaughter = TransitionMatrix(iPinpoint)%Element(iElement)%j

               MotherNatureOkay =      (MotherNature == AnyNature)&
               &                  .OR.((MotherNature == NobleNature   ).AND.(NuclideSpecs(iMother)%NuclideGroup == 1)) &
               &                  .OR.((MotherNature == NonNobleNature).AND.(NuclideSpecs(iMother)%NuclideGroup /= 1))
               DaughterNatureOkay =    (DaughterNature == AnyNature)&
               &                  .OR.((DaughterNature == NobleNature   ).AND.(NuclideSpecs(iDaughter)%NuclideGroup == 1)) &
               &                  .OR.((DaughterNature == NonNobleNature).AND.(NuclideSpecs(iDaughter)%NuclideGroup /= 1))

               NaturesOkay = MotherNatureOkay .AND. DaughterNatureOkay

               if (NaturesOkay) then

                  if (WithDaughters.OR.(iMother == iDaughter)) then
                     x = TransitionMatrix(iPinpoint)%Element(iElement)%x

                     DaughterActivity = x * MyStartCocktail%x(iMother)

                     PinAirDoseRate(iPinpoint,iMother) = PinAirDoseRate(iPinpoint,iMother) &
                     & +   DaughterActivity &
                     &   * DCCEffectiveDose(InAir,Age_Adult)%x(iDaughter) &
                     &   * 1.E-9_Float/3600._Float ! Conversion of nSv/h to Sv/s

                     PinGroundDoseRate(iPinpoint,iMother) = PinGroundDoseRate(iPinpoint,iMother) &
                     & +   DaughterActivity &
                     &   * DCCGround(1,Age_Adult)%x(iDaughter) &
                     &   * 1.E-9_Float/3600._Float ! Conversion of nSv/h to Sv/s

                     PinInhalationDoseRate(iPinpoint,iMother) = PinInhalationDoseRate(iPinpoint,iMother) &
                     & +   DaughterActivity &
                     &   * DCCInhalation(Inhalation_public_adult)%x(iDaughter) &
                     &   * (1._Float/3600._Float) ! Conversion of Sv/(m3/s) to Sv/(m3/h) to facilitate breathing rates in m3/h

                     PinThyroidInhalationDoseRate(iPinpoint,iMother) = PinThyroidInhalationDoseRate(iPinpoint,iMother) &
                     & +   DaughterActivity &
                     &   * DCCThyroidInhalation%x(iDaughter) &
                     &   * (1._Float/3600._Float) ! Conversion of Sv/(m3/s) to Sv/(m3/h) to facilitate breathing rates in m3/h

                     SumAirDoseRate(iMother) = SumAirDoseRate(iMother) + PinAirDoseRate(iPinpoint,iMother)
                  endif ! withdaughters
               endif ! natures okay

            enddo ! loop over sparse matrix at iLower
         endif ! first pinpoint is t=0

      enddo ! loop over all pinpoints
      !
      ! Cut off values < 0.1E-99, as they are printed in an ugly way...
      !
      do iPinpoint = 0,NStartingTimes
         do iNuclide = 1,NNuclides
            if (PinAirDoseRate(iPinpoint,iNuclide) <= 0.1E-99_Float) &
            & PinAirDoseRate(iPinpoint,iNuclide) = 0._Float
            if (PinGroundDoseRate(iPinpoint,iNuclide) <= 0.1E-99_Float) &
            & PinGroundDoseRate(iPinpoint,iNuclide) = 0._Float
            if (PinInhalationDoseRate(iPinpoint,iNuclide) <= 0.1E-99_Float) &
            & PinInhalationDoseRate(iPinpoint,iNuclide) = 0._Float
            if (PinThyroidInhalationDoseRate(iPinpoint,iNuclide) <= 0.1E-99_Float) &
            & PinThyroidInhalationDoseRate(iPinpoint,iNuclide) = 0._Float
         enddo ! loop over nuclides
      enddo ! loop over all pinpoints
      !
      ! Check which nuclides contribute at all on basis of sum of air doserates.
      ! Ground and inhalation rates are not considered to give a different answer.
      !
      NParticipatingNuclides = 0
      do iNuclide = 1,NNuclides
         if (SumAirDoseRate(iNuclide) > 0._Float) then
            NParticipatingNuclides = NParticipatingNuclides + 1
            ParticipatingNuclide(NParticipatingNuclides) = iNuclide
            if (DebugLevel > 0) then
               NProgeny = 0
               do jNuclide = 1,NNuclides
                  if ((MotherDaughterMatrix(jNuclide,iNuclide) /= 0._Float).AND.(jNuclide /= iNuclide)) then
                     NProgeny = NProgeny + 1
                  endif ! found progeny
               enddo ! loop over possible progeny of head of chain
               if (NProgeny == 0) then
                  write(*,'(I4,1X,A7,7X,EN15.5,A)') NParticipatingNuclides,NuclideSpecs(iNuclide)%NuclideName,&
                  & MyStartCocktail%x(iNuclide),' [Bq]'
               else
                  write(*,'(I4,1X,A7,A,EN15.5,A)') NParticipatingNuclides,NuclideSpecs(iNuclide)%NuclideName,' chain:',&
                  &  MyStartCocktail%x(iNuclide),' [Bq]'
                  do jNuclide = 1,NNuclides
                     if ((MotherDaughterMatrix(jNuclide,iNuclide) /= 0._Float).AND.(jNuclide /= iNuclide)) then
                        write(*,'(11X,A,1X,A7)') '+ progeny',NuclideSpecs(jNuclide)%NuclideName
                     endif ! found progeny
                  enddo ! loop over possible progeny of head of chain
               endif ! progeny found
            endif ! participating
         endif ! nuclide is participating
      enddo ! loop over nuclides
      !
      ! option to Write pinpoint doses to file
      !
      if (DebugLevel > 0) then

         if (MotherNature == AnyNature) then
            NaturesName = ' '
         else if (MotherNature == NobleNature) then
            NaturesName = '_MNoble'
         else
            NaturesName = '_MNonNoble'
         endif

         if (DaughterNature == AnyNature) then
            NaturesName = trim(NaturesName)//' '
         else if (DaughterNature == NobleNature) then
            NaturesName = trim(NaturesName)//'_DNoble'
         else
            NaturesName = trim(NaturesName)//'_DNonNoble'
         endif

         TotalAirDoseRate = SUM(PinAirDoseRate,DIM=2)

         FName = trim(MyStartCocktail%MyDirectory)//trim(MyStartCocktail%MyName)//'_PinpointAirDoseRates_'&
         & //trim(WithDaughterString)//trim(NaturesName)//'.txt'

         open(ScratchFile,file = FName,form='FORMATTED',action='WRITE')
         write(ScratchFile,'(A)') 'Values are in (Sv/s)*m3. Multiplication with the air thinning '&
         & //'factor in /m3 gives the submersion dose rate in Sv/s.'
         write(ScratchFile,'(A,5X,F10.1,5X,A,F15.5,5X,A,5X,I0)') 'FirstDelay:',FirstDelay,&
         & 'DelayGrowthFactor:',DelayGrowthFactor,'NStartingTimes:',NStartingTimes
         write(ScratchFile,'(A,2000(A15,1X))') 'Pinpoint      t[s]           SumDoseRate ',&
         & (NuclideSpecs(ParticipatingNuclide(iNuclide))%NuclideName,iNuclide = 1,NParticipatingNuclides)
         do iPinpoint = 0,NStartingTimes
            write(ScratchFile,'(I8,1X,2000(G15.5,1X))') iPinpoint,AvailableDelay(iPinpoint),&
            & TotalAirDoseRate(iPinpoint),&
            & (PinAirDoseRate(iPinpoint,ParticipatingNuclide(iNuclide)),iNuclide = 1,NParticipatingNuclides)
         enddo ! loop over nuclides
         close(ScratchFile)

         TotalGroundDoseRate = SUM(PinGroundDoseRate,DIM=2)
         FName = trim(MyStartCocktail%MyDirectory)//trim(MyStartCocktail%MyName)//'_PinpointGroundDoseRates_'&
         & //trim(WithDaughterString)//trim(NaturesName)//'.txt'
         open(ScratchFile,file = FName,form='FORMATTED',action='WRITE')
         write(ScratchFile,'(A)') 'Values are in (Sv/s)*m2. Multiplication with the deposition thinning '&
         & //'factor in /m2 gives the ground dose rate in Sv/s.'
         write(ScratchFile,'(A,5X,F10.1,5X,A,F15.5,5X,A,5X,I0)') 'FirstDelay:',FirstDelay,&
         & 'DelayGrowthFactor:',DelayGrowthFactor,'NStartingTimes:',NStartingTimes
         write(ScratchFile,'(A,2000(A15,1X))') 'Pinpoint      t[s]          SumDoseRate ',&
         & (NuclideSpecs(ParticipatingNuclide(iNuclide))%NuclideName,iNuclide = 1,NParticipatingNuclides)
         do iPinpoint = 0,NStartingTimes
            write(ScratchFile,'(I8,1X,2000(G15.5,1X))') iPinpoint,AvailableDelay(iPinpoint),&
            & TotalGroundDoseRate(iPinpoint),&
            & (PinGroundDoseRate(iPinpoint,ParticipatingNuclide(iNuclide)),iNuclide = 1,NParticipatingNuclides)
         enddo ! loop over nuclides
         close(ScratchFile)

         TotalImmersionDoseRate = SUM(PinInhalationDoseRate,DIM=2)
         FName = trim(MyStartCocktail%MyDirectory)//trim(MyStartCocktail%MyName)//'_PinpointInhalationDoseRates_'&
         & //trim(WithDaughterString)//trim(NaturesName)//'.txt'
         open(ScratchFile,file = FName,form='FORMATTED',action='WRITE')
         write(ScratchFile,'(A)') 'Values are in Bq * (Sv/Bq) * (h/s). Multiplication with the air thinning '&
         & //'factor in /m3 and the breathing rate in m3/h gives the inhalation dose rate in Sv/s.'
         write(ScratchFile,'(A,5X,F10.1,5X,A,F15.5,5X,A,5X,I0)') 'FirstDelay:',FirstDelay,&
         & 'DelayGrowthFactor:',DelayGrowthFactor,'NStartingTimes:',NStartingTimes
         write(ScratchFile,'(A,2000(A15,1X))') 'Pinpoint      t[s]          SumDoseRate ',&
         & (NuclideSpecs(ParticipatingNuclide(iNuclide))%NuclideName,iNuclide = 1,NParticipatingNuclides)
         do iPinpoint = 0,NStartingTimes
            write(ScratchFile,'(I8,1X,2000(G15.5,1X))') iPinpoint,AvailableDelay(iPinpoint),&
            & TotalImmersionDoseRate(iPinpoint),&
            & (PinInhalationDoseRate(iPinpoint,ParticipatingNuclide(iNuclide)),iNuclide = 1,NParticipatingNuclides)
         enddo ! loop over nuclides
         close(ScratchFile)

         TotalThyroidDoseRate = SUM(PinThyroidInhalationDoseRate,DIM=2)
         FName = trim(MyStartCocktail%MyDirectory)//trim(MyStartCocktail%MyName)//'_PinpointThyroidInhalationDoseRates_'&
         & //trim(WithDaughterString)//trim(NaturesName)//'.txt'
         open(ScratchFile,file = FName,form='FORMATTED',action='WRITE')
         write(ScratchFile,'(A)') 'Values are in Bq * (Sv/Bq) * (h/s). Multiplication with the air thinning '&
         & //'factor in /m3 and the breathing rate in m3/h gives the inhalation dose rate in Sv/s.'
         write(ScratchFile,'(A,5X,F10.1,5X,A,F15.5,5X,A,5X,I0)') 'FirstDelay:',FirstDelay,&
         & 'DelayGrowthFactor:',DelayGrowthFactor,'NStartingTimes:',NStartingTimes
         write(ScratchFile,'(A,2000(A15,1X))') 'Pinpoint      t[s]          SumDoseRate ',&
         & (NuclideSpecs(ParticipatingNuclide(iNuclide))%NuclideName,iNuclide = 1,NParticipatingNuclides)
         do iPinpoint = 0,NStartingTimes
            write(ScratchFile,'(I8,1X,2000(G15.5,1X))') iPinpoint,AvailableDelay(iPinpoint),&
            & TotalThyroidDoseRate(iPinpoint),&
            & (PinThyroidInhalationDoseRate(iPinpoint,ParticipatingNuclide(iNuclide)),iNuclide = 1,NParticipatingNuclides)
         enddo ! loop over nuclides
         close(ScratchFile)
      endif
   end subroutine MakePinpointDoseRates



   subroutine GetPinpointCocktails(MyStartCocktail,WithProgeny)
      !
      ! Construct the decayed cocktail for all pinpoints
      !
      type(CocktailType), intent(in) :: MyStartCocktail
      logical, intent(in) :: WithProgeny

      type(CocktailType), dimension(0:NStartingTimes) :: MyCocktail
      character(DefaultLength) :: FName,WithDaughterString
      integer :: iNuclide,iPinpoint,NParticipatingNuclides
      integer, dimension(MaxNuclides) :: ParticipatingNuclide

      if (WithProgeny) then
         WithDaughterString = 'withprogeny'
      else
         WithDaughterString = 'noprogeny'
      endif
      !
      ! Mature the starting cocktail for all pinpoints
      !
      MyCocktail(0) = MyStartCocktail
      do iPinpoint = 1,NStartingTimes
         call MatureCocktail(MyStartCocktail,TransitionMatrix(iPinpoint),WithProgeny,MyCocktail(iPinpoint))
      enddo ! loop over all pinpoints
      !
      ! Label the contributing nuclides (to reduce output)
      !
      NParticipatingNuclides = 0
      do iNuclide = 1,NNuclides
         if (ANY(MyCocktail%x(iNuclide) > 0._Float)) then
            NParticipatingNuclides = NParticipatingNuclides + 1
            ParticipatingNuclide(NParticipatingNuclides) = iNuclide
         endif ! nuclide is participating
      enddo ! loop over nuclides
      !
      ! Write result to file
      !
      FName = trim(MyStartCocktail%MyDirectory)//trim(MyStartCocktail%MyName)//'_PinpointCocktail_'&
      & //trim(WithDaughterString)//'.txt'

      open(ScratchFile,file=FName,form='FORMATTED',action='WRITE')

      write(ScratchFile,'(A,5X,F10.1,5X,A,F15.5,5X,A,5X,I0)') 'FirstDelay:',FirstDelay,&
         & 'DelayGrowthFactor:',DelayGrowthFactor,'NStartingTimes:',NStartingTimes
      write(ScratchFile,'(I0,A)') NParticipatingNuclides,' = number of nuclides in cocktail'
      write(ScratchFile,'(A,2000(A15,1X))') 'Pinpoint      t[s]        ',&
      & (NuclideSpecs(ParticipatingNuclide(iNuclide))%NuclideName,iNuclide = 1,NParticipatingNuclides)
      do iPinpoint = 0,NStartingTimes
         write(ScratchFile,'(I8,1X,2000(G15.5,1X))') iPinpoint,AvailableDelay(iPinpoint),&
            & (MyCocktail(iPinpoint)%x(ParticipatingNuclide(iNuclide)),iNuclide = 1,NParticipatingNuclides)
      enddo ! loop over all pinpoints

      close(ScratchFile)
   end subroutine GetPinpointCocktails
end module Libpinpoint
