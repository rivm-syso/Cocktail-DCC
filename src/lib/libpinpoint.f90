MODULE LibPinPoint
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
   USE libxmath
   USE libutil
   USE libdcc
   USE libendf

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CocktailType,InitLibPinPoint,NStartingTimes,StartingTimeName,&
   & AvailableDelay,TransitionMatrix,GetPinpointCocktails,&
   & MakePinpointDoseRates,ReadRIVMSourceTermFile,SaveCocktail2,FirstDelay,DelayGrowthFactor,MatureCocktail

   TYPE CocktailType
      CHARACTER(DefaultLength) :: MyName,MyDirectory
      REAL(Float), DIMENSION(MaxNuclides) :: x
   END TYPE CocktailType

   CHARACTER(DefaultLength) :: ProjectPath = './'

   REAL(Float), PARAMETER :: FirstDelay = 60._Float ! seconds
   REAL(Float), PARAMETER :: DelayGrowthFactor = 1.15_Float
   INTEGER, PARAMETER :: NStartingTimes = 200 ! There is also pinpoint 0, which is 0

   CHARACTER(15), DIMENSION(0:NStartingTimes) :: StartingTimeName
   REAL(Float), DIMENSION(0:NStartingTimes) :: AvailableDelay ! In seconds

   TYPE(SparseMatrix), DIMENSION(0:NStartingTimes) :: TransitionMatrix

CONTAINS
   SUBROUTINE InitLibPinPoint(UseICRP)
      !
      ! Prepare some settings before use of this library
      !
      LOGICAL, INTENT(IN) :: UseICRP

      INTEGER :: iStartingTime
      REAL(Float) :: tMin
      TYPE(NuclideType), DIMENSION(0:MaxNuclides) :: RegularizedNuclideSpecs
      REAL(Float), DIMENSION(MaxNuclides,MaxNuclides) :: RegularizedMotherDaughterMatrix

      INTEGER, PARAMETER :: DebugLevel = 0
      !
      ! Get the DCCs
      !
      CALL InitLibDCC(UseICRP) ! includes reading data from ICRP or ENDF
      CALL ReadTissueDCCs()
      CALL ReadDCCEffectiveDose()
      CALL ReadGroundDCCs()
      CALL ReadInhalationDCCs()
      CALL ReadThyroidInhalationDCCs()

      AvailableDelay(0) = 0
      AvailableDelay(1) = 60
      DO iStartingTime = 2,NStartingTimes
         AvailableDelay(iStartingTime) = DelayGrowthFactor*AvailableDelay(iStartingTime-1)
      ENDDO ! loop over starting times

      IF (DebugLevel.GT.0) WRITE(*,*)
      IF (DebugLevel.GT.0) WRITE(*,'(A)') 'Considering the following exponentially distributed delays:'
      DO iStartingTime = 0,NStartingTimes
         WRITE(StartingTimeName(iStartingTime),'(E11.5)') AvailableDelay(iStartingTime)
         IF (DebugLevel.GT.0) WRITE(*,'(A)') StartingTimeName(iStartingTime)//'second'
      ENDDO ! loop over starting times

      CALL GetTransitionMatrices(UseICRP)
   END SUBROUTINE InitLibPinPoint



   SUBROUTINE ReadRIVMSourceTermFile(FName,MyCocktail)
      !
      ! Read sourceterm with RIVM file format
      !
      TYPE(CocktailType), INTENT(OUT) :: MyCocktail
      CHARACTER(*), INTENT(IN) :: FName

      CHARACTER(10) :: MyName
      CHARACTER(DefaultLength) :: ALine
      INTEGER :: iLine,iNuclide,TheIndex,iCharacter
      LOGICAL :: Ready,IsNotABackSlash
      REAL(Float) :: MyActivity,ScalingFactor

      INTEGER, PARAMETER :: DebugLevel = 0
      !
      ! Extract name of the source term from filename
      !
      TheIndex = INDEX(FName,'.RIVMSource')
      IF (TheIndex.EQ.0) THEN
         WRITE(*,'(A)') 'Sourcterm "'//TRIM(FName)//'" is not of type .RIVMSource! Exiting!!'
         CALL EXIT()
      ENDIF
      !
      ! Extract directory
      !
      MyCocktail%MyDirectory = FName
      iCharacter = LEN_TRIM(MyCocktail%MyDirectory)
      IsNotABackSlash = .NOT.(MyCocktail%MyDirectory(iCharacter:iCharacter).EQ.'/')
      DO WHILE (IsNotABackSlash.AND.(iCharacter.GT.0))
         MyCocktail%MyDirectory(iCharacter:iCharacter) = ' '
         iCharacter = iCharacter - 1
         IF (iCharacter.GT.0) IsNotABackSlash = .NOT.(MyCocktail%MyDirectory(iCharacter:iCharacter).EQ.'/')
      ENDDO
      !
      ! Extract name of sourceterm
      !
      MyCocktail%MyName = FName((LEN_TRIM(MyCocktail%MyDirectory)+1):(TheIndex-1))
      !
      ! Read file
      !
      MyCocktail%x = 0._Float

      OPEN(ScratchFile,FILE=FName,FORM='FORMATTED',ACTION='READ')
      !
      ! Read header to find meta-info, e.g. scaling factor for activity
      !
      ScalingFactor = 1._Float

      Ready = .FALSE.
      DO WHILE (.NOT.Ready)
         READ(ScratchFile,'(A)') ALine
         IF (DebugLevel.GT.1) WRITE(*,'(A)') 'ALine = "'//TRIM(ALine)//'"'

         TheIndex = INDEX(ALine,'<scalingfactor>')
         IF (TheIndex.GT.0) THEN
            READ(ALine((TheIndex+15):LEN_TRIM(ALine)),*) ScalingFactor
            IF (DebugLevel.GT.1) WRITE(*,'(A,EN15.5)') 'Found scaling factor ',ScalingFactor
         ENDIF

         Ready = (ALine(1:1).NE.'!')
      ENDDO
      !
      ! Read the nuclide vector and take into account the optional scaling factor
      !
      WRITE(*,'(A)') 'Reading source term at t=0 from '//TRIM(FName)

      Ready = .FALSE.

      DO WHILE (.NOT.Ready)
         IF (DebugLevel.GT.1) WRITE(*,'(A)') 'ALine = "'//TRIM(ALine)//'"'

         IF (ALine(1:1).NE.'!') THEN
            READ(ALine,*) MyName,MyActivity

            CALL MassNuc2NucMass(MyName)
            CALL AllUpCase(MyName)
            CALL EnsureHyphen(MyName)

            iNuclide = GetNuclideNumber(MyName)

            IF (iNuclide.EQ.0) THEN
               WRITE(*,'(A)') 'Could not recognize nuclide "'//TRIM(MyName)//'"! Exiting!'
               CALL EXIT()
            ENDIF
            MyCocktail%x(iNuclide) = ScalingFactor * MyActivity ! [Bq/10kt]
            IF (DebugLevel.GT.0) THEN
               WRITE(*,'(A,A,I4,A,EN20.10,A)') MyName,' = nuclide ',iNuclide,&
               & ' and has activity ',MyCocktail%x(iNuclide),' Bq/10kt'
            ENDIF
         ENDIF ! not a comment line stat starts with !

         READ(ScratchFile,'(A)',END=10) ALine

         Ready = (LEN_TRIM(ALine).EQ.0)
      ENDDO ! loop until ready

      10 CONTINUE
      CLOSE(ScratchFile)
   END SUBROUTINE ReadRIVMSourceTermFile



   SUBROUTINE GetTransitionMatrices(UseICRP)
      !
      ! Get transition matrices for all pinpoints
      !
      LOGICAL, INTENT(IN) :: UseICRP

      CHARACTER(DefaultLength) :: FName,ALine
      INTEGER :: iStartingTime,iLine
      CHARACTER(6) :: PinPointName

      INTEGER, PARAMETER :: DebugLevel = 0

      WRITE(*,'(A)') 'Getting transition matrices...'
      DO iStartingTime = 1,NStartingTimes
         !
         ! Open file
         !
         WRITE(PinPointName,'(I3.3)') iStartingTime

         IF (UseICRP) THEN
            FName = 'SparseMatrix_'//TRIM(StartingTimeName(iStartingTime))//'second.dat'
         ELSE
            FName = 'SparseMatrix_ENDF_'//TRIM(PinPointName)//'.dat'
         ENDIF

         IF (DebugLevel.GT.0) WRITE(*,'(A)') 'Reading '//TRIM(FName)
         CALL FLUSH(6)
         OPEN(ScratchFile, FILE=TransitionMatrixPath() // '/' // TRIM(FName), ACTION='READ')

         IF (.NOT.UseICRP) THEN
            DO iLine=1,7
            READ(ScratchFile,*) ! Skip header
            ENDDO
         ENDIF
         !
         ! Count number of sparse matrix elements to allocate space for
         !
         TransitionMatrix(iStartingTime)%N = 0
         DO
            READ(ScratchFile,'(A)',END=10) ALine
            TransitionMatrix(iStartingTime)%N = TransitionMatrix(iStartingTime)%N + 1
         ENDDO
         10  CONTINUE

         IF (DebugLevel.GT.0) WRITE(*,'(A,I0,A)') TRIM(FName)//' has ',TransitionMatrix(iStartingTime)%N,' elements'

         REWIND(ScratchFile)
         !
         ! Rewind file and read matrix elements
         !
         TransitionMatrix(iStartingTime)%NMax = TransitionMatrix(iStartingTime)%N
         ALLOCATE(TransitionMatrix(iStartingTime)%Element(TransitionMatrix(iStartingTime)%NMax))

         TransitionMatrix(iStartingTime)%N = 0

         IF (.NOT.UseICRP) THEN
            DO iLine=1,7
            READ(ScratchFile,*) ! Skip header
            ENDDO
         ENDIF

         DO
            READ(ScratchFile,'(A)',END=20) ALine
            TransitionMatrix(iStartingTime)%N = TransitionMatrix(iStartingTime)%N + 1
            READ(ALine,*) TransitionMatrix(iStartingTime)%Element(TransitionMatrix(iStartingTime)%N)
         ENDDO
         20  CONTINUE
         CLOSE(ScratchFile)
      ENDDO ! loop over starting times
   END SUBROUTINE GetTransitionMatrices



   SUBROUTINE SaveCocktail2(Cocktail,Tag)
      !
      ! Save a cocktail of radionuclides to file
      !
      TYPE(CocktailType), INTENT(IN) :: Cocktail
      CHARACTER(*), INTENT(IN) :: Tag

      REAL(Float), PARAMETER :: NegligibleActivity = 0._Float ! 10.E9_Float ! Bq/10kt, seems to characterize difference between tables 2 and 5
      CHARACTER(DefaultLength) :: FName
      INTEGER :: iNuclide

      FName = TRIM(Cocktail%MyName)//'_'//TRIM(Tag)//'.dat'
      OPEN(ScratchFile,FILE=TRIM(Cocktail%MyDirectory)//TRIM(FName),ACTION='WRITE')
      WRITE(ScratchFile,'(A)') 'Nuclide    Activity[Bq]'
      DO iNuclide = 1,NNuclides
         IF (Cocktail%x(iNuclide).GT.NegligibleActivity) THEN
            WRITE(ScratchFile,'(A,1X,EN11.2)') NuclideSpecs(iNuclide)%NuclideName,Cocktail%x(iNuclide)
         ENDIF
      ENDDO
      CLOSE(ScratchFile)
   END SUBROUTINE SaveCocktail2




   SUBROUTINE MatureCocktail(InCocktail,TheMatrix,WithProgeny,OutCocktail)
      !
      ! Estimate the cocktail at a given delay from the cocktail at t=0 and the transition matrix for the given delay
      !
      TYPE(CocktailType), INTENT(IN) :: InCocktail
      TYPE(SparseMatrix), INTENT(IN) :: TheMatrix
      LOGICAL, INTENT(IN) :: WithProgeny
      TYPE(CocktailType), INTENT(OUT) :: OutCocktail

      INTEGER :: iElement,iMother,iDaughter,iNuclide
      REAL(Float) :: x

      INTEGER, PARAMETER :: DebugLevel = 0

      IF (DebugLevel.GT.0) THEN
         DO iNuclide = 1,NNuclides
            IF (InCocktail%x(iNuclide).GT.0._Float) THEN
               WRITE(*,'(A,A,A,I4,A,EN20.10,A)') 'MatureCocktail: ',&
            & NuclideSpecs(iNuclide)%NuclideName,' = nuclide ',iNuclide,&
               &   ' and has activity ',InCocktail%x(iNuclide),' Bq'
            ENDIF
         ENDDO
      ENDIF ! debug

      OutCocktail%x = 0._Float
      OutCocktail%MyName = TRIM(InCocktail%MyName)//'_aged'
      OutCocktail%MyDirectory = InCocktail%MyDirectory

      DO iElement = 1,TheMatrix%N
         iMother = TheMatrix%Element(iElement)%i
         iDaughter = TheMatrix%Element(iElement)%j
         x = TheMatrix%Element(iElement)%x
         IF (WithProgeny.OR.(iMother.EQ.iDaughter)) THEN

            IF ((DebugLevel.GT.0).AND.(InCocktail%x(iMother).GT.0._Float)) THEN
               WRITE(*,'(A,A,A,A,A,EN15.5,A,EN15.5,A,EN15.5)') 'Maturization of ',&
               & NuclideSpecs(iDaughter)%NuclideName,&
               & ' because of initial presence of ',&
               & NuclideSpecs(iMother)%NuclideName ,&
               & ' : ',&
               & OutCocktail%x(iDaughter),&
               & ' + ',x*InCocktail%x(iMother),' = ',&
               & OutCocktail%x(iDaughter) + x*InCocktail%x(iMother)
            ENDIF ! debug

            OutCocktail%x(iDaughter) = OutCocktail%x(iDaughter) + x*InCocktail%x(iMother)
            IF (OutCocktail%x(iDaughter).LT.1.E-95_Float) OutCocktail%x(iDaughter) = 0._Float ! to prevent exponents with 3 digits, which print ugly...
         ENDIF
      ENDDO ! loop over elements
   END SUBROUTINE MatureCocktail



   SUBROUTINE MakePinpointDoseRates(MyStartCocktail,WithDaughters,MotherNature,DaughterNature)
      !
      ! Estimate dose contributions at all pinpoints
      ! Add dose contribution to the head of chain, assume ground activity to be fully on top = "planar"
      !
      TYPE(CocktailType), INTENT(IN) :: MyStartCocktail
      LOGICAL, INTENT(IN) :: WithDaughters
      INTEGER, INTENT(IN) :: MotherNature,DaughterNature ! 0 = any, 1 = noble gas, 2 = not a noble gas

      REAL(Float), DIMENSION(0:NStartingTimes,MaxNuclides) :: PinAirDoseRate,PinGroundDoseRate,&
      & PinInhalationDoseRate,PinThyroidInhalationDoseRate
      INTEGER, DIMENSION(MaxNuclides) :: ParticipatingNuclide
      INTEGER :: NParticipatingNuclides
      REAL(Float), DIMENSION(MaxNuclides) :: SumAirDoseRate
      INTEGER :: iPinpoint,iNuclide,iDaughter,iMother,iElement,jNuclide,NProgeny
      REAL(Float) :: DaughterActivity,x
      REAL(Float), DIMENSION(0:NStartingTimes) :: TotalAirDoseRate,TotalGroundDoseRate,TotalImmersionDoseRate,&
      & TotalThyroidDoseRate
      CHARACTER(DefaultLength) :: FName,WithDaughterString,NaturesName
      LOGICAL :: NaturesOkay,MotherNatureOkay,DaughterNatureOkay

      INTEGER, PARAMETER :: AnyNature = 0
      INTEGER, PARAMETER :: NobleNature = 1
      INTEGER, PARAMETER :: NonNobleNature = 2

      INTEGER, PARAMETER :: DebugLevel = 1

      PinAirDoseRate = 0._Float
      PinGroundDoseRate = 0._Float
      PinInhalationDoseRate = 0._Float
      PinThyroidInhalationDoseRate = 0._Float
      SumAirDoseRate = 0._Float

      IF (WithDaughters) THEN
         WithDaughterString = 'withprogeny'
      ELSE
         WithDaughterString = 'noprogeny'
      ENDIF

      DO iPinpoint = 0,NStartingTimes

         IF (iPinpoint.EQ.0) THEN
            ! According to pinpoint 0, there is only the head of chain, which is any fission product available at t=0.
            DO iNuclide = 1,NNuclides

               MotherNatureOkay =      (MotherNature.EQ.AnyNature)&
               &                  .OR.((MotherNature.EQ.NobleNature   ).AND.(NuclideSpecs(iNuclide)%NuclideGroup.EQ.1)) &
               &                  .OR.((MotherNature.EQ.NonNobleNature).AND.(NuclideSpecs(iNuclide)%NuclideGroup.NE.1))
               DaughterNatureOkay =    (DaughterNature.EQ.AnyNature)&
               &                  .OR.((DaughterNature.EQ.NobleNature   ).AND.(NuclideSpecs(iNuclide)%NuclideGroup.EQ.1)) &
               &                  .OR.((DaughterNature.EQ.NonNobleNature).AND.(NuclideSpecs(iNuclide)%NuclideGroup.NE.1))

               NaturesOkay = MotherNatureOkay .AND. DaughterNatureOkay

               IF (NaturesOkay) THEN
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
               ENDIF ! natures okay
            ENDDO ! loop over nuclides
         ELSE
            !
            ! For any pinpoint > 0, you have the full decay chain with progeny under each head of chain (formed at t=0)
            !
            DO iElement = 1,TransitionMatrix(iPinpoint)%N
               iMother = TransitionMatrix(iPinpoint)%Element(iElement)%i
               iDaughter = TransitionMatrix(iPinpoint)%Element(iElement)%j

               MotherNatureOkay =      (MotherNature.EQ.AnyNature)&
               &                  .OR.((MotherNature.EQ.NobleNature   ).AND.(NuclideSpecs(iMother)%NuclideGroup.EQ.1)) &
               &                  .OR.((MotherNature.EQ.NonNobleNature).AND.(NuclideSpecs(iMother)%NuclideGroup.NE.1))
               DaughterNatureOkay =    (DaughterNature.EQ.AnyNature)&
               &                  .OR.((DaughterNature.EQ.NobleNature   ).AND.(NuclideSpecs(iDaughter)%NuclideGroup.EQ.1)) &
               &                  .OR.((DaughterNature.EQ.NonNobleNature).AND.(NuclideSpecs(iDaughter)%NuclideGroup.NE.1))

               NaturesOkay = MotherNatureOkay .AND. DaughterNatureOkay

               IF (NaturesOkay) THEN

                  IF (WithDaughters.OR.(iMother.EQ.iDaughter)) THEN
                     x = TransitionMatrix(iPinpoint)%Element(iElement)%x

                     DaughterActivity = x * MyStartCocktail%x(iMother)

                     PinAirDoseRate(iPinpoint,iMother) = PinAirDoseRate(iPinpoint,iMother) &
                     & +   DaughterActivity &
                     &   * DCCEffectiveDose(InAir,Age_Adult)%x(iDaughter) &
                     &   * 1.E-9/3600._Float ! Conversion of nSv/h to Sv/s

                     PinGroundDoseRate(iPinpoint,iMother) = PinGroundDoseRate(iPinpoint,iMother) &
                     & +   DaughterActivity &
                     &   * DCCGround(1,Age_Adult)%x(iDaughter) &
                     &   * 1.E-9/3600._Float ! Conversion of nSv/h to Sv/s

                     PinInhalationDoseRate(iPinpoint,iMother) = PinInhalationDoseRate(iPinpoint,iMother) &
                     & +   DaughterActivity &
                     &   * DCCInhalation(Inhalation_public_adult)%x(iDaughter) &
                     &   * (1._Float/3600._Float) ! Conversion of Sv/(m3/s) to Sv/(m3/h) to facilitate breathing rates in m3/h

                     PinThyroidInhalationDoseRate(iPinpoint,iMother) = PinThyroidInhalationDoseRate(iPinpoint,iMother) &
                     & +   DaughterActivity &
                     &   * DCCThyroidInhalation%x(iDaughter) &
                     &   * (1._Float/3600._Float) ! Conversion of Sv/(m3/s) to Sv/(m3/h) to facilitate breathing rates in m3/h

                     SumAirDoseRate(iMother) = SumAirDoseRate(iMother) + PinAirDoseRate(iPinpoint,iMother)
                  ENDIF ! withdaughters
               ENDIF ! natures okay

            ENDDO ! loop over sparse matrix at iLower
         ENDIF ! first pinpoint is t=0

      ENDDO ! loop over all pinpoints
      !
      ! Cut off values < 0.1E-99, as they are printed in an ugly way...
      !
      DO iPinpoint = 0,NStartingTimes
         DO iNuclide = 1,NNuclides
            IF (PinAirDoseRate(iPinpoint,iNuclide).LE.0.1E-99_Float) &
            & PinAirDoseRate(iPinpoint,iNuclide) = 0._Float
            IF (PinGroundDoseRate(iPinpoint,iNuclide).LE.0.1E-99_Float) &
            & PinGroundDoseRate(iPinpoint,iNuclide) = 0._Float
            IF (PinInhalationDoseRate(iPinpoint,iNuclide).LE.0.1E-99_Float) &
            & PinInhalationDoseRate(iPinpoint,iNuclide) = 0._Float
            IF (PinThyroidInhalationDoseRate(iPinpoint,iNuclide).LE.0.1E-99_Float) &
            & PinThyroidInhalationDoseRate(iPinpoint,iNuclide) = 0._Float
         ENDDO ! loop over nuclides
      ENDDO ! loop over all pinpoints
      !
      ! Check which nuclides contribute at all on basis of sum of air doserates.
      ! Ground and inhalation rates are not considered to give a different answer.
      !
      NParticipatingNuclides = 0
      DO iNuclide = 1,NNuclides
         IF (SumAirDoseRate(iNuclide).GT.0._Float) THEN
            NParticipatingNuclides = NParticipatingNuclides + 1
            ParticipatingNuclide(NParticipatingNuclides) = iNuclide
            IF (DebugLevel.GT.0) THEN
               NProgeny = 0
               DO jNuclide = 1,NNuclides
                  IF ((MotherDaughterMatrix(jNuclide,iNuclide).NE.0._Float).AND.(jNuclide.NE.iNuclide)) THEN
                     NProgeny = NProgeny + 1
                  ENDIF ! found progeny
               ENDDO ! loop over possible progeny of head of chain
               IF (NProgeny.EQ.0) THEN
                  WRITE(*,'(I4,1X,A7,7X,EN15.5,A)') NParticipatingNuclides,NuclideSpecs(iNuclide)%NuclideName,&
                  & MyStartCocktail%x(iNuclide),' [Bq]'
               ELSE
                  WRITE(*,'(I4,1X,A7,A,EN15.5,A)') NParticipatingNuclides,NuclideSpecs(iNuclide)%NuclideName,' chain:',&
                  &  MyStartCocktail%x(iNuclide),' [Bq]'
                  DO jNuclide = 1,NNuclides
                     IF ((MotherDaughterMatrix(jNuclide,iNuclide).NE.0._Float).AND.(jNuclide.NE.iNuclide)) THEN
                        WRITE(*,'(11X,A,1X,A7)') '+ progeny',NuclideSpecs(jNuclide)%NuclideName
                     ENDIF ! found progeny
                  ENDDO ! loop over possible progeny of head of chain
               ENDIF ! progeny found
            ENDIF ! participating
         ENDIF ! nuclide is participating
      ENDDO ! loop over nuclides
      !
      ! option to Write pinpoint doses to file
      !
      IF (DebugLevel.GT.0) THEN

         IF (MotherNature.EQ.AnyNature) THEN
            NaturesName = ' '
         ELSE IF (MotherNature.EQ.NobleNature) THEN
            NaturesName = '_MNoble'
         ELSE
            NaturesName = '_MNonNoble'
         ENDIF

         IF (DaughterNature.EQ.AnyNature) THEN
            NaturesName = TRIM(NaturesName)//' '
         ELSE IF (DaughterNature.EQ.NobleNature) THEN
            NaturesName = TRIM(NaturesName)//'_DNoble'
         ELSE
            NaturesName = TRIM(NaturesName)//'_DNonNoble'
         ENDIF

         TotalAirDoseRate = SUM(PinAirDoseRate,DIM=2)

         FName = TRIM(MyStartCocktail%MyDirectory)//TRIM(MyStartCocktail%MyName)//'_PinpointAirDoseRates_'&
         & //TRIM(WithDaughterString)//TRIM(NaturesName)//'.txt'

         OPEN(ScratchFile,FILE = FName,FORM='FORMATTED',ACTION='WRITE')
         WRITE(ScratchFile,'(A)') 'Values are in (Sv/s)*m3. Multiplication with the air thinning '&
         & //'factor in /m3 gives the submersion dose rate in Sv/s.'
         WRITE(ScratchFile,'(A,5X,F10.1,5X,A,F15.5,5X,A,5X,I0)') 'FirstDelay:',FirstDelay,&
         & 'DelayGrowthFactor:',DelayGrowthFactor,'NStartingTimes:',NStartingTimes
         WRITE(ScratchFile,'(A,2000(A15,1X))') 'Pinpoint      t[s]           SumDoseRate ',&
         & (NuclideSpecs(ParticipatingNuclide(iNuclide))%NuclideName,iNuclide = 1,NParticipatingNuclides)
         DO iPinpoint = 0,NStartingTimes
            WRITE(ScratchFile,'(I8,1X,2000(G15.5,1X))') iPinpoint,AvailableDelay(iPinpoint),&
            & TotalAirDoseRate(iPinpoint),&
            & (PinAirDoseRate(iPinpoint,ParticipatingNuclide(iNuclide)),iNuclide = 1,NParticipatingNuclides)
         ENDDO ! loop over nuclides
         CLOSE(ScratchFile)

         TotalGroundDoseRate = SUM(PinGroundDoseRate,DIM=2)
         FName = TRIM(MyStartCocktail%MyDirectory)//TRIM(MyStartCocktail%MyName)//'_PinpointGroundDoseRates_'&
         & //TRIM(WithDaughterString)//TRIM(NaturesName)//'.txt'
         OPEN(ScratchFile,FILE = FName,FORM='FORMATTED',ACTION='WRITE')
         WRITE(ScratchFile,'(A)') 'Values are in (Sv/s)*m2. Multiplication with the deposition thinning '&
         & //'factor in /m2 gives the ground dose rate in Sv/s.'
         WRITE(ScratchFile,'(A,5X,F10.1,5X,A,F15.5,5X,A,5X,I0)') 'FirstDelay:',FirstDelay,&
         & 'DelayGrowthFactor:',DelayGrowthFactor,'NStartingTimes:',NStartingTimes
         WRITE(ScratchFile,'(A,2000(A15,1X))') 'Pinpoint      t[s]          SumDoseRate ',&
         & (NuclideSpecs(ParticipatingNuclide(iNuclide))%NuclideName,iNuclide = 1,NParticipatingNuclides)
         DO iPinpoint = 0,NStartingTimes
            WRITE(ScratchFile,'(I8,1X,2000(G15.5,1X))') iPinpoint,AvailableDelay(iPinpoint),&
            & TotalGroundDoseRate(iPinpoint),&
            & (PinGroundDoseRate(iPinpoint,ParticipatingNuclide(iNuclide)),iNuclide = 1,NParticipatingNuclides)
         ENDDO ! loop over nuclides
         CLOSE(ScratchFile)

         TotalImmersionDoseRate = SUM(PinInhalationDoseRate,DIM=2)
         FName = TRIM(MyStartCocktail%MyDirectory)//TRIM(MyStartCocktail%MyName)//'_PinpointInhalationDoseRates_'&
         & //TRIM(WithDaughterString)//TRIM(NaturesName)//'.txt'
         OPEN(ScratchFile,FILE = FName,FORM='FORMATTED',ACTION='WRITE')
         WRITE(ScratchFile,'(A)') 'Values are in Bq * (Sv/Bq) * (h/s). Multiplication with the air thinning '&
         & //'factor in /m3 and the breathing rate in m3/h gives the inhalation dose rate in Sv/s.'
         WRITE(ScratchFile,'(A,5X,F10.1,5X,A,F15.5,5X,A,5X,I0)') 'FirstDelay:',FirstDelay,&
         & 'DelayGrowthFactor:',DelayGrowthFactor,'NStartingTimes:',NStartingTimes
         WRITE(ScratchFile,'(A,2000(A15,1X))') 'Pinpoint      t[s]          SumDoseRate ',&
         & (NuclideSpecs(ParticipatingNuclide(iNuclide))%NuclideName,iNuclide = 1,NParticipatingNuclides)
         DO iPinpoint = 0,NStartingTimes
            WRITE(ScratchFile,'(I8,1X,2000(G15.5,1X))') iPinpoint,AvailableDelay(iPinpoint),&
            & TotalImmersionDoseRate(iPinpoint),&
            & (PinInhalationDoseRate(iPinpoint,ParticipatingNuclide(iNuclide)),iNuclide = 1,NParticipatingNuclides)
         ENDDO ! loop over nuclides
         CLOSE(ScratchFile)

         TotalThyroidDoseRate = SUM(PinThyroidInhalationDoseRate,DIM=2)
         FName = TRIM(MyStartCocktail%MyDirectory)//TRIM(MyStartCocktail%MyName)//'_PinpointThyroidInhalationDoseRates_'&
         & //TRIM(WithDaughterString)//TRIM(NaturesName)//'.txt'
         OPEN(ScratchFile,FILE = FName,FORM='FORMATTED',ACTION='WRITE')
         WRITE(ScratchFile,'(A)') 'Values are in Bq * (Sv/Bq) * (h/s). Multiplication with the air thinning '&
         & //'factor in /m3 and the breathing rate in m3/h gives the inhalation dose rate in Sv/s.'
         WRITE(ScratchFile,'(A,5X,F10.1,5X,A,F15.5,5X,A,5X,I0)') 'FirstDelay:',FirstDelay,&
         & 'DelayGrowthFactor:',DelayGrowthFactor,'NStartingTimes:',NStartingTimes
         WRITE(ScratchFile,'(A,2000(A15,1X))') 'Pinpoint      t[s]          SumDoseRate ',&
         & (NuclideSpecs(ParticipatingNuclide(iNuclide))%NuclideName,iNuclide = 1,NParticipatingNuclides)
         DO iPinpoint = 0,NStartingTimes
            WRITE(ScratchFile,'(I8,1X,2000(G15.5,1X))') iPinpoint,AvailableDelay(iPinpoint),&
            & TotalThyroidDoseRate(iPinpoint),&
            & (PinThyroidInhalationDoseRate(iPinpoint,ParticipatingNuclide(iNuclide)),iNuclide = 1,NParticipatingNuclides)
         ENDDO ! loop over nuclides
         CLOSE(ScratchFile)
      ENDIF
   END SUBROUTINE MakePinpointDoseRates



   SUBROUTINE GetPinpointCocktails(MyStartCocktail,WithProgeny)
      !
      ! Construct the decayed cocktail for all pinpoints
      !
      TYPE(CocktailType), INTENT(IN) :: MyStartCocktail
      LOGICAL, INTENT(IN) :: WithProgeny

      TYPE(CocktailType), DIMENSION(0:NStartingTimes) :: MyCocktail
      CHARACTER(DefaultLength) :: FName,WithDaughterString
      INTEGER :: iNuclide,iPinpoint,NParticipatingNuclides
      INTEGER, DIMENSION(MaxNuclides) :: ParticipatingNuclide

      IF (WithProgeny) THEN
         WithDaughterString = 'withprogeny'
      ELSE
         WithDaughterString = 'noprogeny'
      ENDIF
      !
      ! Mature the starting cocktail for all pinpoints
      !
      MyCocktail(0) = MyStartCocktail
      DO iPinpoint = 1,NStartingTimes
         CALL MatureCocktail(MyStartCocktail,TransitionMatrix(iPinpoint),WithProgeny,MyCocktail(iPinpoint))
      ENDDO ! loop over all pinpoints
      !
      ! Label the contributing nuclides (to reduce output)
      !
      NParticipatingNuclides = 0
      DO iNuclide = 1,NNuclides
         IF (MyStartCocktail%x(iNuclide).GT.0._Float) THEN
            NParticipatingNuclides = NParticipatingNuclides + 1
            ParticipatingNuclide(NParticipatingNuclides) = iNuclide
         ENDIF ! nuclide is participating
      ENDDO ! loop over nuclides
      !
      ! Write result to file
      !
      FName = TRIM(MyStartCocktail%MyDirectory)//TRIM(MyStartCocktail%MyName)//'_PinpointCocktail_'&
      & //TRIM(WithDaughterString)//'.txt'

      OPEN(ScratchFile,FILE=FName,FORM='FORMATTED',ACTION='WRITE')

      WRITE(ScratchFile,'(A,5X,F10.1,5X,A,F15.5,5X,A,5X,I0)') 'FirstDelay:',FirstDelay,&
         & 'DelayGrowthFactor:',DelayGrowthFactor,'NStartingTimes:',NStartingTimes
      WRITE(ScratchFile,'(I0,A)') NParticipatingNuclides,' = number of nuclides in cocktail'
      WRITE(ScratchFile,'(A,2000(A15,1X))') 'Pinpoint      t[s]        ',&
      & (NuclideSpecs(ParticipatingNuclide(iNuclide))%NuclideName,iNuclide = 1,NParticipatingNuclides)
      DO iPinpoint = 0,NStartingTimes
         WRITE(ScratchFile,'(I8,1X,2000(G15.5,1X))') iPinpoint,AvailableDelay(iPinpoint),&
            & (MyCocktail(iPinpoint)%x(ParticipatingNuclide(iNuclide)),iNuclide = 1,NParticipatingNuclides)
      ENDDO ! loop over all pinpoints

      CLOSE(ScratchFile)
   END SUBROUTINE GetPinpointCocktails
END MODULE Libpinpoint
