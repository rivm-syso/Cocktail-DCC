MODULE LibCocktailDCC
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
   ! ____________________________________________________
   !
   ! This module gives access to two functions applicable to the situation after the detonation of
   ! one of the supported types of nuclear device:
   !
   ! One of the functions provides the compound cocktail dose conversion coefficient (DCC) for inhalation,
   ! submersion or groundshine at any given delay after the blast, without a need to follow individual nuclides.
   ! The advantage of this approach is that a lightweight atmospheric dispersion calculation,
   ! involving only 1 substance, suffices to assess the situation. Both the regular DCC and the cumulative
   ! (time integrated) DCC are supported. In a practical application, one should use an atmospheric dispersion
   ! model on a unit release to find the thinning factor or ground deposition fraction at the time/location
   ! of interest. This thinning factor should then be scaled with the yield of the blast, as the DCC given
   ! by libcocktaildcc is per 10kt. Multiplication of the yield-scaled thinning factor (for submersion or
   ! groundshine) with the cocktail DCC gives the dose rate in Sv/s. For inhalation, one also has to multiply
   ! with the inhalation rate in m3/h.
   !
   ! The other function gives the remaining amount [Bq] of any nuclide at any delay after blast, where
   ! all ingrowth of progeny is taken into account. Like the DCC function, to get a meaningful concentration
   ! estimate, this value has to be multiplied with the yield-scaled thinning factor from an atmospheric
   ! transport model. Otherwise, the function gives the total activity of a particular single nuclide, assuming
   ! the 10kt nuclide vector specified at t=0. This function can be used to find concentrations of individual nuclides,
   ! e.g. for the assessment of a necessity to commend iodine profylaxis, closure of greenhouses or putting
   ! in their pens of livestock.
   !
   !
   ! This module offers the following functionality:
   !
   ! You can specify the names of the files specifying your source terms in a text file named RIVMSources.txt.
   ! This file should be located in the COCKTAIL_DCC_SOURCES_DIR directory (environment variable).
   !
   ! The first lines that start with a "!" are considered to be comment lines and are therefore skipped.
   ! All subsequent lines in RIVMSources.txt are assumed to be references to text files containing source terms at t=0.
   ! These source term files must have names with extension ".RIVMSource", which means that
   ! they are 2-column ASCII files specifying names and activities of nuclides at t=0. As with file RIVMSources.txt,
   ! any line starting with '!' is considered to be a comment line and therefore skipped, with 1 exception:
   ! if the tag "<scalingfactor>" is found in a comment line, then the value that is given after it on the same line
   ! is used to scale the activities. In this way a header line (including the exclamation mark!):
   !
   ! <scalingfactor> 1.0E6
   !
   ! can be used to specify that the activities are given in MBq instead of Bq.
   !
   ! There is a utility named "source_term_dose.exe" that, for a given source term at t=0, can generate lookup tables
   ! for the DCCs and for the individual activities. For each source term, four ASCII format lookup tables are used:
   ! three for the different cocktail DCCs and for the nuclide specific concentration.
   ! These files are assumed to be located in the same directory as the source term file.
   ! Upon first use of module libcocktaildcc, the initialization routine assesses availability of the lookup
   ! tables for each source term in the file RIVMSources.txt and reads them for future use. If for any source term the
   ! associated lookup tables cannot be found, then a system call is made to utility source_term_dose.exe with the aim to
   ! add them. This takes about half a minute for 1 source term and needs to be done only once for each new source term.
   ! In this way, one can add support for a new source term by the construction of a correct .RIVMSource file
   ! and specification of its name and location on a new line in file RIVMSources.txt. The lookup tables will be
   ! added automatically.
   !
   ! When you are working with module libcocktaildcc, a source terms in your personal file RIVMSources.txt
   ! is referenced through an index equal to the rank of the source term in the list.
   ! Challenge: find a way to comfortably refer to the different source terms using names that show what they refer to.
   !
   ! After initialization, the names of the source terms can be found via the following string:
   !
   ! CHARACTER(DefaultLength), DIMENSION(MaxNSourceTerms) :: SourceTermName
   !
   !
   ! To facilitate unambiguous specification of the different pathways for DCC,
   ! we have defined the following set of INTEGER tags:
   !
   !   INTEGER, PARAMETER :: NPathways = 3
   !
   !   INTEGER, PARAMETER :: PathwayAir        = 1 : submersion, i.e. external radiation alone
   !   INTEGER, PARAMETER :: PathwayGround     = 2 : external radiation from what is on the ground, assuming no penetration in the ground
   !   INTEGER, PARAMETER :: PathwayInhalation = 3 : inhalation of particles, highest value is taken from available sizes
   !
   ! and their names:
   !
   !   CHARACTER(10), DIMENSION(NPathways), PARAMETER :: PathwayName = &
   !   & (/'Air       ',&!  1
   !   &   'Ground    ',&!  2
   !   &   'Inhalation'/)!  3
   !
   !
   ! The DCCs for submersion have been taken from the supplemental files from the
   ! EDC_Viewer, August 16, 2020, which gives the external dose rate coefficients of ICRP Publication 144.
   !
   ! To facilitate clear reference to either regular DCCs or cumulative DCCs,
   ! we have defined the following set of INTEGER tags:
   !
   !   INTEGER, PARAMETER :: iRegularDCC = 1
   !   INTEGER, PARAMETER :: iCumulativeDCC = 2
   !
   ! and their names:
   !
   !   CHARACTER(13), DIMENSION(2), PARAMETER :: DCCTypeName = &
   !   & (/'DCC          ',&!  1
   !   &   'CumulativeDCC'/)!  2
   !
   !
   ! Two REAL(Float) functions are the actual work horses of this library:
   !
   !   FUNCTION GetCocktailDCC(x,iSourceTerm,iPathway,iType)
   ! This function gives the regular or cumulative cocktail-DCC for a given pathway, cocktail and delay after
   ! blast. The DCC is per 10kt.
   ! The following inputs should be specified:
   !   REAL(Float), INTENT(IN) :: x    = time in seconds since "blast"
   !   INTEGER, INTENT(IN) :: iSourceTerm    = tag specifying the nuclide vector at x=0 (options shown above)
   !   INTEGER, INTENT(IN) :: iPathway = tag specifying the pathway (options shown above)
   !   INTEGER, INTENT(IN) :: iType    = tag specifying if you want to get the regular or the
   !                                     cumulative DCC (options shown above)
   !
   !   FUNCTION GetCocktailNuclide(x,iSourceTerm,MyName)
   ! This function gives the amount of a nuclide for a given cocktail and a given delay after blast.
   ! The amount is per 10kt. For irrelevant nuclides, 0.0 is returned and not an error message.
   ! Where:
   !   REAL(Float), INTENT(IN) :: x       = time in seconds since "blast"
   !   INTEGER, INTENT(IN) :: iSourceTerm       = tag specifying the nuclide vector at x=0 (options shown above)
   !   CHARACTER(*), INTENT(IN) :: MyName = name of the nuclide you want to follow, e.g. 'I-131' or 'Cs137'
   !
   !
   ! An example of how to use libcocktaildcc can be found in test program test_cocktail_dcc.f90.
   !
   USE libxmath
   USE libutil
   USE libinterval
   use libendf

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: GetCocktailDCC,PathwayAir,PathwayGround,PathwayInhalation,NPathways,&
      & NSourceTerms,SourceTermName,PathwayName, DCCTypeName,iRegularDCC,iCumulativeDCC,GetCocktailNuclide

   INTEGER, PARAMETER :: PathwayAir        = 1
   INTEGER, PARAMETER :: PathwayGround     = 2
   INTEGER, PARAMETER :: PathwayInhalation = 3

   INTEGER, PARAMETER :: NPathways = 3

   INTEGER, PARAMETER :: MaxNSourceTerms = 10
   INTEGER :: NSourceTerms = 0

   LOGICAL :: LibCocktailDCCInitialized = .FALSE.

   INTEGER, PARAMETER :: iRegularDCC = 1
   INTEGER, PARAMETER :: iCumulativeDCC = 2

   CHARACTER(13), DIMENSION(2), PARAMETER :: DCCTypeName = &
   & (/'DCC          ',&!  1
   &   'CumulativeDCC'/)!  2

   TYPE(ExponentialIntervalType), DIMENSION(MaxNSourceTerms,NPathways,iRegularDCC:iCumulativeDCC) :: CocktailDCC
   TYPE(ExponentialIntervalType), DIMENSION(MaxNSourceTerms,MaxNuclides) :: LookupCocktail
   CHARACTER(10), DIMENSION(MaxNSourceTerms,MaxNuclides) :: MyNuclideName
   INTEGER, DIMENSION(MaxNSourceTerms) :: NContributingNuclides

   CHARACTER(DefaultLength), DIMENSION(MaxNSourceTerms) :: SourceTermName

   CHARACTER(10), DIMENSION(NPathways), PARAMETER :: PathwayName = &
   & (/'Air       ',&!  1
   &   'Ground    ',&!  2
   &   'Inhalation'/)!  3


CONTAINS
   SUBROUTINE InitLibCocktailDCC()
      !
      ! Read the list of source terms and check availability of lookup tables.
      ! Add missing lookup tables.
      !
      use LibENDF, only: RIVMSourcesPath

      LOGICAL, PARAMETER :: CrashOnError = .FALSE.,DoSilent = .FALSE.
      CHARACTER(DefaultLength) :: ALine,SourceFileName,CheckLookupFileName,Commando,FName,UtilityName
      INTEGER :: TheIndex,Error

      INTEGER, PARAMETER :: DebugLevel = 1

      WRITE(*,'(A)') '--------------------------------------------------------------'
      WRITE(*,*)
      WRITE(*,'(A)') 'Going to initialize libcocktailDCC...'
      WRITE(*,*)

      FName = RIVMSourcesPath() // '/RIVMSources.txt'

      IF (.NOT.FileExists(FName)) THEN
         WRITE(*,'(A)') 'Cannot find list of source term files "'//TRIM(FName)//'"'
      ENDIF

      IF (DebugLevel.GT.0) WRITE(*,'(A)') 'Reading list of source terms from file '//TRIM(FName)

      OPEN(ScrotchFile,FILE=FName,FORM='FORMATTED',ACTION='READ')
      !
      ! Skip header
      !
      READ(ScrotchFile,'(A)') ALine
      DO WHILE (ALine(1:1).EQ.'!')
         READ(ScrotchFile,'(A)') ALine
      ENDDO ! skip header
      !
      ! Read source term names and check lookup tables. Add new lookup tables where necessary.
      !
      NSourceTerms = 0
      DO WHILE (NSourceTerms.LT.MaxNSourceTerms)
         !
         ! Read source file name
         !
         READ(Aline,'(A)',END=10) SourceFileName
         NSourceTerms = NSourceTerms + 1

         IF (DebugLevel.GT.0) THEN
            WRITE(*,'(A,I0,A)') 'Source term ',NSourceTerms,' is '//TRIM(SourceFileName)
            IF (NSourceTerms.EQ.MaxNSourceTerms) WRITE(*,'(A,I0,A)') 'I Reached the maximum number of sources: ',&
            & MaxNSourceTerms,'! Subsequent sources will not be read!'
         ENDIF
         !
         ! Check if name suggests format .RIVMSource
         !
         TheIndex = INDEX(SourceFileName,'.RIVMSource')
         IF (TheIndex.EQ.0) THEN
            WRITE(*,'(A)') 'Source term "'//TRIM(SourceFileName)//'" is not of type .RIVMSource! Exiting!!'
            CALL EXIT()
         ELSE
            IF (DebugLevel.GT.1) WRITE(*,'(A,I0)') 'Extension .RIVMSource starts at position ',TheIndex
         ENDIF

         SourceTermName(NSourceTerms) = SourceFileName(1:(TheIndex-1))
         IF (DebugLevel.GT.1) WRITE(*,'(A)') 'Name of source term is "'//TRIM(SourceTermName(NSourceTerms))//'"'
         DO WHILE (INDEX(SourceTermName(NSourceTerms),'/').GT.0)
            SourceTermName(NSourceTerms) = SourceTermName(NSourceTerms) &
               & ((INDEX(SourceTermName(NSourceTerms),'/')+1):LEN_TRIM(SourceTermName(NSourceTerms)))
            IF (DebugLevel.GT.1) WRITE(*,'(A)') 'Name of source term is "'//TRIM(SourceTermName(NSourceTerms))//'"'
         ENDDO
         IF (DebugLevel.GT.1) WRITE(*,'(A)') 'Name of source term is "'//TRIM(SourceTermName(NSourceTerms))//'"'
         !
         ! Check existence of lookup tables. If not available, add new lookup tables.
         !
         CheckLookupFileName = SourceFileName(1:(TheIndex-1))//'_PinpointAirDoseRates_withprogeny.txt'

         IF (.NOT.FileExists(CheckLookupFileName)) THEN
            IF (DebugLevel.GT.0) WRITE(*,'(A)') 'Cannot find associated lookup tables. Adding new tables...'
            Commando = 'runlog??.txt'
            WRITE(Commando(7:8),'(I2.2)') NSourceTerms

            UtilityName = './build/src/source_term_dose'

            IF (.NOT.FileExists(UtilityName)) THEN
               WRITE(*,'(A)') 'Cannot find utility '//TRIM(UtilityName)//' and its auxiliary files... Exiting!'
               CALL EXIT()
            ENDIF ! utility not found
            Commando = TRIM(UtilityName)//' '//TRIM(SourceFileName)//' yes 0 0 > '//TRIM(Commando)
            Error = SYSTEM(Commando)
            IF (Error.EQ.0) THEN
               WRITE(*,'(A)') 'Execution of external call "'//TRIM(Commando)//'" went well!'
            ELSE
               WRITE(*,'(A)') 'Execution of external call "'//TRIM(Commando)//'" went wrong! Exiting'
               CALL EXIT()
            ENDIF
         ELSE
            IF (DebugLevel.GT.0) WRITE(*,'(A)') 'Found associated lookup tables!'
         ENDIF
         !
         ! Read lookup tables for this source term
         !
         CALL InitializeCocktailDCC(SourceFileName,NSourceTerms)
         CALL InitializeCocktail(SourceFileName,NSourceTerms)

         READ(ScrotchFile,'(A)',END=10) ALine
      ENDDO ! Loop over all source terms declared

      10 CONTINUE
      CLOSE(ScrotchFile)

      LibCocktailDCCInitialized = .TRUE.

      WRITE(*,'(A)') 'Ready initializing libcocktailDCC!'
      WRITE(*,*)
      WRITE(*,'(A)') '--------------------------------------------------------------'
   END SUBROUTINE InitLibCocktailDCC



   SUBROUTINE MakeCumulativeCocktailDCC(iSourceTerm)
      !
      ! Construct cumulative cocktailDCC in Sv (plus the SI dimensions of the volume or surface)
      !
      INTEGER, INTENT(IN) :: iSourceTerm

      INTEGER :: iPathway,iPinpoint
      REAL(Float) :: tLeft,tRight,dDCC,dt

      DO iPathway = 1,NPathways
         CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%IntervalSpecs     = &
         & CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs
         CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%FirstDelay        = &
         & CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%FirstDelay
         CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%DelayGrowthFactor = &
         & CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%DelayGrowthFactor

         ALLOCATE(CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%Values(&
         & 0:CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%IntervalSpecs%N))

         CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%Values = 0._Float

         tLeft = 0._Float
         tRight = CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%FirstDelay
         DO iPinpoint = 1,CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%IntervalSpecs%N
            !
            ! Trapezium rule based on linear interpolation (which should be fair when GrowthFactor is not too far from 1)
            !
            dt = tRight-tLeft
            dDCC = dt * 0.5_Float*(   CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%Values(iPinpoint-1)&
            &                       + CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%Values(iPinpoint  ))

            CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%Values(iPinpoint) = &
            & CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%Values(iPinpoint) + dDCC

            tLeft = tRight
            tRight = tLeft * CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%DelayGrowthFactor
         ENDDO ! loop over pinpoints
      ENDDO ! loop over pathways
   END SUBROUTINE MakeCumulativeCocktailDCC



   SUBROUTINE InitializeCocktailDCC(SourceFileName,iSourceTerm)
      !
      ! Load the lookup tables for the cocktail DCCs in Sv/s (plus the SI dimensions of the volume or surface)
      !
      CHARACTER(*), INTENT(IN) :: SourceFileName
      INTEGER, INTENT(IN) :: iSourceTerm

      INTEGER :: iPathway,iPinpoint,MyPinpoint,TheIndex
      CHARACTER(DefaultLength) :: FName,DumStr
      REAL(Float) :: Dum

      INTEGER, PARAMETER :: DebugLevel = 0

      TheIndex = INDEX(SourceFileName,'.RIVMSource')

      DO iPathway = 1,NPathways
         !
         ! Open file
         !
         FName = SourceFileName(1:(TheIndex-1))//'_Pinpoint'//TRIM(PathwayName(iPathway))//'DoseRates_withprogeny.txt'

         IF (.NOT.FileExists(FName)) THEN
            WRITE(*,'(A)') 'Cannot find lookup table file '//TRIM(FName)//', exiting!'
            CALL EXIT()
         ELSE
            WRITE(*,'(A)') 'Going to read lookup table file '//TRIM(FName)//'...'
         ENDIF

         OPEN(ScratchFile,FILE=FName,FORM='FORMATTED',ACTION='READ')
         !
         ! Skip header
         !
         READ(ScratchFile,*)
         !
         ! Extract interval specs and allocate array
         !
         READ(ScratchFile,*) DumStr,CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%FirstDelay,&
         &                   DumStr,CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%DelayGrowthFactor,&
         &                   DumStr,CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs%N

         CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs%XMin = 0._Float
         CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs%XMax = &
         & CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%FirstDelay &
         & * (CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%DelayGrowthFactor)&
         & **CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs%N

         ALLOCATE(CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%Values(&
         & 0:CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs%N))
         !
         ! Skip header
         !
         READ(ScratchFile,*)
         !
         ! Read actual cocktail DCCs
         !
         DO iPinpoint = 0,CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs%N
            READ(ScratchFile,*) MyPinpoint,Dum,CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%Values(iPinpoint)
            IF (DebugLevel.GT.0) WRITE(*,'(A,I3,A,I3)') 'Pinpoint=',iPinpoint,', MyPinpoint=',MyPinpoint
         ENDDO ! loop over pinpoints
         !
         ! Close file
         !
         CLOSE(ScratchFile)

      ENDDO ! loop over pathways
      !
      ! Construct cumulative cocktailDCC
      !
      CALL MakeCumulativeCocktailDCC(iSourceTerm)
   END SUBROUTINE InitializeCocktailDCC



   FUNCTION GetCocktailDCC(x,iSourceTerm,iPathway,iType)
      !
      ! Give the
      ! regular or cumulative (iType)
      ! cocktail-DCC
      ! for a given pathway (iPathway)
      ! for a given cocktail (iSourceTerm)
      ! and a given delay x [seconds] after start cocktail
      ! The DCC is per 10kt.
      !
      REAL(Float), INTENT(IN) :: x
      INTEGER, INTENT(IN) :: iSourceTerm,iPathway,iType
      REAL(Float) :: GetCocktailDCC

      INTEGER, PARAMETER :: InterpolationWay = 1 ! First order interpolation
      LOGICAL :: Error
      REAL(Float) :: Dum

      IF (.NOT.LibCocktailDCCInitialized) CALL InitLibCocktailDCC()

      IF ((iSourceTerm.LT.1).OR.(iSourceTerm.GT.NSourceTerms)) THEN
         WRITE(*,'(A,I0,A)') 'Invalid source term index: ',iSourceTerm,'! Exiting!!'
         CALL EXIT()
      ENDIF !

      Dum = ExponentialIntervalInterpolate(x,CocktailDCC(iSourceTerm,iPathway,iType),InterpolationWay,Error)
      IF (Error) THEN
         WRITE(*,'(A,G15.5,A,I0,A,I0,A,I0)') 'GetCocktailDCC: Error for x=',x,&
         & ', iSourceTerm=',iSourceTerm,', iPathway=',iPathway,', iType=',iType
         WRITE(*,'(A)') 'Exiting!!'
         CALL EXIT()
      ENDIF
      GetCocktailDCC = Dum
   END FUNCTION GetCocktailDCC



   SUBROUTINE InitializeCocktail(SourceFileName,iSourceTerm)
      !
      ! Load the lookup tables for the cocktails in Bq
      !
      CHARACTER(*), INTENT(IN) :: SourceFileName
      INTEGER, INTENT(IN) :: iSourceTerm

      INTEGER :: iPathway,iPinpoint,MyPinpoint,PinpointN,iNuclide,TheIndex
      CHARACTER(DefaultLength) :: FName,DumStr
      REAL(Float) :: Dum,PinpointFirstDelay,PinpointDelayGrowthFactor

      INTEGER, PARAMETER :: DebugLevel = 0
      !
      ! Open file
      !
      TheIndex = INDEX(SourceFileName,'.RIVMSource')

      FName = SourceFileName(1:(TheIndex-1))//'_PinpointCocktail_withprogeny.txt'

      OPEN(ScratchFile,FILE=FName,FORM='FORMATTED',ACTION='READ')
      !
      ! Extract interval specs and allocate array
      !
      READ(ScratchFile,*) DumStr,PinpointFirstDelay,&
      &                   DumStr,PinpointDelayGrowthFactor,&
      &                   DumStr,PinpointN
      !
      ! Find out which nuclides are contributing this time and initialize lookup tables for them
      !
      READ(ScratchFile,*) NContributingNuclides(iSourceTerm)

      READ(ScratchFile,*) DumStr,DumStr,(MyNuclideName(iSourceTerm,iNuclide),iNuclide=1,NContributingNuclides(iSourceTerm))

      DO iNuclide=1,NContributingNuclides(iSourceTerm)
         CALL AllUpCase(MyNuclideName(iSourceTerm,iNuclide))
         CALL RemoveCharacter(MyNuclideName(iSourceTerm,iNuclide),'-')

         LookupCocktail(iSourceTerm,iNuclide)%FirstDelay         = PinpointFirstDelay
         LookupCocktail(iSourceTerm,iNuclide)%DelayGrowthFactor  = PinpointDelayGrowthFactor
         LookupCocktail(iSourceTerm,iNuclide)%IntervalSpecs%N    = PinpointN
         LookupCocktail(iSourceTerm,iNuclide)%IntervalSpecs%XMin = 0._Float
         LookupCocktail(iSourceTerm,iNuclide)%IntervalSpecs%XMax = &
         & LookupCocktail(iSourceTerm,iNuclide)%FirstDelay &
         & * (LookupCocktail(iSourceTerm,iNuclide)%DelayGrowthFactor)&
         & **LookupCocktail(iSourceTerm,iNuclide)%IntervalSpecs%N

         ALLOCATE(LookupCocktail(iSourceTerm,iNuclide)%Values(0:PinpointN))
      ENDDO ! loop over contributing nuclides
      !
      ! Read actual cocktail
      !
      DO iPinpoint = 0,PinpointN
         READ(ScratchFile,*) MyPinpoint,Dum,&
         & (LookupCocktail(iSourceTerm,iNuclide)%Values(iPinpoint),iNuclide=1,NContributingNuclides(iSourceTerm))
         IF (DebugLevel.GT.0) WRITE(*,'(A,I3,A,I3)') 'Pinpoint=',iPinpoint,', MyPinpoint=',MyPinpoint
      ENDDO ! loop over pinpoints
      !
      ! Close file
      !
      CLOSE(ScratchFile)
   END SUBROUTINE InitializeCocktail



   FUNCTION GetCocktailNuclide(x,iSourceTerm,MyName)
      !
      ! Give the amount of a nuclide for a given cocktail (iSourceTerm)
      ! and a given delay x [seconds] after start cocktail.
      ! The amount is per 10kt.
      !
      REAL(Float), INTENT(IN) :: x
      INTEGER, INTENT(IN) :: iSourceTerm
      CHARACTER(*), INTENT(IN) :: MyName
      REAL(Float) :: GetCocktailNuclide

      INTEGER, PARAMETER :: InterpolationWay = 1 ! First order interpolation
      LOGICAL :: Error,Ready
      REAL(Float) :: Dum
      INTEGER :: iNuclide,TheNuclide
      CHARACTER(10) :: UpName

      IF (.NOT.LibCocktailDCCInitialized) CALL InitLibCocktailDCC()

      IF ((iSourceTerm.LT.1).OR.(iSourceTerm.GT.NSourceTerms)) THEN
         WRITE(*,'(A,I0,A)') 'Invalid source term index: ',iSourceTerm,'! Exiting!!'
         CALL EXIT()
      ENDIF
      !
      ! Identify the nuclide. If found, the give value via lookup,
      ! otherwise return 0 as the nuclide apparently does not participate.
      !
      UpName = MyName
      CALL MassNuc2NucMass(UpName)
      CALL AllUpCase(UpName)
      CALL RemoveCharacter(UpName,'-')
      TheNuclide = 0
      iNuclide = 0
      Ready = .FALSE.
      DO WHILE((.NOT.Ready).AND.(iNuclide.LT.MaxNuclides))
         iNuclide = iNuclide + 1
         IF (TRIM(UpName).EQ.TRIM(MyNuclideName(iSourceTerm,iNuclide))) THEN
            TheNuclide = iNuclide
         ENDIF
      ENDDO ! loop over nuclides

      IF (TheNuclide.EQ.0) THEN
         Dum = 0._Float
      ELSE
         Dum = ExponentialIntervalInterpolate(x,LookupCocktail(iSourceTerm,TheNuclide),InterpolationWay,Error)
         IF (Error) THEN
            WRITE(*,'(A,G15.5,A,I0,A,I0,A,I0,A,A)') 'GetCocktailNuclide: Error for x=',x,&
            & ', iSourceTerm=',iSourceTerm,', iNuclide=',TheNuclide,' for nuclide named ',TRIM(MyName)
            WRITE(*,'(A)') 'Exiting!!'
            CALL EXIT()
         ENDIF ! error
      ENDIF ! nuclide recognized

      GetCocktailNuclide = Dum
   END FUNCTION GetCocktailNuclide
END MODULE LibCocktailDCC
