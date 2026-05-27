module LibCocktailDCC
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
   ! character(DefaultLength), dimension(MaxNSourceTerms) :: SourceTermName
   !
   !
   ! To facilitate unambiguous specification of the different pathways for DCC,
   ! we have defined the following set of integer tags:
   !
   !   integer, parameter :: NPathways = 3
   !
   !   integer, parameter :: PathwayAir        = 1 : submersion, i.e. external radiation alone
   !   integer, parameter :: PathwayGround     = 2 : external radiation from what is on the ground, assuming no penetration in the ground
   !   integer, parameter :: PathwayInhalation = 3 : inhalation of particles, highest value is taken from available sizes
   !
   ! and their names:
   !
   !   character(10), dimension(NPathways), parameter :: PathwayName = &
   !   & (/'Air       ',&!  1
   !   &   'Ground    ',&!  2
   !   &   'Inhalation'/)!  3
   !
   !
   ! The DCCs for submersion have been taken from the supplemental files from the
   ! EDC_Viewer, August 16, 2020, which gives the external dose rate coefficients of ICRP Publication 144.
   !
   ! To facilitate clear reference to either regular DCCs or cumulative DCCs,
   ! we have defined the following set of integer tags:
   !
   !   integer, parameter :: iRegularDCC = 1
   !   integer, parameter :: iCumulativeDCC = 2
   !
   ! and their names:
   !
   !   character(13), dimension(2), parameter :: DCCTypeName = &
   !   & (/'DCC          ',&!  1
   !   &   'CumulativeDCC'/)!  2
   !
   !
   ! Two real(Float) functions are the actual work horses of this library:
   !
   !   function GetCocktailDCC(x,iSourceTerm,iPathway,iType)
   ! This function gives the regular or cumulative cocktail-DCC for a given pathway, cocktail and delay after
   ! blast. The DCC is per 10kt.
   ! The following inputs should be specified:
   !   real(Float), intent(in) :: x    = time in seconds since "blast"
   !   integer, intent(in) :: iSourceTerm    = tag specifying the nuclide vector at x=0 (options shown above)
   !   integer, intent(in) :: iPathway = tag specifying the pathway (options shown above)
   !   integer, intent(in) :: iType    = tag specifying if you want to get the regular or the
   !                                     cumulative DCC (options shown above)
   !
   !   function GetCocktailNuclide(x,iSourceTerm,MyName)
   ! This function gives the amount of a nuclide for a given cocktail and a given delay after blast.
   ! The amount is per 10kt. For irrelevant nuclides, 0.0 is returned and not an error message.
   ! Where:
   !   real(Float), intent(in) :: x       = time in seconds since "blast"
   !   integer, intent(in) :: iSourceTerm       = tag specifying the nuclide vector at x=0 (options shown above)
   !   character(*), intent(in) :: MyName = name of the nuclide you want to follow, e.g. 'I-131' or 'Cs137'
   !
   !
   ! An example of how to use libcocktaildcc can be found in test program test_cocktail_dcc.f90.
   !
   use libxmath
   use libutil
   use libinterval
   use libendf

   implicit none

   private
   public :: GetCocktailDCC,PathwayAir,PathwayGround,PathwayInhalation,NPathways,&
      & NSourceTerms,SourceTermName,PathwayName, DCCTypeName,iRegularDCC,iCumulativeDCC,GetCocktailNuclide

   integer, parameter :: PathwayAir        = 1
   integer, parameter :: PathwayGround     = 2
   integer, parameter :: PathwayInhalation = 3

   integer, parameter :: NPathways = 3

   integer, parameter :: MaxNSourceTerms = 10
   integer :: NSourceTerms = 0

   logical :: LibCocktailDCCInitialized = .FALSE.

   integer, parameter :: iRegularDCC = 1
   integer, parameter :: iCumulativeDCC = 2

   character(13), dimension(2), parameter :: DCCTypeName = &
   & ['DCC          ',&!  1
   &  'CumulativeDCC']!  2

   type(ExponentialIntervalType), dimension(MaxNSourceTerms,NPathways,iRegularDCC:iCumulativeDCC) :: CocktailDCC
   type(ExponentialIntervalType), dimension(MaxNSourceTerms,MaxNuclides) :: LookupCocktail
   character(10), dimension(MaxNSourceTerms,MaxNuclides) :: MyNuclideName
   integer, dimension(MaxNSourceTerms) :: NContributingNuclides

   character(DefaultLength), dimension(MaxNSourceTerms) :: SourceTermName

   character(10), dimension(NPathways), parameter :: PathwayName = &
   & ['Air       ',&!  1
   &  'Ground    ',&!  2
   &  'Inhalation']!  3


contains
   subroutine InitLibCocktailDCC()
      !
      ! Read the list of source terms and check availability of lookup tables.
      ! Add missing lookup tables.
      !
      use LibENDF, only: RIVMSourcesPath

      character(DefaultLength) :: ALine,SourceFileName,CheckLookupFileName,Commando,FName,UtilityName
      integer :: TheIndex,Error

      integer, parameter :: DebugLevel = 1

      write(*,'(A)') '--------------------------------------------------------------'
      write(*,*)
      write(*,'(A)') 'Going to initialize libcocktailDCC...'
      write(*,*)

      FName = RIVMSourcesPath() // '/RIVMSources.txt'

      if (.NOT.FileExists(FName)) then
         write(*,'(A)') 'Cannot find list of source term files "'//trim(FName)//'"'
      endif

      if (DebugLevel > 0) write(*,'(A)') 'Reading list of source terms from file '//trim(FName)

      open(ScrotchFile,file=FName,form='FORMATTED',action='READ')
      !
      ! Skip header
      !
      read(ScrotchFile,'(A)') ALine
      do while (ALine(1:1) == '!')
         read(ScrotchFile,'(A)') ALine
      enddo ! skip header
      !
      ! Read source term names and check lookup tables. Add new lookup tables where necessary.
      !
      NSourceTerms = 0
      do while (NSourceTerms < MaxNSourceTerms)
         !
         ! Read source file name
         !
         read(Aline,'(A)',end=10) SourceFileName
         NSourceTerms = NSourceTerms + 1

         if (DebugLevel > 0) then
            write(*,'(A,I0,A)') 'Source term ',NSourceTerms,' is '//trim(SourceFileName)
            if (NSourceTerms == MaxNSourceTerms) write(*,'(A,I0,A)') 'I Reached the maximum number of sources: ',&
            & MaxNSourceTerms,'! Subsequent sources will not be read!'
         endif
         !
         ! Check if name suggests format .RIVMSource
         !
         TheIndex = INDEX(SourceFileName,'.RIVMSource')
         if (TheIndex == 0) then
            write(*,'(A)') 'Source term "'//trim(SourceFileName)//'" is not of type .RIVMSource! Exiting!!'
            call exit()
         else
            if (DebugLevel > 1) write(*,'(A,I0)') 'Extension .RIVMSource starts at position ',TheIndex
         endif

         SourceTermName(NSourceTerms) = SourceFileName(1:(TheIndex-1))
         if (DebugLevel > 1) write(*,'(A)') 'Name of source term is "'//trim(SourceTermName(NSourceTerms))//'"'
         do while (INDEX(SourceTermName(NSourceTerms),'/') > 0)
            SourceTermName(NSourceTerms) = SourceTermName(NSourceTerms) &
               & ((INDEX(SourceTermName(NSourceTerms),'/')+1):len_trim(SourceTermName(NSourceTerms)))
            if (DebugLevel > 1) write(*,'(A)') 'Name of source term is "'//trim(SourceTermName(NSourceTerms))//'"'
         enddo
         if (DebugLevel > 1) write(*,'(A)') 'Name of source term is "'//trim(SourceTermName(NSourceTerms))//'"'
         !
         ! Check existence of lookup tables. If not available, add new lookup tables.
         !
         CheckLookupFileName = SourceFileName(1:(TheIndex-1))//'_PinpointAirDoseRates_withprogeny.txt'

         if (.NOT.FileExists(CheckLookupFileName)) then
            if (DebugLevel > 0) write(*,'(A)') 'Cannot find associated lookup tables. Adding new tables...'
            Commando = 'runlog??.txt'
            write(Commando(7:8),'(I2.2)') NSourceTerms

            UtilityName = './build/src/source_term_dose'

            if (.NOT.FileExists(UtilityName)) then
               write(*,'(A)') 'Cannot find utility '//trim(UtilityName)//' and its auxiliary files... Exiting!'
               call exit()
            endif ! utility not found
            Commando = trim(UtilityName)//' '//trim(SourceFileName)//' yes 0 0 > ' &
            &          // TransitionMatrixPath() // '/' // trim(Commando)
            Error = SYSTEM(Commando)
            if (Error == 0) then
               write(*,'(A)') 'Execution of external call "'//trim(Commando)//'" went well!'
            else
               write(*,'(A)') 'Execution of external call "'//trim(Commando)//'" went wrong! Exiting'
               call exit()
            endif
         else
            if (DebugLevel > 0) write(*,'(A)') 'Found associated lookup tables!'
         endif
         !
         ! Read lookup tables for this source term
         !
         call InitializeCocktailDCC(SourceFileName,NSourceTerms)
         call InitializeCocktail(SourceFileName,NSourceTerms)

         read(ScrotchFile,'(A)',end=10) ALine
      enddo ! Loop over all source terms declared

      10 continue
      close(ScrotchFile)

      LibCocktailDCCInitialized = .TRUE.

      write(*,'(A)') 'Ready initializing libcocktailDCC!'
      write(*,*)
      write(*,'(A)') '--------------------------------------------------------------'
   end subroutine InitLibCocktailDCC



   subroutine MakeCumulativeCocktailDCC(iSourceTerm)
      !
      ! Construct cumulative cocktailDCC in Sv (plus the SI dimensions of the volume or surface)
      !
      integer, intent(in) :: iSourceTerm

      integer :: iPathway,iPinpoint
      real(Float) :: tLeft,tRight,dDCC,dt

      do iPathway = 1,NPathways
         CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%IntervalSpecs     = &
         & CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs
         CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%FirstDelay        = &
         & CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%FirstDelay
         CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%DelayGrowthFactor = &
         & CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%DelayGrowthFactor

         allocate(CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%Values(&
         & 0:CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%IntervalSpecs%N))

         CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%Values = 0._Float

         tLeft = 0._Float
         tRight = CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%FirstDelay
         do iPinpoint = 1,CocktailDCC(iSourceTerm,iPathway,iCumulativeDCC)%IntervalSpecs%N
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
         enddo ! loop over pinpoints
      enddo ! loop over pathways
   end subroutine MakeCumulativeCocktailDCC



   subroutine InitializeCocktailDCC(SourceFileName,iSourceTerm)
      !
      ! Load the lookup tables for the cocktail DCCs in Sv/s (plus the SI dimensions of the volume or surface)
      !
      character(*), intent(in) :: SourceFileName
      integer, intent(in) :: iSourceTerm

      integer :: iPathway,iPinpoint,MyPinpoint,TheIndex
      character(DefaultLength) :: FName,DumStr
      real(Float) :: Dum

      integer, parameter :: DebugLevel = 0

      TheIndex = INDEX(SourceFileName,'.RIVMSource')

      do iPathway = 1,NPathways
         !
         ! Open file
         !
         FName = SourceFileName(1:(TheIndex-1))//'_Pinpoint'//trim(PathwayName(iPathway))//'DoseRates_withprogeny.txt'

         if (.NOT.FileExists(FName)) then
            write(*,'(A)') 'Cannot find lookup table file '//trim(FName)//', exiting!'
            call exit()
         else
            write(*,'(A)') 'Going to read lookup table file '//trim(FName)//'...'
         endif

         open(ScratchFile,file=FName,form='FORMATTED',action='READ')
         !
         ! Skip header
         !
         read(ScratchFile,*)
         !
         ! Extract interval specs and allocate array
         !
         read(ScratchFile,*) DumStr,CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%FirstDelay,&
         &                   DumStr,CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%DelayGrowthFactor,&
         &                   DumStr,CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs%N

         CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs%XMin = 0._Float
         CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs%XMax = &
         & CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%FirstDelay &
         & * (CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%DelayGrowthFactor)&
         & **CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs%N

         allocate(CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%Values(&
         & 0:CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs%N))
         !
         ! Skip header
         !
         read(ScratchFile,*)
         !
         ! Read actual cocktail DCCs
         !
         do iPinpoint = 0,CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%IntervalSpecs%N
            read(ScratchFile,*) MyPinpoint,Dum,CocktailDCC(iSourceTerm,iPathway,iRegularDCC)%Values(iPinpoint)
            if (DebugLevel > 0) write(*,'(A,I3,A,I3)') 'Pinpoint=',iPinpoint,', MyPinpoint=',MyPinpoint
         enddo ! loop over pinpoints
         !
         ! Close file
         !
         close(ScratchFile)

      enddo ! loop over pathways
      !
      ! Construct cumulative cocktailDCC
      !
      call MakeCumulativeCocktailDCC(iSourceTerm)
   end subroutine InitializeCocktailDCC



   function GetCocktailDCC(x,iSourceTerm,iPathway,iType)
      !
      ! Give the
      ! regular or cumulative (iType)
      ! cocktail-DCC
      ! for a given pathway (iPathway)
      ! for a given cocktail (iSourceTerm)
      ! and a given delay x [seconds] after start cocktail
      ! The DCC is per 10kt.
      !
      real(Float), intent(in) :: x
      integer, intent(in) :: iSourceTerm,iPathway,iType
      real(Float) :: GetCocktailDCC

      integer, parameter :: InterpolationWay = 1 ! First order interpolation
      logical :: Error
      real(Float) :: Dum

      if (.NOT.LibCocktailDCCInitialized) call InitLibCocktailDCC()

      if ((iSourceTerm < 1).OR.(iSourceTerm > NSourceTerms)) then
         write(*,'(A,I0,A)') 'Invalid source term index: ',iSourceTerm,'! Exiting!!'
         call exit()
      endif !

      Dum = ExponentialIntervalInterpolate(x,CocktailDCC(iSourceTerm,iPathway,iType),InterpolationWay,Error)
      if (Error) then
         write(*,'(A,G15.5,A,I0,A,I0,A,I0)') 'GetCocktailDCC: Error for x=',x,&
         & ', iSourceTerm=',iSourceTerm,', iPathway=',iPathway,', iType=',iType
         write(*,'(A)') 'Exiting!!'
         call exit()
      endif
      GetCocktailDCC = Dum
   end function GetCocktailDCC



   subroutine InitializeCocktail(SourceFileName,iSourceTerm)
      !
      ! Load the lookup tables for the cocktails in Bq
      !
      character(*), intent(in) :: SourceFileName
      integer, intent(in) :: iSourceTerm

      integer :: iPinpoint,MyPinpoint,PinpointN,iNuclide,TheIndex
      character(DefaultLength) :: FName,DumStr
      real(Float) :: Dum,PinpointFirstDelay,PinpointDelayGrowthFactor

      integer, parameter :: DebugLevel = 0
      !
      ! Open file
      !
      TheIndex = INDEX(SourceFileName,'.RIVMSource')

      FName = SourceFileName(1:(TheIndex-1))//'_PinpointCocktail_withprogeny.txt'

      open(ScratchFile,file=FName,form='FORMATTED',action='READ')
      !
      ! Extract interval specs and allocate array
      !
      read(ScratchFile,*) DumStr,PinpointFirstDelay,&
      &                   DumStr,PinpointDelayGrowthFactor,&
      &                   DumStr,PinpointN
      !
      ! Find out which nuclides are contributing this time and initialize lookup tables for them
      !
      read(ScratchFile,*) NContributingNuclides(iSourceTerm)

      read(ScratchFile,*) DumStr,DumStr,(MyNuclideName(iSourceTerm,iNuclide),iNuclide=1,NContributingNuclides(iSourceTerm))

      do iNuclide=1,NContributingNuclides(iSourceTerm)
         call AllUpCase(MyNuclideName(iSourceTerm,iNuclide))
         call RemoveCharacter(MyNuclideName(iSourceTerm,iNuclide),'-')

         LookupCocktail(iSourceTerm,iNuclide)%FirstDelay         = PinpointFirstDelay
         LookupCocktail(iSourceTerm,iNuclide)%DelayGrowthFactor  = PinpointDelayGrowthFactor
         LookupCocktail(iSourceTerm,iNuclide)%IntervalSpecs%N    = PinpointN
         LookupCocktail(iSourceTerm,iNuclide)%IntervalSpecs%XMin = 0._Float
         LookupCocktail(iSourceTerm,iNuclide)%IntervalSpecs%XMax = &
         & LookupCocktail(iSourceTerm,iNuclide)%FirstDelay &
         & * (LookupCocktail(iSourceTerm,iNuclide)%DelayGrowthFactor)&
         & **LookupCocktail(iSourceTerm,iNuclide)%IntervalSpecs%N

         allocate(LookupCocktail(iSourceTerm,iNuclide)%Values(0:PinpointN))
      enddo ! loop over contributing nuclides
      !
      ! Read actual cocktail
      !
      do iPinpoint = 0,PinpointN
         read(ScratchFile,*) MyPinpoint,Dum,&
         & (LookupCocktail(iSourceTerm,iNuclide)%Values(iPinpoint),iNuclide=1,NContributingNuclides(iSourceTerm))
         if (DebugLevel > 0) write(*,'(A,I3,A,I3)') 'Pinpoint=',iPinpoint,', MyPinpoint=',MyPinpoint
      enddo ! loop over pinpoints
      !
      ! Close file
      !
      close(ScratchFile)
   end subroutine InitializeCocktail



   function GetCocktailNuclide(x,iSourceTerm,MyName)
      !
      ! Give the amount of a nuclide for a given cocktail (iSourceTerm)
      ! and a given delay x [seconds] after start cocktail.
      ! The amount is per 10kt.
      !
      real(Float), intent(in) :: x
      integer, intent(in) :: iSourceTerm
      character(*), intent(in) :: MyName
      real(Float) :: GetCocktailNuclide

      integer, parameter :: InterpolationWay = 1 ! First order interpolation
      logical :: Error,Ready
      real(Float) :: Dum
      integer :: iNuclide,TheNuclide
      character(10) :: UpName

      if (.NOT.LibCocktailDCCInitialized) call InitLibCocktailDCC()

      if ((iSourceTerm < 1).OR.(iSourceTerm > NSourceTerms)) then
         write(*,'(A,I0,A)') 'Invalid source term index: ',iSourceTerm,'! Exiting!!'
         call exit()
      endif
      !
      ! Identify the nuclide. If found, the give value via lookup,
      ! otherwise return 0 as the nuclide apparently does not participate.
      !
      UpName = MyName
      call MassNuc2NucMass(UpName)
      call AllUpCase(UpName)
      call RemoveCharacter(UpName,'-')
      TheNuclide = 0
      iNuclide = 0
      Ready = .FALSE.
      do while((.NOT.Ready).AND.(iNuclide < MaxNuclides))
         iNuclide = iNuclide + 1
         if (trim(UpName) == trim(MyNuclideName(iSourceTerm,iNuclide))) then
            TheNuclide = iNuclide
         endif
      enddo ! loop over nuclides

      if (TheNuclide == 0) then
         Dum = 0._Float
      else
         Dum = ExponentialIntervalInterpolate(x,LookupCocktail(iSourceTerm,TheNuclide),InterpolationWay,Error)
         if (Error) then
            write(*,'(A,G15.5,A,I0,A,I0,A,I0,A,A)') 'GetCocktailNuclide: Error for x=',x,&
            & ', iSourceTerm=',iSourceTerm,', iNuclide=',TheNuclide,' for nuclide named ',trim(MyName)
            write(*,'(A)') 'Exiting!!'
            call exit()
         endif ! error
      endif ! nuclide recognized

      GetCocktailNuclide = Dum
   end function GetCocktailNuclide
end module LibCocktailDCC
