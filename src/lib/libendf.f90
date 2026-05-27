module LibENDF
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

   use Matrix_Exponential
   use LibUtil, only: DefaultLength, ScratchFile, &
      env_var, FileExists, &
      RemoveCharacter, Capitalize, CharacterIsADigit
   use LibXMath, only: Float, dp, Inc, &
      SparseMatrix, SparseLogicalMatrix, Matrix2SparseMatrix, LMatrix2SparseLMatrix, SparseMatrix2Matrix, SparseLMatrix2LMatrix
   implicit none

   private
   public :: RIVMSourcesPath, TransitionMatrixPath
   public :: NNuclides,NuclideSpecs,MakeMotherDaughterMatrix,&
   & MotherDaughterMatrix,FindOrphans,RegularizeNuclides,&
   & NuclideType,MaxNuclides,Orphanage,NuclideFamily,&
   & CollectProgeny,ListTransitions,MakeEvolutionMatrix,AddOrphansBelow,&
   & GetNuclideNumber,MassNuc2NucMass,EnsureHyphen,ReadNuclideSpecs,AtomName,&
   & IsRelated,RegularizedIsRelated,ReadNProcessENDFNuclideSpecs
   public :: isGas

   integer, parameter :: MaxAtoms = 118 ! largest atom number available

   character(2), dimension(0:MaxAtoms) :: AtomName
   !
   ! Number of nuclide groups
   !
   integer, parameter :: NNuclideGroups = 11
   character(5), dimension(0:NNuclideGroups), parameter :: NuclideGroupName &
      & = ['Other',&
      &     'Noble',&
      &     'I    ',&
      &     'Cs   ',&
      &     'Te   ',&
      &     'Sb   ',&
      &     'Ba   ',&
      &     'Sr   ',&
      &     'Mo   ',&
      &     'La   ',&
      &     'Ce   ',&
      &     'U    ']

   type AtomType
      character(13) :: Name
      character(3) :: Symbol
      integer :: Period,Group,NuclideGroup
      character(24) :: ChemicalSeries
      real(Float) :: Mass
   end type AtomType

   integer :: NAtoms = 0
   type(AtomType), dimension(0:MaxAtoms) :: AtomSpecs
   !
   ! The maximum number of nuclides that we will read from the database:
   !
   integer, parameter :: MaxNuclides = 4000
   !
   !
   integer, parameter :: MaxNDaughters = 15 ! Sufficient for "normal decay", including assignment of very fast decays to their grandmothers
   !
   ! The true number of nuclides found:
   !
   integer :: NNuclides

   type NuclideType
      character(10) :: NuclideName
      character(10), dimension(MaxNDaughters) :: DaughterName
      character(40), dimension(MaxNDaughters) :: DecayName
      integer :: NDecayModes,AtomNumber,NHadrons,NDaughters,MetaStableMode,NuclideGroup
      character(15) :: AtomName
      integer, dimension(MaxNDaughters) :: Daughter
      real(Float) :: HalfTime,AtomMass
      real(Float), dimension(MaxNDaughters) :: DaughterFraction,DecayFraction
      logical :: IsOrphan
   end type NuclideType
   !
   ! The whole lot of information will be stored in:
   !
   type(NuclideType), dimension(0:MaxNuclides) :: NuclideSpecs
   !
   ! A single datafile can have so many lines
   !
   integer, parameter :: MaxNLines = 30000

   integer, dimension(MaxNLines) :: iRecordType
   logical, dimension(MaxNLines) :: IsContinuation
   !
   ! The following types of records can be found:
   !
   integer, parameter :: iUnknownRecord =  0
   integer, parameter :: iTextRecord    =  1
   integer, parameter :: iContRecord    =  2
   integer, parameter :: iListRecord    =  3
   integer, parameter :: iTab1Record    =  4
   integer, parameter :: iTab2Record    =  5
   integer, parameter :: iIntgRecord    =  6
   integer, parameter :: iDecayRecord   =  7

   integer, parameter :: NRecordTypes   =  7 ! should be same as last number above

   character(7), dimension(0:NRecordTypes), parameter :: RecordName = &
      & ['Unknown',&!   0
      &   'Text   ',&!   1
      &   'Cont   ',&!   2
      &   'List   ',&!   3
      &   'Tab1   ',&!   4
      &   'Tab2   ',&!   5
      &   'Intg   ',&!   6
      &   'Decay  ']!   7

   character(1), dimension(0:4), parameter :: MetastabilityName = &
      & [' ',&!   0
      &   'm',&!   1
      &   'n',&!   2
      &   'o',&!   3
      &   'p']!   4

   integer, parameter :: iDecayModeGamma              =  0
   integer, parameter :: iDecayModeBeta               =  1
   integer, parameter :: iDecayModeElectronCapture    =  2
   integer, parameter :: iDecayModeIsomericTransition =  3
   integer, parameter :: iDecayModeAlpha              =  4
   integer, parameter :: iDecayModeNeutron            =  5
   integer, parameter :: iDecayModeSpontaneousFission =  6
   integer, parameter :: iDecayModeProton             =  7
   integer, parameter :: iDecayModeUnknown            = 10

   character(19), dimension(0:10), parameter :: DecayModeName = &
      & ['Gamma              ',&!   0
      &   'Beta               ',&!   1
      &   'Electron capture   ',&!   2
      &   'Isomeric transition',&!   3
      &   'Alpha              ',&!   4
      &   'Neutron            ',&!   5
      &   'Spontaneous fission',&!   6
      &   'Proton             ',&!   7
      &   '-invalid1-         ',&!   8
      &   '-invalid2-         ',&!   9
      &   'Unknown            ']!  10

   integer, dimension(0:10), parameter :: AtomStep = &
      & [ 0,& ! Gamma
      &    1,& ! Beta
      &   -1,& ! Electron capture
      &    0,& ! Isomeric transition
      &   -2,& ! Alpha
      &    0,& ! Neutron
      & 1000,& ! Spontaneous fission
      &   -1,& ! Proton
      & 1000,& ! -invalid1-
      & 1000,& ! -invalid2-
      & 1000]  ! Unknown

   integer, dimension(0:10), parameter :: MassStep = &
      & [ 0,& ! Gamma
      &    0,& ! Beta
      &    0,& ! Electron capture
      &    0,& ! Isomeric transition
      &   -4,& ! Alpha
      &   -1,& ! Neutron
      & 1000,& ! Spontaneous fission
      &   -1,& ! Proton
      & 1000,& ! -invalid1-
      & 1000,& ! -invalid2-
      & 1000]  ! Unknown

   real(Float), dimension(MaxNuclides,MaxNuclides) :: MotherDaughterMatrix

   type OrphanType
      integer :: Mother,Daughter
      real(Float) :: Yield,CompoundYield
      logical :: Active
   end type OrphanType

   integer, parameter :: MaxNOrphans = 1700

   type OrphanageType
      real(Float) :: tMin
      integer :: NOrphans,NTooUnStable
      logical, dimension(MaxNuclides) :: IsTooUnstable
      type(OrphanType), dimension(MaxNOrphans) :: Orphan
      logical, dimension(MaxNuclides) :: OrphanActive
   end type OrphanageType

   type(OrphanageType) :: Orphanage

   logical, dimension(MaxNuclides,MaxNuclides) :: IsRelated,RegularizedIsRelated
   !
   ! Structure for storing relatives stemming from 1 nuclide
   !
   type NuclideFamilyType
      integer :: NFamily,NFarFamily
      ! Equivalent of compound yield: if you start with any (specified) family member with activity 1, how much of the other
      ! nuclides do you find? Only contributions are counted that are linked to the starting nuclide via fast-decaying nuclides.
      real(Float), dimension(MaxNuclides) :: RemainingFraction
      integer, dimension(MaxNuclides) :: FamilyMember,NDaughters,NMothers,FamilyNumber
      integer, dimension(MaxNuclides,MaxNuclides) :: Mother,Daughter
      logical, dimension(MaxNuclides) :: IsFamily,IsFarFamily,IncomingArrowsChecked,OutgoingArrowsChecked,&
      & DescendentOfStartMember
      logical, dimension(MaxNuclides,MaxNuclides) :: IsRelated
   end type NuclideFamilyType

   type(NuclideFamilyType) :: NuclideFamily

contains
   function ENDFPath()
      character(:), allocatable :: ENDFPath

      call env_var('COCKTAIL_DCC_ENDF_DIR', ENDFPath)
      if (.not. allocated(ENDFPath)) error stop 'need to set environment variable: COCKTAIL_DCC_ENDF_DIR'
   end function

   function RIVMSourcesPath()
      character(:), allocatable :: RIVMSourcesPath

      call env_var('COCKTAIL_DCC_SOURCES_DIR', RIVMSourcesPath)
      if (.not. allocated(RIVMSourcesPath)) RIVMSourcesPath = './'
   end function

   function TransitionMatrixPath()
      character(:), allocatable :: TransitionMatrixPath

      call env_var('COCKTAIL_DCC_OUTPUT_DIR', TransitionMatrixPath)
      if (.not. allocated(TransitionMatrixPath)) TransitionMatrixPath = '.'
   end function

   integer function GetAtomNumber(Symbol)
      !
      ! Search for the given name of an atom in the list of known atoms and
      ! return its number
      !
      character(*), intent(in) :: Symbol
      integer :: DumInt,iAtom,TheIndex
      character(10) :: DumSymbol

      DumSymbol = trim(Symbol)
      call Capitalize(DumSymbol)
      TheIndex = INDEX(DumSymbol,'-')
      if (TheIndex > 0) DumSymbol = DumSymbol(1:(TheIndex-1))
      DumInt = 0
      do iAtom = 1,NAtoms
         if (trim(DumSymbol) == trim(AtomSpecs(iAtom)%Symbol)) DumInt = iAtom
      enddo
      GetAtomNumber = DumInt
   end function GetAtomNumber



   function IsGas(nuclide)
      ! Tell me whether the nuclide with a given name is a gas (Kr, Xe, H, Ar)
      type(NuclideType), intent(in) :: nuclide

      logical IsGas, DumLog

      DumLog = (nuclide%AtomNumber == 1) & ! H , hydrogen
         .OR.(nuclide%AtomNumber == 2) &   ! He, Helium
         .OR.(nuclide%AtomNumber == 10) &  ! Ne, Neon
         .OR.(nuclide%AtomNumber == 18) &  ! Ar, Argon
         .OR.(nuclide%AtomNumber == 36) &  ! Kr, Krypton
         .OR.(nuclide%AtomNumber == 54) &  ! Xe, Xenon
         .OR.(nuclide%AtomNumber == 86) &  ! Rn, Radon
         .OR.(nuclide%AtomNumber == 118)   ! Og, Oganesson

      IsGas = DumLog
   end function



   subroutine InitNuclide(iNuclide)
      !
      ! Set all info to 0 or anything else that shows clearly that the nulcide is stil empty
      !
      integer, intent(in) :: iNuclide

      NuclideSpecs(iNuclide)%NuclideName = ' '
      NuclideSpecs(iNuclide)%DaughterName = ' '
      NuclideSpecs(iNuclide)%AtomName = ' '
      NuclideSpecs(iNuclide)%AtomNumber = 0
      NuclideSpecs(iNuclide)%NuclideGroup = 0
      NuclideSpecs(iNuclide)%Daughter = 0
      NuclideSpecs(iNuclide)%HalfTime = 0._Float
      NuclideSpecs(iNuclide)%DaughterFraction = 0._Float
   end subroutine InitNuclide



   subroutine ReadICRPNuclides(InName)
      !
      ! Read the names, decaytimes and daughters for a collection of nuclides
      !
      character(*), intent(in) :: InName

      integer :: iNuclide,jNuclide,DumInt,iDaughter,iAtom,i
      character(DefaultLength) :: ALine,DumStr
      character(2) :: Units
      character(7) :: DumStr7
      character(7), dimension(MaxNuclides) :: NuclideName
      character(11), dimension(4) :: YieldStr
      real(Float) :: C,Dum
      logical, dimension(MaxNuclides) :: Available
      integer, dimension(MaxNuclides) :: SortedIndex,AtomNumber

      integer, parameter :: DebugLevel = 0
      !
      ! Find number of nuclides and sort them before reading the true specs
      !
      open(ScratchFile,file=trim(InName),form='FORMATTED',status='OLD',&
         & position='REWIND')
      read(ScratchFile,*) ! Skip header

      if(DebugLevel > 1) write(*,*)
      if(DebugLevel > 1) write(*,'(A)') 'First pass: the nuclides in the order as in the file:'
      iNuclide = 0
      do
         read(ScratchFile,'(A7)',end=20) DumStr7

         iNuclide = iNuclide + 1
         NuclideName(iNuclide) = trim(DumStr7)
         AtomNumber(iNuclide) = GetAtomNumber(NuclideName(iNuclide))
         if(DebugLevel > 1) write(*,'(I4,1X,A7,5X,I3,1X,A20)') iNuclide,NuclideName(iNuclide),&
            & AtomNumber(iNuclide),AtomSpecs(AtomNumber(iNuclide))%Name
            !
            ! Initialize nuclide
            !
         call InitNuclide(iNuclide)

      enddo
      20 continue ! ready counting nuclides

         NNuclides = iNuclide ! The number of nuclides actually found in the file
      if(DebugLevel > 1) write(*,*)
      !
      ! Sort the nuclides by atom number in stead of alphabetically
      !
      if(DebugLevel > 1) write(*,'(A)') 'Second pass: sorting the nuclides:'
      Available = .TRUE.
      SortedIndex = 0
      i = 0
      do iAtom=1,NAtoms
         do iNuclide = 1,NNuclides
            if ((AtomNumber(iNuclide) == iAtom).AND.(Available(iNuclide))) then
               i = i + 1
               SortedIndex(iNuclide) = i
               Available(iNuclide) = .FALSE.
               if(DebugLevel > 1) write(*,'(I4,1X,I3,1X,A20)') i,iAtom,trim(NuclideName(iNuclide))
            endif
         enddo
      enddo
      if(DebugLevel > 1) write(*,*)
      !
      ! Read the specs
      !
      if(DebugLevel > 1) write(*,'(A)') 'Third pass: reading the specs:'
      REWIND(ScratchFile)
      read(ScratchFile,*) ! Skip header

      do jNuclide = 1,NNuclides
         iNuclide = SortedIndex(jNuclide)
         read(ScratchFile,10,end=20) ALine
         10 format(A)
         if (DebugLevel > 1) write(*,70) jNuclide,iNuclide
         70 format('Line: ',I4,': interpreting specs for sorted nuclide ',I4)
         if (DebugLevel > 2) write(*,60) trim(ALine)
         60 format('Specs from database are: "',A,'"')
         !
         ! Get the nuclide name
         !
         read(ALine,100) NuclideSpecs(iNuclide)%NuclideName,NuclideSpecs(iNuclide)%HalfTime,Units,&
         & DumStr,DumInt,DumInt,DumInt,DumInt,&
         & (NuclideSpecs(iNuclide)%DaughterName(iDaughter),DumStr,&
         &   YieldStr(iDaughter),iDaughter=1,4),Dum,Dum,Dum,DumInt,DumInt,DumInt,DumInt,DumInt,&
         & NuclideSpecs(iNuclide)%AtomMass
         NuclideSpecs(iNuclide)%NHadrons = NINT(NuclideSpecs(iNuclide)%AtomMass)

         100 format(A7,F8.0,A2,A8,3I7,I6,4(1X,A7,A6,A11),f7.0,2f8.0,3i4,i5,i4,e11.0,e10.0,e9.0)
         !
         ! Calibrate halftime
         !
         if (trim(Units) == 'y') then
            C = 3600._Float*24._Float*365.25_Float
            if (DebugLevel > 1) write(*,'(A,EN20.5)') &
            & 'Half time is given in year --> conversion to seconds with factor ',C
         else if (trim(Units) == 'd') then
            C = 3600._Float*24._Float
            if (DebugLevel > 1) write(*,'(A,EN20.5)') &
            & 'Half time is given in day --> conversion to seconds with factor ',C
         else if (trim(Units) == 'h') then
            C = 3600._Float
            if (DebugLevel > 1) write(*,'(A,EN20.5)') &
            & 'Half time is given in hour --> conversion to seconds with factor ',C
         else if (trim(Units) == 'm') then
            C = 60._Float
            if (DebugLevel > 1) write(*,'(A,EN20.5)') &
            & 'Half time is given in minute --> conversion to seconds with factor ',C
         else if (trim(Units) == 's') then
            C = 1._Float
            if (DebugLevel > 1) write(*,'(A,EN20.5)') &
            & 'Half time is given in second --> conversion to seconds with factor ',C
         else if (trim(Units) == 'ms') then
            C = 1.E-3_Float
            if (DebugLevel > 1) write(*,'(A,EN20.5)') &
            & 'Half time is given in millisecond --> conversion to seconds with factor ',C
         else if (trim(Units) == 'us') then
            C = 1.E-6_Float
            if (DebugLevel > 1) write(*,'(A,EN20.5)') &
            & 'Half time is given in microsecond --> conversion to seconds with factor ',C
         else
            write(*,'(A)') 'Strange unit for halftime: '//Units//' !   --> Exiting!'
            call exit()
         endif

         NuclideSpecs(iNuclide)%HalfTime = NuclideSpecs(iNuclide)%HalfTime*C

         NuclideSpecs(iNuclide)%AtomNumber = GetAtomNumber(NuclideSpecs(iNuclide)%NuclideName)
         NuclideSpecs(iNuclide)%AtomName = AtomSpecs(NuclideSpecs(iNuclide)%AtomNumber)%Name
         if (DebugLevel > 1) then
            if (NuclideSpecs(iNuclide)%AtomNumber > 0) then
               write(*,'(A,I3,A,A)') 'Recognized atom ',NuclideSpecs(iNuclide)%AtomNumber,' : ',&
               & trim(NuclideSpecs(iNuclide)%AtomName)
            else
               write(*,'(A)') 'This is not a nuclide related to an atom in my database...'
            endif
         endif

         NuclideSpecs(iNuclide)%NuclideGroup = AtomSpecs(NuclideSpecs(iNuclide)%AtomNumber)%NuclideGroup
         !
         ! Extract daughters
         !
         do iDaughter=1,4
            if (DebugLevel > 2) write(*,'(I1,5X,A)') iDaughter,'"'//trim(YieldStr(iDaughter))//'"'
            if (len_trim(NuclideSpecs(iNuclide)%DaughterName(iDaughter)) > 0) &
            & read(YieldStr(iDaughter),'(E11.0)') NuclideSpecs(iNuclide)%DaughterFraction(iDaughter)
         enddo

         if (DebugLevel > 1) then
            write(*,'(A)') 'Nuclide: '//trim(NuclideSpecs(iNuclide)%NuclideName)
            write(*,'(A,EN20.10,A)') 'Half time: ',NuclideSpecs(iNuclide)%HalfTime,' [s]'
            write(*,'(A,A2)') 'Units of half time: ',Units
            do iDaughter=1,4
               if (len_trim(NuclideSpecs(iNuclide)%DaughterName(iDaughter)) > 0) then
               write(*,'(A,I1,A,A7,A,EN20.5)') 'Daughter(',iDaughter,'): ',&
                  & trim(NuclideSpecs(iNuclide)%DaughterName(iDaughter)),&
                  & '     Yield: ',NuclideSpecs(iNuclide)%DaughterFraction(iDaughter)
               endif
            enddo
            write(*,*)
         endif

         if (DebugLevel > 0) then
            write(*,'(A10,EN20.10)') trim(NuclideSpecs(iNuclide)%NuclideName),NuclideSpecs(iNuclide)%HalfTime
         endif

      enddo ! loop over ICRP-file

      if(DebugLevel > 1) write(*,*)

      close(ScratchFile)
   end subroutine ReadICRPNuclides



   subroutine MatchDaughters
      !
      ! Match the names of the daughters with the other nuclides.
      ! Exclude stable daughters
      !
      integer :: iNuclide,iDaughter
      integer, parameter :: DebugLevel = 0

      IsRelated = .FALSE.

      do iNuclide = 1,NNuclides
         !
         ! By default the aughters is set to be irrelevant...
         !
         do iDaughter = 1,MaxNDaughters
            NuclideSpecs(iNuclide)%Daughter(iDaughter) = &
            & GetNuclideNumber(NuclideSpecs(iNuclide)%DaughterName(iDaughter))
            if (NuclideSpecs(iNuclide)%Daughter(iDaughter) /= 0) then
               IsRelated(NuclideSpecs(iNuclide)%Daughter(iDaughter),iNuclide) = .TRUE.
            endif
         enddo
         !
         ! Print results
         !
         if (DebugLevel > 0) then
            write(*,40) iNuclide,&
               &    '"'//trim(NuclideSpecs(iNuclide)%NuclideName)//'"',&
               &    NuclideSpecs(iNuclide)%AtomNumber,&
               &    '"'//trim(NuclideSpecs(iNuclide)%AtomName)//'"',&
               &    NuclideSpecs(iNuclide)%HalfTime,&
               &    ('"'//trim(NuclideSpecs(NuclideSpecs(iNuclide)%Daughter(iDaughter))%NuclideName)//'"',&
               &    NuclideSpecs(iNuclide)%DaughterFraction(iDaughter),iDaughter=1,4)
            40    format(I4,'=',A12,'; atom(',I3,')=',A17,' T1/2=',G15.10,&
               & '  --> ',4(A12,' frac=',F10.8,'   '))

            write(*,*) 'Ready reading nuclide!'
            write(*,*) '-----------------------------------------------'
            write(*,*)
            write(*,*)
         endif
      enddo ! loop over nuclides
   end subroutine MatchDaughters



   subroutine ReadNuclideSpecs(InName)
      !
      ! Read the names, decaytimes and daughters for a collection of nuclides
      !
      character(*), intent(in) :: InName
      !
      if (NAtoms == 0) call ReadAtoms('atoms.dat')
      !
      ! Read the datafile
      !
      if (.NOT.FileExists(trim(InName))) then
         write(*,'(A)') 'Cannot find ICRP nuclide file '//'/'//trim(InName)
         write(*,'(A)')  'You can download this file via the Supplemental Material associated with '&
            & //'http://www.icrp.org/publication.asp?id=ICRP%20Publication%20107'
         write(*,'(A)')  'Exiting!'
         call exit()
      endif

      write(*,'(A)') 'Going to read nuclide specs from ICRP nuclide file '//trim(InName)

      call ReadICRPNuclides(InName)

      write(*,'(A)') 'Ready reading nuclide specs, going to match daughters...'
      !
      ! Match the names of the daughters with the other nuclides.
      ! Exclude stable daughers
      !
      call MatchDaughters()
      write(*,'(A)') 'Ready matching daughters!'
   end subroutine ReadNuclideSpecs



   subroutine MassNuc2NucMass(MyName)
      !
      ! When applicable, convert the nuclide name from format 137mCs to Cs137m
      !
      character(*), intent(inout) :: MyName
      logical :: HasWrongFormat,IsMetastable,Ready
      character(10) :: MassName,NukeName
      integer :: iCharacter

      HasWrongFormat = CharacterIsADigit(MyName(1:1))

      if (HasWrongFormat) then
      !
      ! Get the mass number
      !
      MassName = ' '
      iCharacter = 1
      Ready = .FALSE.
      do while (.NOT.Ready)
         MassName = trim(MassName)//MyName(iCharacter:iCharacter)
         iCharacter = iCharacter + 1
         Ready = (.NOT.CharacterIsADigit(MyName(iCharacter:iCharacter)))
      enddo ! not ready
      !
      ! Get the optional 'm' for meta-stable nuclides (case sensitive!!)
      !
      IsMetastable = (MyName(iCharacter:iCharacter) == 'm') ! This is case sensitive! Make sure you haven't converted to full uppercase yet!!

      if (IsMetastable) iCharacter = iCharacter + 1
      !
      ! Get the nuclide name
      !
      NukeName = ' '
      Ready = .FALSE.
      do while (.NOT.Ready)
         NukeName = trim(NukeName)//MyName(iCharacter:iCharacter)
         iCharacter = iCharacter + 1
         Ready = (MyName(iCharacter:iCharacter) == ' ')
      enddo ! not ready
      !
      ! Synthesis of elements in order. E.g., Cs137m
      !
      if (IsMetastable) then
         MyName = trim(NukeName)//trim(MassName)//'m'
      else
         MyName = trim(NukeName)//trim(MassName)
      endif

      endif
   end subroutine MassNuc2NucMass



   subroutine EnsureHyphen(MyName)
      !
      ! When applicable, add a hyphen as separator between atom name and mass number
      !
      character(*), intent(inout) :: MyName
      logical :: HasWrongFormat,Ready,HasHyphen
      character(10) :: NukeName
      integer :: iCharacter

      HasWrongFormat = CharacterIsADigit(MyName(1:1))

      if (HasWrongFormat) call MassNuc2NucMass(MyName)

      HasHyphen = (INDEX(MyName,'-') > 0)

      if (.NOT.HasHyphen) then
         !
         ! Get the mass number
         !
         NukeName = ' '
         iCharacter = 1
         Ready = .FALSE.
         do while (.NOT.Ready)
            NukeName = trim(NukeName)//MyName(iCharacter:iCharacter)
            iCharacter = iCharacter + 1
            Ready = CharacterIsADigit(MyName(iCharacter:iCharacter))
         enddo ! not ready

         NukeName = trim(NukeName)//'-'//MyName((iCharacter):len_trim(MyName))
         MyName = trim(NukeName)
      endif
   end subroutine EnsureHyphen



   subroutine ReadAtoms(InName)
      !
      ! Read the specs of all known atoms
      !
      character(*), intent(in) :: InName

      integer :: iAtom,DumInt,iNuclideGroup
      character(DefaultLength) :: ALine,DumStr
      logical :: Ready

      integer, parameter :: DebugLevel = 0

      write(*,'(A)') 'Going to read atom specs...'
      open(ScratchFile,file=(' ./data/')//trim(InName),form='FORMATTED',status='OLD',&
         & position='REWIND')
      read(ScratchFile,*) ! Skip header

      Ready = .FALSE.
      iAtom = 0
      do while (.NOT.Ready)
         read(ScratchFile,'(A)',end=20) ALine
         iAtom = iAtom + 1
         if (DebugLevel > 1) write(*,70) iAtom
         70 format('Interpreting specs for Atom ',I3)
         if (DebugLevel > 1) write(*,60) trim(ALine)
         60 format('Specs from database are: "',A,'"')
         !
         ! Get the atom specs
         !
         read(ALine( 1: 3),'(I3)')    DumInt
         read(ALine( 7:19),'(A13)')   AtomSpecs(iAtom)%Name
         read(ALine(21:23),'(A3)')    AtomSpecs(iAtom)%Symbol
         read(ALine(27:27),'(I1)')    AtomSpecs(iAtom)%Period
         read(ALine(34:35),'(I2)')    AtomSpecs(iAtom)%Group
         read(ALine(41:64),'(A24)')   AtomSpecs(iAtom)%ChemicalSeries
         read(ALine(65:76),'(F12.0)') AtomSpecs(iAtom)%Mass
         DumStr = ALine(78:len(ALine))

         AtomSpecs(iAtom)%NuclideGroup = 0

         if (len_trim(DumStr) /= 0) then
            do iNuclideGroup = 0,NNuclideGroups
            if (trim(DumStr) == trim(NuclideGroupName(iNuclideGroup))) &
               & AtomSpecs(iAtom)%NuclideGroup = iNuclideGroup
            enddo
         endif

         if (DumInt /= iAtom) then
            write(*,'(A)') 'Atom number does not match line number! Exiting!!'
            call exit()
         endif

         if (DebugLevel > 0) then
            write(*,200) AtomSpecs(iAtom),NuclideGroupName(AtomSpecs(iAtom)%NuclideGroup)
            200  format('Specs: "',A13,'"',1X,'"',A3,'"',1X,I2,1X,I2,1X,I2,1X,'"',A24,'"',1X,F12.8,1X,A)
         endif
      enddo
      20 continue ! ready reading data

      NAtoms = iAtom

      write(*,'(A,I3)') '... ready reading specs for this many atoms: ',NAtoms
   end subroutine ReadAtoms



   subroutine InitLibENDF()
      !
      ! Prepare some settings before use of this library
      !
      character(DefaultLength) :: FName,ALine
      logical :: ValidLine,IsException
      integer :: iLine,MyAtomNumber

      integer, parameter :: DebugLevel = 0

      FName = trim(ENDFPath()) // '/decay.list'

      open(ScratchFile,file=FName,form='FORMATTED',action='READ')
      !
      ! First, scan the file to get atom names only
      !
      read(ScratchFile,*)
      read(ScratchFile,*)
      read(ScratchFile,*)

      do
         read(ScratchFile,'(A)',end=10) ALine
         ValidLine = (ALine(1:1) /= '-')
         if (ValidLine) then
            read(ALine(7:9),*) MyAtomNumber
            AtomName(MyAtomNumber) = ALine(11:12)
         endif ! valid line
      enddo

      10 continue
      !
      ! Now, scan again and get the nuclides
      !
      REWIND(ScratchFile)

      read(ScratchFile,*)
      read(ScratchFile,*)
      read(ScratchFile,*)

      NNuclides = 0
      iLine = 3

      do
         read(ScratchFile,'(A)',end=20) ALine

         ValidLine = (ALine(1:1) /= '-')

         if (ValidLine) then
            NNuclides = NNuclides + 1
            iLine = iLine + 1
            !
            ! By default, all nuclides are assumed to be stable. Only if decay data is found, then this is overruled
            !
            NuclideSpecs(NNuclides)%HalfTime = 1.E90_Float ! infinity --> stable
            NuclideSpecs(NNuclides)%IsOrphan = .FALSE.
            !
            ! Interpret line
            !
            if (DebugLevel > 1) then
               write(*,'(A,I4,A)') 'Line ',NNuclides,' = "'//trim(ALine)//'"'
            endif

            read(ALine(7:9),*) NuclideSpecs(NNuclides)%AtomNumber
            NuclideSpecs(NNuclides)%AtomName = ALine(11:12)

            read(ALine(14:16),*) NuclideSpecs(NNuclides)%NHadrons

            if (ALine(17:17) == ' ') then
               NuclideSpecs(NNuclides)%MetaStableMode = 0
            else if (ALine(17:17) == 'M') then
               NuclideSpecs(NNuclides)%MetaStableMode = 1
            else if (ALine(17:17) == 'N') then
               NuclideSpecs(NNuclides)%MetaStableMode = 2
            else if (ALine(17:17) == 'O') then
               NuclideSpecs(NNuclides)%MetaStableMode = 3
            else
               write(*,'(A,I0,A,I0,A,I0,A)') 'Unknown metastability mode "',NuclideSpecs(NNuclides)%MetaStableMode,&
               & '" on line ',iLine,' for nuclide ',NNuclides,'! Exiting!!'
               call exit()
            endif
            !
            ! ICRP-107 knows Pm-137m and ENDF knows Pm-137. These are probably the same nuclide.
            ! To allow the use of both datasets simultaneously, the ENDF nuclide is mapped on the ICRP-nuclide:
            ! The nuclide name gets an "m". All subsequent references in ENDF to Pm-137 are trapped and re-cast to Pm-137m.
            !
            IsException = ((NuclideSpecs(NNuclides)%AtomName == 'Pm').AND.(NuclideSpecs(NNuclides)%NHadrons == 137))
            if (IsException) NuclideSpecs(NNuclides)%MetaStableMode = 1

            NuclideSpecs(NNuclides)%NDaughters = 0

            write(NuclideSpecs(NNuclides)%NuclideName,'(A,A,I0)') &
               & trim(NuclideSpecs(NNuclides)%AtomName),'-',NuclideSpecs(NNuclides)%NHadrons

            NuclideSpecs(NNuclides)%NuclideName = trim(NuclideSpecs(NNuclides)%NuclideName)&
            & //MetastabilityName(NuclideSpecs(NNuclides)%MetaStableMode)

            if (DebugLevel > 0) write(*,'(I4,1X,A)') NNuclides,trim(NuclideSpecs(NNuclides)%NuclideName)
         endif ! valid line
      enddo

      20 continue
      close(ScratchFile)
   end subroutine InitLibENDF



   integer function GetNuclideNumber(NuclideName)
      !
      ! Search for the given name of a nuclide in the list of known nuclides and
      ! return its number
      !
      character(*), intent(in) :: NuclideName
      integer :: DumInt,iNuclide
      character(10) :: DumNuclideName

      DumInt = 0

      if (len_trim(NuclideName) > 0) then
         DumNuclideName = trim(NuclideName)
         call Capitalize(DumNuclideName)
         do iNuclide = 1,NNuclides
            if (trim(DumNuclideName) == trim(NuclideSpecs(iNuclide)%NuclideName)) DumInt = iNuclide
         enddo
      endif ! length of input string <> 0

      GetNuclideNumber = DumInt
   end function GetNuclideNumber



   integer function AtomNMass2Nuclide(AtomNumber,NHadrons,MetastableMode)
      !
      ! Find the nuclide with these properties
      !
      integer, intent(in) :: AtomNumber,NHadrons,MetastableMode
      integer :: DumInt,iNuclide

      DumInt = 0

      do iNuclide = 1,NNuclides
      if (      (AtomNumber     == NuclideSpecs(iNuclide)%AtomNumber)&
         & .AND.(NHadrons       == NuclideSpecs(iNuclide)%NHadrons)&
         & .AND.(MetastableMode == NuclideSpecs(iNuclide)%MetastableMode)) DumInt = iNuclide
      enddo

      AtomNMass2Nuclide = DumInt
   end function AtomNMass2Nuclide



   integer function Zap2Nuclide(Zap,iMetaStableMode)
      !
      ! Convert a ZAP-number to a nuclide index
      !
      integer, intent(in) :: Zap,iMetaStableMode

      integer :: iAtom,iMass,iMother

      iAtom = Zap/1000
      iMass = Zap - 1000*iAtom

      iMother = AtomNMass2Nuclide(iAtom,iMass,iMetastableMode)
      Zap2Nuclide = iMother
   end function Zap2Nuclide



   subroutine CheckRecordType(iLine,ALine)
      !
      ! Assess which type of record this is
      !
      integer, intent(in) :: iLine
      character(*), intent(in) :: ALine

      integer :: MatNumber,MFNumber,MTNumber

      integer, parameter :: DebugLevel = 0
      !
      ! Read the information at the end of the line to see what type of record this is
      !
      read(ALine(67:75),'(I4,I2,I3)') MatNumber,MFNumber,MTNumber

      if ((MFNumber == 1).AND.(MTNumber == 451)) then
         iRecordType(iLine) = iTextRecord
      else if ((MFNumber == 8).AND.(MTNumber == 457)) then
         iRecordType(iLine) = iDecayRecord
      endif

      if (DebugLevel > 0) then
         write(*,'(A,I4,A29,3(1X,I4,1X,I2,1X,I3,4X,A))') 'Line ',iLine,&
         & ' has a record of type '//RecordName(iRecordType(iLine)),MatNumber,MFNumber,MTNumber,'"'//trim(ALine)//'"'
      endif
   end subroutine CheckRecordType


   subroutine Line2Values(ALine,MySubString,MyFloat,MyInteger,SubStringIsInteger,DoingBullocks)
      !
      ! Read 6 values from a line and check if they are real or integer
      !
      character(*), intent(in) :: ALine
      integer, dimension(6), intent(out) :: MyInteger
      real(Float), dimension(6), intent(out) :: MyFloat
      logical, dimension(6), intent(out) :: SubStringIsInteger
      logical, intent(out) :: DoingBullocks
      character(11), dimension(6), intent(inout) :: MySubString

      integer :: iSubString,TheExponent,TheIndex
      character(DefaultLength) :: DumString
      real(Float) :: RelativeDifference

      integer, parameter :: DebugLevel = 0

      read(ALine,'(6A11)') MySubString

      MyFloat = 0._Float
      MyInteger = 0
      DoingBullocks = .FALSE.

      do iSubString = 1,6
         if (DebugLevel > 2) write(*,'(A,I0,A)') 'Substring ',iSubString,' = "'//MySubString(iSubString)//'"'
      enddo
      do iSubString = 1,6
         if (DebugLevel > 1) write(*,'(A,I0,A)') 'Substring ',iSubString,' = "'//MySubString(iSubString)//'"'

         TheIndex = MAX(INDEX(MySubString(iSubString),'.'),&
         & INDEX(MySubString(iSubString),'+'),INDEX(MySubString(iSubString),'-'))
         if (TheIndex > 0) then
            SubStringIsInteger(iSubString) = .FALSE.
            !
            ! Values can be with or without E in the exponential notation, and there can be spaces in the string:
            !
            if ((INDEX(MySubString(iSubString),'E') /= 0)&
            & .OR.((INDEX(MySubString(iSubString),'+') == 0).AND.(INDEX(MySubString(iSubString),'-') == 0)))then
               read(MySubString(iSubString),*) MyFloat(iSubString)
            else
               read(MySubString(iSubString)(1:(TheIndex-1)),*) MyFloat(iSubString)
               DumString = MySubString(iSubString)(TheIndex:11)
               call RemoveCharacter(DumString,' ')
               read(DumString,*) TheExponent
               MyFloat(iSubString) = MyFloat(iSubString) * 10._Float**TheExponent
            endif
            !
            ! If this is an integer in disguise, make it an integer (a real one! ;-) )
            !
            RelativeDifference = 0._Float
            if (abs(MyFloat(iSubString)) > 0._Float) RelativeDifference = &
               & abs(MyFloat(iSubString)-NINT(MyFloat(iSubString)))/abs(MyFloat(iSubString))

            if ((MyFloat(iSubString) <= 1.E6_Float)&
               & .AND.(RelativeDifference < 1.E-10_Float)) then
               SubStringIsInteger(iSubString) = .TRUE.
               MyInteger(iSubString) = NINT(MyFloat(iSubString))
            endif
         else
            SubStringIsInteger(iSubString) = .TRUE.
            DoingBullocks = .TRUE.
            read(MySubString(iSubString),*,end=20) MyInteger(iSubString)
            DoingBullocks = .FALSE.
         endif

         20  continue
         !
         ! Show substring
         !
         if (DebugLevel > 1) then
            if (SubStringIsInteger(iSubString)) then
               write(*,'(A,I0,A,I0)') 'Substring ',iSubString,' is an integer: ',MyInteger(iSubString)
            else
               write(*,'(A,I0,A,F15.5)') 'Substring ',iSubString,' is a real    : ',MyFloat(iSubString)
            endif
         endif ! debug
      enddo ! loop over substrings
   end subroutine Line2Values



   subroutine Parse1DataLine(iLine,iNuclide,ALine,iLineInChapter,ReadyReadingFile,ThisLineHasDecay,iDecay)
      !
      ! Try to find useful info about a known nuclide
      !
      integer, intent(in) :: iLine,iNuclide,iLineInChapter
      character(*), intent(in) :: ALine
      logical, intent(inout) :: ReadyReadingFile,ThisLineHasDecay
      integer, intent(inout) :: iDecay

      character(DefaultLength) :: MyDecayModeLine
      integer :: DaughterAtom,DaughterMass,MyDecayMode,MyMetastableMode,MyAtomStep,MyMassStep,&
         iParticle,iCompoundDecay,iDaughter
      logical, dimension(6) :: SubStringIsInteger
      integer, dimension(6) :: MyInteger
      real(Float), dimension(6) :: MyFloat
      real(Float) :: MyFraction
      logical :: Ready,IsValidDecay,DoingBullocks,IsException
      character(11), dimension(6) :: MySubString

      integer, parameter :: DebugLevel = 0

      DoingBullocks = .FALSE.
      ReadyReadingFile = .FALSE.
      ThisLineHasDecay = .FALSE.

      if (DebugLevel > 1) write(*,'(A,I0,A)') 'Line ',iLine,' : "'//trim(ALine)//'"'
      !
      ! Only interpret decay info
      !
      if (iRecordType(iLine) == iDecayRecord) then
         if (DebugLevel > 1) write(*,'(A,I4,A)') 'Decay information record ',iLineInChapter,' = "'//trim(ALine)//'"'
         !
         ! Read the 6 values and check if they are real/integer
         !
         call Line2Values(ALine,MySubstring,MyFloat,MyInteger,SubStringIsInteger,DoingBullocks)
         ! Check this criterion!!!!!

         ReadyReadingFile =      SubStringIsInteger(1).AND.(MyInteger(1) == 0)&
            &                  .AND.SubStringIsInteger(2).AND.(MyInteger(2) == 0)&
            &                  .AND.SubStringIsInteger(3)&
            &                  .AND.SubStringIsInteger(4).AND.(MyInteger(4) == 0)&
            &                  .AND.SubStringIsInteger(5)&
            &                  .AND.SubStringIsInteger(6)&
            &                  .AND.((iLineInChapter == 2).OR.(iLineInChapter >= 5))

         ReadyReadingFile = ReadyReadingFile.OR.&
         & ((iLineInChapter >= 5).AND.(iDecay >= NuclideSpecs(iNuclide)%NDecayModes))

         if (DoingBullocks) then
            if (DebugLevel > 0) then
               write(*,'(A,I0)') 'Encountering compliance problems in file for '&
               & //NuclideSpecs(iNuclide)%NuclideName//' at line ',iLine
            endif ! debug
         else if (ReadyReadingFile) then
            if (DebugLevel > 1) then
               write(*,'(A,I0)') 'Ready reading file for '&
               & //NuclideSpecs(iNuclide)%NuclideName//' at line ',iLine
            endif ! debug
         else
            if (iLineInChapter == 2) then
               NuclideSpecs(iNuclide)%HalfTime = MyFloat(1) ! Even in cases where you have value e.g. 1.000, this goes well

               if (DebugLevel > 1) write(*,'(A,EN15.5,A)') 'Halflife = ',NuclideSpecs(iNuclide)%HalfTime,' s'
            else if (iLineInChapter == 4) then
               NuclideSpecs(iNuclide)%NDecayModes = MyInteger(6)

               if (DebugLevel > 1) write(*,'(A,I0)') 'Number of decay modes: ',NuclideSpecs(iNuclide)%NDecayModes
            else if (iLineInChapter >= 5) then

               MyMetastableMode = MyInteger(2)
               MyFraction = MyFloat(5)
               !
               ! If you have a non-trivial decay fraction, update nuclide database
               !
               if ((MyFraction > 0._Float).AND.(MyFraction <= 1._Float)) then
                  if (DebugLevel > 1) write(*,'(A,EN15.5)') 'Acceptable fraction = ',MyFraction

                  iParticle = 1
                  MyAtomStep = 0
                  MyMassStep = 0
                  IsValidDecay = .FALSE.
                  MyDecayModeLine = ' '
                  iCompoundDecay = NINT(1000000._Float*MyFloat(1))

                  Ready = (iCompoundDecay < 1000000)

                  do while (.NOT.Ready)
                     MyDecayMode = iCompoundDecay/10**(7-iParticle)

                     iCompoundDecay = iCompoundDecay - MyDecayMode*10**(7-iParticle)

                     Ready = (iCompoundDecay == 0)

                     if (iCompoundDecay < 0) then
                        write(*,'(A,I0)') 'Faulty value for iCompoundDecay: ',iCompoundDecay
                        call exit()
                     endif ! rubbish

                     if (DebugLevel > 1) write(*,'(A,A)') 'Found decay mode ',DecayModeName(MyDecayMode)

                     MyDecayModeLine = trim(MyDecayModeLine)//' + '//DecayModeName(MyDecayMode)

                     if (     (MyDecayMode >= 0)&
                        &   .AND.(MyDecayMode <= 7)) then
                        IsValidDecay = .TRUE.
                        MyAtomStep = MyAtomStep + AtomStep(MyDecayMode)
                        MyMassStep = MyMassStep + MassStep(MyDecayMode)
                     else
                        IsValidDecay = .FALSE.
                        if (DebugLevel > 1) write(*,'(A,A,A)') &
                        & 'Found disqualifying decay mode ',DecayModeName(MyDecayMode),', making this decay invalid!!'
                     endif ! valid decay mode

                     iParticle = iParticle + 1

                  enddo ! ready checking all particles in this single decay

                  MyDecayModeLine = MyDecayModeLine(4:DefaultLength)
                  !
                  ! Account for the full train of released particles in this decay, e.g. NN, in 1 go
                  !
                  if (IsValidDecay) then
                     ThisLineHasDecay = .TRUE.
                     iDecay = iDecay + 1
                     NuclideSpecs(iNuclide)%DecayFraction(iDecay) = MyFraction

                     if (.NOT.(MyDecayMode == iDecayModeSpontaneousFission)) then
                        NuclideSpecs(iNuclide)%NDaughters = NuclideSpecs(iNuclide)%NDaughters + 1
                        iDaughter = NuclideSpecs(iNuclide)%NDaughters

                        DaughterAtom = NuclideSpecs(iNuclide)%AtomNumber + MyAtomStep
                        DaughterMass = NuclideSpecs(iNuclide)%NHadrons   + MyMassStep

                        write(NuclideSpecs(iNuclide)%DaughterName(iDaughter),'(A,A,I0)') &
                           &  trim(AtomName(DaughterAtom)),'-',DaughterMass
                        !
                        ! Trap Pm-137 and recast it to Pm-137m, which is the only isotope in ICRP-107 with this combination.
                        ! To allow for the simultaneous use of ENDF with ICRP-107, Pm-137 is renamed Pm-137m
                        !
                        IsException = ((DaughterAtom == 61).AND.(DaughterMass == 137))
                        if (IsException) MyMetastableMode = 1

                        NuclideSpecs(iNuclide)%DaughterName(iDaughter) = &
                           & trim(NuclideSpecs(iNuclide)%DaughterName(iDaughter))&
                           & //MetaStabilityName(MyMetastableMode)

                        NuclideSpecs(iNuclide)%Daughter(iDaughter) = &
                           & GetNuclideNumber(NuclideSpecs(iNuclide)%DaughterName(iDaughter))

                        NuclideSpecs(iNuclide)%DaughterFraction(iDaughter) = MyFraction
                        NuclideSpecs(iNuclide)%DecayName(iDaughter) = trim(MyDecayModeLine)

                        if (DebugLevel > 0) write(*,'(A,I4,A,A7,A,EN15.5,A,A,A7,A,F20.10,A)') &
                           & 'Nuclide ',iNuclide,': ',&
                           & NuclideSpecs(iNuclide)%NuclideName,&
                           & ' with halflife ',NuclideSpecs(iNuclide)%HalfTime,' s',&
                           & ' to daughter ',NuclideSpecs(iNuclide)%DaughterName(iDaughter),&
                           & ' takes a fraction ',MyFraction,&
                           & ', decay: '//trim(NuclideSpecs(iNuclide)%DecayName(iDaughter))

                     else
                        if (DebugLevel > 0) then
                           write(*,'(A,I4,A,A7,A,EN15.5,A,A,F20.10,A)') &
                              & 'Nuclide ',iNuclide,': ',&
                              & NuclideSpecs(iNuclide)%NuclideName,&
                              & ' with halflife ',NuclideSpecs(iNuclide)%HalfTime,&
                              & ' s to many daughters  ',&
                              & ' takes a fraction ',MyFraction,&
                              & ', decay: Spontaneous fission'
                        endif
                     endif ! SF or other valid decay
                  endif ! valid decay mode
               else
                  if (DebugLevel > 1) write(*,'(A,EN15.5)') 'Unacceptable fraction = ',MyFraction
               endif ! yield > 0
            endif ! line 2 or 5 and later in the decay chapter
         endif ! not doing bullocks...
      endif ! Inside chapter of decay records
   end subroutine Parse1DataLine



   subroutine ParseENDFFile(iNuclide)
      !
      ! Parse 1 of the 3821 ENDF datafiles
      !
      integer, intent(in) :: iNuclide

      character(DefaultLength) :: FName,ALine
      character(20) :: MyTail
      logical :: Ready,IsNewChapter,FoundMeaningfulDecay,ThisLineHasDecay,IsException
      integer :: iLine,iLineInChapter,iDecay

      integer, parameter :: DebugLevel = 0
      !
      ! Construct filename
      !
      FName = 'dec-???_'
      write(FName(5:7),'(I3.3)') NuclideSpecs(iNuclide)%AtomNumber
      FName = trim(FName)//trim(NuclideSpecs(iNuclide)%AtomName)
      MyTail = '_???'
      write(MyTail(2:4),'(I3.3)') NuclideSpecs(iNuclide)%NHadrons
      !
      ! Trap Pm-137 and recast it to Pm-137m, which is the only isotope in ICRP-107 with this combination.
      ! To allow for the simultaneous use of ENDF with ICRP-107, Pm-137 is renamed Pm-137m
      !
      IsException = ((NuclideSpecs(iNuclide)%AtomNumber == 61).AND.(NuclideSpecs(iNuclide)%NHadrons == 137))



      if ((NuclideSpecs(iNuclide)%MetaStableMode /= 0).AND..NOT.IsException) then
         write(MyTail(5:6),'(A1,I1)') 'm',NuclideSpecs(iNuclide)%MetaStableMode
      endif
      FName = trim(FName)//trim(MyTail)

      FName = trim(ENDFPath()) // '/' // trim(FName)//'.endf'

      if (IsException) write(*,'(A)') 'To match ENDF with ICRP-107, instead of '&
      & //trim(NuclideSpecs(iNuclide)%NuclideName)&
      & //' I will read data from '//trim(FName)
      if (.NOT.FileExists(FName)) then
         write(*,'(A,I0,A)') 'Cannot find file "'//trim(FName)//'" for nuclide ',&
         & iNuclide,': '//trim(NuclideSpecs(iNuclide)%NuclideName)//'! Exiting!!'
         call exit()
      endif
      !
      ! Read contents
      !
      open(ScratchFile,file=trim(FName),form='FORMATTED',action='READ')

      Ready = .FALSE.
      iLine = 0
      IsContinuation = .FALSE.
      iRecordType = iUnknownRecord
      iLineInChapter = 0
      IsNewChapter = .TRUE.
      FoundMeaningfulDecay = .FALSE.
      iDecay = 0

      do while (.NOT.Ready)
         read(ScratchFile,'(A)',end=10) ALine
         iLine = iLine + 1
         !    write(*,'(A,I4,A)') 'Line(',iLine,') = "'//trim(ALine)//'"'
         !
         ! Find the type of the record
         !
         call CheckRecordType(iLine,ALine)

         if (iLine > 1) then
            IsNewChapter = (iRecordType(iLine) /= iRecordType(iLine-1))
         endif

         if (IsNewChapter) then
            iLineInChapter = 1
         else
            iLineInChapter = iLineInChapter + 1
         endif
         !
         ! Interpret the record on the line
         !
         call Parse1DataLine(iLine,iNuclide,ALine,iLineInChapter,Ready,ThisLineHasDecay,iDecay)

         FoundMeaningfulDecay = FoundMeaningfulDecay .OR. ThisLineHasDecay

         if (Ready.AND.(.NOT.FoundMeaningfulDecay)) then
            if (NuclideSpecs(iNuclide)%HalfTime < 1.E90_Float) then
               write(*,'(A,I4,A,A,A,EN20.10,A)') 'Nuclide ',iNuclide,': ',NuclideSpecs(iNuclide)%NuclideName,&
               & ' has halflife ',NuclideSpecs(iNuclide)%HalfTime,&
               & ', but I found end of file before I encountered any relevant decay info...'
            else
               if (DebugLevel > 0) write(*,'(A,I4,A,A,A)') 'Nuclide ',iNuclide,': ',&
               & NuclideSpecs(iNuclide)%NuclideName,' is stable!'
            endif
         endif

      enddo ! loop until ready

      10 continue

      close(ScratchFile)
   end subroutine ParseENDFFile



   subroutine ParseSpontaneousFile(FName)
      !
      ! Parse an ENDF file with branching ratio's for spontaneous fission.
      ! Unknown yet is if the branching ratio's account for the yield or not...
      !
      character(*), intent(in) :: FName

      logical :: Ready,IsSFLine,DoingBullocks,PrevLineIsSF,IsNonVoid,LastFissionProductFound
      integer :: iLine,SFLine,iMother,NExpectedFissionProducts,NFissionProducts,iElement
      character(DefaultLength) :: ALine
      character(11), dimension(12) :: MySubString
      logical, dimension(12) :: SubStringIsInteger
      integer, dimension(12) :: MyInteger
      integer, dimension(3) :: SFNuclide
      real(Float), dimension(12) :: MyFloat
      real(Float) :: TotalYield
      real(Float), dimension(3) :: MyBranchingRatio

      character(DefaultLength), parameter :: ENDFSpontaneousPath = './build/external/ENDFB_spontaneous/'

      if (.NOT.FileExists(trim(ENDFSpontaneousPath)//trim(FName))) then
         write(*,'(A)') 'Cannot find file "'//trim(FName)//'" for spontaneous fission! Exiting!!'
         call exit()
      else
         write(*,'(A)') 'Parsing file "'//trim(FName)//'" for spontaneous fission!'
      endif
      !
      ! Read contents
      !
      open(ScratchFile,file=trim(ENDFSpontaneousPath)//trim(FName),form='FORMATTED',action='READ')

      Ready = .FALSE.
      iLine = 0
      SFLine = 0
      TotalYield = 0._Float
      PrevLineIsSF = .FALSE.
      NFissionProducts = 0
      LastFissionProductFound = .FALSE.

      do while (.NOT.Ready)

         read(ScratchFile,'(A)',end=10) ALine
         iLine = iLine + 1

         IsSFLine = (ALine(72:75) == '8454')

         if (.NOT.IsSFLine) then
            if (PrevLineIsSF) then
               write(*,'(A,A10,A,F15.10)') 'Sum of branching ratios for spontaneous fission of ',&
               & NuclideSpecs(iMother)%NuclideName,' is ',TotalYield
               write(*,'(A,A,A,I6,A,I6,A)') 'For ',NuclideSpecs(iMother)%NuclideName,' I expected ',&
               & NExpectedFissionProducts,' fission products, found ',&
               & NFissionProducts,' fission products!'
            endif
            SFLine = 0 ! Reset Spontaneous Fission lines
            TotalYield = 0._Float
            PrevLineIsSF = .FALSE.
            NFissionProducts = 0
            LastFissionProductFound = .FALSE.
         else
            SFLine = SFLine + 1
            !
            ! Read the 6 values and check if they are real/integer
            !
            call Line2Values(ALine,MySubString,MyFloat,MyInteger,SubStringIsInteger,DoingBullocks)

            if (SFLine == 1) then
               iMother = Zap2Nuclide(MyInteger(1),MyInteger(2))

               write(*,'(A,I7,A,I4,A,A)') 'Head of chain for this list of spontaneous '&
                  & //'fissions, starting at line ',iLine,' = nuclide ',&
                  & iMother,': ',NuclideSpecs(iMother)%NuclideName

            else if (SFLine == 2) then
               NExpectedFissionProducts = MyInteger(6)
            else if (SFLine > 2) then
               read(ScratchFile,'(A)',end=10) ALine
               iLine = iLine + 1

               call Line2Values(ALine,MySubstring(7),MyFloat(7),MyInteger(7),SubStringIsInteger(7),DoingBullocks)
               !
               ! Work out 3 nuclides per 2 lines:
               !
               do iElement = 1,3
                  if (.NOT.LastFissionProductFound) then
                     IsNonVoid = (len_trim(MySubstring(4*iElement-3)) > 0)
                     if (IsNonVoid) then
                        SFNuclide(iElement) = Zap2Nuclide(MyInteger(4*iElement-3),MyInteger(4*iElement-2))
                        MyBranchingRatio(iElement) = MyFloat(4*iElement-1)

                        TotalYield = TotalYield + MyBranchingRatio(iElement)
                        NFissionProducts = NFissionProducts + 1

                        if ((SFNuclide(iElement) > 0).AND.(MyBranchingRatio(iElement) > 0._Float)) then

                           write(*,'(A,I8,A,A,A,A,A,EN20.10,1X,I3,1X,I3)') 'Line ',iLine,' SF: ',&
                           & NuclideSpecs(iMother)%NuclideName,' --> ',&
                           & NuclideSpecs(SFNuclide(iElement))%NuclideName,&
                           & ' with branching ratio ',MyBranchingRatio(iElement),&
                           & NuclideSpecs(SFNuclide(iElement))%AtomNumber,NuclideSpecs(SFNuclide(iElement))%NHadrons
                        endif

                        LastFissionProductFound = (NFissionProducts >= NExpectedFissionProducts)
                     endif ! nonvoid
                  endif ! not yet last fission product found
               enddo ! loop over elements on 2 lines

               PrevLineIsSF = .TRUE.

            endif ! SF line 1, 2 or later

         endif ! line with SF data
      enddo ! loop until ready

      10 continue
   end subroutine ParseSpontaneousFile



   subroutine MakeMotherDaughterMatrix(MyNuclideSpecs,MyMotherDaughterMatrix,OutName)
      !
      ! Fill the sparse matrix of which the exponential gives the evolution of the
      ! nuclide strengths (in Bq possibly per unit of volume or surface)
      ! over a macroscopic timestep. The matrix is stored in a global variable,
      ! since it will be needed quite often.
      !
      type(NuclideType), dimension(0:MaxNuclides), intent(in) :: MyNuclideSpecs
      real(Float), dimension(MaxNuclides,MaxNuclides), intent(out) :: MyMotherDaughterMatrix
      character(*), intent(in) :: OutName
      integer :: DaughterNuclide,MotherNuclide,NNonZero,iDaughter
      character(30) :: FormatString
      real(Float) :: Ln2

      integer, parameter :: DebugLevel = 0

      if (DebugLevel > 0) then
         write(*,'(A,I4)') 'MakeMyMotherDaughterMatrix: NNuclides = ',NNuclides
      endif

      ln2 = LOG(2.)
      !
      ! Initialize diagonal with -1, rest at zero
      !
      do MotherNuclide = 1,NNuclides
         do DaughterNuclide = 1,NNuclides
            if (DaughterNuclide == MotherNuclide) then
               MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide) = -ln2&
                  & /MyNuclideSpecs(DaughterNuclide)%HalfTime
            else
               MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide) =  0.
            endif
         enddo
      enddo
      !
      ! Now fill the nonzero off-diagonal elements. Leave open the option that some daughters are accounted for several times.
      !
      do MotherNuclide = 1,NNuclides
         do iDaughter = 1,MaxNDaughters
            DaughterNuclide = MyNuclideSpecs(MotherNuclide)%Daughter(iDaughter)
            if (DaughterNuclide /= 0) then
               MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide) = &
               MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide) &
               & + MyNuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter)*ln2&
               &   /MyNuclideSpecs(DaughterNuclide)%HalfTime
            endif
         enddo
      enddo
      !
      ! Write to output
      !
      if (DebugLevel > 0) write(*,*)
      if (len_trim(OutName) > 0) then
         open(ScratchFile,file = OutName,form='FORMATTED', position='REWIND')
         write(ScratchFile,'(A)') 'Mother/daughter Matrix: horizontal nuclides decay into the vertical nuclides'
         write(ScratchFile,664) NNuclides
         664  format(I4,' nuclides are represented in this matrix. The matrix elements are for the decay of ACTIVITY!!')

         FormatString = '(11X,    (A10,1X))'
         write(FormatString(6:9),'(I4)') NNuclides
         write(ScratchFile,FormatString) (MyNuclideSpecs(MotherNuclide)%NuclideName,MotherNuclide = 1,NNuclides)

         FormatString = '(A10,1X,    (G10.4,1X))'
         write(FormatString(9:12),'(I4)') NNuclides
         do DaughterNuclide = 1,NNuclides
            write(ScratchFile,FormatString) MyNuclideSpecs(DaughterNuclide)%NuclideName,&
         & (MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide),MotherNuclide = 1,NNuclides)
         enddo
         !
         ! Count number of non-zero entries in decay matrix
         !
         NNonZero = 0
         do DaughterNuclide = 1,NNuclides
            do MotherNuclide = 1,NNuclides
               if (abs(MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide)) > 1.E-24) call Inc(NNonZero)
            enddo
         enddo
         !
         ! Write only non-zero entries to output
         !
         write(ScratchFile,*)
         write(ScratchFile,'(A)') 'Here come the non-zero entries of the mother/daughter matrix:'
         write(ScratchFile,668) NNuclides,NNonZero
         668  format(I10,1X,I10,'  = NNuclides and NNonZero')
         do DaughterNuclide = 1,NNuclides
            do MotherNuclide = 1,NNuclides
               if (abs(MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide)) > 1.E-24) then
                  write(ScratchFile,667) &
                     & DaughterNuclide,MotherNuclide,MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide)
                     667 format(I4,1X,I4,1X,G10.4)
               endif
            enddo
         enddo
         close(ScratchFile)
      endif
   end subroutine MakeMotherDaughterMatrix



   subroutine ReadENDFNuclideSpecs()
      !
      ! Get all nuclide stuff you want to have from the ENDF database
      !
      integer :: iNuclide

      call InitLibENDF()

      do iNuclide = 1,NNuclides
         call ParseENDFFile(iNuclide)
      enddo ! loop over nuclides
      !
      ! Spontaneous fission is skippped, as for none of these nuclides there are meaningful dose conversion coefficients.
      ! You can give it a go by uncommenting the commands below. No guarantee that it works.
      !
      !   call ParseSpontaneousFile('nt501')
      !   call ParseSpontaneousFile('nt502')
      !   call ParseSpontaneousFile('nt503')
      !   call ParseSpontaneousFile('nt505')
      !   call ParseSpontaneousFile('nt506')
      !   call ParseSpontaneousFile('nt507')
      !   call ParseSpontaneousFile('nt508')
      !   call ParseSpontaneousFile('nt509')
      !   call ParseSpontaneousFile('nt510')
      !   call ParseSpontaneousFile('nt511')
      !   call ParseSpontaneousFile('nt512')
      !   call ParseSpontaneousFile('nt513')
      !   call ParseSpontaneousFile('nt514')
      !   call ParseSpontaneousFile('nt515')
      !   call ParseSpontaneousFile('nt516')
      !   call ParseSpontaneousFile('nt517')
      !   call ParseSpontaneousFile('nt518')
      !   call ParseSpontaneousFile('nt666a')
      !   call ParseSpontaneousFile('nt666b')
      !   call ParseSpontaneousFile('nt998')
      !   call ParseSpontaneousFile('nt999')
      !
      ! Match the names of the daughters with the other nuclides.
      ! Exclude stable daughers
      !
      call MatchDaughters()
      write(*,'(A)') 'Ready matching daughters!'
   end subroutine ReadENDFNuclideSpecs



   subroutine ReadNProcessENDFNuclideSpecs(tMin,RegularizedNuclideSpecs,RegularizedMotherDaughterMatrix)
      !
      ! Read ENDF nuclide database, find orphans and construct mother-daughter matrices
      !
      real(Float), intent(inout) :: tMin
      type(NuclideType), dimension(0:MaxNuclides), intent(inout) :: RegularizedNuclideSpecs
      real(Float), dimension(MaxNuclides,MaxNuclides), intent(inout) :: RegularizedMotherDaughterMatrix

      type(SparseMatrix) :: SparseMMMatrix,SparseRegularizedMMMatrix
      real(Float) :: tMinSug
      logical :: ENDFBinFileExists,IsDifferenttMin
      type(SparseLogicalMatrix) :: SparseIsRelated,SparseRegularizedIsRelated
      character(DefaultLength) :: FName

      integer, parameter :: DebugLevel = 0

      write(*,'(A)') 'Reading the ENDF database of nuclides...'

      if (.NOT.FileExists(RIVMSourcesPath() // '/ICRP-07.NDX')) then
         write(*,'(A)') 'Make sure you have file ICRP-07.NDX in the directory where you call this utility! Exiting!!'
         call exit()
      endif
      !
      ! Read all data on nuclides from an earlier run or construct them afresh
      !
      ENDFBinFileExists = FileExists(TransitionMatrixPath() // '/ENDFBinFile.dat')

      if (ENDFBinFileExists) then

         write(*,'(A)') 'I found the binary file ENDFBinFile.dat with all parsed ENDF data readily available!'
         write(*,'(A)') 'This is not the first time that you run this library!'
         write(*,'(A)') ''
         write(*,'(A)') 'Initialization will go smooth and fast!'
         write(*,'(A)') 'Going to read this preprocessed file ...'

         open(ScratchFile,file=TransitionMatrixPath() // '/ENDFBinFile.dat',form='UNFORMATTED',position='REWIND',action='READ')

         read(ScratchFile) AtomSpecs
         read(ScratchFile) AtomName
         read(ScratchFile) NNuclides
         read(ScratchFile) NuclideSpecs

         if (ALLOCATED(SparseMMMatrix%Element)) deallocate(SparseMMMatrix%Element)
         read(ScratchFile) SparseMMMatrix%N,SparseMMMatrix%NMax
         allocate(SparseMMMatrix%Element(SparseMMMatrix%N))
         read(ScratchFile) SparseMMMatrix%Element
         call SparseMatrix2Matrix(SparseMMMatrix,MotherDaughterMatrix)

         read(ScratchFile) SparseIsRelated%N
         allocate(SparseIsRelated%Element(SparseIsRelated%N))
         read(ScratchFile) SparseIsRelated%Element
         call SparseLMatrix2LMatrix(SparseIsRelated,IsRelated)

         read(ScratchFile) tMinSug
         IsDifferenttMin = (abs(tMinSug/tMin-1._Float) > 0.001_Float)
         if (IsDifferenttMin) then
            write(*,'(A,EN20.10,A,EN20.10,A)') 'The value of tMin in the file of ',tMinSug,&
            & ' seconds differs from the one in your subroutine call: ',tMin,' seconds!'
            write(*,'(A)') 'Going to replace your value by the one in tne file!'
            tMin = tMinSug
         endif

         read(ScratchFile) Orphanage

         read(ScratchFile) RegularizedNuclideSpecs

         if (ALLOCATED(SparseRegularizedMMMatrix%Element)) deallocate(SparseRegularizedMMMatrix%Element)
         read(ScratchFile) SparseRegularizedMMMatrix%N,SparseRegularizedMMMatrix%NMax
         allocate(SparseRegularizedMMMatrix%Element(SparseRegularizedMMMatrix%N))
         read(ScratchFile) SparseRegularizedMMMatrix%Element
         call SparseMatrix2Matrix(SparseRegularizedMMMatrix,RegularizedMotherDaughterMatrix)

         read(ScratchFile) SparseRegularizedIsRelated%N
         allocate(SparseRegularizedIsRelated%Element(SparseRegularizedIsRelated%N))
         read(ScratchFile) SparseRegularizedIsRelated%Element
         call SparseLMatrix2LMatrix(SparseRegularizedIsRelated,RegularizedIsRelated)

         close(ScratchFile)

         write(*,'(A)') 'Ready reading preprocessed ENDF file ...'
      else
         write(*,'(A)') 'I could not find the binary file ENDFBinFile.dat with all parsed ENDF data readily available...'
         write(*,'(A)') 'This is probably the first time that you run this library...'
         write(*,'(A)') ''
         write(*,'(A)') 'Going to read all the ENDF data files one by one and parse them ...'
         write(*,'(A)') 'This may take some time ...'

         call ReadENDFNuclideSpecs()
         !
         ! Find all orphan decay steps, i.e. those involving a halflife < tMin, and put them in a list
         !
         call FindOrphans(tMin)
         !
         ! Make a copy of the nuclide specs with which we can play
         ! The original is kept such that we will be able to trace the true family relations
         !
         RegularizedNuclideSpecs = NuclideSpecs
         !
         ! Regularize nuclide specs by elimination of too fast nuclides
         !
         call RegularizeNuclides(RegularizedNuclideSpecs)

         write(*,*)
         !================================================================================================================
         !
         ! Construct matrix with time-rate-of-change of nuclide vector, based on the regularized transition characteristics
         !
         write(*,*)
         write(*,'(A)') 'Going to construct mother-daughter matrix for time-rate-of-change!'
         !
         ! Once for the original full set of decay options
         !
         FName = ' '
         call MakeMotherDaughterMatrix(NuclideSpecs,MotherDaughterMatrix,FName)
         !
         ! And once for the regularized set of decay relations, where all orphans have been taken out
         !
         FName = ' '
         call MakeMotherDaughterMatrix(RegularizedNuclideSpecs,RegularizedMotherDaughterMatrix,FName)

         write(*,'(A)') 'Ready reading and parsting all the ENDF data files!'

         open(ScratchFile,file=TransitionMatrixPath() // '/ENDFBinFile.dat',form='UNFORMATTED',position='REWIND',action='WRITE')

         write(ScratchFile) AtomSpecs
         write(ScratchFile) AtomName
         write(ScratchFile) NNuclides
         write(ScratchFile) NuclideSpecs
         SparseMMMatrix = Matrix2SparseMatrix(MotherDaughterMatrix)
         write(ScratchFile) SparseMMMatrix%N,SparseMMMatrix%NMax
         write(ScratchFile) SparseMMMatrix%Element
         SparseIsRelated = LMatrix2SparseLMatrix(IsRelated)
         write(ScratchFile) SparseIsRelated%N
         write(ScratchFile) SparseIsRelated%Element
         write(ScratchFile) tMin
         write(ScratchFile) Orphanage
         write(ScratchFile) RegularizedNuclideSpecs
         SparseRegularizedMMMatrix = Matrix2SparseMatrix(RegularizedMotherDaughterMatrix)
         write(ScratchFile) SparseRegularizedMMMatrix%N,SparseRegularizedMMMatrix%NMax
         write(ScratchFile) SparseRegularizedMMMatrix%Element
         SparseRegularizedIsRelated = LMatrix2SparseLMatrix(RegularizedIsRelated)
         write(ScratchFile) SparseRegularizedIsRelated%N
         write(ScratchFile) SparseRegularizedIsRelated%Element

         close(ScratchFile)

         write(*,'(A)') 'Created binary file ENDFBinFile.dat with all parsed ENDF data for future use!'

      endif ! Binary data already available

      if (DebugLevel > 0) then
         write(*,'(A,I0,A)') 'AtomSpecs                  takes : ',SIZEOF(AtomSpecs),' bytes'
         write(*,'(A,I0,A)') 'AtomName                   takes : ',SIZEOF(AtomName),' bytes'
         write(*,'(A,I0,A)') 'NNuclides                  takes : ',SIZEOF(NNuclides),' bytes'
         write(*,'(A,I0,A)') 'NuclideSpecs               takes : ',SIZEOF(NuclideSpecs),' bytes'
         write(*,'(A,I0,A)') 'SparseMMMatrix             takes : ',SIZEOF(SparseMMMatrix%Element),' bytes'
         write(*,'(A,I0,A)') 'SparseIsRelated            takes : ',SIZEOF(SparseIsRelated%Element),' bytes'
         write(*,'(A,I0,A)') 'tMinSug                    takes : ',SIZEOF(tMinSug),' bytes'
         write(*,'(A,I0,A)') 'Orphanage                  takes : ',SIZEOF(Orphanage),' bytes'
         write(*,'(A,I0,A)') 'RegularizedNuclideSpecs    takes : ',SIZEOF(RegularizedNuclideSpecs),' bytes'
         write(*,'(A,I0,A)') 'SparseRegularizedMMMatrix  takes : ',SIZEOF(SparseRegularizedMMMatrix%Element),' bytes'
         write(*,'(A,I0,A)') 'SparseRegularizedIsRelated takes : ',SIZEOF(SparseRegularizedIsRelated%Element),' bytes'
      endif
   end subroutine ReadNProcessENDFNuclideSpecs



   subroutine FindOrphans(tMin)
      !================================================================================================================
      !
      ! Find the orphans: nuclides with halftime shorter than Orphanage%tMin
      !
      !================================================================================================================
      real(Float), intent(in) :: tMin

      integer :: DaughterNuclide,MotherNuclide,NMothers,iDaughter,HerDaughterNuclide,iOrphan

      integer, parameter :: DebugLevel = 0

      write(*,*)
      write(*,'(A)') 'Starting to look for orphan nuclides!'

      Orphanage%tMin = tMin
      !
      ! Find orphans in the original nuclide specs dataset (to prevent missing combos of Mother --> Orphanage%Orphan-daughter --> Orphanage%Orphan-granddaughter)
      !
      Orphanage%NOrphans = 0
      Orphanage%NTooUnStable = 0

      do DaughterNuclide = 1,NNuclides
         Orphanage%IsTooUnstable(DaughterNuclide) = (NuclideSpecs(DaughterNuclide)%HalfTime < Orphanage%tMin)
         if (Orphanage%IsTooUnstable(DaughterNuclide)) then
            if (DebugLevel > 0) write(*,'(A,EN20.10,A,EN20.10,A)') 'Handling too fast nuclide '&
               & //trim(NuclideSpecs(DaughterNuclide)%NuclideName)//' with halflife ',NuclideSpecs(DaughterNuclide)%HalfTime,&
               & ' [s] < ',Orphanage%tMin,'[s]'

            NuclideSpecs(DaughterNuclide)%IsOrphan = .TRUE.

            Orphanage%NTooUnstable = Orphanage%NTooUnStable + 1
            !
            ! Check all nuclides if they can be the mother of this daughter
            !
            NMothers = 0

            do MotherNuclide = 1,NNuclides
               do iDaughter = 1,MaxNDaughters
                  HerDaughterNuclide = NuclideSpecs(MotherNuclide)%Daughter(iDaughter)
                  !
                  ! If a possible mother has been found, then define a new orphan
                  !
               if (DaughterNuclide == HerDaughterNuclide) then
                  NMothers = NMothers + 1

                  if (DebugLevel > 0) write(*,'(A,I0,A,F15.10)') &
                     & 'Found mother '//trim(NuclideSpecs(MotherNuclide)%NuclideName)//&
                     & ' and '//trim(NuclideSpecs(DaughterNuclide)%NuclideName)//' is her daughter number ',iDaughter,&
                     & ' with branching ratio ',NuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter)

                  if (Orphanage%NOrphans == MaxNOrphans) then
                     write(*,'(A,I0,A)') 'Cannot add new orphan after maximum number of ',&
                        & Orphanage%NOrphans,' orphans! Exiting!!'
                     call exit()
                  endif ! maximum dimension reached

                  Orphanage%NOrphans = Orphanage%NOrphans + 1

                  Orphanage%Orphan(Orphanage%NOrphans)%Mother   = MotherNuclide
                  Orphanage%Orphan(Orphanage%NOrphans)%Daughter = DaughterNuclide
                  Orphanage%Orphan(Orphanage%NOrphans)%yield    = NuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter)
               endif ! Recognition of daughter by mother
            enddo ! loop over daughters
         enddo ! Loop over candidate mothers

         if (DebugLevel > 0) then
            if (NMothers > 0) then
               write(*,'(A,I0,A)') 'Found ',NMothers,' mothers for '//trim(NuclideSpecs(DaughterNuclide)%NuclideName)//'!'
            else
               write(*,'(A)') 'No mothers were found, '//trim(NuclideSpecs(DaughterNuclide)%NuclideName)&
                  & //' is an artificial nuclide and not a decay product!'
            endif

            write(*,'(A)') 'This is all I wanted to say about too fast nuclide '&
               & //trim(NuclideSpecs(DaughterNuclide)%NuclideName)
            endif ! DebugLevel > 0
         endif ! Halflife of daughter too short
      enddo ! Loop over daughters
      !
      ! Show list of orphans
      !
      write(*,'(A)') 'List of orphans: '

      do iOrphan = 1,Orphanage%NOrphans
         write(*,'(A,I4,A,A,A,A,5X,A,EN20.10,A,F20.10,A)') &
            & 'Orphanage%Orphan ',iOrphan,' is daughter ',&
            & NuclideSpecs(Orphanage%Orphan(iOrphan)%Daughter)%NuclideName,&
            & ' from mother ',&
            & NuclideSpecs(Orphanage%Orphan(iOrphan)%Mother  )%NuclideName,&
            & 'because halflife ',&
            & NuclideSpecs(Orphanage%Orphan(iOrphan)%Daughter)%HalfTime,&
            & '[s] < limit ',Orphanage%tMin,'[s]'
      enddo ! loop over orphans

      write(*,'(A)') 'Ready looking for orphan nuclides!'
   end subroutine FindOrphans



   subroutine RegularizeNuclides(MyNuclideSpecs)
      !
      ! All orphan nuclides with halftime shorter than Orphanage%tMin are eliminated.
      ! Their contributions are attributed to the mother nuclide decaying directly to the grand-daughters.
      !
      type(NuclideType), dimension(0:MaxNuclides), intent(inout) :: MyNuclideSpecs

      integer :: DaughterNuclide,NGrandDaughters,iGrandDaughter,GrandDaughterNuclide,NMothers,MotherNuclide,iDaughter,&
      & HerDaughterNuclide,FreeDaughter,iNuclide
      logical :: IsAlreadyADaughter,IsFirstFreeSlot

      integer, parameter :: DebugLevel = 0

      write(*,*)
      write(*,'(A,EN20.10,A)') 'Going to eliminate all nuclides with halflife faster than ',Orphanage%tMin,' s'
      write(*,*)
      !
      ! "DaughterNuclide" is the possibly too fast nuclide, i.e. an orphan nuclide, that has to be brushed under the carpet.
      !
      do DaughterNuclide = 1,NNuclides
         if (Orphanage%IsTooUnstable(DaughterNuclide)) then
            if (DebugLevel > 0) write(*,'(A,EN20.10,A)') 'Handling too fast nuclide '&
               & //trim(MyNuclideSpecs(DaughterNuclide)%NuclideName)//' with halflife ',&
               & MyNuclideSpecs(DaughterNuclide)%HalfTime,&
               & ' s'
            !
            ! Count number of granddaughters that is formed and that we have to re-arrange
            !
            write(*,'(A)') 'Going to search for radioactive granddaughters...'
            NGrandDaughters = 0
            do iGrandDaughter = 1,MaxNDaughters
               GrandDaughterNuclide = MyNuclideSpecs(DaughterNuclide)%Daughter(iGrandDaughter)
               if (GrandDaughterNuclide > 0) then
                  write(*,'(A)') 'Spotted radioactive granddaughter '&
                     & //MyNuclideSpecs(GrandDaughterNuclide)%NuclideName
                  NGrandDaughters = NGrandDaughters + 1
               endif
            enddo
            if (NGrandDaughters > 0) then
               write(*,'(A,I0,A)') 'This daughter nuclide forms ',NGrandDaughters,' radioactive granddaughters'
            else
               write(*,'(A,I0,A)') 'This daughter nuclide forms no radioactive granddaughters, '&
               & //'I will only eliminate the formation of this nuclide'
            endif
            !
            ! Check all nuclides if they can be the mother of this orphan daughter
            !
            ! "MotherNuclide" is the mother of the orphan, whose decay to the orphan will be shifted to her dranddaughters
            !
            NMothers = 0

            do MotherNuclide = 1,NNuclides
               do iDaughter = 1,MaxNDaughters
                  HerDaughterNuclide = MyNuclideSpecs(MotherNuclide)%Daughter(iDaughter)
                  !
                  !#############################################################################################################
                  !
                  ! If a possible mother has been found, eliminate the family relation and pass all info on to the next generation
                  !
                  if (DaughterNuclide == HerDaughterNuclide) then

                     NMothers = NMothers + 1

                     write(*,'(A,I0,A,F15.10)') 'Found mother '//trim(MyNuclideSpecs(MotherNuclide)%NuclideName)//&
                        & ' and '//trim(MyNuclideSpecs(DaughterNuclide)%NuclideName)//' is her daughter number ',&
                        & iDaughter,&
                        & ' with branching ratio ',MyNuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter)
                     !
                     ! Pass heritage to granddaughters via daughter
                     !
                     write(*,'(A)') 'Going to pass heritage of mother '&
                        & //trim(MyNuclideSpecs(MotherNuclide)%NuclideName)&
                        & //' to her radioactive granddaughters...'

                     do iGrandDaughter = 1,MyNuclideSpecs(DaughterNuclide)%NDaughters
                        GrandDaughterNuclide = MyNuclideSpecs(DaughterNuclide)%Daughter(iGrandDaughter)
                        if (GrandDaughterNuclide > 0) then
                           !
                           ! Found a grand-daughter
                           !
                           write(*,'(A,I0,A,F15.10)') 'Found granddaughter '&
                              & //trim(MyNuclideSpecs(GrandDaughterNuclide)%NuclideName)&
                              &   //' with index ',iGrandDaughter,' and with branching ratio ',&
                              & MyNuclideSpecs(DaughterNuclide)%DaughterFraction(iGrandDaughter)
                           !
                           ! Check if this granddaughter is already a daughter of the grandmother
                           !
                           IsAlreadyADaughter = .FALSE.
                           FreeDaughter = 0
                           do iNuclide = 1,MaxNDaughters
                              if (MyNuclideSpecs(MotherNuclide)%Daughter(iNuclide) == GrandDaughterNuclide) then
                                 FreeDaughter = iNuclide
                                 IsAlreadyADaughter = .TRUE.
                              endif ! Granddaughter is already also a daughter
                           enddo
                           !
                           ! If not already a daughter, find the first available new daughter slot
                           !
                           if (FreeDaughter == 0) then
                              do iNuclide = 1,MaxNDaughters
                                 IsFirstFreeSlot=((FreeDaughter == 0).AND.(MyNuclideSpecs(MotherNuclide)%Daughter(iNuclide) == 0))
                                 if (IsFirstFreeSlot) FreeDaughter = iNuclide
                              enddo
                           endif ! granddaughter was not also a daughter

                           if (FreeDaughter > 0) then
                              if (IsAlreadyADaughter) then
                                 write(*,'(A,I0,A)') 'Daughter slot ',FreeDaughter,&
                                 & ' of mother '//trim(MyNuclideSpecs(MotherNuclide)%NuclideName)&
                                 & //' already contained granddaughter '//trim(MyNuclideSpecs(GrandDaughterNuclide)%NuclideName)&
                                 & //'! Adding contribution!'

                                 write(*,'(A,I0,A)') 'Augmenting branching ratio for old daughter '&
                                 & //trim(MyNuclideSpecs(GrandDaughterNuclide)%NuclideName)&
                                 & //' of mother '//trim(MyNuclideSpecs(MotherNuclide)%NuclideName)//' in daughter slot number ',&
                                 & FreeDaughter,':'
                                 write(*,'(5(A,F15.10))') &
                                 & 'Adding ',MyNuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter),&
                                 & ' x ',MyNuclideSpecs(DaughterNuclide)%DaughterFraction(iGrandDaughter),&
                                 & ' = ',MyNuclideSpecs(MotherNuclide  )%DaughterFraction(iDaughter)&
                                 & *MyNuclideSpecs(DaughterNuclide)%DaughterFraction(iGrandDaughter),&
                                 & ' to ',MyNuclideSpecs(MotherNuclide)%DaughterFraction(FreeDaughter),&
                                 & ' which gives ',MyNuclideSpecs(MotherNuclide)%DaughterFraction(FreeDaughter) + &
                                 &  MyNuclideSpecs(MotherNuclide  )%DaughterFraction(iDaughter     )&
                                 & *MyNuclideSpecs(DaughterNuclide)%DaughterFraction(iGrandDaughter)

                                 MyNuclideSpecs(MotherNuclide)%DaughterFraction(FreeDaughter) = &
                                 & MyNuclideSpecs(MotherNuclide)%DaughterFraction(FreeDaughter) + &
                                 &  MyNuclideSpecs(MotherNuclide  )%DaughterFraction(iDaughter     )&
                                 & *MyNuclideSpecs(DaughterNuclide)%DaughterFraction(iGrandDaughter)

                              else
                                 write(*,'(A,I0,A)') 'Daughter slot ',FreeDaughter,&
                                 & ' of mother is still free, adopting it for granddaughter!'

                                 write(*,'(A,I0)') 'Installing new daughter '&
                                 & //trim(MyNuclideSpecs(GrandDaughterNuclide)%NuclideName)&
                                 & //' of mother '//trim(MyNuclideSpecs(MotherNuclide)%NuclideName)//' in daughter slot number ',&
                                 & FreeDaughter

                                 MyNuclideSpecs(MotherNuclide)%Daughter(FreeDaughter) = &
                                 & GrandDaughterNuclide

                                 MyNuclideSpecs(MotherNuclide)%DaughterName(FreeDaughter) = &
                                 & MyNuclideSpecs(DaughterNuclide)%DaughterName(iGrandDaughter)

                                 MyNuclideSpecs(MotherNuclide)%DaughterFraction(FreeDaughter) = &
                                 &  MyNuclideSpecs(MotherNuclide  )%DaughterFraction(iDaughter     )&
                                 & *MyNuclideSpecs(DaughterNuclide)%DaughterFraction(iGrandDaughter)

                              endif
                           else
                              write(*,'(A)') 'Somehow I cannot find a free daughter...'
                              call exit()
                           endif

                        endif
                     enddo ! loop over granddaughters
                     !
                     ! Eliminate the decay from this mother to the daughter at hand in the nuclide database:
                     !
                     write(*,'(A,I0,A)') 'Eliminating daughter slot ',iDaughter,' of mother!'

                     MyNuclideSpecs(MotherNuclide)%Daughter(iDaughter) = 0
                     MyNuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter) = 0._Float
                     MyNuclideSpecs(MotherNuclide)%DaughterName(iDaughter) = ' '

                  endif ! Recognition of daughter by mother
                  !
                  ! End of bookkeeping of the removal of 1 orphan nuclide by adding its contributions to 1 of its mothers
                  !
                  !#############################################################################################################
                  !
               enddo ! loop over 20 daughters
            enddo ! Loop over candidate mothers

            if (NMothers > 0) then
               write(*,'(A,I0,A)') 'Found ',NMothers,' mothers!'
            else
               write(*,'(A)') 'No mothers were found, this is an artificial nuclide and not a decay product!'
            endif

            write(*,'(A)')'Ready handling too fast nuclide '//trim(MyNuclideSpecs(DaughterNuclide)%NuclideName)
            write(*,*)
         endif ! Halflife of daughter too short
      enddo ! Loop over daughters
      !
      ! Make IsRelated(Daughter,Mother) matrix for the regularized set
      !
      RegularizedIsRelated = .FALSE.

      do MotherNuclide = 1,NNuclides
         do iDaughter = 1,MaxNDaughters
            HerDaughterNuclide = MyNuclideSpecs(MotherNuclide)%Daughter(iDaughter)
            if (HerDaughterNuclide /= 0) then
               RegularizedIsRelated(HerDaughterNuclide,MotherNuclide) = .TRUE.
            endif
         enddo
      enddo

      write(*,'(A)') 'Ready regularizing mother-daughter matrix!'
   end subroutine RegularizeNuclides



   subroutine CollectProgeny(iNuclide,DoPrint)
      !================================================================================================================
      !
      ! Collect all progeny of this nuclide, first the regular ones, then add the short-living ones.
      !
      !================================================================================================================
      integer, intent(in) :: iNuclide
      logical, intent(in) :: DoPrint

      integer :: NAdded,MotherNuclide,DaughterNuclide,OrphanCounter,SwiftNuclide,NFarAdded
      logical :: Ready,IsNewFamily,IsNewFarFamily

      integer, parameter :: DebugLevel = 0

      if (DoPrint) write(*,'(A,I0,A,EN20.10,A)') 'Collecting the progeny family for nuclide ',iNuclide,&
      & ': '//NuclideSpecs(iNuclide)%NuclideName//' with halflife ',NuclideSpecs(iNuclide)%HalfTime,' [s]:'

      NuclideFamily%NFamily = 1 ! i.e.: the mother nuclide
      NuclideFamily%NFarFamily = 1 ! i.e.: the mother nuclide
      NuclideFamily%FamilyMember = 0
      NuclideFamily%FamilyNumber = 0
      NuclideFamily%NMothers = 0
      NuclideFamily%NDaughters = 0
      NuclideFamily%Mother = 0
      NuclideFamily%Daughter = 0
      NuclideFamily%FamilyMember(1) = iNuclide
      NuclideFamily%FamilyNumber(iNuclide) = 1
      NuclideFamily%IsFamily = .FALSE.
      NuclideFamily%IsFamily(iNuclide) = .TRUE.
      NuclideFamily%IsFarFamily = .FALSE.
      NuclideFamily%IsFarFamily(iNuclide) = .TRUE.
      !
      ! Iterate until no new long living family member is added
      !
      Ready = .FALSE.

      do while (.NOT.Ready)
         NAdded = 0
         NFarAdded = 0
         !
         ! Loop over possible mothers
         !
         do MotherNuclide = 1,NNuclides
            !
            ! Loop over possible daughters
            !
            do DaughterNuclide = 1,NNuclides
               !
               ! First check if this nuclide is part of the full family tree, including orphans
               ! The search is downward in the tree, so new members are identified if the mother is already a member,
               ! but the daughter not yet.
               !
               IsNewFarFamily =  (IsRelated(DaughterNuclide,MotherNuclide) &
               &    .AND. NuclideFamily%IsFarFamily(MotherNuclide) &
               &    .AND. (.NOT.NuclideFamily%IsFarFamily(DaughterNuclide))&
               &    .AND. (NuclideSpecs(DaughterNuclide)%HalfTime < 1.E80_Float))

               if (IsNewFarFamily) then
                  NuclideFamily%IsFarFamily(DaughterNuclide) = .TRUE.
                  NFarAdded = NFarAdded + 1
                  NuclideFamily%NFarFamily = NuclideFamily%NFarFamily + 1
               endif
               !
               ! Now check if it is part of the narrow family of non-orphans and update bookkeeping if affirmative
               !
               IsNewFamily =  (RegularizedIsRelated(DaughterNuclide,MotherNuclide) &
               &    .AND. NuclideFamily%IsFamily(MotherNuclide) &
               &    .AND. (.NOT.NuclideFamily%IsFamily(DaughterNuclide))&
               &    .AND. (NuclideSpecs(DaughterNuclide)%HalfTime < 1.E80_Float))
               !
               ! Add new family member
               !
               if (IsNewFamily) then
                  NuclideFamily%IsFamily(DaughterNuclide) = .TRUE.
                  NAdded = NAdded + 1
                  NuclideFamily%NFamily = NuclideFamily%NFamily + 1
                  NuclideFamily%FamilyMember(NuclideFamily%NFamily) = DaughterNuclide
                  NuclideFamily%FamilyNumber(DaughterNuclide) = NuclideFamily%NFamily

                  NuclideFamily%NDaughters(MotherNuclide) = NuclideFamily%NDaughters(MotherNuclide) + 1
                  NuclideFamily%Daughter(MotherNuclide,NuclideFamily%NDaughters(MotherNuclide)) = &
                  & DaughterNuclide

                  NuclideFamily%NMothers(DaughterNuclide) = NuclideFamily%NMothers(DaughterNuclide) + 1
                  NuclideFamily%Mother(DaughterNuclide,NuclideFamily%NMothers(DaughterNuclide)) = &
                  & MotherNuclide

                  if (DebugLevel > 2) then
                     write(*,'(A)') 'Adding daughter '//NuclideSpecs(DaughterNuclide)%NuclideName//' for mother '//&
                     & NuclideSpecs(MotherNuclide)%NuclideName
                  endif
               else if (IsNewFarFamily) then
                  if (DebugLevel > 2) then
                     write(*,'(A)') 'Adding orphan daughter '//NuclideSpecs(DaughterNuclide)%NuclideName//' for mother '//&
                     & NuclideSpecs(MotherNuclide)%NuclideName
                  endif
               endif
            enddo ! loop over all other nuclides
         enddo ! loop over nuclides already in the family

         Ready = (NFarAdded == 0)

         if (DebugLevel > 2) then
            if (Ready) then
               write(*,'(A)') 'Found no new members, halting iterations!'
            else
               write(*,'(A,I0,A,I0,A)') 'Found ',NFarAdded,' new members, of which ',NFarAdded-NAdded,' core members!'&
               & //' Proceeding with next iteration!'
            endif
         endif
      enddo ! while loop
      !
      ! Now add short living family members
      !
      if (DebugLevel > 2) then
         write(*,'(A)') 'The regular family has been formed, now going to add possible fast nuclides!'
      endif

      OrphanCounter = NuclideFamily%NFamily

      NAdded = 0

      do SwiftNuclide = 1,NNuclides
         if (NuclideFamily%IsFarFamily(SwiftNuclide).AND..NOT.(NuclideFamily%IsFamily(SwiftNuclide))) then
            if (DebugLevel > 2) then
               write(*,'(A)') 'Nuclide '//trim(NuclideSpecs(SwiftNuclide)%NuclideName)//' seems interesting! Adding it!'
            endif

            NAdded = NAdded + 1

            OrphanCounter = OrphanCounter + 1

            NuclideFamily%FamilyMember(OrphanCounter) = SwiftNuclide
            NuclideFamily%FamilyNumber(SwiftNuclide) = OrphanCounter
            !
            ! Register its daughters
            !
            do DaughterNuclide = 1,NNuclides
               if (IsRelated(DaughterNuclide,SwiftNuclide)) then
               NuclideFamily%NDaughters(SwiftNuclide) = NuclideFamily%NDaughters(SwiftNuclide) + 1
               NuclideFamily%Daughter(SwiftNuclide,NuclideFamily%NDaughters(SwiftNuclide)) = &
               &   DaughterNuclide
               endif
            enddo ! loop over possible mothers
            !
            ! Register its mothers
            !
            do MotherNuclide = 1,NNuclides
               if (IsRelated(SwiftNuclide,MotherNuclide)) then
               NuclideFamily%NMothers(SwiftNuclide) = NuclideFamily%NMothers(SwiftNuclide) + 1
               NuclideFamily%Mother(SwiftNuclide,NuclideFamily%NMothers(SwiftNuclide)) = &
               & MotherNuclide
               endif
            enddo ! loop over possible mothers
         endif ! This is such a swift nuclide
      enddo ! loop over possibly fast decaying nuclides that have to be added to the family, but with high index numbers

      Ready = (NAdded == 0)

      if (DebugLevel > 2) then
         write(*,'(A,I0,A)') 'Found ',NAdded,' new too fast (orphan) members!'
      endif
      !
      ! Show all relatives for this nuclide
      !
      if (DebugLevel == 1) then
         write(*,'(A)') NuclideSpecs(iNuclide)%NuclideName//' chain:'
         do DaughterNuclide = 1,NNuclides
            if(DaughterNuclide /= iNuclide) then
               if(NuclideFamily%IsFamily(DaughterNuclide)) then
                  write(*,'(5X,A)') NuclideSpecs(DaughterNuclide)%NuclideName
               else if (NuclideFamily%IsFarFamily(DaughterNuclide)) then
                  write(*,'(5X,A)') NuclideSpecs(DaughterNuclide)%NuclideName//' <-- orphan'
               endif
            endif ! Daughter not iNuclide
         enddo
      elseif ((DebugLevel > 1).AND. DoPrint) then
         write(*,'(A,I0,A)') 'This nuclide has a family with ',NuclideFamily%NFarFamily-1,' more members:'
         do DaughterNuclide = 1,NNuclides
            if (DebugLevel > 2) then
               write(*,'(4A,2L5)') 'CollectProgeny: show all relatives: ',NuclideSpecs(iNuclide)%NuclideName,&
               & ' --> ? ',NuclideSpecs(DaughterNuclide)%NuclideName,&
               & NuclideFamily%IsFamily(DaughterNuclide),NuclideFamily%IsFarFamily(DaughterNuclide)
            endif

            if(DaughterNuclide /= iNuclide) then
               if(NuclideFamily%IsFamily(DaughterNuclide)) then
                  write(*,'(A,EN20.10,A)') NuclideSpecs(DaughterNuclide)%NuclideName//' with halflife ',&
                  & NuclideSpecs(DaughterNuclide)%HalfTime,' [s]'
               else if (NuclideFamily%IsFarFamily(DaughterNuclide)) then
                  write(*,'(A,EN20.10,A)') NuclideSpecs(DaughterNuclide)%NuclideName//' with halflife ',&
                  & NuclideSpecs(DaughterNuclide)%HalfTime,' [s]  <-- orphan'
               endif
            endif ! Daughter not iNuclide
         enddo
      endif
   end subroutine CollectProgeny



   subroutine ListTransitions()
      !
      ! List all transitions that take place within 1 family, including orphans
      !
      integer :: MotherNuclide,DaughterNuclide,iDaughter,MAtom,DAtom

      write(*,*)
      write(*,'(A)') 'Transitions in family with head of chain '&
      & //NuclideSpecs(NuclideFamily%FamilyMember(1))%NuclideName//':'

      do MotherNuclide = 1,NNuclides
         if (NuclideFamily%IsFarFamily(MotherNuclide)) then
            do iDaughter = 1,NuclideSpecs(MotherNuclide)%NDaughters
               DaughterNuclide = NuclideSpecs(MotherNuclide)%Daughter(iDaughter)
               if (DaughterNuclide /= 0) then
               !
               ! Find type of transition
               !
               MAtom = NuclideSpecs(MotherNuclide)%AtomNumber
               DAtom = NuclideSpecs(DaughterNuclide)%AtomNumber

               write(*,'(A,A,A30,A,A,A,F15.10)') &
               & NuclideSpecs(MotherNuclide  )%NuclideName,&
               & ' --- ',trim(NuclideSpecs(MotherNuclide)%DecayName(iDaughter)),' --> ',&
               & NuclideSpecs(DaughterNuclide)%NuclideName,&
               & ' fraction ',&
               & NuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter)
               else if ((len_trim(NuclideSpecs(MotherNuclide)%DaughterName(iDaughter)) > 0)&
               & .AND.(trim(NuclideSpecs(MotherNuclide)%DaughterName(iDaughter)) /= 'SF')) then
               write(*,'(4A,F15.10,A)') &
               & NuclideSpecs(MotherNuclide)%NuclideName,&
               & ' ------------> ',&
               & NuclideSpecs(MotherNuclide)%DaughterName(iDaughter),&
               & ' fraction ',&
               & NuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter),&
               & ' STABLE'

               endif
            enddo
         endif ! mother is family
      enddo
   end subroutine ListTransitions



   subroutine MakeEvolutionMatrix(RegularizedMotherDaughterMatrix,x,NSecondsDelay,zDelay)
      !================================================================================================================
      !
      ! Construct exponential of mother/daughter matrix of this family
      !
      !================================================================================================================
      real(Float), intent(in) :: NSecondsDelay
      real(Float), dimension(MaxNuclides,MaxNuclides), intent(in) :: RegularizedMotherDaughterMatrix
      real(dp), dimension(:,:), allocatable, intent(out) :: zDelay
      real(dp), dimension(:,:), allocatable, intent(out) :: x

      integer :: MotherNuclide,DaughterNuclide,ErrorCode
      real(dp), dimension(:,:), allocatable :: yDelay

      if (ALLOCATED(x)) deallocate(x)
      allocate(x(NuclideFamily%NFamily,NuclideFamily%NFamily))
      x = 0._Float

      if (ALLOCATED(yDelay)) deallocate(yDelay)
      allocate(yDelay(NuclideFamily%NFamily,NuclideFamily%NFamily))
      yDelay = 0._Float

      if (ALLOCATED(zDelay)) deallocate(zDelay)
      allocate(zDelay(NuclideFamily%NFamily,NuclideFamily%NFamily))
      zDelay = 0._Float
      !
      ! Copy time-rate of change from large matrix for all nuclides to smaller matrix
      !
      do MotherNuclide = 1,NuclideFamily%NFamily
         do DaughterNuclide = 1,NuclideFamily%NFamily
            x(DaughterNuclide,MotherNuclide) = &
            & RegularizedMotherDaughterMatrix(NuclideFamily%FamilyMember(DaughterNuclide),&
            & NuclideFamily%FamilyMember(MotherNuclide))
         enddo
      enddo
      !
      ! Macroscopic time-step over 1 minute, 1 hour and 1 day
      !
      yDelay = NSecondsDelay*x
      !
      ! Exponential
      !
      call dmexp(yDelay,NuclideFamily%NFamily,zDelay,ErrorCode)
   end subroutine MakeEvolutionMatrix



   subroutine AddOrphansBelow(iNuclide,TheDecayMatrix,zDelay,NFirstMotherOrphans)
      !================================================================================================================
      !
      ! Modify a transition matrix to include the orphan transitions below this mother nuclide that have been short-circuited
      ! from the database plus orphans of orphans etcetera. Only the set of decay products for the head of chain is repaired,
      ! i.e. the line in the transition matrix for element iNuclide.
      !
      !================================================================================================================
      integer, intent(in) :: iNuclide
      real(Float), dimension(MaxNuclides,MaxNuclides), intent(inout) :: TheDecayMatrix
      real(dp), dimension(:,:), allocatable, intent(inout) :: zDelay
      integer, intent(out) :: NFirstMotherOrphans

      integer :: iOrphan,NAdded,jOrphan,iIteration,iiNuclide,jjNuclide,MotherNuclide,&
      & NAddedRegular,NAddedOrphans
      logical :: Ready,NotAlreadyAMember,ItsMotherIs_ii
      logical, dimension(MaxNuclides) :: TooFastFullyChecked

      integer, parameter :: DebugLevel = 0

      if (DebugLevel > 0) write(*,*)
      if (DebugLevel > 0) write(*,'(A)') '-------------------- Start of search for orphan chain ------------------------'
      if (DebugLevel > 0) write(*,*)
      if (DebugLevel > 0) write(*,'(A,I4,A,A)') 'Checking orphan chain below nuclide ',iNuclide,&
      & ' : ',trim(NuclideSpecs(iNuclide)%NuclideName)

      NFirstMotherOrphans = 0
      Orphanage%Orphan(:)%Active = .FALSE. ! None of the orphan transitions is yet recognized as relevant
      !
      ! Search for too fast progeny in iterations until no new too fast nuclides are found
      !
      Ready = .FALSE.
      iIteration = 0

      do while (.NOT.Ready)
         iIteration = iIteration + 1
         if (DebugLevel > 0) write(*,*)
         if (DebugLevel > 0) write(*,'(A,I0)') '############# Iteration ',iIteration
         !
         ! If any arrow into a too fast nuclide, either from regular progeny of the head of chain or from another
         ! too fast nuclide in the family, is inactive, switch off this daughter nuclide from TooFastFullyChecked
         !
         call MarkFullyCheckedOrphans()
         !
         ! Now fill in the gaps in the transition matrix for too fast nuclide
         !
         call AddOrphansFromRegularFamily(NAddedRegular)
         call AddOrphansFromCheckedOrphans(NAddedOrphans)

         NAdded = NAddedRegular + NAddedOrphans

         NFirstMotherOrphans = NFirstMotherOrphans + NAdded

         if (DebugLevel > 0) write(*,'(A,I0)') 'Number of orphans added in this iteration: ',NAdded

         Ready = (NAdded == 0)

      enddo ! loop until no new orphans added

      if (DebugLevel > 0) write(*,*)
      if (DebugLevel > 0) write(*,'(A)') '-------------------- End of search for orphan chain --------------------------'
      if (DebugLevel > 0) write(*,*)

      contains

         subroutine MarkFullyCheckedOrphans
            TooFastFullyChecked = .TRUE. ! De-selection of non-fully checked too fast nuclides will be done immediately below

            do iOrphan = 1,Orphanage%NOrphans ! This loop includes all orphans, also the irrelevant ones not in this chain
               !
               ! See if its mother is either regular progeny in the chain at hand (which needs to be addressed)...
               !
               do MotherNuclide = 1,NuclideFamily%NFamily
                  iiNuclide = NuclideFamily%FamilyMember(MotherNuclide)
                  if (Orphanage%Orphan(iOrphan)%Mother == iiNuclide) then ! So now the mother of the orphan is regular family

                     if (.NOT.Orphanage%Orphan(iOrphan)%Active) then ! but the transition to a too fast nuclide in not yet active
                        TooFastFullyChecked(Orphanage%Orphan(iOrphan)%Daughter) = .FALSE.
                        ! Then the too fast daughter is not yet a good basis for searching for new too fast progeny below

                        if (DebugLevel > 0) write(*,'(A,I0,A)') 'MarkFullyCheckedOrphans: In iteration ',iIteration,&
                        & ' switched off orphan '&
                        & //trim(NuclideSpecs(Orphanage%Orphan(iOrphan)%Daughter)%NuclideName)&
                        & //' for lack of progeny transition '&
                        & //trim(NuclideSpecs(Orphanage%Orphan(iOrphan)%Mother)%NuclideName)//' --> '&
                        & //trim(NuclideSpecs(Orphanage%Orphan(iOrphan)%Daughter)%NuclideName)
                     endif

                  endif ! We found the mother of this orphan in the regular progeny of the head of chain
               enddo ! loop over mother nuclides
               !
               ! ... or a fully checked orphan that itself is in the big family
               !
               do jOrphan = 1,Orphanage%NOrphans
                  iiNuclide = Orphanage%Orphan(jOrphan)%Daughter
                  if (NuclideFamily%IsFarFamily(iiNuclide)) then
                     if (Orphanage%Orphan(iOrphan)%Mother == iiNuclide) then ! So now the mother of the orphan is an orphan

                        if (.NOT.Orphanage%Orphan(iOrphan)%Active) then ! but this mother orphan transition is not yet active..
                           TooFastFullyChecked(Orphanage%Orphan(iOrphan)%Daughter) = .FALSE. ! then, de-select this nuclide

                           if (DebugLevel > 0) write(*,'(A,I0,A)') 'MarkFullyCheckedOrphans: In iteration ',&
                           & iIteration,' switched off orphan '&
                           & //trim(NuclideSpecs(Orphanage%Orphan(iOrphan)%Daughter)%NuclideName)&
                           & //' for lack of orphan-orphan transition '&
                           & //trim(NuclideSpecs(Orphanage%Orphan(iOrphan)%Mother)%NuclideName)//' --> '&
                           & //trim(NuclideSpecs(Orphanage%Orphan(iOrphan)%Daughter)%NuclideName)
                        endif

                     endif ! We found the mother of this orphan in the regular progeny of the head of chain
                  endif ! Check if the mother of this orphan belongs to this (big) family
               enddo ! loop over jOrphan
            enddo ! loop over iOrphan
         end subroutine MarkFullyCheckedOrphans

         subroutine AddOrphansFromRegularFamily(MyNAdded)
            integer, intent(out) :: MyNAdded
            !
            ! Step 1: Loop over all regular nuclides in this family to see if they have a yet inactive orphan in their progeny
            !
            MyNAdded = 0

            do MotherNuclide = 1,NuclideFamily%NFamily
               iiNuclide = NuclideFamily%FamilyMember(MotherNuclide)
               !
               ! Check all non-active orphans if they can originate from this mother
               !
               do jOrphan = 1,Orphanage%NOrphans
                  !
                  ! You found a new family member if any new orphan is not already a member and if the orphan at hand is the mother of the new orphan
                  !
                  NotAlreadyAMember = (.NOT.Orphanage%Orphan(jOrphan)%Active)

                  ItsMotherIs_ii = (Orphanage%Orphan(jOrphan)%Mother == iiNuclide)

                  if (NotAlreadyAMember.AND.ItsMotherIs_ii) then
                     Orphanage%Orphan(jOrphan)%Active = .TRUE.
                     !
                     ! Modify transition matrix by adding the yield x activity of the mother of the orphan to the activity of the orphan
                     !
                     jjNuclide = Orphanage%Orphan(jOrphan)%Daughter
                     Orphanage%OrphanActive(jjNuclide) = .TRUE. ! Set this nuclide to active orphan nuclide (possibly already active...)
                     !
                     ! zDelay(MotherNuclide,1) gives the amount of regular family member MotherNuclide
                     ! that results from the head of chain after the given delay.
                     ! Orphanage%Orphan(jOrphan)%yield gives the fraction of this MotherNuclide that is passed on to its too fast daughter
                     !
                     TheDecayMatrix(jjNuclide,iNuclide) = TheDecayMatrix(jjNuclide,iNuclide) &
                     & + Orphanage%Orphan(jOrphan)%yield * zDelay(MotherNuclide,1)

                     MyNAdded = MyNAdded + 1

                     if (DebugLevel > 0) write(*,'(A,F20.15)') 'AddOrphansFromRegularFamily: Found new orphan transition: '&
                     & //NuclideSpecs(jjNuclide)%NuclideName &
                     & //' from mother '&
                     & //NuclideSpecs(iiNuclide)%NuclideName&
                     & //' with yield ',Orphanage%Orphan(jOrphan)%Yield
                  endif
               enddo ! loop over non-active orphans
            enddo ! loop over all regular nuclides in this family
         end subroutine AddOrphansFromRegularFamily

         subroutine AddOrphansFromCheckedOrphans(MyNAdded)
            integer, intent(out) :: MyNAdded
            !
            ! Step 2: Loop over all fully checked orphan nuclides in this big family to see if they have a yet inactive orphan in their progeny
            !
            MyNAdded = 0

            do iOrphan = 1,Orphanage%NOrphans
               iiNuclide = Orphanage%Orphan(iOrphan)%Daughter ! i.e. the endpoint of an orphan arrow, the too fast nuclide itself
               if (NuclideFamily%IsFarFamily(iiNuclide)) then ! the too fast nuclide is far family!
                  if (TooFastFullyChecked(iiNuclide)) then ! and all its incoming arrow have been accounted for, so not it can pass its activity to even lower decks!
                     !
                     ! Check all non-active orphans if they can originate from this (too fast) mother
                     !
                     do jOrphan = 1,Orphanage%NOrphans
                        !
                        ! You found a new family member if any new orphan is not already a member and if the orphan at hand is the mother of the new orphan
                        !
                        NotAlreadyAMember = (.NOT.Orphanage%Orphan(jOrphan)%Active)

                        ItsMotherIs_ii = (Orphanage%Orphan(jOrphan)%Mother == iiNuclide)

                        if (NotAlreadyAMember.AND.ItsMotherIs_ii) then
                           Orphanage%Orphan(jOrphan)%Active = .TRUE.
                           !
                           ! Modify transition matrix by adding the yield x activity of the mother of the orphan to the activity of the orphan.
                           ! As the orphan nuclides do not form part of the family, we take the associated element
                           ! of the full matrix TheDecayMatrix instead of from family matrix zDelay.
                           !
                           jjNuclide = Orphanage%Orphan(jOrphan)%Daughter
                           Orphanage%OrphanActive(jjNuclide) = .TRUE. ! Set this nuclide to active orphan nuclide (possibly already active...)

                           TheDecayMatrix(jjNuclide,iNuclide) = TheDecayMatrix(jjNuclide,iNuclide) &
                           & + Orphanage%Orphan(jOrphan)%yield * TheDecayMatrix(iiNuclide,iNuclide)

                           MyNAdded = MyNAdded + 1

                           if (DebugLevel > 0) write(*,'(A)') 'AddOrphansFromCheckedOrphans: Found new orphan transition: '&
                           & //NuclideSpecs(jjNuclide)%NuclideName &
                           & //' from mother '&
                           & //NuclideSpecs(iiNuclide)%NuclideName
                        endif
                     enddo ! loop over non-active orphans
                  endif ! orphan active
               endif ! possible mother is an orphan herself that belongs to the big family
            enddo ! loop over all active orphans
         end subroutine AddOrphansFromCheckedOrphans
   end subroutine AddOrphansBelow
end module LibENDF
