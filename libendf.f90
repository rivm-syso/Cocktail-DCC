MODULE LibENDF
   ! ____________________________________________________
   !
   ! Version: 11 July 2023
   !
   ! Developed For:
   !
   ! VLH   Centre for Environmental Safety and Security
   ! RIVM - National Institute for Public Health and the Environment
   !
   ! PO Box 1,
   ! NL - 3720 BA Bilthoven
   ! The Netherlands

   USE libxmath
   USE libutil
   USE libexponential

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: NNuclides,NuclideSpecs,MakeMotherDaughterMatrix,&
   & MotherDaughterMatrix,FindOrphans,RegularizeNuclides,&
   & NuclideType,MaxNuclides,Orphanage,NuclideFamily,&
   & CollectProgeny,ListTransitions,MakeEvolutionMatrix,AddOrphansBelow,&
   & GetNuclideNumber,MassNuc2NucMass,EnsureHyphen,ReadNuclideSpecs,AtomName,&
   & IsRelated,RegularizedIsRelated,ReadNProcessENDFNuclideSpecs

   CHARACTER(DefaultLength) :: ENDFPath = './../ENDF-B-VIII.0_decay/'
   CHARACTER(DefaultLength) :: ENDFSpontaneousPath = './../ENDFB_spontaneous/'

   INTEGER, PARAMETER :: MaxAtoms = 118 ! largest atom number available

   CHARACTER(2), DIMENSION(0:MaxAtoms) :: AtomName
   !
   ! Number of nuclide groups
   !
   INTEGER, PARAMETER :: NNuclideGroups = 11
   CHARACTER(5), DIMENSION(0:NNuclideGroups), PARAMETER :: NuclideGroupName &
      & = (/'Other',&
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
      &     'U    '/)

   TYPE AtomType
      CHARACTER(13) :: Name
      CHARACTER(3) :: Symbol
      INTEGER :: Period,Group,NuclideGroup
      CHARACTER(24) :: ChemicalSeries
      REAL(Float) :: Mass
   END TYPE AtomType

   INTEGER :: NAtoms = 0
   TYPE(AtomType), DIMENSION(0:MaxAtoms) :: AtomSpecs
   !
   ! The maximum number of nuclides that we will read from the database:
   !
   INTEGER, PARAMETER :: MaxNuclides = 4000
   !
   !
   INTEGER, PARAMETER :: MaxNDaughters = 15 ! Sufficient for "normal decay", including assignment of very fast decays to their grandmothers
   !
   ! The true number of nuclides found:
   !
   INTEGER :: NNuclides

   TYPE NuclideType
      CHARACTER(10) :: NuclideName
      CHARACTER(10), DIMENSION(MaxNDaughters) :: DaughterName
      CHARACTER(40), DIMENSION(MaxNDaughters) :: DecayName
      INTEGER :: NDecayModes,AtomNumber,NHadrons,NDaughters,MetaStableMode,NuclideGroup
      CHARACTER(15) :: AtomName
      INTEGER, DIMENSION(MaxNDaughters) :: Daughter
      REAL(Float) :: HalfTime,AtomMass
      REAL(Float), DIMENSION(MaxNDaughters) :: DaughterFraction,DecayFraction
      LOGICAL :: IsOrphan
   END TYPE NuclideType
   !
   ! The whole lot of information will be stored in:
   !
   TYPE(NuclideType), DIMENSION(0:MaxNuclides) :: NuclideSpecs
   !
   ! A single datafile can have so many lines
   !
   INTEGER, PARAMETER :: MaxNLines = 30000

   INTEGER, DIMENSION(MaxNLines) :: iRecordType
   LOGICAL, DIMENSION(MaxNLines) :: IsContinuation
   !
   ! The following types of records can be found:
   !
   INTEGER, PARAMETER :: iUnknownRecord =  0
   INTEGER, PARAMETER :: iTextRecord    =  1
   INTEGER, PARAMETER :: iContRecord    =  2
   INTEGER, PARAMETER :: iListRecord    =  3
   INTEGER, PARAMETER :: iTab1Record    =  4
   INTEGER, PARAMETER :: iTab2Record    =  5
   INTEGER, PARAMETER :: iIntgRecord    =  6
   INTEGER, PARAMETER :: iDecayRecord   =  7

   INTEGER, PARAMETER :: NRecordTypes   =  7 ! should be same as last number above

   CHARACTER(7), DIMENSION(0:NRecordTypes), PARAMETER :: RecordName = &
      & (/'Unknown',&!   0
      &   'Text   ',&!   1
      &   'Cont   ',&!   2
      &   'List   ',&!   3
      &   'Tab1   ',&!   4
      &   'Tab2   ',&!   5
      &   'Intg   ',&!   6
      &   'Decay  '/)!   7

   CHARACTER(1), DIMENSION(0:4), PARAMETER :: MetastabilityName = &
      & (/' ',&!   0
      &   'm',&!   1
      &   'n',&!   2
      &   'o',&!   3
      &   'p'/)!   4

   INTEGER, PARAMETER :: iDecayModeGamma              =  0
   INTEGER, PARAMETER :: iDecayModeBeta               =  1
   INTEGER, PARAMETER :: iDecayModeElectronCapture    =  2
   INTEGER, PARAMETER :: iDecayModeIsomericTransition =  3
   INTEGER, PARAMETER :: iDecayModeAlpha              =  4
   INTEGER, PARAMETER :: iDecayModeNeutron            =  5
   INTEGER, PARAMETER :: iDecayModeSpontaneousFission =  6
   INTEGER, PARAMETER :: iDecayModeProton             =  7
   INTEGER, PARAMETER :: iDecayModeUnknown            = 10

   CHARACTER(19), DIMENSION(0:10), PARAMETER :: DecayModeName = &
      & (/'Gamma              ',&!   0
      &   'Beta               ',&!   1
      &   'Electron capture   ',&!   2
      &   'Isomeric transition',&!   3
      &   'Alpha              ',&!   4
      &   'Neutron            ',&!   5
      &   'Spontaneous fission',&!   6
      &   'Proton             ',&!   7
      &   '-invalid1-         ',&!   8
      &   '-invalid2-         ',&!   9
      &   'Unknown            '/)!  10

   INTEGER, DIMENSION(0:10), PARAMETER :: AtomStep = &
      & (/ 0,& ! Gamma
      &    1,& ! Beta
      &   -1,& ! Electron capture
      &    0,& ! Isomeric transition
      &   -2,& ! Alpha
      &    0,& ! Neutron
      & 1000,& ! Spontaneous fission
      &   -1,& ! Proton
      & 1000,& ! -invalid1-
      & 1000,& ! -invalid2-
      & 1000/)  ! Unknown

   INTEGER, DIMENSION(0:10), PARAMETER :: MassStep = &
      & (/ 0,& ! Gamma
      &    0,& ! Beta
      &    0,& ! Electron capture
      &    0,& ! Isomeric transition
      &   -4,& ! Alpha
      &   -1,& ! Neutron
      & 1000,& ! Spontaneous fission
      &   -1,& ! Proton
      & 1000,& ! -invalid1-
      & 1000,& ! -invalid2-
      & 1000/)  ! Unknown

   REAL(Float), DIMENSION(MaxNuclides,MaxNuclides) :: MotherDaughterMatrix,DecayMatrix

   TYPE OrphanType
      INTEGER :: Mother,Daughter
      REAL(Float) :: Yield,CompoundYield
      LOGICAL :: Active
   END TYPE OrphanType

   INTEGER, PARAMETER :: MaxNOrphans = 1700

   TYPE OrphanageType
      REAL(Float) :: tMin
      INTEGER :: NOrphans,NTooUnStable
      LOGICAL, DIMENSION(MaxNuclides) :: IsTooUnstable
      TYPE(OrphanType), DIMENSION(MaxNOrphans) :: Orphan
      LOGICAL, DIMENSION(MaxNuclides) :: OrphanActive
   END TYPE OrphanageType

   TYPE(OrphanageType) :: Orphanage

   LOGICAL, DIMENSION(MaxNuclides,MaxNuclides) :: IsRelated,RegularizedIsRelated
   !
   ! Structure for storing relatives stemming from 1 nuclide
   !
   TYPE NuclideFamilyType
      INTEGER :: NFamily,NFarFamily
      ! Equivalent of compound yield: if you start with any (specified) family member with activity 1, how much of the other
      ! nuclides do you find? Only contributions are counted that are linked to the starting nuclide via fast-decaying nuclides.
      REAL(Float), DIMENSION(MaxNuclides) :: RemainingFraction
      INTEGER, DIMENSION(MaxNuclides) :: FamilyMember,NDaughters,NMothers,FamilyNumber
      INTEGER, DIMENSION(MaxNuclides,MaxNuclides) :: Mother,Daughter
      LOGICAL, DIMENSION(MaxNuclides) :: IsFamily,IsFarFamily,IncomingArrowsChecked,OutgoingArrowsChecked,&
      & DescendentOfStartMember
      LOGICAL, DIMENSION(MaxNuclides,MaxNuclides) :: IsRelated
   END TYPE NuclideFamilyType

   TYPE(NuclideFamilyType) :: NuclideFamily

CONTAINS
   INTEGER FUNCTION GetAtomNumber(Symbol)
      !
      ! Search for the given name of an atom in the list of known atoms and
      ! return its number
      !
      CHARACTER(*), INTENT(IN) :: Symbol
      INTEGER :: DumInt,iAtom,TheIndex
      CHARACTER(10) :: DumSymbol

      DumSymbol = TRIM(Symbol)
      CALL Capitalize(DumSymbol)
      TheIndex = INDEX(DumSymbol,'-')
      IF (TheIndex.GT.0) DumSymbol = DumSymbol(1:(TheIndex-1))
      DumInt = 0
      DO iAtom = 1,NAtoms
         IF (TRIM(DumSymbol).EQ.TRIM(AtomSpecs(iAtom)%Symbol)) DumInt = iAtom
      ENDDO
      GetAtomNumber = DumInt
   END FUNCTION GetAtomNumber



   SUBROUTINE InitNuclide(iNuclide)
      !
      ! Set all info to 0 or anything else that shows clearly that the nulcide is stil empty
      !
      INTEGER, INTENT(IN) :: iNuclide

      NuclideSpecs(iNuclide)%NuclideName = ' '
      NuclideSpecs(iNuclide)%DaughterName = ' '
      NuclideSpecs(iNuclide)%AtomName = ' '
      NuclideSpecs(iNuclide)%AtomNumber = 0
      NuclideSpecs(iNuclide)%NuclideGroup = 0
      NuclideSpecs(iNuclide)%Daughter = 0
      NuclideSpecs(iNuclide)%HalfTime = 0._Float
      NuclideSpecs(iNuclide)%DaughterFraction = 0._Float
   END SUBROUTINE InitNuclide



   SUBROUTINE ReadICRPNuclides(InName)
      !
      ! Read the names, decaytimes and daughters for a collection of nuclides
      !
      CHARACTER(*), INTENT(IN) :: InName

      INTEGER :: iNuclide,jNuclide,DumInt,iDaughter,iAtom,i
      CHARACTER(DefaultLength) :: ALine,DumStr
      CHARACTER(2) :: Units
      CHARACTER(7) :: DumStr7
      CHARACTER(7), DIMENSION(MaxNuclides) :: NuclideName
      CHARACTER(11), DIMENSION(4) :: YieldStr
      REAL(Float) :: C,Dum
      LOGICAL, DIMENSION(MaxNuclides) :: Available
      INTEGER, DIMENSION(MaxNuclides) :: SortedIndex,AtomNumber

      INTEGER, PARAMETER :: DebugLevel = 0
      !
      ! Find number of nuclides and sort them before reading the true specs
      !
      OPEN(ScratchFile,FILE=TRIM(InName),FORM='FORMATTED',STATUS='OLD',&
         & POSITION='REWIND')
      READ(ScratchFile,*) ! Skip header

      IF(DebugLevel.GT.1) WRITE(*,*)
      IF(DebugLevel.GT.1) WRITE(*,'(A)') 'First pass: the nuclides in the order as in the file:'
      iNuclide = 0
      DO
         READ(ScratchFile,'(A7)',END=20) DumStr7

         iNuclide = iNuclide + 1
         NuclideName(iNuclide) = TRIM(DumStr7)
         AtomNumber(iNuclide) = GetAtomNumber(NuclideName(iNuclide))
         IF(DebugLevel.GT.1) WRITE(*,'(I4,1X,A7,5X,I3,1X,A20)') iNuclide,NuclideName(iNuclide),&
            & AtomNumber(iNuclide),AtomSpecs(AtomNumber(iNuclide))%Name
            !
            ! Initialize nuclide
            !
         CALL InitNuclide(iNuclide)

      ENDDO
      20 CONTINUE ! ready counting nuclides

         NNuclides = iNuclide ! The number of nuclides actually found in the file
      IF(DebugLevel.GT.1) WRITE(*,*)
      !
      ! Sort the nuclides by atom number in stead of alphabetically
      !
      IF(DebugLevel.GT.1) WRITE(*,'(A)') 'Second pass: sorting the nuclides:'
      Available = .TRUE.
      SortedIndex = 0
      i = 0
      DO iAtom=1,NAtoms
         DO iNuclide = 1,NNuclides
            IF ((AtomNumber(iNuclide).EQ.iAtom).AND.(Available(iNuclide))) THEN
               i = i + 1
               SortedIndex(iNuclide) = i
               Available(iNuclide) = .FALSE.
               IF(DebugLevel.GT.1) WRITE(*,'(I4,1X,I3,1X,A20)') i,iAtom,TRIM(NuclideName(iNuclide))
            ENDIF
         ENDDO
      ENDDO
      IF(DebugLevel.GT.1) WRITE(*,*)
      !
      ! Read the specs
      !
      IF(DebugLevel.GT.1) WRITE(*,'(A)') 'Third pass: reading the specs:'
      REWIND(ScratchFile)
      READ(ScratchFile,*) ! Skip header

      DO jNuclide = 1,NNuclides
         iNuclide = SortedIndex(jNuclide)
         READ(ScratchFile,10,END=20) ALine
         10 FORMAT(A)
         IF (DebugLevel.GT.1) WRITE(*,70) jNuclide,iNuclide
         70 FORMAT('Line: ',I4,': interpreting specs for sorted nuclide ',I4)
         IF (DebugLevel.GT.2) WRITE(*,60) TRIM(ALine)
         60 FORMAT('Specs from database are: "',A,'"')
         !
         ! Get the nuclide name
         !
         READ(ALine,100) NuclideSpecs(iNuclide)%NuclideName,NuclideSpecs(iNuclide)%HalfTime,Units,&
         & DumStr,DumInt,DumInt,DumInt,DumInt,&
         & (NuclideSpecs(iNuclide)%DaughterName(iDaughter),DumStr,&
         &   YieldStr(iDaughter),iDaughter=1,4),Dum,Dum,Dum,DumInt,DumInt,DumInt,DumInt,DumInt,&
         & NuclideSpecs(iNuclide)%AtomMass
         NuclideSpecs(iNuclide)%NHadrons = NINT(NuclideSpecs(iNuclide)%AtomMass)

         100 FORMAT(A7,F8.0,A2,A8,3I7,I6,4(1X,A7,A6,A11),f7.0,2f8.0,3i4,i5,i4,e11.0,e10.0,e9.0)
         !
         ! Calibrate halftime
         !
         IF (TRIM(Units).EQ.'y') THEN
            C = 3600._Float*24._Float*365.25_Float
            IF (DebugLevel.GT.1) WRITE(*,'(A,EN20.5)') &
            & 'Half time is given in year --> conversion to seconds with factor ',C
         ELSE IF (TRIM(Units).EQ.'d') THEN
            C = 3600._Float*24._Float
            IF (DebugLevel.GT.1) WRITE(*,'(A,EN20.5)') &
            & 'Half time is given in day --> conversion to seconds with factor ',C
         ELSE IF (TRIM(Units).EQ.'h') THEN
            C = 3600._Float
            IF (DebugLevel.GT.1) WRITE(*,'(A,EN20.5)') &
            & 'Half time is given in hour --> conversion to seconds with factor ',C
         ELSE IF (TRIM(Units).EQ.'m') THEN
            C = 60._Float
            IF (DebugLevel.GT.1) WRITE(*,'(A,EN20.5)') &
            & 'Half time is given in minute --> conversion to seconds with factor ',C
         ELSE IF (TRIM(Units).EQ.'s') THEN
            C = 1._Float
            IF (DebugLevel.GT.1) WRITE(*,'(A,EN20.5)') &
            & 'Half time is given in second --> conversion to seconds with factor ',C
         ELSE IF (TRIM(Units).EQ.'ms') THEN
            C = 1.E-3_Float
            IF (DebugLevel.GT.1) WRITE(*,'(A,EN20.5)') &
            & 'Half time is given in millisecond --> conversion to seconds with factor ',C
         ELSE IF (TRIM(Units).EQ.'us') THEN
            C = 1.E-6_Float
            IF (DebugLevel.GT.1) WRITE(*,'(A,EN20.5)') &
            & 'Half time is given in microsecond --> conversion to seconds with factor ',C
         ELSE
            WRITE(*,'(A)') 'Strange unit for halftime: '//Units//' !   --> Exiting!'
            CALL EXIT()
         ENDIF

         NuclideSpecs(iNuclide)%HalfTime = NuclideSpecs(iNuclide)%HalfTime*C

         NuclideSpecs(iNuclide)%AtomNumber = GetAtomNumber(NuclideSpecs(iNuclide)%NuclideName)
         NuclideSpecs(iNuclide)%AtomName = AtomSpecs(NuclideSpecs(iNuclide)%AtomNumber)%Name
         IF (DebugLevel.GT.1) THEN
            IF (NuclideSpecs(iNuclide)%AtomNumber.GT.0) THEN
               WRITE(*,'(A,I3,A,A)') 'Recognized atom ',NuclideSpecs(iNuclide)%AtomNumber,' : ',&
               & TRIM(NuclideSpecs(iNuclide)%AtomName)
            ELSE
               WRITE(*,'(A)') 'This is not a nuclide related to an atom in my database...'
            ENDIF
         ENDIF

         NuclideSpecs(iNuclide)%NuclideGroup = AtomSpecs(NuclideSpecs(iNuclide)%AtomNumber)%NuclideGroup
         !
         ! Extract daughters
         !
         DO iDaughter=1,4
            IF (DebugLevel.GT.2) WRITE(*,'(I1,5X,A)') iDaughter,'"'//TRIM(YieldStr(iDaughter))//'"'
            IF (LEN_TRIM(NuclideSpecs(iNuclide)%DaughterName(iDaughter)).GT.0) &
            & READ(YieldStr(iDaughter),'(E11.0)') NuclideSpecs(iNuclide)%DaughterFraction(iDaughter)
         ENDDO

         IF (DebugLevel.GT.1) THEN
            WRITE(*,'(A)') 'Nuclide: '//TRIM(NuclideSpecs(iNuclide)%NuclideName)
            WRITE(*,'(A,EN20.10,A)') 'Half time: ',NuclideSpecs(iNuclide)%HalfTime,' [s]'
            WRITE(*,'(A,A2)') 'Units of half time: ',Units
            DO iDaughter=1,4
               IF (LEN_TRIM(NuclideSpecs(iNuclide)%DaughterName(iDaughter)).GT.0) THEN
               WRITE(*,'(A,I1,A,A7,A,EN20.5)') 'Daughter(',iDaughter,'): ',&
                  & TRIM(NuclideSpecs(iNuclide)%DaughterName(iDaughter)),&
                  & '     Yield: ',NuclideSpecs(iNuclide)%DaughterFraction(iDaughter)
               ENDIF
            ENDDO
            WRITE(*,*)
         ENDIF

         IF (DebugLevel.GT.0) THEN
            WRITE(*,'(A10,EN20.10)') TRIM(NuclideSpecs(iNuclide)%NuclideName),NuclideSpecs(iNuclide)%HalfTime
         ENDIF

      ENDDO ! loop over ICRP-file

      IF(DebugLevel.GT.1) WRITE(*,*)

      CLOSE(ScratchFile)
   END SUBROUTINE ReadICRPNuclides



   SUBROUTINE MatchDaughters
      !
      ! Match the names of the daughters with the other nuclides.
      ! Exclude stable daughters
      !
      INTEGER :: iNuclide,iDaughter
      INTEGER, PARAMETER :: DebugLevel = 0

      IsRelated = .FALSE.

      DO iNuclide = 1,NNuclides
         !
         ! By default the aughters is set to be irrelevant...
         !
         DO iDaughter = 1,MaxNDaughters
            NuclideSpecs(iNuclide)%Daughter(iDaughter) = &
            & GetNuclideNumber(NuclideSpecs(iNuclide)%DaughterName(iDaughter))
            IF (NuclideSpecs(iNuclide)%Daughter(iDaughter).NE.0) THEN
               IsRelated(NuclideSpecs(iNuclide)%Daughter(iDaughter),iNuclide) = .TRUE.
            ENDIF
         ENDDO
         !
         ! Print results
         !
         IF (DebugLevel.GT.0) THEN
            WRITE(*,40) iNuclide,&
               &    '"'//TRIM(NuclideSpecs(iNuclide)%NuclideName)//'"',&
               &    NuclideSpecs(iNuclide)%AtomNumber,&
               &    '"'//TRIM(NuclideSpecs(iNuclide)%AtomName)//'"',&
               &    NuclideSpecs(iNuclide)%HalfTime,&
               &    ('"'//TRIM(NuclideSpecs(NuclideSpecs(iNuclide)%Daughter(iDaughter))%NuclideName)//'"',&
               &    NuclideSpecs(iNuclide)%DaughterFraction(iDaughter),iDaughter=1,4)
            40    FORMAT(I4,'=',A12,'; atom(',I3,')=',A17,' T1/2=',G15.10,&
               & '  --> ',4(A12,' frac=',F10.8,'   '))

            WRITE(*,*) 'Ready reading nuclide!'
            WRITE(*,*) '-----------------------------------------------'
            WRITE(*,*)
            WRITE(*,*)
         ENDIF
      ENDDO ! loop over nuclides
   END SUBROUTINE MatchDaughters



   SUBROUTINE ReadNuclideSpecs(InName)
      !
      ! Read the names, decaytimes and daughters for a collection of nuclides
      !
      CHARACTER(*), INTENT(IN) :: InName
      !
      IF (NAtoms.EQ.0) CALL ReadAtoms('atoms.dat')
      !
      ! Read the datafile
      !
      IF (.NOT.FileExists(TRIM(InName))) THEN
         WRITE(*,'(A)') 'Cannot find ICRP nuclide file '//'/'//TRIM(InName)
         WRITE(*,'(A)')  'You can download this file via the Supplemental Material associated with '&
            & //'http://www.icrp.org/publication.asp?id=ICRP%20Publication%20107'
         WRITE(*,'(A)')  'Exiting!'
         CALL EXIT()
      ENDIF

      WRITE(*,'(A)') 'Going to read nuclide specs from ICRP nuclide file '//TRIM(InName)

      CALL ReadICRPNuclides(InName)

      WRITE(*,'(A)') 'Ready reading nuclide specs, going to match daughters...'
      !
      ! Match the names of the daughters with the other nuclides.
      ! Exclude stable daughers
      !
      CALL MatchDaughters()
      WRITE(*,'(A)') 'Ready matching daughters!'
   END SUBROUTINE ReadNuclideSpecs



   SUBROUTINE MassNuc2NucMass(MyName)
      !
      ! When applicable, convert the nuclide name from format 137mCs to Cs137m
      !
      CHARACTER(*), INTENT(INOUT) :: MyName
      LOGICAL :: HasWrongFormat,IsMetastable,Ready
      CHARACTER(10) :: MassName,NukeName
      INTEGER :: iCharacter

      HasWrongFormat = CharacterIsADigit(MyName(1:1))

      IF (HasWrongFormat) THEN
      !
      ! Get the mass number
      !
      MassName = ' '
      iCharacter = 1
      Ready = .FALSE.
      DO WHILE (.NOT.Ready)
         MassName = TRIM(MassName)//MyName(iCharacter:iCharacter)
         iCharacter = iCharacter + 1
         Ready = (.NOT.CharacterIsADigit(MyName(iCharacter:iCharacter)))
      ENDDO ! not ready
      !
      ! Get the optional 'm' for meta-stable nuclides (case sensitive!!)
      !
      IsMetastable = (MyName(iCharacter:iCharacter).EQ.'m') ! This is case sensitive! Make sure you haven't converted to full uppercase yet!!

      IF (IsMetastable) iCharacter = iCharacter + 1
      !
      ! Get the nuclide name
      !
      NukeName = ' '
      Ready = .FALSE.
      DO WHILE (.NOT.Ready)
         NukeName = TRIM(NukeName)//MyName(iCharacter:iCharacter)
         iCharacter = iCharacter + 1
         Ready = (MyName(iCharacter:iCharacter).EQ.' ')
      ENDDO ! not ready
      !
      ! Synthesis of elements in order. E.g., Cs137m
      !
      IF (IsMetastable) THEN
         MyName = TRIM(NukeName)//TRIM(MassName)//'m'
      ELSE
         MyName = TRIM(NukeName)//TRIM(MassName)
      ENDIF

      ENDIF
   END SUBROUTINE MassNuc2NucMass



   SUBROUTINE EnsureHyphen(MyName)
      !
      ! When applicable, add a hyphen as separator between atom name and mass number
      !
      CHARACTER(*), INTENT(INOUT) :: MyName
      LOGICAL :: HasWrongFormat,IsMetastable,Ready,HasHyphen
      CHARACTER(10) :: MassName,NukeName
      INTEGER :: iCharacter

      HasWrongFormat = CharacterIsADigit(MyName(1:1))

      IF (HasWrongFormat) CALL MassNuc2NucMass(MyName)

      HasHyphen = (INDEX(MyName,'-').GT.0)

      IF (.NOT.HasHyphen) THEN
         !
         ! Get the mass number
         !
         NukeName = ' '
         iCharacter = 1
         Ready = .FALSE.
         DO WHILE (.NOT.Ready)
            NukeName = TRIM(NukeName)//MyName(iCharacter:iCharacter)
            iCharacter = iCharacter + 1
            Ready = CharacterIsADigit(MyName(iCharacter:iCharacter))
         ENDDO ! not ready

         NukeName = TRIM(NukeName)//'-'//MyName((iCharacter):LEN_TRIM(MyName))
         MyName = TRIM(NukeName)
      ENDIF
   END SUBROUTINE EnsureHyphen



   SUBROUTINE ReadAtoms(InName)
      !
      ! Read the specs of all known atoms
      !
      CHARACTER(*), INTENT(IN) :: InName

      INTEGER :: iAtom,DumInt,iNuclideGroup
      CHARACTER(DefaultLength) :: ALine,DumStr
      LOGICAL :: Ready

      INTEGER, PARAMETER :: DebugLevel = 0

      WRITE(*,'(A)') 'Going to read atom specs...'
      OPEN(ScratchFile,FILE=(' ./data/')//TRIM(InName),FORM='FORMATTED',STATUS='OLD',&
         & POSITION='REWIND')
      READ(ScratchFile,*) ! Skip header

      Ready = .FALSE.
      iAtom = 0
      DO WHILE (.NOT.Ready)
         READ(ScratchFile,'(A)',END=20) ALine
         iAtom = iAtom + 1
         IF (DebugLevel.GT.1) WRITE(*,70) iAtom
         70 FORMAT('Interpreting specs for Atom ',I3)
         IF (DebugLevel.GT.1) WRITE(*,60) TRIM(ALine)
         60 FORMAT('Specs from database are: "',A,'"')
         !
         ! Get the atom specs
         !
         READ(ALine( 1: 3),'(I3)')    DumInt
         READ(ALine( 7:19),'(A13)')   AtomSpecs(iAtom)%Name
         READ(ALine(21:23),'(A3)')    AtomSpecs(iAtom)%Symbol
         READ(ALine(27:27),'(I1)')    AtomSpecs(iAtom)%Period
         READ(ALine(34:35),'(I2)')    AtomSpecs(iAtom)%Group
         READ(ALine(41:64),'(A24)')   AtomSpecs(iAtom)%ChemicalSeries
         READ(ALine(65:76),'(F12.0)') AtomSpecs(iAtom)%Mass
         DumStr = ALine(78:LEN(ALine))

         AtomSpecs(iAtom)%NuclideGroup = 0

         IF (LEN_TRIM(DumStr).NE.0) THEN
            DO iNuclideGroup = 0,NNuclideGroups
            IF (TRIM(DumStr).EQ.TRIM(NuclideGroupName(iNuclideGroup))) &
               & AtomSpecs(iAtom)%NuclideGroup = iNuclideGroup
            ENDDO
         ENDIF

         IF (DumInt.NE.iAtom) THEN
            WRITE(*,'(A)') 'Atom number does not match line number! Exiting!!'
            CALL EXIT()
         ENDIF

         IF (DebugLevel.GT.0) THEN
            WRITE(*,200) AtomSpecs(iAtom),NuclideGroupName(AtomSpecs(iAtom)%NuclideGroup)
            200  FORMAT('Specs: "',A13,'"',1X,'"',A3,'"',1X,I2,1X,I2,1X,I2,1X,'"',A24,'"',1X,F12.8,1X,A)
         ENDIF
      ENDDO
      20 CONTINUE ! ready reading data

      NAtoms = iAtom

      WRITE(*,'(A,I3)') '... ready reading specs for this many atoms: ',NAtoms
   END SUBROUTINE ReadAtoms



   SUBROUTINE InitLibENDF()
      !
      ! Prepare some settings before use of this library
      !
      CHARACTER(DefaultLength) :: FName,ALine
      LOGICAL :: IsMetaStable,ValidLine,IsException
      INTEGER :: iLine,MyAtomNumber

      INTEGER, PARAMETER :: DebugLevel = 0

      FName = TRIM(ENDFPath)//'decay.list'

      OPEN(ScratchFile,FILE=FName,FORM='FORMATTED',ACTION='READ')
      !
      ! First, scan the file to get atom names only
      !
      READ(ScratchFile,*)
      READ(ScratchFile,*)
      READ(ScratchFile,*)

      DO
         READ(ScratchFile,'(A)',END=10) ALine
         ValidLine = (ALine(1:1).NE.'-')
         IF (ValidLine) THEN
            READ(ALine(7:9),*) MyAtomNumber
            AtomName(MyAtomNumber) = ALine(11:12)
         ENDIF ! valid line
      ENDDO

      10 CONTINUE
      !
      ! Now, scan again and get the nuclides
      !
      REWIND(ScratchFile)

      READ(ScratchFile,*)
      READ(ScratchFile,*)
      READ(ScratchFile,*)

      NNuclides = 0
      iLine = 3

      DO
         READ(ScratchFile,'(A)',END=20) ALine

         ValidLine = (ALine(1:1).NE.'-')

         IF (ValidLine) THEN
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
            IF (DebugLevel.GT.1) THEN
               WRITE(*,'(A,I4,A)') 'Line ',NNuclides,' = "'//TRIM(ALine)//'"'
            ENDIF

            READ(ALine(7:9),*) NuclideSpecs(NNuclides)%AtomNumber
            NuclideSpecs(NNuclides)%AtomName = ALine(11:12)

            READ(ALine(14:16),*) NuclideSpecs(NNuclides)%NHadrons

            IF (ALine(17:17).EQ.' ') THEN
               NuclideSpecs(NNuclides)%MetaStableMode = 0
            ELSE IF (ALine(17:17).EQ.'M') THEN
               NuclideSpecs(NNuclides)%MetaStableMode = 1
            ELSE IF (ALine(17:17).EQ.'N') THEN
               NuclideSpecs(NNuclides)%MetaStableMode = 2
            ELSE IF (ALine(17:17).EQ.'O') THEN
               NuclideSpecs(NNuclides)%MetaStableMode = 3
            ELSE
               WRITE(*,'(A,I0,A,I0,A,I0,A)') 'Unknown metastability mode "',NuclideSpecs(NNuclides)%MetaStableMode,&
               & '" on line ',iLine,' for nuclide ',NNuclides,'! Exiting!!'
               CALL EXIT()
            ENDIF
            !
            ! ICRP-107 knows Pm-137m and ENDF knows Pm-137. These are probably the same nuclide.
            ! To allow the use of both datasets simultaneously, the ENDF nuclide is mapped on the ICRP-nuclide:
            ! The nuclide name gets an "m". All subsequent references in ENDF to Pm-137 are trapped and re-cast to Pm-137m.
            !
            IsException = ((NuclideSpecs(NNuclides)%AtomName.EQ.'Pm').AND.(NuclideSpecs(NNuclides)%NHadrons.EQ.137))
            IF (IsException) NuclideSpecs(NNuclides)%MetaStableMode = 1

            NuclideSpecs(NNuclides)%NDaughters = 0

            WRITE(NuclideSpecs(NNuclides)%NuclideName,'(A,A,I0)') &
               & TRIM(NuclideSpecs(NNuclides)%AtomName),'-',NuclideSpecs(NNuclides)%NHadrons

            NuclideSpecs(NNuclides)%NuclideName = TRIM(NuclideSpecs(NNuclides)%NuclideName)&
            & //MetastabilityName(NuclideSpecs(NNuclides)%MetaStableMode)

            IF (DebugLevel.GT.0) WRITE(*,'(I4,1X,A)') NNuclides,TRIM(NuclideSpecs(NNuclides)%NuclideName)
         ENDIF ! valid line
      ENDDO

      20 CONTINUE
      CLOSE(ScratchFile)
   END SUBROUTINE InitLibENDF



   INTEGER FUNCTION GetNuclideNumber(NuclideName)
      !
      ! Search for the given name of a nuclide in the list of known nuclides and
      ! return its number
      !
      CHARACTER(*), INTENT(IN) :: NuclideName
      INTEGER :: DumInt,iNuclide
      CHARACTER(10) :: DumNuclideName

      DumInt = 0

      IF (LEN_TRIM(NuclideName).GT.0) THEN
         DumNuclideName = TRIM(NuclideName)
         CALL Capitalize(DumNuclideName)
         DO iNuclide = 1,NNuclides
            IF (TRIM(DumNuclideName).EQ.TRIM(NuclideSpecs(iNuclide)%NuclideName)) DumInt = iNuclide
         ENDDO
      ENDIF ! length of input string <> 0

      GetNuclideNumber = DumInt
   END FUNCTION GetNuclideNumber



   INTEGER FUNCTION AtomNMass2Nuclide(AtomNumber,NHadrons,MetastableMode)
      !
      ! Find the nuclide with these properties
      !
      INTEGER, INTENT(IN) :: AtomNumber,NHadrons,MetastableMode
      INTEGER :: DumInt,iNuclide

      DumInt = 0

      DO iNuclide = 1,NNuclides
      IF (      (AtomNumber    .EQ.NuclideSpecs(iNuclide)%AtomNumber)&
         & .AND.(NHadrons      .EQ.NuclideSpecs(iNuclide)%NHadrons)&
         & .AND.(MetastableMode.EQ.NuclideSpecs(iNuclide)%MetastableMode)) DumInt = iNuclide
      ENDDO

      AtomNMass2Nuclide = DumInt
   END FUNCTION AtomNMass2Nuclide



   INTEGER FUNCTION Zap2Nuclide(Zap,iMetaStableMode)
      !
      ! Convert a ZAP-number to a nuclide index
      !
      INTEGER, INTENT(IN) :: Zap,iMetaStableMode

      INTEGER :: iAtom,iMass,iMother

      iAtom = Zap/1000
      iMass = Zap - 1000*iAtom

      iMother = AtomNMass2Nuclide(iAtom,iMass,iMetastableMode)
      Zap2Nuclide = iMother
   END FUNCTION Zap2Nuclide



   SUBROUTINE CheckRecordType(iNuclide,iLine,ALine)
      !
      ! Assess which type of record this is
      !
      INTEGER, INTENT(IN) :: iNuclide,iLine
      CHARACTER(*), INTENT(IN) :: ALine

      INTEGER :: MatNumber,MFNumber,MTNumber

      INTEGER, PARAMETER :: DebugLevel = 0
      !
      ! Read the information at the end of the line to see what type of record this is
      !
      READ(ALine(67:75),'(I4,I2,I3)') MatNumber,MFNumber,MTNumber

      IF ((MFNumber.EQ.1).AND.(MTNumber.EQ.451)) THEN
         iRecordType(iLine) = iTextRecord
      ELSE IF ((MFNumber.EQ.8).AND.(MTNumber.EQ.457)) THEN
         iRecordType(iLine) = iDecayRecord
      ENDIF

      IF (DebugLevel.GT.0) THEN
         WRITE(*,'(A,I4,A29,3(1X,I4,1X,I2,1X,I3,4X,A))') 'Line ',iLine,&
         & ' has a record of type '//RecordName(iRecordType(iLine)),MatNumber,MFNumber,MTNumber,'"'//TRIM(ALine)//'"'
      ENDIF
   END SUBROUTINE CheckRecordType


   SUBROUTINE Line2Values(ALine,MySubString,MyFloat,MyInteger,SubStringIsInteger,DoingBullocks)
      !
      ! Read 6 values from a line and check if they are real or integer
      !
      CHARACTER(*), INTENT(IN) :: ALine
      INTEGER, DIMENSION(6), INTENT(OUT) :: MyInteger
      REAL(Float), DIMENSION(6), INTENT(OUT) :: MyFloat
      LOGICAL, DIMENSION(6), INTENT(OUT) :: SubStringIsInteger
      LOGICAL, INTENT(OUT) :: DoingBullocks
      CHARACTER(11), DIMENSION(6), INTENT(INOUT) :: MySubString

      INTEGER :: iSubString,TheExponent,TheIndex
      CHARACTER(DefaultLength) :: DumString
      REAL(Float) :: RelativeDifference

      INTEGER, PARAMETER :: DebugLevel = 0

      READ(ALine,'(6A11)') MySubString

      MyFloat = 0._Float
      MyInteger = 0
      DoingBullocks = .FALSE.

      DO iSubString = 1,6
         IF (DebugLevel.GT.2) WRITE(*,'(A,I0,A)') 'Substring ',iSubString,' = "'//MySubString(iSubString)//'"'
      ENDDO
      DO iSubString = 1,6
         IF (DebugLevel.GT.1) WRITE(*,'(A,I0,A)') 'Substring ',iSubString,' = "'//MySubString(iSubString)//'"'

         TheIndex = MAX(INDEX(MySubString(iSubString),'.'),&
         & INDEX(MySubString(iSubString),'+'),INDEX(MySubString(iSubString),'-'))
         IF (TheIndex.GT.0) THEN
            SubStringIsInteger(iSubString) = .FALSE.
            !
            ! Values can be with or without E in the exponential notation, and there can be spaces in the string:
            !
            IF ((INDEX(MySubString(iSubString),'E').NE.0)&
            & .OR.((INDEX(MySubString(iSubString),'+').EQ.0).AND.(INDEX(MySubString(iSubString),'-').EQ.0)))THEN
               READ(MySubString(iSubString),*) MyFloat(iSubString)
            ELSE
               READ(MySubString(iSubString)(1:(TheIndex-1)),*) MyFloat(iSubString)
               DumString = MySubString(iSubString)(TheIndex:11)
               CALL RemoveCharacter(DumString,' ')
               READ(DumString,*) TheExponent
               MyFloat(iSubString) = MyFloat(iSubString) * 10._Float**TheExponent
            ENDIF
            !
            ! If this is an integer in disguise, make it an integer (a real one! ;-) )
            !
            RelativeDifference = 0._Float
            IF (ABS(MyFloat(iSubString)).GT.0._Float) RelativeDifference = &
               & ABS(MyFloat(iSubString)-NINT(MyFloat(iSubString)))/ABS(MyFloat(iSubString))

            IF ((MyFloat(iSubString).LE.1.E6_Float)&
               & .AND.(RelativeDifference.LT.1.E-10_Float)) THEN
               SubStringIsInteger(iSubString) = .TRUE.
               MyInteger(iSubString) = NINT(MyFloat(iSubString))
            ENDIF
         ELSE
            SubStringIsInteger(iSubString) = .TRUE.
            DoingBullocks = .TRUE.
            READ(MySubString(iSubString),*,END=20) MyInteger(iSubString)
            DoingBullocks = .FALSE.
         ENDIF

         20  CONTINUE
         !
         ! Show substring
         !
         IF (DebugLevel.GT.1) THEN
            IF (SubStringIsInteger(iSubString)) THEN
               WRITE(*,'(A,I0,A,I0)') 'Substring ',iSubString,' is an integer: ',MyInteger(iSubString)
            ELSE
               WRITE(*,'(A,I0,A,F15.5)') 'Substring ',iSubString,' is a real    : ',MyFloat(iSubString)
            ENDIF
         ENDIF ! debug
      ENDDO ! loop over substrings
   END SUBROUTINE Line2Values



   SUBROUTINE Parse1DataLine(iLine,iNuclide,ALine,iLineInChapter,ReadyReadingFile,ThisLineHasDecay,iDecay)
      !
      ! Try to find useful info about a known nuclide
      !
      INTEGER, INTENT(IN) :: iLine,iNuclide,iLineInChapter
      CHARACTER(*), INTENT(IN) :: ALine
      LOGICAL, INTENT(INOUT) :: ReadyReadingFile,ThisLineHasDecay
      INTEGER, INTENT(INOUT) :: iDecay

      CHARACTER(DefaultLength) :: MyDecayModeLine
      CHARACTER(80) :: PartOfLine
      INTEGER :: TheIndex,iSubString,TheExponent,DaughterAtom,DaughterMass,MyDecayMode,MyMetastableMode,&
      & MyAtomStep,MyMassStep,iParticle,iCompoundDecay,DaughterNuclide,iDaughter
      LOGICAL, DIMENSION(6) :: SubStringIsInteger
      INTEGER, DIMENSION(6) :: MyInteger
      REAL(Float), DIMENSION(6) :: MyFloat
      REAL(Float) :: MyFraction
      LOGICAL :: Ready,IsValidDecay,DoingBullocks,IsException
      CHARACTER(11), DIMENSION(6) :: MySubString

      INTEGER, PARAMETER :: DebugLevel = 0

      DoingBullocks = .FALSE.
      ReadyReadingFile = .FALSE.
      ThisLineHasDecay = .FALSE.

      IF (DebugLevel.GT.1) WRITE(*,'(A,I0,A)') 'Line ',iLine,' : "'//TRIM(ALine)//'"'
      !
      ! Only interpret decay info
      !
      IF (iRecordType(iLine).EQ.iDecayRecord) THEN
         IF (DebugLevel.GT.1) WRITE(*,'(A,I4,A)') 'Decay information record ',iLineInChapter,' = "'//TRIM(ALine)//'"'
         !
         ! Read the 6 values and check if they are real/integer
         !
         CALL Line2Values(ALine,MySubstring,MyFloat,MyInteger,SubStringIsInteger,DoingBullocks)
         ! Check this criterion!!!!!

         ReadyReadingFile =      SubStringIsInteger(1).AND.(MyInteger(1).EQ.0)&
            &                  .AND.SubStringIsInteger(2).AND.(MyInteger(2).EQ.0)&
            &                  .AND.SubStringIsInteger(3)&
            &                  .AND.SubStringIsInteger(4).AND.(MyInteger(4).EQ.0)&
            &                  .AND.SubStringIsInteger(5)&
            &                  .AND.SubStringIsInteger(6)&
            &                  .AND.((iLineInChapter.EQ.2).OR.(iLineInChapter.GE.5))

         ReadyReadingFile = ReadyReadingFile.OR.&
         & ((iLineInChapter.GE.5).AND.(iDecay.GE.NuclideSpecs(iNuclide)%NDecayModes))

         IF (DoingBullocks) THEN
            IF (DebugLevel.GT.0) THEN
               WRITE(*,'(A,I0)') 'Encountering compliance problems in file for '&
               & //NuclideSpecs(iNuclide)%NuclideName//' at line ',iLine
            ENDIF ! debug
         ELSE IF (ReadyReadingFile) THEN
            IF (DebugLevel.GT.1) THEN
               WRITE(*,'(A,I0)') 'Ready reading file for '&
               & //NuclideSpecs(iNuclide)%NuclideName//' at line ',iLine
            ENDIF ! debug
         ELSE
            IF (iLineInChapter.EQ.2) THEN
               NuclideSpecs(iNuclide)%HalfTime = MyFloat(1) ! Even in cases where you have value e.g. 1.000, this goes well

               IF (DebugLevel.GT.1) WRITE(*,'(A,EN15.5,A)') 'Halflife = ',NuclideSpecs(iNuclide)%HalfTime,' s'
            ELSE IF (iLineInChapter.EQ.4) THEN
               NuclideSpecs(iNuclide)%NDecayModes = MyInteger(6)

               IF (DebugLevel.GT.1) WRITE(*,'(A,I0)') 'Number of decay modes: ',NuclideSpecs(iNuclide)%NDecayModes
            ELSE IF (iLineInChapter.GE.5) THEN

               MyMetastableMode = MyInteger(2)
               MyFraction = MyFloat(5)
               !
               ! If you have a non-trivial decay fraction, update nuclide database
               !
               IF ((MyFraction.GT.0._Float).AND.(MyFraction.LE.1._Float)) THEN
                  IF (DebugLevel.GT.1) WRITE(*,'(A,EN15.5)') 'Acceptable fraction = ',MyFraction

                  iParticle = 1
                  MyAtomStep = 0
                  MyMassStep = 0
                  IsValidDecay = .FALSE.
                  MyDecayModeLine = ' '
                  iCompoundDecay = NINT(1000000._Float*MyFloat(1))

                  Ready = (iCompoundDecay.LT.1000000)

                  DO WHILE (.NOT.Ready)
                     MyDecayMode = iCompoundDecay/10**(7-iParticle)

                     iCompoundDecay = iCompoundDecay - MyDecayMode*10**(7-iParticle)

                     Ready = (iCompoundDecay.EQ.0)

                     IF (iCompoundDecay.LT.0) THEN
                        WRITE(*,'(A,I0)') 'Faulty value for iCompoundDecay: ',iCompoundDecay
                        CALL EXIT()
                     ENDIF ! rubbish

                     IF (DebugLevel.GT.1) WRITE(*,'(A,A)') 'Found decay mode ',DecayModeName(MyDecayMode)

                     MyDecayModeLine = TRIM(MyDecayModeLine)//' + '//DecayModeName(MyDecayMode)

                     IF (     (MyDecayMode.GE.0)&
                        &   .AND.(MyDecayMode.LE.7)) THEN
                        IsValidDecay = .TRUE.
                        MyAtomStep = MyAtomStep + AtomStep(MyDecayMode)
                        MyMassStep = MyMassStep + MassStep(MyDecayMode)
                     ELSE
                        IsValidDecay = .FALSE.
                        IF (DebugLevel.GT.1) WRITE(*,'(A,A,A)') &
                        & 'Found disqualifying decay mode ',DecayModeName(MyDecayMode),', making this decay invalid!!'
                     ENDIF ! valid decay mode

                     iParticle = iParticle + 1

                  ENDDO ! ready checking all particles in this single decay

                  MyDecayModeLine = MyDecayModeLine(4:DefaultLength)
                  !
                  ! Account for the full train of released particles in this decay, e.g. NN, in 1 go
                  !
                  IF (IsValidDecay) THEN
                     ThisLineHasDecay = .TRUE.
                     iDecay = iDecay + 1
                     NuclideSpecs(iNuclide)%DecayFraction(iDecay) = MyFraction

                     IF (.NOT.(MyDecayMode.EQ.iDecayModeSpontaneousFission)) THEN
                        NuclideSpecs(iNuclide)%NDaughters = NuclideSpecs(iNuclide)%NDaughters + 1
                        iDaughter = NuclideSpecs(iNuclide)%NDaughters

                        DaughterAtom = NuclideSpecs(iNuclide)%AtomNumber + MyAtomStep
                        DaughterMass = NuclideSpecs(iNuclide)%NHadrons   + MyMassStep

                        WRITE(NuclideSpecs(iNuclide)%DaughterName(iDaughter),'(A,A,I0)') &
                           &  TRIM(AtomName(DaughterAtom)),'-',DaughterMass
                        !
                        ! Trap Pm-137 and recast it to Pm-137m, which is the only isotope in ICRP-107 with this combination.
                        ! To allow for the simultaneous use of ENDF with ICRP-107, Pm-137 is renamed Pm-137m
                        !
                        IsException = ((DaughterAtom.EQ.61).AND.(DaughterMass.EQ.137))
                        IF (IsException) MyMetastableMode = 1

                        NuclideSpecs(iNuclide)%DaughterName(iDaughter) = &
                           & TRIM(NuclideSpecs(iNuclide)%DaughterName(iDaughter))&
                           & //MetaStabilityName(MyMetastableMode)

                        NuclideSpecs(iNuclide)%Daughter(iDaughter) = &
                           & GetNuclideNumber(NuclideSpecs(iNuclide)%DaughterName(iDaughter))

                        NuclideSpecs(iNuclide)%DaughterFraction(iDaughter) = MyFraction
                        NuclideSpecs(iNuclide)%DecayName(iDaughter) = TRIM(MyDecayModeLine)

                        IF (DebugLevel.GT.0) WRITE(*,'(A,I4,A,A7,A,EN15.5,A,A,A7,A,F20.10,A)') &
                           & 'Nuclide ',iNuclide,': ',&
                           & NuclideSpecs(iNuclide)%NuclideName,&
                           & ' with halflife ',NuclideSpecs(iNuclide)%HalfTime,' s',&
                           & ' to daughter ',NuclideSpecs(iNuclide)%DaughterName(iDaughter),&
                           & ' takes a fraction ',MyFraction,&
                           & ', decay: '//TRIM(NuclideSpecs(iNuclide)%DecayName(iDaughter))

                     ELSE
                        IF (DebugLevel.GT.0) THEN
                           WRITE(*,'(A,I4,A,A7,A,EN15.5,A,A,F20.10,A)') &
                              & 'Nuclide ',iNuclide,': ',&
                              & NuclideSpecs(iNuclide)%NuclideName,&
                              & ' with halflife ',NuclideSpecs(iNuclide)%HalfTime,&
                              & ' s to many daughters  ',&
                              & ' takes a fraction ',MyFraction,&
                              & ', decay: Spontaneous fission'
                        ENDIF
                     ENDIF ! SF or other valid decay
                  ENDIF ! valid decay mode
               ELSE
                  IF (DebugLevel.GT.1) WRITE(*,'(A,EN15.5)') 'Unacceptable fraction = ',MyFraction
               ENDIF ! yield > 0
            ENDIF ! line 2 or 5 and later in the decay chapter
         ENDIF ! not doing bullocks...
      ENDIF ! Inside chapter of decay records
   END SUBROUTINE Parse1DataLine



   SUBROUTINE ParseENDFFile(iNuclide)
      !
      ! Parse 1 of the 3821 ENDF datafiles
      !
      INTEGER, INTENT(IN) :: iNuclide

      CHARACTER(DefaultLength) :: FName,ALine
      CHARACTER(20) :: MyTail
      LOGICAL :: Ready,IsNewChapter,FoundMeaningfulDecay,ThisLineHasDecay,IsException
      INTEGER :: iLine,iLineInChapter,iDecay

      INTEGER, PARAMETER :: DebugLevel = 0
      !
      ! Construct filename
      !
      FName = 'dec-???_'
      WRITE(FName(5:7),'(I3.3)') NuclideSpecs(iNuclide)%AtomNumber
      FName = TRIM(FName)//TRIM(NuclideSpecs(iNuclide)%AtomName)
      MyTail = '_???'
      WRITE(MyTail(2:4),'(I3.3)') NuclideSpecs(iNuclide)%NHadrons
      !
      ! Trap Pm-137 and recast it to Pm-137m, which is the only isotope in ICRP-107 with this combination.
      ! To allow for the simultaneous use of ENDF with ICRP-107, Pm-137 is renamed Pm-137m
      !
      IsException = ((NuclideSpecs(iNuclide)%AtomNumber.EQ.61).AND.(NuclideSpecs(iNuclide)%NHadrons.EQ.137))



      IF ((NuclideSpecs(iNuclide)%MetaStableMode.NE.0).AND..NOT.IsException) THEN
         WRITE(MyTail(5:6),'(A1,I1)') 'm',NuclideSpecs(iNuclide)%MetaStableMode
      ENDIF
      FName = TRIM(FName)//TRIM(MyTail)

      FName = TRIM(ENDFPath)//TRIM(FName)//'.endf'

      IF (IsException) WRITE(*,'(A)') 'To match ENDF with ICRP-107, instead of '&
      & //TRIM(NuclideSpecs(iNuclide)%NuclideName)&
      & //' I will read data from '//TRIM(FName)
      !   IF (DebugLevel.GT.0) THEN
      !     WRITE(*,'(A,I0,A)') 'Parsing file '//TRIM(FName)//' for nuclide ',iNuclide,'...'
      !   ENDIF ! Debug
      IF (.NOT.FileExists(FName)) THEN
         WRITE(*,'(A,I0,A)') 'Cannot find file "'//TRIM(FName)//'" for nuclide ',&
         & iNuclide,': '//TRIM(NuclideSpecs(iNuclide)%NuclideName)//'! Exiting!!'
         CALL EXIT()
      ENDIF
      !
      ! Read contents
      !
      OPEN(ScratchFile,FILE=TRIM(FName),FORM='FORMATTED',ACTION='READ')

      Ready = .FALSE.
      iLine = 0
      IsContinuation = .FALSE.
      iRecordType = iUnknownRecord
      iLineInChapter = 0
      IsNewChapter = .TRUE.
      FoundMeaningfulDecay = .FALSE.
      iDecay = 0

      DO WHILE (.NOT.Ready)
         READ(ScratchFile,'(A)',END=10) ALine
         iLine = iLine + 1
         !    WRITE(*,'(A,I4,A)') 'Line(',iLine,') = "'//TRIM(ALine)//'"'
         !
         ! Find the type of the record
         !
         CALL CheckRecordType(iNuclide,iLine,ALine)

         IF (iLine.GT.1) THEN
            IsNewChapter = (iRecordType(iLine).NE.iRecordType(iLine-1))
         ENDIF

         IF (IsNewChapter) THEN
            iLineInChapter = 1
         ELSE
            iLineInChapter = iLineInChapter + 1
         ENDIF
         !
         ! Interpret the record on the line
         !
         CALL Parse1DataLine(iLine,iNuclide,ALine,iLineInChapter,Ready,ThisLineHasDecay,iDecay)

         FoundMeaningfulDecay = FoundMeaningfulDecay .OR. ThisLineHasDecay

         IF (Ready.AND.(.NOT.FoundMeaningfulDecay)) THEN
            IF (NuclideSpecs(iNuclide)%HalfTime.LT.1.E90_Float) THEN
               WRITE(*,'(A,I4,A,A,A,EN20.10,A)') 'Nuclide ',iNuclide,': ',NuclideSpecs(iNuclide)%NuclideName,&
               & ' has halflife ',NuclideSpecs(iNuclide)%HalfTime,&
               & ', but I found end of file before I encountered any relevant decay info...'
            ELSE
               IF (DebugLevel.GT.0) WRITE(*,'(A,I4,A,A,A)') 'Nuclide ',iNuclide,': ',&
               & NuclideSpecs(iNuclide)%NuclideName,' is stable!'
            ENDIF
         ENDIF

      ENDDO ! loop until ready

      10 CONTINUE

      CLOSE(ScratchFile)
   END SUBROUTINE ParseENDFFile



   SUBROUTINE ParseSpontaneousFile(FName)
      !
      ! Parse an ENDF file with branching ratio's for spontaneous fission.
      ! Unknown yet is if the branching ratio's account for the yield or not...
      !
      CHARACTER(*), INTENT(IN) :: FName

      LOGICAL :: Ready,IsSFLine,DoingBullocks,PrevLineIsSF,IsNonVoid,LastFissionProductFound
      INTEGER :: iLine,SFLine,iMother,iDaughter,iParameter,NExpectedFissionProducts,NFissionProducts,iElement
      CHARACTER(DefaultLength) :: ALine,DumString
      CHARACTER(11), DIMENSION(12) :: MySubString
      LOGICAL, DIMENSION(12) :: SubStringIsInteger
      INTEGER, DIMENSION(12) :: MyInteger
      INTEGER, DIMENSION(3) :: SFNuclide
      REAL(Float), DIMENSION(12) :: MyFloat
      REAL(Float) :: TotalYield
      REAL(Float), DIMENSION(3) :: MyBranchingRatio

      INTEGER, PARAMETER :: DebugLevel = 0

      IF (.NOT.FileExists(TRIM(ENDFSpontaneousPath)//TRIM(FName))) THEN
         WRITE(*,'(A)') 'Cannot find file "'//TRIM(FName)//'" for spontaneous fission! Exiting!!'
         CALL EXIT()
      ELSE
         WRITE(*,'(A)') 'Parsing file "'//TRIM(FName)//'" for spontaneous fission!'
      ENDIF
      !
      ! Read contents
      !
      OPEN(ScratchFile,FILE=TRIM(ENDFSpontaneousPath)//TRIM(FName),FORM='FORMATTED',ACTION='READ')

      Ready = .FALSE.
      iLine = 0
      SFLine = 0
      TotalYield = 0._Float
      PrevLineIsSF = .FALSE.
      NFissionProducts = 0
      LastFissionProductFound = .FALSE.

      DO WHILE (.NOT.Ready)

         READ(ScratchFile,'(A)',END=10) ALine
         iLine = iLine + 1

         IsSFLine = (ALine(72:75).EQ.'8454')

         IF (.NOT.IsSFLine) THEN
            IF (PrevLineIsSF) THEN
               WRITE(*,'(A,A10,A,F15.10)') 'Sum of branching ratios for spontaneous fission of ',&
               & NuclideSpecs(iMother)%NuclideName,' is ',TotalYield
               WRITE(*,'(A,A,A,I6,A,I6,A)') 'For ',NuclideSpecs(iMother)%NuclideName,' I expected ',&
               & NExpectedFissionProducts,' fission products, found ',&
               & NFissionProducts,' fission products!'
            ENDIF
            SFLine = 0 ! Reset Spontaneous Fission lines
            TotalYield = 0._Float
            PrevLineIsSF = .FALSE.
            NFissionProducts = 0
            LastFissionProductFound = .FALSE.
         ELSE
            SFLine = SFLine + 1
            !
            ! Read the 6 values and check if they are real/integer
            !
            CALL Line2Values(ALine,MySubString,MyFloat,MyInteger,SubStringIsInteger,DoingBullocks)

            IF (SFLine.EQ.1) THEN
               iMother = Zap2Nuclide(MyInteger(1),MyInteger(2))

               WRITE(*,'(A,I7,A,I4,A,A)') 'Head of chain for this list of spontaneous '&
                  & //'fissions, starting at line ',iLine,' = nuclide ',&
                  & iMother,': ',NuclideSpecs(iMother)%NuclideName

            ELSE IF (SFLine.EQ.2) THEN
               NExpectedFissionProducts = MyInteger(6)
            ELSE IF (SFLine.GT.2) THEN
               READ(ScratchFile,'(A)',END=10) ALine
               iLine = iLine + 1
               ! IF (DebugLevel.GT.2) WRITE(*,'(A,I6,A)') 'Line ',iLine,' : "'//TRIM(ALine)//'"'

               CALL Line2Values(ALine,MySubstring(7),MyFloat(7),MyInteger(7),SubStringIsInteger(7),DoingBullocks)
               !
               ! Work out 3 nuclides per 2 lines:
               !
               DO iElement = 1,3
                  IF (.NOT.LastFissionProductFound) THEN
                     IsNonVoid = (LEN_TRIM(MySubstring(4*iElement-3)).GT.0)
                     IF (IsNonVoid) THEN
                        SFNuclide(iElement) = Zap2Nuclide(MyInteger(4*iElement-3),MyInteger(4*iElement-2))
                        MyBranchingRatio(iElement) = MyFloat(4*iElement-1)

                        TotalYield = TotalYield + MyBranchingRatio(iElement)
                        NFissionProducts = NFissionProducts + 1

                        IF ((SFNuclide(iElement).GT.0).AND.(MyBranchingRatio(iElement).GT.0._Float)) THEN

                           WRITE(*,'(A,I8,A,A,A,A,A,EN20.10,1X,I3,1X,I3)') 'Line ',iLine,' SF: ',&
                           & NuclideSpecs(iMother)%NuclideName,' --> ',&
                           & NuclideSpecs(SFNuclide(iElement))%NuclideName,&
                           & ' with branching ratio ',MyBranchingRatio(iElement),&
                           & NuclideSpecs(SFNuclide(iElement))%AtomNumber,NuclideSpecs(SFNuclide(iElement))%NHadrons
                        ENDIF

                        LastFissionProductFound = (NFissionProducts.GE.NExpectedFissionProducts)
                     ENDIF ! nonvoid
                  ENDIF ! not yet last fission product found
               ENDDO ! loop over elements on 2 lines

               PrevLineIsSF = .TRUE.

            ENDIF ! SF line 1, 2 or later

         ENDIF ! line with SF data
      ENDDO ! loop until ready

      10 CONTINUE
   END SUBROUTINE ParseSpontaneousFile



   SUBROUTINE MakeMotherDaughterMatrix(MyNuclideSpecs,MyMotherDaughterMatrix,OutName)
      !
      ! Fill the sparse matrix of which the exponential gives the evolution of the
      ! nuclide strengths (in Bq possibly per unit of volume or surface)
      ! over a macroscopic timestep. The matrix is stored in a global variable,
      ! since it will be needed quite often.
      !
      TYPE(NuclideType), DIMENSION(0:MaxNuclides), INTENT(IN) :: MyNuclideSpecs
      REAL(Float), DIMENSION(MaxNuclides,MaxNuclides), INTENT(OUT) :: MyMotherDaughterMatrix
      CHARACTER(*), INTENT(IN) :: OutName
      INTEGER :: DaughterNuclide,MotherNuclide,NNonZero,HerDaughterNuclide,iDaughter,iGrandDaughter
      CHARACTER(30) :: FormatString
      REAL(Float) :: Ln2

      INTEGER, PARAMETER :: DebugLevel = 0

      IF (DebugLevel.GT.0) THEN
         WRITE(*,'(A,I4)') 'MakeMyMotherDaughterMatrix: NNuclides = ',NNuclides
      ENDIF

      ln2 = LOG(2.)
      !
      ! Initialize diagonal with -1, rest at zero
      !
      DO MotherNuclide = 1,NNuclides
         DO DaughterNuclide = 1,NNuclides
            IF (DaughterNuclide.EQ.MotherNuclide) THEN
               MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide) = -ln2&
                  & /MyNuclideSpecs(DaughterNuclide)%HalfTime
            ELSE
               MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide) =  0.
            ENDIF
         ENDDO
      ENDDO
      !
      ! Now fill the nonzero off-diagonal elements. Leave open the option that some daughters are accounted for several times.
      !
      DO MotherNuclide = 1,NNuclides
         DO iDaughter = 1,MaxNDaughters
            DaughterNuclide = MyNuclideSpecs(MotherNuclide)%Daughter(iDaughter)
            IF (DaughterNuclide.NE.0) THEN
               MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide) = &
               MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide) &
               & + MyNuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter)*ln2&
               &   /MyNuclideSpecs(DaughterNuclide)%HalfTime
            ENDIF
         ENDDO
      ENDDO
      !
      ! Write to output
      !
      IF (DebugLevel.GT.0) WRITE(*,*)
      IF (LEN_TRIM(OutName).GT.0) THEN
         OPEN(ScratchFile,FILE = OutName,FORM='FORMATTED', POSITION='REWIND')
         WRITE(ScratchFile,'(A)') 'Mother/daughter Matrix: horizontal nuclides decay into the vertical nuclides'
         WRITE(ScratchFile,664) NNuclides
         664  FORMAT(I4,' nuclides are represented in this matrix. The matrix elements are for the decay of ACTIVITY!!')

         FormatString = '(11X,    (A10,1X))'
         WRITE(FormatString(6:9),'(I4)') NNuclides
         WRITE(ScratchFile,FormatString) (MyNuclideSpecs(MotherNuclide)%NuclideName,MotherNuclide = 1,NNuclides)

         FormatString = '(A10,1X,    (G10.4,1X))'
         WRITE(FormatString(9:12),'(I4)') NNuclides
         DO DaughterNuclide = 1,NNuclides
            WRITE(ScratchFile,FormatString) MyNuclideSpecs(DaughterNuclide)%NuclideName,&
         & (MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide),MotherNuclide = 1,NNuclides)
         ENDDO
         !
         ! Count number of non-zero entries in decay matrix
         !
         NNonZero = 0
         DO DaughterNuclide = 1,NNuclides
            DO MotherNuclide = 1,NNuclides
               IF (ABS(MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide)).GT.1.E-24) CALL Inc(NNonZero)
            ENDDO
         ENDDO
         !
         ! Write only non-zero entries to output
         !
         WRITE(ScratchFile,*)
         WRITE(ScratchFile,'(A)') 'Here come the non-zero entries of the mother/daughter matrix:'
         WRITE(ScratchFile,668) NNuclides,NNonZero
         668  FORMAT(I10,1X,I10,'  = NNuclides and NNonZero')
         DO DaughterNuclide = 1,NNuclides
            DO MotherNuclide = 1,NNuclides
               IF (ABS(MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide)).GT.1.E-24) THEN
                  WRITE(ScratchFile,667) &
                     & DaughterNuclide,MotherNuclide,MyMotherDaughterMatrix(DaughterNuclide,MotherNuclide)
                     667 FORMAT(I4,1X,I4,1X,G10.4)
               ENDIF
            ENDDO
         ENDDO
         CLOSE(ScratchFile)
      ENDIF
   END SUBROUTINE MakeMotherDaughterMatrix



   SUBROUTINE ReadENDFNuclideSpecs()
      !
      ! Get all nuclide stuff you want to have from the ENDF database
      !
      INTEGER :: iNuclide

      CALL InitLibENDF()

      DO iNuclide = 1,NNuclides
         CALL ParseENDFFile(iNuclide)
      ENDDO ! loop over nuclides
      !
      ! Spontaneous fission is skippped, as for none of these nuclides there are meaningful dose conversion coefficients.
      ! You can give it a go by uncommenting the commands below. No guarantee that it works.
      !
      !   CALL ParseSpontaneousFile('nt501')
      !   CALL ParseSpontaneousFile('nt502')
      !   CALL ParseSpontaneousFile('nt503')
      !   CALL ParseSpontaneousFile('nt505')
      !   CALL ParseSpontaneousFile('nt506')
      !   CALL ParseSpontaneousFile('nt507')
      !   CALL ParseSpontaneousFile('nt508')
      !   CALL ParseSpontaneousFile('nt509')
      !   CALL ParseSpontaneousFile('nt510')
      !   CALL ParseSpontaneousFile('nt511')
      !   CALL ParseSpontaneousFile('nt512')
      !   CALL ParseSpontaneousFile('nt513')
      !   CALL ParseSpontaneousFile('nt514')
      !   CALL ParseSpontaneousFile('nt515')
      !   CALL ParseSpontaneousFile('nt516')
      !   CALL ParseSpontaneousFile('nt517')
      !   CALL ParseSpontaneousFile('nt518')
      !   CALL ParseSpontaneousFile('nt666a')
      !   CALL ParseSpontaneousFile('nt666b')
      !   CALL ParseSpontaneousFile('nt998')
      !   CALL ParseSpontaneousFile('nt999')
      !
      ! Match the names of the daughters with the other nuclides.
      ! Exclude stable daughers
      !
      CALL MatchDaughters()
      WRITE(*,'(A)') 'Ready matching daughters!'
   END SUBROUTINE ReadENDFNuclideSpecs



   SUBROUTINE ReadNProcessENDFNuclideSpecs(tMin,RegularizedNuclideSpecs,RegularizedMotherDaughterMatrix)
      !
      ! Read ENDF nuclide database, find orphans and construct mother-daughter matrices
      !
      REAL(Float), INTENT(INOUT) :: tMin
      TYPE(NuclideType), DIMENSION(0:MaxNuclides), INTENT(INOUT) :: RegularizedNuclideSpecs
      REAL(Float), DIMENSION(MaxNuclides,MaxNuclides), INTENT(INOUT) :: RegularizedMotherDaughterMatrix

      TYPE(SparseMatrix) :: SparseMMMatrix,SparseRegularizedMMMatrix
      REAL(Float) :: tMinSug
      LOGICAL :: ENDFBinFileExists,IsDifferenttMin
      TYPE(SparseLogicalMatrix) :: SparseIsRelated,SparseRegularizedIsRelated
      INTEGER :: MostDaughters,iNuclide
      CHARACTER(DefaultLength) :: FName

      INTEGER, PARAMETER :: DebugLevel = 0

      WRITE(*,'(A)') 'Reading the ENDF database of nuclides...'

      IF (.NOT.FileExists('ICRP-07.NDX')) THEN
         WRITE(*,'(A)') 'Make sure you have file ICRP-07.NDX in the directory where you call this utility! Exiting!!'
         CALL EXIT()
      ENDIF
      !
      ! Read all data on nuclides from an earlier run or construct them afresh
      !
      ENDFBinFileExists = FileExists('ENDFBinFile.dat')

      IF (ENDFBinFileExists) THEN

         WRITE(*,'(A)') 'I found the binary file ENDFBinFile.dat with all parsed ENDF data readily available!'
         WRITE(*,'(A)') 'This is not the first time that you run this library!'
         WRITE(*,'(A)') ''
         WRITE(*,'(A)') 'Initialization will go smooth and fast!'
         WRITE(*,'(A)') 'Going to read this preprocessed file ...'

         OPEN(ScratchFile,FILE='ENDFBinFile.dat',FORM='UNFORMATTED',POSITION='REWIND',ACTION='READ')

         READ(ScratchFile) AtomSpecs
         READ(ScratchFile) AtomName
         READ(ScratchFile) NNuclides
         READ(ScratchFile) NuclideSpecs

         IF (ALLOCATED(SparseMMMatrix%Element)) DEALLOCATE(SparseMMMatrix%Element)
         READ(ScratchFile) SparseMMMatrix%N,SparseMMMatrix%NMax
         ALLOCATE(SparseMMMatrix%Element(SparseMMMatrix%N))
         READ(ScratchFile) SparseMMMatrix%Element
         CALL SparseMatrix2Matrix(SparseMMMatrix,MotherDaughterMatrix)

         READ(ScratchFile) SparseIsRelated%N
         ALLOCATE(SparseIsRelated%Element(SparseIsRelated%N))
         READ(ScratchFile) SparseIsRelated%Element
         CALL SparseLMatrix2LMatrix(SparseIsRelated,IsRelated)

         READ(ScratchFile) tMinSug
         IsDifferenttMin = (ABS(tMinSug/tMin-1._Float).GT.0.001_Float)
         IF (IsDifferenttMin) THEN
            WRITE(*,'(A,EN20.10,A,EN20.10,A)') 'The value of tMin in the file of ',tMinSug,&
            & ' seconds differs from the one in your subroutine call: ',tMin,' seconds!'
            WRITE(*,'(A)') 'Going to replace your value by the one in tne file!'
            tMin = tMinSug
         ENDIF

         READ(ScratchFile) RegularizedNuclideSpecs

         IF (ALLOCATED(SparseRegularizedMMMatrix%Element)) DEALLOCATE(SparseRegularizedMMMatrix%Element)
         READ(ScratchFile) SparseRegularizedMMMatrix%N,SparseRegularizedMMMatrix%NMax
         ALLOCATE(SparseRegularizedMMMatrix%Element(SparseRegularizedMMMatrix%N))
         READ(ScratchFile) SparseRegularizedMMMatrix%Element
         CALL SparseMatrix2Matrix(SparseRegularizedMMMatrix,RegularizedMotherDaughterMatrix)

         READ(ScratchFile) SparseRegularizedIsRelated%N
         ALLOCATE(SparseRegularizedIsRelated%Element(SparseRegularizedIsRelated%N))
         READ(ScratchFile) SparseRegularizedIsRelated%Element
         CALL SparseLMatrix2LMatrix(SparseRegularizedIsRelated,RegularizedIsRelated)

         CLOSE(ScratchFile)

         WRITE(*,'(A)') 'Ready reading preprocessed ENDF file ...'
      ELSE
         WRITE(*,'(A)') 'I could not find the binary file ENDFBinFile.dat with all parsed ENDF data readily available...'
         WRITE(*,'(A)') 'This is probably the first time that you run this library...'
         WRITE(*,'(A)') ''
         WRITE(*,'(A)') 'Going to read all the ENDF data files one by one and parse them ...'
         WRITE(*,'(A)') 'This may take some time ...'

         CALL ReadENDFNuclideSpecs()
         !
         ! Find all orphan decay steps, i.e. those involving a halflife < tMin, and put them in a list
         !
         CALL FindOrphans(tMin)
         !
         ! Make a copy of the nuclide specs with which we can play
         ! The original is kept such that we will be able to trace the true family relations
         !
         RegularizedNuclideSpecs = NuclideSpecs
         !
         ! Regularize nuclide specs by elimination of too fast nuclides
         !
         CALL RegularizeNuclides(RegularizedNuclideSpecs)

         WRITE(*,*)
         !================================================================================================================
         !
         ! Construct matrix with time-rate-of-change of nuclide vector, based on the regularized transition characteristics
         !
         WRITE(*,*)
         WRITE(*,'(A)') 'Going to construct mother-daughter matrix for time-rate-of-change!'
         !
         ! Once for the original full set of decay options
         !
         FName = ' '
         CALL MakeMotherDaughterMatrix(NuclideSpecs,MotherDaughterMatrix,FName)
         !
         ! And once for the regularized set of decay relations, where all orphans have been taken out
         !
         FName = ' '
         CALL MakeMotherDaughterMatrix(RegularizedNuclideSpecs,RegularizedMotherDaughterMatrix,FName)

         WRITE(*,'(A)') 'Ready reading and parsting all the ENDF data files!'

         OPEN(ScratchFile,FILE='ENDFBinFile.dat',FORM='UNFORMATTED',POSITION='REWIND',ACTION='WRITE')

         WRITE(ScratchFile) AtomSpecs
         WRITE(ScratchFile) AtomName
         WRITE(ScratchFile) NNuclides
         WRITE(ScratchFile) NuclideSpecs
         SparseMMMatrix = Matrix2SparseMatrix(MotherDaughterMatrix)
         WRITE(ScratchFile) SparseMMMatrix%N,SparseMMMatrix%NMax
         WRITE(ScratchFile) SparseMMMatrix%Element
         SparseIsRelated = LMatrix2SparseLMatrix(IsRelated)
         WRITE(ScratchFile) SparseIsRelated%N
         WRITE(ScratchFile) SparseIsRelated%Element
         WRITE(ScratchFile) tMin
         WRITE(ScratchFile) RegularizedNuclideSpecs
         SparseRegularizedMMMatrix = Matrix2SparseMatrix(RegularizedMotherDaughterMatrix)
         WRITE(ScratchFile) SparseRegularizedMMMatrix%N,SparseRegularizedMMMatrix%NMax
         WRITE(ScratchFile) SparseRegularizedMMMatrix%Element
         SparseRegularizedIsRelated = LMatrix2SparseLMatrix(RegularizedIsRelated)
         WRITE(ScratchFile) SparseRegularizedIsRelated%N
         WRITE(ScratchFile) SparseRegularizedIsRelated%Element

         CLOSE(ScratchFile)

         WRITE(*,'(A)') 'Created binary file ENDFBinFile.dat with all parsed ENDF data for future use!'

      ENDIF ! Binary data already available

      IF (DebugLevel.GT.0) THEN
         WRITE(*,'(A,I0,A)') 'AtomSpecs                  takes : ',SIZEOF(AtomSpecs),' bytes'
         WRITE(*,'(A,I0,A)') 'AtomName                   takes : ',SIZEOF(AtomName),' bytes'
         WRITE(*,'(A,I0,A)') 'NNuclides                  takes : ',SIZEOF(NNuclides),' bytes'
         WRITE(*,'(A,I0,A)') 'NuclideSpecs               takes : ',SIZEOF(NuclideSpecs),' bytes'
         WRITE(*,'(A,I0,A)') 'SparseMMMatrix             takes : ',SIZEOF(SparseMMMatrix%Element),' bytes'
         WRITE(*,'(A,I0,A)') 'SparseIsRelated            takes : ',SIZEOF(SparseIsRelated%Element),' bytes'
         WRITE(*,'(A,I0,A)') 'tMinSug                    takes : ',SIZEOF(tMinSug),' bytes'
         WRITE(*,'(A,I0,A)') 'RegularizedNuclideSpecs    takes : ',SIZEOF(RegularizedNuclideSpecs),' bytes'
         WRITE(*,'(A,I0,A)') 'SparseRegularizedMMMatrix  takes : ',SIZEOF(SparseRegularizedMMMatrix%Element),' bytes'
         WRITE(*,'(A,I0,A)') 'SparseRegularizedIsRelated takes : ',SIZEOF(SparseRegularizedIsRelated%Element),' bytes'
      ENDIF
   END SUBROUTINE ReadNProcessENDFNuclideSpecs



   SUBROUTINE FindOrphans(tMin)
      !================================================================================================================
      !
      ! Find the orphans: nuclides with halftime shorter than Orphanage%tMin
      !
      !================================================================================================================
      REAL(Float), INTENT(IN) :: tMin

      INTEGER :: DaughterNuclide,MotherNuclide,NMothers,iDaughter,HerDaughterNuclide,iOrphan

      INTEGER, PARAMETER :: DebugLevel = 0

      WRITE(*,*)
      WRITE(*,'(A)') 'Starting to look for orphan nuclides!'

      Orphanage%tMin = tMin
      !
      ! Find orphans in the original nuclide specs dataset (to prevent missing combos of Mother --> Orphanage%Orphan-daughter --> Orphanage%Orphan-granddaughter)
      !
      Orphanage%NOrphans = 0
      Orphanage%NTooUnStable = 0

      DO DaughterNuclide = 1,NNuclides
         Orphanage%IsTooUnstable(DaughterNuclide) = (NuclideSpecs(DaughterNuclide)%HalfTime.LT.Orphanage%tMin)
         IF (Orphanage%IsTooUnstable(DaughterNuclide)) THEN
            IF (DebugLevel.GT.0) WRITE(*,'(A,EN20.10,A,EN20.10,A)') 'Handling too fast nuclide '&
               & //TRIM(NuclideSpecs(DaughterNuclide)%NuclideName)//' with halflife ',NuclideSpecs(DaughterNuclide)%HalfTime,&
               & ' [s] < ',Orphanage%tMin,'[s]'

            NuclideSpecs(DaughterNuclide)%IsOrphan = .TRUE.

            Orphanage%NTooUnstable = Orphanage%NTooUnStable + 1
            !
            ! Check all nuclides if they can be the mother of this daughter
            !
            NMothers = 0

            DO MotherNuclide = 1,NNuclides
               DO iDaughter = 1,MaxNDaughters
                  HerDaughterNuclide = NuclideSpecs(MotherNuclide)%Daughter(iDaughter)
                  !
                  ! If a possible mother has been found, then define a new orphan
                  !
               IF (DaughterNuclide.EQ.HerDaughterNuclide) THEN
                  NMothers = NMothers + 1

                  IF (DebugLevel.GT.0) WRITE(*,'(A,I0,A,F15.10)') &
                     & 'Found mother '//TRIM(NuclideSpecs(MotherNuclide)%NuclideName)//&
                     & ' and '//TRIM(NuclideSpecs(DaughterNuclide)%NuclideName)//' is her daughter number ',iDaughter,&
                     & ' with branching ratio ',NuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter)

                  IF (Orphanage%NOrphans.EQ.MaxNOrphans) THEN
                     WRITE(*,'(A,I0,A)') 'Cannot add new orphan after maximum number of ',&
                        & Orphanage%NOrphans,' orphans! Exiting!!'
                     CALL EXIT()
                  ENDIF ! maximum dimension reached

                  Orphanage%NOrphans = Orphanage%NOrphans + 1

                  Orphanage%Orphan(Orphanage%NOrphans)%Mother   = MotherNuclide
                  Orphanage%Orphan(Orphanage%NOrphans)%Daughter = DaughterNuclide
                  Orphanage%Orphan(Orphanage%NOrphans)%yield    = NuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter)
               ENDIF ! Recognition of daughter by mother
            ENDDO ! loop over daughters
         ENDDO ! Loop over candidate mothers

         IF (DebugLevel.GT.0) THEN
            IF (NMothers.GT.0) THEN
               WRITE(*,'(A,I0,A)') 'Found ',NMothers,' mothers for '//TRIM(NuclideSpecs(DaughterNuclide)%NuclideName)//'!'
            ELSE
               WRITE(*,'(A)') 'No mothers were found, '//TRIM(NuclideSpecs(DaughterNuclide)%NuclideName)&
                  & //' is an artificial nuclide and not a decay product!'
            ENDIF

            WRITE(*,'(A)') 'This is all I wanted to say about too fast nuclide '&
               & //TRIM(NuclideSpecs(DaughterNuclide)%NuclideName)
            ENDIF ! DebugLevel > 0
         ENDIF ! Halflife of daughter too short
      ENDDO ! Loop over daughters
      !
      ! Show list of orphans
      !
      WRITE(*,'(A)') 'List of orphans: '

      DO iOrphan = 1,Orphanage%NOrphans
         WRITE(*,'(A,I4,A,A,A,A,5X,A,EN20.10,A,F20.10,A)') &
            & 'Orphanage%Orphan ',iOrphan,' is daughter ',&
            & NuclideSpecs(Orphanage%Orphan(iOrphan)%Daughter)%NuclideName,&
            & ' from mother ',&
            & NuclideSpecs(Orphanage%Orphan(iOrphan)%Mother  )%NuclideName,&
            & 'because halflife ',&
            & NuclideSpecs(Orphanage%Orphan(iOrphan)%Daughter)%HalfTime,&
            & '[s] < limit ',Orphanage%tMin,'[s]'
      ENDDO ! loop over orphans

      WRITE(*,'(A)') 'Ready looking for orphan nuclides!'
   END SUBROUTINE FindOrphans



   SUBROUTINE RegularizeNuclides(MyNuclideSpecs)
      !
      ! All orphan nuclides with halftime shorter than Orphanage%tMin are eliminated.
      ! Their contributions are attributed to the mother nuclide decaying directly to the grand-daughters.
      !
      TYPE(NuclideType), DIMENSION(0:MaxNuclides), INTENT(INOUT) :: MyNuclideSpecs

      INTEGER :: DaughterNuclide,NGrandDaughters,iGrandDaughter,GrandDaughterNuclide,NMothers,MotherNuclide,iDaughter,&
      & HerDaughterNuclide,FreeDaughter,iNuclide
      LOGICAL :: IsAlreadyADaughter,IsFirstFreeSlot

      INTEGER, PARAMETER :: DebugLevel = 0

      WRITE(*,*)
      WRITE(*,'(A,EN20.10,A)') 'Going to eliminate all nuclides with halflife faster than ',Orphanage%tMin,' s'
      WRITE(*,*)
      !
      ! "DaughterNuclide" is the possibly too fast nuclide, i.e. an orphan nuclide, that has to be brushed under the carpet.
      !
      DO DaughterNuclide = 1,NNuclides
         IF (Orphanage%IsTooUnstable(DaughterNuclide)) THEN
            IF (DebugLevel.GT.0) WRITE(*,'(A,EN20.10,A)') 'Handling too fast nuclide '&
               & //TRIM(MyNuclideSpecs(DaughterNuclide)%NuclideName)//' with halflife ',&
               & MyNuclideSpecs(DaughterNuclide)%HalfTime,&
               & ' s'
            !
            ! Count number of granddaughters that is formed and that we have to re-arrange
            !
            WRITE(*,'(A)') 'Going to search for radioactive granddaughters...'
            NGrandDaughters = 0
            DO iGrandDaughter = 1,MaxNDaughters
               GrandDaughterNuclide = MyNuclideSpecs(DaughterNuclide)%Daughter(iGrandDaughter)
               IF (GrandDaughterNuclide.GT.0) THEN
                  WRITE(*,'(A)') 'Spotted radioactive granddaughter '&
                     & //MyNuclideSpecs(GrandDaughterNuclide)%NuclideName
                  NGrandDaughters = NGrandDaughters + 1
               ENDIF
            ENDDO
            IF (NGrandDaughters.GT.0) THEN
               WRITE(*,'(A,I0,A)') 'This daughter nuclide forms ',NGrandDaughters,' radioactive granddaughters'
            ELSE
               WRITE(*,'(A,I0,A)') 'This daughter nuclide forms no radioactive granddaughters, '&
               & //'I will only eliminate the formation of this nuclide'
            ENDIF
            !
            ! Check all nuclides if they can be the mother of this orphan daughter
            !
            ! "MotherNuclide" is the mother of the orphan, whose decay to the orphan will be shifted to her dranddaughters
            !
            NMothers = 0

            DO MotherNuclide = 1,NNuclides
               DO iDaughter = 1,MaxNDaughters
                  HerDaughterNuclide = MyNuclideSpecs(MotherNuclide)%Daughter(iDaughter)
                  !
                  !#############################################################################################################
                  !
                  ! If a possible mother has been found, eliminate the family relation and pass all info on to the next generation
                  !
                  IF (DaughterNuclide.EQ.HerDaughterNuclide) THEN

                     NMothers = NMothers + 1

                     WRITE(*,'(A,I0,A,F15.10)') 'Found mother '//TRIM(MyNuclideSpecs(MotherNuclide)%NuclideName)//&
                        & ' and '//TRIM(MyNuclideSpecs(DaughterNuclide)%NuclideName)//' is her daughter number ',&
                        & iDaughter,&
                        & ' with branching ratio ',MyNuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter)
                     !
                     ! Pass heritage to granddaughters via daughter
                     !
                     WRITE(*,'(A)') 'Going to pass heritage of mother '&
                        & //TRIM(MyNuclideSpecs(MotherNuclide)%NuclideName)&
                        & //' to her radioactive granddaughters...'

                     DO iGrandDaughter = 1,MyNuclideSpecs(DaughterNuclide)%NDaughters
                        GrandDaughterNuclide = MyNuclideSpecs(DaughterNuclide)%Daughter(iGrandDaughter)
                        IF (GrandDaughterNuclide.GT.0) THEN
                           !
                           ! Found a grand-daughter
                           !
                           WRITE(*,'(A,I0,A,F15.10)') 'Found granddaughter '&
                              & //TRIM(MyNuclideSpecs(GrandDaughterNuclide)%NuclideName)&
                              &   //' with index ',iGrandDaughter,' and with branching ratio ',&
                              & MyNuclideSpecs(DaughterNuclide)%DaughterFraction(iGrandDaughter)
                           !
                           ! Check if this granddaughter is already a daughter of the grandmother
                           !
                           IsAlreadyADaughter = .FALSE.
                           FreeDaughter = 0
                           DO iNuclide = 1,MaxNDaughters
                              IF (MyNuclideSpecs(MotherNuclide)%Daughter(iNuclide).EQ.GrandDaughterNuclide) THEN
                                 FreeDaughter = iNuclide
                                 IsAlreadyADaughter = .TRUE.
                              ENDIF ! Granddaughter is already also a daughter
                           ENDDO
                           !
                           ! If not already a daughter, find the first available new daughter slot
                           !
                           IF (FreeDaughter.EQ.0) THEN
                              DO iNuclide = 1,MaxNDaughters
                                 IsFirstFreeSlot=((FreeDaughter.EQ.0).AND.(MyNuclideSpecs(MotherNuclide)%Daughter(iNuclide).EQ.0))
                                 IF (IsFirstFreeSlot) FreeDaughter = iNuclide
                              ENDDO
                           ENDIF ! granddaughter was not also a daughter

                           IF (FreeDaughter.GT.0) THEN
                              IF (IsAlreadyADaughter) THEN
                                 WRITE(*,'(A,I0,A)') 'Daughter slot ',FreeDaughter,&
                                 & ' of mother '//TRIM(MyNuclideSpecs(MotherNuclide)%NuclideName)&
                                 & //' already contained granddaughter '//TRIM(MyNuclideSpecs(GrandDaughterNuclide)%NuclideName)&
                                 & //'! Adding contribution!'

                                 WRITE(*,'(A,I0,A)') 'Augmenting branching ratio for old daughter '&
                                 & //TRIM(MyNuclideSpecs(GrandDaughterNuclide)%NuclideName)&
                                 & //' of mother '//TRIM(MyNuclideSpecs(MotherNuclide)%NuclideName)//' in daughter slot number ',&
                                 & FreeDaughter,':'
                                 WRITE(*,'(5(A,F15.10))') &
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

                              ELSE
                                 WRITE(*,'(A,I0,A)') 'Daughter slot ',FreeDaughter,&
                                 & ' of mother is still free, adopting it for granddaughter!'

                                 WRITE(*,'(A,I0)') 'Installing new daughter '&
                                 & //TRIM(MyNuclideSpecs(GrandDaughterNuclide)%NuclideName)&
                                 & //' of mother '//TRIM(MyNuclideSpecs(MotherNuclide)%NuclideName)//' in daughter slot number ',&
                                 & FreeDaughter

                                 MyNuclideSpecs(MotherNuclide)%Daughter(FreeDaughter) = &
                                 & GrandDaughterNuclide

                                 MyNuclideSpecs(MotherNuclide)%DaughterName(FreeDaughter) = &
                                 & MyNuclideSpecs(DaughterNuclide)%DaughterName(iGrandDaughter)

                                 MyNuclideSpecs(MotherNuclide)%DaughterFraction(FreeDaughter) = &
                                 &  MyNuclideSpecs(MotherNuclide  )%DaughterFraction(iDaughter     )&
                                 & *MyNuclideSpecs(DaughterNuclide)%DaughterFraction(iGrandDaughter)

                              ENDIF
                           ELSE
                              WRITE(*,'(A)') 'Somehow I cannot find a free daughter...'
                              CALL EXIT()
                           ENDIF

                        ENDIF
                     ENDDO ! loop over granddaughters
                     !
                     ! Eliminate the decay from this mother to the daughter at hand in the nuclide database:
                     !
                     WRITE(*,'(A,I0,A)') 'Eliminating daughter slot ',iDaughter,' of mother!'

                     MyNuclideSpecs(MotherNuclide)%Daughter(iDaughter) = 0
                     MyNuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter) = 0._Float
                     MyNuclideSpecs(MotherNuclide)%DaughterName(iDaughter) = ' '

                  ENDIF ! Recognition of daughter by mother
                  !
                  ! End of bookkeeping of the removal of 1 orphan nuclide by adding its contributions to 1 of its mothers
                  !
                  !#############################################################################################################
                  !
               ENDDO ! loop over 20 daughters
            ENDDO ! Loop over candidate mothers

            IF (NMothers.GT.0) THEN
               WRITE(*,'(A,I0,A)') 'Found ',NMothers,' mothers!'
            ELSE
               WRITE(*,'(A)') 'No mothers were found, this is an artificial nuclide and not a decay product!'
            ENDIF

            WRITE(*,'(A)')'Ready handling too fast nuclide '//TRIM(MyNuclideSpecs(DaughterNuclide)%NuclideName)
            WRITE(*,*)
         ENDIF ! Halflife of daughter too short
      ENDDO ! Loop over daughters
      !
      ! Make IsRelated(Daughter,Mother) matrix for the regularized set
      !
      RegularizedIsRelated = .FALSE.

      DO MotherNuclide = 1,NNuclides
         DO iDaughter = 1,MaxNDaughters
            HerDaughterNuclide = MyNuclideSpecs(MotherNuclide)%Daughter(iDaughter)
            IF (HerDaughterNuclide.NE.0) THEN
               RegularizedIsRelated(HerDaughterNuclide,MotherNuclide) = .TRUE.
            ENDIF
         ENDDO
      ENDDO

      WRITE(*,'(A)') 'Ready regularizing mother-daughter matrix!'
   END SUBROUTINE RegularizeNuclides



   SUBROUTINE CollectProgeny(iNuclide,DoPrint)
      !================================================================================================================
      !
      ! Collect all progeny of this nuclide, first the regular ones, then add the short-living ones.
      !
      !================================================================================================================
      INTEGER, INTENT(IN) :: iNuclide
      LOGICAL, INTENT(IN) :: DoPrint

      INTEGER :: NAdded,MotherNuclide,DaughterNuclide,OrphanCounter,SwiftNuclide,NFarAdded
      LOGICAL :: Ready,IsNewFamily,IsNewFarFamily

      INTEGER, PARAMETER :: DebugLevel = 0

      IF (DoPrint) WRITE(*,'(A,I0,A,EN20.10,A)') 'Collecting the progeny family for nuclide ',iNuclide,&
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

      DO WHILE (.NOT.Ready)
         NAdded = 0
         NFarAdded = 0
         !
         ! Loop over possible mothers
         !
         DO MotherNuclide = 1,NNuclides
            !
            ! Loop over possible daughters
            !
            DO DaughterNuclide = 1,NNuclides
               !
               ! First check if this nuclide is part of the full family tree, including orphans
               ! The search is downward in the tree, so new members are identified if the mother is already a member,
               ! but the daughter not yet.
               !
               IsNewFarFamily =  (IsRelated(DaughterNuclide,MotherNuclide) &
               &    .AND. NuclideFamily%IsFarFamily(MotherNuclide) &
               &    .AND. (.NOT.NuclideFamily%IsFarFamily(DaughterNuclide))&
               &    .AND. (NuclideSpecs(DaughterNuclide)%HalfTime.LT.1.E80_Float))

               IF (IsNewFarFamily) THEN
                  NuclideFamily%IsFarFamily(DaughterNuclide) = .TRUE.
                  NFarAdded = NFarAdded + 1
                  NuclideFamily%NFarFamily = NuclideFamily%NFarFamily + 1
               ENDIF
               !
               ! Now check if it is part of the narrow family of non-orphans and update bookkeeping if affirmative
               !
               IsNewFamily =  (RegularizedIsRelated(DaughterNuclide,MotherNuclide) &
               &    .AND. NuclideFamily%IsFamily(MotherNuclide) &
               &    .AND. (.NOT.NuclideFamily%IsFamily(DaughterNuclide))&
               &    .AND. (NuclideSpecs(DaughterNuclide)%HalfTime.LT.1.E80_Float))
               !
               ! Add new family member
               !
               IF (IsNewFamily) THEN
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

                  IF (DebugLevel.GT.2) THEN
                     WRITE(*,'(A)') 'Adding daughter '//NuclideSpecs(DaughterNuclide)%NuclideName//' for mother '//&
                     & NuclideSpecs(MotherNuclide)%NuclideName
                  ENDIF
               ELSE IF (IsNewFarFamily) THEN
                  IF (DebugLevel.GT.2) THEN
                     WRITE(*,'(A)') 'Adding orphan daughter '//NuclideSpecs(DaughterNuclide)%NuclideName//' for mother '//&
                     & NuclideSpecs(MotherNuclide)%NuclideName
                  ENDIF
               ENDIF
            ENDDO ! loop over all other nuclides
         ENDDO ! loop over nuclides already in the family

         Ready = (NFarAdded .EQ. 0)

         IF (DebugLevel.GT.2) THEN
            IF (Ready) THEN
               WRITE(*,'(A)') 'Found no new members, halting iterations!'
            ELSE
               WRITE(*,'(A,I0,A,I0,A)') 'Found ',NFarAdded,' new members, of which ',NFarAdded-NAdded,' core members!'&
               & //' Proceeding with next iteration!'
            ENDIF
         ENDIF
      ENDDO ! while loop
      !
      ! Now add short living family members
      !
      IF (DebugLevel.GT.2) THEN
         WRITE(*,'(A)') 'The regular family has been formed, now going to add possible fast nuclides!'
      ENDIF

      OrphanCounter = NuclideFamily%NFamily

      NAdded = 0

      DO SwiftNuclide = 1,NNuclides
         IF (NuclideFamily%IsFarFamily(SwiftNuclide).AND..NOT.(NuclideFamily%IsFamily(SwiftNuclide))) THEN
            IF (DebugLevel.GT.2) THEN
               WRITE(*,'(A)') 'Nuclide '//TRIM(NuclideSpecs(SwiftNuclide)%NuclideName)//' seems interesting! Adding it!'
            ENDIF

            NAdded = NAdded + 1

            OrphanCounter = OrphanCounter + 1

            NuclideFamily%FamilyMember(OrphanCounter) = SwiftNuclide
            NuclideFamily%FamilyNumber(SwiftNuclide) = OrphanCounter
            !
            ! Register its daughters
            !
            DO DaughterNuclide = 1,NNuclides
               IF (IsRelated(DaughterNuclide,SwiftNuclide)) THEN
               NuclideFamily%NDaughters(SwiftNuclide) = NuclideFamily%NDaughters(SwiftNuclide) + 1
               NuclideFamily%Daughter(SwiftNuclide,NuclideFamily%NDaughters(SwiftNuclide)) = &
               &   DaughterNuclide
               ENDIF
            ENDDO ! loop over possible mothers
            !
            ! Register its mothers
            !
            DO MotherNuclide = 1,NNuclides
               IF (IsRelated(SwiftNuclide,MotherNuclide)) THEN
               NuclideFamily%NMothers(SwiftNuclide) = NuclideFamily%NMothers(SwiftNuclide) + 1
               NuclideFamily%Mother(SwiftNuclide,NuclideFamily%NMothers(SwiftNuclide)) = &
               & MotherNuclide
               ENDIF
            ENDDO ! loop over possible mothers
         ENDIF ! This is such a swift nuclide
      ENDDO ! loop over possibly fast decaying nuclides that have to be added to the family, but with high index numbers

      Ready = (NAdded .EQ. 0)

      IF (DebugLevel.GT.2) THEN
         WRITE(*,'(A,I0,A)') 'Found ',NAdded,' new too fast (orphan) members!'
      ENDIF
      !
      ! Show all relatives for this nuclide
      !
      IF (DebugLevel.EQ.1) THEN
         WRITE(*,'(A)') NuclideSpecs(iNuclide)%NuclideName//' chain:'
         DO DaughterNuclide = 1,NNuclides
            IF(DaughterNuclide.NE.iNuclide) THEN
               IF(NuclideFamily%IsFamily(DaughterNuclide)) THEN
                  WRITE(*,'(5X,A)') NuclideSpecs(DaughterNuclide)%NuclideName
               ELSE IF (NuclideFamily%IsFarFamily(DaughterNuclide)) THEN
                  WRITE(*,'(5X,A)') NuclideSpecs(DaughterNuclide)%NuclideName//' <-- orphan'
               ENDIF
            ENDIF ! Daughter not iNuclide
         ENDDO
      ELSEIF ((DebugLevel.GT.1).AND. DoPrint) THEN
         WRITE(*,'(A,I0,A)') 'This nuclide has a family with ',NuclideFamily%NFarFamily-1,' more members:'
         DO DaughterNuclide = 1,NNuclides
            IF (DebugLevel.GT.2) THEN
               WRITE(*,'(4A,2L5)') 'CollectProgeny: show all relatives: ',NuclideSpecs(iNuclide)%NuclideName,&
               & ' --> ? ',NuclideSpecs(DaughterNuclide)%NuclideName,&
               & NuclideFamily%IsFamily(DaughterNuclide),NuclideFamily%IsFarFamily(DaughterNuclide)
            ENDIF

            IF(DaughterNuclide.NE.iNuclide) THEN
               IF(NuclideFamily%IsFamily(DaughterNuclide)) THEN
                  WRITE(*,'(A,EN20.10,A)') NuclideSpecs(DaughterNuclide)%NuclideName//' with halflife ',&
                  & NuclideSpecs(DaughterNuclide)%HalfTime,' [s]'
               ELSE IF (NuclideFamily%IsFarFamily(DaughterNuclide)) THEN
                  WRITE(*,'(A,EN20.10,A)') NuclideSpecs(DaughterNuclide)%NuclideName//' with halflife ',&
                  & NuclideSpecs(DaughterNuclide)%HalfTime,' [s]  <-- orphan'
               ENDIF
            ENDIF ! Daughter not iNuclide
         ENDDO
      ENDIF
   END SUBROUTINE CollectProgeny



   SUBROUTINE ListTransitions()
      !
      ! List all transitions that take place within 1 family, including orphans
      !
      INTEGER :: MotherNuclide,DaughterNuclide,iDaughter,MAtom,DAtom

      WRITE(*,*)
      WRITE(*,'(A)') 'Transitions in family with head of chain '&
      & //NuclideSpecs(NuclideFamily%FamilyMember(1))%NuclideName//':'

      DO MotherNuclide = 1,NNuclides
         IF (NuclideFamily%IsFarFamily(MotherNuclide)) THEN
            DO iDaughter = 1,NuclideSpecs(MotherNuclide)%NDaughters
               DaughterNuclide = NuclideSpecs(MotherNuclide)%Daughter(iDaughter)
               IF (DaughterNuclide.NE.0) THEN
               !
               ! Find type of transition
               !
               MAtom = NuclideSpecs(MotherNuclide)%AtomNumber
               DAtom = NuclideSpecs(DaughterNuclide)%AtomNumber

               WRITE(*,'(A,A,A30,A,A,A,F15.10)') &
               & NuclideSpecs(MotherNuclide  )%NuclideName,&
               & ' --- ',TRIM(NuclideSpecs(MotherNuclide)%DecayName(iDaughter)),' --> ',&
               & NuclideSpecs(DaughterNuclide)%NuclideName,&
               & ' fraction ',&
               & NuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter)
               ELSE IF ((LEN_TRIM(NuclideSpecs(MotherNuclide)%DaughterName(iDaughter)).GT.0)&
               & .AND.(TRIM(NuclideSpecs(MotherNuclide)%DaughterName(iDaughter)).NE.'SF')) THEN
               WRITE(*,'(4A,F15.10,A)') &
               & NuclideSpecs(MotherNuclide)%NuclideName,&
               & ' ------------> ',&
               & NuclideSpecs(MotherNuclide)%DaughterName(iDaughter),&
               & ' fraction ',&
               & NuclideSpecs(MotherNuclide)%DaughterFraction(iDaughter),&
               & ' STABLE'

               ENDIF
            ENDDO
         ENDIF ! mother is family
      ENDDO
   END SUBROUTINE ListTransitions



   SUBROUTINE MakeEvolutionMatrix(iNuclide,RegularizedMotherDaughterMatrix,x,NSecondsDelay,zDelay)
      !================================================================================================================
      !
      ! Construct exponential of mother/daughter matrix of this family
      !
      !================================================================================================================
      INTEGER, INTENT(IN) :: iNuclide
      REAL(Float), INTENT(IN) :: NSecondsDelay
      REAL(Float), DIMENSION(MaxNuclides,MaxNuclides), INTENT(IN) :: RegularizedMotherDaughterMatrix
      REAL(dp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: zDelay
      REAL(dp), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: x

      INTEGER :: MotherNuclide,DaughterNuclide,ErrorCode,iSecond,iMinute,iHour,iDay,jNuclide
      REAL(dp), DIMENSION(:,:), ALLOCATABLE :: yDelay
      CHARACTER(DefaultLength) :: MyFormat,FName

      INTEGER, PARAMETER :: DebugLevel = 0

      IF (ALLOCATED(x)) DEALLOCATE(x)
      ALLOCATE(x(NuclideFamily%NFamily,NuclideFamily%NFamily))
      x = 0._Float

      IF (ALLOCATED(yDelay)) DEALLOCATE(yDelay)
      ALLOCATE(yDelay(NuclideFamily%NFamily,NuclideFamily%NFamily))
      yDelay = 0._Float

      IF (ALLOCATED(zDelay)) DEALLOCATE(zDelay)
      ALLOCATE(zDelay(NuclideFamily%NFamily,NuclideFamily%NFamily))
      zDelay = 0._Float
      !
      ! Copy time-rate of change from large matrix for all nuclides to smaller matrix
      !
      DO MotherNuclide = 1,NuclideFamily%NFamily
         DO DaughterNuclide = 1,NuclideFamily%NFamily
            x(DaughterNuclide,MotherNuclide) = &
            & RegularizedMotherDaughterMatrix(NuclideFamily%FamilyMember(DaughterNuclide),&
            & NuclideFamily%FamilyMember(MotherNuclide))
         ENDDO
      ENDDO
      !
      ! Macroscopic time-step over 1 minute, 1 hour and 1 day
      !
      yDelay = NSecondsDelay*x
      !
      ! Exponential
      !
      CALL MatrixExponential(yDelay,NuclideFamily%NFamily,zDelay,ErrorCode)
   END SUBROUTINE MakeEvolutionMatrix



   SUBROUTINE AddOrphansBelow(iNuclide,TheDecayMatrix,zDelay,NFirstMotherOrphans)
      !================================================================================================================
      !
      ! Modify a transition matrix to include the orphan transitions below this mother nuclide that have been short-circuited
      ! from the database plus orphans of orphans etcetera. Only the set of decay products for the head of chain is repaired,
      ! i.e. the line in the transition matrix for element iNuclide.
      !
      !================================================================================================================
      INTEGER, INTENT(IN) :: iNuclide
      REAL(Float), DIMENSION(MaxNuclides,MaxNuclides), INTENT(INOUT) :: TheDecayMatrix
      REAL(dp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: zDelay
      INTEGER, INTENT(OUT) :: NFirstMotherOrphans

      INTEGER :: iOrphan,NAdded,jOrphan,iIteration,iiNuclide,jjNuclide,MotherNuclide,&
      & GhostNuclide,NAddedRegular,NAddedOrphans
      LOGICAL :: Ready,NotAlreadyAMember,IsMotherOfNewOrphan,HaveCandidateMother,ItsMotherIs_ii
      LOGICAL, DIMENSION(MaxNuclides) :: TooFastFullyChecked

      INTEGER, PARAMETER :: DebugLevel = 0

      IF (DebugLevel.GT.0) WRITE(*,*)
      IF (DebugLevel.GT.0) WRITE(*,'(A)') '-------------------- Start of search for orphan chain ------------------------'
      IF (DebugLevel.GT.0) WRITE(*,*)
      IF (DebugLevel.GT.0) WRITE(*,'(A,I4,A,A)') 'Checking orphan chain below nuclide ',iNuclide,&
      & ' : ',TRIM(NuclideSpecs(iNuclide)%NuclideName)

      NFirstMotherOrphans = 0
      Orphanage%Orphan(:)%Active = .FALSE. ! None of the orphan transitions is yet recognized as relevant
      !
      ! Search for too fast progeny in iterations until no new too fast nuclides are found
      !
      Ready = .FALSE.
      iIteration = 0

      DO WHILE (.NOT.Ready)
         iIteration = iIteration + 1
         IF (DebugLevel.GT.0) WRITE(*,*)
         IF (DebugLevel.GT.0) WRITE(*,'(A,I0)') '############# Iteration ',iIteration
         !
         ! If any arrow into a too fast nuclide, either from regular progeny of the head of chain or from another
         ! too fast nuclide in the family, is inactive, switch off this daughter nuclide from TooFastFullyChecked
         !
         CALL MarkFullyCheckedOrphans()
         !
         ! Now fill in the gaps in the transition matrix for too fast nuclide
         !
         CALL AddOrphansFromRegularFamily(NAddedRegular)
         CALL AddOrphansFromCheckedOrphans(NAddedOrphans)

         NAdded = NAddedRegular + NAddedOrphans

         NFirstMotherOrphans = NFirstMotherOrphans + NAdded

         IF (DebugLevel.GT.0) WRITE(*,'(A,I0)') 'Number of orphans added in this iteration: ',NAdded

         Ready = (NAdded.EQ.0)

      ENDDO ! loop until no new orphans added

      IF (DebugLevel.GT.0) WRITE(*,*)
      IF (DebugLevel.GT.0) WRITE(*,'(A)') '-------------------- End of search for orphan chain --------------------------'
      IF (DebugLevel.GT.0) WRITE(*,*)

      CONTAINS

         SUBROUTINE MarkFullyCheckedOrphans
            TooFastFullyChecked = .TRUE. ! De-selection of non-fully checked too fast nuclides will be done immediately below

            DO iOrphan = 1,Orphanage%NOrphans ! This loop includes all orphans, also the irrelevant ones not in this chain
               !
               ! See if its mother is either regular progeny in the chain at hand (which needs to be addressed)...
               !
               DO MotherNuclide = 1,NuclideFamily%NFamily
                  iiNuclide = NuclideFamily%FamilyMember(MotherNuclide)
                  IF (Orphanage%Orphan(iOrphan)%Mother.EQ.iiNuclide) THEN ! So now the mother of the orphan is regular family

                     IF (.NOT.Orphanage%Orphan(iOrphan)%Active) THEN ! but the transition to a too fast nuclide in not yet active
                        TooFastFullyChecked(Orphanage%Orphan(iOrphan)%Daughter) = .FALSE.
                        ! Then the too fast daughter is not yet a good basis for searching for new too fast progeny below

                        IF (DebugLevel.GT.0) WRITE(*,'(A,I0,A)') 'MarkFullyCheckedOrphans: In iteration ',iIteration,&
                        & ' switched off orphan '&
                        & //TRIM(NuclideSpecs(Orphanage%Orphan(iOrphan)%Daughter)%NuclideName)&
                        & //' for lack of progeny transition '&
                        & //TRIM(NuclideSpecs(Orphanage%Orphan(iOrphan)%Mother)%NuclideName)//' --> '&
                        & //TRIM(NuclideSpecs(Orphanage%Orphan(iOrphan)%Daughter)%NuclideName)
                     ENDIF

                  ENDIF ! We found the mother of this orphan in the regular progeny of the head of chain
               ENDDO ! loop over mother nuclides
               !
               ! ... or a fully checked orphan that itself is in the big family
               !
               DO jOrphan = 1,Orphanage%NOrphans
                  iiNuclide = Orphanage%Orphan(jOrphan)%Daughter
                  IF (NuclideFamily%IsFarFamily(iiNuclide)) THEN
                     IF (Orphanage%Orphan(iOrphan)%Mother.EQ.iiNuclide) THEN ! So now the mother of the orphan is an orphan

                        IF (.NOT.Orphanage%Orphan(iOrphan)%Active) THEN ! but this mother orphan transition is not yet active..
                           TooFastFullyChecked(Orphanage%Orphan(iOrphan)%Daughter) = .FALSE. ! then, de-select this nuclide

                           IF (DebugLevel.GT.0) WRITE(*,'(A,I0,A)') 'MarkFullyCheckedOrphans: In iteration ',&
                           & iIteration,' switched off orphan '&
                           & //TRIM(NuclideSpecs(Orphanage%Orphan(iOrphan)%Daughter)%NuclideName)&
                           & //' for lack of orphan-orphan transition '&
                           & //TRIM(NuclideSpecs(Orphanage%Orphan(iOrphan)%Mother)%NuclideName)//' --> '&
                           & //TRIM(NuclideSpecs(Orphanage%Orphan(iOrphan)%Daughter)%NuclideName)
                        ENDIF

                     ENDIF ! We found the mother of this orphan in the regular progeny of the head of chain
                  ENDIF ! Check if the mother of this orphan belongs to this (big) family
               ENDDO ! loop over jOrphan
            ENDDO ! loop over iOrphan
         END SUBROUTINE MarkFullyCheckedOrphans

         SUBROUTINE AddOrphansFromRegularFamily(MyNAdded)
            INTEGER, INTENT(OUT) :: MyNAdded
            !
            ! Step 1: Loop over all regular nuclides in this family to see if they have a yet inactive orphan in their progeny
            !
            MyNAdded = 0

            DO MotherNuclide = 1,NuclideFamily%NFamily
               iiNuclide = NuclideFamily%FamilyMember(MotherNuclide)
               !
               ! Check all non-active orphans if they can originate from this mother
               !
               DO jOrphan = 1,Orphanage%NOrphans
                  !
                  ! You found a new family member if any new orphan is not already a member and if the orphan at hand is the mother of the new orphan
                  !
                  NotAlreadyAMember = (.NOT.Orphanage%Orphan(jOrphan)%Active)

                  ItsMotherIs_ii = (Orphanage%Orphan(jOrphan)%Mother.EQ.iiNuclide)

                  IF (NotAlreadyAMember.AND.ItsMotherIs_ii) THEN
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

                     IF (DebugLevel.GT.0) WRITE(*,'(A,F20.15)') 'AddOrphansFromRegularFamily: Found new orphan transition: '&
                     & //NuclideSpecs(jjNuclide)%NuclideName &
                     & //' from mother '&
                     & //NuclideSpecs(iiNuclide)%NuclideName&
                     & //' with yield ',Orphanage%Orphan(jOrphan)%Yield
                  ENDIF
               ENDDO ! loop over non-active orphans
            ENDDO ! loop over all regular nuclides in this family
         END SUBROUTINE AddOrphansFromRegularFamily

         SUBROUTINE AddOrphansFromCheckedOrphans(MyNAdded)
            INTEGER, INTENT(OUT) :: MyNAdded
            !
            ! Step 2: Loop over all fully checked orphan nuclides in this big family to see if they have a yet inactive orphan in their progeny
            !
            MyNAdded = 0

            DO iOrphan = 1,Orphanage%NOrphans
               iiNuclide = Orphanage%Orphan(iOrphan)%Daughter ! i.e. the endpoint of an orphan arrow, the too fast nuclide itself
               IF (NuclideFamily%IsFarFamily(iiNuclide)) THEN ! the too fast nuclide is far family!
                  IF (TooFastFullyChecked(iiNuclide)) THEN ! and all its incoming arrow have been accounted for, so not it can pass its activity to even lower decks!
                     !
                     ! Check all non-active orphans if they can originate from this (too fast) mother
                     !
                     DO jOrphan = 1,Orphanage%NOrphans
                        !
                        ! You found a new family member if any new orphan is not already a member and if the orphan at hand is the mother of the new orphan
                        !
                        NotAlreadyAMember = (.NOT.Orphanage%Orphan(jOrphan)%Active)

                        ItsMotherIs_ii = (Orphanage%Orphan(jOrphan)%Mother.EQ.iiNuclide)

                        IF (NotAlreadyAMember.AND.ItsMotherIs_ii) THEN
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

                           IF (DebugLevel.GT.0) WRITE(*,'(A)') 'AddOrphansFromCheckedOrphans: Found new orphan transition: '&
                           & //NuclideSpecs(jjNuclide)%NuclideName &
                           & //' from mother '&
                           & //NuclideSpecs(iiNuclide)%NuclideName
                        ENDIF
                     ENDDO ! loop over non-active orphans
                  ENDIF ! orphan active
               ENDIF ! possible mother is an orphan herself that belongs to the big family
            ENDDO ! loop over all active orphans
         END SUBROUTINE AddOrphansFromCheckedOrphans
   END SUBROUTINE AddOrphansBelow
END MODULE LibENDF
