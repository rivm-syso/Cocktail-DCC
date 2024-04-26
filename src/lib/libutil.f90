MODULE LibUtil
   !
   ! General purpose utility routines
   !
   ! ____________________________________________________
   !
   ! Developed For:
   ! RIVM - National Institute for Public Health and the Environment
   !
   ! PO Box 1,
   ! NL - 3720 BA Bilthoven
   ! The Netherlands
   ! ____________________________________________________
   !
   USE libxmath
   IMPLICIT NONE

   PRIVATE

   PUBLIC :: RunShellCommand,FileExists,ListModules,PathAndName,MSYSDir,&
      & UpCase,LowCase,Capitalize,DeCapitalize,CutStringInTwo,&
      & AllUpCase,AllLowCase,DefaultLength,ScratchFile,ScrotchFile,&
      & Male,Female,GenderName,RemoveCharacter,CharacterIsADigit
   !
   ! The following file-units are reserved for scratch/quickservice
   ! etcetera. Only use it if you are SURE that you will close
   ! it immediatelly!
   !
   INTEGER, PARAMETER :: ScratchFile = 17
   INTEGER, PARAMETER :: ScrotchFile = 18
   !
   ! Gender
   !
   INTEGER, PARAMETER :: Male = 1
   INTEGER, PARAMETER :: Female = 2

   CHARACTER(6), DIMENSION(2), PARAMETER :: GenderName = ((/'male  ','female'/))
   !
   ! To accommodate for long file- and directory names, by default such strings
   ! have the following length:
   !
   INTEGER, PARAMETER :: DefaultLength = 512 ! With 256 you get problems when wgrib is asked to take action...
   !
   ! If you call these sources from MinGW/MSYS, then SYSTEM-calls require the path to MSYS.
   ! The default path is set here. With a setting in a settings-file the default can be overruled.
   !
   CHARACTER(DefaultLength) :: MSYSDir = 'c:/MinGW/msys/1.0/bin'
   !!
   ! There are always problems with directory/filename combinations:
   ! Do you have to end the path with a (back-)slash or not?
   ! To solve this problem, a new datatype is created, plus a routine
   ! that glues them together.
   TYPE PathNameType
      CHARACTER(DefaultLength) :: Path,Name
   END TYPE PathNameType
   !
   ! To discriminate between forward-slash and backslash, you can refer to the following:
   !
   INTEGER, PARAMETER :: System_Linux = 1
   INTEGER, PARAMETER :: System_Windows = 2
   INTEGER :: TheSystem = System_Linux
   !
   ! The following parameter can be set for debugging this MODULE
   !
   INTEGER, PARAMETER :: DebugLevel = 0

CONTAINS
   FUNCTION CharacterIsADigit(Ch)
      CHARACTER(1), INTENT(IN) :: Ch
      LOGICAL :: CharacterIsADigit
      CharacterIsADigit = (INDEX('0123456789',Ch).NE.0)
   END FUNCTION CharacterIsADigit



   FUNCTION PathAndName(PathName)
      !
      ! Return a single string resulting from glueing together directory and filename.
      ! Take a very close look at the slash between them!
      !
      TYPE(PathNameType), INTENT(IN) :: PathName

      INTEGER :: TheEnd
      CHARACTER(DefaultLength) :: PathAndName

      TheEnd = LEN_TRIM(PathName%Path)
      IF (PathName%Path(TheEnd:TheEnd).EQ.'/') THEN
         PathAndName = TRIM(PathName%Path)//TRIM(PathName%Name)
      ELSE
         PathAndName = TRIM(PathName%Path)//'/'//TRIM(PathName%Name)
      ENDIF
   END FUNCTION PathAndName



   SUBROUTINE RunShellCommand(Commando,CrashOnError,DoSilent)
      !
      ! Perform an external shell-command and if necessary halt on error
      !
      CHARACTER(*), INTENT(IN) :: Commando
      INTEGER :: Error,SYSTEM !,ThePlace
      LOGICAL, INTENT(IN) :: CrashOnError,DoSilent
      CHARACTER(DefaultLength) :: ReversedCommand

      INTEGER, PARAMETER :: DebugLevel = 0
      !
      ! Make a copy
      !
      ReversedCommand = TRIM(MSYSDir)//'/'//TRIM(Commando)
      !
      ! Try the command
      !
      IF (TheSystem.EQ.System_Windows) THEN
         IF (.NOT.DoSilent) WRITE(*,10) TRIM(ReversedCommand)
         Error = SYSTEM(TRIM(ReversedCommand))
      ELSE
         IF (.NOT.DoSilent) WRITE(*,10) TRIM(Commando)
         Error = SYSTEM(TRIM(Commando))
      ENDIF
      10 FORMAT('RunShellCommand : "',A,'"')
      !
      ! If failure, try with reversed slashes
      !
      IF (Error.NE.0) THEN
         IF (.NOT.DoSilent) WRITE(*,'(A)') 'Failure in shell-command. Retrying with path c:/msys/1.0/bin/ prepended...'

         IF (TheSystem.EQ.System_Windows) THEN
            IF (.NOT.DoSilent) WRITE(*,10) TRIM(Commando)
            Error = SYSTEM(TRIM(Commando))
            IF (Error.EQ.0) THEN
               TheSystem = System_Linux ! Switch system mode to Linux
               WRITE(*,*)
               WRITE(*,'(A)') '--> System reacts better to original way of calling: changing to Linux-mode!'
               WRITE(*,*)
            ENDIF
         ELSE
            IF (.NOT.DoSilent) WRITE(*,10) TRIM(ReversedCommand)
            Error = SYSTEM(TRIM(ReversedCommand))
            IF (Error.EQ.0) THEN
               TheSystem = System_Windows ! Switch system mode to Windows
               WRITE(*,*)
               WRITE(*,'(A)') '--> System reacts better to MSYS-specific calling; changing to Windows-mode!'
               WRITE(*,*)
            ENDIF
         ENDIF

         IF (Error.NE.0) THEN
            IF (.NOT.DoSilent) WRITE(*,*) '... sorry, second attempt of shell-command failed either!'
            IF (.NOT.DoSilent) WRITE(*,20) Error
            20    FORMAT('Error ',I7,' in execution of command')
            IF (CrashOnError) STOP
         ENDIF
      ENDIF

      IF ((Error.EQ.0).AND.(.NOT.DoSilent)) WRITE(*,30)
      30  FORMAT('Shell-command succeeded!')

   END SUBROUTINE RunShellCommand



   FUNCTION FileExists(InName)
      !
      ! Check if a certain file exists
      !
      LOGICAL ThisFileExists,FileExists
      CHARACTER(*), INTENT(IN) :: InName
      INQUIRE(FILE=TRIM(InName),EXIST=ThisFileExists)
      FileExists = ThisFileExists
   END FUNCTION FileExists



   FUNCTION PathNameExists(PathName)
      !
      ! Check if a certain file in a given directory exists
      !
      LOGICAL ThisPathNameExists,PathNameExists
      TYPE(PathNameType), INTENT(IN) :: PathName
      INQUIRE(FILE=TRIM(PathAndName(PathName)),EXIST=ThisPathNameExists)
      PathNameExists = ThisPathNameExists
   END FUNCTION PathNameExists



   SUBROUTINE RemoveCharacter(MyString,MyCharacter)
      !
      ! Remove 1 character completely from a string
      !
      CHARACTER(*), INTENT(INOUT) :: MyString
      CHARACTER, INTENT(IN) :: MyCharacter
      INTEGER :: MyCharacterPos

      MyCharacterPos = INDEX(MyString,MyCharacter)
      DO WHILE ((MyCharacterPos.GT.0).AND.(MyCharacterPos.LE.LEN_TRIM(MyString)))
         MyString = MyString(1:(MyCharacterPos-1))//MyString((MyCharacterPos+1):LEN_TRIM(MyString))
         MyCharacterPos = INDEX(MyString,MyCharacter)
      ENDDO
   END SUBROUTINE RemoveCharacter



   SUBROUTINE CutStringInTwo(StringIn,String1,String2,NoSecondString)
      !
      ! Cut a string in two pieces. The first part is the first non-trivial
      ! substring until either a space or a tab character
      !
      CHARACTER(*), INTENT(INOUT) :: StringIn
      CHARACTER(*), INTENT(OUT) :: String1,String2
      LOGICAL, INTENT(OUT) :: NoSecondString
      !
      ! Local variables
      !
      INTEGER :: SpacePosition,TabPosition,CutPosition
      CHARACTER(500) :: DumStr1,DumStr2
      !
      ! Remove leading spaces or tabs
      !
      IF (DebugLevel.GT.0) WRITE(*,10) '"'//TRIM(StringIn)//'"'
      10 FORMAT('Cut : StringIn = ',A)
      DumStr1 = StringIn
      DO WHILE (((DumStr1(1:1).EQ.CHAR(32)).OR.(DumStr1(1:1).EQ.CHAR(9))).AND.(LEN_TRIM(DumStr1).GT.1))
         DumStr1 = DumStr1(2:LEN_TRIM(DumStr1))
      ENDDO
      IF (DebugLevel.GT.0) WRITE(*,11) TRIM(DumStr1)
      11 FORMAT(' Cut : DumStr1 zonder spaties aan het begin = "',A,'"')
      !
      ! Remove trailing tabs
      !
      DO WHILE ((DumStr1(LEN_TRIM(DumStr1):LEN_TRIM(DumStr1)).EQ.CHAR(9)).AND.(LEN_TRIM(DumStr1).GT.1))
         DumStr1 = DumStr1(1:(LEN_TRIM(DumStr1)-1))
      ENDDO
      IF (DebugLevel.GT.0) WRITE(*,12) TRIM(DumStr1)
      12 FORMAT(' Cut : DumStr1 zonder spaties aan het eind = "',A,'"')
      !
      ! Find a separator: space or tab
      !
      SpacePosition = INDEX(DumStr1,CHAR(32))
      TabPosition = INDEX(DumStr1,CHAR(9))
      CutPosition = SpacePosition
      IF ((TabPosition.GT.0).AND.(TabPosition.LT.SpacePosition)) CutPosition = TabPosition

      NoSecondString = ((CutPosition.EQ.0).OR.(CutPosition.EQ.(LEN_TRIM(DumStr1)+1)))

      IF (DebugLevel.GT.0) WRITE(*,20) TRIM(DumStr1),SpacePosition,TabPosition,CutPosition
      20 FORMAT(' "',A,'" heeft spatie op ',I3,' en tab op ',I3,' dus we gaan knippen op ',I3)
      IF (.NOT.NoSecondString) THEN
         DumStr2 = DumStr1((CutPosition+1):LEN_TRIM(DumStr1))
         DumStr1 = DumStr1(1:(CutPosition-1))
      ENDIF
      IF (DebugLevel.GT.0) THEN
         IF (NoSecondString) THEN
            WRITE(*,31) TRIM(DumStr1)
            31 FORMAT(' Aan het eind: "',A,'"  en verder niets!')
         ELSE
            WRITE(*,30) TRIM(DumStr1),TRIM(DumStr2)
            30 FORMAT(' Aan het eind: "',A,'"  en  "',A,'"')
         ENDIF
      ENDIF
      String1 = TRIM(DumStr1)
      IF (.NOT.NoSecondString) String2 = TRIM(DumStr2)
      IF (DebugLevel.GT.0) WRITE(*,*)
   END SUBROUTINE CutStringInTwo



   CHARACTER FUNCTION UpCase(Ch)
      !
      ! Make a letter upper case
      !
      CHARACTER, INTENT(IN) :: Ch
      CHARACTER :: DumChar
      IF ((ICHAR(Ch).GE.97).AND.(ICHAR(Ch).LE.122)) THEN
         DumChar = CHAR(ICHAR(Ch)-32)
      ELSE
         DumChar = Ch
      ENDIF
      UpCase = DumChar
   END FUNCTION UpCase



   CHARACTER FUNCTION LowCase(Ch)
      !
      ! Make a letter lower case
      !
      CHARACTER, INTENT(IN) :: Ch
      CHARACTER :: DumChar
      IF ((ICHAR(Ch).GE.65).AND.(ICHAR(Ch).LE.90)) THEN
         DumChar = CHAR(ICHAR(Ch)+32)
      ELSE
         DumChar = Ch
      ENDIF
      LowCase = DumChar
   END FUNCTION LowCase



   SUBROUTINE AllUpCase(String)
      !
      ! Set all letters to capital
      !
      CHARACTER(*), INTENT(INOUT) :: String
      INTEGER :: i,StringLength
      StringLength = LEN_TRIM(String)
      DO i=1,StringLength
         String(i:i) = UpCase(String(i:i))
      ENDDO
   END SUBROUTINE AllUpCase



   SUBROUTINE AllLowCase(String)
      !
      ! Set all letters to lower case
      !
      CHARACTER(*), INTENT(INOUT) :: String
      INTEGER :: i,StringLength
      StringLength = LEN_TRIM(String)
      DO i=1,StringLength
         String(i:i) = LowCase(String(i:i))
      ENDDO
   END SUBROUTINE AllLowCase



   SUBROUTINE Capitalize(String)
      !
      ! Set first letter to capital, rest to lower case
      !
      CHARACTER(*), INTENT(INOUT) :: String
      INTEGER :: i,StringLength
      StringLength = LEN_TRIM(String)
      String(1:1) = UpCase(String(1:1))
      DO i=2,StringLength
         String(i:i) = LowCase(String(i:i))
      ENDDO
   END SUBROUTINE Capitalize



   SUBROUTINE DeCapitalize(String)
      !
      ! Set all letters to lower case
      !
      CHARACTER(*), INTENT(INOUT) :: String
      INTEGER :: i,StringLength
      StringLength = LEN_TRIM(String)
      DO i=1,StringLength
         String(i:i) = LowCase(String(i:i))
      ENDDO
   END SUBROUTINE DeCapitalize



   FUNCTION ListModules(SourceFileName)
      !
      ! Make a list of the modulenames used in a sourcefile
      !
      CHARACTER(*), INTENT(IN) :: SourceFileName
      CHARACTER(DefaultLength) :: Commando
      LOGICAL, PARAMETER :: CrashOnError = .FALSE., DoSilent = .TRUE.
      LOGICAL :: ListModules,Okay
      !
      ! Make a list of the MODULES used by the sourcefile
      !
      Okay = FileExists(SourceFileName)
      IF (Okay) THEN
         Commando = 'grep -i -w use '//TRIM(SourceFileName)//' > listmodules.tmp'
         CALL RunShellCommand(Commando,CrashOnError,DoSilent)
      ENDIF

      ListModules = Okay
   END FUNCTION ListModules
END MODULE LibUtil
