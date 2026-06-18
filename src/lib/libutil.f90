module LibUtil
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
   implicit none

   private
   public :: env_var
   public :: RunShellCommand,FileExists,ListModules,PathAndName,MSYSDir,&
      & UpCase,LowCase,Capitalize,DeCapitalize,CutStringInTwo,&
      & AllUpCase,AllLowCase,DefaultLength,ScratchFile,ScrotchFile,&
      & Male,Female,GenderName,RemoveCharacter,CharacterIsADigit
   !
   ! The following file-units are reserved for scratch/quickservice
   ! etcetera. Only use it if you are SURE that you will close
   ! it immediatelly!
   !
   integer, parameter :: ScratchFile = 17
   integer, parameter :: ScrotchFile = 18
   !
   ! Gender
   !
   integer, parameter :: Male = 1
   integer, parameter :: Female = 2

   character(6), dimension(2), parameter :: GenderName = (['male  ','female'])
   !
   ! To accommodate for long file- and directory names, by default such strings
   ! have the following length:
   !
   integer, parameter :: DefaultLength = 512 ! With 256 you get problems when wgrib is asked to take action...
   !
   ! If you call these sources from MinGW/MSYS, then SYSTEM-calls require the path to MSYS.
   ! The default path is set here. With a setting in a settings-file the default can be overruled.
   !
   character(DefaultLength) :: MSYSDir = 'c:/MinGW/msys/1.0/bin'
   !!
   ! There are always problems with directory/filename combinations:
   ! Do you have to end the path with a (back-)slash or not?
   ! To solve this problem, a new datatype is created, plus a routine
   ! that glues them together.
   type PathNameType
      character(DefaultLength) :: Path,Name
   end type PathNameType
   !
   ! To discriminate between forward-slash and backslash, you can refer to the following:
   !
   integer, parameter :: System_Linux = 1
   integer, parameter :: System_Windows = 2
   integer :: TheSystem = System_Linux
   !
   ! The following parameter can be set for debugging this module
   !
   integer, parameter :: DebugLevel = 0

contains
   subroutine env_var(name, val)
      character(*), intent(in) :: name
      character(:), allocatable, intent(out) :: val

      integer :: n
      integer :: status
      logical :: is_in_env

      call get_environment_variable(name, length=n, status=status)

      select case(status)
       case(1);      is_in_env = .false.  ! environment variable has not been set
       case(2);      is_in_env = .false.  ! machine does not support environment variables
       case(0);      is_in_env = .true.
       case default; error stop 'unknown STATUS from GET_ENVIRONMENT_VARIABLE'
      end select

      if (is_in_env) then
         allocate(character(n) :: val)
         call get_environment_variable(name, value=val, status=status)
         if (status /= 0) error stop 'inconsistent STATUS from intrinsic GET_ENVIRONMENT_VARIABLE'
      endif
   end subroutine

   function CharacterIsADigit(Ch)
      character(1), intent(in) :: Ch
      logical :: CharacterIsADigit
      CharacterIsADigit = (INDEX('0123456789',Ch) /= 0)
   end function CharacterIsADigit



   function PathAndName(PathName)
      !
      ! Return a single string resulting from glueing together directory and filename.
      ! Take a very close look at the slash between them!
      !
      type(PathNameType), intent(in) :: PathName

      integer :: TheEnd
      character(DefaultLength) :: PathAndName

      TheEnd = len_trim(PathName%Path)
      if (PathName%Path(TheEnd:TheEnd) == '/') then
         PathAndName = trim(PathName%Path)//trim(PathName%Name)
      else
         PathAndName = trim(PathName%Path)//'/'//trim(PathName%Name)
      endif
   end function PathAndName



   subroutine RunShellCommand(Commando,CrashOnError,DoSilent)
      !
      ! Perform an external shell-command and if necessary halt on error
      !
      character(*), intent(in) :: Commando
      integer :: Error,SYSTEM !,ThePlace
      logical, intent(in) :: CrashOnError,DoSilent
      character(DefaultLength) :: ReversedCommand
      !
      ! Make a copy
      !
      ReversedCommand = trim(MSYSDir)//'/'//trim(Commando)
      !
      ! Try the command
      !
      if (TheSystem == System_Windows) then
         if (.NOT.DoSilent) write(*,10) trim(ReversedCommand)
         Error = SYSTEM(trim(ReversedCommand))
      else
         if (.NOT.DoSilent) write(*,10) trim(Commando)
         Error = SYSTEM(trim(Commando))
      endif
      10 format('RunShellCommand : "',A,'"')
      !
      ! If failure, try with reversed slashes
      !
      if (Error /= 0) then
         if (.NOT.DoSilent) write(*,'(A)') 'Failure in shell-command. Retrying with path c:/msys/1.0/bin/ prepended...'

         if (TheSystem == System_Windows) then
            if (.NOT.DoSilent) write(*,10) trim(Commando)
            Error = SYSTEM(trim(Commando))
            if (Error == 0) then
               TheSystem = System_Linux ! Switch system mode to Linux
               write(*,*)
               write(*,'(A)') '--> System reacts better to original way of calling: changing to Linux-mode!'
               write(*,*)
            endif
         else
            if (.NOT.DoSilent) write(*,10) trim(ReversedCommand)
            Error = SYSTEM(trim(ReversedCommand))
            if (Error == 0) then
               TheSystem = System_Windows ! Switch system mode to Windows
               write(*,*)
               write(*,'(A)') '--> System reacts better to MSYS-specific calling; changing to Windows-mode!'
               write(*,*)
            endif
         endif

         if (Error /= 0) then
            if (.NOT.DoSilent) write(*,*) '... sorry, second attempt of shell-command failed either!'
            if (.NOT.DoSilent) write(*,20) Error
            20    format('Error ',I7,' in execution of command')
            if (CrashOnError) stop
         endif
      endif

      if ((Error == 0).AND.(.NOT.DoSilent)) write(*,30)
      30  format('Shell-command succeeded!')

   end subroutine RunShellCommand



   function FileExists(InName)
      !
      ! Check if a certain file exists
      !
      logical ThisFileExists,FileExists
      character(*), intent(in) :: InName
      inquire(file=trim(InName),exist=ThisFileExists)
      FileExists = ThisFileExists
   end function FileExists



   function PathNameExists(PathName)
      !
      ! Check if a certain file in a given directory exists
      !
      logical ThisPathNameExists,PathNameExists
      type(PathNameType), intent(in) :: PathName
      inquire(file=trim(PathAndName(PathName)),exist=ThisPathNameExists)
      PathNameExists = ThisPathNameExists
   end function PathNameExists



   subroutine RemoveCharacter(MyString,MyCharacter)
      !
      ! Remove 1 character completely from a string
      !
      character(*), intent(inout) :: MyString
      character, intent(in) :: MyCharacter
      integer :: MyCharacterPos

      MyCharacterPos = INDEX(MyString,MyCharacter)
      do while ((MyCharacterPos > 0).AND.(MyCharacterPos <= len_trim(MyString)))
         MyString = MyString(1:(MyCharacterPos-1))//MyString((MyCharacterPos+1):len_trim(MyString))
         MyCharacterPos = INDEX(MyString,MyCharacter)
      enddo
   end subroutine RemoveCharacter



   subroutine CutStringInTwo(StringIn,String1,String2,NoSecondString)
      !
      ! Cut a string in two pieces. The first part is the first non-trivial
      ! substring until either a space or a tab character
      !
      character(*), intent(inout) :: StringIn
      character(*), intent(out) :: String1,String2
      logical, intent(out) :: NoSecondString
      !
      ! Local variables
      !
      integer :: SpacePosition,TabPosition,CutPosition
      character(500) :: DumStr1,DumStr2
      !
      ! Remove leading spaces or tabs
      !
      if (DebugLevel > 0) write(*,10) '"'//trim(StringIn)//'"'
      10 format('Cut : StringIn = ',A)
      DumStr1 = StringIn
      do while (((DumStr1(1:1) == CHAR(32)).OR.(DumStr1(1:1) == CHAR(9))).AND.(len_trim(DumStr1) > 1))
         DumStr1 = DumStr1(2:len_trim(DumStr1))
      enddo
      if (DebugLevel > 0) write(*,11) trim(DumStr1)
      11 format(' Cut : DumStr1 zonder spaties aan het begin = "',A,'"')
      !
      ! Remove trailing tabs
      !
      do while ((DumStr1(len_trim(DumStr1):len_trim(DumStr1)) == CHAR(9)).AND.(len_trim(DumStr1) > 1))
         DumStr1 = DumStr1(1:(len_trim(DumStr1)-1))
      enddo
      if (DebugLevel > 0) write(*,12) trim(DumStr1)
      12 format(' Cut : DumStr1 zonder spaties aan het eind = "',A,'"')
      !
      ! Find a separator: space or tab
      !
      SpacePosition = INDEX(DumStr1,CHAR(32))
      TabPosition = INDEX(DumStr1,CHAR(9))
      CutPosition = SpacePosition
      if ((TabPosition > 0).AND.(TabPosition < SpacePosition)) CutPosition = TabPosition

      NoSecondString = ((CutPosition == 0).OR.(CutPosition == (len_trim(DumStr1)+1)))

      if (DebugLevel > 0) write(*,20) trim(DumStr1),SpacePosition,TabPosition,CutPosition
      20 format(' "',A,'" heeft spatie op ',I3,' en tab op ',I3,' dus we gaan knippen op ',I3)
      if (.NOT.NoSecondString) then
         DumStr2 = DumStr1((CutPosition+1):len_trim(DumStr1))
         DumStr1 = DumStr1(1:(CutPosition-1))
      endif
      if (DebugLevel > 0) then
         if (NoSecondString) then
            write(*,31) trim(DumStr1)
            31 format(' Aan het eind: "',A,'"  en verder niets!')
         else
            write(*,30) trim(DumStr1),trim(DumStr2)
            30 format(' Aan het eind: "',A,'"  en  "',A,'"')
         endif
      endif
      String1 = trim(DumStr1)
      if (.NOT.NoSecondString) String2 = trim(DumStr2)
      if (DebugLevel > 0) write(*,*)
   end subroutine CutStringInTwo



   character function UpCase(Ch)
      !
      ! Make a letter upper case
      !
      character, intent(in) :: Ch
      character :: DumChar
      if ((ICHAR(Ch) >= 97).AND.(ICHAR(Ch) <= 122)) then
         DumChar = CHAR(ICHAR(Ch)-32)
      else
         DumChar = Ch
      endif
      UpCase = DumChar
   end function UpCase



   character function LowCase(Ch)
      !
      ! Make a letter lower case
      !
      character, intent(in) :: Ch
      character :: DumChar
      if ((ICHAR(Ch) >= 65).AND.(ICHAR(Ch) <= 90)) then
         DumChar = CHAR(ICHAR(Ch)+32)
      else
         DumChar = Ch
      endif
      LowCase = DumChar
   end function LowCase



   subroutine AllUpCase(String)
      !
      ! Set all letters to capital
      !
      character(*), intent(inout) :: String
      integer :: i,StringLength
      StringLength = len_trim(String)
      do i=1,StringLength
         String(i:i) = UpCase(String(i:i))
      enddo
   end subroutine AllUpCase



   subroutine AllLowCase(String)
      !
      ! Set all letters to lower case
      !
      character(*), intent(inout) :: String
      integer :: i,StringLength
      StringLength = len_trim(String)
      do i=1,StringLength
         String(i:i) = LowCase(String(i:i))
      enddo
   end subroutine AllLowCase



   subroutine Capitalize(String)
      !
      ! Set first letter to capital, rest to lower case
      !
      character(*), intent(inout) :: String
      integer :: i,StringLength
      StringLength = len_trim(String)
      String(1:1) = UpCase(String(1:1))
      do i=2,StringLength
         String(i:i) = LowCase(String(i:i))
      enddo
   end subroutine Capitalize



   subroutine DeCapitalize(String)
      !
      ! Set all letters to lower case
      !
      character(*), intent(inout) :: String
      integer :: i,StringLength
      StringLength = len_trim(String)
      do i=1,StringLength
         String(i:i) = LowCase(String(i:i))
      enddo
   end subroutine DeCapitalize



   function ListModules(SourceFileName)
      !
      ! Make a list of the modulenames used in a sourcefile
      !
      character(*), intent(in) :: SourceFileName
      character(DefaultLength) :: Commando
      logical, parameter :: CrashOnError = .FALSE., DoSilent = .TRUE.
      logical :: ListModules,Okay
      !
      ! Make a list of the MODULES used by the sourcefile
      !
      Okay = FileExists(SourceFileName)
      if (Okay) then
         Commando = 'grep -i -w use '//trim(SourceFileName)//' > listmodules.tmp'
         call RunShellCommand(Commando,CrashOnError,DoSilent)
      endif

      ListModules = Okay
   end function ListModules
end module LibUtil
