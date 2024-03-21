PROGRAM GetDeps
   !
   ! Make a list of all sources necessary to compile a given F90 input-file.
   ! All USE statements are checked and recursively followed
   ! until no more USE statements are found.
   ! The source-files are placed in order such that no sourcefile
   ! uses modules that have not yet been specified.
   !
   ! Convention:
   !
   ! - All modules must be referred to by a USE statement, where the word USE is in capital letters
   ! - All modules are stored in files having the same name as the name of the module, but then with extension .f90.
   ! - Module-names are case sensitive
   !
   ! ____________________________________________________
   !
   ! Version: 4 July 2022
   !
   ! Developed For:
   !
   ! RIVM - National Institute for Public Health and the Environment
   !
   ! PO Box 1,
   ! NL - 3720 BA Bilthoven
   ! The Netherlands
   ! Website: http://www.rivm.nl/
   !
   ! ____________________________________________________
   !
   USE libutil

   IMPLICIT NONE

   INTEGER, PARAMETER :: MaxNSources = 100

   INTEGER :: IARGC,NSources,NArguements,iLow,iHigh,iSource,jSource,&
   & iLine,NewiLow,NewiHigh,SourceCounter,PrevSourceCounter,i,NCharacters,iIntrinsicModule
   INTEGER, DIMENSION(MaxNSources) :: iMother
   CHARACTER(DefaultLength), DIMENSION(MaxNSources) :: FName
   CHARACTER(DefaultLength) :: ThisLine,FirstString,SecondString,ThirdString,&
   & FormatString,DumName,Commando
   CHARACTER(10000) :: OutLine
   LOGICAL, DIMENSION(MaxNSources,MaxNSources) :: UsesModule,Unresolved
   LOGICAL :: Ready,IsOldSource,NoSecondString,NoThirdString,ItIsThisFile,ValidLine,DaughterExists
   INTEGER, DIMENSION(MaxNSources) :: SourceOrder,NUnresolved
   LOGICAL, DIMENSION(MaxNSources) :: Resolved
   LOGICAL, PARAMETER :: CrashOnError = .TRUE.
   LOGICAL, PARAMETER :: DoSilent = .TRUE.

   INTEGER, PARAMETER :: NIntrinsicModules = 1
   CHARACTER(20), DIMENSION(NIntrinsicModules) :: IntrinsicModuleName = &
   & ('iso_c_binding.f90   ')

   INTEGER :: DebugLevel
   !
   ! Parse commandline parameter: which file to process?
   !
   NArguements = IARGC()
   !
   ! Debugmode is activated by giving a third commandline parameter
   !
   IF (NArguements.LE.1) THEN
      DebugLevel = 0
   ELSE
      DebugLevel = 1
   ENDIF
   !
   ! Give info on utility if no parameter is passed
   !
   IF ((NArguements.EQ.0).OR.(DebugLevel.GT.0)) THEN
      WRITE(*,*)
      WRITE(*,10)
      10  FORMAT('call: getdeps <f90sourcefilename> [<debuglevel>]')
      WRITE(*,*)
      WRITE(*,20)
      20  FORMAT('creates a string of all sourcefiles with all recursive module dependencies in order of dependence')
      WRITE(*,21)
      21  FORMAT('Second parameter <debuglevel> is optional.')
      WRITE(*,22)
      22  FORMAT('Giving any second parameter generates lots of intermediate output.')
   ENDIF

   IF (NArguements.GT.0) THEN
      !
      ! Read name of file to process
      !
      CALL GETARG(1,FName(1))
      IF (DebugLevel.GT.0) WRITE(*,50) TRIM(FName(1))
      50  FORMAT('Basic sourcefile =                                         "',A,'"')
      iMother(1) = 0
      !
      ! Recursive loop over sources
      !
      NSources = 1 ! Until now only the sourcefile specified on the commandline...
      UsesModule = .FALSE. ! No dependencies yet...
      iLow = 1
      iHigh = 1

      Ready = .FALSE.

      DO WHILE (.NOT.Ready)
         IF (DebugLevel.GT.0) WRITE(*,*) '###########################################################################'
         IF (DebugLevel.GT.0) WRITE(*,*) 'New iteration started for finding new sourcefiles that are used'
         IF (DebugLevel.GT.0) WRITE(*,*) 'Loop over file ',iLow,' to ',iHigh
         !
         ! Already keep track of the next iterations boundaries
         !
         NewiLow = iHigh+1
         NewiHigh = iHigh
         !
         ! Loop over sources found in last iteration (1 file if just started)
         !
         DO iSource = iLow, iHigh
            IF (DebugLevel.GT.0) THEN
               WRITE(*,*)
               WRITE(*,60) iSource,TRIM(FName(iSource))
               60        FORMAT('Checking dependencies for source file(',I0,'): "',A,'"')
            ENDIF
            !
            ! Convert to unix-style ASCII, just to make sure...
            !
            !         Commando = 'dos2unix '//TRIM(FName(iSource))
            !         CALL RunShellCommand(Commando,CrashOnError,DoSilent)
            !
            ! Dump all USE statements in scratchfile named 'listmodules.tmp'
            !
            DaughterExists = ListModules(FName(iSource))
            IF (.NOT.DaughterExists) THEN
               IF (iMother(iSource).GT.0) THEN
                  WRITE(*,'(A)') 'Unresolved dependence for '//TRIM(FName(iSource))//&
                  & ' in sourcefile '//TRIM(FName(iMother(iSource)))
               ELSE
                  WRITE(*,'(A)') 'File '//TRIM(FName(1))//' does not exist!'
               ENDIF
            !
            ! Check if new sourcefiles are encountered and update dependence matrix
            !
            ELSE IF (FileExists('listmodules.tmp')) THEN
               OPEN(ScratchFile,FILE='listmodules.tmp',FORM='FORMATTED',&
                  & POSITION='REWIND')
               iLine = 0
               DO
                  READ(ScratchFile,'(A255)',END=911) ThisLine
                  iLine = iLine + 1
                  IF (DebugLevel.GT.0) WRITE(*,80) iLine,TRIM(ThisLine)
                  80 FORMAT('     Line(',I3,') = "',A,'"')
                  !
                  ! Check if this is a valid line, which means that "USE" should be the first
                  ! text on this line other than spaces and tabs
                  !
                  CALL DeCapitalize(ThisLine)
                  ValidLine = .TRUE.
                  NCharacters = INDEX(ThisLine,'use')-1
                  DO i=1,NCharacters
                     ValidLine = (ValidLine.AND.&
                     &             (&
                     &               (ThisLine(i:i).EQ.CHAR(32)).OR.(ThisLine(i:i).EQ.CHAR(9))&
                     &              )&
                     &            )
                  ENDDO

                  IF (.NOT.ValidLine) THEN
                     IF (DebugLevel.GT.0) WRITE(*,'(A)') 'Invalid line, proceeding to next...'
                  ELSE
                     !
                     ! Distil name of sourcefile from this line
                     !
                     CALL CutStringInTwo(ThisLine,FirstString,SecondString,NoSecondString)
                     CALL CutStringInTwo(SecondString,DumName,ThirdString,NoThirdString)
                     DumName = TRIM(DumName)//'.f90'
                     IF (DebugLevel.GT.0) WRITE(*,70) TRIM(DumName)
                     70 FORMAT('     Dependence found on = "',A,'"')
                     !
                     ! Did we already have this sourcefile?
                     ! One possibility of an "old" source file is an intrinsic module (for which there is no source file...)
                     !
                     IsOldSource = .FALSE.

                     DO iIntrinsicModule = 1,NIntrinsicModules
                        ItIsThisFile = (TRIM(DumName).EQ.(TRIM(IntrinsicModuleName(iIntrinsicModule))))
                        IF (DebugLevel.GT.1) WRITE(*,*) 'Suggested module:"'//TRIM(DumName)//'"'
                        IF (DebugLevel.GT.1) WRITE(*,*) 'Internal module: "'//TRIM(IntrinsicModuleName(iIntrinsicModule))
                        IF (DebugLevel.GT.1) WRITE(*,*) 'Same? ',ItIsThisFile
                        IsOldSource = (IsOldSource.OR.ItIsThisFile)
                        IF ((DebugLevel.GT.0).AND.ItIsThisFile) WRITE(*,71) TRIM(IntrinsicModuleName(iIntrinsicModule))
                        71 FORMAT('This is an internal module: ',A,', so I will NOT try to read its sourcefile!')
                     ENDDO
                     !
                     ! Continue checking if the module is already known
                     !
                     IF (.NOT.IsOldSource) THEN
                        DO SourceCounter = 1,NSources
                           ItIsThisFile = (TRIM(DumName).EQ.TRIM(FName(SourceCounter)))
                           IsOldSource = (IsOldSource.OR.ItIsThisFile)
                           IF (ItIsThisFile) jSource = SourceCounter
                        ENDDO
                     ENDIF
                     IF (DebugLevel.GT.0) WRITE(*,*) '     OldFile? = ',IsOldSource
                     !
                     ! If this is a new sourcefile, update the filelist and the new loop boundaries
                     !
                     IF (.NOT.IsOldSource) THEN
                        NSources = NSources + 1
                        FName(NSources) = TRIM(DumName)
                        iMother(NSources) = iSource
                        IF (DebugLevel.GT.0) WRITE(*,90) TRIM(FName(NSources))
                        90 FORMAT('     Adding new sourcefile:                                "',A,'"')
                        NewiHigh = NewiHigh + 1
                        jSource = NSources
                     ENDIF
                     IF (DebugLevel.GT.0) WRITE(*,*) '     Dependent file has number = ',jSource
                     !
                     ! Update dependence matrix
                     !
                     UsesModule(iSource,jSource) = .TRUE.
                  ENDIF ! valid line
               ENDDO ! loop over contents of file 'listmodules.tmp', with dependencies for 1 daughter sourcefile
               911   CLOSE(ScratchFile)
            ENDIF ! found file listmodules.tmp
         ENDDO ! loop over all new sourcefiles found in the previous iteration
         !
         ! Check for convergence of iteration: no more USE statements found
         !
         IF (DebugLevel.GT.0) WRITE(*,*) 'NewiHigh,iHigh = ',NewiHigh,iHigh
         Ready = (NewiHigh.EQ.iHigh)
         !
         iLow = NewiLow
         iHigh = NewiHigh
      ENDDO ! Recursive loop
      IF (DebugLevel.GT.0) WRITE(*,*)
      IF (DebugLevel.GT.0) WRITE(*,100)
      IF (DebugLevel.GT.0) WRITE(*,*)
      100 FORMAT('-------------------------------- ready checking dependencies! ---------------------------------------')
      !
      ! Show the unsorted dependence matrix
      !
      IF (DebugLevel.GT.0) THEN
         WRITE(*,'(A)') 'The following sources are found. Dependencies refer to files in the same order:'
         FormatString = '(A30,1X,     (L1,1X))'
         WRITE(FormatString(9:12),'(I4)') NSources
         DO iSource = 1,NSources
            WRITE(*,FormatString) TRIM(FName(iSource)),(UsesModule(iSource,jSource),jSource = 1,NSources)
         ENDDO
      ENDIF ! debuglevel > 0
      !
      ! Sort the dependence matrix in order of dependence
      !
      IF (DebugLevel.GT.0) WRITE(*,*)
      200 FORMAT('-------------------------------- Ordering dependencies ----------------------------------------------')
      IF (DebugLevel.GT.0) WRITE(*,200)
      IF (DebugLevel.GT.0) WRITE(*,*)
      !
      ! Make a matrix of unresolved dependencies in the list of files
      !
      Unresolved = UsesModule
      DO iSource = 1,NSources
         Resolved(iSource) = .FALSE.
      ENDDO
      Ready = .FALSE.
      SourceCounter = 0
      PrevSourceCounter = 0
      DO WHILE (.NOT.Ready)
         !
         ! Show the unresolved part of the unsorted dependence matrix
         !
         IF (DebugLevel.GT.0) THEN
            WRITE(*,'(A)') 'The following sources have not yet been resolved:'
            FormatString = '(A30,1X,     (L1,1X))'
            WRITE(FormatString(9:12),'(I4)') NSources-SourceCounter
            DO iSource = 1,NSources
               IF (.NOT.Resolved(iSource)) THEN
                  ThisLine = TRIM(FName(iSource))
                  i=0
                  DO jSource = 1,NSources
                     IF (.NOT.Resolved(jSource)) THEN
                        i = i + 1
                        IF (UsesModule(iSource,jSource)) THEN
                           WRITE(ThisLine((32+(i-1)*2):(32+(i-1)*2)),'(A1)') 'X'
                        ELSE
                           WRITE(ThisLine((32+(i-1)*2):(32+(i-1)*2)),'(A1)') '.'
                        ENDIF
                     ENDIF
                  ENDDO
                  WRITE(*,'(A)') TRIM(ThisLine)
               ENDIF
            ENDDO
            WRITE(*,*)
         ENDIF ! debuglevel > 0

         IF (DebugLevel.GT.0) WRITE(*,'(A)') 'Entering another iteration in an attempt to order the sources...'
         !
         ! Count all unresolved dependencies
         !
         DO iSource = 1,NSources
            NUnresolved(iSource) = 0
            DO jSource = 1,NSources
               IF (Unresolved(iSource,jSource)) NUnresolved(iSource) = NUnresolved(iSource) + 1
            ENDDO
         ENDDO
         FormatString = '(        I3)'
         WRITE(FormatString(2:5),'(I4)') NSources
         IF (DebugLevel.GT.0) WRITE(*,FormatString) (NUnresolved(iSource),iSource = 1,NSources)
         !
         ! Put all resolved sources in the output list, as long as they have not already been added
         !
         DO iSource = 1,NSources
            IF ((NUnresolved(iSource).EQ.0).AND.(.NOT.(Resolved(iSource)))) THEN
               IF (DebugLevel.GT.0) WRITE(*,300) TRIM(FName(iSource))
               300 FORMAT('Adding file to the list:                                               ',A)
               SourceCounter = SourceCounter + 1
               SourceOrder(SourceCounter) = iSource
               !
               ! Check this file as resolved
               !
               Resolved(iSource) = .TRUE.
               !
               ! Remove all unresolved dependencies on this file
               !
               DO jSource = 1,NSources
                  Unresolved(jSource,iSource) = .FALSE.
               ENDDO
            ENDIF ! New file added to list
         ENDDO ! loop over all sources
         !
         ! Check if we have added all sourcefiles to the ordered list
         !
         Ready = (SourceCounter.EQ.NSources)
         !
         ! Check if this iteration has resolved at least 1 dependence. If not: terminate furhter execution!
         !
         IF (SourceCounter.EQ.PrevSourceCounter) THEN
            WRITE(*,'(A)') 'This iteration did not resolve any dependence --> terminating further execution!'
            WRITE(*,'(A)') 'Check your MODULE sources for circular references!!'
            CALL EXIT()
         ENDIF
         PrevSourceCounter = SourceCounter
      ENDDO ! iterative loop until all files are added to the list
      !
      ! Show the list:
      !
      IF (DebugLevel.GT.0) THEN
         WRITE(*,*)
         WRITE(*,*)
         WRITE(*,400)
         400 FORMAT('-------------------------------- Final result: ------------------------------------------------------')
         WRITE(*,*)
         WRITE(*,110)
         110 FORMAT('The final list of sources in order of dependence is:')
         WRITE(*,*)
         DO iSource = 1,SourceCounter !NSources
            WRITE(*,'(A)') TRIM(FName(SourceOrder(iSource)))
         ENDDO
         WRITE(*,*)

         WRITE(*,'(A)') &
         & 'Source-dependence matrix in order of dependence. Dependencies refer to files in the same order:'
         FormatString = '(A30,1X,     (A1,1X))'
         WRITE(FormatString(9:12),'(I4)') NSources
         DO iSource = 1,NSources
            ThisLine = TRIM(FName(SourceOrder(iSource)))
            DO jSource = 1,NSources
               IF (UsesModule(SourceOrder(iSource),SourceOrder(jSource))) THEN
                  WRITE(ThisLine((32+(jSource-1)*2):(32+(jSource-1)*2)),'(A1)') 'X'
               ELSE
                  WRITE(ThisLine((32+(jSource-1)*2):(32+(jSource-1)*2)),'(A1)') '.'
               ENDIF
            ENDDO
            WRITE(*,'(A)') TRIM(ThisLine)
         ENDDO

         WRITE(*,*)
         WRITE(*,'(A)') 'Output string = '
      ENDIF

      OutLine = ' '
      DO iSource = 1,SourceCounter !NSources
         OutLine = TRIM(OutLine)//' '//TRIM(FName(SourceOrder(iSource)))
      ENDDO
      WRITE(*,'(A)') TRIM(OutLine)

   ENDIF ! More than 0 commandline parameters specified
END PROGRAM GetDeps


