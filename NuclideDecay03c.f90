PROGRAM NuclideDecay
   !
   ! This utility constructs the decay matrix for (a vector of) all ENDF nuclides for a given arbitrary delay.
   ! All different delay pinpoints used in the cocktail DCC library are handled in 1 go!
   !
   USE libxmath
   USE libutil
   USE libinterval
   USE libendf
   USE libexponential

   IMPLICIT NONE
   !
   ! Declarations
   !
   INTEGER, PARAMETER :: NIntervals = 16
   TYPE(NuclideType), DIMENSION(0:MaxNuclides) :: RegularizedNuclideSpecs
   REAL(Float), DIMENSION(MaxNuclides,MaxNuclides) :: RegularizedMotherDaughterMatrix,TheDecayMatrix
   CHARACTER(6), DIMENSION(4), PARAMETER :: TimeUnitName = ((/'second','minute','hour  ','day   '/))
   REAL(Float) :: tMin,Lambda,Yield,Dum,tMax
   INTEGER :: MotherNuclide,DaughterNuclide,iNuclide,jNuclide,NArguements,NFirstMotherOrphans,iStartingTime,MyNuclide
   INTEGER, DIMENSION(NIntervals) :: tStarType
   LOGICAL :: IsUnbranchedFirstDaughter,HasIAEAEntry,HasRP122Entry,DoPrint,DoFull,JustStarted,&
   & ENDFBinFileExists
   REAL(dp), DIMENSION(:,:), ALLOCATABLE :: x,zDelay
   CHARACTER(DefaultLength) :: MyTime,ALine,FullSparseString
   CHARACTER(DefaultLength), DIMENSION(20) :: RunlogText
   !
   ! Specify how the pinpoints i=1,N are made: t_i = t0*G^(i-1)
   ! There is also pinpoint 0, which is 0
   !
   REAL(Float), PARAMETER :: FirstDelay = 60._Float ! seconds
   REAL(Float), PARAMETER :: DelayGrowthFactor = 1.15_Float
   INTEGER, PARAMETER :: NStartingTimes = 200
   REAL(Float), PARAMETER :: VerySmallNumber = 1.E-95_Float

   CHARACTER(6), DIMENSION(0:NStartingTimes) :: PinPointName
   CHARACTER(DefaultLength), DIMENSION(0:NStartingTimes) :: FName
   REAL(Float), DIMENSION(0:NStartingTimes) :: AvailableDelay ! In seconds

   INTEGER, PARAMETER :: DebugLevel = 0 ! 0: no comments, 1: some comments, 2: verbose
   !
   ! Specify locations of auxiliary files and of output via commandline parameters
   !
   WRITE(*,'(A)') '==================================================================================================='
   WRITE(*,*)
   WRITE(*,'(A)') 'This is utility NuclideDecay, version 14 March 2024,'
   WRITE(*,'(A)') 'based on the ~3800 nuclides from the ENDF dataset.'
   WRITE(*,*)
   WRITE(*,'(A)') 'developed for:'
   WRITE(*,*)
   WRITE(*,'(A)') 'VLH (Centre for Environmental Safety and Security)'
   WRITE(*,'(A)') 'RIVM (National Institute for Public Health and the Environment)'
   WRITE(*,'(A)') 'Visiting address: A. van Leeuwenhoeklaan 9, Bilthoven, The Netherlands'
   WRITE(*,'(A)') 'Mail address:     Interne Postbak 21, Postbus 1, 3720 BA Bilthoven, The Netherlands'
   WRITE(*,'(A)') 'Website:          http://www.rivm.com/'
   WRITE(*,*)
   WRITE(*,'(A)') '==================================================================================================='

   NArguements = IARGC()

   IF (NArguements.NE.1) THEN
      WRITE(*,'(A)') 'Call:'
      WRITE(*,*)
      WRITE(*,'(A)') 'NuclideDecay <full/sparse>'
      WRITE(*,*)
      WRITE(*,'(A)') 'where:'
      WRITE(*,'(A)') '<full/sparse> can be used to select full (files 31MB) or sparse (files 115kB) format.'
      CALL EXIT()
   ENDIF

   CALL GETARG(1,FullSparseString)
   WRITE(*,'(A)') 'You specified the following FullSparseString : '//TRIM(FullSparseString)
   DoFull = (INDEX(FullSparseString,'full').GT.0)
   !
   ! Initialize the set of delays
   !
   AvailableDelay(0) = 0
   AvailableDelay(1) = 60
   DO iStartingTime = 2,NStartingTimes
      AvailableDelay(iStartingTime) = DelayGrowthFactor*AvailableDelay(iStartingTime-1)
   ENDDO ! loop over starting times

   DO iStartingTime = 0,NStartingTimes
      WRITE(PinPointName(iStartingTime),'(I3.3)') iStartingTime
      IF (DebugLevel.GT.0) WRITE(*,'(A)') PinPointName(iStartingTime)//'second'

      IF (DoFull) THEN
         FName(iStartingTime) = 'DecayMatrix_ENDF_'//TRIM(PinPointName(iStartingTime))//'.dat'
         WRITE(*,'(A)') 'Going to write fill matrix output to file "'//TRIM(FName(iStartingTime))//'"'
      ELSE
         FName(iStartingTime) = 'SparseMatrix_ENDF_'//TRIM(PinPointName(iStartingTime))//'.dat'
         WRITE(*,'(A)') 'Going to write sparse matrix output to file "'//TRIM(FName(iStartingTime))//'"'
      ENDIF
   ENDDO ! loop over starting times
   !
   ! Set shortest halflife that is admitted to the numerical scheme for solving the Bateman equations.
   ! Shorter ones are skipped: the grandmother is directly transformed into the grand-daughter.
   ! The missing decay steps are added afterwards.
   !
   tMin = 10._Float ! [s]

   WRITE(*,'(A)') 'Set shortest halflife that is admitted to the numerical scheme for solving the Bateman equations.'
   WRITE(*,'(A)') 'Shorter ones are skipped: the grandmother is directly transformed into the grand-daughter.'
   WRITE(*,'(A)') 'The missing decay steps are added afterwards.'
   WRITE(*,'(A,F15.5,A)') 'For this run, we choose tMin = ',tMin,' [s]'
   !
   ! Read ENDF nuclide database
   !
   CALL ReadNProcessENDFNuclideSpecs(tMin,RegularizedNuclideSpecs,RegularizedMotherDaughterMatrix)

   !################################################################################################################
   !################################################################################################################
   !################################################################################################################
   !################################################################################################################
   !################################################################################################################
   !################################################################################################################
   !
   ! Loop over all nuclides as "mater familias" = "head of chain"
   !
   MyNuclide = 0 ! Independent counter that always starts with 0, irrespective of any debugging with iNuclide

   DO iStartingTime = 0,NStartingTimes
      OPEN(ScrotchFile,FILE=FName(iStartingTime),FORM='FORMATTED',ACTION='WRITE',POSITION='REWIND')

      WRITE(ScrotchFile,'(I0,A)') NStartingTimes,' NStartingTimes'
      WRITE(ScrotchFile,'(F15.5,A)') FirstDelay,' First delay [s]'
      WRITE(ScrotchFile,'(F15.10,A)') DelayGrowthFactor,' DelayGrowthFactor'
      WRITE(ScrotchFile,'(I0,A)') iStartingTime,' Power of the delay growth factor (0 = First delay)'
      WRITE(ScrotchFile,'(EN20.10,A)') AvailableDelay(iStartingTime),' Delay [s] associated with this file'
      WRITE(ScrotchFile,*)

      IF (DoFull) THEN
         WRITE(ScrotchFile,'(A15,4000A20)') 'Mother Daughter ',&
         & (TRIM(NuclideSpecs(jNuclide)%NuclideName),jNuclide = 1,NNuclides)
      ELSE
         WRITE(ScrotchFile,'(A)') 'Mother Daughter Factor[Bq/Bq]'
      ENDIF ! full or sparse
      CLOSE(ScrotchFile)
   ENDDO ! loop over pinpoints to write header


   DO iNuclide = 1,NNuclides
      MyNuclide = MyNuclide + 1

      WRITE(*,*)
      WRITE(*,'(A)') '######################################################################'
      WRITE(*,*)
      WRITE(*,'(A,I4,A,EN15.5,A)') 'Head of Chain nuclide ',iNuclide,&
         & ': '//TRIM(NuclideSpecs(iNuclide)%NuclideName)
      WRITE(*,*)

      IF (DebugLevel.GT.3) WRITE(*,'(A,I0)') 'DEBUG: iNuclide = ',iNuclide
      Lambda = LOG(2._Float)/NuclideSpecs(iNuclide)%HalfTime
      !
      ! Only "long-lived" nuclides can be the mater familias, contamination with fast nuclides is not considered to be a contamination, since it resolves itself fast
      !
      TheDecayMatrix = 0._Float

      IF (NuclideSpecs(iNuclide)%HalfTime.GE.Orphanage%tMin) THEN
         WRITE(*,'(A)') '----------------------------------------------------------------------------------------------'
         !
         ! Construct a list of all progeny from this head of chain
         !
         Orphanage%OrphanActive = .FALSE. ! No orphans in this family yet
         DoPrint = .TRUE.
         CALL CollectProgeny(iNuclide,DoPrint)
         !
         ! List all transitions in this family, including progeny
         !
         IF (NuclideSpecs(iNuclide)%HalfTime.GE.1.E85_Float) THEN
            WRITE(*,'(A)') 'This is a stable nuclide, so no transitions!'
         ELSE
            CALL ListTransitions()
         ENDIF

         !#######################################################################
         !
         ! Construct an evolution matrix for this family to account for all decay transformations in the given delay time
         !
         DO iStartingTime = 0,NStartingTimes

            TheDecayMatrix = 0._Float

            CALL MakeEvolutionMatrix(iNuclide,RegularizedMotherDaughterMatrix,x,AvailableDelay(iStartingTime),zDelay)
            !
            ! Add decay contribution of this mother to output decay matrix
            !
            MotherNuclide = 1 ! Head of chain only!

            DO DaughterNuclide = 1,NuclideFamily%NFamily
               TheDecayMatrix(NuclideFamily%FamilyMember(DaughterNuclide),iNuclide) &
               & = zDelay(DaughterNuclide,MotherNuclide)
            ENDDO
            !
            ! Add activities for orphan nuclides
            !
            CALL AddOrphansBelow(iNuclide,TheDecayMatrix,zDelay,NFirstMotherOrphans)
            !
            ! Save results in regular full or sparse format
            !
            DO jNuclide = 1,NNuclides
               IF (ABS(TheDecayMatrix(jNuclide,iNuclide)).LT.VerySmallNumber) TheDecayMatrix(jNuclide,iNuclide) = 0._Float
            ENDDO

            OPEN(ScrotchFile,FILE=FName(iStartingTime),FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')

            IF (DoFull) THEN
               WRITE(ScrotchFile,'(A20,1300(G19.10,1X))') TRIM(NuclideSpecs(iNuclide)%NuclideName),&
               & (TheDecayMatrix(jNuclide,iNuclide),jNuclide = 1,NNuclides)
            ELSE
               DO jNuclide = 1,NNuclides
                  IF (TheDecayMatrix(jNuclide,iNuclide).NE.0._Float) THEN
                     WRITE(ScrotchFile,'(I0,1X,I0,1X,G17.10)') iNuclide,jNuclide,TheDecayMatrix(jNuclide,iNuclide)
                  ENDIF
               ENDDO
            ENDIF ! full or sparse

            CLOSE(ScrotchFile)

         ENDDO ! loop over pinpoints in time
         !#######################################################################
      ELSE
         WRITE(*,'(A)') 'This is an orphan, so I am not considering it as a head of chain!'
         !#######################################################################
         !
         ! Write zeros for this nuclide, assuming that progeny of orphan nuclides is negligible (which may not be the case...)
         !
         DO iStartingTime = 0,NStartingTimes
            DO jNuclide = 1,NNuclides
               IF (ABS(TheDecayMatrix(jNuclide,iNuclide)).LT.VerySmallNumber) TheDecayMatrix(jNuclide,iNuclide) = 0._Float
            ENDDO

            OPEN(ScrotchFile,FILE=FName(iStartingTime),FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')

            IF (DoFull) THEN
               WRITE(ScrotchFile,'(A20,1300(G19.10,1X))') TRIM(NuclideSpecs(iNuclide)%NuclideName),&
               & (TheDecayMatrix(jNuclide,iNuclide),jNuclide = 1,NNuclides)
            ELSE
               DO jNuclide = 1,NNuclides
                  IF (TheDecayMatrix(jNuclide,iNuclide).NE.0._Float) THEN
                     WRITE(ScrotchFile,'(I0,1X,I0,1X,G17.10)') iNuclide,jNuclide,TheDecayMatrix(jNuclide,iNuclide)
                  ENDIF
               ENDDO
            ENDIF ! full or sparse
            CLOSE(ScrotchFile)

         ENDDO ! loop over pinpoints in time
         !#######################################################################
      ENDIF ! halflife > Orphanage%tMin

   ENDDO ! loop over nuclides
END PROGRAM NuclideDecay
