program nuclide_decay
   !
   ! This utility constructs the decay matrix for (a vector of) all ENDF nuclides for a given arbitrary delay.
   ! All different delay pinpoints used in the cocktail DCC library are handled in 1 go!
   !
   use libxmath
   use libutil
   use libinterval
   use libendf
   use Matrix_Exponential

   implicit none
   !
   ! Declarations
   !
   integer, parameter :: NIntervals = 16
   type(NuclideType), dimension(0:MaxNuclides) :: RegularizedNuclideSpecs
   real(Float), dimension(MaxNuclides,MaxNuclides) :: RegularizedMotherDaughterMatrix,TheDecayMatrix
   character(6), dimension(4), parameter :: TimeUnitName = (['second','minute','hour  ','day   '])
   real(Float) :: tMin,Lambda,Yield,Dum,tMax
   integer :: MotherNuclide,DaughterNuclide,iNuclide,jNuclide,NArguements,NFirstMotherOrphans,iStartingTime,MyNuclide
   integer, dimension(NIntervals) :: tStarType
   logical :: IsUnbranchedFirstDaughter,HasIAEAEntry,HasRP122Entry,DoPrint,DoFull,JustStarted,&
   & ENDFBinFileExists
   real(dp), dimension(:,:), allocatable :: x,zDelay
   character(DefaultLength) :: MyTime,ALine,FullSparseString
   character(DefaultLength), dimension(20) :: RunlogText
   !
   ! Specify how the pinpoints i=1,N are made: t_i = t0*G^(i-1)
   ! There is also pinpoint 0, which is 0
   !
   real(Float), parameter :: FirstDelay = 60._Float ! seconds
   real(Float), parameter :: DelayGrowthFactor = 1.15_Float
   integer, parameter :: NStartingTimes = 200
   real(Float), parameter :: VerySmallNumber = 1.E-95_Float

   character(6), dimension(0:NStartingTimes) :: PinPointName
   character(DefaultLength), dimension(0:NStartingTimes) :: FName
   real(Float), dimension(0:NStartingTimes) :: AvailableDelay ! In seconds

   integer, parameter :: DebugLevel = 0 ! 0: no comments, 1: some comments, 2: verbose
   !
   ! Specify locations of auxiliary files and of output via commandline parameters
   !
   write(*,'(A)') '==================================================================================================='
   write(*,*)
   write(*,'(A)') 'This utility was developed for:'
   write(*,*)
   write(*,'(A)') 'VLH (Centre for Environmental Safety and Security)'
   write(*,'(A)') 'RIVM (National Institute for Public Health and the Environment)'
   write(*,'(A)') 'Visiting address: A. van Leeuwenhoeklaan 9, Bilthoven, The Netherlands'
   write(*,'(A)') 'Mail address:     Interne Postbak 21, Postbus 1, 3720 BA Bilthoven, The Netherlands'
   write(*,'(A)') 'Website:          http://www.rivm.com/'
   write(*,*)
   write(*,'(A)') '==================================================================================================='

   NArguements = IARGC()

   if (NArguements /= 1) then
      write(*,'(A)') 'Call:'
      write(*,*)
      write(*,'(A)') 'nuclide_decay <full/sparse>'
      write(*,*)
      write(*,'(A)') 'where:'
      write(*,'(A)') '<full/sparse> can be used to select full (files 31MB) or sparse (files 115kB) format.'
      call exit()
   endif

   call getarg(1,FullSparseString)
   write(*,'(A)') 'You specified the following FullSparseString : '//trim(FullSparseString)
   DoFull = (INDEX(FullSparseString,'full') > 0)
   !
   ! Initialize the set of delays
   !
   AvailableDelay(0) = 0
   AvailableDelay(1) = FirstDelay
   do iStartingTime = 2,NStartingTimes
      AvailableDelay(iStartingTime) = DelayGrowthFactor*AvailableDelay(iStartingTime-1)
   enddo ! loop over starting times

   do iStartingTime = 0,NStartingTimes
      write(PinPointName(iStartingTime),'(I3.3)') iStartingTime
      if (DebugLevel > 0) write(*,'(A)') PinPointName(iStartingTime)//'second'

      if (DoFull) then
         FName(iStartingTime) = TransitionMatrixPath() // '/' // 'DecayMatrix_ENDF_' // trim(PinPointName(iStartingTime)) // '.dat'
         write(*,'(A)') 'Going to write fill matrix output to file "'//trim(FName(iStartingTime))//'"'
      else
         FName(iStartingTime) = TransitionMatrixPath() // '/' // 'SparseMatrix_ENDF_' // trim(PinPointName(iStartingTime)) // '.dat'
         write(*,'(A)') 'Going to write sparse matrix output to file "'//trim(FName(iStartingTime))//'"'
      endif
   enddo ! loop over starting times
   !
   ! Set shortest halflife that is admitted to the numerical scheme for solving the Bateman equations.
   ! Shorter ones are skipped: the grandmother is directly transformed into the grand-daughter.
   ! The missing decay steps are added afterwards.
   !
   tMin = 10._Float ! [s]

   write(*,'(A)') 'Set shortest halflife that is admitted to the numerical scheme for solving the Bateman equations.'
   write(*,'(A)') 'Shorter ones are skipped: the grandmother is directly transformed into the grand-daughter.'
   write(*,'(A)') 'The missing decay steps are added afterwards.'
   write(*,'(A,F15.5,A)') 'For this run, we choose tMin = ',tMin,' [s]'
   !
   ! Read ENDF nuclide database
   !
   call ReadNProcessENDFNuclideSpecs(tMin,RegularizedNuclideSpecs,RegularizedMotherDaughterMatrix)

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

   do iStartingTime = 0,NStartingTimes
      open(ScrotchFile,file=FName(iStartingTime),form='FORMATTED',action='WRITE',position='REWIND')

      write(ScrotchFile,'(I0,A)') NStartingTimes,' NStartingTimes'
      write(ScrotchFile,'(F15.5,A)') FirstDelay,' First delay [s]'
      write(ScrotchFile,'(F15.10,A)') DelayGrowthFactor,' DelayGrowthFactor'
      write(ScrotchFile,'(I0,A)') iStartingTime,' Power of the delay growth factor (0 = First delay)'
      write(ScrotchFile,'(EN20.10,A)') AvailableDelay(iStartingTime),' Delay [s] associated with this file'
      write(ScrotchFile,*)

      if (DoFull) then
         write(ScrotchFile,'(A15,4000A20)') 'Mother Daughter ',&
         & (trim(NuclideSpecs(jNuclide)%NuclideName),jNuclide = 1,NNuclides)
      else
         write(ScrotchFile,'(A)') 'Mother Daughter Factor[Bq/Bq]'
      endif ! full or sparse
      close(ScrotchFile)
   enddo ! loop over pinpoints to write header


   do iNuclide = 1,NNuclides
      MyNuclide = MyNuclide + 1

      write(*,*)
      write(*,'(A)') '######################################################################'
      write(*,*)
      write(*,'(A,I4,A,EN15.5,A)') 'Head of Chain nuclide ',iNuclide,&
         & ': '//trim(NuclideSpecs(iNuclide)%NuclideName)
      write(*,*)

      if (DebugLevel > 3) write(*,'(A,I0)') 'DEBUG: iNuclide = ',iNuclide
      Lambda = LOG(2._Float)/NuclideSpecs(iNuclide)%HalfTime
      !
      ! Only "long-lived" nuclides can be the mater familias, contamination with fast nuclides is not considered to be a contamination, since it resolves itself fast
      !
      TheDecayMatrix = 0._Float

      if (NuclideSpecs(iNuclide)%HalfTime >= Orphanage%tMin) then
         write(*,'(A)') '----------------------------------------------------------------------------------------------'
         !
         ! Construct a list of all progeny from this head of chain
         !
         Orphanage%OrphanActive = .FALSE. ! No orphans in this family yet
         DoPrint = .TRUE.
         call CollectProgeny(iNuclide,DoPrint)
         !
         ! List all transitions in this family, including progeny
         !
         if (NuclideSpecs(iNuclide)%HalfTime >= 1.E85_Float) then
            write(*,'(A)') 'This is a stable nuclide, so no transitions!'
         else
            call ListTransitions()
         endif

         !#######################################################################
         !
         ! Construct an evolution matrix for this family to account for all decay transformations in the given delay time
         !
         do iStartingTime = 0,NStartingTimes

            TheDecayMatrix = 0._Float

            call MakeEvolutionMatrix(RegularizedMotherDaughterMatrix,x,AvailableDelay(iStartingTime),zDelay)
            !
            ! Add decay contribution of this mother to output decay matrix
            !
            MotherNuclide = 1 ! Head of chain only!

            do DaughterNuclide = 1,NuclideFamily%NFamily
               TheDecayMatrix(NuclideFamily%FamilyMember(DaughterNuclide),iNuclide) &
               & = zDelay(DaughterNuclide,MotherNuclide)
            enddo
            !
            ! Add activities for orphan nuclides
            !
            call AddOrphansBelow(iNuclide,TheDecayMatrix,zDelay,NFirstMotherOrphans)
            !
            ! Save results in regular full or sparse format
            !
            do jNuclide = 1,NNuclides
               if (abs(TheDecayMatrix(jNuclide,iNuclide)) < VerySmallNumber) TheDecayMatrix(jNuclide,iNuclide) = 0._Float
            enddo

            open(ScrotchFile,file=FName(iStartingTime),form='FORMATTED',action='WRITE',position='APPEND')

            if (DoFull) then
               write(ScrotchFile,'(A20,1300(G19.10,1X))') trim(NuclideSpecs(iNuclide)%NuclideName),&
               & (TheDecayMatrix(jNuclide,iNuclide),jNuclide = 1,NNuclides)
            else
               do jNuclide = 1,NNuclides
                  if (TheDecayMatrix(jNuclide,iNuclide) /= 0._Float) then
                     write(ScrotchFile,'(I0,1X,I0,1X,G17.10)') iNuclide,jNuclide,TheDecayMatrix(jNuclide,iNuclide)
                  endif
               enddo
            endif ! full or sparse

            close(ScrotchFile)

         enddo ! loop over pinpoints in time
         !#######################################################################
      else
         write(*,'(A)') 'This is an orphan, so I am not considering it as a head of chain!'
         !#######################################################################
         !
         ! Write zeros for this nuclide, assuming that progeny of orphan nuclides is negligible (which may not be the case...)
         !
         do iStartingTime = 0,NStartingTimes
            do jNuclide = 1,NNuclides
               if (abs(TheDecayMatrix(jNuclide,iNuclide)) < VerySmallNumber) TheDecayMatrix(jNuclide,iNuclide) = 0._Float
            enddo

            open(ScrotchFile,file=FName(iStartingTime),form='FORMATTED',action='WRITE',position='APPEND')

            if (DoFull) then
               write(ScrotchFile,'(A20,1300(G19.10,1X))') trim(NuclideSpecs(iNuclide)%NuclideName),&
               & (TheDecayMatrix(jNuclide,iNuclide),jNuclide = 1,NNuclides)
            else
               do jNuclide = 1,NNuclides
                  if (TheDecayMatrix(jNuclide,iNuclide) /= 0._Float) then
                     write(ScrotchFile,'(I0,1X,I0,1X,G17.10)') iNuclide,jNuclide,TheDecayMatrix(jNuclide,iNuclide)
                  endif
               enddo
            endif ! full or sparse
            close(ScrotchFile)

         enddo ! loop over pinpoints in time
         !#######################################################################
      endif ! halflife > Orphanage%tMin

   enddo ! loop over nuclides
end program nuclide_decay
