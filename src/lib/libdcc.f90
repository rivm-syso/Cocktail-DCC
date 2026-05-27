module LibDCC
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
   use LibENDF, only: RIVMSourcesPath, isGas, &
      ReadNProcessENDFNuclideSpecs, ReadNuclideSpecs, NuclideSpecs, GetNuclideNumber, NuclideType, MaxNuclides
   use LibUtil, only: DefaultLength, ScratchFile, &
      Male, Female, GenderName
   use LibXMath, only: Float
   implicit none

   private

   public :: InitLibDCC,ReadTissueDCCs,ReadDCCEffectiveDose,ReadGroundDCCs,DCCEffectiveDose,DCCGround,&
   & InAir,Age_Adult,DCCInhalation,ReadInhalationDCCs,Inhalation_public_adult,&
   & DCCThyroidInhalation,ReadThyroidInhalationDCCs

   character(DefaultLength) :: ProjectPath = './'

   character(DefaultLength) :: DCCAirSubmersionPath,DCCSoilContaminationPath,DCCWaterImmersionPath,&
   & DCCInhalationPath

   logical :: OPTION_NobleGroundShine = .False.

   integer, parameter :: InAir = 1
   integer, parameter :: InWater = 2

   integer, parameter :: Tissue_R_marrow   =  1
   integer, parameter :: Tissue_Colon      =  2
   integer, parameter :: Tissue_Lungs      =  3
   integer, parameter :: Tissue_ST_wall    =  4
   integer, parameter :: Tissue_Breast     =  5
   integer, parameter :: Tissue_Ovaries    =  6
   integer, parameter :: Tissue_Testes     =  7
   integer, parameter :: Tissue_UB_wall    =  8
   integer, parameter :: Tissue_Oesophagus =  9
   integer, parameter :: Tissue_Liver      = 10
   integer, parameter :: Tissue_Thyroid    = 11
   integer, parameter :: Tissue_Endost_BS  = 12
   integer, parameter :: Tissue_Brain      = 13
   integer, parameter :: Tissue_S_glands   = 14
   integer, parameter :: Tissue_Skin       = 15
   integer, parameter :: Tissue_Adrenals   = 16
   integer, parameter :: Tissue_ET         = 17
   integer, parameter :: Tissue_GB_wall    = 18
   integer, parameter :: Tissue_Ht_wall    = 19
   integer, parameter :: Tissue_Kidneys    = 20
   integer, parameter :: Tissue_Lymph      = 21
   integer, parameter :: Tissue_Muscle     = 22
   integer, parameter :: Tissue_O_mucosa   = 23
   integer, parameter :: Tissue_Pancreas   = 24
   integer, parameter :: Tissue_Prostate   = 25
   integer, parameter :: Tissue_SI_wall    = 26
   integer, parameter :: Tissue_Spleen     = 27
   integer, parameter :: Tissue_Thymus     = 28
   integer, parameter :: Tissue_Uterus     = 29

   integer, parameter :: NTissues          = 29

   character(10), parameter, dimension(NTissues) :: TissueName = &
      & ['R-marrow  ',&
      &   'Colon     ',&
      &   'Lungs     ',&
      &   'ST-wall   ',&
      &   'Breast    ',&
      &   'Ovaries   ',&
      &   'Testes    ',&
      &   'UB-wall   ',&
      &   'Oesophagus',&
      &   'Liver     ',&
      &   'Thyroid   ',&
      &   'Endost-BS ',&
      &   'Brain     ',&
      &   'S-glands  ',&
      &   'Skin      ',&
      &   'Adrenals  ',&
      &   'ET        ',&
      &   'GB-wall   ',&
      &   'Ht-wall   ',&
      &   'Kidneys   ',&
      &   'Lymph     ',&
      &   'Muscle    ',&
      &   'O-mucosa  ',&
      &   'Pancreas  ',&
      &   'Prostate  ',&
      &   'SI-wall   ',&
      &   'Spleen    ',&
      &   'Thymus    ',&
      &   'Uterus    ']

   integer, parameter :: Age_Adult   = 1
   integer, parameter :: Age_15yr    = 2
   integer, parameter :: Age_10yr    = 3
   integer, parameter :: Age_5yr     = 4
   integer, parameter :: Age_1yr     = 5
   integer, parameter :: Age_Newborn = 6

   integer, parameter :: NDCCAges = 6

   integer, parameter :: NGroundDepths = 9

   character(45), parameter, dimension(NGroundDepths) :: GroundDepthName = &
      & ['Planar sources at specific depths/0.5 g.cm-2/',&
      &   'Exponential sources/0.5 g.cm-2/              ',&
      &   'Exponential sources/1.0 g.cm-2/              ',&
      &   'Exponential sources/2.5 g.cm-2/              ',&
      &   'Exponential sources/5.0 g.cm-2/              ',&
      &   'Exponential sources/10.0 g.cm-2/             ',&
      &   'Exponential sources/20.0 g.cm-2/             ',&
      &   'Exponential sources/50.0 g.cm-2/             ',&
      &   'Exponential sources/100.0 g.cm-2/            ']

   integer, parameter :: NInhalationDCCs = 8

   integer, parameter :: Inhalation_worker_1um     = 1
   integer, parameter :: Inhalation_worker_5um     = 2
   integer, parameter :: Inhalation_public_infant  = 3
   integer, parameter :: Inhalation_public_1_year  = 4
   integer, parameter :: Inhalation_public_5_year  = 5
   integer, parameter :: Inhalation_public_10_year = 6
   integer, parameter :: Inhalation_public_15_year = 7
   integer, parameter :: Inhalation_public_adult   = 8

   type DCCVector
      real(Float), dimension(MaxNuclides) :: x
   end type DCCVector

   type(DCCVector), dimension(InAir:InWater,Male:Female,NTissues,NDCCAges) :: DCCTissue
   type(DCCVector), dimension(NGroundDepths,NDCCAges) :: DCCGround
   type(DCCVector), dimension(NGroundDepths) :: DCCEquivalentDoseGround,DCCKermaRateGround
   type(DCCVector), dimension(InAir:InWater,NDCCAges) :: DCCEffectiveDose
   type(DCCVector) :: DCCEquivalentDoseInAir,DCCKermaRateInAir,DCCThyroidInhalation
   type(DCCVector), dimension(NInhalationDCCs) :: DCCInhalation

   type(NuclideType), dimension(0:MaxNuclides) :: RegularizedNuclideSpecs
   real(Float), dimension(MaxNuclides,MaxNuclides) :: RegularizedMotherDaughterMatrix

contains

   subroutine InitLibDCC(UseICRP)
      !
      ! Initialize this library
      !
      use LibUtil, only: env_var

      logical, intent(in) :: UseICRP

      character(:), allocatable :: DCCCalcPath
      real(Float) :: tMin
      integer, parameter :: DebugLevel = 0

      if (UseICRP) then
         call ReadNuclideSpecs(RIVMSourcesPath() // '/ICRP-07.NDX')
      else
         !
         ! Set shortest halflife that is admitted to the numerical scheme for solving the Bateman equations.
         ! Shorter ones are skipped: the grandmother is directly transformed into the grand-daughter.
         ! The missing decay staps are added afterwards.
         !
         tMin = 10._Float ! [s]

         if (DebugLevel > 0) then
            write(*,'(A)') 'Set shortest halflife admitted to the numerical scheme for solving the Bateman equations.'
            write(*,'(A)') 'Shorter ones are skipped: the grandmother is directly transformed into the grand-daughter.'
            write(*,'(A)') 'The missing decay steps are added afterwards.'
            write(*,'(A,F15.5,A)') 'For this run, we choose tMin = ',tMin,' [s]'
         endif

         call ReadNProcessENDFNuclideSpecs(tMin,RegularizedNuclideSpecs,RegularizedMotherDaughterMatrix)
      endif

      call env_var('COCKTAIL_DCC_ICRP_SJ_DIR', DCCCalcPath)
      if (.not. allocated(DCCCalcPath)) error stop 'need to set environment variable: COCKTAIL_DCC_ICRP_SJ_DIR'
      DCCAirSubmersionPath = trim(DCCCalcPath) // '/Air submersion/'
      DCCSoilContaminationPath = trim(DCCCalcPath) // '/Soil contamination/'
      DCCWaterImmersionPath = trim(DCCCalcPath) // '/Water immersion/'
      DCCInhalationPath = trim(ProjectPath)//'resources/'
   end subroutine InitLibDCC

   subroutine ReadTissueDCCs()
      !
      ! Read tissue- and gender specific DCCs
      ! Attention: This file has only the ICRP nuclides, but the reference in this model may be ENDF,
      ! which has ~ 3 x more nuclides!
      !
      character(DefaultLength) :: FName,ALine
      integer :: iGender,iLine,jNuclide,iTissue,iAge,iImmersion,NDashLinesFound
      character(10) :: MyTissueName,MyNuclideName
      logical :: IsDashLine

      do iImmersion = InAir,InWater
         do iGender = Male,Female
            FName = 'Nuclide-specific_EquivalentDose_'//trim(GenderName(iGender))//'.txt'
            if (iImmersion == InAir) then
               FName = 'Air_'//trim(FName)
               FName = trim(DCCAirSubmersionPath)//trim(FName)
            else
               FName = 'Water_'//trim(FName)
               FName = trim(DCCWaterImmersionPath)//trim(FName)
            endif
            open(ScratchFile,file=FName,form='FORMATTED',action='READ')
            !
            ! Skip header
            !
            NDashLinesFound = 0
            iLine = 0

            do while (NDashLinesFound < 2)
               read(ScratchFile,'(A)') ALine
               iLine = iLine + 1
               if (ALine(1:1) == '-') NDashLinesFound = NDashLinesFound + 1
            enddo
            !
            ! Read data for all nuclides
            !
            IsDashLine = .FALSE.

            do while (.NOT.IsDashLine)
               do iTissue = 1,NTissues
                  read(ScratchFile,'(A)',end=10) ALine
                  iLine = iLine + 1

                  IsDashLine = (INDEX(ALine,'-----') > 0)

                  if (.NOT.IsDashLine) then
                     !
                     ! Get nuclide and check name
                     !
                     if (iTissue == 1) then
                        MyNuclideName = ALine(1:7)
                        jNuclide = GetNuclideNumber(MyNuclideName)
                        if (jNuclide == 0) then
                           write(*,'(A,I0,A,A)') 'At line ',iLine,' in file '//trim(FName),&
                           & ' : found nuclide '//MyNuclideName//', but I cannot give it a number... Exiting!!'
                           call exit()
                        endif
                     endif ! First organ
                     !
                     ! Get and check tissue name
                     !
                     MyTissueName = ALine(9:18)
                     if (MyTissueName /= TissueName(iTissue)) then
                        write(*,'(A,I0,A)') 'At line ',iLine,&
                        & ' : expected tissue '//TissueName(iTissue)//', but found '//MyTissueName//'! Exiting!!'
                        call exit()
                     endif
                     !
                     ! Read values
                     !
                     read(ALine(25:DefaultLength),*) (DCCTissue(iImmersion,iGender,iTissue,iAge)%x(jNuclide),iAge=1,NDCCAges)
                  endif ! this is not a dashline
               enddo ! loop over tissues
               read(ScratchFile,*) ! Skip line with dashes
               iLine = iLine + 1
            enddo ! loop over nuclides

            10    continue

            close(ScratchFile)
         enddo ! loop over gender
      enddo ! loop over in air and in water
   end subroutine ReadTissueDCCs



   subroutine ReadDCCEffectiveDose()
      !
      ! Read effective dose DCCs
      ! Attention: This file has only the ICRP nuclides, but the reference in this model may be ENDF,
      ! which has ~ 3 x more nuclides!
      !
      character(DefaultLength) :: FName,ALine
      integer :: iLine,jNuclide,iAge,NDashLinesFound,iImmersion
      character(10) :: MyNuclideName
      logical :: IsDashLine

      do iImmersion = InAir,InWater
         FName = 'Nuclide-specific_EffectiveDose.txt'
         if (iImmersion == InAir) then
            FName = 'Air_'//trim(FName)
            FName = trim(DCCAirSubmersionPath)//trim(FName)
         else
            FName = 'Water_'//trim(FName)
            FName = trim(DCCWaterImmersionPath)//trim(FName)
         endif

         open(ScratchFile,file=FName,form='FORMATTED',action='READ')
         !
         ! Skip header
         !
         NDashLinesFound = 0
         iLine = 0

         do while (NDashLinesFound < 2)
            read(ScratchFile,'(A)') ALine
            iLine = iLine + 1
            if (ALine(1:1) == '-') NDashLinesFound = NDashLinesFound + 1
         enddo
         !
         ! Read data for all nuclides
         !
         do
            read(ScratchFile,'(A)',end=10) ALine
            iLine = iLine + 1

            IsDashLine = (INDEX(ALine,'-----') > 0)

            if (.NOT.IsDashLine) then
               !
               ! Get nuclide and check name
               !
               MyNuclideName = ALine(1:7)
               jNuclide = GetNuclideNumber(MyNuclideName)
               if (jNuclide == 0) then
                  write(*,'(A,I0,A,A)') 'At line ',iLine,' in file '//trim(FName),&
                  & ' : found nuclide '//MyNuclideName//', but I cannot give it a number... Exiting!!'
                  call exit()
               endif
               !
               ! Read values
               !
               if (iImmersion == InAir) then
                  read(ALine(18:DefaultLength),*) (DCCEffectiveDose(iImmersion,iAge)%x(jNuclide),iAge=1,NDCCAges),&
                     & DCCEquivalentDoseInAir%x(jNuclide),DCCKermaRateInAir%x(jNuclide)
               else
                  read(ALine(18:DefaultLength),*) (DCCEffectiveDose(iImmersion,iAge)%x(jNuclide),iAge=1,NDCCAges)
               endif
            endif ! not a dash line
         enddo ! loop over nuclides

         10  continue

         close(ScratchFile)
      enddo ! loop over air and water
   end subroutine ReadDCCEffectiveDose



   subroutine ReadGroundDCCs()
      !
      ! Read ground DCCs
      !
      character(DefaultLength) :: FName,ALine
      integer :: iLine,jNuclide,iGroundDepth,iAge,iSkip
      character(10) :: MyNuclideName
      logical :: IsDashLine

      do iGroundDepth = 1,NGroundDepths
         FName = trim(GroundDepthName(iGroundDepth))//'Soil_Nuclide-specific_EffectiveDose.txt'
         FName = trim(DCCSoilContaminationPath)//trim(FName)
         open(ScratchFile,file=FName,form='FORMATTED',action='READ')
         !
         ! Skip header 10 lines
         !
         do iSkip = 1,10
            read(ScratchFile,*)
         enddo
         iLine = 10
         !
         ! Read data for all nuclides
         !
         do
            read(ScratchFile,'(A)',end=10) ALine
            iLine = iLine + 1

            IsDashLine = (INDEX(ALine,'-----') > 0)

            if (.NOT.IsDashLine) then
               !
               ! Get nuclide and check name
               !
               MyNuclideName = ALine(1:7)
               jNuclide = GetNuclideNumber(MyNuclideName)
               if (jNuclide == 0) then
                  write(*,'(A,I0,A,A)') 'At line ',iLine,' in file '//trim(FName),&
                  & ' : found nuclide '//MyNuclideName//', but I cannot give it a number... Exiting!!'
                  call exit()
               endif
               !
               ! Read values
               !
               read(ALine(9:DefaultLength),*) (DCCGround(iGroundDepth,iAge)%x(jNuclide),iAge=1,NDCCAges),&
               & DCCEquivalentDoseGround(iGroundDepth)%x(jNuclide),DCCKermaRateGround(iGroundDepth)%x(jNuclide)
               !
               ! check if it is a noble gas and if they are supposed to have groundshine
               !
               if (isGas(NuclideSpecs(jNuclide)) .and. (.not. OPTION_NobleGroundShine)) then
                  DCCGround(:,:)%x(jNuclide) = 0._Float
                  DCCEquivalentDoseGround(:)%x(jNuclide) = 0._Float
                  DCCKermaRateGround(:)%x(jNuclide) = 0._Float
               endif

            endif ! not a dash line
         enddo ! loop over nuclides

         10  continue

         close(ScratchFile)
      enddo ! loop over in air and in water
   end subroutine ReadGroundDCCs



   subroutine ReadInhalationDCCs()
      !
      ! Read inhalation DCCs
      !
      character(DefaultLength) :: FName,ALine
      integer :: iLine,iNuclide,iInhalationDCC
      character(10) :: MyNuclideName

      FName = 'CoVeGa_v63_voorRapport_Inhalation.prn'
      FName = trim(DCCInhalationPath)//trim(FName)
      open(ScratchFile,file=FName,form='FORMATTED',action='READ')
      !
      ! Skip header
      !
      read(ScratchFile,*)
      !
      ! Read data for all nuclides
      !
      do iLine = 2,1253
         read(ScratchFile,'(A)') ALine
         !
         ! Get nuclide and check name
         !
         read(ALine,*) MyNuclideName
         iNuclide = GetNuclideNumber(MyNuclideName)
         if (iNuclide == 0) then
            write(*,'(A,I0,A,A)') 'At line ',iLine,' in file '//trim(FName),&
            & ' : found nuclide '//MyNuclideName//', but I cannot give it a number... Exiting!!'
            call exit()
         endif
         !
         ! Read values in Sv/Bq
         !
         read(ALine(58:DefaultLength),*) (DCCInhalation(iInhalationDCC)%x(iNuclide),iInhalationDCC=1,NInhalationDCCs)

      enddo ! loop over lines in file

      close(ScratchFile)
   end subroutine ReadInhalationDCCs



   subroutine ReadThyroidInhalationDCCs()
      !
      ! Read inhalation DCCs
      !
      character(DefaultLength) :: FName,ALine
      integer :: iLine,iNuclide
      character(10) :: MyNuclideName
      real(Float) :: Dum

      FName = RIVMSourcesPath() // '/dccunix.dat'
      open(ScratchFile,file=FName,form='FORMATTED',action='READ')
      !
      ! Skip header
      !
      read(ScratchFile,*)
      !
      ! Read data for all nuclides
      !
      do iLine = 2,247
         read(ScratchFile,'(A)') ALine
         !
         ! Get nuclide and check name
         !
         read(ALine,*) MyNuclideName
         iNuclide = GetNuclideNumber(MyNuclideName)
         if (iNuclide == 0) then
            write(*,'(A,I0,A,A)') 'At line ',iLine,' in file '//trim(FName),&
            & ' : found nuclide '//MyNuclideName//', but I cannot give it a number... Skipping nuclide!!'
         else
            !
            ! Read values in Sv/Bq
            !
            read(ALine(9:DefaultLength),*) Dum,DCCThyroidInhalation%x(iNuclide)
         endif

      enddo ! loop over lines in file

      close(ScratchFile)
   end subroutine ReadThyroidInhalationDCCs
end module LibDCC
