MODULE LibDCC
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
   USE libxmath
   USE libutil
   USE libendf

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: InitLibDCC,ReadTissueDCCs,ReadDCCEffectiveDose,ReadGroundDCCs,DCCEffectiveDose,DCCGround,&
   & InAir,Age_Adult,DCCInhalation,ReadInhalationDCCs,Inhalation_public_adult,&
   & DCCThyroidInhalation,ReadThyroidInhalationDCCs

   CHARACTER(DefaultLength) :: ProjectPath = './'

   CHARACTER(DefaultLength) :: DCCCalcPath,DCCAirSubmersionPath,DCCSoilContaminationPath,DCCWaterImmersionPath,&
   & DCCInhalationPath

   INTEGER, PARAMETER :: InAir = 1
   INTEGER, PARAMETER :: InWater = 2

   INTEGER, PARAMETER :: Tissue_R_marrow   =  1
   INTEGER, PARAMETER :: Tissue_Colon      =  2
   INTEGER, PARAMETER :: Tissue_Lungs      =  3
   INTEGER, PARAMETER :: Tissue_ST_wall    =  4
   INTEGER, PARAMETER :: Tissue_Breast     =  5
   INTEGER, PARAMETER :: Tissue_Ovaries    =  6
   INTEGER, PARAMETER :: Tissue_Testes     =  7
   INTEGER, PARAMETER :: Tissue_UB_wall    =  8
   INTEGER, PARAMETER :: Tissue_Oesophagus =  9
   INTEGER, PARAMETER :: Tissue_Liver      = 10
   INTEGER, PARAMETER :: Tissue_Thyroid    = 11
   INTEGER, PARAMETER :: Tissue_Endost_BS  = 12
   INTEGER, PARAMETER :: Tissue_Brain      = 13
   INTEGER, PARAMETER :: Tissue_S_glands   = 14
   INTEGER, PARAMETER :: Tissue_Skin       = 15
   INTEGER, PARAMETER :: Tissue_Adrenals   = 16
   INTEGER, PARAMETER :: Tissue_ET         = 17
   INTEGER, PARAMETER :: Tissue_GB_wall    = 18
   INTEGER, PARAMETER :: Tissue_Ht_wall    = 19
   INTEGER, PARAMETER :: Tissue_Kidneys    = 20
   INTEGER, PARAMETER :: Tissue_Lymph      = 21
   INTEGER, PARAMETER :: Tissue_Muscle     = 22
   INTEGER, PARAMETER :: Tissue_O_mucosa   = 23
   INTEGER, PARAMETER :: Tissue_Pancreas   = 24
   INTEGER, PARAMETER :: Tissue_Prostate   = 25
   INTEGER, PARAMETER :: Tissue_SI_wall    = 26
   INTEGER, PARAMETER :: Tissue_Spleen     = 27
   INTEGER, PARAMETER :: Tissue_Thymus     = 28
   INTEGER, PARAMETER :: Tissue_Uterus     = 29

   INTEGER, PARAMETER :: NTissues          = 29

   CHARACTER(10), PARAMETER, DIMENSION(NTissues) :: TissueName = &
      & (/'R-marrow  ',&
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
      &   'Uterus    '/)

   INTEGER, PARAMETER :: Age_Adult   = 1
   INTEGER, PARAMETER :: Age_15yr    = 2
   INTEGER, PARAMETER :: Age_10yr    = 3
   INTEGER, PARAMETER :: Age_5yr     = 4
   INTEGER, PARAMETER :: Age_1yr     = 5
   INTEGER, PARAMETER :: Age_Newborn = 6

   INTEGER, PARAMETER :: NDCCAges = 6

   CHARACTER(7), PARAMETER, DIMENSION(NDCCAges) :: DCCAgeName = &
      & (/'Adult  ',&
      &   '15yr   ',&
      &   '10yr   ',&
      &   '5yr    ',&
      &   '1yr    ',&
      &   'Newborn'/)

   INTEGER, PARAMETER :: NGroundDepths = 9

   CHARACTER(45), PARAMETER, DIMENSION(NGroundDepths) :: GroundDepthName = &
      & (/'Planar sources at specific depths/0.5 g.cm-2/',&
      &   'Exponential sources/0.5 g.cm-2/              ',&
      &   'Exponential sources/1.0 g.cm-2/              ',&
      &   'Exponential sources/2.5 g.cm-2/              ',&
      &   'Exponential sources/5.0 g.cm-2/              ',&
      &   'Exponential sources/10.0 g.cm-2/             ',&
      &   'Exponential sources/20.0 g.cm-2/             ',&
      &   'Exponential sources/50.0 g.cm-2/             ',&
      &   'Exponential sources/100.0 g.cm-2/            '/)

   INTEGER, PARAMETER :: NInhalationDCCs = 8

   INTEGER, PARAMETER :: Inhalation_worker_1um     = 1
   INTEGER, PARAMETER :: Inhalation_worker_5um     = 2
   INTEGER, PARAMETER :: Inhalation_public_infant  = 3
   INTEGER, PARAMETER :: Inhalation_public_1_year  = 4
   INTEGER, PARAMETER :: Inhalation_public_5_year  = 5
   INTEGER, PARAMETER :: Inhalation_public_10_year = 6
   INTEGER, PARAMETER :: Inhalation_public_15_year = 7
   INTEGER, PARAMETER :: Inhalation_public_adult   = 8

   CHARACTER(14), PARAMETER, DIMENSION(NInhalationDCCs) :: InhalationDCCName = &
      & (/'worker 1um    ',&
      &   'worker 5um    ',&
      &   'public infant ',&
      &   'public 1 year ',&
      &   'public 5 year ',&
      &   'public 10 year',&
      &   'public 15 year',&
      &   'public adult  '/)

   TYPE DCCVector
      REAL(Float), DIMENSION(MaxNuclides) :: x
   END TYPE DCCVector

   TYPE(DCCVector), DIMENSION(InAir:InWater,Male:Female,NTissues,NDCCAges) :: DCCTissue
   TYPE(DCCVector), DIMENSION(NGroundDepths,NDCCAges) :: DCCGround
   TYPE(DCCVector), DIMENSION(NGroundDepths) :: DCCEquivalentDoseGround,DCCKermaRateGround
   TYPE(DCCVector), DIMENSION(InAir:InWater,NDCCAges) :: DCCEffectiveDose
   TYPE(DCCVector) :: DCCEquivalentDoseInAir,DCCKermaRateInAir,DCCThyroidInhalation
   TYPE(DCCVector), DIMENSION(NInhalationDCCs) :: DCCInhalation

   REAL(Float) :: BreathingRate = 1.2_Float ! [m3/h] Value can be overruled later

   TYPE(NuclideType), DIMENSION(0:MaxNuclides) :: RegularizedNuclideSpecs
   REAL(Float), DIMENSION(MaxNuclides,MaxNuclides) :: RegularizedMotherDaughterMatrix

CONTAINS

   SUBROUTINE InitLibDCC(UseICRP)
      !
      ! Initialize this library
      !
      LOGICAL, INTENT(IN) :: UseICRP
      REAL(Float) :: tMin

      INTEGER, PARAMETER :: DebugLevel = 0

      IF (UseICRP) THEN
         CALL ReadNuclideSpecs('ICRP-07.NDX')
      ELSE
         !
         ! Set shortest halflife that is admitted to the numerical scheme for solving the Bateman equations.
         ! Shorter ones are skipped: the grandmother is directly transformed into the grand-daughter.
         ! The missing decay staps are added afterwards.
         !
         tMin = 10._Float ! [s]

         IF (DebugLevel.GT.0) THEN
            WRITE(*,'(A)') 'Set shortest halflife admitted to the numerical scheme for solving the Bateman equations.'
            WRITE(*,'(A)') 'Shorter ones are skipped: the grandmother is directly transformed into the grand-daughter.'
            WRITE(*,'(A)') 'The missing decay steps are added afterwards.'
            WRITE(*,'(A,F15.5,A)') 'For this run, we choose tMin = ',tMin,' [s]'
         ENDIF

         CALL ReadNProcessENDFNuclideSpecs(tMin,RegularizedNuclideSpecs,RegularizedMotherDaughterMatrix)
      ENDIF

      DCCCalcPath = TRIM(ProjectPath)//'build/'
      DCCAirSubmersionPath = TRIM(DCCCalcPath)//'external/Supplemental Files -v4/Air submersion/'
      DCCSoilContaminationPath = TRIM(DCCCalcPath)//'external/Supplemental Files -v4/Soil contamination/'
      DCCWaterImmersionPath = TRIM(DCCCalcPath)//'external/Supplemental Files -v4/Water immersion/'
      DCCInhalationPath = TRIM(ProjectPath)//'resources/'
   END SUBROUTINE InitLibDCC



   SUBROUTINE ReadTissueDCCs()
      !
      ! Read tissue- and gender specific DCCs
      ! Attention: This file has only the ICRP nuclides, but the reference in this model may be ENDF,
      ! which has ~ 3 x more nuclides!
      !
      CHARACTER(DefaultLength) :: FName,ALine
      INTEGER :: iGender,iLine,jNuclide,iTissue,iAge,iImmersion,NDashLinesFound
      CHARACTER(10) :: MyTissueName,MyNuclideName
      LOGICAL :: IsDashLine

      DO iImmersion = InAir,InWater
         DO iGender = Male,Female
            FName = 'Nuclide-specific_EquivalentDose_'//TRIM(GenderName(iGender))//'.txt'
            IF (iImmersion.EQ.InAir) THEN
               FName = 'Air_'//TRIM(FName)
               FName = TRIM(DCCAirSubmersionPath)//TRIM(FName)
            ELSE
               FName = 'Water_'//TRIM(FName)
               FName = TRIM(DCCWaterImmersionPath)//TRIM(FName)
            ENDIF
            OPEN(ScratchFile,FILE=FName,FORM='FORMATTED',ACTION='READ')
            !
            ! Skip header
            !
            NDashLinesFound = 0
            iLine = 0

            DO WHILE (NDashLinesFound.LT.2)
               READ(ScratchFile,'(A)') ALine
               iLine = iLine + 1
               IF (ALine(1:1).EQ.'-') NDashLinesFound = NDashLinesFound + 1
            ENDDO
            !
            ! Read data for all nuclides
            !
            IsDashLine = .FALSE.

            DO WHILE (.NOT.IsDashLine)
               DO iTissue = 1,NTissues
                  READ(ScratchFile,'(A)',END=10) ALine
                  iLine = iLine + 1

                  IsDashLine = (INDEX(ALine,'-----').GT.0)

                  IF (.NOT.IsDashLine) THEN
                     !
                     ! Get nuclide and check name
                     !
                     IF (iTissue.EQ.1) THEN
                        MyNuclideName = ALine(1:7)
                        jNuclide = GetNuclideNumber(MyNuclideName)
                        IF (jNuclide.EQ.0) THEN
                           WRITE(*,'(A,I0,A,A)') 'At line ',iLine,' in file '//TRIM(FName),&
                           & ' : found nuclide '//MyNuclideName//', but I cannot give it a number... Exiting!!'
                           CALL EXIT()
                        ENDIF
                     ENDIF ! First organ
                     !
                     ! Get and check tissue name
                     !
                     MyTissueName = ALine(9:18)
                     IF (MyTissueName.NE.TissueName(iTissue)) THEN
                        WRITE(*,'(A,I0,A)') 'At line ',iLine,&
                        & ' : expected tissue '//TissueName(iTissue)//', but found '//MyTissueName//'! Exiting!!'
                        CALL EXIT()
                     ENDIF
                     !
                     ! Read values
                     !
                     READ(ALine(25:DefaultLength),*) (DCCTissue(iImmersion,iGender,iTissue,iAge)%x(jNuclide),iAge=1,NDCCAges)
                  ENDIF ! this is not a dashline
               ENDDO ! loop over tissues
               READ(ScratchFile,*) ! Skip line with dashes
               iLine = iLine + 1
            ENDDO ! loop over nuclides

            10    CONTINUE

            CLOSE(ScratchFile)
         ENDDO ! loop over gender
      ENDDO ! loop over in air and in water
   END SUBROUTINE ReadTissueDCCs



   SUBROUTINE ReadDCCEffectiveDose()
      !
      ! Read effective dose DCCs
      ! Attention: This file has only the ICRP nuclides, but the reference in this model may be ENDF,
      ! which has ~ 3 x more nuclides!
      !
      CHARACTER(DefaultLength) :: FName,ALine
      INTEGER :: iLine,jNuclide,iAge,NDashLinesFound,iImmersion
      CHARACTER(10) :: MyNuclideName
      LOGICAL :: IsDashLine

      INTEGER, PARAMETER :: DebugLevel = 0

      DO iImmersion = InAir,InWater
         FName = 'Nuclide-specific_EffectiveDose.txt'
         IF (iImmersion.EQ.InAir) THEN
            FName = 'Air_'//TRIM(FName)
            FName = TRIM(DCCAirSubmersionPath)//TRIM(FName)
         ELSE
            FName = 'Water_'//TRIM(FName)
            FName = TRIM(DCCWaterImmersionPath)//TRIM(FName)
         ENDIF

         OPEN(ScratchFile,FILE=FName,FORM='FORMATTED',ACTION='READ')
         !
         ! Skip header
         !
         NDashLinesFound = 0
         iLine = 0

         DO WHILE (NDashLinesFound.LT.2)
            READ(ScratchFile,'(A)') ALine
            iLine = iLine + 1
            IF (ALine(1:1).EQ.'-') NDashLinesFound = NDashLinesFound + 1
         ENDDO
         !
         ! Read data for all nuclides
         !
         DO
            READ(ScratchFile,'(A)',END=10) ALine
            iLine = iLine + 1

            IsDashLine = (INDEX(ALine,'-----').GT.0)

            IF (.NOT.IsDashLine) THEN
               !
               ! Get nuclide and check name
               !
               MyNuclideName = ALine(1:7)
               jNuclide = GetNuclideNumber(MyNuclideName)
               IF (jNuclide.EQ.0) THEN
                  WRITE(*,'(A,I0,A,A)') 'At line ',iLine,' in file '//TRIM(FName),&
                  & ' : found nuclide '//MyNuclideName//', but I cannot give it a number... Exiting!!'
                  CALL EXIT()
               ENDIF
               !
               ! Read values
               !
               IF (iImmersion.EQ.InAir) THEN
                  READ(ALine(18:DefaultLength),*) (DCCEffectiveDose(iImmersion,iAge)%x(jNuclide),iAge=1,NDCCAges),&
                     & DCCEquivalentDoseInAir%x(jNuclide),DCCKermaRateInAir%x(jNuclide)
               ELSE
                  READ(ALine(18:DefaultLength),*) (DCCEffectiveDose(iImmersion,iAge)%x(jNuclide),iAge=1,NDCCAges)
               ENDIF
            ENDIF ! not a dash line
         ENDDO ! loop over nuclides

         10  CONTINUE

         CLOSE(ScratchFile)
      ENDDO ! loop over air and water
   END SUBROUTINE ReadDCCEffectiveDose



   SUBROUTINE ReadGroundDCCs()
      !
      ! Read ground DCCs
      !
      CHARACTER(DefaultLength) :: FName,ALine
      INTEGER :: iLine,jNuclide,iGroundDepth,iAge,DumLine, iSkip
      CHARACTER(10) :: MyNuclideName
      LOGICAL :: IsDashLine

      DO iGroundDepth = 1,NGroundDepths
         FName = TRIM(GroundDepthName(iGroundDepth))//'Soil_Nuclide-specific_EffectiveDose.txt'
         FName = TRIM(DCCSoilContaminationPath)//TRIM(FName)
         OPEN(ScratchFile,FILE=FName,FORM='FORMATTED',ACTION='READ')
         !
         ! Skip header 10 lines
         !
         do iSkip = 1,10
            READ(ScratchFile,*)
         enddo
         iLine = 10
         !
         ! Read data for all nuclides
         !
         DO
            READ(ScratchFile,'(A)',END=10) ALine
            iLine = iLine + 1

            IsDashLine = (INDEX(ALine,'-----').GT.0)

            IF (.NOT.IsDashLine) THEN
               !
               ! Get nuclide and check name
               !
               MyNuclideName = ALine(1:7)
               jNuclide = GetNuclideNumber(MyNuclideName)
               IF (jNuclide.EQ.0) THEN
                  WRITE(*,'(A,I0,A,A)') 'At line ',iLine,' in file '//TRIM(FName),&
                  & ' : found nuclide '//MyNuclideName//', but I cannot give it a number... Exiting!!'
                  CALL EXIT()
               ENDIF
               !
               ! Read values
               !
               READ(ALine(9:DefaultLength),*) (DCCGround(iGroundDepth,iAge)%x(jNuclide),iAge=1,NDCCAges),&
               & DCCEquivalentDoseGround(iGroundDepth)%x(jNuclide),DCCKermaRateGround(iGroundDepth)%x(jNuclide)

            ENDIF ! not a dash line
         ENDDO ! loop over nuclides

         10  CONTINUE

         CLOSE(ScratchFile)
      ENDDO ! loop over in air and in water
   END SUBROUTINE ReadGroundDCCs



   SUBROUTINE ReadInhalationDCCs()
      !
      ! Read inhalation DCCs
      !
      CHARACTER(DefaultLength) :: FName,ALine
      INTEGER :: iLine,iNuclide,DumLine,iInhalationDCC
      CHARACTER(10) :: MyNuclideName

      FName = 'CoVeGa_v63_voorRapport_Inhalation.prn'
      FName = TRIM(DCCInhalationPath)//TRIM(FName)
      OPEN(ScratchFile,FILE=FName,FORM='FORMATTED',ACTION='READ')
      !
      ! Skip header
      !
      READ(ScratchFile,*)
      !
      ! Read data for all nuclides
      !
      DO iLine = 2,1253
         READ(ScratchFile,'(A)') ALine
         !
         ! Get nuclide and check name
         !
         READ(ALine,*) MyNuclideName
         iNuclide = GetNuclideNumber(MyNuclideName)
         IF (iNuclide.EQ.0) THEN
            WRITE(*,'(A,I0,A,A)') 'At line ',iLine,' in file '//TRIM(FName),&
            & ' : found nuclide '//MyNuclideName//', but I cannot give it a number... Exiting!!'
            CALL EXIT()
         ENDIF
         !
         ! Read values in Sv/Bq
         !
         READ(ALine(58:DefaultLength),*) (DCCInhalation(iInhalationDCC)%x(iNuclide),iInhalationDCC=1,NInhalationDCCs)

      ENDDO ! loop over lines in file

      CLOSE(ScratchFile)
   END SUBROUTINE ReadInhalationDCCs



   SUBROUTINE ReadThyroidInhalationDCCs()
      !
      ! Read inhalation DCCs
      !
      CHARACTER(DefaultLength) :: FName,ALine
      INTEGER :: iLine,iNuclide,DumLine
      CHARACTER(10) :: MyNuclideName
      REAL(Float) :: Dum

      FName = 'dccunix.dat'
      OPEN(ScratchFile,FILE=FName,FORM='FORMATTED',ACTION='READ')
      !
      ! Skip header
      !
      READ(ScratchFile,*)
      !
      ! Read data for all nuclides
      !
      DO iLine = 2,247
         READ(ScratchFile,'(A)') ALine
         !
         ! Get nuclide and check name
         !
         READ(ALine,*) MyNuclideName
         iNuclide = GetNuclideNumber(MyNuclideName)
         IF (iNuclide.EQ.0) THEN
            WRITE(*,'(A,I0,A,A)') 'At line ',iLine,' in file '//TRIM(FName),&
            & ' : found nuclide '//MyNuclideName//', but I cannot give it a number... Skipping nuclide!!'
         ELSE
            !
            ! Read values in Sv/Bq
            !
            READ(ALine(9:DefaultLength),*) Dum,DCCThyroidInhalation%x(iNuclide)
         ENDIF

      ENDDO ! loop over lines in file

      CLOSE(ScratchFile)
   END SUBROUTINE ReadThyroidInhalationDCCs
END MODULE LibDCC
