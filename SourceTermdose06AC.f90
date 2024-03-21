PROGRAM SourceTermDose
!
! This utility reads a nuclide vector at t=0 (blast) from file
! and returns lookup tables for (regular and cumulative) dose rates at a set of pinpoints in time
! and lookup tables for the total activitues at these pinpoints.
! You can make estimates for particular natures of either mother or daughter: "any", "noble gas" or "other than noble gas"
!
! Groundshine DCCs: ICRP-144 (supplementary information)
!
! AC difference with version AB: only ENDF data are used for decay, ICRP data only used for DCC values.
!
   USE libxmath
   USE libutil
   USE libendf
   USE libpinpoint
   USE libdcc

   IMPLICIT NONE

   TYPE(CocktailType) :: StartCocktail
   CHARACTER(DefaultLength) :: FName,TheDirectory,TheFile,WithDaughterString,MotherNatureString,DaughterNatureString
   REAL(Float), DIMENSION(0:NStartingTimes,MaxNuclides) :: PinAirDoseRate,PinGroundDoseRate,PinInhalationDoseRate
   INTEGER :: NArguements,MotherNature,DaughterNature
   LOGICAL :: WithDaughters,UseICRP
   !
   ! Get commmand line arguement:
   !
   NArguements = IARGC()

   IF (NArguements.NE.4) THEN
      WRITE(*,'(A)') 'Call:'
      WRITE(*,*)
      WRITE(*,'(A)') 'SourceTermdose06AB.exe <fname> <withdaughters> <Mother nature> <daughter nature>'
      WRITE(*,*)
      WRITE(*,'(A)') 'where:'
      WRITE(*,'(A)') '<fname> is the name of a file with the nuclide vector at t=0'
      WRITE(*,'(A)') '<withdaughers> is "yes" or "no"'
      WRITE(*,'(A)') '<mother nature> is 0 for "any", 1 for "noble gases", 2 for "other than noble gases"'
      WRITE(*,'(A)') '<daughter nature> is 0 for "any", 1 for "noble gases", 2 for "other than noble gases"'
      CALL EXIT()
   ENDIF

   CALL GETARG(1,FName)
   CALL GETARG(2,WithDaughterString)

   CALL AllLowCase(WithDaughterString)

   WithDaughters = (INDEX(WithDaughterString,'yes').GT.0)
   IF (WithDaughters) THEN
      WRITE(*,'(A)') 'Going to include full ingrowth of progeny!'
   ELSE
      WRITE(*,'(A)') 'Assuming clean decay: any decay makes the substance vanish...'
   ENDIF
   WRITE(*,*)

   CALL GETARG(3,MotherNatureString)
   READ(MotherNatureString,*) MotherNature
   CALL GETARG(4,DaughterNatureString)
   READ(DaughterNatureString,*) DaughterNature
   !
   ! Initialize
   !
   UseICRP = .FALSE.
   CALL InitLibPinpoint(UseICRP)
   !
   ! Specify source term at t=0:
   !
   CALL ReadRIVMSourceTermFile(FName,StartCocktail)

   WRITE(*,'(A)') 'Constructing lookup tables for sourceterm "'//TRIM(StartCocktail%MyName)&
   & //'" in directory "'//TRIM(StartCocktail%MyDirectory)//'"'
   !
   ! Construct dose-rate contributions at all pinpoints in time, including thinning scenario
   !
   CALL MakePinpointDoseRates(StartCocktail,WithDaughters,MotherNature,DaughterNature)
   !
   ! Construct cocktails at all pinpoints in time
   !
   CALL GetPinpointCocktails(StartCocktail,WithDaughters)
END PROGRAM SourceTermDose