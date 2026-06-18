program source_term_dose
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
   use LibPinpoint, only: InitLibPinPoint, ReadRIVMSourceTermFile, CocktailType, &
      GetPinpointCocktails, MakePinpointDoseRates
   use LibUtil, only: AllLowCase, DefaultLength

   implicit none

   type(CocktailType) :: StartCocktail
   character(DefaultLength) :: FName,WithDaughterString,MotherNatureString,DaughterNatureString
   integer :: NArguements,MotherNature,DaughterNature
   logical :: WithDaughters,UseICRP

   logical, parameter :: OPTION_NobleGroundShine = .True.
   !
   ! Get commmand line arguement:
   !
   NArguements = IARGC()

   if (NArguements /= 4) then
      write(*,'(A)') 'Call:'
      write(*,*)
      write(*,'(A)') 'source_term_dose.exe <fname> <withdaughters> <Mother nature> <daughter nature>'
      write(*,*)
      write(*,'(A)') 'where:'
      write(*,'(A)') '<fname> is the name of a file with the nuclide vector at t=0'
      write(*,'(A)') '<withdaughers> is "yes" or "no"'
      write(*,'(A)') '<mother nature> is 0 for "any", 1 for "noble gases", 2 for "other than noble gases"'
      write(*,'(A)') '<daughter nature> is 0 for "any", 1 for "noble gases", 2 for "other than noble gases"'
      call exit()
   endif

   call getarg(1,FName)
   call getarg(2,WithDaughterString)

   call AllLowCase(WithDaughterString)

   WithDaughters = (INDEX(WithDaughterString,'yes') >0)
   if (WithDaughters) then
      write(*,'(A)') 'Going to include full ingrowth of progeny!'
   else
      write(*,'(A)') 'Assuming clean decay: any decay makes the substance vanish...'
   endif
   write(*,*)

   call getarg(3,MotherNatureString)
   read(MotherNatureString,*) MotherNature
   call getarg(4,DaughterNatureString)
   read(DaughterNatureString,*) DaughterNature
   !
   ! Initialize
   !
   UseICRP = .FALSE.
   call InitLibPinpoint(UseICRP)
   !
   ! Specify source term at t=0:
   !
   call ReadRIVMSourceTermFile(FName,StartCocktail)

   write(*,'(A)') 'Constructing lookup tables for sourceterm "'//trim(StartCocktail%MyName)&
   & //'" in directory "'//trim(StartCocktail%MyDirectory)//'"'
   !
   ! Construct dose-rate contributions at all pinpoints in time, including thinning scenario
   !
   call MakePinpointDoseRates(StartCocktail,WithDaughters,MotherNature,DaughterNature)
   !
   ! Construct cocktails at all pinpoints in time
   !
   call GetPinpointCocktails(StartCocktail,WithDaughters)
end program source_term_dose