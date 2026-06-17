program test_cocktail_dcc
   !
   ! Developed For:
   !
   ! VLH   Centre for Environmental Safety and Security
   ! RIVM - National Institute for Public Health and the Environment
   !
   ! PO Box 1,
   ! NL - 3720 BA Bilthoven
   ! The Netherlands
   !
   ! Sample program showing how to use the functionality made available by module libcocktaildcc
   !
   use libcocktaildcc, only: GetCocktailDCC, GetCocktailNuclide, &
      PathwayAir, PathwayInhalation, PathwayName, &
      iRegularDCC,iCumulativeDCC, DCCTypeName, &
      NSourceTerms, NPathways, &
      SourceTermName
   use libxmath, only: Float

   implicit none

   real(Float), parameter :: GrowthFactor = 1.1_Float
   real(Float) :: t,MyEstimate,MyCumulativeEstimate
   integer :: iTime,iSourceTerm,iPathway,iType,iNuclide
   real(Float), dimension(10,3,2) :: MyDCC
   real(Float), dimension(10) :: MyAmountI133,MyAmountKr85m,MyAmountCs137

   integer, parameter :: SourceTermUKrausFoster2 = 1
   integer, parameter :: SourceTermAxelsson = 2
   !
   ! Example 1: Show some specific regular and cumulative DCCs for specific times
   !
   write(*,'(A)') 'Some specific cocktail DCCs and cumulative cocktail DCCs:'

   t = 24*3600._Float
   MyEstimate = GetCocktailDCC(t,SourceTermUKrausFoster2,PathwayAir,iRegularDCC)
   MyCumulativeEstimate = GetCocktailDCC(t,SourceTermUKrausFoster2,PathwayAir,iCumulativeDCC)
   write(*,'(A,2EN15.5)') 'DCC(Uranium-KrausFoster,Air) after 1 day  = ',MyEstimate,MyCumulativeEstimate

   !t = 7*24*3600._Float
   !MyEstimate = GetCocktailDCC(t,SourceTermPuSTUK,PathwayGround,iRegularDCC)
   !MyCumulativeEstimate = GetCocktailDCC(t,SourceTermPuSTUK,PathwayGround,iCumulativeDCC)

   t = 365.24_Float*24*3600._Float
   MyEstimate = GetCocktailDCC(t,SourceTermAxelsson,PathwayInhalation,iRegularDCC)
   MyCumulativeEstimate = GetCocktailDCC(t,SourceTermAxelsson,PathwayInhalation,iCumulativeDCC)
   write(*,'(A,2EN15.5)') 'DCC(Uranium-RIVM,Inhalation) after 1 year = ',MyEstimate,MyCumulativeEstimate
   !
   ! Example 2: Estimate all types of DCC and cumulative DCC for all known types of SourceTerm
   ! for a set of times starting at delay 10 seconds and each next delay is 10% longer.
   !
   write(*,*)
   write(*,'(A)') 'For all types of SourceTerm and all pathways:'
   write(*,'(15X,50A16)') (((SourceTermName(iSourceTerm)(1:12),iType=1,2),iPathway=1,NPathways),iSourceTerm=1,NSourceTerms)
   write(*,'(13X,50A16)') (((PathwayName(iPathway),iType=1,2),iPathway=1,NPathways),iSourceTerm=1,NSourceTerms)
   write(*,'(3X,A,9X,50A16)') 't[s]',(((DCCTypeName(iType),iType=1,2),iPathway=1,NPathways),iSourceTerm=1,NSourceTerms)

   t = 10._Float ! seconds

   do iTime = 1,300

      do iSourceTerm = 1,NSourceTerms
         do iPathway = 1,NPathways
            do iType = 1,2
               MyDCC(iSourceTerm,iPathway,iType) = GetCocktailDCC(t,iSourceTerm,iPathway,iType)
            enddo
         enddo
      enddo

      write(*,'(50(EN15.5,1X))') t,(((MyDCC(iSourceTerm,iPathway,iType),&
      & iType=1,2),iPathway=1,NPathways),iSourceTerm=1,NSourceTerms)

      t = t*GrowthFactor
   enddo
   !
   ! Test 3: For each SourceTerm, some individual nuclides are followed in the first year:
   !
   write(*,*)
   write(*,'(A)') 'Some specific nuclides for all types of SourceTerm:'
   write(*,'(15X,50A16)') ((SourceTermName(iSourceTerm)(1:12),iNuclide=1,3),iSourceTerm=1,NSourceTerms)
   write(*,'(3X,A,9X,10A)') 't[day] ',('I-133           Kr-85m          Cs-137          ',iSourceTerm=1,NSourceTerms)

   do iTime = 0,365

      t = iTime*24._Float*3600._Float ! seconds

      do iSourceTerm = 1,NSourceTerms
         MyAmountI133( iSourceTerm) = GetCocktailNuclide(t,iSourceTerm,'I-133' )
         MyAmountKr85m(iSourceTerm) = GetCocktailNuclide(t,iSourceTerm,'Kr-85m')
         MyAmountCs137(iSourceTerm) = GetCocktailNuclide(t,iSourceTerm,'Cs137' ) ! The '-' in the name of the nuclide can be left out
      enddo

      write(*,'(I3,12X,50(EN15.5,1X))') iTime,(MyAmountI133( iSourceTerm),&
      &                                       MyAmountKr85m(iSourceTerm),&
      &                                       MyAmountCs137(iSourceTerm),iSourceTerm=1,NSourceTerms)

   enddo

end program test_cocktail_dcc