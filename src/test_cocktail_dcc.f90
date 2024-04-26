PROGRAM test_cocktail_dcc
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
   USE libxmath
   USE libutil
   USE libcocktaildcc

   IMPLICIT NONE

   REAL(Float), PARAMETER :: GrowthFactor = 1.1_Float
   REAL(Float) :: t,MyEstimate,MyCumulativeEstimate
   INTEGER :: iTime,iSourceTerm,iPathway,iType,iNuclide
   REAL(Float), DIMENSION(10,3,2) :: MyDCC
   REAL(Float), DIMENSION(10) :: MyAmountI133,MyAmountKr85m,MyAmountCs137

   INTEGER, PARAMETER :: SourceTermUKrausFoster2 = 1
   INTEGER, PARAMETER :: SourceTermUKrausFoster5 = 2
   INTEGER, PARAMETER :: SourceTermPuSTUK        = 3
   INTEGER, PARAMETER :: SourceTermUSTUK         = 4
   INTEGER, PARAMETER :: SourceTermSTCCON1       = 5
   !
   ! Specify the location of the DCC files
   !
   RIVMSourcesPath = './../brontermen/'
   !
   ! Example 1: Show some specific regular and cumulative DCCs for specific times
   !
   WRITE(*,'(A)') 'Some specific cocktail DCCs and cumulative cocktail DCCs:'

   t = 24*3600._Float
   MyEstimate = GetCocktailDCC(t,SourceTermUKrausFoster2,PathwayAir,iRegularDCC)
   MyCumulativeEstimate = GetCocktailDCC(t,SourceTermUKrausFoster2,PathwayAir,iCumulativeDCC)
   WRITE(*,'(A,2EN15.5)') 'DCC(Uranium-KrausFoster,Air) after 1 day  = ',MyEstimate,MyCumulativeEstimate

   t = 7*24*3600._Float
   MyEstimate = GetCocktailDCC(t,SourceTermPuSTUK,PathwayGround,iRegularDCC)
   MyCumulativeEstimate = GetCocktailDCC(t,SourceTermPuSTUK,PathwayGround,iCumulativeDCC)
   WRITE(*,'(A,2EN15.5)') 'DCC(Plutonium-STUK,Ground)   after 1 week = ',MyEstimate,MyCumulativeEstimate

   t = 365.24_Float*24*3600._Float
   MyEstimate = GetCocktailDCC(t,SourceTermSTCCON1,PathwayInhalation,iRegularDCC)
   MyCumulativeEstimate = GetCocktailDCC(t,SourceTermSTCCON1,PathwayInhalation,iCumulativeDCC)
   WRITE(*,'(A,2EN15.5)') 'DCC(Uranium-RIVM,Inhalation) after 1 year = ',MyEstimate,MyCumulativeEstimate
   !
   ! Example 2: Estimate all types of DCC and cumulative DCC for all known types of SourceTerm
   ! for a set of times starting at delay 10 seconds and each next delay is 10% longer.
   !
   WRITE(*,*)
   WRITE(*,'(A)') 'For all types of SourceTerm and all pathways:'
   WRITE(*,'(15X,50A16)') (((SourceTermName(iSourceTerm)(1:12),iType=1,2),iPathway=1,NPathways),iSourceTerm=1,NSourceTerms)
   WRITE(*,'(13X,50A16)') (((PathwayName(iPathway),iType=1,2),iPathway=1,NPathways),iSourceTerm=1,NSourceTerms)
   WRITE(*,'(3X,A,9X,50A16)') 't[s]',(((DCCTypeName(iType),iType=1,2),iPathway=1,NPathways),iSourceTerm=1,NSourceTerms)

   t = 10._Float ! seconds

   DO iTime = 1,300

      DO iSourceTerm = 1,NSourceTerms
         DO iPathway = 1,NPathways
            DO iType = 1,2
               MyDCC(iSourceTerm,iPathway,iType) = GetCocktailDCC(t,iSourceTerm,iPathway,iType)
            ENDDO
         ENDDO
      ENDDO

      WRITE(*,'(50(EN15.5,1X))') t,(((MyDCC(iSourceTerm,iPathway,iType),&
      & iType=1,2),iPathway=1,NPathways),iSourceTerm=1,NSourceTerms)

      t = t*GrowthFactor
   ENDDO
   !
   ! Test 3: For each SourceTerm, some individual nuclides are followed in the first year:
   !
   WRITE(*,*)
   WRITE(*,'(A)') 'Some specific nuclides for all types of SourceTerm:'
   WRITE(*,'(15X,50A16)') ((SourceTermName(iSourceTerm)(1:12),iNuclide=1,3),iSourceTerm=1,NSourceTerms)
   WRITE(*,'(3X,A,9X,10A)') 't[day] ',('I-133           Kr-85m          Cs-137          ',iSourceTerm=1,NSourceTerms)

   DO iTime = 0,365

      t = iTime*24._Float*3600._Float ! seconds

      DO iSourceTerm = 1,NSourceTerms
         MyAmountI133( iSourceTerm) = GetCocktailNuclide(t,iSourceTerm,'I-133' )
         MyAmountKr85m(iSourceTerm) = GetCocktailNuclide(t,iSourceTerm,'Kr-85m')
         MyAmountCs137(iSourceTerm) = GetCocktailNuclide(t,iSourceTerm,'Cs137' ) ! The '-' in the name of the nuclide can be left out
      ENDDO

      WRITE(*,'(I3,12X,50(EN15.5,1X))') iTime,(MyAmountI133( iSourceTerm),&
      &                                       MyAmountKr85m(iSourceTerm),&
      &                                       MyAmountCs137(iSourceTerm),iSourceTerm=1,NSourceTerms)

   ENDDO

END PROGRAM test_cocktail_dcc