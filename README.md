# Explanation of files #


**Libendf.f90**

  This library connects to the ENDF database of nuclides. This database is quite complex to decipher. Therefore, the resulting datastructures are stored in a binary file that can be found after it is run for the first time. 

**NuclideDecay03b.f90**

For all different delay pinpoints in time in 1 go, this utility constructs the decay matrices for (a vector of) all nuclides. If the ENDF data are already available in a pre-digested form, i.e. in a binary file that has been generated by this same very library, then this binary file is used, otherwise, the long road is taken and the original ASCII ENDF data are parsed. The output is a set of matrices with filenames SparseMatrix_ENDF_XXX.dat, where XXX is the number of the pinpoint, which runs from 0 to 200. Store these matrices in a location that later in libpinpoint.f90 will be called “TransitionMatrixPath”.
The ready-made files are available upon request, with the authors of the repository.
	
**Mature_nuclides.f90**

This utility reads a nuclide vector at t=0 from file and a delay. It returns the aged vector.

**Sourcetermdose06AB.f90**

This utility reads a nuclide vector at t=0 (blast) from file and returns lookup tables for (regular and cumulative) dose rates at a set of pinpoints in time and lookup tables for the total activities at these pinpoints. You can make estimates for particular natures of either mother or daughter: "any", "noble gas" or "other than noble gas".

**Getdeps.f90**

This utility creates a string with the needed dependencies in the correct order. Easing the build of different utilities.

**testCocktailDCC02.f90**

Sample program showing how to use the functionality made available by module libcocktaildcc. Cocktail DCCs for different pathways and activities for 3 nuclides are estimated.

**Libxmath**

Supporting library with calculation options
  
**Libutil**

Supporting library with useful function
  
**LibInterval**

Supporting library handling the interval calculations

**LibExponential**

Supporting library originally from J Blevins, handling matrix exponential calculations

**Libpinpoint**

Supporting library handling reading and saving the pinpoint tables and transitionmatrices.

# Supporting library #
handling reading and saving the pinpoint tables and transitionmatrices.

# Supporting files #
We use dataset of the ENDF-B-VIII, link: [ENDF/B-VIII.0 Evaluated Nuclear Data Library (bnl.gov)](https://www.nndc.bnl.gov/endf-b8.0/)

ICRP nuclides, link: http://www.icrp.org/publication.asp?id=ICRP%20Publication%20107

DCC from the ICRP-144, doi: https://doi.org/10.1177/0146645320906277

Supplemented with data from: ICRP, supplemental files, 2020. 

# How to build the utilities #

These are the commands that we used but might be different on another system.

+	gfortran -o getdeps.exe libxmath.f90 libutil.f90 getdeps.f90
+	gfortran -Ofast -o NuclideDecay03c.exe `./getdeps.exe NuclideDecay03c.f90`
+	gfortran -o Sourcetermdose06AC.exe `./getdeps.exe sourcetermdose06AC.f90`
+	gfortran -o testcocktailDCC02.exe `./getdeps.exe testcocktailDCC02.f90`
+	gfortran -o mature_nuclides.exe `./getdeps.exe mature_nuclides.f90`

# Remarks #

The NuclideDecay gives an overview of decay chains, within these chains the order of Mother nuclide to stable is top to bottom in case of beta- decay. And bottom to top in case of beta+ decay. 

When making different Pinpoint tables with the cocktails, no orphans (Dillen et al. 2019, appendix B) are allowed. Any radioisotope that has a halflife time of less than 10 seconds is considered an orphan and is not considered an appropriate head of chain. NOTE: In case you have the computing power to allow for more calculations, you could set the tmin (in both libdcc.f90 and nuclideDecay03c.f90) to lower values.

The paper was made without the use of the binary file, which currently misses the orphan nuclides. 

# References #
Dillen T van, Dijk A van, Kloosterman A, Russo F, Mommaert C. Accounting for ingrowth of radioactive progeny in dose assessments: generic weighting factors for dose coefficients. J Radiol Prot 40:83; 2019.
