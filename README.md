# Explanation of files #


**testCocktailDCC02.f90**
Sample program showing how to use the functionality made available by module libcocktaildcc. Cocktail DCCs for different pathways and activities for 3 nuclides are estimated.

**libcocktaildcc.f90**
Library to support the use of cocktail DCCs. 

**SourcetermDose06AB.f90**
This utility reads a source term/nuclide vector at t=0 (blast) from file and returns lookup tables for (regular and cumulative) cocktail DCCs at the agreed set of pinpoints in time. It also generates lookup tables at these pinpoints for the total individual nuclide activities. You can construct cocktail DCCs for particular natures of either the mother nuclide or the daughter nuclide: "any", "noble gas" or "other than noble gas".

**Mature_Nuclides.f90**
This utility reads a nuclide vector at t=0 from file and a delay time interval. It returns the aged vector. This utility can be used to e.g. calculate the inventory of a nuclear reactor that has been shut down for some time.

**libpinpoint.f90**
Supporting library handling reading and saving the pinpoint tables and transition matrices.

**NuclideDecay03c.f90**
Utility. For all different delay pinpoints in time in 1 go, this utility constructs the decay matrix for (any combination of) all nuclides available in the ENDF database. The output is a set files with names SparseMatrix_ENDF_XXX.dat, where XXX is the number of the pinpoint, which runs from 0 to 200. These ready-made matrix files are available as a zip-file for all users in the repository of the cocktail DCC model. Alternatively, one can construct these matrix files oneselves with this utility NuclideDecay03c.f90. Please store these matrix files in a location that later in libpinpoint.f90 will be called “TransitionMatrixPath”.

**libendf.f90**
This library connects to the ENDF database of nuclides, which has a complex format. During the first invocation of this library, after reading the decay data from the ENDF database, the decay data are stored in an unformatted file. This is done to speed up future use of this library.

**libinterval.f90**
Supporting library for support of interval type variables and interpolation.

**libexponential.f90**
Supporting library for numerical estimation of the matrix exponential, based on the original subroutines from A. Miller as made available by J Blevins.

**libutil.f90**
Supporting library with general purpose material

**libxmath.f90**
Supporting library with extra mathematics not available in FORTRAN intrinsic functions.

**getdeps.f90**
Utility. For a given source file (module or main program), this utility returns the list of source files on which this source file is dependent. The list is returned as a text string, where the order of the fine names is put in order of dependency. The aim of the utility is to facilitate compilation of FORTRAN-90 code without a need for the user to specify the names of the source files, let alone order them.

# Supporting library #
handling reading and saving the pinpoint tables and transition matrices.

# Supporting files #
We use dataset of the ENDF-B-VIII, link: [ENDF/B-VIII.0 Evaluated Nuclear Data Library (bnl.gov)](https://www.nndc.bnl.gov/endf-b8.0/)

ICRP nuclides, link: http://www.icrp.org/publication.asp?id=ICRP%20Publication%20107

DCC from the ICRP-144, doi: https://doi.org/10.1177/0146645320906277

DCC from the ICRP-119, doi: http://dx.doi.org/10.1016/j.ympev.2012.04.018

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

The paper was made without the use of the binary file.

Several long-lived nuclides are currently considered as stable in the code because the ENDF-B database has no value in their numerical spontaneous decay description (MT = 457). these radionuclides are:
+	Cr-50   halflive: 1.3    e18	[year]
+	Zn-70   halflive: >= 3.8 e18	[year]
+	Se-82   halflive: 9.6    e19	[year]
+	Te-123  halflive: >9.2   e16	[year]
+	Te-130  halflive: 7.9    e20	[year]
+	Xe-124  halflive: 1.6    e14	[year]
+	Xe-134  halflive: 5.8    e22	[year]
+	Xe-136  halflive: 2.165  e21	[year]
+	Ba-132  halflive: 3.0    e21	[year]
+	Ce-138  halflive: 4.4    e16	[year]
+	Ce-142  halflive: 5      e16	[year]
+	Eu-151  halflive: 1.7    e18	[year]
+	Ta-180m halflive: 1.2    e15	[year]
+	W-183   halflive: 6.7    e20	[year]
+	Os-184  halflive: 5.6    e13	[year]

Er-145, a very exotic nuclide did not yet get a halflife, thus has a value that classifies it as stable. The code now gives it a half-live time of 1.0e-6 [s].

# References #
Dillen T van, Dijk A van, Kloosterman A, Russo F, Mommaert C. Accounting for ingrowth of radioactive progeny in dose assessments: generic weighting factors for dose coefficients. J Radiol Prot 40:83; 2019.

