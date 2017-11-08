gap> START_TEST("d8.tst");
gap> testdatadir := DirectoriesPackageLibrary( "cohomolo", "testdata" );;
gap> CHMLINFO:=InfoLevel(InfoCohomolo);;
gap> SetInfoLevel(InfoCohomolo,1);

#
gap> Read( Filename ( testdatadir, "a4" ) );

#
gap> chr:= CHR( G, 2, F, m2 );;
gap> M:= SchurMultiplier( chr );
#I  Indices in the subgroup chain are:  3 
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
[ 2 ]
gap> D:= CoveringGroup( chr );;
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
gap> Size(D);
24
gap> FirstCohomologyDimension( chr );
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
2
gap> SecondCohomologyDimension( chr );
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
2

#
gap> chr:= CHR( G, 3, F, m3 );;
gap> M:= SchurMultiplier( chr );
The Sylow p-subgroup of the group is cyclic - so the multiplier is trivial.
[  ]
gap> D:= CoveringGroup( chr );;
gap> Size(D);
12
gap> FirstCohomologyDimension( chr );
#I  Indices in the subgroup chain are:  4 
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
2
gap> SecondCohomologyDimension( chr );
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
2

#
gap> SetInfoLevel(InfoCohomolo,CHMLINFO);
gap> STOP_TEST( "a4.tst", 10000 );
