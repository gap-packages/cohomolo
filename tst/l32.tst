gap> START_TEST("l32.tst");
gap> testdatadir := DirectoriesPackageLibrary( "cohomolo", "testdata" );;
gap> CHMLINFO:=InfoLevel(InfoCohomolo);;
gap> SetInfoLevel(InfoCohomolo,1);

#
gap> Read( Filename ( testdatadir, "l32" ) );

#
gap> chr:= CHR( G, 2, F, m2 );;
gap> M:= SchurMultiplier( chr );
#I  Indices in the subgroup chain are:  7 3 
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
[ 2 ]
gap> D:= CoveringGroup( chr );;
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
gap> Size(D);
336
gap> FirstCohomologyDimension( chr );
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
1
gap> SecondCohomologyDimension( chr );
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
1

#
gap> chr:= CHR( G, 7, F, m7 );;
gap> M:= SchurMultiplier( chr );
The Sylow p-subgroup of the group is cyclic - so the multiplier is trivial.
[  ]
gap> D:= CoveringGroup( chr );;
gap> Size(D);
168
gap> FirstCohomologyDimension( chr );
#I  Indices in the subgroup chain are:  8 3 
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
0
gap> SecondCohomologyDimension( chr );
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
1

#
gap> SetInfoLevel(InfoCohomolo,CHMLINFO);
gap> STOP_TEST( "l32.tst", 10000 );
