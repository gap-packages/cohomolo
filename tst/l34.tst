gap> START_TEST("l34.tst");
gap> testdatadir := DirectoriesPackageLibrary( "cohomolo", "testdata" );;
gap> CHMLINFO:=InfoLevel(InfoCohomolo);;
gap> SetInfoLevel(InfoCohomolo,1);

#
gap> Read( Filename ( testdatadir, "l34" ) );

#
gap> chr:= CHR( G, 2 );;
gap> M:= SchurMultiplier( chr );
#I  Indices in the subgroup chain are:  21 5 3 
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
[ 4, 4 ]

#
gap> chr:= CHR( G, 3 );;
gap> M:= SchurMultiplier( chr );
#I  Indices in the subgroup chain are:  280 8 
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
[ 3 ]

#
gap> SetInfoLevel(InfoCohomolo,CHMLINFO);
gap> STOP_TEST( "l34.tst", 10000 );
