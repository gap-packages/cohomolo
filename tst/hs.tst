gap> START_TEST("hs.tst");
gap> testdatadir := DirectoriesPackageLibrary( "cohomolo", "testdata" );;
gap> CHMLINFO:=InfoLevel(InfoCohomolo);;
gap> SetInfoLevel(InfoCohomolo,1);

#
gap> Read( Filename ( testdatadir, "hs" ) );

#
gap> chr:= CHR( G, 2 );;
gap> M:= SchurMultiplier( chr );
#I  Indices in the subgroup chain are:  4125 7 3 
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
[ 2 ]

#
gap> SetInfoLevel(InfoCohomolo,CHMLINFO);
gap> STOP_TEST( "hs.tst", 10000 );
