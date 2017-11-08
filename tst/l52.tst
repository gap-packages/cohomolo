gap> START_TEST("l52.tst");
gap> testdatadir := DirectoriesPackageLibrary( "cohomolo", "testdata" );;
gap> CHMLINFO:=InfoLevel(InfoCohomolo);;
gap> SetInfoLevel(InfoCohomolo,1);

#
gap> Read( Filename ( testdatadir, "l52" ) );

#
gap> chr:= CHR( G, 2, 0, m2 );;
gap> M:= SchurMultiplier( chr );
#I  Indices in the subgroup chain are:  31 15 7 3 
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
[  ]

#
gap> SetInfoLevel(InfoCohomolo,CHMLINFO);
gap> STOP_TEST( "l52.tst", 10000 );
