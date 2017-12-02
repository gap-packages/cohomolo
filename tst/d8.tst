gap> START_TEST("d8.tst");
gap> testdatadir := DirectoriesPackageLibrary( "cohomolo", "testdata" );;
gap> CHMLINFO:=InfoLevel(InfoCohomolo);;
gap> SetInfoLevel(InfoCohomolo,1);

#
gap> Read( Filename ( testdatadir, "d8" ) );

#
gap> chr:= CHR(G,2,F,m2);;
gap> M:= SchurMultiplier( chr );
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
[ 2 ]
gap> D:= CoveringGroup( chr );;
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
gap> Size( D );
16
gap> FirstCohomologyDimension( chr );
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
2
gap> SecondCohomologyDimension( chr );
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
3
gap> Ex:= SplitExtensionCHR( chr );;
gap> Size(Ex);
32
gap> Ex:= NonsplitExtension( chr );;
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
gap> Size(Ex);
32
gap> Ex:= NonsplitExtension( chr, [0,1,1] );;
gap> Size(Ex);
32

#
gap> SetInfoLevel(InfoCohomolo,CHMLINFO);
gap> STOP_TEST( "d8.tst", 10000 );
