gap> START_TEST("testall.tst");
gap> testdatadir := DirectoriesPackageLibrary( "cohomolo", "testdata" );;
gap> CHMLINFO:=InfoLevel(InfoCohomolo);;
gap> SetInfoLevel(InfoCohomolo,1);
gap> Read( Filename ( testdatadir, "d8" ) );
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
gap> F:= FirstCohomologyDimension( chr );
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
2
gap> F:= SecondCohomologyDimension( chr );
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
gap> Read( Filename ( testdatadir, "a6" ) );
gap> chr:= CHR( G, 3, F, m3 );;
gap> M:= SchurMultiplier( chr );
#I  Indices in the subgroup chain are:  10 4 
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
[ 3 ]
gap> D:= CoveringGroup( chr );;
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
gap> Size(D);
1080
gap> F:= SecondCohomologyDimension( chr );
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
2
gap> Ex:= SplitExtensionCHR( chr );;
gap> SE:= Subgroup( Ex, [Ex.4,Ex.5,Ex.6,Ex.7,Ex.8] );;
gap> Index( Ex, SE );
1080
gap> Ex:= NonsplitExtension( chr, [1,2] );;
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
gap> SE:= Subgroup( Ex, [Ex.4,Ex.5,Ex.6,Ex.7,Ex.8] );;
gap> Index( Ex, SE );
1080
gap> Read( Filename( TestDir, "sz8" ) );
gap> chr:= CHR( G, 2, 0, m2 );;
gap> orig_size:= SizeScreen();;
gap> SizeScreen( [ 72 ] );;
gap> SecondCohomologyDimension( chr );
#I  Indices in the subgroup chain are:  65 7 
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
0
gap> SizeScreen( orig_size );;
gap> SetInfoLevel(InfoCohomolo,CHMLINFO);
gap> STOP_TEST( "testall.tst", 10000 );
