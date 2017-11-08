gap> START_TEST("sz8.tst");
gap> testdatadir := DirectoriesPackageLibrary( "cohomolo", "testdata" );;
gap> CHMLINFO:=InfoLevel(InfoCohomolo);;
gap> SetInfoLevel(InfoCohomolo,1);

#
gap> Read( Filename( testdatadir, "sz8" ) );
gap> chr:= CHR( G, 2, 0, m2 );;

# change screen size to verify line wrapping bug is gone
gap> orig_size:= SizeScreen();;
gap> SizeScreen( [ 72 ] );;

#
gap> FirstCohomologyDimension( chr );
#I  Indices in the subgroup chain are:  65 7 
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
3
gap> SecondCohomologyDimension( chr );
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
0

#
gap> SizeScreen( orig_size );;

#
gap> SetInfoLevel(InfoCohomolo,CHMLINFO);
gap> STOP_TEST( "sz8.tst", 10000 );
