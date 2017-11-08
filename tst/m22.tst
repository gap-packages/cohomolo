gap> START_TEST("m22.tst");
gap> testdatadir := DirectoriesPackageLibrary( "cohomolo", "testdata" );;
gap> CHMLINFO:=InfoLevel(InfoCohomolo);;
gap> SetInfoLevel(InfoCohomolo,1);

#
gap> Read( Filename ( testdatadir, "m22" ) );

#
gap> chr:= CHR( G, 2, F );;
gap> M:= SchurMultiplier( chr );
#I  Indices in the subgroup chain are:  77 15 3 
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
[ 4 ]
gap> FirstCohomologyDimension( chr );
Error, Cohomolo: if mult is false, matrices field must be defined
gap> SecondCohomologyDimension( chr );
Error, Cohomolo: if mult is false, matrices field must be defined

#
gap> chr:= CHR( G, 3, F );;
gap> M:= SchurMultiplier( chr );
#I  Indices in the subgroup chain are:  22 280 8 
#I   Cohomolo package: Calling external program.
#I   External program complete.
#I   Removing temporary files.
[ 3 ]
gap> FirstCohomologyDimension( chr );
Error, Cohomolo: if mult is false, matrices field must be defined
gap> SecondCohomologyDimension( chr );
Error, Cohomolo: if mult is false, matrices field must be defined

#
gap> SetInfoLevel(InfoCohomolo,CHMLINFO);
gap> STOP_TEST( "m22.tst", 10000 );
