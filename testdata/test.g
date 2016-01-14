TestDir:= "../testdata";

SetInfoLevel(InfoCohomolo,1);
Read( Concatenation( TestDir, "d8" ) );
chr:= CHR(G,2,F,m2);;
M:= SchurMultiplier( chr );
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
#>[ 2 ]

D:= CoveringGroup( chr );;
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
Size( D );
#>16

F:= FirstCohomologyDimension( chr );
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
#>2

F:= SecondCohomologyDimension( chr );
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
#>3

Ex:= SplitExtensionCHR( chr );;
Size(Ex);
#>32

Ex:= NonsplitExtension( chr );;
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
Size(Ex);
#>32

Ex:= NonsplitExtension( chr, [0,1,1] );;
Size(Ex);
#>32

Read( Concatenation( TestDir, "a6" ) );
chr:= CHR( G, 3, F, m3 );;
M:= SchurMultiplier( chr );
#>#Indices in the subgroup chain are:  10 4 
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
#>[ 3 ]

D:= CoveringGroup( chr );;
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
Size(D);
#>1080

F:= SecondCohomologyDimension( chr );
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
#>2

Ex:= SplitExtensionCHR( chr );;
SE:= Subgroup( Ex, [Ex.4,Ex.5,Ex.6,Ex.7,Ex.8] );;
Index( Ex, SE );
#>1080

Ex:= NonsplitExtension( chr, [1,2] );;
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Cohomology package: Calling external program.
#>#I  External program complete.
#>#I  Removing temporary files.
SE:= Subgroup( Ex, [Ex.4,Ex.5,Ex.6,Ex.7,Ex.8] );;
Index( Ex, SE );
#>1080

