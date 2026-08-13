gap> START_TEST("exec.tst");

#
gap> CHMLINFO:=InfoLevel(InfoCohomolo);;
gap> SetInfoLevel(InfoCohomolo,0);

# a program that fails must not be mistaken for one that has run (issue #37);
# it prints its usage message to stdout, which is not part of the test output
gap> COHOMOLO.Run("cohomology.gap", ["-q"]);
Error, the cohomolo program cohomology.gap failed with exit code 1

#
gap> name := TmpName();;
gap> PrintTo(Concatenation(name, ".one"), "");
gap> PrintTo(Concatenation(name, ".two"), "");
gap> COHOMOLO.RemoveTempFiles(name);
gap> List([".one", ".two"], s -> IsExistingFile(Concatenation(name, s)));
[ false, false ]

#
gap> SetInfoLevel(InfoCohomolo,CHMLINFO);
gap> STOP_TEST( "exec.tst", 10000 );
