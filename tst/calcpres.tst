gap> START_TEST("calcpres.tst");
gap> CHMLINFO:=InfoLevel(InfoCohomolo);;
gap> SetInfoLevel(InfoCohomolo,0);

# CalcPres presents the group on its own generators
gap> G := Group((1,2,3,4,5),(1,2,3));;
gap> chr := CHR( G, 2 );;
gap> CalcPres( chr );
gap> F := chr.fpgp;;
gap> Length( GeneratorsOfGroup( F ) ) = Length( GeneratorsOfGroup( G ) );
true
gap> Size( F ) = Size( G );
true
gap> ForAll( RelatorsOfFpGroup( F ), r -> MappedWord( r,
>        FreeGeneratorsOfFpGroup( F ), GeneratorsOfGroup( G ) ) = () );
true
gap> CalcPres( chr );
CalcPres: presentation is already known.
gap> CalcPres( 3 );
Error, <chr> must be a cohomology record

# PermRep acts on the cosets of a subgroup
gap> K := Subgroup( F, [ F.1 ] );;
gap> P := PermRep( F, K );;
gap> NrMovedPoints( P ) = Index( F, K );
true
gap> Size( P ) = Size( G );
true
gap> PermRep( G, K );
Error, PermRep(F,K): K must be a subgroup of the fp-group F

#
gap> SetInfoLevel(InfoCohomolo,CHMLINFO);
gap> STOP_TEST( "calcpres.tst", 10000 );
