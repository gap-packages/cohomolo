##############################################################################
##
#A  coho4.g                     March 2000                      Derek Holt
##
##  Converted from the GAP 3.4.4 verion, Derek Holt 2000/03/10.
##  Changes are minimal - SplitExtension renamed SplitExtensionCHR
##
##  This file contains the interface to my cohomology `C' programs.
##
DeclareInfoClass("InfoCohomolo");

#############################################################################
##
#V  COHOMOLO . . . . . . . . . . . . . . . . . . . . . . . .  global variable
##
COHOMOLO := rec(
                 Multiplier := [],
                 CoDim1     := 0,
                 CoDim2     := 0,
                 RelVals    := [],
                 PermRels   := [],
                 RM_F       := "rm -f ",
                 CALL       := DirectoriesPackagePrograms("cohomolo")
                 #This is the directory of the external executables
                );

#############################################################################
##
#F  COHOMOLO.WritePermGroup( <G>, <deg>, <base>, <filename> )
##
##  <G> should be permutation group of degree <deg> and with base <base>.
##  The group is written to the file <filename> in the format of the
##  cohomology programs.
##
COHOMOLO.WritePermGroup := function ( G, deg, base, filename)

   local  gens, ng, nb, perm, i, j, stream;

   gens := GeneratorsOfGroup(G);
   ng := Length(gens);
   nb := Length(base);

   stream := OutputTextFile(filename, false);
   SetPrintFormattingStatus(stream, false);

   PrintTo(stream,deg," ",ng," ",nb," 0\n");
   for i in [1..nb] do
     PrintTo(stream," ",base[i]);
   od;
   PrintTo(stream,"\n");
   for i in [1..ng] do
     perm := gens[i];
     for j in [1..deg] do
         PrintTo(stream," ",j^perm);
     od;
     PrintTo(stream,"\n");
   od;
   CloseStream(stream);
end;

#############################################################################
##
#F  COHOMOLO.WriteFPGroup( <G>, <filename> )
##
##  <G> should be a finitely presented group.
##  The group is written to the file <filename> in the format of the
##  cohomology programs.
##
COHOMOLO.WriteFPGroup := function ( G, filename)
   local  gens, genstring, igens, allgens, ng, rels, nr, LL, LG, ILG, ALG,
          w, i, j, ct, gno, lgno, stream;

   gens := FreeGeneratorsOfFpGroup (G);
   ng := Length(gens);
   rels := RelatorsOfFpGroup(G);
   nr := Length(rels);
   igens := List(gens, x -> x^-1);
   allgens := Concatenation(gens,igens);
   LL := "abcdefghijklmnopqrstuvwxyz";
   LG := List( [ 1 .. ng ],  x -> [LL[x]] );
   ILG := List(LG, x -> Concatenation(x,"-"));
   ALG := Concatenation(LG,ILG);
   #lower case letters will be used as the generator names in the printout.

   stream := OutputTextFile(filename, false);
   SetPrintFormattingStatus(stream, false);

   PrintTo(stream,ng," 0 ",nr,"\n");
   genstring := List( [ 1 .. ng ],  x -> LL[x] );
   PrintTo(stream,genstring,"\n");
   for i in [1..nr] do
      w := rels[i];
      lgno := 0; ct := 0;
      for j in [1..Length(w)] do
        gno := Position(allgens,Subword(w,j,j));
        if gno = lgno then
           ct := ct+1; #current power of this generator
        else
           if ct>1 then PrintTo(stream,ct); fi;
           ct := 1;
           PrintTo(stream,ALG[gno]);
           lgno := gno;
        fi;
      od;
      if ct>1 then PrintTo(stream,ct); fi;
      PrintTo(stream,"\n");
   od;
   PrintTo(stream,"n\n");
   CloseStream(stream);
end;

#############################################################################
##
#F  COHOMOLO.WriteMatrices( <mats>, <filename> )
##
##  <mats> should be a list of square matrices, all of same dimension
##  over the field of order p, for some prime p.
##  They are written to the file <filename> in the format of the
##  cohomology programs.
##
COHOMOLO.WriteMatrices := function ( mats, filename)
   local  f, w, p, nm, dim, mat, row, entry, table, n, i, j, k, stream;

   f := Field(Flat(mats));
   p := Size( f );
   w := Z(p);

   #We first make ourselves a little table to convert field entries to integers.
   table := [];
   for i in [1..p-1] do
     table[LogFFE(i*w^0,w)+1] := i;
   od;
   nm := Length(mats);
   dim := Length(mats[1]);

   stream := OutputTextFile(filename, false);
   SetPrintFormattingStatus(stream, false);

   PrintTo(stream,p," ",dim," ",nm,"\n");
   for i in [1..nm] do
     PrintTo(stream,"\n");
     mat := mats[i];
     for j in [1..dim] do
        row := mat[j];
        for k in [1..dim] do
          entry := row[k];
          if entry = 0*w then n := 0; else n := table[LogFFE(entry,w)+1]; fi;
          PrintTo(stream," ",n);
        od;
        PrintTo(stream,"\n");
     od;
   od;
   CloseStream(stream);
end;


#############################################################################
##
#F  CHR( <G>, <p>, [<F>], [<mats>] ) . . . . . . make a cohomology record
##
##  This function makes a cohomology record. It should be called before
##  using any of the main cohomology functions.
##  <G> must be a permutation group, and <p> a prime (ususally dividing the
##  order of <G>).
##  <F> must either be 0 (if not required, but <mats> is), or a finitely
##  presented group with the same number of generators as <G> and mapping
##  epimorphically onto <G>. (This will be checked.)
##  (To obtain useful results, <F> would normally be isomorphic to <G>.)
##  <mats> (if present) should be a list of matrices, all of the same
##  dimension, over the field of order <p>.
##  The length of the list must be the same as the number of generators of
##  <G>, and the matrices must define a module for <G> over the field of
##  order <p>.
##  Apart from checking validity, a Sylow <p>-subgroup of <G> is calculated.
##  The record (with various components set) is returned.
##
CHR := function ( arg )
   local  chr, na, G, p, F, mats, Ggens, Fgens, ng, deg, dim, error, mat, rels,
          w, i, d;
   na := Number(arg);

   if na<2 or na>4 then
      Error("Number of arguments to CHR must be 2, 3 or 4");
   fi;
   G:=arg[1]; p:=arg[2]; F:=0; mats:=0;
   if na >= 3 then F:=arg[3]; fi;
   if na = 4 then mats:=arg[4]; fi;
   
   chr := rec();
   chr.isCHR := true;

   #First check arguments have correct types
   if not IsPermGroup(G) then
      Error("First argument of CHR must be a permutation group");
   fi;
   chr.permgp := G;
   Ggens := GeneratorsOfGroup(G);
   ng := Length(Ggens);
   if not IsInt(p) or not IsPrime(p) then
      Error("Second argument of CHR must be a prime number");
   fi;
   chr.prime := p;
   if not IsGroup(F) and F<>0 then
      Error("Third argument of CHR must be a finitely presented group or 0");
   fi;
   if IsGroup(F) then
      if not IsFpGroup(F) then
        Error("Third argument of CHR must be a finitely presented group");
      fi;
      Fgens := FreeGeneratorsOfFpGroup(F);
      if Length(Fgens)<>ng then 
        Error(
          "Arguments <G> and <F> of CHR must have same numbers of generators");
      fi;
   fi;
   chr.fpgp := F;
   error := false;
   if mats<>0 then
     if not IsList(mats) or Length(mats)<>ng then
       error := true;
     else
       for i in [1..ng] do
         mat := mats[i];
         if not IsMatrix(mat) or Field(Flat(mat))<>GF(p) then
            error := true;
         else
            d := DimensionsMat(mat);
            if i=1 then dim := d[1]; fi;
            if d[1]<>dim or d[2]<>dim or RankMat(mat)<>dim then
               error := true;
            fi;
         fi;
       od;
     fi;
     if error then
       Error(
     "Fourth argument of CHR must be a list of invertible matrices over GF(p)");
     fi;
   fi;
   chr.matrices := mats;

   #If F is present we now check that the permutations satisfy the relations,
   #and if mats is also present, we check that the matrices also satisfy them.
   if  IsGroup(F) then
      rels := RelatorsOfFpGroup(F);
      for w in rels do
         if MappedWord(w,Fgens,Ggens) <> () then
            Print("Relator: ",w,"\n");
            Error("This relator of F is not satisfied by the generators of G");
         fi;
         if mats<>0 and MappedWord(w,Fgens,mats) <> IdentityMat(dim,GF(p)) then
           Error("A relator of F is not satisfied by the matrices");
         fi;
      od;
   fi;

   #Now calculate some other fields - first degree of group.
   deg := 0;
   for i in Ggens do
      if i<>() and LargestMovedPointPerm(i)>deg then
         deg := LargestMovedPointPerm(i);
      fi;
   od;
   if deg > 4096 then
     Error(
  "Sorry - the cohomology programs are restricted to groups of degree <= 4096");
   fi;
   chr.degree := deg;
   #Next a base, and then a Sylow p-subgroup - and that's it for now!
   chr.base := BaseOfGroup(G);
   chr.sylow := SylowSubgroup(G,p);

   chr.verbose := false;
   
   return chr;
end;

#############################################################################
##
#F  IsCHR( x ) . . . . . . . . checks if the argument is a cohomology record.
##
IsCHR := function(x)
   return IsRecord(x) and IsBound(x.isCHR) and x.isCHR=true;
end;

#############################################################################
##
#F  FindSubSeq( <chr>, <n> ) . . . . . . find a subgroup chain
##
##  <chr> should be a cohomology record, and n a positive integer.
##  Let G=chr.permgp and P=chr.sylow.
##  FindSubSeq calculates a sequence [S1,...,St] of subgroups of G, where
##  Si > S(i+1) S1 = <G>,  S2 = <P> and the indices are a ssmall as possible.
##  It works, by trying to find a subgroup S2 of G with small index, and then
##  repeating with S2 in place of G.
##  The parameter <n> indicates how hard it should try.
##  It selects S2 as soon as it has found an S2 with index <= <n>.
##  The sequence is stored as chr.chain.
##
FindSubSeq := function ( chr, n )
   local  G, P, p, deg, seq, H, N, foundsub, orb, norbs, orbdone,
          sylind, bestind, ind, badind,
          orbsize, bestsub, orbno, orbstab, Z, OZ, OZG, centels, cent,
          el, elno, i, a, b, c;

   if not IsCHR(chr) then
     Error("First argument of FindSubSeq must be a cohomology record");
   fi;
   if IsBound(chr.chain) then
      return; #Must have been called already.
   fi;
   G := chr.permgp;  P := chr.sylow;  deg := chr.degree; p := chr.prime;
   seq := [G];
   sylind := Index(G,P);
   if sylind = 1 then
      chr.chain := seq;
      return;
   fi;

   #We will be trying to find subgroups of P as stabilizers of orbits of P,
   # so we first calculate the orbits.
   if sylind >= n then  
     orb := List(Orbits(P,[1..deg]),AsSet);
     norbs := Length(orb);
     orbdone := [];
     for i in [1..norbs] do
        orbdone[i] := false; # this means orbit number i is not yet processed.
     od;
   fi;

   #We also might try centralizers of central elements, so let's find
   #the centre of P, and then find a small random set of elements
   #of order p.
   Z := Centre(P);
   OZG := [];
   for el in GeneratorsOfGroup(Z) do
     Add(OZG,el^(Order(el)/p));
   od;
   OZ := Subgroup(Z,OZG);
   if IsCyclic(OZ) then
      centels := [OZG[1]];
   elif Length(OZG) = 2 then
     a := OZG[1]; b := OZG[2];
     if p=2 then centels:=[a,b,a*b]; elif p=3 then centels:=[a,b,a*b,a*b^2];
       elif p=5 then centels:=[a,b,a*b,a*b^2,a*b^3,a*b^4];
       else centels:=[a,b,a*b,a*b^2,a*b^3,a*b^4,a*b^5,a*b^6];
     fi;
   else
     a := OZG[1]; b := OZG[2]; c := OZG[3];
     centels := [a,b,c,a*b,a*c,b*c,a*b*c];
   fi;

   #Now we start looking for subgroups.
   H := G;
   foundsub := true;
   while sylind >= n and foundsub do
      foundsub := false;
      bestind := sylind; #The smallest index we have so far.

      #First try the normalizer of P
      N := Normalizer(H,P);
      if H <> N and P <> N then
        ind :=  Index(H,N);
        if ind < n then
           #Found a suitable subgroup!!
           H := N; Add(seq,H);
           foundsub := true; sylind := sylind/ind;
        else
           bestind := ind;  bestsub := N;
        fi;
      fi;

      #We next look for a suitable subgroup as a stabilizer of an orbit
      #of P - we consider the orbits in order of increasing size - we
      #go only up to size 16, since it takes too long to calculate set
      #stabilizers for large sets.
      orbsize := 1;
      while not foundsub and orbsize <= 16 do
         orbno := 1;
         while not foundsub and orbno <= norbs do
            if not orbdone[orbno] and Length(orb[orbno]) = orbsize then
               if orb[orbno] = Set(Orbit(H,orb[orbno][1])) then
                  orbstab := H;
               else
                  orbstab := Stabilizer(H,orb[orbno],OnSets);
               fi;
               if orbstab = H or orbstab = P then
                  orbdone[orbno] := true;
               else
                  ind := Index(H,orbstab);
                  if ind < n then
                     #Found a suitable subgroup!!
                     H := orbstab; Add(seq,H);
                     foundsub := true; sylind := sylind/ind;
                     orbdone[orbno] := true;
                  else   
                     #We could try the nomalizer!
                     N := Normalizer(H,orbstab);
                     if N<>H and N<>orbstab then
                        ind := Index(H,N);
                        orbstab := N;
                        if ind < n then
                           #Found a suitable subgroup!!
                           H := N; Add(seq,H);
                           foundsub := true; sylind := sylind/ind;
                        fi;
                     fi;
                     if ind < bestind then
                       #Not ideal, but the best we have so far, so remember it!
                       bestind := ind;  bestsub := orbstab;
                     fi;
                  fi;
               fi;
            fi;
            orbno := orbno + 1;
         od;
         orbsize := orbsize+1;
      od;

      if not foundsub then
         #Orbit stabilizers wasn't satisfactory, so we'll try
         #centralizers of central elements of P - 
         #We have already collected a few such elements of order p.
         elno := 1;
         while  not foundsub and elno <= Length(centels) do
            el := centels[elno];
            cent := Centralizer(H,el);
            if cent <> H and cent <> P then
              ind := Index(H,cent);
              if ind < n then
                 #Found a suitable subgroup!!
                 H := cent; Add(seq,H);
                 foundsub := true; sylind := sylind/ind;
              else   
                  #We could try the normalizer!
                  N := Normalizer(H,cent);
                  if N<>H and N<>cent then
                     ind := Index(H,N);
                     cent := N;
                     if ind < n then
                        #Found a suitable subgroup!!
                        H := N; Add(seq,H);
                        foundsub := true; sylind := sylind/ind;
                     fi;
                  fi;
                  if ind < bestind then
                 #Not ideal, but the best we have so far, so remember it!
                     bestind := ind;  bestsub := cent;
                  fi;
              fi;
            fi;
            elno := elno+1;
         od;
      fi;

      if not foundsub and bestind < sylind then
         #We make do with the best one we found!
         H := bestsub; Add(seq,H);
         foundsub := true; sylind := sylind/bestind;
      fi;
   od;

   N := Normalizer(H,P);
   if H <> N and P <> N then
     Add(seq,N);
   fi;

   Add(seq,P);

   if InfoLevel(InfoCohomolo)>=1 then
     Print("#I  Indices in the subgroup chain are:  ");
   fi;
   badind := false;
   for i in [1..Length(seq)-1] do
     ind := Index(seq[i],seq[i+1]);
     if ind>50000 then badind := true; fi;
     if InfoLevel(InfoCohomolo)>=1 then
       Print(Index(seq[i],seq[i+1])," ");
     fi;
   od;
   if InfoLevel(InfoCohomolo)>=1 then
     Print("\n");
   fi;

   if badind then
     Print(
  "#WARNING: An index in the subgroup chain found is larger than 50000.\n");
     Print(
            "#This calculation may fail. See manual for possible remedies.\n");
   fi;

   chr.chain := seq;
end;

#############################################################################
##
#F  Cohomolo( <chr>, <mult>, <pres>, <first>, <filename> )
##             . . . . perform cohomology calculation using external package.
##
##  This is the routine that calls the external routines to carry out the
##  cohomology calculations. It is not intended to be called directly, but
##  it is called by the access functions which follow.
##  <chr> must be a cohomology record, created with CHR.
##  <mult> is true for multiplier calculations, otherwise false
##         (and when false, chr.matrices must be defined).
##  <pres> is true if presentations of extensions are to be defined,
##         (and when true, chr.fpgp must be defined).
##  <first> is true for first cohomology group computation, in which case
##         mult and pres must be false.
##  <filename> is the name of the base of the files to be created for use of
##       the external package. The files have names of form <filename>.<suffix>.
##
Cohomolo := function( chr, mult, pres, first, filename )
   local deg, G, P, p, base, lc, nint, ch, F, mats, i, optstring, callstring,
         ct, ok;

   if not IsCHR(chr) then
     Error("First argument of Cohomolo must be a cohomology record");
   fi;

   F := chr.fpgp; mats := chr.matrices;
   if not mult and mats=0 then
     Error("Cohomolo: if mult is false, matrices field must be defined");
   fi;
   if first and mats=0 then
     Error("Cohomolo: if first is true, matrices field must be defined");
   fi;
   if first and (mult or pres) then
     Error("Cohomolo: if first is true, mult and pres must be false");
   fi;
   if pres and not IsFpGroup(F) then
     Error("Cohomolo: if pres is true, fpgp field must be defined");
   fi;

   G := chr.permgp; P := chr.sylow; deg := chr.degree; p := chr.prime;
   #Check that P really is a Sylow-subgroup for safety.
   if Set(FactorsInt(Size(P))) <> [p] or p in FactorsInt(Index(G,P)) then
      Error("Cohomolo: sylow field of <chr> is not a Sylow-subgroup!");
   fi;
   #First check for trivial cases.
   if IsTrivial(P) then
     Print(
"The Sylow p-subgroup of the group is trivial - so all cohomology is trivial.\n"
     );
     if mats<>0 then chr.codim1 := 0; chr.codim2 := 0; fi;
     chr.multiplier := [];
     if IsFpGroup(F) then chr.multrelvals := []; fi;
     if IsFpGroup(F) and mats<>0 then chr.modrelvals:=[]; fi;
     return;
   fi;
   if mult and IsCyclic(P) then
     Print(
"The Sylow p-subgroup of the group is cyclic - so the multiplier is trivial.\n"
     );
     chr.multiplier := []; chr.multrelvals := [];
     return;
   fi;
   
   FindSubSeq(chr,20);

   #We will check the validity of the chain for safety - 
   #just in case the user may have supplied it!

   ct := 1; ok := true;
   ch := chr.chain; lc := Length(ch);
   while ok and ct <= lc do
     if ct=1 and ch[ct]<>G  then ok := false;
     elif ct=lc and ch[ct]<>P then ok := false;
     elif ct<lc and (not IsSubgroup(ch[ct],ch[ct+1]) or ch[ct]=ch[ct+1]) then
        ok := false;
     fi;
     ct := ct+1;
   od;
   if not ok then
      Error(
    "Cohomolo: <chr>.chain must be a strictly decreasing chain from G to P.");
   fi;

   ## chr.norm is true iff the penultimate member of the chain normalises P.
   ## chr.neqg is true iff the chain has length<=2 and P is normal in G.
   chr.norm := lc > 1 and IsNormal(ch[lc-1],ch[lc]);
   chr.neqg := lc=1 or (lc=2 and chr.norm);
   #Now write things to files -
   #but first make sure there is no rubbish left from a previous run.
   Exec( Concatenation( COHOMOLO.RM_F, filename, ".*" ) );
   base := chr.base;
   COHOMOLO.WritePermGroup(G,deg,base,Concatenation(filename,".inperm"));
   if lc > 1 then
     COHOMOLO.WritePermGroup(P,deg,base,Concatenation(filename,".gp"));
     if chr.norm then
       COHOMOLO.WritePermGroup(ch[lc-1],deg,base,Concatenation(filename,".gn"));
       nint := lc-3;
     else
       nint := lc-2;
     fi;
     if nint > 0 then
        for i in [1..nint] do
           COHOMOLO.WritePermGroup(ch[nint+2-i],deg,base,
                                      Concatenation(filename,".g",String(i)));
        od;
     fi;
   fi;

   if pres then
      COHOMOLO.WriteFPGroup(F,Concatenation(filename,".tc"));
      PrintTo(Concatenation(filename,".sg.rel"),Length(base),"\n");
      #The last is a technicality to avoid a silly problem.
   fi;
   if pres and not mult then
      #Write input for the external program nqip. This is the number of
      #the basis of H^2(G,M) that is to be used for corestriction and
      #calculation of an extension. This number always starts at 1.
      #If H^2(G,M) has dimension>1,then we will repeat later with higher values.
      PrintTo(Concatenation(filename,".nqip"),1,"\n");
   fi;

   if not mult then
      COHOMOLO.WriteMatrices(mats,Concatenation(filename,".inmat"));
   fi;
   
   #Now work out the command for the external program.
   callstring := Filename(COHOMOLO.CALL, "cohomology.gap");
   optstring := "-";
   if mult        then Add( optstring, 'm'); fi;
   if pres        then Add( optstring, 'c'); fi;
   if first       then Add( optstring, '1'); fi;
   if chr.norm    then Add( optstring, 'n'); fi;
   if chr.neqg    then Add( optstring, 'e'); fi;
   if chr.verbose then Add( optstring, 'v'); fi;
   if 1 < Length( optstring ) then
      callstring := Concatenation(callstring," ", optstring );
   fi;
   callstring:= Concatenation( callstring,
                               " ", String( chr.prime ), " ", filename );
   if chr.verbose then Print(callstring,"\n"); fi;

   Info(InfoCohomolo,1," Cohomolo package: Calling external program." );
   if InfoLevel(InfoCohomolo)>=2 then
     Print("#I", callstring, "\n");
   fi;
   Exec( callstring );
   Info(InfoCohomolo,1," External program complete." );

   if mult then

      if not READ( Concatenation( filename, ".mult" ) ) then
        Error( "'Cohomolo' failed for some reason.\n" );
      fi;
      chr.multiplier := COHOMOLO.Multiplier;

   elif READ(Concatenation(filename,".cdim")) then
      if first then chr.codim1 := COHOMOLO.CoDim1;
      else  chr.codim2 := COHOMOLO.CoDim2;
      fi;
   else
     Error( "'Cohomolo' failed for some reason.\n" );
   fi;

   if pres then
      if mult then
        if chr.multiplier=[] then chr.multrelvals := [];
        elif READ(Concatenation(filename,".rvals")) then
           chr.multrelvals := COHOMOLO.RelVals;
        else
           Error( "'Cohomolo' failed for some reason.\n" );
        fi;
      else
         if chr.codim2=0 then chr.modrelvals := []; 
         elif READ(Concatenation(filename,".rvals")) then
           chr.modrelvals := [COHOMOLO.RelVals];
         else
           Error( "'Cohomolo' failed for some reason.\n" );
         fi;
      fi;
   fi;

   if pres and not mult and chr.codim2>1 then

     #We have to calculate the other extensions for the basis of H^2(G,M).
      Add( optstring, 'r' );
      callstring := Filename(COHOMOLO.CALL, "cohomology.gap");
      callstring := Concatenation(callstring," ",
                      optstring, " ", String(chr.prime), " ", filename );
            
      for i in [2..chr.codim2] do
         PrintTo(Concatenation(filename,".nqip"),i,"\n");
         if chr.verbose then Print(callstring,"\n"); fi;

         Info(InfoCohomolo,1,
                    " Cohomolo package: Calling external program." );
         Exec( callstring );
         Info(InfoCohomolo,1," External program complete." );

         if not READ(Concatenation(filename,".rvals")) then
           Error( "'Cohomolo' failed for some reason.\n" );
         fi;
         chr.modrelvals[i] := COHOMOLO.RelVals;
      od;
   fi;

   Info(InfoCohomolo,1," Removing temporary files." );
   Exec( Concatenation( COHOMOLO.RM_F, filename, ".*" ) );
end;

#############################################################################
##
#F  CalcExtPres( <chr>, <mult>, [<vec>] ) .  calculate presentation of extension
##
##  <chr> must be a cohomology record.
##  Let <chr>.permgp=G, <chr>.prime=p  and  <chr>.fpgp =F. 
##  CalcExtPres calculates a presentation of an extension of M by F,
##  where M is either the p-part of the Schur multiplier of G, when <mult> is
##  true, or the F-module over GF(p) defined by <chr>.matrices, when <mult>
##  is false.
##  In the first case, Cohomolo must already have been called on <chr> with
##  <mult> and <pres> true, and the presentation is of a p-cover of G
##  (assuming F and G are isomorphic).
##  In the second case, <vec> must be a list of integers, and either <vec>=[],
##  in which case the split extension is returned, or Cohomolo must already
##  been called already with <mult> false and <pres> true, and the length of
##  <vec> must be equal to <chr>.codim2. The extension returned is the one
##  corresponding to the element <vec> of H^2(G,M).
##  This function is not intended to be called directly - it is accessed by
##  the individual functions, and by ExtPresentation.
##
CalcExtPres := function( arg )
    local chr, mult, vec, E, EG, p, F, ng, codim, mp, dim, rels, rel, relval, w,
          mats, row, modrv, multrv, f, table, rt, i, j, k, Erels;

   chr := arg[1]; mult := arg[2];
   if (mult and Number(arg)<>2) or (not mult and Number(arg)<>3) then
      Error("Wrong number of arguments to CalcExtPres");
   fi; 
   if not mult then
      vec := arg[3];
      if not IsList(vec) then
         Error("Third argument of CalcExtPres must be a list");
      fi;
   fi;
   if not IsCHR(chr) then
     Error("First argument of CalcExtPres must be a cohomology record");
   fi;

   F := chr.fpgp; mats := chr.matrices; p := chr.prime;
   if not mult and mats=0 then
     Error("CalcExtPres: if mult is false, matrices field must be defined");
   fi;
   if not IsFpGroup(F) then
     Error("CalcExtPres: if pres is true, fpgp field must be defined");
   fi;
   if mult then
      if not IsBound(chr.multiplier) or not IsBound(chr.multrelvals) then
         Error("CalcExtPres: Cohomolo must be run before CalcExtPres");
      fi;
      multrv := chr.multrelvals;
      mp := chr.multiplier;
      dim := Length(mp);
   else
      if Length(vec)<>0 then
         if not IsBound(chr.codim2) or not IsBound(chr.modrelvals) then
            Error("CalcExtPres: Cohomolo must be run before CalcExtPres");
         fi;
         codim := chr.codim2;
         if  Length(vec)<>codim then
            Error("CalcExtPres: Length of third argument wrong");
         fi;
         modrv := chr.modrelvals;
      fi;
      #We make ourselves a little table to convert field entries to integers.
      f := Field(Flat(mats));
      rt := Z(Size(f));
      table := [];
      for i in [1..p-1] do
        table[LogFFE(i*rt^0,rt)+1] := i;
      od;
      dim := Length(mats[1]);
   fi;

   #Now we can start to build up the extension E
   ng := Length(FreeGeneratorsOfFpGroup(F));
   E := FreeGroup(ng+dim);
   EG := GeneratorsOfGroup(E);
   Erels := [];
   #First the relators that say the normal subgroup is finite abelian.
   for i in [1..dim-1] do for j in [i+1..dim] do
      Add(Erels,Comm(EG[j+ng],EG[i+ng]));
   od; od;
   for i in [1..dim] do
      if mult then j:=mp[i]; else j:=p; fi;
      Add(Erels,EG[i+ng]^j);
   od;
   #Now the module action relators
   for i in [1..ng] do for j in [1..dim] do
      w := EG[i]^-1*EG[ng+j]*EG[i];
      if mult then w := w/EG[ng+j];
      else
         row := mats[i][j];
         for k in [1..dim] do
           if row[k]<>0*rt then
             w := w/EG[ng+k]^table[LogFFE(row[k],rt)+1];
           fi;
         od;
      fi;
      Add(Erels,w);
   od; od;
   #And finally the relators, that tell us the values of the relators of
   #F in the normal subgroup.
   rels := RelatorsOfFpGroup(F);
   for i in [1..Length(rels)] do
      rel := rels[i];
      w := MappedWord(rel,FreeGeneratorsOfFpGroup(F),EG{[1..ng]});
      if mult then
         if mp=[] then relval := [];
         else relval :=  multrv[i];
         fi;
      else 
         relval := [];
         for j in [1..dim] do relval[j]:=0; od;
         if vec<>[] then
            for j in [1..dim] do
               for k in [1..codim] do
                  relval[j] := relval[j] + vec[k]*modrv[k][i][j];
               od;
               relval[j] := relval[j] mod p;
            od;
         fi;
      fi;
      for j in [1..dim] do
        if relval[j]<>0 then
          w := w/EG[ng+j]^relval[j];
        fi;
      od;
      Add(Erels,w);
   od; 
   return E/Erels;
end;

#############################################################################
##
#F  SchurMultiplier( <chr> )
##      . . . . calculate the p-part of the Schur multiplier of <chr>.permgp
##
##  <chr> is a cohomology record.
##  The result is returned as a list of the abelian invariants.
##
SchurMultiplier := function( chr )
   if not IsCHR(chr) then
     Error("First argument of Cohomolo must be a cohomology record");
   fi;
   if not IsBound(chr.multiplier) then
      Cohomolo(chr,true,false,false, TmpName() );
   fi;
   return chr.multiplier;
end;

#############################################################################
##
#F  CoveringGroup( <chr> )  . . . calculate p-covering group of <chr>.permgp 
##
##  To use this <chr>.fpgp must be defined, and the answer will be accurate only
##  if <chr>.fpgp and <chr>.permgp are isomorphic.
##  The p-cover is returned as an fp-group.
##
CoveringGroup := function( chr )
   if not IsCHR(chr) then
     Error("First argument of Cohomolo must be a cohomology record");
   fi;
   if not IsBound(chr.multrelvals) then
      Cohomolo(chr,true,true,false, TmpName() );
   fi;
   return CalcExtPres(chr,true);
end;

#############################################################################
##
#F  FirstCohomologyDimension( <chr> ) 
##  . . dimension of H^1(G,M), with G=<chr>.permgp, M=module of <chr>.matrices 
##
##  <chr> is a cohomology record with <chr>.matrices defined.
##
FirstCohomologyDimension := function( chr )
   if not IsCHR(chr) then
     Error("First argument of Cohomolo must be a cohomology record");
   fi;
   if not IsBound(chr.codim1) then
      Cohomolo(chr,false,false,true, TmpName() );
   fi;
   return chr.codim1;
end;

#############################################################################
##
#F  SecondCohomologyDimension( <chr> ) 
##  . . dimension of H^2(G,M), with G=<chr>.permgp, M=module of <chr>.matrices 
##
##  <chr> is a cohomology record with <chr>.matrices defined.
##
SecondCohomologyDimension := function( chr )
   if not IsCHR(chr) then
     Error("First argument of Cohomolo must be a cohomology record");
   fi;
   if not IsBound(chr.codim2) then
     Cohomolo(chr,false,false,false, TmpName() );
   fi;
   return chr.codim2;
end;

#############################################################################
##
#F  SplitExtensionCHR( <chr> )
## . . presentation of split extension of <chr>.matrices by <chr>.fpgp
##
##  To use this <chr>.fpgp and <chr>.matrices must both be define
##  This is a routine calculation and does not call Cohomolo.
##
SplitExtensionCHR := function( chr )
   if not IsCHR(chr) then
     Error("First argument of Cohomolo must be a cohomology record");
   fi;
   return CalcExtPres(chr,false,[]);
end;

#############################################################################
##
#F  NonsplitExtension( <chr> , [<vec>] )
##  . . presentation of a nonsplit extension of <chr>.matrices by <chr>.fpgp
##
##  To use this <chr>.fpgp and <chr>.matrices must both be define
##  The equivalence classes of nonsplit extensions of M by G are in one-one
##  correspondence with the elements of H^2(G,M). Therefore an element of
##  H^2(G,M) needs to be specified. If the optional argument <vec> is
##  present, it must be a list of integers of length Dim(H^2(G,M)), which
##  specifies the element of H^2(G,M) as a vector (the integers may as well
##  be in range [0...p-1] but that is no compulsory).
##  By default, <vec> is taken to be [1,0,...,0], where the length is equal
##  to Dim(H^2(G,M)).
##
NonsplitExtension := function( arg )
   local chr, vec, i;
   if Number(arg)<>1 and Number(arg)<>2 then
     Error("Number of arguments of NonSplitExtension wrong");
   fi;
   chr := arg[1];
   vec := [];
   if Number(arg)=2 then vec := arg[2]; fi;
   if not IsCHR(chr) then
     Error("First argument of Cohomolo must be a cohomology record");
   fi;
   if not IsBound(chr.modrelvals) then
      Cohomolo(chr,false,true,false, TmpName() );
   fi;
   if chr.codim2 = 0 then
      Error("There is no nonsplit extension.");
   fi;
   if vec=[] then
      vec[1]:=1;
      for i in [2..chr.codim2] do vec[i]:=0; od;
   fi;
   return CalcExtPres(chr,false,vec);
end;

#############################################################################
##
#F  PermRep( <F>, <K> ). . . calculate permutation represenation of fp-group
##
##  <F> should be a finitely presented group and <K> a subgroup of finite
##  index.
##  This function return the permutation group giving the action of <F> on
##  right cosets of <K>.
##  Of course, there is no guarantee in general that the perm rep of <F> will
##  be faithful.
##  It is a straightforward application of standard GAP functions.
##
PermRep := function(F,K)
   local  ng, Ggens, GGi, CT, i;
   if not IsFpGroup(F) or not IsSubgroup(F,K) then
     Error("PermRep(F,K): K must be a subgroup of the fp-group F");
   fi;
   ng := Length(FreeGeneratorsOfFpGroup(F));
   CT := CosetTable(F,K);
   GGi := List(CT,PermList); #list of generators and inverses.
   Ggens := [];
   for i in [1..ng] do Ggens[i] := GGi[2*i-1]; od;
   return Group(Ggens,());
end;

#############################################################################
##
#F  CalcPres( <chr> ). . . calculate a presentation of <chr>.permgp
##
##  <chr> must be a cohomology record, created with CHR.
##  This routine can be used to calculate a presentation of <chr>.permgp,
##  when none is known to begin with (i.e. when <chr>.fpgp is zero).
##  The computed presentation is stored as <chr>.fpgp.
##  It calls an external program to do this.
##  It is only really useful for moderately small groups - currently, its
##  use is restricted to groups of order at most 32767.
##
CalcPres := function( chr )
   local deg, G, base, F, Fg, Fr, optstring, callstring, ng, ordgen, o, rel,
   i, l, w, g, h, ct, neg, filename;

   if not IsCHR(chr) then
     Error( "<chr> must be a cohomology record");
   fi;

   F := chr.fpgp;
   if IsFpGroup(F) then
      Print("CalcPres: presentation is already known.\n");
      return;
   fi;

   G := chr.permgp;  deg := chr.degree;

   #Now write group to file -
   base := chr.base;
   filename:= TmpName();
   COHOMOLO.WritePermGroup(G,deg,base,Concatenation( filename, ".inperm" ) );
   
   #Now work out the command for the external program.
   callstring := Filename( COHOMOLO.CALL, "calcpres.gap" );
   optstring := "-";
   if chr.verbose then Add( optstring, 'v' ); fi;
   if 1 < Length( optstring ) then
      callstring := Concatenation(callstring," ", optstring );
   fi;
   callstring := Concatenation(callstring," ", filename );
   if chr.verbose then Print(callstring,"\n"); fi;

   Info(InfoCohomolo,1," Calling external program." );
   Exec(callstring);
   Info(InfoCohomolo,1," External program complete." );

   if not READ(Concatenation( filename, ".reg.relg" ) ) then
      Error( "'CalcPres' failed for some reason.\n" );
   fi;
   chr.permrels := COHOMOLO.PermRels;

   Info(InfoCohomolo,1," Removing temporary files." );
   Exec( Concatenation( COHOMOLO.RM_F, filename, ".*" ) );

   # Now we build the finitely presented group.
   ng := Length(GeneratorsOfGroup(G));
   F := FreeGroup(ng);
   Fg := GeneratorsOfGroup(F);
   Fr := [];
   #first look for those relators that give the orders of generators.
   ordgen := [];
   for i in [1..ng] do ordgen[i] := 0; od;
   for rel in chr.permrels do
      if Size(Set(rel))=1 then
        g := rel[1];
        if g<0 then g:= -g; fi;
        o := Length(rel);
        ordgen[g] := o;
        Add(Fr,Fg[g]^o);
      fi;
   od;
   #Now process the other relators.
   for rel in chr.permrels do
      if Size(Set(rel))>1 then
        w := One(F);
        ct := 0; l := Length(rel); g := 0;
        for i in [1..l+1] do
           if i<=l then h := rel[i]; else h:=0; fi;
           if h<>g and ct>0 then
              neg := g<0;
              if neg then g := -g; fi;
              o := ordgen[g];
              if o=0 then
                 if neg then w:=w*Fg[g]^-ct; else w:=w*Fg[g]^ct; fi;
              elif neg and 2*ct>=o then
                 w:=w*Fg[g]^(o-ct);
              elif not neg and 2*ct>=o+1 then
                 w:=w*Fg[g]^(ct-o);
              else
                 if neg then w:=w*Fg[g]^-ct; else w:=w*Fg[g]^ct; fi;
              fi;
              ct := 0;
           fi;
           g := h; ct := ct+1;
        od;
        if w<> One(F) then
          Add(Fr,w);
        fi;
      fi;
   od;

   Unbind(chr.permrels);
   chr.fpgp := F/Fr;
end;

