#############################################################################
##
#W    init.g               share package 'cohomolo'            Derek Holt
##
##    @(#)$Id: init.g,v 1.0 2000/04/19 16:50:50 andrews Exp $
##

# announce the package version and test for the existence of the binary
DeclarePackage("cohomolo","1.0",
  function()
  local path,file;
    # test for existence of the compiled binary
    path:=DirectoriesPackagePrograms("cohomolo");
    file:=Filename(path,"extprun");
    if file=fail then
      Info(InfoWarning,1,
 "Package ``cohomolo'': The program `extprun' (for example) is not compiled");
    fi;
    return file<>fail;
  end);

# install the documentation
DeclarePackageAutoDocumentation( "cohomolo", "doc" );

