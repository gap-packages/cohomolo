#! /bin/sh
#Gap version of shell-script to find a presnetation of a permutation group
#on the given generators.
#Tests and checks that are done within GAP have been removed.
#First set the path directory
DIR=`echo $0 | sed -e 's/[^\/]*$//'`
RM=/bin/rm
verbose=false
gpname=
for i
do case $i in
  -*)
   flags="`echo ' '$i|awk '{for (j=2;j<=length($1);j++) print substr($1,j,1)}'`"
    for j in $flags
    do case $j in
         v) verbose=true;;
         *) echo Usage:  calcpres [-v] gpname; exit 1;;
       esac
    done;;
  *) gpname=$i
  esac
done

if test $verbose = true
then cmddest=
else cmddest=" > /dev/null"
fi

export gpname verbose cmdsource cmddest cmd
    
cmdsource=
cmd="${DIR}gprun -b $gpname"; ${DIR}execcmd.gap || exit 1
cmd="${DIR}conrun -h $gpname outperm triv inperm";
     ${DIR}execcmd.gap || exit 1
mv ${gpname}.inperm.nr ${gpname}.reg.inperm
cmd="${DIR}gprun -b ${gpname}.reg"; ${DIR}execcmd.gap || exit 1

echo "1.5 50 30 1" > ${gpname}.grip
cmdsource="< ${gpname}.grip"
cmd="${DIR}grrun -g ${gpname}.reg outperm"; ${DIR}execcmd.gap || exit 1
