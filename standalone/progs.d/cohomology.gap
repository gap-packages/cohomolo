#! /bin/sh
#Gap version of cohomology shell-script
#Tests and checks that are done within GAP have been removed.
#First set the path directory
DIR=`echo $0 | sed -e 's/[^\/]*$//'`
RM=/bin/rm
newjob=true co=false mult=false check=true first=false neqg=false norm=false
verbose=false
gpname=
for i
do case $i in
  -*) 
   flags="`echo ' '$i|awk '{for (j=2;j<=length($1);j++) print substr($1,j,1)}'`"
    for j in $flags
    do case $j in
         r) newjob=false co=true;;
         c) co=true;;
         n) norm=true;;
         e) neqg=true;;
         m) mult=true;;
         x) check=false;;
         1) first=true;;
         s) step=true;;
         v) verbose=true;;
         *) echo Usage:  cohomology [-rcnemx1sv] prime gpname; exit 1;;
       esac
    done;;
  [1-9]*) prime=$i;;
  *) gpname=$i
  esac
done

if test $verbose = true
then cmddest=
else cmddest=" > /dev/null"
fi

export gpname cmd cmddest cmdsource step verbose

if test -r ${gpname}.tc
then tc=true
else tc=false
fi

case $newjob in
true)
    if test $mult = true
    then flag="-b"
    else flag="-b -w"
    fi

    echo $prime > ${gpname}.sylip
    
    cmdsource=
    cmd="${DIR}gprun -n $flag $gpname"; ${DIR}execcmd.gap || exit 1
    
    if test $co = true -a $tc = true
    then flag=-f
    else flag=
    fi
    
    cmd="${DIR}egrun $flag $gpname"; ${DIR}execcmd.gap || exit 1

    if test $neqg = true -a $norm = false 
    then cp ${gpname}.sg ${gpname}.psg
    else cmd="${DIR}gprun -b $gpname gp ogp"; ${DIR}execcmd.gap || exit 1
         cmd="${DIR}egrun $gpname ogp psg"; ${DIR}execcmd.gap || exit 1
         ${RM} -f ${gpname}.*gp
    fi
         
    cmdsource="< ${gpname}.sylip"
    cmd="${DIR}sylrun $gpname psg"; ${DIR}execcmd.gap || exit 1
    cmdsource=

    if test $norm = true
    then  if test $neqg = true
          then cmd="${DIR}normrun -n $gpname sg psg"; ${DIR}execcmd.gap || exit 1
               cp ${gpname}.sg ${gpname}.nsg
          else cmd="${DIR}gprun -b $gpname gn ogn"; ${DIR}execcmd.gap || exit 1
               cmd="${DIR}egrun $gpname ogn nsg"; ${DIR}execcmd.gap || exit 1
               ${RM} -f ${gpname}.*gn
               cmd="${DIR}normrun -n $gpname nsg psg"; ${DIR}execcmd.gap || exit 1
          fi
    fi
    
    case $mult in
      true) flag1=-m;;
      false) flag1=;;
    esac
    case $co in
      true) flag2=-c;;
      false) flag2=;;
    esac
    
    cmd="${DIR}pcrun $flag1 $flag2 $gpname"; ${DIR}execcmd.gap || exit 1
    cmd="${DIR}selgen -w $gpname"; ${DIR}execcmd.gap || exit 1
    
    if test $co = true -o '(' $check = true -a $mult = false ')'
    then  echo "1.5 50 30 1" > ${gpname}.grip
          cmdsource="< ${gpname}.grip"
    fi
    
    case $co in
      true) case $mult in
            true) cflag1=-g; repfile=;;
            false) cflag1="-g -c"; repfile=cr0;;
            esac
            cmd="${DIR}grrun $gpname psg"; ${DIR}execcmd.gap || exit 1;;
      false) cflag1=; repfile=
        if test $neqg = true -a $mult = false -a $norm = false -a $check = true
        then
        cmd="${DIR}grrun $gpname psg"; ${DIR}execcmd.gap || exit 1
        fi;;
    esac
    cmdsource=
    
    case $mult in
    true) scflag1=-m;;
    false) scflag1=;;
    esac
    
    list=
    case $mult in
    false) case $neqg in
             true) mc3=pg;;
             false) mc3="pg dcr";
           esac
           case $co in
           true)  case $neqg in
                  true) mc1=psg; mc2=;;
                  false) mc1="psg sg"; mc2=-cr;;
                  esac;;
           false) mc1=;  mc2=;;
           esac;;
    esac
    
    case $norm in
      true)
        case $co in
        true)
          cmd="${DIR}conrun $cflag1 $gpname nsg psg $repfile"
          ${DIR}execcmd.gap || exit 1
          cmdsource="< ${gpname}.grip"
          cmd="${DIR}grrun $gpname nsg"; ${DIR}execcmd.gap || exit 1;;
        false)
          if test $mult = false -a $neqg = true -a $check = true
          then
            cmdsource="< ${gpname}.grip"
            cmd="${DIR}grrun $gpname nsg"; ${DIR}execcmd.gap || exit 1
          fi;;
        esac
        cmdsource=
        cmd="${DIR}scrun $scflag1 $gpname sc0 ng"; ${DIR}execcmd.gap || exit 1
        suba=nsg; cflag2=-d0
    
        case $mult in
        false) mc3=${mc3}" ng"
               case $co in
               true) mc1=${mc1}" nsg"; mc2=${mc2}" -cr0";;
               esac;;
        esac;;
    
      false)
        suba=psg; cflag2=-d;;
    esac
    
    case $neqg in
    false)
        ocflag1=$cflag1
        if test ! "$repfile"
        then cflag1=${cflag1}" -c"
        fi
        
        oldno=;  scflag2=
        for no in 1 2 3 4 5 6 7 8 9
        do
           subb=g$no
           if test -r ${gpname}.$subb
           then
              cmd="${DIR}gprun -b $gpname $subb o$subb"; ${DIR}execcmd.gap || exit 1
              cmd="${DIR}egrun $gpname o$subb sg$no"; ${DIR}execcmd.gap || exit 1
              ${RM} -f ${gpname}.$subb; ${RM} -f ${gpname}.o$subb
              subb=sg$no
                
              cmd="${DIR}conrun\
                            $cflag1 $cflag2 $gpname $subb $suba pg dcr$no cr$no"
              ${DIR}execcmd.gap || exit 1
              cmd="${DIR}scrun $scflag1 $scflag2 $gpname sc$no dcr$no";
              ${DIR}execcmd.gap || exit 1
              oldno=$no; suba=$subb; cflag2=-d$oldno; scflag2=-s$oldno
             
              list=${list}" "${no}
              case $mult in
              false) mc3=${mc3}" dcr"$no
                     case $co in
                     true) mc1=${mc1}" sg"$no; mc2=${mc2}" -cr"$no;;
                     esac;;
              esac
              case $co in
              true) cmdsource="< ${gpname}.grip"
                cmd="${DIR}grrun $gpname $subb"; ${DIR}execcmd.gap || exit 1
                cmdsource=;;
              esac
           else break
           fi
        done
        
        cflag1=$ocflag1
        
        cmd="${DIR}conrun $cflag1 $cflag2 $gpname sg $suba";
        ${DIR}execcmd.gap || exit 1
        cmd="${DIR}scrun $scflag1 $scflag2 $gpname"; ${DIR}execcmd.gap || exit 1
        
        if test $check = true -a $mult = false
        then cmdsource="< ${gpname}.grip"
          cmd="${DIR}grrun $gpname sg"; ${DIR}execcmd.gap || exit 1
          cmdsource=
        fi;;
    esac
    
    case $mult in
    false) case $check in
           true) mflag=-t
                 case $neqg in
                 true) case $norm in
                       true) cp ${gpname}.nsg.rel ${gpname}.sg.rel;;
                       false)  cp ${gpname}.psg.rel ${gpname}.sg.rel;;
                       esac;;
                 esac;;
           false) mflag=;;
           esac
           cmd="${DIR}matcalc $mflag $gpname $mc1 $mc2 $mc3";
           ${DIR}execcmd.gap || exit 1
           case $first in
              true) nqcall="nqrun -g -1";;
              false) nqcall="nqrun -g";;
           esac;;
    true)  nqcall="nqmrun -g";;
    esac
    
    case $tc in
    true) case $co in
          true) case $neqg in
                true) case $norm in
                      true) rrarg=nsg;;
                      false) rrarg=psg;;
                      esac;;
                false) rrarg=;;
                esac
                cmd="${DIR}readrels -a $gpname $rrarg"; ${DIR}execcmd.gap || exit 1;;
          esac;;
    esac

    cmd="${DIR}$nqcall $gpname"; ${DIR}execcmd.gap || exit 1
    
    if test $mult = false -a $co = true
    then cmdsource=" < ${gpname}.nqip"
    fi
    case $norm in
    true) case $mult in
          true) nqarg=;;
          false) nqarg=ngmat;;
          esac
          case $neqg in
            true) case $co in
                    true) nqflag=-c;;
                    false) nqflag=;;
                  esac;;
            false) nqflag=;;
          esac
          cmd="${DIR}$nqcall $nqflag -a $gpname sc0 $nqarg";
          ${DIR}execcmd.gap || exit 1;;
    false) case $co in
           true) case $neqg in
                 true)  case $mult in
                          true) cmd="${DIR}$nqcall -c $gpname";
                                ${DIR}execcmd.gap || exit 1;;
                         false) cmd="${DIR}$nqcall -a -c $gpname $$"
                                ${DIR}execcmd.gap || exit 1;;
                        esac;;
                 esac;;
           esac;;
    esac
    
    case $neqg in
      false)
        for no in $list
        do case $mult in
           true) nqarg=;;
           false) nqarg=dcr${no}mat;;
           esac
           cmd="${DIR}$nqcall -a $gpname sc$no $nqarg"; ${DIR}execcmd.gap || exit 1
        done

        case $co in
        true) nqflag=-c;;
        false) nqflag=;;
        esac
        cmd="${DIR}$nqcall -a $nqflag $gpname"; ${DIR}execcmd.gap || exit 1;;
      esac
    cmdsource=;;
false)
     case $neqg in
     false)
        for no in 1 2 3 4 5 6 7 8 9
        do
         subb=sg$no
         if test -r ${gpname}.$subb
         then list=${list}" "$no
         else break
         fi
        done;;
    esac
    cmdsource=" < ${gpname}.nqip"
    cmd="${DIR}nqrun -g -a -c $gpname $$"; ${DIR}execcmd.gap || exit 1
    cmdsource=;;
esac
    
case $co in
true) case $mult in
      true) flag=-m; crarg=;;
      false) flag=;;
      esac
      echo "1.5 50 -10" > ${gpname}.erip

      if test $neqg = false -o $norm = true
      then
        cmdsource="< ${gpname}.erip"
        cmd="${DIR}extprun $flag $gpname psg"; ${DIR}execcmd.gap || exit 1
        cmdsource=
      fi

      case $norm in
      true) case $mult in
            false) crarg=cr0;;
            esac
            cmd="${DIR}crrun $flag $gpname nsg psg $crarg";
            ${DIR}execcmd.gap || exit 1
            case $neqg in
            false)
              cmdsource="< ${gpname}.erip"
              cmd="${DIR}extprun $flag $gpname nsg"; ${DIR}execcmd.gap || exit 1
              cmdsource=
              suba=nsg;;
            esac;;
      false) suba=psg;;
      esac

      case $neqg in
        false)
         for no in $list
         do case $mult in
            false) crarg=cr$no;;
            esac
            subb=sg$no
            cmd="${DIR}crrun $flag $gpname $subb $suba $crarg";
            ${DIR}execcmd.gap || exit 1
            cmdsource="< ${gpname}.erip"
            cmd="${DIR}extprun $flag $gpname $subb"; ${DIR}execcmd.gap || exit 1
            cmdsource=
            suba=$subb
         done

         case $mult in
         false) crarg=cr;;
         esac
         cmd="${DIR}crrun $flag $gpname sg $suba $crarg";
         ${DIR}execcmd.gap || exit 1;;
      esac

      case $tc in
      true) case $mult in
            true) rflag=-m;;
            false) rflag=;;
            esac
            case $neqg in
                true) case $norm in
                      true) rrarg=nsg;;
                      false) rrarg=psg;;
                      esac;;
                false) rrarg=;;
            esac
            cmd="${DIR}readrels -g $rflag $gpname $rrarg";
            ${DIR}execcmd.gap || exit 1;;
      esac;;
esac
