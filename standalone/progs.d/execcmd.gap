case $step in
  true) echo -n $cmd $cmdsource $cmddest"   ?"
        read ans
        case $ans in
          n*) exit 0;;
          q*) exit 1;;
        esac;;
  *)    case $verbose in
            true) echo RUNNING: $cmd $cmdsource $cmddest
        esac;;
esac

eval $cmd $cmdsource $cmddest || exit 1

exit 0
