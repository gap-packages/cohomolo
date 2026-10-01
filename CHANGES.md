This file describes changes in the cohomolo package.

1.7.1 (2026-08-31)
  - In SplitExtensionCHR, NonsplitExtension, and CoveringGroup, the
    order of the result is known if the order of the given finitely
    presented group is known

1.7.0 (2026-08-13)
  - Report an error when an external program fails, instead of returning the
    result of an earlier stage of the computation
  - Raise the space `nqrun` reserves for relations tenfold
  - Replace the recursive build system by a single Makefile, with new
    `standalone`, `doc` and `check` targets
  - Fix the type mismatches between the C files, and make the
    standalone-only programs compile again
  - Convert the documentation to GAPDoc/AutoDoc, and fix stale parts of it
  - Various janitorial changes

1.6.12 (2025-10-17)
  - Various janitorial changes

1.6.11 (2023-01-03)
  - Enhance the `configure` script to accept `--with-gaproot=PATH`
  - Various janitorial changes

1.6.10 (2022-03-30)
  - Fix more compiler warnings
  - Fix building the manual in latest GAP

1.6.9 (2021-02-19)
  - Fix linker error when using GCC 10

1.6.8 (2019-07-07)
  - Allow CC and CFLAGS to be overridden
  - Fix a bunch of compiler warning
  - Simplify parts of the build system
  - PackageInfo.g: specify license in machine readable form

1.6.7 (2018-09-14)
  - Set the GAP team as new package maintainer

1.6.6 (2017-12-02)
  - Internal changes (make test suite compatible again with GAP 4.8)

1.6.5 (2017-11-08)
  - Fixed a bug were using GAP in a narrow terminal could lead
    to a communication problem between the GAP package and the
    standalone C code.

1.6.4 (2016-01-20)
  - No user visible changes (some test files were added)

1.6.3 (2016-01-16)
  - Update README
  - Clarify that the package is licensed under the GPL v2 or later
  - Workaround a miscompilation bug in Apple clang-700.1.81

1.6.2 (2016-01-15)
  - Add Markus Pfeiffer and Max Horn as maintainers
  - Convert C functions from K&R to ANSI C

1.6.1 (2016-01-14)
  - Moved package homepage to https://gap-packages.github.io/cohomolo/
  - Added GPL v2 license
  - Improve build system
  - Fix many warnings in the C code
  - Add HTML version of the manual
