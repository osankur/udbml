# AC_INIT (package, version, [bug-report], [tarname])
AC_INIT(udbml,0.2,[Maximilien.Colange@lsv.fr])

# We use automake to build makefiles
AM_INIT_AUTOMAKE([foreign -Wall, subdir-objects])

AC_CONFIG_SUBDIRS([uppaal-dbm/modules])

OCAML_INCLUDEDIR=`ocamlc -where`

AC_SUBST([OCAML_INCLUDEDIR])

# If the user didn't supply a CFLAGS value,
# set an empty one to prevent autoconf to stick -O2 -g here.
test -z "$CFLAGS" && CFLAGS=
test -z "$CXXFLAGS" && CXXFLAGS=

# We compile c11/c++11
CFLAGS="-std=c11 $CFLAGS"
CXXFLAGS="-std=c++11 $CXXFLAGS"

AC_ARG_ENABLE([debug],
              [AC_HELP_STRING([--enable-debug], [turn on debug build])])
if test "$enable_debug" = yes; then
  # enable debug symbols, and disable optimizations
  CFLAGS="-g -O0 $CFLAGS"
  CXXFLAGS="-g -O0 $CXXFLAGS"
else
  CFLAGS="-O3 $CFLAGS"
  CXXFLAGS="-O3 $CXXFLAGS"
fi

AC_ARG_ENABLE([assert],
              [AC_HELP_STRING([--enable-assert], [turn on assertions])])
if test "$enable_assert" != yes; then
  CPPFLAGS="$CPPFLAGS -DNDEBUG"
fi

AC_PROG_CC
AC_PROG_CXX

# the files to be produced
AC_CONFIG_FILES([Makefile uppaal-dbm/Makefile src/Makefile])

AC_PREFIX_DEFAULT(`ocamlfind printconf destdir`)

# Do it !! (mandatory)
AC_OUTPUT
