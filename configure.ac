# AC_INIT (package, version, [bug-report], [tarname])
AC_INIT(udbml,0.1,[Maximilien.Colange@lsv.fr])

# We use automake to build makefiles
AM_INIT_AUTOMAKE([foreign -Wall, subdir-objects])

UPPAAL_INCLUDEDIR=`udbm-config --inc`
OCAML_INCLUDEDIR=`ocamlc -where`
UPPAAL_LIBDIR=`udbm-config --libs`

AC_SUBST([UPPAAL_INCLUDEDIR])
AC_SUBST([OCAML_INCLUDEDIR])
AC_SUBST([UPPAAL_LIBDIR])

# the files to be produced
AC_CONFIG_FILES([Makefile src/Makefile])

AC_PREFIX_DEFAULT(`ocamlfind printconf destdir`)

# Do it !! (mandatory)
AC_OUTPUT
