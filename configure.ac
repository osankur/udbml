# AC_INIT (package, version, [bug-report], [tarname])
AC_INIT(udbml,0.2,[Maximilien.Colange@lsv.fr])

# We use automake to build makefiles
AM_INIT_AUTOMAKE([foreign -Wall, subdir-objects])

AC_CONFIG_SUBDIRS([uppaal-dbm/modules])

OCAML_INCLUDEDIR=`ocamlc -where`

AC_SUBST([OCAML_INCLUDEDIR])

AC_PROG_CC
AC_PROG_CXX

# Use libtool
AM_PROG_AR
LT_INIT

# the files to be produced
AC_CONFIG_FILES([Makefile uppaal-dbm/Makefile uppaal-dbm/modules/lib/Makefile src/Makefile])

AC_PREFIX_DEFAULT(`ocamlfind printconf destdir`)

# Do it !! (mandatory)
AC_OUTPUT
