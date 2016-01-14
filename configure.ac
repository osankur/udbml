# AC_INIT (package, version, [bug-report], [tarname])
AC_INIT(udbml,0.1,[Maximilien.Colange@lsv.fr])

# We use automake to build makefiles
AM_INIT_AUTOMAKE([foreign -Wall, subdir-objects])

AC_CONFIG_SUBDIRS([uppaal-dbm/modules])

UPPAAL_INCLUDEDIR=-I${prefix}/uppaal/include
OCAML_INCLUDEDIR=`ocamlc -where`
UPPAAL_LIBDIR="-L${prefix}/uppaal/lib -ludbm"

AC_SUBST([UPPAAL_INCLUDEDIR])
AC_SUBST([OCAML_INCLUDEDIR])
AC_SUBST([UPPAAL_LIBDIR])

# the files to be produced
AC_CONFIG_FILES([Makefile uppaal-dbm/Makefile src/Makefile])

AC_PREFIX_DEFAULT(`ocamlfind printconf destdir`)

# Do it !! (mandatory)
AC_OUTPUT
