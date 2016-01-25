# -*- mode: autoconf; -*-
###############################################################################
#
# Filename: aclocal.m4
#
# This file is a part of the UPPAAL toolkit.
# Copyright (c) 1995 - 2000, Uppsala University and Aalborg University.
# All right reserved.
#
# Local configure macros for UPPAAL.
#
# $Id: acinclude.m4,v 1.3 2004/08/26 10:38:54 behrmann Exp $
#
###############################################################################


AC_DEFUN([UA_CXX_STREAMBUF],
[
  AC_CACHE_CHECK([flavour of std::streambuf], ua_cv_cxx_streambuf, 
  [
    AC_LANG_SAVE
    AC_LANG_CPLUSPLUS
    cat > conftest.$ac_ext <<EOF
#include <iostream>
using std::streambuf;

class Test : public std::streambuf
{
   void test(std::streambuf *a) { a->sync(); }
};
EOF

    if AC_TRY_EVAL(ac_compile); then
      ua_cv_cxx_streambuf="old"
    else
      ua_cv_cxx_streambuf="new"
    fi
    AC_LANG_RESTORE
  ])

  if test $ua_cv_cxx_streambuf = old; then
    AC_DEFINE(CXX_OLD_STREAMBUF,1,[Define if the library implementation of streambuf have public interface to override])
  fi
])
