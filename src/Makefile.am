UPPAAL_INCLUDES=@UPPAAL_INCLUDEDIR@
OCAML_INCLUDES=-I@OCAML_INCLUDEDIR@

UDBML_INCLUDES=$(foreach i,$(UPPAAL_INCLUDES) $(OCAML_INCLUDEDIR),-ccopt $i)

UPPAALLIB=@UPPAAL_LIBDIR@

PACKAGE_FILES=libudbml.a udbml.a udbml.cmi udbml.cmxa
INSTALL_FILES=$(foreach i,$(PACKAGE_FILES),$(BUILDDIR)/$i)

BUILDDIR=_build

OCAMLBUILD=ocamlbuild -use-ocamlfind -classic-display -build-dir $(BUILDDIR)

UDBML_PREFIX=${prefix}

CXX_FLAGS=-fPIC -O3 -DNDEBUG
OCAML_CXX_FLAGS=$(foreach i,$(CXX_FLAGS),-ccopt $i)

all-local:
	$(OCAMLBUILD) -cflags "-cc g++ $(UDBML_INCLUDES) $(OCAML_CXX_FLAGS)" udbm_stubs.o
	$(OCAMLBUILD) udbml.cmx
	ocamlmklib -custom -o $(BUILDDIR)/udbml $(BUILDDIR)/udbml.cmx $(BUILDDIR)/udbm_stubs.o $(UPPAALLIB) -lstdc++

install-data-local:
	ocamlfind install -destdir $(UDBML_PREFIX) udbml META $(INSTALL_FILES)
uninstall-local:
	ocamlfind remove -destdir $(UDBML_PREFIX) udbml
clean-local:
	ocamlbuild -clean
reinstall:
	@$(MAKE) uninstall
	@$(MAKE) install
