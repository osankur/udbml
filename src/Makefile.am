UPPAAL_INCLUDES=-I@abs_top_srcdir@/uppaal-dbm/modules/include
UDBML_INCLUDES=$(foreach i,$(UPPAAL_INCLUDES) $(OCAML_INCLUDEDIR),-ccopt $i)
OCAML_INCLUDES=-I@OCAML_INCLUDEDIR@

UPPAALLIB=-L${prefix}/uppaal/lib -ludbm

MODULES=udbml_unpriced udbml_priced udbml
MODULES_FILES=$(MODULES:=.cmi) $(MODULES:=.cmx)
PACKAGE_FILES=libudbml.a udbml.a udbml.cmxa $(MODULES_FILES)

INSTALL_FILES=$(foreach i,$(PACKAGE_FILES),$(BUILDDIR)/$i) udbm_stubs.h

BUILDDIR=_build

OCAMLBUILD=ocamlbuild -use-ocamlfind -classic-display -build-dir $(BUILDDIR)

UDBML_PREFIX=${prefix}

CXX_FLAGS=-xc++ -fPIC -I@abs_top_srcdir@/src $(CPPFLAGS) $(CXXFLAGS)
OCAML_CXX_FLAGS=$(foreach i,$(CXX_FLAGS),-ccopt $i)

all-local:
	$(OCAMLBUILD) -cflags "-cc g++ $(UDBML_INCLUDES) $(OCAML_CXX_FLAGS)" udbm_stubs.o
	$(OCAMLBUILD) -cflags "-cc g++ $(UDBML_INCLUDES) $(OCAML_CXX_FLAGS)" updbm_stubs.o
	$(OCAMLBUILD) udbml_unpriced.cmx
	$(OCAMLBUILD) udbml_priced.cmx
	ocamlmklib -custom -o $(BUILDDIR)/udbml $(BUILDDIR)/udbml.cmx $(BUILDDIR)/udbml_unpriced.cmx $(BUILDDIR)/udbm_stubs.o $(BUILDDIR)/udbml_priced.cmx $(BUILDDIR)/updbm_stubs.o $(UPPAALLIB) -lstdc++

install-data-local:
	ocamlfind install -destdir $(UDBML_PREFIX) udbml META $(INSTALL_FILES)
uninstall-local:
	ocamlfind remove -destdir $(UDBML_PREFIX) udbml
clean-local:
	ocamlbuild -clean
reinstall:
	@$(MAKE) uninstall
	@$(MAKE) install
