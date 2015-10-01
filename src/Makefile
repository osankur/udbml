UPPAAL_INCLUDEDIR=$(shell udbm-config --inc)
OCAML_INCLUDEDIR=-I$(shell ocamlc -where)

INCLUDEDIR=$(foreach i,$(UPPAAL_INCLUDEDIR) $(OCAML_INCLUDEDIR),-ccopt $i)

UPPAALLIB=$(shell udbm-config --libs)

PACKAGE_FILES=libudbml.a udbml.a udbml.cmi udbml.cmxa
INSTALL_FILES=$(foreach i,$(PACKAGE_FILES),$(BUILDDIR)/$i)

BUILDDIR=_build

OCAMLBUILD=ocamlbuild -use-ocamlfind -classic-display -build-dir $(BUILDDIR)

all:
	$(OCAMLBUILD) -cflags "-cc g++ $(INCLUDEDIR)" udbm_stubs.o
	$(OCAMLBUILD) udbml.cmx
	ocamlmklib -o $(BUILDDIR)/udbml $(BUILDDIR)/udbml.cmx $(BUILDDIR)/udbm_stubs.o $(UPPAALLIB) -lstdc++

install:
	ocamlfind install udbml META $(INSTALL_FILES)
uninstall:
	ocamlfind remove udbml
clean:
	ocamlbuild -clean
reinstall:
	@ocamlfind remove udbml
	@ocamlfind install udbml META $(PACKAGE_FILES)
