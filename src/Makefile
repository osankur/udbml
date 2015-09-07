INCLUDEDIR=$(shell udbm-config --inc) -I$(shell ocamlc -where)
UPPAALDIR=$(shell udbml-config --dir)
PACKAGE_FILES=libudbml.a udbml.a udbml.cmi udbml.cmxa
main: udbml.cmxa main.ml
	ocamlopt -I . -cclib -lstdc++ -cclib -std=c++11 -cclib -L. -cclib $(UPPAALDIR)/lib -cclib -ludbm udbml.cmxa -o main main.ml
udbm_stubs.o: udbm_stubs.c
	g++ -o udbm_stubs.o -std=c++11 $(INCLUDEDIR) -fPIC -c udbm_stubs.c
udbml.mli: udbml.ml
	ocamlc -i $< > $@
udbml.cmx: udbml.ml
	ocamlopt -c $<
libudbml.a dlludbml.so: udbm_stubs.o
	ocamlmklib -v -o udbml $<
udbml.cmxa: udbml.cmx libudbml.a
	ocamlmklib -v -o udbml udbml.cmx
clean:
	rm -f *.[oa] *.so *.cm[ixoa] *.cmxa
install:
	@ocamlfind install udbml META $(PACKAGE_FILES)
uninstall:
	@ocamlfind remove udbml
reinstall:
	@ocamlfind remove udbml
	@ocamlfind install udbml $(PACKAGE_FILES)
