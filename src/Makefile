INCLUDEDIR=$(shell udbm-config --inc) -I$(shell ocamlc -where)
UPPAALLIB=$(shell udbm-config --libs)
PACKAGE_FILES=libudbml.a udbml.a udbml.cmi udbml.cmxa udbml.cma
main: main.ml udbml.cmxa
	ocamlopt -I . -cclib -lstdc++ -cclib -std=c++11 -cclib -L. -cclib "$(UPPAALLIB)" udbml.cmxa -o main main.ml
udbml.cmxa udbml.cma: udbml.cmx udbml.cmo libudbml.a
	ocamlmklib -v -o udbml udbml.cmx
	ocamlmklib -v -o udbml udbml.cmo
udbm_stubs.o: udbm_stubs.cc
	g++ -o udbm_stubs.o -std=c++11 $(INCLUDEDIR) -fPIC -c udbm_stubs.cc
udbml.mli: udbml.ml
	ocamlc -i $< > $@
udbml.cmx: udbml.ml
	ocamlopt -c $<
udbml.cmo: udbml.ml
	ocamlc -c $<
libudbml.a dlludbml.so: udbm_stubs.o
	ocamlmklib -v -o udbml $<
clean:
	rm -f *.[oa] *.so *.cm[ixoa] *.cmxa
install:
	@ocamlfind install udbml META $(PACKAGE_FILES)
uninstall:
	@ocamlfind remove udbml
reinstall:
	@ocamlfind remove udbml
	@ocamlfind install udbml META $(PACKAGE_FILES)
