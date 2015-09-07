INCLUDEDIR=$(shell udbm-config --inc) -I$(shell ocamlc -where)
UPPAALDIR=$(shell udbm-config --dir)
main: main.ml
	ocamlfind ocamlopt -linkpkg -package udbml -I . -cclib -lstdc++ -cclib -std=c++11 -cclib -L. -cclib -L$(UPPAALDIR)/lib -cclib -ludbm -o main main.ml
