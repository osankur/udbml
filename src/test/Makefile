INCLUDEDIR=$(shell udbm-config --inc) -I$(shell ocamlc -where)
UPPAALLIB=$(shell udbm-config --libs)
main: main.ml
	OCAMLPATH=/home/osankur/lrde/tiamo/local ocamlfind ocamlopt -linkpkg -package udbml -I . -I $(INCLUDEDIR) -cclib -lstdc++ -cclib -std=c++11 -cclib -L. -cclib "$(UPPAALLIB)" -o main main.ml
