UPPAALINCLUDEDIR=$(shell udbm-config --inc)
UPPAALLIB=$(foreach i,$(shell udbm-config --libs),-cclib $i)

main.native:
	ocamlbuild -use-ocamlfind -cflags '-cclib $(UPPAALINCLUDEDIR)' -lflags '-cclib -lstdc++ -cclib -std=c++11 $(UPPAALLIB)' main.native 
main.byte:
	ocamlbuild -use-ocamlfind -cflags '-cclib $(UPPAALINCLUDEDIR)' -lflags '-cclib -lstdc++ -cclib -std=c++11 -cclib $(UPPAALLIB)' main.byte
clean:
	ocamlbuild -clean

