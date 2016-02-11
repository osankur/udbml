libudbm.a: libbase.a libdbm.a libdebug.a libhash.a
	mkdir tmp-lib && \
	cd tmp-lib && \
	ar x ../libbase.a && \
	ar x ../libdbm.a && \
	ar x ../libdebug.a && \
	ar x ../libhash.a && \
	cd .. && \
	ar rucs libudbm.a tmp-lib/*.o && \
	rm -rf tmp-lib
	
