include $(BASILISK)/Makefile.defs

CFLAGS += -O2 -I../src

centripetal.tst: CC = mpicc -D_MPI=4
icentripetal.tst: CC = mpicc -D_MPI=4
opensmoke.tst: CFLAGS += -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build
opensmoke-properties.tst: CFLAGS += -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build
