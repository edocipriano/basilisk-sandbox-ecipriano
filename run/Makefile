include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -events -I../src -I$(BREW_INCLUDE_PATHS) -L$(BREW_LIBRARY_PATHS) -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build

c7pathak-varprop-binary.tst: CC = mpicc -D_MPI=2
