include $(BASILISK)/Makefile.defs

CFLAGS += -O2 \
	-disable-dimensions \
	-DMOVIE=1 \
	-DUSE_GSL=1 \
	-DUSE_LIBCONFIG=0 \
	-I$(BASILISKSANDBOX)/../dev \
	-I$(BREW_INCLUDE_PATHS) \
	-L$(BREW_LIBRARY_PATHS) \
	-I$(OPENSMOKE_INTERFACE_NEW)/src \
	-rpath $(OPENSMOKE_INTERFACE_NEW)/build

expansion-T375.c: expansion.c
	ln -s expansion.c expansion-T375.c
expansion-T375.s: CFLAGS += -DTEMPERATURE=375
expansion-T375.tst: CFLAGS += -DTEMPERATURE=375

expansion-T400.c: expansion.c
	ln -s expansion.c expansion-T400.c
expansion-T400.s: CFLAGS += -DTEMPERATURE=400
expansion-T400.tst: CFLAGS += -DTEMPERATURE=400

expansion.tst: expansion-T375.tst expansion-T400.tst

microgravity.tst: CFLAGS += -DTEMPERATURE=773 -DPRESSURE=10 \
	-DKINFOLDER=evaporation/n-heptane-in-nitrogen

normalgravity.tst: CFLAGS += -DTEMPERATURE=773 -DPRESSURE=10 \
	-DKINFOLDER=evaporation/n-heptane-in-nitrogen -DGRAVITY=-9.81

