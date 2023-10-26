include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -dimensions -I../src -I$(BASILISKWIKI) -I$(BREW_INCLUDE_PATHS) -L$(BREW_LIBRARY_PATHS)

dropdilation.tst: CFLAGS += -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build
bubblethermalexpansion.tst: CFLAGS += -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build
bubblethermalshrinking.tst: CFLAGS += -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build

fixedbubblevelocity-extrapolation.c: fixedbubblevelocity.c
	ln -s fixedbubblevelocity.c fixedbubblevelocity-extrapolation.c
fixedbubblevelocity-extrapolation.s: CFLAGS += -DEXTRAPOLATION=1
fixedbubblevelocity-extrapolation.tst: CFLAGS += -DEXTRAPOLATION=1

fixedbubblevelocity.tst: fixedbubblevelocity-extrapolation.tst
