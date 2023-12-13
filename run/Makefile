include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -disable-dimensions -I../src -I$(BASILISKWIKI) -I$(BREW_INCLUDE_PATHS) -L$(BREW_LIBRARY_PATHS) -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build

dropdilation.tst: CFLAGS += -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build
bubblethermalexpansion.tst: CFLAGS += -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build
bubblethermalshrinking.tst: CFLAGS += -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build

fixedbubblevelocity-jump.c: fixedbubblevelocity.c
	ln -s fixedbubblevelocity.c fixedbubblevelocity-jump.c
fixedbubblevelocity-jump.s: CFLAGS += -DJUMP=1
fixedbubblevelocity-jump.tst: CFLAGS += -DJUMP=1

fixedbubblevelocity.tst: fixedbubblevelocity-jump.tst
