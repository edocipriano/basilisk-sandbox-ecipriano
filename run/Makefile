include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -I../src -I$(BASILISKWIKI) -I$(BREW_INCLUDE_PATHS) -L$(BREW_LIBRARY_PATHS)

fixedbubblevelocity-extrapolation.c: fixedbubblevelocity.c
	ln -s fixedbubblevelocity.c fixedbubblevelocity-extrapolation.c
fixedbubblevelocity-extrapolation.s: CFLAGS += -DEXTRAPOLATION=1
fixedbubblevelocity-extrapolation.tst: CFLAGS += -DEXTRAPOLATION=1

fixedbubblevelocity.tst: fixedbubblevelocity-extrapolation.tst
