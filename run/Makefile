include $(BASILISK)/Makefile.defs

BASILISK_INCLUDES ?= -I../src
BASILISK_LIBS ?= 

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -disable-dimensions $(BASILISK_INCLUDES) $(BASILISK_LIBS)

fixedbubblevelocity-jump.c: fixedbubblevelocity.c
	ln -s fixedbubblevelocity.c fixedbubblevelocity-jump.c
fixedbubblevelocity-jump.s: CFLAGS += -DJUMP=1
fixedbubblevelocity-jump.tst: CFLAGS += -DJUMP=1

fixedbubblevelocity.tst: fixedbubblevelocity-jump.tst
