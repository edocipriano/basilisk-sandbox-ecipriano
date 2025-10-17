include $(BASILISK)/Makefile.defs

CFLAGS += -O2 \
	-disable-dimensions \
	-DMOVIE=1 \
	-DUSE_GSL=1 \
	-DMTRACE=2 \
	$(BASILISK_INCLUDES) \
	$(BASILISK_LIBS)

poisson-gfm.tst: CFLAGS += -DGFM=1
capwave.tst: CFLAGS += -DGFM=1
spurious.tst: CFLAGS += -DGFM=1
oscillation.tst: CFLAGS += -DGFM=1

pinning-axi.c: pinning.c
	ln -s pinning.c pinning-axi.c
pinning-axi.s: CFLAGS += -DAXI=1
pinning-axi.tst: CFLAGS += -DAXI=1

pinning.tst: pinning-axi.tst

low-mach.tst: CFLAGS += -DUSE_VOFI=1 -I$(VOFI_LIBRARY_PATH)/include
