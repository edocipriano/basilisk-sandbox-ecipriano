include $(BASILISK)/Makefile.defs

CFLAGS += -O2 -I../src -disable-dimensions

aslamvof.tst: CFLAGS += -I$(BASILISKWIKI)
opensmoke.tst: CFLAGS += -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build
opensmoke-properties.tst: CFLAGS += -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build
batchreactor.tst: CFLAGS += -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build
fsolve-gsl.tst: CFLAGS += -I/opt/homebrew/include -L/opt/homebrew/lib
boilingtemperature.tst: CFLAGS += -I/opt/homebrew/include -L/opt/homebrew/lib
poisson-gfm.tst: CFLAGS += -DGFM=1
capwave.tst: CFLAGS += -DGFM=1
spurious.tst: CFLAGS += -DGFM=1
oscillation.tst: CFLAGS += -DGFM=1

pinning-axi.c: pinning.c
	ln -s pinning.c pinning-axi.c
pinning-axi.s: CFLAGS += -DAXI=1
pinning-axi.tst: CFLAGS += -DAXI=1

pinning.tst: pinning-axi.tst

# New targets
GFM = poisson-gfm.tst capwave.tst spurious.tst oscillation.tst
quick: aslam.tst aslamvof.tst fsolve-gsl.tst interfaceregression.tst \
  mapregion.tst radiation.tst spark.tst fickcorrected.tst \
  boilingtemperature.tst $(GFM)
gfm: $(GFM)
slow: icentripetal.tst centripetal.tst pinning.tst
all: quick slow
