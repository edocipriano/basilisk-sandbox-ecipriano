include $(BASILISK)/Makefile.defs

CFLAGS += -O2 -I../src

aslamvof.tst: CFLAGS += -I$(BASILISKWIKI)
opensmoke.tst: CFLAGS += -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build
opensmoke-properties.tst: CFLAGS += -I$(OPENSMOKE_INTERFACE)/src -rpath $(OPENSMOKE_INTERFACE)/build
fsolve-gsl.tst: CFLAGS += -I/opt/homebrew/include -L/opt/homebrew/lib
boilingtemperature.tst: CFLAGS += -I/opt/homebrew/include -L/opt/homebrew/lib

# New targets
quick: aslam.tst aslamvof.tst fsolve-gsl.tst interfaceregression.tst \
  mapregion.tst radiation.tst spark.tst fickcorrected.tst \
  boilingtemperature.tst
slow: icentripetal.tst centripetal.tst
all: quick slow
