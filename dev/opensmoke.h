/**
# OpenSMOKE++

Basilisk interface for the [OpenSMOKE++ library](#cuoci2015opensmoke++),
which is a C++ framework for numerical simulations of reacting systems
with detailed kinetic mechanisms, including thousands of chemical
species and reactions. It provides the functionalities for managing
the calculation of thermodynamics and transport properties of mixtures
and for handling the chemical kinetics.

## OpenSMOKEInterface
This module requires an external library, which contains both the
OpenSMOKE++ library itself, and the C interface which contains the
functions that can be called from Basilisk. It can be downloaded
from [github](https://github.com/edocipriano/OpenSMOKEppInterface)
and it can be easily compiled using cmake, provided the following
external dependencies: 
[OpenSMOKE++](#cuoci2015opensmoke++), 
[boost](https://www.boost.org/), 
[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page), 
[libconfig++](https://github.com/hyperrealm/libconfig).
*/

#include "OpenSMOKE_Interface.h"
#define OPENSMOKE 1

#pragma autolink -L$OPENSMOKE_INTERFACE/build -lopensmoke

/**
## User Data

The variable *kinfolder*, with the path of the kinetic scheme folder
relative to the path $OPENSMOKE_INTERFACE/kinetics, must be set by the
user.
*/

char * kinfolder, * liqfolder;

void molecular_weights (size_t ns, double * MW) {
  for (int i = 0; i < ns; i++)
    MW[i] = OpenSMOKE_MW (i);
}

int index_species (char * name) {
  return OpenSMOKE_IndexOfSpecies (name);
}

void species_names (size_t ns, char ** names) {
  for (int i = 0; i < ns; i++)
    names[i] = strdup (OpenSMOKE_NamesOfSpecies (i));
}

void species_names_liquid (size_t ns, char ** names) {
  for (int i = 0; i < ns; i++) {
    const char * species = OpenSMOKE_NamesOfLiquidSpecies (i);
    int len = strlen (species);
    char corrname[len+1];
    strcpy (corrname, species);
    corrname[3 <= len ? len-3 : 0] = '\0';
    names[i] = strdup (corrname);
  }
}

#define MAX_KINFOLDER_LEN 1200

void kinetics (char * kinfolder) {
  char kinfolder_root[MAX_KINFOLDER_LEN];
  sprintf (kinfolder_root, "%s/kinetics/%s/kinetics",
      getenv ("OPENSMOKE_INTERFACE"), kinfolder);

  OpenSMOKE_Init();
  OpenSMOKE_ReadKinetics (kinfolder_root);
}

void kinetics_liquid (char * kinfolder) {
  char kinfolder_root[MAX_KINFOLDER_LEN];
  sprintf (kinfolder_root, "%s/kinetics/%s/kinetics",
      getenv ("OPENSMOKE_INTERFACE"), kinfolder);

  OpenSMOKE_ReadLiquidKinetics (kinfolder_root);
}

void properties_liquid (char * liqfolder) {
  char liqfolder_root[MAX_KINFOLDER_LEN];
  sprintf (liqfolder_root, "%s/kinetics/LiquidProperties/%s",
      getenv ("OPENSMOKE_INTERFACE"), liqfolder);

  OpenSMOKE_ReadLiquidProperties (liqfolder_root);
}

event cleanup (t = end)
{
  OpenSMOKE_Clean ();
}

/**
## References

~~~bib
@article{cuoci2015opensmoke++,
  title={OpenSMOKE++: An object-oriented framework for the numerical modeling of reactive systems with detailed kinetic mechanisms},
  author={Cuoci, Alberto and Frassoldati, Alessio and Faravelli, Tiziano and Ranzi, Eliseo},
  journal={Computer Physics Communications},
  volume={192},
  pages={237--264},
  year={2015},
  publisher={Elsevier}
}
~~~
*/

