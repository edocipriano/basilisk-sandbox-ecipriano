#define NGS 20
#define NLS 1

char* gas_species[NGS] = { 
  "CO2",
  "CH4",
  "CH3",
  "CH2",
  "CH2(S)",
  "CH3OH",
  "CH3O",
  "CH2OH",
  "CH2O",
  "HCO",
  "N2",
  "H2",
  "H",
  "O2",
  "O",
  "H2O",
  "OH",
  "H2O2",
  "HO2",
  "CO"
};

char* liq_species[NLS] = { 
  "CH3OH"
};

char* inert_species[1] = { 
  "N2"
};

double gas_start[NGS] = { 
  0., 
  0., 
  0., 
  0., 
  0., 
  0., 
  0., 
  0., 
  0., 
  0., 
  0.7670907862,
  0., 
  0., 
  0.2329092138,
  0., 
  0., 
  0., 
  0., 
  0., 
  0.  
};

double liq_start[NLS] = {
  1.
};

double inDmix2[NGS] = {
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5,
  1.64e-5
};

double inDmix1[NLS] = {0.};

double inKeq[NLS] = {0.7};

double Tboil[NLS] = {337.85};

double lambda1 = 0.;

double lambda2 = 0.;

double dhev = 0.;

double cp1 = 0.;

double cp2 = 0.;

