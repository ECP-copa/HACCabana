#include "Parameters.h"

namespace HACCabana
{

Parameters::Parameters() 
{
  ng = 256;
  np = 256;
  rL = 256.0;
  nsteps = 10;
  nsub = 3;
  rmax = 3.116326355;
  rsm = 0.01;
  z_in = 200.0;
  z_fin = 50.0;
  a_in = 1.0 / (1.0+z_in);
  a_fin = 1.0 / (1.0+z_fin);
  hubble = 0.6766;
  deut = 0.02242;
  Tcmb = 2.726; 
  omega_cdm = 0.26067;
  omega_baryon = deut / hubble / hubble;
  omega_nu = 0.0;
  omega_cb = omega_cdm + omega_baryon;
  omega_matter = omega_cb + omega_nu;
  omega_radiation = 2.471e-5*pow(Tcmb/2.725f,4.0f)/pow(hubble,2.0f);
  alpha = 1.0;
  neff_massless = 3.04;
  neff_massive = 0.0;
  f_nu_massless = neff_massless*7.0/8.0*pow(4.0/11.0,4.0/3.0);
  f_nu_massive = neff_massive*7.0/8.0*pow(4.0/11.0,4.0/3.0);
  w_de = -1.0;
  wa_de = 0.0;
  gpscal = static_cast< float >(ng) / static_cast< float >(np);
}

Parameters::~Parameters() 
{
  ;
}

void Parameters::load_from_file(std::string file_name)
{
  ;
}

}
