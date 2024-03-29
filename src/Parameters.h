#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <map>

#include <Cabana_Core.hpp>

#include "Definitions.h"

namespace HACCabana
{

  class Parameters
  {
    public:

      // simulation and cosmology params
      int ng;
      int np;
      float rL;
      float oL;
      int nsteps;
      int nsub;
      float rmax;
      float rsm;
      float z_in;
      float z_fin;
      float a_in;
      float a_fin;
      float hubble;
      float deut;
      float Tcmb; 
      float omega_cdm;
      float omega_baryon;
      float omega_nu;
      float omega_cb;
      float omega_matter;
      float omega_radiation;
      float alpha;
      float neff_massless;
      float neff_massive;
      float f_nu_massless;
      float f_nu_massive;
      float w_de;
      float wa_de;
      float gpscal;
      float cm_size;

      Parameters();
      ~Parameters();

      void load_from_file(std::string file_name);

    private:
      std::map<std::string, std::string> m_params;
  };

}
#endif
