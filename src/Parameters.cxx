#include <fstream>
#include <sstream> 
#include "Parameters.h"

namespace HACCabana
{

Parameters::Parameters() :
  ng(0),
  np(0),
  rL(0.0f),
  oL(0.0f),
  nsteps(0),
  nsub(0),
  rmax(3.116326355),
  rsm(0.0),
  z_in(0.0),
  z_fin(0.0),
  hubble(0.0),
  deut(0.0),
  Tcmb(0.0),
  omega_cdm(0.0),
  omega_nu(0.0),
  alpha(0.0),
  neff_massless(0.0),
  neff_massive(0.0),
  w_de(0.0),
  wa_de(0.0),
  cm_size(0.0)
{
  ;
}

Parameters::~Parameters() 
{
  ;
}

void Parameters::load_from_file(std::string file_name)
{
  m_params.clear();

  std::ifstream file(file_name);
  std::stringstream filestream;

  // don't do this when using MPI
  if (file) 
  {
    filestream << file.rdbuf();
    file.close();
  }

  while (filestream.good()) 
  {
    std::string line;
    getline(filestream, line);
    size_t pos = line.find_first_not_of(" \t");
    if (pos < line.size())
      line = line.substr(pos);

    if (line.empty()) 
    {
      continue;
    } else if (line[0] == '#') {
      continue;
    }

    pos = line.find_first_of(" \t");
    std::string key = line.substr(0, pos);

    if (pos < line.size())
      line = line.substr(pos);
    else
      line = "";

    pos = line.find_first_not_of(" \t");
    if (pos < line.size())
      line = line.substr(pos);
    else
      line = "";

    std::string value = line;
    m_params.insert(std::make_pair(key, value));
  }

  if (!m_params.count("HACC_HEADER_VERSION")) 
  {
    throw std::runtime_error( "HACC_HEADER_VERSION not found in parameters file!" );
  } 
  else if (m_params["HACC_HEADER_VERSION"] != "1.0.0") 
  {
     throw std::runtime_error( "Only HACC_HEADER_VERSION 1.0.0 is understood!" );
  }

  ng = atoi(m_params["NG"].c_str());
  np = atoi(m_params["NP"].c_str());
  rL = atof(m_params["RL"].c_str());
  oL = atof(m_params["OL"].c_str());
  nsteps = atoi(m_params["N_STEPS"].c_str());
  nsub = atoi(m_params["N_SUB"].c_str());
  //rmax(3.116326355),
  rsm = atof(m_params["RSM"].c_str());
  z_in = atof(m_params["Z_IN"].c_str());
  z_fin = atof(m_params["Z_FIN"].c_str());
  hubble = atof(m_params["HUBBLE"].c_str());
  deut = atof(m_params["DEUT"].c_str());
  Tcmb = atof(m_params["T_CMB"].c_str());
  omega_cdm = atof(m_params["OMEGA_CDM"].c_str());
  omega_nu = atof(m_params["OMEGA_NU"].c_str());
  alpha = atof(m_params["ALPHA"].c_str()); 
  neff_massless = atof(m_params["N_EFF_MASSLESS"].c_str());
  neff_massive = atof(m_params["N_EFF_MASSIVE"].c_str());
  w_de = atof(m_params["W_DE"].c_str());
  wa_de = atof(m_params["WA_DE"].c_str());
  cm_size = atof(m_params["CM_SIZE"].c_str());

  a_in = 1.0 / (1.0+z_in);
  a_fin = 1.0 / (1.0+z_fin);
  omega_baryon = deut / hubble / hubble;
  omega_cb = omega_cdm + omega_baryon;
  omega_matter = omega_cb + omega_nu;
  omega_radiation = 2.471e-5*pow(Tcmb/2.725f,4.0f)/pow(hubble,2.0f);
  f_nu_massless = neff_massless*7.0/8.0*pow(4.0/11.0,4.0/3.0);
  f_nu_massive = neff_massive*7.0/8.0*pow(4.0/11.0,4.0/3.0);
  gpscal = static_cast< float >(ng) / static_cast< float >(np);
}

}
