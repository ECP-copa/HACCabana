#include <iostream>
#include <Cabana_Core.hpp>
#include <Cabana_AoSoA.hpp>
#include <Cabana_LinkedCellList.hpp>

#include "Definitions.h"
#include "TimeStepper.h"
#include "Particles.h"
#include "ParticleActions.h"

using namespace std;

inline bool floatCompare(float f1, float f2) {
  static constexpr auto epsilon = MAX_ERR;
  if (fabs(f1 - f2) <= epsilon)
    return true;
  return fabs(f1 - f2) <= epsilon * fmax(fabs(f1), fabs(f2));
}

int main( int argc, char* argv[] )
{

  // Kokkos::ScopeGuard initializes Kokkos and guarantees it is finalized.
  Kokkos::ScopeGuard scope_guard(argc, argv);

  if (argc != 4)
  {
    cout << "Usage: ./driver_sr <particles before subcycle>.bin <particles after subcycle>.bin step_number" << endl;
    exit(1);
  }

  // simulation and cosmology params
  const int ng = 256;
  const int np = 256;
  const int nsteps = 10;
  const int nsub = 3;
  const float rmax = 3.116326355;
  const float rsm = 0.01;
  const float z_in = 200.0;
  const float z_fin = 50.0;
  const float a_in = 1.0 / (1.0+z_in);
  const float a_fin = 1.0 / (1.0+z_fin);
  const float hubble = 0.6766;
  const float deut = 0.02242;
  const float Tcmb = 2.726; 
  const float omega_cdm = 0.26067;
  const float omega_baryon = deut / hubble / hubble;
  const float omega_nu = 0.0;
  const float omega_cb = omega_cdm + omega_baryon;
  const float omega_matter = omega_cb + omega_nu;
  const float omega_radiation = 2.471e-5*pow(Tcmb/2.725f,4.0f)/pow(hubble,2.0f);
  const float alpha = 1.0;
  const float neff_massless = 3.04;
  const float neff_massive = 0.0;
  const float f_nu_massless = neff_massless*7.0/8.0*pow(4.0/11.0,4.0/3.0);
  const float f_nu_massive = neff_massive*7.0/8.0*pow(4.0/11.0,4.0/3.0);
  const float w_de = -1.0;
  const float wa_de = 0.0;
  const float gpscal = static_cast< float >(ng) / static_cast< float >(np);

  TimeStepper ts(
     alpha,
		 a_in,
		 a_fin,
		 nsteps,
		 omega_matter,
		 omega_cdm,
		 omega_baryon,
		 omega_cb,
		 omega_nu,
		 omega_radiation,
		 f_nu_massless,
		 f_nu_massive,
		 w_de,
     wa_de);
  
  const int step0 = atoi(argv[3]);  // the full simulation step number we are subcycling

  cout << "Advancing timestepper to step " << step0 << endl;

  //get timestepper up to speed
  for(int step = 0; step < step0; step++)
    ts.advanceFullStep();

  // we're starting to subcycle after a PM kick
  ts.advanceHalfStep();

  HACCabana::Particles P;
  P.readRawData(argv[1]);
  cout << "Finished reading file: " << argv[1] << endl;
  
  HACCabana::ParticleActions PA(&P);
  PA.subCycle(ts, nsub, gpscal, rmax*rmax, rsm*rsm);

  // verify against the answer from the simulation
  // --------------------------------------------------------------------------------------------------------------------------

  if (true)
  {
    cout << "\nVerifying result." << endl;
    HACCabana::Particles P_ans;
    P_ans.readRawData(argv[2]);
    cout << "Finished reading file: " << argv[2] << endl;

    auto particle_id = Cabana::slice<HACCabana::Particles::Fields::ParticleID>( P.aosoa_host, "particle_id" );
    auto sort_data = Cabana::sortByKey( particle_id );
    Cabana::permute( sort_data, P.aosoa_host );

    auto particle_id_ans = Cabana::slice<HACCabana::Particles::Fields::ParticleID>( P_ans.aosoa_host, "particle_id_ans" );
    auto sort_data_ans = Cabana::sortByKey( particle_id_ans );
    Cabana::permute( sort_data_ans, P_ans.aosoa_host );

    auto position = Cabana::slice<HACCabana::Particles::Fields::Position>( P.aosoa_host, "position" );
    auto position_ans = Cabana::slice<HACCabana::Particles::Fields::Position>( P_ans.aosoa_host, "position_ans" );

    cout << "Checking " << P.num_p << " particles against " << P_ans.num_p << " answer particles." << endl;
    assert(P.num_p == P_ans.num_p);

    // don't check near the boundary
    const float dx_boundary = 4.0;

    cout << "\tChceking particles within [" << MIN_POS+dx_boundary << "," << MAX_POS-dx_boundary << ")" << endl;

    int count = 0;
    int err_n = 0;
    for (int i=0; i<P_ans.num_p; ++i)
    {
      assert(particle_id(i) == particle_id_ans(i));
      bool is_inside = false;
      if (position(i,0) >= MIN_POS+dx_boundary &&\
          position(i,1) >= MIN_POS+dx_boundary &&\
          position(i,2) >= MIN_POS+dx_boundary &&\
          position(i,0) <  MAX_POS-dx_boundary &&\
          position(i,1) <  MAX_POS-dx_boundary &&\
          position(i,2) <  MAX_POS-dx_boundary)
      {
        is_inside = true;
        ++count;
      }
      if (is_inside && (!floatCompare(position(i,0),position_ans(i,0)) ||\
                        !floatCompare(position(i,1),position_ans(i,1)) ||\
                        !floatCompare(position(i,2),position_ans(i,2))))
      {
        ++err_n;
      }
    }
    cout << "\t" << err_n << " particles (out of " << count << ") have position relative error greater than " << MAX_ERR << endl;
  }

  return 0;
}
