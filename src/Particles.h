#ifndef PARTICLES_H
#define PARTICLES_H

#include <string>
#include <Cabana_Core.hpp>
#include <Cabana_AoSoA.hpp>
#include "Definitions.h"

namespace HACCabana
{

  class Particles
  {
    public:
      enum Fields
      {
        ParticleID = 0,
        Position = 1,
        Velocity = 2
      };

      using data_types = Cabana::MemberTypes<int64_t, float[3], float[3]>;
      using aosoa_host_type = Cabana::AoSoA<data_types, Kokkos::Device<Kokkos::Serial, Kokkos::HostSpace>>;

      size_t num_p = 0;
      size_t begin = 0;
      size_t end = 0;
      aosoa_host_type aosoa_host;

      Particles();
      ~Particles();
      void readRawData(std::string file_name);
  };

}
#endif
