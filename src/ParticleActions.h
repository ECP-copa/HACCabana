#ifndef PARTICLE_ACTIONS_H
#define PARTICLE_ACTIONS_H

#include <string>
#include <Cabana_Core.hpp>
#include <Cabana_AoSoA.hpp>
#include <Cabana_LinkedCellList.hpp>
#include "Definitions.h"
#include "TimeStepper.h"
#include "Particles.h"

namespace HACCabana 
{
  class ParticleActions
  {
  private:
    Particles *P;

  public:
    using device_exec = Kokkos::Cuda;
    using device_type = Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>;

    ParticleActions();
    ParticleActions(Particles *P_);
    ~ParticleActions();
    void setParticles(Particles *P_);
    void subCycle(TimeStepper &ts, const int nsub, const float gpscal, const float rmax2, const float rsm2);
    void updatePos(Cabana::AoSoA<HACCabana::Particles::data_types, device_type> &aosoa_device,\
        float prefactor);
    void updateVel(Cabana::AoSoA<HACCabana::Particles::data_types, device_type> &aosoa_device,\
        Cabana::LinkedCellList<device_type> cell_list,\
        const float c, const float rmax2, const float rsm2);
  };
}
#endif
