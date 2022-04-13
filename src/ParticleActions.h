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
    using device_exec = Kokkos::DefaultExecutionSpace::execution_space;
    using device_mem = Kokkos::DefaultExecutionSpace::memory_space;
    using device_type = Kokkos::Device<device_exec, device_mem>;

    ParticleActions();
    ParticleActions(Particles *P_);
    ~ParticleActions();
    void setParticles(Particles *P_);
    void subCycle(TimeStepper &ts, const int nsub, const float gpscal, const float rmax2, const float rsm2,\
        const float cm_size, const float min_pos, const float max_pos);
    void updatePos(Cabana::AoSoA<HACCabana::Particles::data_types, device_type> &aosoa_device,\
        float prefactor);
    void updateVel(Cabana::AoSoA<HACCabana::Particles::data_types, device_type> &aosoa_device,\
        Cabana::LinkedCellList<device_type> cell_list,\
        const float c, const float rmax2, const float rsm2);
  };
}
#endif
