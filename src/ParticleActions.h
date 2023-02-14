#ifndef PARTICLE_ACTIONS_H
#define PARTICLE_ACTIONS_H

#include <string>
#include <Cabana_Core.hpp>
#include <Cabana_AoSoA.hpp>
#include <Cabana_LinkedCellList.hpp>
#include "Definitions.h"
#include "TimeStepper.h"
#include "Particles.h"
#include "Neighbor.h"
#include "HACCabana_Config.h"

// This might need to be changed for AMD and Intel GPUs. Nvidia warp size is 32. 
#define VECTOR_LENGTH 64

#include <sys/time.h>
inline double mytime() {
  timeval tv;
  gettimeofday(&tv, NULL);
  double time = 1.0*tv.tv_sec;
  time += tv.tv_usec*1.e-6;
  return time;
}

namespace HACCabana 
{

// Polynomial long range force calculation
KOKKOS_INLINE_FUNCTION
float FGridEvalPoly(float r2)
{
#if POLY_ORDER == 6
  return (0.271431f + r2*(-0.0783394f + r2*(0.0133122f + r2*(-0.00159485f + r2*(0.000132336f + r2*(-0.00000663394f + r2*0.000000147305f))))));
#elif POLY_ORDER == 5
  return (0.269327f + r2*(-0.0750978f + r2*(0.0114808f + r2*(-0.00109313f + r2*(0.0000605491f + r2*-0.00000147177f)))));
#elif POLY_ORDER == 4
  return (0.263729f + r2*(-0.0686285f + r2*(0.00882248f + r2*(-0.000592487f + r2*0.0000164622f))));
#else
  return 0.0f;
#endif
}

  template <typename Derived>
  class ParticleActionsBase
  {
  protected:
    Particles *P;

  public:
    using device_exec = Kokkos::DefaultExecutionSpace::execution_space;
    using device_mem = Kokkos::DefaultExecutionSpace::memory_space;
    using device_type = Kokkos::Device<device_exec, device_mem>;
    //using device_scratch = Kokkos::ScratchMemorySpace<device_exec>;

    using derived_type = Derived;

    ParticleActionsBase() {};
    ParticleActionsBase(Particles *P_) : P(P_)
    {
    };
    ~ParticleActionsBase() {};
    void setParticles(Particles *P_)
    {
      P = P_;
    }

    void updatePos(\
        Cabana::AoSoA<HACCabana::Particles::data_types, device_type, VECTOR_LENGTH> aosoa_device,\
        float prefactor)
    {
      auto position = Cabana::slice<HACCabana::Particles::Fields::Position>(aosoa_device, "position");
      auto velocity = Cabana::slice<HACCabana::Particles::Fields::Velocity>(aosoa_device, "velocity");

      Kokkos::parallel_for("stream", Kokkos::RangePolicy<device_exec>(0, P->num_p),
      KOKKOS_LAMBDA(const int i) {
        position(i,0) = position(i,0) + prefactor * velocity(i,0);
        position(i,1) = position(i,1) + prefactor * velocity(i,1);
        position(i,2) = position(i,2) + prefactor * velocity(i,2);
      });
      Kokkos::fence();
    }

    KOKKOS_INLINE_FUNCTION void calcVel(float force[3], const float dx,
                                        const float dy, const float dz,
                                        const float rmax2, const float rsm2) {
        const float dist2 = dx * dx + dy * dy + dz * dz;
        if (dist2 < rmax2) {
            const float dist2Err = dist2 + rsm2;
            const float tmp =
                1.0f / Kokkos::Experimental::sqrt(dist2Err * dist2Err * dist2Err) -
                FGridEvalPoly(dist2);
            force[0] += dx * tmp;
            force[1] += dy * tmp;
            force[2] += dz * tmp;
        }
    }

  void subCycle(TimeStepper &ts, const int nsub, const float gpscal, const float rmax2, const float rsm2,
  const float cm_size, const float min_pos, const float max_pos)
  {
    // copy particles to GPU
    Cabana::AoSoA<HACCabana::Particles::data_types, device_type, VECTOR_LENGTH> aosoa_device("aosoa_device", P->num_p);
    Cabana::deep_copy(aosoa_device, P->aosoa_host);
    auto neighbors = createNeighbors<typename derived_type::neighbor_tag>(aosoa_device, rmax2, cm_size, min_pos, max_pos);
    Kokkos::fence();

    // ------------------------------------------------------------------------------------
    const double stepFraction = 1.0/nsub;

    // mass convertion from IC np units to ng grid units. Required since mass = 1 in Gravity only HACC. 
    const float divscal = gpscal*gpscal*gpscal;
    const float pi = 4.0*atanf(1.0);
    // c = G*dt/dy*1/a in code units. 
    // a = Gm/x^2. 
    // G = (3/2)Omega_m *(1/4pi) 
    // in code units, dt/dy for time transformation, 1/a since dp/dy = -gradphi/a.
    const float c = divscal/4.0/pi*ts.fscal()*ts.tau()*stepFraction;

    const float pf = powf(ts.pp(), (1.0 + 1.0 / ts.alpha()));
    const float prefactor = 1.0 / (ts.alpha() * ts.adot() * pf);
    float tau = ts.tau()*stepFraction;

    double kick_time = 0.0f;
    // SKS subcycles
    for(int step=0; step < nsub; ++step) 
    {
    std::cout << "Doing substep " << step << std::endl;

    //half stream
    updatePos(aosoa_device, prefactor*tau*0.5);

    // kick
    double tmp = mytime();
    static_cast<derived_type*>( this )->updateVel(aosoa_device, neighbors, c, rmax2, rsm2);
    kick_time += mytime() - tmp;

    //half stream
    updatePos(aosoa_device, prefactor*tau*0.5);
    }

    std::cout << "kick time " << kick_time << std::endl;

    // copy GPU particles back to host
    Cabana::deep_copy(P->aosoa_host, aosoa_device);
  }

  };

  
template <typename NeighborTag>
class ParticleActions;
  
template <>
class ParticleActions<LinkedCellTag> : public ParticleActionsBase<ParticleActions<LinkedCellTag>>
{
  public:
  using base_type = ParticleActionsBase<ParticleActions<LinkedCellTag>>;
  // Use base class constructors.
  using base_type::base_type;
  using neighbor_tag = LinkedCellTag;

  void updateVel(\
    Cabana::AoSoA<HACCabana::Particles::data_types, device_type, VECTOR_LENGTH> aosoa_device,\
    Neighbors<Cabana::AoSoA<HACCabana::Particles::data_types, device_type, VECTOR_LENGTH>, LinkedCellTag> neighbors,\
    const float c, const float rmax2, const float rsm2)
  {
    auto& cell_list = neighbors.cell_list;
    auto position = Cabana::slice<HACCabana::Particles::Fields::Position>(aosoa_device, "position");
    auto velocity = Cabana::slice<HACCabana::Particles::Fields::Velocity>(aosoa_device, "velocity");
    auto bin_index = Cabana::slice<HACCabana::Particles::Fields::BinIndex>(aosoa_device, "bin_index");

    Kokkos::parallel_for("copy_bin_index", Kokkos::RangePolicy<device_exec>(0, cell_list.totalBins()),
    KOKKOS_LAMBDA(const int i)
    {
      int bin_ijk[3];
      cell_list.ijkBinIndex(i, bin_ijk[0], bin_ijk[1], bin_ijk[2]);
      for (size_t ii = cell_list.binOffset(bin_ijk[0], bin_ijk[1], bin_ijk[2]); 
          ii < cell_list.binOffset(bin_ijk[0], bin_ijk[1], bin_ijk[2]) +
          cell_list.binSize(bin_ijk[0], bin_ijk[1], bin_ijk[2]); 
          ++ii)
        bin_index(ii) = i;
    });
    Kokkos::fence();

    auto vector_kick = KOKKOS_LAMBDA(const int s, const int a)
    {
      int bin_ijk[3];
      cell_list.ijkBinIndex(bin_index.access(s,a), bin_ijk[0], bin_ijk[1], bin_ijk[2]);

      float force[3] = {0.0, 0.0, 0.0};
      for (int ii=-1; ii<2; ++ii) 
      {
        if (bin_ijk[0] + ii < 0 || bin_ijk[0] + ii >= cell_list.numBin(0))
          continue;
        for (int jj=-1; jj<2; ++jj) 
        {
          if (bin_ijk[1] + jj < 0 || bin_ijk[1] + jj >= cell_list.numBin(1))
            continue;
          for (int kk=-1; kk<2; ++kk) 
          {
            if (bin_ijk[2] + kk < 0 || bin_ijk[2] + kk >= cell_list.numBin(2))
              continue;

            const size_t binOffset = cell_list.binOffset(bin_ijk[0] + ii, bin_ijk[1] + jj, bin_ijk[2] + kk);
            const int binSize = cell_list.binSize(bin_ijk[0] + ii, bin_ijk[1] + jj, bin_ijk[2] + kk);

            for (int j = binOffset; j < binOffset+binSize; ++j) 
            {
              const float dx = position(j,0)-position.access(s,a,0);
              const float dy = position(j,1)-position.access(s,a,1);
              const float dz = position(j,2)-position.access(s,a,2);
              calcVel(force, dx, dy, dz, rmax2, rsm2);
            }

          }
        }
      }
      velocity.access(s,a,0) += force[0] * c;
      velocity.access(s,a,1) += force[1] * c;
      velocity.access(s,a,2) += force[2] * c;
    };

    Cabana::SimdPolicy<VECTOR_LENGTH, device_exec> simd_policy(P->begin, P->end);
    Cabana::simd_parallel_for( simd_policy, vector_kick, "kick" ); 

    Kokkos::fence();
  }

  protected:
  using base_type::P;
};

#ifdef HACCabana_ENABLE_ARBORX
template <>
class ParticleActions<ArborXTag> : public ParticleActionsBase<ParticleActions<LinkedCellTag>>
{
  public:
  using base_type = ParticleActionsBase<ParticleActions<LinkedCellTag>>;
  // Use base class constructors.
  using base_type::base_type;
  using neighbor_tag = ArborXTag;

  void updateVel(\
  Cabana::AoSoA<HACCabana::Particles::data_types, device_type, VECTOR_LENGTH> aosoa_device,\
  Neighbors<Cabana::AoSoA<HACCabana::Particles::data_types, device_type, VECTOR_LENGTH>, ArborXTag> neighbors,\
  const float c, const float rmax2, const float rsm2)
  {
    auto& neighbor_list = neighbors.neighbor_list;
    auto position = Cabana::slice<HACCabana::Particles::Fields::Position>(aosoa_device, "position");
    auto velocity = Cabana::slice<HACCabana::Particles::Fields::Velocity>(aosoa_device, "velocity");

    auto kick = KOKKOS_LAMBDA(const int i, const int j)
    {
      float force[3] = {0.0, 0.0, 0.0};
      const float dx = position(j,0)-position(i,0);
      const float dy = position(j,1)-position(i,1);
      const float dz = position(j,2)-position(i,2);
      calcVel(force, dx, dy, dz, rmax2, rsm2);

      velocity(i,0) += force[0] * c;
      velocity(i,1) += force[1] * c;
      velocity(i,2) += force[2] * c;
    };
    Kokkos::RangePolicy<device_exec> policy( P->begin, P->end );
    Cabana::neighbor_parallel_for( policy, kick, neighbor_list,
                                  Cabana::FirstNeighborsTag(),
                                  Cabana::SerialOpTag(), "kick" );
    Kokkos::fence();
  }

  protected:
  using base_type::P;
};
#endif

}
#endif
