#include "ParticleActions.h"

#include <sys/time.h>
double mytime() {
  timeval tv;
  gettimeofday(&tv, NULL);
  double time = 1.0*tv.tv_sec;
  time += tv.tv_usec*1.e-6;
  return time;
}

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

namespace HACCabana
{

ParticleActions::ParticleActions() {};

ParticleActions::ParticleActions(Particles *P_) : P(P_)
{
  ;
};

ParticleActions::~ParticleActions()
{
  ;
};

void ParticleActions::setParticles(Particles *P_)
{
  P = P_;
}

// Stream
void ParticleActions::updatePos(\
    Cabana::AoSoA<HACCabana::Particles::data_types, device_mem, VECTOR_LENGTH> aosoa_device,\
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

// Kick
void ParticleActions::updateVel(\
    Cabana::AoSoA<HACCabana::Particles::data_types, device_mem, VECTOR_LENGTH> aosoa_device,\
    Cabana::LinkedCellList<device_mem> cell_list,\
    const float c, const float rmax2, const float rsm2)
{
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
            const float dist2 = dx * dx + dy * dy + dz * dz;
            if (dist2 < rmax2) 
            {
              const float dist2Err = dist2 + rsm2;
              const float tmp =  1.0f/Kokkos::sqrt(dist2Err*dist2Err*dist2Err) - FGridEvalPoly(dist2);
              force[0] += dx * tmp;
              force[1] += dy * tmp;
              force[2] += dz * tmp;
            }
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

void ParticleActions::subCycle(TimeStepper &ts, const int nsub, const float gpscal, const float rmax2, const float rsm2, 
    const float cm_size, const float min_pos, const float max_pos)
{
  // copy particles to GPU
  Cabana::AoSoA<HACCabana::Particles::data_types, device_mem, VECTOR_LENGTH> aosoa_device("aosoa_device", P->num_p);
  Cabana::deep_copy(aosoa_device, P->aosoa_host);

  // create the cell list on the GPU
  // NOTE: fuzz particles (outside of overload) are not included
  float dx = cm_size;
  float x_min = min_pos;
  float x_max = max_pos;

  float grid_delta[3] = {dx, dx, dx};
  float grid_min[3] = {x_min, x_min, x_min};
  float grid_max[3] = {x_max, x_max, x_max};

  auto position = Cabana::slice<HACCabana::Particles::Fields::Position>(aosoa_device, "position");
  Cabana::LinkedCellList<device_mem> cell_list(position, P->begin, P->end, grid_delta, grid_min, grid_max);
  Cabana::permute(cell_list, aosoa_device);
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
    this->updatePos(aosoa_device, prefactor*tau*0.5);

    // kick
    double tmp = mytime();
    this->updateVel(aosoa_device, cell_list, c, rmax2, rsm2);
    kick_time += mytime() - tmp;

    //half stream
    this->updatePos(aosoa_device, prefactor*tau*0.5);
  }

  std::cout << "kick time " << kick_time << std::endl;

  // copy GPU particles back to host
  Cabana::deep_copy(P->aosoa_host, aosoa_device);
}

}
