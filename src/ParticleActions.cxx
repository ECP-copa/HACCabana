#include "ParticleActions.h"

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
    Cabana::AoSoA<HACCabana::Particles::data_types, device_type> &aosoa_device,\
    float prefactor)
{
  auto position = Cabana::slice<HACCabana::Particles::Fields::Position>(aosoa_device, "position");
  auto velocity = Cabana::slice<HACCabana::Particles::Fields::Velocity>(aosoa_device, "velocity");

  Kokkos::parallel_for("initialize", Kokkos::RangePolicy<device_exec>(0, P->num_p),
  KOKKOS_LAMBDA(const int i) {
    position(i,0) = position(i,0) + prefactor * velocity(i,0);
    position(i,1) = position(i,1) + prefactor * velocity(i,1);
    position(i,2) = position(i,2) + prefactor * velocity(i,2);
  });
  Kokkos::fence();
}

// Kick
void ParticleActions::updateVel(\
    Cabana::AoSoA<HACCabana::Particles::data_types, device_type> &aosoa_device,\
    Cabana::LinkedCellList<device_type> cell_list,\
    const float c, const float rmax2, const float rsm2,
    const float cell_size, const float min_pos)
{
  auto position = Cabana::slice<HACCabana::Particles::Fields::Position>(aosoa_device, "position");
  auto velocity = Cabana::slice<HACCabana::Particles::Fields::Velocity>(aosoa_device, "velocity");

  Kokkos::parallel_for("initialize", Kokkos::RangePolicy<device_exec>(P->begin, P->end),
  KOKKOS_LAMBDA(const int i)
  {
    // the cell that contains this particle
    int my_bin_pos[3];

    // AoSoA needs to be permuted for this to work 
    bool found = false;
    for (int ii=0; ii<cell_list.totalBins(); ++ii) 
    {
      cell_list.ijkBinIndex(ii, my_bin_pos[0], my_bin_pos[1], my_bin_pos[2]);
      size_t binOffset = cell_list.binOffset(my_bin_pos[0], my_bin_pos[1], my_bin_pos[2]);
      size_t binSize = cell_list.binSize(my_bin_pos[0], my_bin_pos[1], my_bin_pos[2]);
      if (i >= binOffset && i < binOffset + binSize) 
      {
        found = true;
        break;
      }
    }

#ifdef DEBUG_ME
    // verifying that the particle is indeed within the cell
    // NOTE: padding the cell boundary since the particle may have moved outside 
    // during a previous subcycle
    float pad_val = 0.5;
    assert( found );
    assert( position(i,0) >= (float)my_bin_pos[0] * cell_size + min_pos - pad_val );
    assert( position(i,1) >= (float)my_bin_pos[1] * cell_size + min_pos - pad_val );
    assert( position(i,2) >= (float)my_bin_pos[2] * cell_size + min_pos - pad_val );
    assert( position(i,0) <= my_bin_pos[0] * cell_size + cell_size + min_pos + pad_val );
    assert( position(i,1) <= my_bin_pos[1] * cell_size + cell_size + min_pos + pad_val );
    assert( position(i,2) <= my_bin_pos[2] * cell_size + cell_size + min_pos + pad_val );
#endif
    
    float force[3] = {0.0, 0.0, 0.0};
    for (int ii=-1; ii<2; ++ii) 
    {
      if (my_bin_pos[0] + ii < 0 || my_bin_pos[0] + ii >= cell_list.numBin(0))
        continue;
      for (int jj=-1; jj<2; ++jj) 
      {
        if (my_bin_pos[1] + jj < 0 || my_bin_pos[1] + jj >= cell_list.numBin(1))
          continue;
        for (int kk=-1; kk<2; ++kk) 
        {
          if (my_bin_pos[2] + kk < 0 || my_bin_pos[2] + kk >= cell_list.numBin(2))
            continue;
          size_t binOffset = cell_list.binOffset(my_bin_pos[0] + ii, my_bin_pos[1] + jj, my_bin_pos[2] + kk);
          size_t binSize = cell_list.binSize(my_bin_pos[0] + ii, my_bin_pos[1] + jj, my_bin_pos[2] + kk);
          for (int j = binOffset; j < binOffset+binSize; ++j) 
          {
            if (i != j) // ??? 
            {
              const float dx = position(j,0)-position(i,0);
              const float dy = position(j,1)-position(i,1);
              const float dz = position(j,2)-position(i,2);
              const float dist2 = dx * dx + dy * dy + dz * dz;
              if (dist2 < rmax2) 
              {
                const float dist2Err = dist2 + rsm2;
                const float tmp =  rsqrtf(dist2Err*dist2Err*dist2Err) - FGridEvalPoly(dist2);
                force[0] += dx * tmp;
                force[1] += dy * tmp;
                force[2] += dz * tmp;
              }
            }
          }
        }
      }
    }
    velocity(i,0) += force[0] * c;
    velocity(i,1) += force[1] * c;
    velocity(i,2) += force[2] * c;
  });
  Kokkos::fence();
}

void ParticleActions::subCycle(TimeStepper &ts, const int nsub, const float gpscal, const float rmax2, const float rsm2, 
    const float cm_size, const float min_pos, const float max_pos)
{
  // copy particles to GPU
  Cabana::AoSoA<HACCabana::Particles::data_types, device_type> aosoa_device("aosoa_device", P->num_p);
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
  Cabana::LinkedCellList<device_type> cell_list(position, P->begin, P->end, grid_delta, grid_min, grid_max);
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
    double tmp = MPI_Wtime();
    this->updateVel(aosoa_device, cell_list, c, rmax2, rsm2, cm_size, min_pos);
    kick_time += MPI_Wtime() - tmp;

    //half stream
    this->updatePos(aosoa_device, prefactor*tau*0.5);
  }

  std::cout << "kick time " << kick_time << std::endl;

  // copy GPU particles back to host
  Cabana::deep_copy(P->aosoa_host, aosoa_device);
}

}
