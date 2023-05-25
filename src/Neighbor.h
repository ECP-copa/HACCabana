#ifndef NEIGHBORS_H
#define NEIGHBORS_H

#include "HACCabana_Config.h"

namespace HACCabana
{
  struct LinkedCellTag
  {
  };
  struct ArborXTag
  {
  };

  template <typename AoSoADevice, typename NeighborType>
  struct Neighbors;

  template <typename AoSoADevice>
  struct Neighbors<AoSoADevice, LinkedCellTag>
  {
    using device_type = typename AoSoADevice::device_type;
    using neighbor_type = Cabana::LinkedCellList<device_type>;
    using tag_type = LinkedCellTag;
    // Store to avoid reallocation during rebuild.
    neighbor_type cell_list;

    Neighbors(AoSoADevice aosoa_device,
              const float, const float cm_size, const float min_pos,
              const float max_pos)
    {
      auto position = Cabana::slice<HACCabana::Particles::Fields::Position>(
          aosoa_device, "position");

      // create the cell list on the GPU
      // NOTE: fuzz particles (outside of overload) are not included
      float dx = cm_size;
      float x_min = min_pos;
      float x_max = max_pos;

      float grid_delta[3] = {dx, dx, dx};
      float grid_min[3] = {x_min, x_min, x_min};
      float grid_max[3] = {x_max, x_max, x_max};

      cell_list = neighbor_type(position, grid_delta, grid_min, grid_max);
    }

    template <typename ParticlesType>
    void update(ParticlesType aosoa_device)
    {
      cell_list.build(aosoa_device);
      Cabana::permute(cell_list, aosoa_device);
    }
  };

#ifdef HACCabana_ENABLE_ARBORX
  template <typename AoSoADevice>
  struct Neighbors<AoSoADevice, ArborXTag>
  {
    using device_type = typename AoSoADevice::device_type;
    using neighbor_type =
        Cabana::Experimental::CrsGraph<typename device_type::memory_space,
                                       Cabana::FullNeighborTag>;
    using tag_type = ArborXTag;
    neighbor_type neighbor_list;

    Neighbors(AoSoADevice aosoa_device, const float rmax2, const float,
              const float, const float)
    {
      auto position = Cabana::slice<HACCabana::Particles::Fields::Position>(
          aosoa_device, "position");
      // Interface in Cabana to ArborX neighbor search.
      neighbor_list = Cabana::Experimental::makeNeighborList<device_type>(
          Cabana::FullNeighborTag{}, position, 0, position.size(), rmax2);
    }
  };
#endif

  template <typename NeighborType, typename AoSoADevice>
  auto createNeighbors(AoSoADevice &aosoa_device, const float rmax2, const float cm_size, const float min_pos,
                       const float max_pos)
  {
    return Neighbors<AoSoADevice, NeighborType>(aosoa_device, rmax2, cm_size, min_pos, max_pos);
  }
}
#endif