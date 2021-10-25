# HACCabana
is a miniapp for an N-body cosmology code, based on HACC and the CoPA [Cabana Particle Toolkit](https://github.com/ECP-copa/Cabana). The driver bootstraps a direct particle-particle short-range solver designed for GPUs. The solver follows the method used by HACC's P3M implementation for computing the short-range component of the gravitational force.

## Build
The only required input to build is the Cabana install location. The default Kokkos device is used to run HACCabana.
```
mkdir build
cd build
cmake -DCMAKE_PREFIX_PATH=$CABANA_INSTALL ..
make
```

## Run
Using input and verification files (generated from HACC).

``./build/driver_short-range -i ../PRE_rank0_particle_data.bin -v ../POST_rank0_particle_data.bin -t 0 -c ./256.indat.params``

or use synthetically generated data

``./build/driver_short-range -s -t 0 -c ./256.indat.params``
