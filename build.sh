module use /soft/modulefiles
module load cmake gcc/6.5.0 cuda openmpi

mkdir build
cd build

export KOKKOS_INSTALL=/gpfs/jlse-fs0/users/erangel/HACCabana/kokkos/build/install
export CABANA_INSTALL=/gpfs/jlse-fs0/users/erangel/HACCabana/Cabana/build/install

cmake -DCMAKE_PREFIX_PATH=$CABANA_INSTALL -DCMAKE_CXX_COMPILER=$KOKKOS_INSTALL/bin/nvcc_wrapper -DENABLE_GPU=ON ..
make
cd ..
