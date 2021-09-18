# example run using real data taken at step 7
#OMP_PROC_BIND=false ./build/driver_short-range -i ../PRE_rank0_step7_particle_data.bin -v ../POST_rank0_step7_particle_data.bin -t 7
OMP_PROC_BIND=false ./build/driver_short-range -i ../PRE_rank0_particle_data.bin -v ../POST_rank0_particle_data.bin -t 0
#OMP_PROC_BIND=false ./build/driver_short-range -s -t 0
