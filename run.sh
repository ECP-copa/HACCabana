# example run using real data taken at step 7
OMP_PROC_BIND=false ./build/driver_short-range -i ../1-rank-data/PRE_rank0_step7_particle_data.bin -v ../1-rank-data/POST_rank0_step7_particle_data.bin -t 7 -c 256.indat.params
#OMP_PROC_BIND=false ./build/driver_short-range -s -t 0 -c 256.indat.params
