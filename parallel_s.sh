#!/bin/sh
#$ -cwd
#$ -l f_node=12
#$ -l h_rt=05:00:00
#$ -p -5

## The initialization part and mpi environment vary with different cluster systems
## Setting here applies to TSUBAME3.0 at the Tokyo Institute of Technology

. /etc/profile.d/modules.sh

module load cuda intel intel-mpi
#export PSM2_MEMORY=large
#export PSM2_MQ_RECVREQS_MAX=10000000

mpiexec.hydra -genv PSM2_MEMORY large -ppn 28 -n 336 ./svd_sonly
# mpiexec.hydra -genv PSM2_MEMORY large -ppn 28 -n 336 ./mdinvs_sonly6000
# mpiexec.hydra -genv PSM2_MEMORY large -ppn 28 -n 336 ./rconst_sonly_low6000
mpiexec.hydra -genv PSM2_MEMORY large -ppn 28 -n 336 ./resopr_sonly_low_full6000
