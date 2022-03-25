#!/bin/bash

#SBATCH -p sched_mit_darwin 
#SBATCH -J r3dSAMPLE 
#SBATCH -N 2                    # number of nodes to use (each has 16 cores)
#SBATCH --ntasks-per-node=16   #1/18/17: from engaging email, since some nodes have > 16 cpus 
##SBATCH -n 32                   # total number of cores to run on
#SBATCH --mem-per-cpu 4000      # MB memory per core needed (16-core nodes have 64GB)
#SBATCH --time 12:00:00         # run no longer than 16 hours (HH:MM:SS)
#SBATCH --error stderr
#SBATCH --output stdout
#SBATCH --mail-type FAIL,END
#SBATCH --mail-user ezakem@mit.edu
##SBATCH --constraint=centos6
#SBATCH --constraint=centos7

module add gcc
##module add mvapich2/gcc
module add mvapich2/gcc/64/2.0b
module add slurm
module add harvard/centos6/hdf5-1.8.11_gcc-4.4.7
module add harvard/centos6/netcdf-4.3.0_gcc-4.4.7

echo `date`

pickupD="/nobackup1/ezakem/gudruns/threed/runNitrifyModel_NatComm_newyields_eqaff_withnitrtemp_lowy_L_hihi3"
runD="$pickupD/SAMPLE"

if [ ! -d "$runD" ]; then
  mkdir -p $runD
fi

cp $pickupD/pickup*720000* $runD

cp SAMPLEdata* $runD
cp data* $runD
cp eedata $runD
cp POLY3.COEFFS $runD
cp mitgcmuv $runD
cp Makefile $runD

cd $runD

mv SAMPLEdata data
mv SAMPLEdata.diagnostics data.diagnostics

mpirun -n 32 ./mitgcmuv

echo `date`
