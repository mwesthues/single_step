#!/bin/sh 
########## Begin MOAB/Slurm header ##########
#
# Give job a reasonable name
#MOAB -N ssBLUP
#
# Request number of nodes and CPU cores per node for job
#MOAB -l nodes=1:ppn=16
#
# Request memory per process
#MOAB -l pmem=3800mb
#
# Estimated wallclock time for job
#MOAB -l walltime=00:22:00:00
#
# Job submission directory
#MOAB -d /pfs/work2/workspace/scratch/ho_westhues-single_step-0/gamazon
#
# Standard output naming
#MOAB -o $(JOBNAME)_$(JOBID)
#
# Queue: Use a single node 
#MOAB -q singlenode
#
# Send email when job begins (b), aborts (a) and ends (e)
#MOAB -m a
#
# Write standard output (o) and errors (e) to the same file
#MOAB -j oe
#
#MOAB -v TRAIT=
#
#MOAB -v ITER=30000
#
#MOAB -v PRED1=ped100
#
#MOAB -v PRED2=
#
#MOAB -v PRED3=
#
#MOAB -v REPLICATION=
#
#MOAB -v RUNS=
#
##### **********************************************************************
########### End MOAB header ##########

echo "Working Directory:                    $PWD"
echo "Running on host                       $HOSTNAME"
echo "Job id:                               $MOAB_JOBID"
echo "Job name:                             $MOAB_JOBNAME"
echo "Number of nodes allocated to job:     $MOAB_NODECOUNT"
echo "Number of cores allocated to job:     $MOAB_PROCCOUNT"

# Setup R Environment
module load math/R/3.3.1

# Echo input variables
echo "Trait=${TRAIT}\
      Iter=${ITER}\
      Pred1=${PRED1}\
      Pred2=${PRED2}\
      Pred3=${PRED3}\
      Replication=${REPLICATION}\
      Runs=${RUNS}"

# Set-up program
startprog="Rscript --no-save --no-restore --slave\
           ./analysis/repeated_snp77_prediction.R\
           ${TRAIT} ${ITER} ${PRED1} ${PRED2} ${PRED3}\
           ${REPLICATION} ${RUNS}"

# Start program
echo $startprog
echo $date
exec $startprog
echo $date
