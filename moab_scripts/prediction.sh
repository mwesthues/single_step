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
#MOAB -l pmem=3000mb
#
# Estimated wallclock time for job
#MOAB -l walltime=01:10:00:00
#
# Job submission directory
#MOAB -d /pfs/work2/workspace/scratch/ho_westhues-gamazon-0/gamazon
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
#MOAB -v TRAIT=GTM
#
#MOAB -v ITER=50000
#
#MOAB -v MODEL=BRR
#
#MOAB -v VCOV=RadenII
#
#MOAB -v PI=0.5
#
#MOAB -v PRIOR_PI_COUNT=10
#
#MOAB -v PREDICTOR=mrna
#
#MOAB -v DENT_NA_FRACTION=0.05
#
#MOAB -v FLINT_NA_FRACTION=0.05
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
module load math/R/3.2.1

# Echo input variables
echo "Trait=${TRAIT}\
      Iter=${ITER}\
      Model=${MODEL}\
      VCOV=${VCOV}\
      Pi=${PI}\
      PriorPiCount=${PRIOR_PI_COUNT}\
      Predictor=${PREDICTOR}\
      Dent_NA_Fraction=${DENT_NA_FRACTION}\
      Flint_NA_Fraction=${FLINT_NA_FRACTION}"

# Set-up program
startprog="Rscript --no-save --no-restore --slave\
           ./analysis/fernando_ssBLUP.R\
           ${TRAIT} ${ITER} ${MODEL} ${VCOV} ${PI} ${PRIOR_PI_COUNT}\
           ${PREDICTOR} ${DENT_NA_FRACTION} ${FLINT_NA_FRACTION}"

# Start program
echo $startprog
echo $date
exec $startprog
echo $date
