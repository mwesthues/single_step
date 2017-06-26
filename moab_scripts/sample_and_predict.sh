#!/bin/sh 
########## Begin MOAB/Slurm header ##########
#
# Give job a reasonable name
#MOAB -N sample_predict
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
#MOAB -d /pfs/work2/workspace/scratch/ho_westhues-single_step-0
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
# Use either 'Hybrid' (UHOH data) or 'Inbred' (Yan-lab data)
#MOAB -v DATA_TYPE=
#
#MOAB -v TRAIT=
#
#MOAB -v ITER=30000
#
#MOAB -v MODEL=BRR
#
#MOAB -v VCOV=RadenII
#
#MOAB -v PI=0.5
#
#MOAB -v PRIOR_PI_COUNT=10
#
# Predictor: 'ped', 'snp' or 'mrna'
#MOAB -v PRED1=
#
#MOAB -v PRED2=
#
# Specify the core set set size, as a fraction of genotypes covered by mrna
# data
#MOAB -v CORE_SET=
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
echo "Data_Type=${DATA_TYPE}\
      Trait=${TRAIT}\
      Iter=${ITER}\
      Model=${MODEL}\
      VCOV=${VCOV}\
      Pi=${PI}\
      PriorPiCount=${PRIOR_PI_COUNT}\
      Pred1=${PRED1}\
      Pred2=${PRED2}\
      Core_Set=${CORE_SET}\
      Runs=${RUNS}"

# Set-up program
startprog="Rscript --no-save --no-restore --slave\
           ./analysis/sample_and_predict.R\
           ${DATA_TYPE} ${TRAIT} ${ITER} ${MODEL} ${VCOV} ${PI}\
           ${PRIOR_PI_COUNT} ${PRED1} ${PRED2} ${CORE_SET}\
           ${RUNS}"

# Start program
echo $startprog
echo $date
exec $startprog
echo $date
