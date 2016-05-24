#!/bin/sh 
########## Begin MOAB/Slurm header ##########
#
# Give job a reasonable name
#MOAB -N mrna_imputation
#
# Request number of nodes and CPU cores per node for job
#MOAB -l nodes=1:ppn=16
#
# Request memory per process
#MOAB -l pmem=3500mb
#
# Estimated wallclock time for job
#MOAB -l walltime=00:07:00:00
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
##### USER-SPECIFIC INPUT *************************************************
### Add some parameters that will be used inside the R script.
#
#MOAB -v ITER=50000
#
#MOAB -v MODEL=BRR_Kernel
#
#MOAB -v SNP_FILTER=FALSE
#
#MOAB -v VCOV=RadenII
#
#MOAB -v POOL=Dent

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
echo "Iter=${ITER} \
      Model=${MODEL}\
      SNP_FILTER=${SNP_FILTER}\
      VCOV=${VCOV}\
      POOL=${POOL}"
# Set-up program
startprog="Rscript --no-save --no-restore --slave\
           ./analysis/mrna_imputation.R\
           ${ITER} ${MODEL} ${SNP_FILTER} ${VCOV} ${POOL}"

# Start program
echo $startprog
echo $(date)
exec $startprog
echo $(date)
