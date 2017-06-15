#!/bin/sh 
########## Begin MOAB/Slurm header ##########
#
# Give job a reasonable name
#MOAB -N STRUCTURE
#
# Request number of nodes and CPU cores per node for job
#MOAB -l nodes=1:ppn=1
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
#MOAB -v K=

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
module load bio/structure/2.3.4

# Echo input variables
echo "K=${K}"

# Set-up program
startprog="Rscript --no-save --no-restore --slave\
           ./analysis/maizego_structure.R\
           ${K}"

# Start program
echo $startprog
echo $date
exec $startprog
echo $date