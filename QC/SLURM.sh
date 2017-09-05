#!/bin/bash
#SBATCH --time=99:00:00          # walltime
#SBATCH --nodes=1                # number of cluster nodes
#SBATCH -o z-%j.out-%N	 # name of the stdout
#SBATCH -e z-%j.err-%N	 # name of the stderr, using job and first node values
#SBATCH --ntasks=1               # number of MPI tasks
# additional information for allocated clusters
#SBATCH --partition=usubio-kp
#SBATCH --account=usubio-kp
#SBATCH --mail-type=END
#SBATCH --mail-user=jenessa.lemon@gmail.com
#
# set data and working directories
WORKDIR=/uufs/chpc.utah.edu/common/home/wolf-group2/jlemon/ #set working
SCRDIR=/scratch/local/$USER         #set scratch
mkdir -p SCRDIR                        #make scratch
cp -r $WORKDIR/* $SCRDIR/.             #transfer everyting in working to scratch
cd $SCRDIR                             #go to scratch
#
# run the program
/uufs/chpc.utah.edu/common/home/u6009817/miniconda2/bin/ipyrad -p params-51_individuals.txt -s 234567
#
# transfer output and remove scratch
rsync -avzh ./ $WORKDIR
cd $WORKDIR
rm -rf $SCRATCH

