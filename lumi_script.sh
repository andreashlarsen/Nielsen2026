#!/bin/bash
#SBATCH --job-name=d22_rep2
#SBATCH --account=project_465001110
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=small
#SBATCH --array=
module purge
module load LUMI/23.09
module load partition/L
module load GROMACS/2021.7-cpeGNU-23.09-CPU
gmx mdrun -deffnm umbrella${SLURM_ARRAY_TASK_ID} -v -nt ${SLURM_CPUS_PER_TASK}
#!/bin/bash
#SBATCH --job-name=d22_rep2
#SBATCH --account=project_465001110
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=small
#SBATCH --array=
module purge
module load LUMI/23.09
module load partition/L
module load GROMACS/2021.7-cpeGNU-23.09-CPU
gmx mdrun -deffnm umbrella${SLURM_ARRAY_TASK_ID} -v -nt ${SLURM_CPUS_PER_TASK}
