## (de)select part of script
DOWNLOAD=0
PREPARE=0
NPT=0
PULL=0
PUSH=0
SETUP=1
LOCAL=0
ROME=0
ANALYZE=0

## parallelisation
pinoffset=0
CPU=8

name=kor2
rep=2

######################################################################################################
### MAKE WORKING DIR #################################################################################
######################################################################################################

umb_dir=umb_${name}_bE_rep${rep}
mkdir -p ${name}/rep${rep}/${umb_dir}
cd ${name}/rep${rep}

######################################################################################################
### DOWNLOAD #########################################################################################
######################################################################################################

if [ $DOWNLOAD -eq 1 ]
then
  ## download martini sctipts and itp files
  wget http://cgmartini.nl/images/tools/insane/insane.py
  wget http://cgmartini.nl/images/parameters/ITP/martini_v2.2refP.itp
  wget http://cgmartini.nl/images/parameters/ITP/martini_v2.0_ions.itp
  wget http://cgmartini.nl/images/parameters/lipids/Collections/martini_v2.0_lipids_all_201506.itp
  wget http://cgmartini.nl/images/parameters/lipids/PC/POPC/martini_v2.0_POPC_02.itp
  wget http://cgmartini.nl/images/tools/martinize/python3/martinize-2.6/martinize.py
fi

cp ../../kor2/rep1/insane.py .
cp ../../kor2/rep1/martini_v2.2refP.itp .
cp ../../kor2/rep1/martini_v2.0_ions.itp .
cp ../../kor2/rep1/martini_v2.0_lipids_all_201506.itp .
cp ../../kor2/rep1/martini_v2.0_POPC_02.itp .
cp ../../kor2/rep1/martinize.py .

######################################################################################################
### PREPARE ##########################################################################################
######################################################################################################

if [ $PREPARE -eq 1 ]
then

## generate protein with alphafold and place in folder ${name}
protein=../${name}.pdb

## martinize
python3 martinize.py -v -f $protein -o topol.top -dssp dssp -ff martini22p -x cg.pdb -elastic

## convert pdb to gro
echo 1 | gmx trjconv -f cg.pdb -s cg.pdb -o cg.gro -quiet

## solvate with insane
#python2 insane.py -f cg.pdb -o bilayer.gro -dm 2.5 -p topol.top -x 7 -y 7 -z 20 -l POPC:60  -l POPS:20 -l CHOL:20 -sol PW -salt 0.05
python2 insane.py -f cg.pdb -o bilayer.gro -dm 2.5 -p topol.top -x 10 -y 10 -z 25 -l POPC:60  -l POPS:20 -l CHOL:20 -sol PW -salt 0.05

## change topology file
#sed -i -e '0,/POPC /{s/POPC /POPR /}' topol.top
cat << EOF > tmp1
#include "martini_v2.2refP.itp"
#include "martini_v2.0_ions.itp"
#include "martini_v2.0_lipids_all_201506.itp"
#include "Protein_A.itp"

EOF
tail -n +2 topol.top > tmp2
cat tmp1 tmp2 > topol.top
sed -i -e 's/Protein   /Protein_A /g' topol.top

## minimization
mkdir -p mdp_files
cat << EOF > mdp_files/min.mdp
integrator               = steep
nsteps                   = 20000
nstxout                  = 0
nstfout                  = 0
nstlog                   = 100 
cutoff-scheme            = Verlet
nstlist                  = 20
pbc                      = xyz
verlet-buffer-tolerance  = 0.005
coulombtype              = reaction-field 
rcoulomb                 = 1.1
epsilon_r                = 2.5
epsilon_rf               = 0
vdw_type                 = cutoff  
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1
EOF
gmx grompp -f mdp_files/min.mdp -c bilayer.gro -p topol.top -o min.tpr -quiet -maxwarn 1 # non-matching names in protein (Martini2 -> Martini2P)
gmx mdrun -deffnm min -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset

# if too many lincs warnings then uncomment below snippet to rerun with, e.g., -nsteps 200, iteratively, i.e. grompp with min.gro instead of bilayer.gro
LINCS=1
if [ $LINCS -eq 1 ]
then
    steps=20 # highest amount of steps before stopping because of too many LINCS warnings
    steps_ite=$((20000/steps))
    gmx mdrun -deffnm min -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps $steps
    for i in $(seq 2 $steps_ite)
    do
        echo ---------------------------
        echo ite $i out ouf $steps_ite
        echo ---------------------------
        gmx grompp -f mdp_files/min.mdp -c min.gro -p topol.top -o min.tpr -quiet
        gmx mdrun -deffnm min -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps $steps
        rm -f \#*
    done
fi

echo Potential | gmx energy -f min.edr -o potential.xvg -quiet
xmgrace potential.xvg

# vizualization script for pymol
cat << EOF > min.pml
remove all
load min.gro
remove resname PW
remove resname ION
extract PC, resname POPC
extract PS, resname POPS
extract CH, resname CHOL
color orange, min
set orthoscopic, on
show spheres
show cell
EOF
pymol min.pml

## make index file with LIP and SOL
cat << EOF > index_input
r POPC | r POPS | r CHOL
name 18 LIP
r PW | r ION
name 19 SOL
1 | 18 
q
EOF
gmx make_ndx -f min.gro -quiet < index_input
rm index_input

fi

######################################################################################################
### NPT ##############################################################################################
######################################################################################################

if [ $NPT -eq 1 ]
then

echo -----------------------------------------------------------
echo NPT run, with posres
echo -----------------------------------------------------------

cat << EOF > mdp_files/npt.mdp
refcoord_scaling         = com
integrator               = md
dt                       = 0.02 
nsteps                   = 5000000 ; 100 ns
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 50000 ; every 1 ns 
nstenergy                = 50000
nstxout-compressed       = 50000
cutoff-scheme            = Verlet
nstlist                  = 20
pbc                      = xyz
verlet-buffer-tolerance  = 0.005
coulombtype              = reaction-field 
rcoulomb                 = 1.1
epsilon_r                = 2.5
epsilon_rf               = 0
vdw_type                 = cutoff  
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1
tcoupl                   = v-rescale 
tc-grps                  = LIP SOL Protein
tau_t                    = 1.0 1.0 1.0
ref_t                    = 320 320 320
Pcoupl                   = parrinello-rahman
Pcoupltype               = semiisotropic
tau_p                    = 12.0
compressibility          = 3e-4 3e-4
ref_p                    = 1.0 1.0
gen_vel                  = no
constraints              = none 
constraint_algorithm     = Lincs
continuation             = yes
EOF
gmx grompp -f mdp_files/npt.mdp -c min.gro -r min.gro -p topol.top -n index.ndx -o npt.tpr -quiet
gmx mdrun -v -deffnm npt -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 10000 # short run to get .gro file 
gmx mdrun -v -deffnm npt -cpi npt.cpt -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset

## vizualization scripts for pymol
cat << EOF > npt.pml
remove all
load npt.gro
load_traj npt.xtc
remove resname PW
remove resname ION
extract PC, resname POPC
extract PS, resname POPS
extract CH, resname CHOL
color orange, npt
set orthoscopic, on
show spheres
show cell
EOF
pymol npt.pml

fi

######################################################################################################
### PULL #############################################################################################
######################################################################################################

if [ $PULL -eq 1 ]
then

## PULL

echo -----------------------------------------------------------
echo PULL peptide along z
echo -----------------------------------------------------------

cat << EOF > mdp_files/pull.mdp
integrator               = md
dt                       = 0.02 
nsteps                   = 40000000 ; 800 ns
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 5000 ; every 100 ps 
nstenergy                = 5000
nstxout-compressed       = 5000
cutoff-scheme            = Verlet
nstlist                  = 20
pbc                      = xyz
verlet-buffer-tolerance  = 0.005
coulombtype              = reaction-field 
rcoulomb                 = 1.1
epsilon_r                = 2.5
epsilon_rf               = 0
vdw_type                 = cutoff  
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1
tcoupl                   = v-rescale 
tc-grps                  = LIP SOL Protein
tau_t                    = 1.0 1.0 1.0
ref_t                    = 320 320 320
Pcoupl                   = no
gen_vel                  = no
constraints              = none 
constraint_algorithm     = Lincs
continuation             = yes
pull                     = yes
pull_ngroups             = 2
pull_ncoords             = 1
pull_group1_name         = LIP
pull_group2_name         = Protein
pull_coord1_type         = umbrella
pull_coord1_start        = yes
pull_coord1_rate         = 0.0001 
pull_coord1_k            = 500
pull_coord1_geometry     = direction
pull_coord1_dim          = N N Y ; pull along z
pull_coord1_vec          = 0 0 1
pull_coord1_groups       = 1 2
pull_nstfout             = 5000 ; should be equal to nstxout-compressed
pull_nstxout             = 5000 ; should be equal to nstxout-compressed
refcoord_scaling         = com
pull-group1-pbcatom      = -1 
EOF
gmx grompp -f mdp_files/pull.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o pull.tpr -quiet
gmx mdrun -v -deffnm pull -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 10000 # short run to get .gro file 
gmx mdrun -v -deffnm pull -cpi pull.cpt -px pull_pullx.xvg -pf pull_pullf.xvg -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset

# vizualization scripts for pymol
cat << EOF > pull.pml
remove all
load pull.gro
load_traj pull.xtc
remove resname PW
remove resname ION
extract PR, resname POPR
extract PC, resname POPC
extract PS, resname POPS
extract CH, resname CHOL
color red, PR
color orange, pull
set orthoscopic, on
show spheres
show cell
EOF
pymol pull.pml

fi

#####################################################################################################
## PUSH #############################################################################################
#####################################################################################################

if [ $PUSH -eq 1 ]
then

cp mdp_files/pull.mdp mdp_files/push.mdp
sed -i -e "s/pull_coord1_rate         = 0.0001/pull_coord1_rate         = -0.0001/g" mdp_files/push.mdp
gmx grompp -f mdp_files/push.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o push.tpr -quiet
gmx mdrun -v -deffnm push -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 10000 # short run to get .gro file 
gmx mdrun -v -deffnm push -cpi push.cpt -px push_pullx.xvg -pf push_pullf.xvg -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 800000

cp pull.pml push.pml
sed -i -e "s/pull/push/g" push.pml
pymol push.pml

fi

#####################################################################################################
## SETUP UMBRELLA ###################################################################################
#####################################################################################################


if [ $SETUP -eq 1 ]
then

cat << EOF > extract_frames.py
import numpy as np

file_in = 'pull_pullx.xvg'
file_out = 'extract_frames.ndx'

# generate output file
with open(file_out,'w') as f:
        f.write(' [ frames ] \n\n')

# find and add frames to output file
time,dist = np.genfromtxt(file_in,skip_header=17,skip_footer=3,unpack=True)   
# skip the last few lines - might be corrupted as pull is stopped with error when dist reach half the box size
min_dist = dist[0]
max_dist = 7.0       # protein pulled far enough from membrane to be converged
step_size_nm = 0.05 # stepsize in nm
dist_steps = np.arange(min_dist,max_dist,step_size_nm)

index_prev = 0
for d in dist_steps:
    indices = np.where(dist>d)
    index = indices[0][0] 
    frame = index + 1
    if index != index_prev:
        with open(file_out,'a') as f:
            f.write('%d\n' % frame)
    index_prev = index

with open(file_out,'a') as f:
        f.write('\n')
EOF

python extract_frames.py
echo 0 | gmx trjconv -f pull.xtc -s pull.tpr -n  -fr extract_frames.ndx -o extract_frames.xtc -quiet

cat << EOF > extract_frames_push.py
import numpy as np

file_in = 'push_pullx.xvg'
file_out = 'extract_frames_push.ndx'

# generate output file
with open(file_out,'w') as f:
        f.write(' [ frames ] \n\n')

# find and add frames to output file
time,dist = np.genfromtxt(file_in,skip_header=17,skip_footer=3,unpack=True)   
# skip the last few lines - might be corrupted as pull is stopped with error when dist reach half the box size
min_dist = dist[-1]
max_dist = dist[0]
step_size_nm = 0.05 # stepsize in nm
dist_steps = np.arange(min_dist,max_dist,step_size_nm)

index_prev = 0
for d in dist_steps:
    indices = np.where(dist<=d)
    index = indices[0][0] 
    frame = index + 1
    if index != index_prev:
        with open(file_out,'a') as f:
            f.write('%d\n' % frame)
    index_prev = index

with open(file_out,'a') as f:
        f.write('\n')
EOF

python extract_frames_push.py
echo 0 | gmx trjconv -f push.xtc -s push.tpr -n  -fr extract_frames_push.ndx -o extract_frames_push.xtc -quiet

Nframes_pull=-3; while read -r LINE; do (( Nframes_pull++ )); done < extract_frames.ndx
Nframes_push=-3; while read -r LINE; do (( Nframes_push++ )); done < extract_frames_push.ndx
Nframes=$((Nframes_push+Nframes_pull))
echo $Nframes_pull
echo $Nframes_push
echo $Nframes

cat << EOF > mdp_files/umbrella.mdp
integrator               = md
dt                       = 0.02 
nsteps                   = 50000000 ; 1000 ns
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 500000 ; every 10 ns 
nstenergy                = 500000
nstxout-compressed       = 500000
cutoff-scheme            = Verlet
nstlist                  = 20
pbc                      = xyz
verlet-buffer-tolerance  = 0.005
coulombtype              = reaction-field 
rcoulomb                 = 1.1
epsilon_r                = 2.5
epsilon_rf               = 0
vdw_type                 = cutoff  
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1
tcoupl                   = v-rescale 
tc-grps                  = LIP SOL Protein
tau_t                    = 1.0 1.0 1.0
ref_t                    = 320 320 320
Pcoupl                   = no
gen_vel                  = no
constraints              = none 
constraint_algorithm     = Lincs
continuation             = yes
pull                     = yes
pull_ngroups             = 2
pull_ncoords             = 1
pull_group1_name         = LIP
pull_group2_name         = Protein
pull_coord1_type         = umbrella
pull_coord1_start        = yes
pull_coord1_rate         = 0.0 
pull_coord1_k            = 500
pull_coord1_geometry     = direction
pull_coord1_dim          = N N Y ; pull along z
pull_coord1_vec          = 0 0 1
pull_coord1_groups       = 1 2
pull_nstfout             = 500000 ; should be equal to nstxout-compressed
pull_nstxout             = 500000 
refcoord_scaling         = com
pull-group1-pbcatom      = -1 
EOF

# GROMPP
mkdir -p $umb_dir
for frame in $(seq 1 $Nframes)
do

if [ $frame -le $Nframes_push ]
then

cat << EOF > frames_step.ndx
 [ frames ]
$((frame))
EOF

echo 0 | gmx trjconv -f extract_frames_push.xtc -s push.tpr -fr frames_step.ndx -o step_${frame}.gro -quiet

else

cat << EOF > frames_step.ndx
 [ frames ]
$((frame-Nframes_push))
EOF

echo 0 | gmx trjconv -f extract_frames.xtc -s pull.tpr -fr frames_step.ndx -o step_${frame}.gro -quiet

fi

rm frames_step.ndx

gmx grompp -f mdp_files/umbrella.mdp -c step_${frame}.gro -r step_${frame}.gro -n index.ndx -p topol.top -o $umb_dir/umbrella_step${frame}.tpr -maxwarn 1 -quiet
rm step_${frame}.gro

done

#make files for WHAM
cd $umb_dir
for frame in $(seq 1 $Nframes)
do
echo "umbrella_step${frame}.tpr" >> tpr_files_umbrella.dat
echo "umbrella_step${frame}_pullf.xvg" >> pullf_files_umbrella.dat
done
cd ..

# make script for computerome
rome_nodes=8
cat << EOF > $umb_dir/script.sh
#!/bin/sh
#PBS -W group_list=ku_00142 -A ku_00142
#PBS -N ${name}_bE_rep${rep}
#PBS -m n
#PBS -l nodes=1:thinnode:ppn=${rome_nodes}
#PBS -l mem=35gb
#PBS -l walltime=7:00:00:00
#PBS -t 1-$Nframes
echo Working dir is \$PBS_O_WORKDIR
cd \$PBS_O_WORKDIR
NPROCS="wc -l < \$PBS_NODEFILE"
echo This job has allocated \$NPROCS nodes
module purge
module load tools computerome_utils/2.0
module load gcc/12.1.0
module load cuda/toolkit/10.2.89
module load gromacs/2021.3-plumed
gmx_mpi mdrun -deffnm umbrella_step\${PBS_ARRAYID} -ntomp ${rome_nodes}
EOF

fi

#####################################################################################################
## RUN UMBRELLA #####################################################################################
#####################################################################################################


if [ $LOCAL -eq 1 ]
then

cd $umb_dir
for frame in $(seq 1 $Nframes)
do
echo -----------------------------------------------------------------------
echo Running umbrella test, frame no $frame out of $Nframes
echo -----------------------------------------------------------------------
# testrun for 2500000 steps = 50 ns 
#gmx mdrun -v -deffnm umbrella_step${frame} -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 2500000
gmx mdrun -v -deffnm umbrella_step${frame} -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset
done
cd ..

fi

if [ $ROME -eq 1 ]
then

# copy to computerome
scp -r $umb_dir andhaa@ssh.computerome.dk:/home/projects/ku_00142

fi

######################################################################################################
## ANALYZE UMBRELLA ##################################################################################
######################################################################################################

if [ $ANALYZE -eq 1 ]
then

# download from computerome
#scp -r andhaa@ssh.computerome.dk:/home/projects/ku_00142/$umb_dir .
scp -r andhaa@ssh.computerome.dk:/home/projects/ku_00142/$umb_dir/*_pull*.xvg $umb_dir

cd $umb_dir
#taskset -c 8 gmx wham -it tpr_files_umbrella.dat -if pullf_files_umbrella.dat -hist -temp 320 -b 500 -quiet -nBootstrap 200
taskset -c 8 gmx wham -it tpr_files_umbrella.dat -if pullf_files_umbrella.dat -hist -temp 320 -b 500 -bins 80 -min 1.5 -max 6 -quiet
cd ..

fi

######################################################################################################
### CLEAN UP #########################################################################################
######################################################################################################

## clean up
cd $umb_dir
rm -f \#*
cd ..
rm -f \#*
rm -f step*.pdb

cd ../..


