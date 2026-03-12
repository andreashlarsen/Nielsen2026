
## (de)select part of script
DOWNLOAD=1
PREPARE=1
NVT=1
PULL=1
SETUP=1
ROME=0
ANALYZE=0

## parallelisation
pinoffset=0
CPU=8
GPU=0

name=d22
rep=1
IDP=0

######################################################################################################
### WORKING DIR ######################################################################################
######################################################################################################

# prepare working dir
wdir=${name}/rep${rep}
mkdir -p $wdir/UMBRELLA
cp *.ipynb $wdir
cp umbrella_new-rf_base.mdp $wdir/UMBRELLA
cd $wdir
sed -i -e "s/NAME/$name/g" setup_umbrella.ipynb
sed -i -e "s/NAME/$name/g" umbrella_analysis.ipynb
sed -i -e "s/REP/$rep/g" setup_umbrella.ipynb
sed -i -e "s/REP/$rep/g" umbrella_analysis.ipynb
rome_folder=${name}_rep${rep}
echo rome folder: $rome_folder
mkdir -p $rome_folder

######################################################################################################
### DOWNLOAD #########################################################################################
######################################################################################################

if [ $DOWNLOAD -eq 1 ]
then
  ## download (and change) martini sctipts and itp files
  wget www.cgmartini.nl/images/tools/insane/insane.py
  #wget http://cgmartini.nl/images/parameters/ITP/martini_v2.2.itp
  wget http://cgmartini.nl/images/parameters/ITP/martini_v2.2refP.itp
  wget http://cgmartini.nl/images/parameters/ITP/martini_v2.0_ions.itp
  wget http://cgmartini.nl/images/parameters/lipids/Collections/martini_v2.0_lipids_all_201506.itp
  wget http://cgmartini.nl/images/parameters/lipids/PC/POPC/martini_v2.0_POPC_02.itp
  wget http://cgmartini.nl/images/tools/martinize/python3/martinize-2.6/martinize.py
  mv martini_v2.0_POPC_02.itp martini_v2.0_POPR.itp
  sed -i -e "s/POPC/POPR/g" martini_v2.0_POPR.itp
fi

######################################################################################################
### PREPARE ##########################################################################################
######################################################################################################

if [ $PREPARE -eq 1 ]
then

echo ############################ PREPARE ############################################################

# copy nvt.gro from buckle_hecate_x 
cp ../../Martini2/buckle_hecate_x/nvt.gro .
cp ../../Martini2/buckle_hecate_x/nvt.tpr .
cp ../../Martini2/buckle_hecate_x/index.ndx .

# extract membrane
echo 20 | gmx trjconv -f nvt.gro -s nvt.tpr -n index.ndx -o nvt_lip.gro -quiet

## protein
protein=../../Structures/$name.pdb

## martinize
if [ $IDP -eq 1 ]
then
  echo MARTINI - NO ELASTIC NETWORK
  if [ $name == d2_idp ]
  then
    echo 1 1 | gmx trjconv -f $protein -s $protein -pbc whole -o prot_pbc.pdb -quiet
    python3 martinize.py -v -f $protein -o topol.top -dssp dssp -ff martini22p -x cg.pdb 
  else
    python3 martinize.py -v -f $protein -o topol.top -dssp dssp -ff martini22p -x cg.pdb 
  fi
else
  python3 martinize.py -v -f $protein -o topol.top -dssp dssp -ff martini22p -x cg.pdb -elastic
fi
## convert pdb to gro
#echo 1 | gmx trjconv -f cg.pdb -s cg.pdb -o cg.gro -quiet
echo 1 0 | gmx trjconv -f cg.pdb -s cg.pdb -center -o cg.gro -quiet

## move protein close to membrane (which it will be merged with in the next step), and close to the correct leaflet (opposite POPRs)
echo 0 | gmx trjconv -f cg.gro -s cg.gro -o cg_trans.gro -trans 12 3 10 -quiet

## add lip gro file and prot gro file
head -n -1 cg_trans.gro > tmp0
tail -n +3 tmp0 > prot 
head -1 nvt_lip.gro > header
tail -1 nvt_lip.gro > footer
head -n -1 nvt_lip.gro > tmp1
tail -n +3 tmp1 > lip       
N_prot=$(sed -n '2p' cg_trans.gro)
echo N_prot $N_prot 
N_lip=$(sed -n '2p' nvt_lip.gro)
echo N_lip $N_lip
N_prot_lip=$((N_prot+N_lip))
echo N_prot_lip $N_prot_lip 
cat << EOF > N
$N_prot_lip
EOF
cat header N prot lip footer > nvt_lip_prot.gro
rm header N prot lip footer tmp0 tmp1

## renumber
gmx genconf -f nvt_lip_prot.gro -renumber -o nvt_lip_prot_renumber.gro -quiet

## get box size
tail -1 nvt_lip_prot_renumber.gro 

## solvate with insane (use same box size)
## peptide needs to be in box, check with pymol nvt_lip_prot_renumber.gro
python2 insane.py -f nvt_lip_prot_renumber.gro -o bilayer.gro -p topol.top -x 21 -y 10.0 -z 22.20881 -sol PW -salt 0.05

## change topology file
if [ $name == alps ] || [ $IDP -eq 1 ]
then
   prot_name=Protein
else
    prot_name=Protein_A
fi
cat << EOF > tmp1
#include "martini_v2.2refP.itp"
#include "martini_v2.0_ions.itp"
#include "martini_v2.0_lipids_all_201506.itp"
#include "martini_v2.0_POPR.itp"
#ifdef POSRES_POPR
#include "posres_popr.itp"
#endif
#include "$prot_name.itp"

EOF
tail -n +2 topol.top > tmp2
cat tmp1 tmp2 > topol.top
sed -i -e "s/Protein   /$prot_name /g" topol.top

## add lips to topol.top (from ../buckle_heacete_x/topol.top) - before Protein and water and ions
cat << EOF > tmp3
POPR           304
POPS           101
CHOL           101
POPC           300
POPS           100
CHOL           100
EOF
head -n +17 topol.top > tmp4
tail -n +18 topol.top > tmp5
cat tmp4 tmp3 tmp5 > topol.top
rm tmp1 tmp2 tmp3 tmp4 tmp5

## generate position restraint on PO4's of POPCs in lower leaflet (named POPR)
cat << EOF > posres_popr.itp
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   2    1         20         20         20
EOF

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
epsilon_r                = 15
epsilon_rf               = 0
vdw_type                 = cutoff  
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1
EOF
gmx grompp -f mdp_files/min.mdp -c bilayer.gro -p topol.top -o min.tpr -quiet -maxwarn 1 # non-matching names
gmx mdrun -deffnm min -v -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset

echo Potential | gmx energy -f min.edr -o potential.xvg -quiet
#xmgrace potential.xvg

## vizualization script for pymol: ensure that the peptide is placed on the right side of the membrane
cat << EOF > min.pml
remove all
load min.gro
remove resname PW
remove resname ION
extract PR, resname POPR
extract PC, resname POPC
extract PS, resname POPS
extract CH, resname CHOL
color red, PR
color orange, min
set orthoscopic, on
show spheres
show cell
EOF
#pymol min.pml

## make index file with LIP and SOL
cat << EOF > index_input
13 | 14 | 15 | 16
name 19 LIP
17 | 18
name 20 SOL
1 | 19
a PO4 & r POPC
name 22 PO4_POPC
r POPR | r POPC
q
EOF
gmx make_ndx -f min.gro -quiet < index_input
rm index_input

fi

######################################################################################################
### NVT ##############################################################################################
######################################################################################################

if [ $NVT -eq 1 ]
then

echo ############################ NVT ################################################################

## equilibrate with lower leaflet POPC's (named POPR) contrained 

echo -----------------------------------------------------------
echo NVT run, with posres
echo -----------------------------------------------------------

cat << EOF > mdp_files/nvt.mdp
define                   = -DPOSRES_POPR
refcoord_scaling         = com
integrator               = md
dt                       = 0.02 
nsteps                   = 5000000000 ; 100 us
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
ref_t                    = 310 310 310
Pcoupl                   = no
gen_vel                  = no
constraints              = none 
constraint_algorithm     = Lincs
continuation             = yes
EOF
gmx grompp -f mdp_files/nvt.mdp -c min.gro -r min.gro -p topol.top -n index.ndx -o nvt.tpr -quiet
gmx mdrun -v -deffnm nvt -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -nsteps 2000 # short run to get .gro file 
gmx mdrun -v -deffnm nvt -cpi nvt.cpt -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -nsteps 5000000 # 100 ns

## vizualization scripts for pymol
cat << EOF > nvt.pml
remove all
load nvt.gro
load_traj nvt.xtc
remove resname PW
remove resname ION
extract PR, resname POPR
extract PC, resname POPC
extract PS, resname POPS
extract CH, resname CHOL
color red, PR
color orange, nvt
set orthoscopic, on
show spheres
show cell
EOF
#pymol nvt.pml

## reduce trajectory size and remove solvent
echo 21 | gmx trjconv -s nvt.tpr -f nvt.gro -n index.ndx -o nvt_prot_lip.gro -quiet
echo 21 | gmx trjconv -s nvt.tpr -f nvt.xtc -n index.ndx -dt 100 -o nvt_prot_lip.xtc -quiet

## vizualization scripts for pymol
cat << EOF > nvt_prot_lip.pml
remove all
load nvt_prot_lip.gro
load_traj nvt_prot_lip.xtc
extract PR, resname POPR
extract PC, resname POPC
extract PS, resname POPS
extract CH, resname CHOL
color red, PR
color orange, nvt
set orthoscopic, on
show spheres
show cell
EOF
#pymol nvt_prot_lip.pml

fi

######################################################################################################
### PULL #############################################################################################
######################################################################################################

if [ $PULL -eq 1 ]
then

echo ############################ PULL ###############################################################

## rotate nvt to be like stroh2021 setup: peptide above membrane
# get box dim
tail -1 nvt.gro
gmx editconf -f nvt.gro -rotate 0 180 0 -box 21 10 22.20881 -o nvt_rot.gro -quiet

## PULL

echo -----------------------------------------------------------
echo PULL peptide along x
echo -----------------------------------------------------------

cat << EOF > mdp_files/pull.mdp
define                   = -DPOSRES_POPR
refcoord_scaling         = com
integrator               = md
dt                       = 0.02 
nsteps                   = 60000000 ; 1200 ns
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
ref_t                    = 310 310 310
Pcoupl                   = no
gen_vel                  = no
constraints              = none 
constraint_algorithm     = Lincs
continuation             = yes
pull                     = yes
pull_ngroups             = 1
pull_ncoords             = 1
pull_group1_name         = Protein
pull_coord1_type         = umbrella
pull_coord1_start        = yes
pull_coord1_rate         = -0.00003 ; negative sign as in Kai Stroh example
pull_coord1_k            = 500
pull_coord1_geometry     = direction-periodic
pull_coord1_dim          = Y N N ; pull along x
pull_coord1_vec          = 1 0 0
pull_coord1_groups       = 0 1
pull_coord1_origin       = 0 0 0
pull_nstfout             = 30300 ; should be equal to nstxout-compressed
pull_nstxout             = 30300 ; should be equal to nstxout-compressed
EOF
gmx grompp -f mdp_files/pull.mdp -c nvt_rot.gro -r nvt_rot.gro -p topol.top -n index.ndx -o pull.tpr -quiet
gmx mdrun -v -deffnm pull -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -nsteps 2000 # short run to get .gro file 
gmx mdrun -v -deffnm pull -cpi pull.cpt -px pull_pullx.xvg -pf pull_pullf.xvg -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset

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
#pymol pull.pml

# visual inspection: ensure that the peptide is placed on the right side of the membrane
# reduce pull sim to cover the curvatures you are intereted in:
# frame 
echo 0 | gmx trjconv -s pull.tpr -f pull.xtc -b 290000 -e 660000 -o pull_reduced.xtc -quiet

# vizualization script for pymol
cp pull.pml pull_reduced.pml
sed -i -e "s/load_traj pull.xtc/load_traj pull_reduced.xtc/g" pull_reduced.pml
#pymol pull_reduced.pml

fi

#####################################################################################################
## SETUP UMBRELLA ###################################################################################
#####################################################################################################

if [ $SETUP -eq 1 ]
then

echo ############################ SETUP UMB #########################################################

if [ $name == alps ] || [ $IDP -eq 1 ]
then
    echo Protein instead of Protein_A in setup_umbrella.ipynb
    sed -i -e "s/Protein_A/Protein/g" setup_umbrella.ipynb
fi
#jupyter-lab setup_umbrella.ipynb
jupyter nbconvert --to notebook --inplace --execute setup_umbrella.ipynb

## copy index file to UMBRELLA
cp index.ndx UMBRELLA

## get frames
cat << EOF > get_frames.py
import numpy as np

frames = np.genfromtxt('UMBRELLA/caught-output.dat',skip_header=1,usecols=[0],unpack=True)

f00 = open('tmp00','w')
f01 = open('tmp01','w')
f10 = open('tmp10','w')
f11 = open('tmp11','w')
fcr = open('tmpcr','w') # for running on computerome
f00.write('for i in ')
f01.write('for i in ')
f10.write('for i in ')
f11.write('for i in ')
fcr.write('#PBS -t ')

frame_prev = -1
N = 0
with open('tmp0','w') as f:
    f.write('for i in ')
    for frame in frames:
        if frame != frame_prev:
            f.write('%3d ' % frame)
            fcr.write('%d,' % frame)
            N += 1
            R = N%4 # modulos 4
            if R == 0:
               f00.write('%3d ' % frame) 
            elif R == 1:
               f01.write('%3d ' % frame) 
            elif R == 2:
               f10.write('%3d ' % frame) 
            elif R == 3:
               f11.write('%3d ' % frame)  
            else:
               print('ERROR: somethings wrong!!')
               exit()
        frame_prev = frame

f00.close()
f01.close()
f10.close()
f11.close()
fcr.close()

print('Total number of frames = %d' % N)

EOF
python3 get_frames.py

## make file for GROMPPing and execute it

cat << EOF > tmp1

cd UMBRELLA

EOF

cat << EOF > tmp2

do
    echo --------------------------------------
    echo GROMPP: generate umbrella\${i}.tpr 
    echo --------------------------------------

    gmx grompp -f umbrella\${i}.mdp -c conf\${i}.gro -r conf\${i}.gro -p ../topol.top -n index.ndx -o umbrella\${i}.tpr -quiet
done

cd ..

EOF

cat tmp1 tmp0 tmp2 > do_grompp.sh

echo ------------------
echo GROMPPing...
echo ------------------
bash do_grompp.sh

fi 

#####################################################################################################
## RUN UMBRELLA COMPUTEROME #########################################################################
#####################################################################################################

if [ $ROME -eq 1 ]
then

echo ############################ ROME ##############################################################

## make computerome script
rome_nodes=8
cat << EOF > tmpcr0
#!/bin/sh
#PBS -W group_list=ku_00142 -A ku_00142
#PBS -N ${name}_rep${rep}
#PBS -m n
#PBS -l nodes=1:thinnode:ppn=${rome_nodes}
#PBS -l mem=35gb
#PBS -l walltime=14:00:00:00
EOF

cat << EOF > tmpcr1

echo Working dir is \$PBS_O_WORKDIR
cd \$PBS_O_WORKDIR
NPROCS="wc -l < \$PBS_NODEFILE"
echo This job has allocated \$NPROCS nodes
module purge
module load tools computerome_utils/2.0
module load gcc/12.1.0
module load cuda/toolkit/10.2.89
module load gromacs/2021.3-plumed
gmx_mpi mdrun -deffnm umbrella\${PBS_ARRAYID} -ntomp ${rome_nodes}
EOF

cat tmpcr0 tmpcr tmpcr1 > script.sh
mv script.sh UMBRELLA

## make lumi script
cd $umb_dir
lumi_nodes=4
frames=$(ls umbrella*.tpr | sed 's/[^0-9]*\([0-9]*\).*/\1/' | tr '\n' ',' | sed 's/,$//')
cat << EOF >> lumi_script.sh
#!/bin/bash
#SBATCH --job-name=${protein}_rep${rep}
#SBATCH --account=project_465001110
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=$lumi_nodes
#SBATCH --mem=16G
#SBATCH --partition=small
#SBATCH --array=${frames}
module purge
module load LUMI/23.09
module load partition/L
module load GROMACS/2021.7-cpeGNU-23.09-CPU
gmx mdrun -deffnm umbrella\${SLURM_ARRAY_TASK_ID} -v -nt \${SLURM_CPUS_PER_TASK}
EOF
cd ..
scp -r $umb_dir lumi:/project/project_465001110

## run locally
for filename in umbrella*.tpr 
do
filename_without_extension="${filename%.tpr}"
if [ ! -f "${filename_without_extension}.xtc" ]
then
echo =============================================
echo gmx mdrun -v -deffnm $filename_without_extension -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset
echo =============================================
gmx mdrun -v -deffnm $filename_without_extension -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset
fi
done

## clean up
rm -f tmp* do_grompp.sh change_mdp.sh get_frames.py

## move files to computerome
cp UMBRELLA/*tpr $rome_folder
cp UMBRELLA/script.sh $rome_folder
#scp -r $rome_folder andhaa@ssh.computerome.dk:/home/projects/ku_00142
scp -r $rome_folder rome:/home/projects/ku_00142

## download results from computerome
cd UMBRELLA
#scp -r andhaa@ssh.computerome.dk:/home/projects/ku_00142/$rome_folder/* .
scp -r rome:/home/projects/ku_00142/$rome_folder/* .
cd ..

## download results from LUMI
cd UMBRELLA
scp -r lumi:/project/project_465001110/$rome_folder/* .
cd ..

## clean up
rm -r $rome_folder

fi 
######################################################################################################
## ANALYZE UM/home/projects/ku_00142/BRELLA ##################################################################################
######################################################################################################

if [ $ANALYZE -eq 1 ]
then

echo ############################ ANALYZE ############################################################

if [ $name == alps ] || [ $IDP -eq 1 ]
then
    sed -i -e "s/Protein_A/Protein/g" setup_umbrella.ipynb
fi
#jupyter-lab umbrella_analysis.ipynb
jupyter nbconvert --to notebook --inplace --execute umbrella_analysis.ipynb

# plot result
cat << EOF > plot_pmf.py
import numpy as np
import matplotlib.pyplot as plt
x,y,z = np.genfromtxt('pmf_100_1000.dat',unpack=True)
plt.plot(x,y)
plt.show()
EOF
python plot_pmf.py
fi

######################################################################################################
### CLEAN UP #########################################################################################
######################################################################################################

## clean up
cd UMBRELLA
rm -f \#*
rm slurm*.out
cd ..
rm -f \#*
rm -f step*.pdb

cd ../..
