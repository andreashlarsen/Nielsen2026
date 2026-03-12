
## (de)select part of script
DOWNLOAD=0
PREPARE=0
NPT=0
RESCALE=0
NVT=0
PULL=1

## parallelisation
pinoffset=0
CPU=8
GPU=1

######################################################################################################
### DOWNLOAD #########################################################################################
######################################################################################################

if [ $DOWNLOAD -eq 1 ]
then
  ## download (and change) martini sctipts and itp files
  wget www.cgmartini.nl/images/tools/insane/insane.py
  wget http://cgmartini.nl/images/parameters/ITP/martini_v2.2.itp
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

## protein
protein=../Structures/hecate.pdb

## box
x=30
y=10
z=15

## martinize
python3 martinize.py -v -f $protein -o topol.top -dssp dssp -ff martini22 -x cg.pdb -elastic

## random translation in x and y
cat << EOF > xy_trans.py
from random import random
print($x*(random()-0.5))
EOF
x_trans=$(python3 xy_trans.py)
cat << EOF > xy_trans.py
from random import random
print($y*(random()-0.5))
EOF
y_trans=$(python3 xy_trans.py)
rm xy_trans.py

## make bilayer with insane
python2 insane_trans.py -f cg.pdb -dm -2.0 -o bilayer.gro -p topol.top -x $x -y $y -z $z -l POPC:60 -l POPS:20 -l CHOL:20 -sol W:90 -sol WF:10 -salt 0.05 -rotate random -transx $x_trans -transy $y_trans

## change topology file
sed -i -e '0,/POPC /{s/POPC /POPR /}' topol.top
cat << EOF > tmp1
#include "martini_v2.2.itp"
#include "martini_v2.0_ions.itp"
#include "martini_v2.0_lipids_all_201506.itp"
#include "martini_v2.0_POPR.itp"
#ifdef POSRES_POPR
#include "posres_popr.itp"
#endif
#include "Protein_A.itp"

EOF
tail -n +2 topol.top > tmp2
cat tmp1 tmp2 > topol.top
rm tmp1 tmp2
sed -i -e 's/Protein   /Protein_A /g' topol.top

# position restraint on PO4's of POPCs in lower leaflet (named POPR)
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
gmx grompp -f mdp_files/min.mdp -c bilayer.gro -p topol.top -o min.tpr -quiet -maxwarn 1 # non-matching names: CNO - CN0
gmx mdrun -deffnm min -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset

echo Potential | gmx energy -f min.edr -o potential.xvg -quiet
xmgrace potential.xvg

# vizualization script for pymol
cat << EOF > min.pml
remove all
load min.gro
remove resname W
remove resname WF
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

## make index file with LIP and SOL
cat << EOF > index_input
13 | 14 | 15 | 16
name 20 LIP
17 | 18 | 19
name 21 SOL
1 | 20
a PO4 & r POPC
name 23 PO4_POPC
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

# equilibrate, NPT
cat << EOF > mdp_files/npt.mdp
integrator               = md
dt                       = 0.033 
nsteps                   = 3030303 ; 100 ns
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 303000 ; ca every 10 ns 
nstenergy                = 303000
nstxout-compressed       = 303000
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
tcoupl                   = v-rescale 
tc-grps                  = LIP SOL Protein
tau_t                    = 1.0  1.0 1.0
ref_t                    = 320 320 320
Pcoupl                   = berendsen
Pcoupltype               = anisotropic
tau_p                    = 12.0
compressibility          = 0.0  0.0  3e-4  0  0  0  
ref_p                    = 0.0  0.0  1.0   0  0  0
gen_vel                  = yes
gen_temp                 = 320
gen_seed                 = 473529
constraints              = none 
constraint_algorithm     = Lincs
EOF
gmx grompp -f mdp_files/npt.mdp -c min.gro -p topol.top -n index.ndx -o npt_0.tpr -quiet
gmx mdrun -v -deffnm npt_0 -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset

# vizualization scripts for pymol
cat << EOF > npt_0.pml
remove all
load npt_0.gro
load_traj npt_0.xtc
remove resname W
remove resname WF
remove resname ION
extract PR, resname POPR
extract PC, resname POPC
extract PS, resname POPS
extract CH, resname CHOL
color red, PR
color orange, npt_0
set orthoscopic, on
show spheres
show cell
EOF

fi

######################################################################################################
### RESCALE ##########################################################################################
######################################################################################################

n_ite=10 # needed by NVT, therefore outside if statement

if [ $RESCALE -eq 1 ]
then

vel=1 # velocities in gro file? (1: yes, 0= no)
scale_final=0.7
cat << EOF > calc_scale.py
print(${scale_final}**(1/$n_ite))
EOF
scale_x=$(python3 calc_scale.py)
rm calc_scale.py
echo scale_final = $scale_final
echo scale_x = $scale_x

## python script for rescale x
cat << EOF > scale_gro_file.py
import numpy as np

filename_in = input('enter name of input .gro file:  ')
filename_out = input('enter name of output .gro file:  ')
scale_x = float(input('enter scale factor x:  '))
scale_y = float(input('enter scale factor y:  '))
scale_z = float(input('enter scale factor z:  '))
vel = int(input('gro file include velocities? (0 or 1):  '))

print('\ninput file: %s' % filename_in)
print('scaling factor x: %2.2f' % scale_x)
print('scaling factor y: %2.2f' % scale_y)
print('scaling factor z: %2.2f' % scale_z)
print('output file: %s' % filename_out)

n_vel = 0
if vel:
    n_vel = 3

f_in = open(filename_in,'r')
f_out = open(filename_out,'w')

line = f_in.readline()
while line:
    line_split = line.split()
    n = len(line_split)
    if n in [3,5+n_vel,6+n_vel,7+n_vel]:
        idx = n-3-n_vel
        x,y,z = float(line_split[idx]),float(line_split[idx+1]),float(line_split[idx+2])
        x_scale,y_scale,z_scale = x*scale_x,y*scale_y,z*scale_z
        if n == 3:
            new_line = '%10.5f%10.5f%10.5f\n' % (x_scale,y_scale,z_scale)
        else:
            if vel:
                new_line = '%s%8.3f%8.3f%8.3f%s' %  (line[0:20],x_scale,y_scale,z_scale,line[44:])
            else:
                new_line = '%s%8.3f%8.3f%8.3f\n' %  (line[0:20],x_scale,y_scale,z_scale)
    else:
        new_line = line
    f_out.write(new_line)
    line = f_in.readline()
f_in.close()
f_out.close()
EOF
           
## rescale y-coordinate iteratively, equilibration at each step
for i in $(seq 1 $n_ite)
do
  echo -------------------------
  echo iteration $i
  echo -------------------------
  
  i_prev=$((i-1))
  
  ## input file for rescaling y-coordinate
  cat << EOF > input_scale
npt_$i_prev.gro
npt_${i_prev}_scale.gro
$scale_x
1.0
1.0
$vel
EOF
 
  # rescale
  python3 scale_gro_file.py < input_scale
  
  # equilibrate
  gmx grompp -f mdp_files/npt.mdp -c npt_${i_prev}_scale.gro -p topol.top -n index.ndx -o npt_$i.tpr -quiet
  gmx mdrun -v -deffnm npt_$i -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -nsteps 1010101 # 33 ns 
  
  # vizualization scripts for pymol
  cat << EOF > npt_$i.pml
remove all
load npt_$i.gro
load_traj npt_$i.xtc
remove resname W
remove resname WF
remove resname ION
extract PR, resname POPR
extract PC, resname POPC
extract PS, resname POPS
extract CH, resname CHOL
color red, PR
color orange, npt_$i
set orthoscopic, on
show spheres
show cell
EOF

done

fi

######################################################################################################
### NVT unconstrained?? ##############################################################################
######################################################################################################

# make extra NVT step here - without  POPR constrained - and extract frame for next step: NVT with constrained POPR
# in order to select frame with smallest deviation between membrane and the analytical shape (stroh2021, pageF)
#

######################################################################################################
### NVT ##############################################################################################
######################################################################################################

if [ $NVT -eq 1 ]
then

## equilibrate with lower leaflet POPC's (named POPR) contrained 

echo -----------------------------------------------------------
echo NVT run, with posres
echo -----------------------------------------------------------

cat << EOF > mdp_files/nvt.mdp
define                   = -DPOSRES_POPR
refcoord_scaling         = com
integrator               = md
dt                       = 0.033 
nsteps                   = 3030303030 ; 100 us
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 303000 ; every 10 ns 
nstenergy                = 303000
nstxout-compressed       = 303000
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
tcoupl                   = v-rescale 
tc-grps                  = LIP SOL Protein
tau_t                    = 1.0 1.0 1.0
ref_t                    = 320 320 320
Pcoupl                   = no
gen_vel                  = no
constraints              = none 
constraint_algorithm     = Lincs
continuation             = yes
EOF
gmx grompp -f mdp_files/nvt.mdp -c npt_$n_ite.gro -r npt_$n_ite.gro -p topol.top -n index.ndx -o nvt.tpr -quiet
gmx mdrun -v -deffnm nvt -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -nsteps 10000 # short run to get .gro file 
gmx mdrun -v -deffnm nvt -cpi nvt.cpt -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset 

## vizualization scripts for pymol
cat << EOF > nvt.pml
remove all
load nvt.gro
load_traj nvt.xtc
remove resname W
remove resname WF
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

## reduce trajectory size and remove solvent
echo 22 | gmx trjconv -s nvt.tpr -f nvt.gro -n index.ndx -o nvt_prot_lip.gro -quiet
echo 22 | gmx trjconv -s nvt.tpr -f nvt.xtc -n index.ndx -dt 100 -o nvt_prot_lip.xtc -quiet

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

## density map: protein
echo 1 | gmx densmap -f nvt.xtc -n index.ndx -b 2000000 -n1 90 -n2 1 -od densmap_prot -quiet # exclude first 2 us (biased by initial position)
echo 1 | gmx densmap -f nvt.xtc -n index.ndx -b 2000000 -n1 90 -n2 30 -o densmap_prot -quiet # exclude first 2 us (biased by initial position)
xv densmap_prot.xpm
## control density map: po4 from popc
echo 23 | gmx densmap -f nvt.xtc -n index.ndx -b 2000000 -n1 90 -n2 1 -od densmap_popc -quiet # exclude first 2 us (biased by initial position)
echo 23 | gmx densmap -f nvt.xtc -n index.ndx -b 2000000 -n1 90 -n2 30 -o densmap_popc -quiet # exclude first 2 us (biased by initial position)
xv densmap_popc.xpm
## plot ratio
cat << EOF > densmap.py
import numpy as np
import matplotlib.pyplot as plt

skip=1
x,dens_prot = np.genfromtxt('densmap_prot.dat',skip_header=skip,skip_footer=0,usecols=[0,1],unpack=True)
dens_popc = np.genfromtxt('densmap_popc.dat',skip_header=skip,skip_footer=0,usecols=[1],unpack=True)

freq_prot = dens_prot/np.mean(dens_prot)
freq_popc = dens_popc/np.mean(dens_popc)

n = len(freq_prot)
ratio = np.zeros(n)
for i in range(n):
    if freq_popc[i] > 0:   
        ratio[i] = freq_prot[i]/freq_popc[i]

with open('ratio.dat','w') as f:
    for i in range(n):
        f.write('%f %f\n' % (x[i],ratio[i]))

plt.plot(x,ratio,linestyle='none',marker='.',color='black',label='ratio')
plt.plot(x,freq_popc,linestyle='none',marker='.',color='blue',label='frequency, popc')
plt.plot(x,freq_prot,linestyle='none',marker='.',color='orange',label='frequence, protein')
plt.plot(x,np.ones(n),color='grey',linestyle='--')
plt.plot(x,np.zeros(n),color='grey',linestyle='--')

plt.xlim(0,np.amax(x))
plt.legend(frameon=False)

G = -np.log(ratio)
plt.plot(x,G)
#plt.show()
plt.show()

EOF
python densmap.py

fi 

######################################################################################################
### PULL #############################################################################################
######################################################################################################

if [ $PULL -eq 1 ]
then

## rotate nvt to be like stroh2021 setup: peptide above membrane
# get box dim
tail -1 nvt.gro
gmx editconf -f nvt.gro -rotate 0 180 0 -box 21 10 22.20881 -o nvt_rot.gro -quiet
#gmx editconf -f nvt_rot.gro -translate -6 0 0 -o nvt_rot_trans.gro -quiet

## PULL

echo -----------------------------------------------------------
echo PULL peptide along x
echo -----------------------------------------------------------

cat << EOF > mdp_files/pull.mdp
define                   = -DPOSRES_POPR
refcoord_scaling         = com
integrator               = md
dt                       = 0.033 
; nsteps                   = 15151515 ; 500 ns
nsteps                   = 24242424 ; 800 ns
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 30300 ; every 1 ns 
nstenergy                = 30300
nstxout-compressed       = 30300
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
gmx mdrun -v -deffnm pull -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset -nsteps 10000 # short run to get .gro file 
gmx mdrun -v -deffnm pull -cpi pull.cpt -px pull_pullx.xvg -pf pull_pullf.xvg -quiet -pin on -ntomp $CPU -ntmpi 1 -gpu_id $GPU -pinoffset $pinoffset

# vizualization scripts for pymol
cat << EOF >pull.pml
remove all
load pull.gro
load_traj pull.xtc
remove resname W
remove resname WF
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

# reduce pull sim to cover the curvatures you are intereted in
# frame 
echo 0 | gmx trjconv -s pull.tpr -f pull.xtc -b 120000 -e 520000 -o pull_reduced.xtc -quiet

# vizualization script for pymol
cp pull.pml pull_reduced.pml
sed -i -e "s/load_traj pull.xtc/load_traj pull_reduced.xtc/g" pull_reduced.pml

fi 

#####################################################################################################
## SETUP UMBRELLA ###################################################################################
#####################################################################################################

jupyter-lab setup_umbrella_xdir.ipynb

#####################################################################################################
## RUN UMBRELLA #####################################################################################
#####################################################################################################

for i in 00 01 10 11
do
rm nohup_mdrun$i.out
nohup bash do_mdrun$i.sh > nohup_mdrun$i.out &
done

# vizualize
bash pymol.sh
 
######################################################################################################
## ANALYZE UMBRELLA ##################################################################################
######################################################################################################

jupyter-lab umbrella_analysis.ipynb

bash wham.sh
 
######################################################################################################
### CLEAN UP #########################################################################################
######################################################################################################

## clean up
rm -f \#*
rm -f step*.pdb

