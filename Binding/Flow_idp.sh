
protein=kor2 # hecate kor d1 d2 b2 atr div4a OR kor2 or d22
rep=3
pinoffset=0
CPU=8

DOWNLOAD=1 # download martini files
PUSH=1

#### PREPARE ####
## make idp ff availble in gromacs
# cp -r charmm36IDPSFF.ff /usr/local/gromacs/share/gromacs/top
# rm -r charmm36IDPSFF.ff

mkdir -p idp_$protein
cd idp_$protein

## make gro file
cp ../${protein}/${protein}.pdb .
printf "9\n1" | gmx pdb2gmx -f ${protein}.pdb -o ${protein}.gro -ignh -quiet
gmx editconf -f ${protein}.gro -o ${protein}_box.gro -c -d 1.0 -bt cubic -quiet

## solvate
gmx solvate -cp ${protein}_box.gro -cs spc216.gro -o ${protein}_solv.gro -p topol.top -quiet

## ions
mkdir -p mdp_files
cat << EOF > mdp_files/ions.mdp
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
EOF
gmx grompp -f mdp_files/ions.mdp -c ${protein}_solv.gro -p topol.top -o ions.tpr -quiet
echo SOL | gmx genion -s ions.tpr -o ${protein}_solv_ions.gro -p topol.top -pname NA -nname CL -neutral -quiet

## minimization, AA
cat << EOF > mdp_files/min.mdp
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 200000         ; Maximum number of (minimization) steps to perform
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
EOF
gmx grompp -f mdp_files/min.mdp -c ${protein}_solv_ions.gro -p topol.top -o min.tpr -quiet
gmx mdrun -v -deffnm min -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -quiet
echo Potential | gmx energy -f min.edr -o potential.xvg -quiet
xmgrace potential.xvg

## NVT, AA
cat << EOF > mdp_files/nvt.mdp
define                  = -DPOSRES  
integrator              = md        
nsteps                  = 500000     ; 1 ns
dt                      = 0.002     ; 2 fs
nstxout                 = 0       
nstvout                 = 0       
nstfout                 = 0       
nstenergy               = 5000       
nstlog                  = 5000   
continuation            = no     
constraint_algorithm    = lincs   
constraints             = h-bonds  
lincs_iter              = 1    
lincs_order             = 4     
cutoff-scheme           = Verlet 
nstlist                 = 10      
rcoulomb                = 1.0      
rvdw                    = 1.0       
coulombtype             = PME       
pme_order               = 4         
fourierspacing          = 0.16      
tcoupl                  = V-rescale            
tc-grps                 = Protein Non-Protein   
tau_t                   = 0.1     0.1
ref_t                   = 300     300 
pcoupl                  = no        
pbc                     = xyz       
gen_vel                 = yes       
gen_temp                = 300       
gen_seed                = -1        
EOF
gmx grompp -f mdp_files/nvt.mdp -c min.gro -r min.gro -p topol.top -o nvt.tpr -quiet
gmx mdrun -v -deffnm nvt -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -quiet
echo Temperature | gmx energy -f nvt.edr -o temperature.xvg -quiet
xmgrace temperature.xvg

# NPT, AA
cat << EOF > mdp_files/npt.mdp
define                  = -DPOSRES 
integrator              = md        
nsteps                  = 500000     ; 2 * 500000 = 1 ns
dt                      = 0.002     ; 2 fs
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstenergy               = 5000
nstlog                  = 5000     
nstxout-compressed      = 5000
continuation            = yes        
constraint_algorithm    = lincs      
constraints             = h-bonds   
lincs_iter              = 1         
lincs_order             = 4         
cutoff-scheme           = Verlet    
nstlist                 = 10        
rcoulomb                = 1.0       
rvdw                    = 1.0       
coulombtype             = PME       
pme_order               = 4         
fourierspacing          = 0.16      
tcoupl                  = V-rescale             
tc-grps                 = Protein Non-Protein   
tau_t                   = 0.1     0.1           
ref_t                   = 300     300           
pcoupl                  = Parrinello-Rahman     
pcoupltype              = isotropic             
tau_p                   = 2.0                   
ref_p                   = 1.0                   
compressibility         = 4.5e-5                
refcoord_scaling        = com
pbc                     = xyz       
gen_vel                 = no       
EOF
gmx grompp -f mdp_files/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -quiet
gmx mdrun -v -deffnm npt -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -quiet
echo Pressure | gmx energy -f npt.edr -o pressure.xvg -quiet
xmgrace pressure.xvg
echo Density | gmx energy -f npt.edr -o density.xvg -quiet
xmgrace density.xvg

## production, AA
cat << EOF > mdp_files/prod.mdp
integrator              = md        
nsteps                  = 500000000     ; 1 us
dt                      = 0.002     ; 2 fs
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstenergy               = 500000
nstlog                  = 500000       
nstxout-compressed      = 500000 ; 1 ns
compressed-x-grps       = System
continuation            = yes                             
constraint_algorithm    = lincs                            
constraints             = h-bonds   
lincs_iter              = 1
lincs_order             = 4
cutoff-scheme           = Verlet    
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16      
tcoupl                  = V-rescale             
tc-grps                 = Protein Non-Protein   
tau_t                   = 0.1     0.1
ref_t                   = 300     300
pcoupl                  = Parrinello-Rahman     
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
gen_vel                 = no
EOF
gmx grompp -f mdp_files/prod.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o prod.tpr -quiet
gmx mdrun -v -deffnm prod -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -quiet -nsteps 5000
gmx mdrun -v -deffnm prod -cpi prod.cpt -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -quiet

## vizualization PyMOL
echo Protein Protein | gmx trjconv -f prod.gro -s prod.tpr -pbc whole -center -o prod_pbc.gro -quiet
echo Protein Protein | gmx trjconv -f prod.xtc -s prod.tpr -pbc whole -center -o prod_pbc.xtc -quiet

cat << EOF > prod.pml
load prod_pbc.gro
load_traj prod_pbc.xtc
set orthoscopic, on
show cell
EOF
pymol prod.pml

## monitor alpha helix structure with plumed
echo 1 1 | gmx trjconv -f prod_pbc.gro -s prod_pbc.gro -pbc mol -center -o prod_prot.pdb -quiet
cp prod_prot.pdb prod_prot_OX.pdb
sed -i -e 's/OT1/O  /g' prod_prot_OX.pdb
sed -i -e 's/ OT2/OXT1/g' prod_prot_OX.pdb
mkdir -p plumed
if [ $protein = hecate ]
then
Natoms=286
Nres=15
elif [ $protein = kor ]
then
Natoms=212
Nres=12
elif [ $protein = kor2 ]
then
Natoms=349
Nres=20
elif [ $protein = d2 ]
then
Natoms=254
Nres=14
elif [ $protein = d22 ]
then
Natoms=274
Nres=15
elif [ $protein = b2 ]
then
Natoms=303
Nres=18
elif [ $protein = atr ]
then
Natoms=369
Nres=20
elif [ $protein = div4a ]
then
Natoms=206
Nres=12
elif [ $protein = d1 ]
then
Natoms=345
Nres=22
elif [ $protein = hec23 ]
then
Natoms=414
Nres=23
else
echo ERROR: protein is not in list! 
fi
cat << EOF > plumed/alpha_helix.dat
## uncomment if you want to restart (continue) a previous run
RESTART
UNITS LENGTH=nm ENERGY=kj/mol TIME=ns
## make sure peptide is whole
MOLINFO STRUCTURE=prod_prot_OX.pdb
WHOLEMOLECULES ENTITY0=1-$Natoms
## define secondary structure collective variable
ah: ALPHARMSD RESIDUES=1-$Nres TYPE=OPTIMAL
PRINT STRIDE=1 ARG=ah FILE=ALPHA_HELIX
EOF
rm ALPHA_HELIX
plumed driver --mf_xtc prod_pbc.xtc --plumed plumed/alpha_helix.dat

## plot AH structure
cat << EOF > plot_AH_content.py
import numpy as np
import matplotlib.pyplot as plt
PLOT=int(input())
filename = 'ALPHA_HELIX'
time,AH = np.genfromtxt(filename,skip_header=2,usecols=[0,1],unpack=True)
time *= 100 # error in plumed units - I don't understand why...?!
# CV measures helical content for each 6-residue seqments. max value is 1 for each seqment. There are N-5 6res seqments in a Nres peptide, so:
AH_ideal = $Nres-5
AH_prc = AH/AH_ideal*100
AH_ini = AH_prc[0]
idx= np.argmin(AH_prc)
print(idx)
if PLOT:
    print('min AH content at time %1.1f ns' % time[idx])
    print('min AH content is %1.2f' % AH_prc[idx])
    plt.plot(time,AH_prc,color='black')
    plt.xlabel('time,ns')
    plt.ylabel('%s alpha-helix compared to ideal AH' % '%')
    plt.xlim(0,time[-1])
    plt.ylim(0,100)
    plt.savefig('AH_content.png')
    plt.show()
EOF
echo 1 | python plot_AH_content.py

## extract frames with lowest AH content
min_AH=$(echo 0 | python plot_AH_content.py)
cat << EOF > frame.ndx 
 [ frames ]
$min_AH
EOF
#echo 1 | gmx trjconv -f prod.xtc -s prod.tpr -fr frame.ndx -pbc nojump -o min_AH.pdb -quiet 
echo System | gmx trjconv -f prod_pbc.xtc -s prod_pbc.gro -fr frame.ndx -o min_AH.pdb -quiet 

############################################################################################
################################# CG #######################################################
############################################################################################

##prepare CG simulations
mkdir -p CG_rep$rep
cd CG_rep$rep
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

cp ../CG_rep1/insane.py .
cp ../CG_rep1/martini_v2.2refP.itp .
cp ../CG_rep1/martini_v2.0_ions.itp .
cp ../CG_rep1/martini_v2.0_lipids_all_201506.itp .
cp ../CG_rep1/martini_v2.0_POPC_02.itp .
cp ../CG_rep1/martinize.py .

## coarse-grain pdb
python3 martinize.py -v -f ../min_AH.pdb -o topol.top -dssp dssp -ff martini22p -x cg.pdb
echo 1 1 | gmx trjconv -f cg.pdb -s cg.pdb -center -o cg_cent.pdb -quiet
python2 insane.py -f cg_cent.pdb -o bilayer.gro -dm 3.3 -p topol.top -x 10 -y 10 -z 25 -l POPC:60  -l POPS:20 -l CHOL:20 -sol PW -salt 0.05
## change topology file
#sed -i -e '0,/POPC /{s/POPC /POPR /}' topol.top
cat << EOF > tmp1
#include "martini_v2.2refP.itp"
#include "martini_v2.0_ions.itp"
#include "martini_v2.0_lipids_all_201506.itp"
#include "Protein.itp"
EOF
tail -n +2 topol.top > tmp2
cat tmp1 tmp2 > topol.top
rm tmp1 tmp2

## minimization, CG
mkdir -p mdp_files
cat << EOF > mdp_files/min.mdp
integrator               = steep
nsteps                   = 5000
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
for i in {0..50}
do
echo iteration $i of 50
gmx mdrun -deffnm min -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 20 -nobackup
gmx grompp -f mdp_files/min.mdp -c min.gro -p topol.top -o min.tpr -quiet -nobackup
done
rm -f \#*
rm step*pdb
for i in {0..50}
do
echo iteration $i of 50
gmx mdrun -deffnm min -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 100 -nobackup
gmx grompp -f mdp_files/min.mdp -c min.gro -p topol.top -o min.tpr -quiet -nobackup
done
rm -f \#*
rm step*pdb
for i in {0..20}
do
echo iteration $i of 20
gmx mdrun -deffnm min -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 500 -nobackup
gmx grompp -f mdp_files/min.mdp -c min.gro -p topol.top -o min.tpr -quiet -nobackup
done
rm -f \#*
rm step*pdb
gmx mdrun -deffnm min -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nobackup
echo Potential | gmx energy -f min.edr -o potential.xvg -quiet -nobackup
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

echo -----------------------------------------------------------
echo NPT, CG, with posres and slowly increasing timesteps
echo -----------------------------------------------------------

cat << EOF > mdp_files/npt.mdp
integrator               = md
dt                       = 0.002
nsteps                   = 50000000 ; 1000 ns with final timestep 0.02
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 50000 
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
gen_vel                  = yes
gen_temp                 = 320
constraints              = none
constraint_algorithm     = Lincs
continuation             = no
EOF

gmx grompp -f mdp_files/npt.mdp -c min.gro -r min.gro -p topol.top -n index.ndx -o npt.tpr -quiet -nobackup -maxwarn 1
gmx mdrun -v -deffnm npt -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 20000 -nobackup
sed -i -e 's/continuation             = no/continuation             = yes/g' mdp_files/npt.mdp
sed -i -e 's/gen_vel                  = yes/gen_vel                  = no/g' mdp_files/npt.mdp
sed -i -e 's/dt                       = 0.002/dt                       = 0.003/g' mdp_files/npt.mdp
gmx grompp -f mdp_files/npt.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o npt.tpr -quiet -nobackup
gmx mdrun -v -deffnm npt -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 50000 -nobackup
sed -i -e 's/dt                       = 0.003/dt                       = 0.004/g' mdp_files/npt.mdp
gmx grompp -f mdp_files/npt.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o npt.tpr -quiet -nobackup
gmx mdrun -v -deffnm npt -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 50000 -nobackup
sed -i -e 's/dt                       = 0.004/dt                       = 0.005/g' mdp_files/npt.mdp
gmx grompp -f mdp_files/npt.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o npt.tpr -quiet -nobackup
gmx mdrun -v -deffnm npt -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 50000 -nobackup
sed -i -e 's/dt                       = 0.005/dt                       = 0.007/g' mdp_files/npt.mdp
gmx grompp -f mdp_files/npt.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o npt.tpr -quiet -nobackup
gmx mdrun -v -deffnm npt -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 50000 -nobackup
sed -i -e 's/dt                       = 0.007/dt                       = 0.009/g' mdp_files/npt.mdp
gmx grompp -f mdp_files/npt.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o npt.tpr -quiet -nobackup
gmx mdrun -v -deffnm npt -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 50000 -nobackup
sed -i -e 's/dt                       = 0.009/dt                       = 0.011/g' mdp_files/npt.mdp
gmx grompp -f mdp_files/npt.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o npt.tpr -quiet -nobackup
gmx mdrun -v -deffnm npt -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 50000 -nobackup
sed -i -e 's/dt                       = 0.011/dt                       = 0.013/g' mdp_files/npt.mdp
gmx grompp -f mdp_files/npt.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o npt.tpr -quiet -nobackup
gmx mdrun -v -deffnm npt -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 50000 -nobackup
sed -i -e 's/dt                       = 0.013/dt                       = 0.015/g' mdp_files/npt.mdp
gmx grompp -f mdp_files/npt.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o npt.tpr -quiet -nobackup
gmx mdrun -v -deffnm npt -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 50000 -nobackup
sed -i -e 's/dt                       = 0.015/dt                       = 0.020/g' mdp_files/npt.mdp -nobackup
gmx grompp -f mdp_files/npt.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o npt.tpr -quiet -nobackup
gmx mdrun -v -deffnm npt -cpi npt.cpt -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nobackup

## apply upper roof with plumed
#if [ $protein = hecate ]
#then
#NCGatoms=42
#NCGres=14
#Nliplast=3786
#elif [ $protein = kor ]
#then
#NCGatoms=35
#NGCres=11
#Nliplast=3779
#else
#echo ERROR: protein is not in list!
#fi
#mkdir -p plumed_files
#cat << EOF > plumed_files/plumed_${protein}.dat
#UNITS LENGTH=nm ENERGY=kj/mol TIME=ns
#protein: GROUP ATOMS=1-$NCGatoms
#com_prot: COM ATOMS=protein
#bilayer: GROUP ATOMS=$((NCGatoms+1))-$Nliplast
#com_lip: COM ATOMS=bilayer
#dist_prot_lip: DISTANCE ATOMS=com_lip,com_prot COMPONENTS
#roof: UPPER_WALLS ARG=dist_prot_lip.z AT=7.00 KAPPA=50 EXP=2 EPS=1 OFFSET=0
#PRINT STRIDE=10000 ARG=roof.bias FILE=COLVAR
#EOF
#gmx mdrun -v -deffnm npt -cpi npt.cpt -plumed plumed_files/plumed_${protein}.dat -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset

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
nstlog                   = 10000 ; every 200 ps
nstenergy                = 10000
nstxout-compressed       = 10000
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
pull_nstfout             = 10000 ; should be equal to nstxout-compressed
pull_nstxout             = 10000 ; should be equal to nstxout-compressed
refcoord_scaling         = com
pull-group1-pbcatom      = -1
EOF
gmx grompp -f mdp_files/pull.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o pull.tpr -quiet
gmx mdrun -v -deffnm pull -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 2000 # short run to get .gro file
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

if [ $PUSH -eq 1 ]
then
cp mdp_files/pull.mdp mdp_files/push.mdp
sed -i -e "s/pull_coord1_rate         = 0.0001/pull_coord1_rate         = -0.0001/g" mdp_files/push.mdp
gmx grompp -f mdp_files/push.mdp -c npt.gro -r npt.gro -p topol.top -n index.ndx -o push.tpr -quiet
gmx mdrun -v -deffnm push -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 2000 # short run to get .gro file
gmx mdrun -v -deffnm push -cpi push.cpt -px push_pullx.xvg -pf push_pullf.xvg -quiet -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -nsteps 500000
cp pull.pml push.pml
sed -i -e "s/pull/push/g" push.pml
pymol push.pml
fi

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
umb_dir=umb_idp_${protein}_bE_rep${rep}
mkdir -p $umb_dir
for frame in $(seq 1 $Nframes)
do
echo ===================
echo FRAME $frame
echo ===================
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

# make script for LUMI
lumi_nodes=6
cat << EOF > $umb_dir/script.sh
#!/bin/bash
#SBATCH --job-name=${protein}_idp
#SBATCH --account=project_465000460
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=$lumi_nodes
#SBATCH --mem=16G
#SBATCH --partition=small
#SBATCH --array=1-${Nframes}%100
module purge
module load LUMI/22.08
module load partition/L
module load GROMACS/2021.4-cpeCray-22.08-PLUMED-2.7.4-cray-python-3.9.12.1-CPU
gmx mdrun -deffnm umbrella_step\${SLURM_ARRAY_TASK_ID} -v -nt \${SLURM_CPUS_PER_TASK} -nsteps 25000000 # 500 ns
EOF
scp -r $umb_dir lumi:/project/project_465000460

# make script for computerome
rome_nodes=6
cat << EOF > $umb_dir/script.sh
#!/bin/sh
#PBS -W group_list=ku_00142 -A ku_00142
#PBS -N ${protein}_idp_bE_rep${rep}
#PBS -m n
#PBS -l nodes=1:thinnode:ppn=${rome_nodes}
#PBS -l mem=35gb
#PBS -l walltime=5:00:00:00
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
gmx_mpi mdrun -deffnm umbrella_step\${PBS_ARRAYID} -ntomp ${rome_nodes} -nsteps 25000000 # 500 ns
EOF
scp -r $umb_dir andhaa@ssh.computerome.dk:/home/projects/ku_00142

# download from lumi
scp -r lumi:/project/project_465000460/$umb_dir .

# download from computerome
scp -r andhaa@ssh.computerome.dk:/home/projects/ku_00142/$umb_dir .

# run locally
for step in {38..75}
do
echo ---------------------------------
echo umbrella step $step out of 75
echo ---------------------------------
gmx mdrun -v -deffnm umbrella_step$step -pin on -ntomp $CPU -ntmpi 1 -pinoffset $pinoffset -quiet -nsteps 25000000 # 500 ns 
done

if [ $ANALYZE -eq 1 ]
then
cd $umb_dir
taskset -c 10 gmx wham -it tpr_files_umbrella.dat -if pullf_files_umbrella.dat -hist -temp 320 -b 500 -quiet -nBootstrap 200
xmgrace bsResult.xvg
cd ..

fi

## clean up
cd $umb_dir
rm -f \#*
cd ..
rm -f \#*
rm step*.pdb

## parent dir
cd ..
