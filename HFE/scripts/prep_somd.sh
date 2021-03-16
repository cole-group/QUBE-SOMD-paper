# this util script will prep the somd systems to run qube
# this should be ran from the amvber folder

cp ../../scripts/amber/* .

# make the calculation folder
mkdir fep

# make the main files
for i in free vacuum; do
mkdir fep/$i
# make the sub files
for j in discharge input vanish; do
mkdir fep/$i/$j
# end subfile creation
done
# end main file creation
done 

# now we need to set up the input file folder
# first move the top and crd files

# free files first
cp equilibration/system.rst7 fep/free/input/SYSTEM.crd
cp equilibration/MOL_sol.prm7 fep/free/input/SYSTEM.top
cp free_sim.cfg fep/free/input/sim.cfg


# now the vacuum files
cp gas/MOL.rst7 fep/vacuum/input/SYSTEM.crd
cp gas/MOL.prm7 fep/vacuum/input/SYSTEM.top
cp vacuum_sim.cfg fep/vacuum/input/sim.cfg
# move the morph script to gas phase
mv morph_step*.py fep/vacuum/input

cd fep/vacuum/input
# run the morph scripts
python morph_step1.py
python morph_step2.py
# now cp the morph files to the free calc
cp MORPH.pert.discharge ../../free/input
cp MORPH.pert.vanish ../../free/input

# we can move the vanish/discharge files 
cd ../../../
cp fep_discharge.sh fep/free/discharge/
cp fep_discharge.sh fep/vacuum/discharge

cp fep_vanish.sh fep/free/vanish/
cp fep_vanish.sh fep/vacuum/vanish 

# now we want to make a master job control file
cp run_fep.sh fep/

