# a util scrip to prep the gromacs file inputs
# this should be ran from the gromacs folder

cp ../../scripts/gromacs/* .

# make the calculation folder
mkdir fep
# now for each of the windows prep a folder and write the run files
for i in $(seq 0 1 19); do
# make the window folder
mkdir fep/lambda$i

# cp the relaxed input files
cp equilibration/MOL_sol.grotop fep/lambda$i/system.top
cp equilibration/system.gro87 fep/lambda$i/system.gro

# now prep all files
for filename in minimize equil_npt equil_npt2 equil_nvt prod; do
cp $filename.X.mdp $filename.$i.mdp
# now sed
sed -i "" "s/XXX/$i/g" $filename.$i.mdp
# now mv over 
mv $filename.$i.mdp fep/lambda$i
done

# now make the main control file for the window
echo "#FEP JOB SCRIPT" >> fep_window$i.sh
echo "echo 'STARTING MINIMIZATION'" >> fep_window$i.sh
echo "gmx grompp -f minimize.$i.mdp -p system.top -c system.gro -o minimize.$i.tpr" >> fep_window$i.sh
echo "gmx mdrun -nt 4 -v -deffnm minimize.$i" >> fep_window$i.sh
echo "echo 'STARTING NVT EQUILIBRATION'" >> fep_window$i.sh
echo "gmx grompp -f equil_nvt.$i.mdp -c minimize.$i.gro -p system.top -o nvt_equil.$i.tpr" >> fep_window$i.sh
echo "gmx mdrun -nt 4 -v -deffnm nvt_equil.$i" >> fep_window$i.sh
echo "echo 'STARTING NPT RUN 1'" >> fep_window$i.sh
echo "gmx grompp -f equil_npt.$i.mdp -c nvt_equil.$i.gro -t nvt_equil.$i.cpt -p system.top -o npt_equil1.$i.tpr" >> fep_window$i.sh
echo "gmx mdrun -nt 4 -v -deffnm npt_equil1.$i" >> fep_window$i.sh
echo "echo 'STARTING NPT RUN 2'" >> fep_window$i.sh
echo "gmx grompp -f equil_npt2.$i.mdp -c npt_equil1.$i.gro -t npt_equil1.$i.cpt -p system.top -o npt_equil2.$i.tpr" >> fep_window$i.sh
echo "gmx mdrun -nt 4 -v -deffnm npt_equil2.$i" >> fep_window$i.sh
echo "echo 'STARTING PRODUCTION'" >> fep_window$i.sh
echo "gmx grompp -f prod.$i.mdp -c npt_equil2.$i.gro -t npt_equil2.$i.cpt -p system.top -o prod.$i.tpr" >> fep_window$i.sh
echo "gmx mdrun -nt 4 -v -deffnm prod.$i" >> fep_window$i.sh
echo "echo 'Done!'" >> fep_window$i.sh

# move the file into the window
mv fep_window$i.sh fep/lambda$i
done

# now write out a master script which will run all of the fep windows
echo "echo 'STARTING FEP RUN'" >> fep_run.sh
echo 'for i in $(seq 0 1 19); do' >> fep_run.sh
echo "echo 'Staring window'" >> fep_run.sh
echo 'cd lambda$i' >> fep_run.sh
echo 'bash fep_window$i.sh' >> fep_run.sh
echo 'cd ../' >> fep_run.sh
echo 'done' >> fep_run.sh
# put the master file in the fep dir 
mv fep_run.sh fep/
