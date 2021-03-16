# script to build up the solvated systems
# get all the helper scripts
cp ../scripts/solvated_prep/* .

# first convert the qube input to prmst amber format
python qube_to_prmRst.py -x MOL.xml -p MOL.pdb

# now convert to python
python create_gromacs_input.py

mkdir amber
mkdir gromacs

# make the amber and gromacs files
for state in gas solvated equilibration; do
mkdir amber/$state
mkdir gromacs/$state
done


# now solvate the system this will also produce gromacs files
python solvate.py --input MOL.prm7 MOL.rst7 --output MOL_sol --water tip3p --box_dim 26

# now run the amber equil
python amberequilibration.py --input MOL_sol.prm7 MOL_sol.rst7 --output system


# clean all gromacs files
python clean_gromacs.py


# now move all gromacs files
cp MOL_sol.grotop gromacs/equilibration
mv MOL_sol.gro* gromacs/solvated
mv system.gro87 gromacs/equilibration
mv MOL.gro* gromacs/gas

# now move all amber files
mv MOL.prm7 amber/gas
mv MOL.rst7 amber/gas
cp MOL_sol.prm7 amber/equilibration
mv MOL_sol.* amber/solvated
mv system.rst7 amber/equilibration

# clean up the python files
rm *.py 

echo "Done!"

