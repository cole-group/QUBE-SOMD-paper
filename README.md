## BioSimSpace protocols and scripts used to generate the input files for FEP calculations. 

These data/tutorials accompany the paper: [Nelson L, Bariami S, Ringrose C, Horton JT, Kurdekar V, Mey ASJS, Michel J, Cole DJ. Implementation of the QUBE force field in SOMD for high-throughput alchemical free energy calculations, Journal of Chemical Information and Modeling, 2021, 61, 2124-2130.](https://doi.org/10.1021/acs.jcim.1c00328)

The following details the procedure for protein-ligand binding free energy calculations. A separate tutorial on hydration free energy calculations is provided in the directory HFE/.

## File Setup for Free Energy Calculations

### 01) Protein setup: 

  *For use with **AMBER**:*
  - Given the four fragments of the protein in pdb format, we are going to use the [parameterise.py](https://github.com/michellab/BioSimSpace/blob/devel/nodes/playground/parameterise.py) script of BioSimSpace to parameterise them with the ff14SB amber forcefield.
  - Install BioSimSpace (https://github.com/michellab/BioSimSpace/) and from now on use the BioSimSpace python (~/biosimspace.app/bin/ipython) to run the following commands: 
  - Produce AMBER parameterised .rst7 and .prm7 files for each protein fragment: ```run ./parameterise.py --input FILE.pdb --forcefield ff14SB --output FILE```
  - Combine the protein fragments to get the whole protein: ```run ./combine.py --system1 FILE_1.prm7 FILE_1.rst7 --system2 FILE_2.prm7 FILE_2.rst7 --output FILE_12``` (The combining process is done three times to combine the four fragments together.)
  
  *For use with **QUBE**:*
  - Given the four fragments of the protein in xml and pdb format, we are going to use the [qube_to_prmRst.py](https://github.com/cole-group/QUBE-SOMD-paper/blob/master/qube_to_prmRst.py) to read the xml/pdb files and generate the corresponding amber files for each fragment:
  - Use the Sire python to run the following command: ```~/sire.app/bin/ipython qube_to_prmRst.py -x fragX.xml -p fragX.pdb```
  - This will produce AMBER parameterised .rst7 and .prm7 files for each protein fragment.
  - The combining process to get the whole protein is the same as the one we did for the amber parameterisation:
  - With ~/biosimspace.app/bin/ipython run: ```run ./combine.py --system1 FILE_1.prm7 FILE_1.rst7 --system2 FILE_2.prm7 FILE_2.rst7 --output FILE_12``` (The combining process is done three times to combine the four fragments together)
  
  
### 02) Ligand setup:

  *For use with **AMBER**:*

  - Make a "Ligands" folder in which you will have the pdb files for each ligand and the parameterise.py script.
  - Use the BioSimSpace.app python (e.g. ~/biosimspace.app/bin/ipython) to run the following command: ```run ./parameterise.py --input FILE.pdb --forcefield gaff2 --output FILE```
  - This will produce AMBER parameterised .rst7 and .prm7 files for each ligand.

  *For use with **QUBE**:*

  - Make a "Ligands" folder into which you will copy the QUBE parameterised pdb and xml files for each ligand and the qube_to_prmRst.py script. An example can be found in the 02_Ligand_Setup directory, and a full set of ligand inputs can be found in the "QUBE-FFfiles" folder above.
  - Use the Sire.app python (e.g. ~/sire.app/bin/ipython) to run the following command: ```run ./qube_to_prmRst.py -p FILE.pdb -x FILE.xml```
  - This will produce .rst7 and .prm7 files for each ligand.
  

**From here the setup is the same for either force field.**
  
### 03) Solvate the ligands:

  - In the "Ligands" folder you should now copy in the solvate.py script.
  - Use the BioSimSpace.app python (e.g. ~/biosimspace.app/bin/ipython) to run the following command: ```run ./solvate.py --input FILE.prm7 FILE.rst7 --output FILE_sol --water tip3p --extent 26```
  - This will create solvated, unbound ligand files (e.g. FILE_sol.prm7). The above is the box size used for our system.
  
### 04) Combine the ligands and protein:

  - Make a "Complex" folder in which you will need the (unsolvated) ligand and protein prm7 and rst7 files, and the combine.py script.
  - Nativage to the ipython in BioSimSpace.app (e.g. biosimspace.app/bin/ipython)
  - Run the following command in ipython: ```run ./combine.py --system1 LigFILE.prm7 LigFILE.rst7 --system2 PROTEIN.prm7 PROTEIN.rst7 --output PROT_LigFILE```
  - This will create unsolvated prm7 and rst7 files for the ligand in complex with the protein.
  
### 05) Solvate the complex:

  - Still within the "Complex" folder you will now also need the solvate.py script.
  - Nativage to the ipython in BioSimSpace.app (e.g. biosimspace.app/bin/ipython)
  - Run the following command in ipython: ```run ./solvate.py --input PROT_LigFILE.prm7 PROT_LigFILE.rst7 --output Complex_sol --water tip3p --box_dim 88```
  - This will create solvated prm7 and rst7 files of the ligand in complex with the protein. The above is the box size used for our system.
  
### 06) Equilibrate the solvated systems:

  - Both the solvated complexes and ligands need to equilibrated, this will be conducted in both the "Complex" and "Ligands" folders respectively. You will now also need the amberequilibration.py script in both of these folders.
  - Nativage to the ipython in BioSimSpace.app (e.g. biosimspace.app/bin/ipython)
  - Run the following command for both bound and unbound ligand systems: ```run ./amberequilibration.py --input FILE_sol.prm7 FILE_sol.rst7 --output FILE_sol_eq```
  - This will create equilibrated rst7 files for the bound and unbound systems (e.g. Lig_sol_eq.rst7 or Complex_sol_eq.rst7)
  
### 07) Generate the files for free energy calculations:

  - Make a "Perturbations" folder in which you will need the final prm7 and equilibrated rst7 files for both environments (e.g. ligands in complex with the protein, and ligands unbound in solution). You will also need the prepareFEP.py script.
  - Nativage to the ipython in BioSimSpace.app (e.g. biosimspace.app/bin/ipython)
  - Run the following command for both bound and unbound environments: ```run ./prepareFEP.py --input1 FILE1_sol.prm7 FILE1_sol_eq.rst7 --input2 FILE2_sol.prm7 FILE2_sol_eq.rst7 --output FILE1_to_FILE2```
  - This will create the perturbation files for your free energy calculations (e.g. for a transition of Lig1 to Lig2 the above command will create .mapping, .mergeat0.pdb. .pert, .prm7 and .rst7 files for this pertubration, both bound and unbound). 

### 08) Free Energy Calculations

1) Create the folder setup:
  - You will need a main folder from which to run the free energy scripts, ```ligand_lambdarun-comb.sh``` and ```complex_lambdarun-comb.sh```, the "pertlist", and where you will also have:
    - A "Perturbations" folder in which will be the pertubations you wish to run, followed by bound and unbound simulation folders containing the relevant files (*i.e.* Lig1_to_Lig2/bound/L1_to_L2_bound.* ) 
    - A "Parameters" folder where you will have the lambda.cfg file

2) The "pertlist"
  - This list details the perturbations you wish to simulate and should include only numbers (*e.g.* 1-2, 2-1 not, Lig1-Lig2 etc.)

3) The lambda.cfg file
  - Configuration files for both our AMBER and QUBE runs can be found in the Parameters/ folder. 
  - There are various parameters which can be altered in these files, namely the number of moves and cycles, the timestep, the type of constraints, the lambda windows used and the platform on which to run the calculation. 

4) The ```ligand_lambdarun-comb.sh``` and ```complex_lambdarun-comb.sh``` scripts
- Script ```ligand_lambdarun-comb.sh``` runs the command for the unbound perturbations, whilst ```complex_lambdarun-comb.sh``` runs the bound perturbations. 
- You will need to specify in both scripts where to read the "somd-freenrg" file from, but if you have used the file structure suggested you should not need to change anything else in these scripts.

5) Once the above is ready, you can start you free energy calculations by running: ```./ligand_lambdarun-comb.sh``` and ```./complex_lambdarun-comb.sh```


## Results

In addition to the above instructions for generating the input files for this paper, the Results/ directory also contains the raw data (free energy calculations and single point energy validations of our SOMD implementation.


