Here are the instructions to run automated Hydration Free Energy calculations using the QUBE force field.
The following steps are performed automatically using [build_solvated_systems.sh](https://github.com/cole-group/QUBE-SOMD-paper/blob/master/HFE/scripts/build_solvated_systems.sh) and will produce simulation ready input files for both SOMD and Gromacs.

## File Preparation 
- copy the scripts folder to your desired workplace
- make a new folder for your target molecule in the same place as the scripts folder and add the QUBE pdb and xml file, these should be named `MOL.pdb` and `MOL.xml`
- copy the [build_solvated_systems.sh](https://github.com/cole-group/QUBE-SOMD-paper/blob/master/HFE/scripts/build_solvated_systems.sh) script into this folder and run via `bash build_solvated_systems.sh`

Running `tree` should show a directory structure like this

```
.
├── MOL.pdb
├── MOL.xml
├── amber
│   ├── equilibration
│   │   ├── MOL_sol.prm7
│   │   └── system.rst7
│   ├── gas
│   │   ├── MOL.prm7
│   │   └── MOL.rst7
│   └── solvated
│       ├── MOL_sol.pdb
│       ├── MOL_sol.prm7
│       └── MOL_sol.rst7
├── build_solvated_systems.sh
├── gromacs
│   ├── equilibration
│   │   ├── MOL_sol.grotop
│   │   └── system.gro87
│   ├── gas
│   │   ├── MOL.gro87
│   │   └── MOL.grotop
│   └── solvated
│       ├── MOL_sol.gro87
│       └── MOL_sol.grotop
├── output.yaml
└── system.pdb
```
Now we have Amber and Gromacs style coordinate and system files for the molecule in vacuum and solvated both before (solvated) and after (equilibration) a short equilibration simulation.

The simulation control input files for Amber and Gromacs can then be generated using the `prep_somd.sh` and `prep_gromacs.sh` scripts respectively, to run them first copy them to their respective folder and run as `bash prep_somd(gromacs).sh`.
Running tree again should now give the following structure with the fep directories and control scripts (`run_fep.sh`) which will launch the calculations.

```
.
├── MOL.pdb
├── MOL.sdf
├── MOL.xml
├── amber
│   ├── equilibration
│   │   ├── MOL_sol.prm7
│   │   └── system.rst7
│   ├── fep
│   │   ├── free
│   │   │   ├── discharge
│   │   │   │   └── fep_discharge.sh
│   │   │   ├── input
│   │   │   │   ├── MORPH.pert.discharge
│   │   │   │   ├── MORPH.pert.vanish
│   │   │   │   ├── SYSTEM.crd
│   │   │   │   ├── SYSTEM.top
│   │   │   │   └── sim.cfg
│   │   │   └── vanish
│   │   │       └── fep_vanish.sh
│   │   ├── run_fep.sh
│   │   └── vacuum
│   │       ├── discharge
│   │       │   └── fep_discharge.sh
│   │       ├── input
│   │       │   ├── MORPH.pert.discharge
│   │       │   ├── MORPH.pert.vanish
│   │       │   ├── SYSTEM.crd
│   │       │   ├── SYSTEM.top
│   │       │   ├── morph_step1.py
│   │       │   ├── morph_step2.py
│   │       │   └── sim.cfg
│   │       └── vanish
│   │           └── fep_vanish.sh
│   ├── fep_discharge.sh
│   ├── fep_vanish.sh
│   ├── free_sim.cfg
│   ├── gas
│   │   ├── MOL.prm7
│   │   └── MOL.rst7
│   ├── prep_somd.sh
│   ├── run_fep.sh
│   ├── solvated
│   │   ├── MOL_sol.pdb
│   │   ├── MOL_sol.prm7
│   │   └── MOL_sol.rst7
│   └── vacuum_sim.cfg
├── build_solvated_systems.sh
├── gromacs
│   ├── equil_npt.X.mdp
│   ├── equil_npt2.X.mdp
│   ├── equil_nvt.X.mdp
│   ├── equilibration
│   │   ├── MOL_sol.grotop
│   │   └── system.gro87
│   ├── fep
│   │   ├── fep_run.sh
│   │   ├── lambda0
│   │   │   ├── equil_npt.0.mdp
│   │   │   ├── equil_npt2.0.mdp
│   │   │   ├── equil_nvt.0.mdp
│   │   │   ├── fep_window0.sh
│   │   │   ├── minimize.0.mdp
│   │   │   ├── prod.0.mdp
│   │   │   ├── system.gro
│   │   │   └── system.top

```

##Script details##
Here we explain the steps the scripts follow to prepare the input files in more detail for each phase of the hydration free energy calculations.

- For the vacuum simulations: 
1. Starting with the pdb/xml files of the ligands, generate the corresponding amber files (prm7/rst7) using [qube_to_prmRst.py](https://github.com/cole-group/QUBE-SOMD-paper/blob/master/qube_to_prmRst.py)
3. Generate the MORPH.pert files using the scripts [morph_step1.py](https://github.com/cole-group/QUBE-SOMD-paper/blob/master/HFE/scripts/morph_step1.py) and [morph_step2.py](https://github.com/cole-group/QUBE-SOMD-paper/blob/master/HFE/scripts/morph_step2.py). These scripts read the parameters of the molecules and return the pert files for the discharge and the vanish steps.
4. Create the folder architecture and submit the simulations. 

- For the solvated simulations: 
1. Starting with the pdb/xml files of the ligands, generate the corresponding amber files (prm7/rst7) using [qube_to_prmRst.py](https://github.com/cole-group/QUBE-SOMD-paper/blob/master/qube_to_prmRst.py)
2. Use the BioSimSpace.app python (e.g. ~/biosimspace.app/bin/ipython) to solvate the molecule with [solvate.py](https://github.com/cole-group/QUBE-SOMD-paper/blob/master/solvate.py):
```
run ./solvate.py --input MOL.prm7 MOL.rst7 --output MOL_sol --water tip3p --box_dim 26
```
The dimension of the side of the box is in Angstrom. (This might need to be changed)

3. Equilibrate the system with [amberequilibration.py](https://github.com/cole-group/QUBE-SOMD-paper/blob/master/amberequilibration.py) also with BioSimSpace.
```
run ./amberequilibration.py --input MOL_sol.prm7 MOL_sol.rst7 --output MOL
```
From this point on, the process is the same as the one for the molecules in vacuum:
5. Generate the MORPH.pert files using the scripts [morph_step1.py](https://github.com/cole-group/QUBE-SOMD-paper/blob/master/HFE/scripts/morph_step1.py) and [morph_step2.py](https://github.com/cole-group/QUBE-SOMD-paper/blob/master/HFE/scripts/morph_step2.py).
6. Create the folder architecture and submit the simulations.

Analysis: 

`~/sire.app/bin/analyse_freenrg mbar -i lambda-*/simfile.dat -o out.dat -p 90` generates a dat file with the free energies for each step (discharge, vanish) for both legs. 

Corrections: 
- `FUNC.py`: Evaluates the electrostatic correction for the free energy change: FUNC_corr. This is run for lambda= 0 at the discharge leg.
- `~/sire.app/bin/lj-tailcorrection -C sim.cfg -l <lambda> -b 1.00 -r traj000000001.dcd -s 20` Evaluates the end-point correction for the truncated vdW potentials. This is run for lambda= 0 and lambda= 1 of the vanish leg. 

DG_LJCOR = (LJ correction at lambda 1.0 ) - (LJ correction at lambda 0.0)

To derive the hydration free energy, we use the following formula:

`DDG_hyd = (DG_Vac_Discharge + DG_Vac_Vanish) - (DG_Solv_Discharge + DG_Solv_Vanish) + FUNC_corr - DG_LJCOR`
