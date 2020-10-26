Here are all the files (scripts and results) with the single point energies of proteins and other molecules.
These are used to validate the correct implementation of new methods/features in SOMD and the correct parsing of input files.
- The "protein-fragments" folder contains files with the single point energies of proteins, and the protein fragments.
- The "molecules" folder contains input files of small molecules to calculate single point energies.
- test_combRules.py: Compares Sire and SOMD SPEs using both the geometric and arithmetic combining rules 
- fragment_energetics.py, protein_energetics.py: Computes the OpenMM energy of the protein fragments and the protein
- sire_nrg_calc.py: Computes single point energies with Sire. Usage: ~/sire.app/bin/ipython sire_nrg_calc.py -x MOL.xml -p MOL.pdb


A combination of OpenMM and ParmEd is used to analyse the files. This is done by creating an OpenMM system in the usual way for the first fragment, then each successive fragment is loaded into this same system. Since QUBEKit uses geometric combination rules for the Lennard-Jones parameters (rather than Amber’s arithmetic combination rules), a small change is made to the OpenMM system object to allow for this. Having built the system with OpenMM, the energetics are calculated for the whole system, then using ParmEd, a full breakdown of each contributor is obtained (energy from bonds, angles, torsions and non-bonded). These energy results can then be compared with the energies computed using our SOMD interface as a check.
