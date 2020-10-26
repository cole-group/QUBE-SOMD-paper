To run FEP calculations with the QuBe forcefield, fisrtly we need to convert the protein files (xml/pdb) to amber type format(prm7/rst7).
We use [qube_to_prmRst.py](https://github.com/cole-group/qube_project/blob/master/QuBe-SOMD_paper/FEP_preparation/qube_to_prmRst.py) to read the xml/pdb files of the protein fragments and generate the corresponding amber files:
Usage: ~/sire.app/bin/ipython qube_to_prmRst.py -x fragX.xml -p fragX.pdb

Then we combine the fragments to get the whole protein using [combine.py](https://github.com/cole-group/qube_project/blob/master/QuBe-SOMD_paper/FEP_preparation/combine.py).
Usage: run ./combine.py --system1 fragX.prm7 fragX.rst7 --system2 fragY.prm7 fragY.rst7 --output fragX_Y

Folders "Group1" and "Group2" contain the qube files (pdb/xml) of the protein fragments that can be used to generate the amber files for the proteins for the 
respective set of ligands. 
