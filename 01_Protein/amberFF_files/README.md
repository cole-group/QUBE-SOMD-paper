Both proteins of groups 1 and 2 consist of four fragments. 
We need to parameterise each fragment separately and then combine them to form the proteins. 
To parametrise the protein fragments with the ff14SB forcefield, we are using the [parameterise.py](https://github.com/michellab/BioSimSpace/blob/devel/nodes/playground/parameterise.py) script of BioSimSpace. 
The script takes as an input the protein pdb file and returns the topology and coordinates amber files (prm7/rst7)

- usage: run ./parameterise.py --input FILE.pdb --forcefield ff14SB --output amber_FILE

To combine the newly parameterised fragments we use the [combine.py](https://github.com/michellab/BioSimSpace/blob/devel/nodes/playground/combine.py) script of BioSimSpace, that reads the amber files (2 moleculear systems at a time) 
and combines them to form a complex. Then the coordinates and topology files of the complex are returned. 
This process is done three times, in order to gradually construct the whole protein. 

- usage: run ./combine.py --system1 amber_FILEx.prm7 amber_FILEx.rst7 --system2 amber_FILEy.prm7 amber_FILEy.rst7 --output amber_FILEx_y

