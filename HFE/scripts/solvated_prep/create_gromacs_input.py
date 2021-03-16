"Use biosimspace to write gromacs input files."
import BioSimSpace as BSS

system = BSS.IO.readMolecules(["MOL.rst7", "MOL.prm7"])
BSS.IO.saveMolecules("MOL", system, ["gro87", "grotop"])
