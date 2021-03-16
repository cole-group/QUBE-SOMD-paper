"""
Convert the pdb file to sdf and retype the molecule transfer the torsion parameters from the new molecule to the old and serialise again.
"""
from QUBEKit.ligand import Ligand
from QUBEKit.parametrisation import OpenFF, XML
from openff.toolkit.topology import Molecule as OFFMolecule
import parmed
from simtk.openmm import app, XmlSerializer

# convert to sdf first
off_mol = OFFMolecule.from_file("MOL.pdb")
off_mol.to_file("MOL.sdf", "sdf")

# make the qube mol
qb_mol1 = Ligand.from_file("MOL.sdf")
# apply params
OpenFF(qb_mol1)
# write out the params
qb_mol1.write_parameters("temp")

# now load a second molecule
qb_mol2 = Ligand.from_file("MOL.sdf")
# apply params from xml
XML(qb_mol2, input_file="MOL.xml")

# transfer the torsions
qb_mol2.PeriodicTorsionForce = qb_mol1.PeriodicTorsionForce
qb_mol2.combination = "opls"
qb_mol2.write_parameters(name="new")

# now build the parmed system
pdb = app.PDBFile("MOL.pdb")
ff = app.ForceField("new.xml")
system = ff.createSystem(pdb.topology)
p_sys = parmed.openmm.load_topology(pdb.topology, system, xyz=pdb.positions)
p_sys.save("test.prmtop")
