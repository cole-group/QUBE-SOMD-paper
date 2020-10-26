### Creating the files

QUBEKit-pro is an extension to [QUBEKit](https://github.com/qubekit/QUBEKit) for analysing proteins.
It comes bundled with the QUBEKit install and is accessed with the `QUBEKit-pro` command.

Practically, creating the necessary files via QUBEKit is done simply through the command line interface. 
Navigating to a directory containing the ONETEP output file and the relevant pdb files, the command 
`QUBEKit-pro -build <name of pdb>` will perform the necessary steps described below.

The requisite non-bonded parameter information is extracted from the ONETEP file. 
QUBEKit-pro searches the ddec.onetep file top to bottom until all charge data is populated.
As such, if using QUBEKit-pro for a fragment, 
only the charge data for that fragment should be left in the ddec file (note that this process is not currently automated and users are responsible for generating and preparing the required ONETEP files).
For each atom, this is the partial charge and the effective volume. 
These parameters are then used to calculate the Lennard-Jones parameters via atom-in-molecule electron density partitioning. 
Following this calculation, the atoms previously marked for symmetrisation are symmetrised. 
These parameters are then stored for use in producing the final xml file.

The files needed to reproduce the data in this paper are provided in the directory above.
For example, in group1/frag2/ the input files are `frag2.pdb` and the (truncated) ONETEP output file corresponding to the atoms in fragment 2. The output files (`QUBE_pro_frag2.pdb` and `QUBE_pro_frag2.xml`) are generated using the command `QUBEKit-pro -build frag2.pdb`.


Note that since each atom in the protein is in a unique environment, 
and therefore has unique charge and Lennard-Jones parameters, 
each atom in QUBEKit-pro is assigned a unique type. 
This means the entire protein is to be treated as a single molecule rather than a collection of residues with fixed atom types.


