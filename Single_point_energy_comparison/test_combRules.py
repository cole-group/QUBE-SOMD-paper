#
# Sire script to compare Sire and SOMD single point energies
#

import os
import re
import sys

from nose.tools import assert_almost_equal
from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.Qt import *
from Sire.ID import *
from Sire.Config import *
from Sire.Analysis import *
from Sire.Tools.DCDFile import *
from Sire.Tools import Parameter, resolveParameters
import Sire.Stream
import time
import numpy as np

####################################################################################################
#
#   Config file parameters
#
####################################################################################################
reference = Parameter("reference", 0, """Where to run the simulation.""")

rf_dielectric = Parameter("reaction field dielectric", 78.3,
                          """Dielectric constant to use if the reaction field cutoff method is used.""")

temperature = Parameter("temperature", 25 * celsius, """Simulation temperature""")

pressure = Parameter("pressure", 1 * atm, """Simulation pressure""")

topfile = Parameter("topfile", "SYSTEM.top",
                    """File name of the topology file containing the system to be simulated.""")

crdfile = Parameter("crdfile", "SYSTEM.crd",
                    """File name of the coordinate file containing the coordinates of the
                       system to be simulated.""")

s3file = Parameter("s3file", "SYSTEM.s3",
                   """Filename for the system state file. The system state after topology and and coordinates
                   were loaded are saved in this file.""")

restart_file = Parameter("restart file", "sim_restart.s3",
                         """Filename of the restart file to use to save progress during the simulation.""")

dcd_root = Parameter("dcd root", "traj", """Root of the filename of the output DCD trajectory files.""")

nmoves = Parameter("nmoves", 1000, """Number of Molecular Dynamics moves to be performed during the simulation.""")

random_seed = Parameter("random seed", None, """Random number seed. Set this if you
                         want to have reproducible simulations.""")

ncycles = Parameter("ncycles", 1,
                    """The number of MD cycles. The total elapsed time will be nmoves*ncycles*timestep""")

ncycles_per_snap = Parameter("ncycles_per_snap", 1, """Number of cycles between saving snapshots""")

save_coords = Parameter("save coordinates", True, """Whether or not to save coordinates.""")

buffered_coords_freq = Parameter("buffered coordinates frequency", 1,
                                 """The number of time steps between saving of coordinates during
                                 a cycle of MD. 0 disables buffering.""")
minimal_coordinate_saving = Parameter("minimal coordinate saving", False, "Reduce the number of coordiantes writing for states"
                                                                    "with lambda in ]0,1[")

time_to_skip = Parameter("time to skip", 0 * picosecond, """Time to skip in picoseconds""")

minimise = Parameter("minimise", False, """Whether or not to perform minimization before the simulation.""")

minimise_tol = Parameter("minimise tolerance", 1, """Tolerance used to know when minimization is complete.""")

minimise_max_iter = Parameter("minimise maximum iterations", 1000, """Maximum number of iterations for minimization.""")

equilibrate = Parameter("equilibrate", False , """Whether or not to perform equilibration before dynamics.""")

equil_iterations = Parameter("equilibration iterations", 2000, """Number of equilibration steps to perform.""")

equil_timestep = Parameter("equilibration timestep", 0.5 * femtosecond, """Timestep to use during equilibration.""")

timestep = Parameter("timestep", 2 * femtosecond, """Timestep for the dynamics simulation.""")

platform = Parameter("platform", "Reference", """Which OpenMM platform should be used to perform the dynamics.""")

precision = Parameter("precision", "mixed", """The floating point precision model to use during dynamics.""")

constraint = Parameter("constraint", "none", """The constraint model to use during dynamics.""")

#cutoff_type = Parameter("cutoff type", "cutoffperiodic", """The cutoff method to use during the simulation.""")
cutoff_type = Parameter("cutoff type", "nocutoff", """The cutoff method to use during the simulation.""")

cutoff_dist = Parameter("cutoff distance", 10 * angstrom,
                        """The cutoff distance to use for the non-bonded interactions.""")

integrator_type = Parameter("integrator", "leapfrogverlet", """The integrator to use for dynamics.""")

inverse_friction = Parameter("inverse friction", 0.1 * picosecond,
                             """Inverse friction time for the Langevin thermostat.""")

andersen = Parameter("thermostat", False,
                     """Whether or not to use the Andersen thermostat (needed for NVT or NPT simulation).""")

barostat = Parameter("barostat", False, """Whether or not to use a barostat (needed for NPT simulation).""")

andersen_frequency = Parameter("andersen frequency", 10.0, """Collision frequency in units of (1/ps)""")

barostat_frequency = Parameter("barostat frequency", 25,
                               """Number of steps before attempting box changes if using the barostat.""")

lj_dispersion = Parameter("lj dispersion", False, """Whether or not to calculate and include the LJ dispersion term.""")

cmm_removal = Parameter("center of mass frequency", 10,
                        "Frequency of which the system center of mass motion is removed.""")

center_solute = Parameter("center solute", False,
                          """Whether or not to centre the centre of geometry of the solute in the box.""")

use_restraints = Parameter("use restraints", False, """Whether or not to use harmonic restaints on the solute atoms.""")

k_restraint = Parameter("restraint force constant", 100.0, """Force constant to use for the harmonic restaints.""")

heavy_mass_restraint = Parameter("heavy mass restraint", 1.10,
                                 """Only restrain solute atoms whose mass is greater than this value.""")

unrestrained_residues = Parameter("unrestrained residues", ["WAT", "HOH"],
                                  """Names of residues that are never restrained.""")

freeze_residues = Parameter("freeze residues", False, """Whether or not to freeze certain residues.""")

frozen_residues = Parameter("frozen residues", ["LGR", "SIT", "NEG", "POS"],
                            """List of residues to freeze if 'freeze residues' is True.""")


use_distance_restraints = Parameter("use distance restraints", False,
                                    """Whether or not to use restraints distances between pairs of atoms.""")

distance_restraints_dict = Parameter("distance restraints dictionary", {},
                                     """Dictionary of pair of atoms whose distance is restrained, and restraint
                                     parameters. Syntax is {(atom0,atom1):(reql, kl, Dl)} where atom0, atom1 are atomic
                                     indices. reql the equilibrium distance. Kl the force constant of the restraint.
                                     D the flat bottom radius. WARNING: PBC distance checks not implemented, avoid
                                     restraining pair of atoms that may diffuse out of the box.""")

hydrogen_mass_repartitioning_factor = Parameter("hydrogen mass repartitioning factor",None,
                                     """If not None and is a number, all hydrogen atoms in the molecule will
                                        have their mass increased by the input factor. The atomic mass of the heavy atom
                                        bonded to the hydrogen is decreased to keep the mass constant.""")

## Free energy specific keywords

lambda_val = Parameter("lambda_val", 0.0,
                       """Value of the lambda parameter at which to evaluate free energy gradients.""")

delta_lambda = Parameter("delta_lambda", 0.001,
                         """Value of the lambda interval used to evaluate free energy gradients by finite difference.""")

lambda_array = Parameter("lambda array",[] ,
                        """Array with all lambda values lambda_val needs to be part of the array. """)

shift_delta = Parameter("shift delta", 2.0,
                        """Value of the Lennard-Jones soft-core parameter.""")

coulomb_power = Parameter("coulomb power", 0,
                          """Value of the Coulombic soft-core parameter.""")

energy_frequency = Parameter("energy frequency", 1,
                             """The number of time steps between evaluation of free energy gradients.""")

simfile = Parameter("outdata_file", "simfile.dat", """Filename that records all output needed for the free energy analysis""")

perturbed_resnum = Parameter("perturbed residue number",1,"""The residue number of the molecule to morph.""")

verbose = Parameter("verbose", False, """Print debug output""")


def getDummies(molecule):
 
    natoms = molecule.nAtoms()
    atoms = molecule.atoms()

    from_dummies = None
    to_dummies = None

    for x in range(0, natoms):
        atom = atoms[x]
        if atom.property("initial_ambertype") == "du":
            if from_dummies is None:
                from_dummies = molecule.selectAll(atom.index())
            else:
                from_dummies += molecule.selectAll(atom.index())
        elif atom.property("final_ambertype") == "du":
            if to_dummies is None:
                to_dummies = molecule.selectAll(atom.index())
            else:
                to_dummies += molecule.selectAll(atom.index())

    return to_dummies, from_dummies



def createSystemFreeEnergy(molecules, PERT):
    r"""creates the system for free energy calculation
    Parameters
    ----------
    molecules : Sire.molecules
        Sire object that contains a lot of information about molecules
    Returns
    -------
    system : Sire.system
    """
   

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber)[0].molecule()
        moleculeList.append(molecule)

    # Scan input to find a molecule with passed residue number 
    # The residue name of the first residue in this molecule is
    # used to name the solute. This is used later to match
    # templates in the flex/pert files.

    solute = None
    for molecule in moleculeList:
        if ( molecule.residue(ResIdx(0)).number() == ResNum(perturbed_resnum.val) ):
            solute = molecule
            moleculeList.remove(molecule)
            break

    if solute is None:
        print ("FATAL ! Could not find a solute to perturb with residue number %s in the input ! Check the value of your cfg keyword 'perturbed residue number'" % perturbed_resnum.val)
        sys.exit(-1)

    #solute = moleculeList[0]

    lig_name = solute.residue(ResIdx(0)).name().value()

    solute = solute.edit().rename(lig_name).commit()

    perturbations_lib = PerturbationsLibrary(PERT)
    solute = perturbations_lib.applyTemplate(solute)

    perturbations = solute.property("perturbations")

    lam = Symbol("lambda")

    initial = Perturbation.symbols().initial()
    final = Perturbation.symbols().final()

    solute = solute.edit().setProperty("perturbations",
                                       perturbations.recreate((1 - lam) * initial + lam * final)).commit()

    # We put atoms in three groups depending on what happens in the perturbation
    # non dummy to non dummy --> the hard group, uses a normal intermolecular FF
    # non dummy to dummy --> the todummy group, uses SoftFF with alpha = Lambda
    # dummy to non dummy --> the fromdummy group, uses SoftFF with alpha = 1 - Lambda
    # We start assuming all atoms are hard atoms. Then we call getDummies to find which atoms
    # start/end as dummies and update the hard, todummy and fromdummy groups accordingly

    solute_grp_ref = MoleculeGroup("solute_ref", solute)
    solute_grp_ref_hard = MoleculeGroup("solute_ref_hard")
    solute_grp_ref_todummy = MoleculeGroup("solute_ref_todummy")
    solute_grp_ref_fromdummy = MoleculeGroup("solute_ref_fromdummy")

    solute_ref_hard = solute.selectAllAtoms()
    solute_ref_todummy = solute_ref_hard.invert()
    solute_ref_fromdummy = solute_ref_hard.invert()

    to_dummies, from_dummies = getDummies(solute)

    if to_dummies is not None:
        ndummies = to_dummies.count()
        dummies = to_dummies.atoms()

        for x in range(0, ndummies):
            dummy_index = dummies[x].index()
            solute_ref_hard = solute_ref_hard.subtract(solute.select(dummy_index))
            solute_ref_todummy = solute_ref_todummy.add(solute.select(dummy_index))

    if from_dummies is not None:
        ndummies = from_dummies.count()
        dummies = from_dummies.atoms()

        for x in range(0, ndummies):
            dummy_index = dummies[x].index()
            solute_ref_hard = solute_ref_hard.subtract(solute.select(dummy_index))
            solute_ref_fromdummy = solute_ref_fromdummy.add(solute.select(dummy_index))

    solute_grp_ref_hard.add(solute_ref_hard)
    solute_grp_ref_todummy.add(solute_ref_todummy)
    solute_grp_ref_fromdummy.add(solute_ref_fromdummy)

    solutes = MoleculeGroup("solutes")
    solutes.add(solute)

    molecules = MoleculeGroup("molecules")
    molecules.add(solute)

    solvent = MoleculeGroup("solvent")

    #for molecule in moleculeList[1:]:
    for molecule in moleculeList:
        molecules.add(molecule)
        solvent.add(molecule)

    all = MoleculeGroup("all")

    all.add(molecules)
    all.add(solvent)

    all.add(solutes)
    all.add(solute_grp_ref)
    all.add(solute_grp_ref_hard)
    all.add(solute_grp_ref_todummy)
    all.add(solute_grp_ref_fromdummy)

    # Add these groups to the System
    system = System()

    system.add(solutes)
    system.add(solute_grp_ref)
    system.add(solute_grp_ref_hard)
    system.add(solute_grp_ref_todummy)
    system.add(solute_grp_ref_fromdummy)

    system.add(molecules)

    system.add(solvent)

    system.add(all)

    return system


def setupForceFieldsFreeEnergy(system, space, RULE):
    r"""sets up the force field for the free energy calculation
    Parameters
    ----------
    system : Sire.system
    space : Sire.space
    Returns
    -------
    system : Sire.system
    """


    solutes = system[MGName("solutes")]

    solute = system[MGName("solute_ref")]
    solute_hard = system[MGName("solute_ref_hard")]
    solute_todummy = system[MGName("solute_ref_todummy")]
    solute_fromdummy = system[MGName("solute_ref_fromdummy")]

    solvent = system[MGName("solvent")]

    all = system[MGName("all")]

    # ''solvent'' is actually every molecule that isn't perturbed !
    solvent_intraff = InternalFF("solvent_intraff")
    solvent_intraff.add(solvent)

    # Solute bond, angle, dihedral energy
    solute_intraff = InternalFF("solute_intraff")
    solute_intraff.add(solute)

    # Solvent-solvent coulomb/LJ (CLJ) energy
    solventff = InterCLJFF("solvent:solvent")
    if (cutoff_type.val != "nocutoff"):
        solventff.setUseReactionField(True)
        solventff.setReactionFieldDielectric(rf_dielectric.val)
    solventff.add(solvent)

    #Solvent intramolecular CLJ energy
    solvent_intraclj = IntraCLJFF("solvent_intraclj")
    if (cutoff_type.val != "nocutoff"):
        solvent_intraclj.setUseReactionField(True)
        solvent_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solvent_intraclj.add(solvent)

    # Solute intramolecular CLJ energy
    solute_hard_intraclj = IntraCLJFF("solute_hard_intraclj")
    if (cutoff_type.val != "nocutoff"):
        solute_hard_intraclj.setUseReactionField(True)
        solute_hard_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_intraclj.add(solute_hard)

    solute_todummy_intraclj = IntraSoftCLJFF("solute_todummy_intraclj")
    solute_todummy_intraclj.setShiftDelta(shift_delta.val)
    solute_todummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_todummy_intraclj.setUseReactionField(True)
        solute_todummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_todummy_intraclj.add(solute_todummy)

    solute_fromdummy_intraclj = IntraSoftCLJFF("solute_fromdummy_intraclj")
    solute_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    solute_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_fromdummy_intraclj.setUseReactionField(True)
        solute_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_fromdummy_intraclj.add(solute_fromdummy)

    solute_hard_todummy_intraclj = IntraGroupSoftCLJFF("solute_hard:todummy_intraclj")
    solute_hard_todummy_intraclj.setShiftDelta(shift_delta.val)
    solute_hard_todummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_hard_todummy_intraclj.setUseReactionField(True)
        solute_hard_todummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_todummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_todummy_intraclj.add(solute_todummy, MGIdx(1))

    solute_hard_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_hard:fromdummy_intraclj")
    solute_hard_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    solute_hard_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_hard_fromdummy_intraclj.setUseReactionField(True)
        solute_hard_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_fromdummy_intraclj.add(solute_hard, MGIdx(0))
    solute_hard_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))

    solute_todummy_fromdummy_intraclj = IntraGroupSoftCLJFF("solute_todummy:fromdummy_intraclj")
    solute_todummy_fromdummy_intraclj.setShiftDelta(shift_delta.val)
    solute_todummy_fromdummy_intraclj.setCoulombPower(coulomb_power.val)
    if (cutoff_type.val != "nocutoff"):
        solute_todummy_fromdummy_intraclj.setUseReactionField(True)
        solute_todummy_fromdummy_intraclj.setReactionFieldDielectric(rf_dielectric.val)
    solute_todummy_fromdummy_intraclj.add(solute_todummy, MGIdx(0))
    solute_todummy_fromdummy_intraclj.add(solute_fromdummy, MGIdx(1))

    #Solute-solvent CLJ energy
    solute_hard_solventff = InterGroupCLJFF("solute_hard:solvent")
    if (cutoff_type.val != "nocutoff"):
        solute_hard_solventff.setUseReactionField(True)
        solute_hard_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_hard_solventff.add(solute_hard, MGIdx(0))
    solute_hard_solventff.add(solvent, MGIdx(1))

    solute_todummy_solventff = InterGroupSoftCLJFF("solute_todummy:solvent")
    if (cutoff_type.val != "nocutoff"):
        solute_todummy_solventff.setUseReactionField(True)
        solute_todummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_todummy_solventff.add(solute_todummy, MGIdx(0))
    solute_todummy_solventff.add(solvent, MGIdx(1))

    solute_fromdummy_solventff = InterGroupSoftCLJFF("solute_fromdummy:solvent")
    if (cutoff_type.val != "nocutoff"):
        solute_fromdummy_solventff.setUseReactionField(True)
        solute_fromdummy_solventff.setReactionFieldDielectric(rf_dielectric.val)
    solute_fromdummy_solventff.add(solute_fromdummy, MGIdx(0))
    solute_fromdummy_solventff.add(solvent, MGIdx(1))


    # TOTAL
    forcefields = [solute_intraff,
                   solute_hard_intraclj, solute_todummy_intraclj, solute_fromdummy_intraclj,
                   solute_hard_todummy_intraclj, solute_hard_fromdummy_intraclj,
                   solute_todummy_fromdummy_intraclj,
                   solvent_intraff,
                   solventff, solvent_intraclj,
                   solute_hard_solventff, solute_todummy_solventff, solute_fromdummy_solventff]


    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)

    if (cutoff_type.val != "nocutoff"):
        system.setProperty("switchingFunction", CHARMMSwitchingFunction(cutoff_dist.val))
    else:
        system.setProperty("switchingFunction", NoCutoff())

    system.setProperty("combiningRules", VariantProperty(RULE))
    system.setProperty("coulombPower", VariantProperty(coulomb_power.val))
    system.setProperty("shiftDelta", VariantProperty(shift_delta.val))

    # TOTAL
    total_nrg = solute_intraff.components().total() + solute_hard_intraclj.components().total() + \
                solute_todummy_intraclj.components().total(0) + solute_fromdummy_intraclj.components().total(0) + \
                solute_hard_todummy_intraclj.components().total(
                    0) + solute_hard_fromdummy_intraclj.components().total(0) + \
                solute_todummy_fromdummy_intraclj.components().total(0) + \
                solvent_intraff.components().total() + solventff.components().total() + \
                solvent_intraclj.components().total() + \
                solute_hard_solventff.components().total() + \
                solute_todummy_solventff.components().total(0) + \
                solute_fromdummy_solventff.components().total(0)

    e_total = system.totalComponent()

    lam = Symbol("lambda")

    system.setComponent(e_total, total_nrg)

    system.setConstant(lam, 0.0)

    system.add(PerturbationConstraint(solutes))

    # NON BONDED Alpha constraints for the soft force fields

    system.add(PropertyConstraint("alpha0", FFName("solute_todummy_intraclj"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_fromdummy_intraclj"), 1 - lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_hard:todummy_intraclj"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_hard:fromdummy_intraclj"), 1 - lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_todummy:fromdummy_intraclj"), Max(lam, 1 - lam)))
    system.add(PropertyConstraint("alpha0", FFName("solute_todummy:solvent"), lam))
    system.add(PropertyConstraint("alpha0", FFName("solute_fromdummy:solvent"), 1 - lam))

    system.setComponent(lam, lambda_val.val)

    # printEnergies( system.componentValues() )

    return system


def setupMovesFreeEnergy(system, random_seed, REFERENCE, lam_val, RULE):



    molecules = system[MGName("molecules")]
    solute = system[MGName("solute_ref")]
    solute_hard = system[MGName("solute_ref_hard")]
    solute_todummy = system[MGName("solute_ref_todummy")]
    solute_fromdummy = system[MGName("solute_ref_fromdummy")]

    Integrator_OpenMM = OpenMMFrEnergyST(molecules, solute, solute_hard, solute_todummy, solute_fromdummy)
    Integrator_OpenMM.setRandomSeed(random_seed)
    Integrator_OpenMM.setIntegrator(integrator_type.val)
    Integrator_OpenMM.setFriction(inverse_friction.val)  # Only meaningful for Langevin/Brownian integrators
    Integrator_OpenMM.setPlatform(platform.val)
    Integrator_OpenMM.setCombiningRules(RULE)
    Integrator_OpenMM.setConstraintType(constraint.val)
    Integrator_OpenMM.setCutoffType(cutoff_type.val)
    Integrator_OpenMM.setFieldDielectric(rf_dielectric.val)
    Integrator_OpenMM.setAlchemicalValue(lambda_val.val)
    Integrator_OpenMM.setAlchemicalArray(lambda_array.val)
    Integrator_OpenMM.setDeviceIndex(str(REFERENCE))
    Integrator_OpenMM.setCoulombPower(coulomb_power.val)
    Integrator_OpenMM.setShiftDelta(shift_delta.val)
    Integrator_OpenMM.setDeltatAlchemical(delta_lambda.val)
    Integrator_OpenMM.setPrecision(precision.val)
    Integrator_OpenMM.setTimetoSkip(time_to_skip.val)
    Integrator_OpenMM.setBufferFrequency(buffered_coords_freq.val)

    if cutoff_type != "nocutoff":
        Integrator_OpenMM.setCutoffDistance(cutoff_dist.val)

    Integrator_OpenMM.setCMMremovalFrequency(cmm_removal.val)

    Integrator_OpenMM.setEnergyFrequency(energy_frequency.val)

    if use_restraints.val:
        Integrator_OpenMM.setRestraint(True)

    if andersen.val:
        Integrator_OpenMM.setTemperature(temperature.val)
        Integrator_OpenMM.setAndersen(andersen.val)
        Integrator_OpenMM.setAndersenFrequency(andersen_frequency.val)

    if barostat.val:
        Integrator_OpenMM.setPressure(pressure.val)
        Integrator_OpenMM.setMCBarostat(barostat.val)
        Integrator_OpenMM.setMCBarostatFrequency(barostat_frequency.val)

    #This calls the OpenMMFrEnergyST initialise function
    Integrator_OpenMM.initialise()
    velocity_generator = MaxwellBoltzmann(temperature.val)
    velocity_generator.setGenerator(RanGenerator(random_seed))

    mdmove = MolecularDynamics(molecules, Integrator_OpenMM, timestep.val,
                              {"velocity generator":velocity_generator})

  

    moves = WeightedMoves()
    moves.add(mdmove, 1)

    if (not random_seed):
        random_seed = RanGenerator().randInt(100000, 1000000)

 

    moves.setGenerator(RanGenerator(random_seed))

    return moves
def test_cr(TOP, CRD, PERT, RULE):
    combining_rules = Parameter("combining rules", RULE,
                            """Combining rules to use for the non-bonded interactions.""")

    amber = Amber()
    
    (molecules, space) = amber.readCrdTop(CRD, TOP)
    system = createSystemFreeEnergy(molecules, PERT)
    
    system = setupForceFieldsFreeEnergy(system, space, RULE)
    if random_seed.val:
        ranseed = random_seed.val
    else:
        ranseed = RanGenerator().randInt(100000, 1000000)
    
    

    moves = setupMovesFreeEnergy(system, ranseed, reference.val, lambda_val.val, RULE)

    
    # Get energy from Sire
    nrg_sire = system.energy()
    # Get energy from SOMD
    mdmoves = moves.moves()[0]
    integrator = mdmoves.integrator()
    nrg_somd = integrator.getPotentialEnergy(system)

    if (RULE == 'arithmetic'):
        global diff_arith
        diff_arith  = (nrg_sire - nrg_somd).value()
        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print ("TESTING ARITHMETIC COMBINING RULES")
        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        print ("* The energy from Sire is: %s" %nrg_sire)
        print ("* The energy from SOMD is: %s" % nrg_somd)

        print (" For the arithmetic combining rules the single point energy difference between sire and somd at lambda %s is %s " % (lambda_val.val,diff_arith) )
    elif (RULE == 'geometric'):
        global diff_geom
        diff_geom = (nrg_sire - nrg_somd).value()
        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print ("TESTING GEOMETRIC COMBINING RULES")
        print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

        print ("* The energy from Sire is: %s" %nrg_sire)
        print ("* The energy from SOMD is: %s" % nrg_somd)

        print (" For the geometric combining rules the single point energy difference between sire and somd at lambda %s is %s " % (lambda_val.val,diff_geom) )
if __name__ == '__main__':
    TOP = 'MOL.prm7'
    CRD = 'MOL.rst7'
    PERT = '../MORPH.pert'

    test_cr(TOP, CRD, PERT, 'arithmetic')
    test_cr(TOP, CRD, PERT, 'geometric')
    d= diff_arith - diff_geom
    print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print ("The difference of the single point energies is %s kcal/mol when using geometric and arithmetic combining rules" %d)
    print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    if abs(d) > 0.0001:
        if verbose:
            print("WARNING! GEOMETRIC AND ARITHMETIC COMBINING RULES DO NOT GIVE THE SAME RESULT")
    assert_almost_equal(diff_arith, 0.0, 4)
    assert_almost_equal(diff_geom , 0.0, 4)
