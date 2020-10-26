#!/usr/bin/env python3

from simtk.openmm import CustomNonbondedForce, LangevinIntegrator, Platform
from simtk.openmm.app import ForceField, Modeller, PDBFile, Simulation, NoCutoff

import parmed
from parmed import unit

import os
import sys

import numpy as np


def apply_opls_combo(system, switching_distance=None):
    """
    Apply the OPLS combination rules (geometric) to an OpenMM system.
    :param system: OpenMM system
    :param switching_distance: Distance at which the switching function begins to reduce the interaction.
    :return: New, altered OpenMM system.
    """

    # Get the system information from the openmm system
    forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in
              range(system.getNumForces())}
    # Use the nondonded_force to get the same rules
    nonbonded_force = forces['NonbondedForce']
    lorentz = CustomNonbondedForce(
        'epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)*4.0')
    lorentz.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setUseLongRangeCorrection(True)
    if switching_distance is not None:
        lorentz.setUseSwitchingFunction(True)
        lorentz.setSwitchingDistance(switching_distance)
    system.addForce(lorentz)

    l_j_set = {}
    # For each particle, calculate the combination list again
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        l_j_set[index] = (sigma, epsilon, charge)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(index, charge, 0, 0)

    for i in range(nonbonded_force.getNumExceptions()):
        p1, p2, q, sig, eps = nonbonded_force.getExceptionParameters(i)
        # ALL THE 12, 13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED FORCE
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            charge = 0.5 * (l_j_set[p1][2] * l_j_set[p2][2])
            sig14 = np.sqrt(l_j_set[p1][0] * l_j_set[p2][0])
            nonbonded_force.setExceptionParameters(i, p1, p2, charge, sig14, eps)

    return system


def calculate_protein_energetics():
    """
    * Create an OpenMM system using the first fragment.
    * Add each fragment into the system.
    * Calculate the energy of the system and print.
    """

    os.chdir('group2')
    # Necessary due to size of calculation
    sys.setrecursionlimit(15000)

    frag1 = PDBFile('frag1/no_QUP_frag1.pdb')
    forcefield = ForceField(
        'frag1/QUBE_pro_frag1.xml',
        'frag2/QUBE_pro_frag2_plus.xml',
        'frag3/QUBE_pro_frag3_plus.xml',
        'frag4/QUBE_pro_frag4_plus.xml',
    )

    modeller = Modeller(frag1.topology, frag1.positions)

    frag2 = PDBFile('frag2/no_QUP_frag2.pdb')
    modeller.add(frag2.topology, frag2.positions)

    frag3 = PDBFile('frag3/no_QUP_frag3.pdb')
    modeller.add(frag3.topology, frag3.positions)

    frag4 = PDBFile('frag4/no_QUP_frag4.pdb')
    modeller.add(frag4.topology, frag4.positions)

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
    )

    system = apply_opls_combo(system)

    integrator = LangevinIntegrator(
        298.15 * unit.kelvin,        # Temperature of heat bath
        1.0 / unit.picoseconds,      # Friction coefficient
        2.0 * unit.femtoseconds,     # Time step
    )

    platform = Platform.getPlatformByName('CPU')
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    print('energy from openmm library')
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

    positions = simulation.context.getState(getPositions=True).getPositions()

    with open('output.pdb', 'w') as out_file:
        PDBFile.writeFile(simulation.topology, positions, out_file)

    structure = parmed.load_file('output.pdb')
    energy_comps = parmed.openmm.energy_decomposition_system(structure, system)

    total_energy = 0.0
    for comp in energy_comps:
        total_energy += comp[1]
        print(*comp)

    print(f'Total energy {total_energy: 6.6f}')


if __name__ == '__main__':
    calculate_protein_energetics()
