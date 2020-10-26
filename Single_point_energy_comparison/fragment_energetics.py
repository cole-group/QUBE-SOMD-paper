#!/usr/bin/env python3

from simtk.openmm import CustomNonbondedForce, LangevinIntegrator, Platform, XmlSerializer
from simtk.openmm.app import ForceField, NoCutoff, PDBFile, Simulation

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


def calculate_fragment_energetics(frag_no=1):
    """
    * Create an OpenMM system with a fragment.
    * Calculate the energy of the system and print.
    :param frag_no: The number of the fragment being analysed (used to access files).
    """
    os.chdir(f'group2/frag{frag_no}')
    # Necessary due to size of calculation
    sys.setrecursionlimit(15000)

    pdb = PDBFile(f'QUBE_pro_frag{frag_no}.pdb')
    forcefield = ForceField(f'QUBE_pro_frag{frag_no}_plus.xml')

    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=NoCutoff,
    )

    system = apply_opls_combo(system)

    with open(f'QUBE_pro_frag{frag_no}_out.xml', 'w') as outfile:
        serialized_system = XmlSerializer.serialize(system)
        outfile.write(serialized_system)

    # Create the integrator to do Langevin dynamics
    integrator = LangevinIntegrator(
        298.15 * unit.kelvin,  # Temperature of heat bath
        1.0 / unit.picoseconds,  # Friction coefficient
        2.0 * unit.femtoseconds,  # Time step
    )

    platform = Platform.getPlatformByName('CPU')
    simulation = Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)
    print('energy from openmm library')
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())

    structure = parmed.load_file(f'QUBE_pro_frag{frag_no}.pdb')
    energy_comps = parmed.openmm.energy_decomposition_system(structure, system)

    total_energy = 0.0
    for comp in energy_comps:
        total_energy += comp[1]
        print(*comp)

    print(f'Total energy {total_energy: 6.6f}')


if __name__ == '__main__':
    calculate_fragment_energetics()
