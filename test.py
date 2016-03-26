from ase.lattice.surface import fcc100, add_adsorbate
from ase.constraints import FixAtoms, FixedPlane
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton
from ase.optimize import BFGS
from ase.optimize.sciopt import SciPyFminBFGS, SciPyFminCG
from ase.optimize.bfgslinesearch import BFGSLineSearch
from ase.calculators.emt import EMT
from ase.structure import molecule
from ase import Atoms
from ase import Atom
from ase.units import *

import sys

from ElectronicMinimize import ElectronicMinimize, ElectronicMinimizeGPU
import numpy as np
from mpi4py import MPI

# 2x2-Al(001) surface with 3 layers and an
# Au atom adsorbed in a hollow site:
def slabCalculation():
    slab = fcc100('Al', size=(2, 2, 3))
    add_adsorbate(slab, 'Au', 1.7, 'hollow')
    slab.center(axis=2, vacuum=4.0)

    # Make sure the structure is correct:
    #from ase.visualize import view
    #view(slab)

    # Fix second and third layers:
    mask = [atom.tag > 1 for atom in slab]
    #print mask
    fixlayers = FixAtoms(mask=mask)

    # Constrain the last atom (Au atom) to move only in the yz-plane:
    plane = FixedPlane(-1, (1, 0, 0))

    slab.set_constraint([fixlayers, plane])

    print(slab.cell / Bohr)


    # calc = ElectronicMinimize(atoms = slab)
    # slab.set_calculator(calc)
    #
    # slab.get_potential_energy()

def waterMolecule():
    atoms = molecule('H2O')
    print (atoms)
    atoms.center(vacuum=3.5)
    atoms.calc = ElectronicMinimize(atoms = atoms)
    print(atoms.get_potential_energy())
    print(atoms.get_forces())

def optimizeWater():
    d = 0.9575
    t = np.pi / 180 * 104.51
    water = Atoms('H2O',
                  positions=[(d, 0, 0),
                             (d * np.cos(t), d * np.sin(t), 0),
                             (0, 0, 0)],
                  cell = np.eye(3) * 7.0)
    water.calc = ElectronicMinimize(atoms = water, log = True)
    water.calc.dragWavefunctions = False
    # dyn = SciPyFminBFGS(water, callback_always=True)
    # dyn = BFGSLineSearch(water)
    # dyn = BFGS(water)
    # dyn.run(fmax=0.000005, steps = 10)
    print(water.get_potential_energy() / Hartree)

def compareForces():
    d = 0.9575
    t = np.pi / 180 * 104.51
    water = Atoms('H2O',
                  positions=[(d, 0, 0),
                             (d * np.cos(t), d * np.sin(t), 0),
                             (0, 0, 0)],
                  cell = np.eye(3) * 5.0)
    water.calc = ElectronicMinimize(atoms = water, log = False)
    print(water.get_forces())
    print(water.calc.calculate_numerical_forces(water))

def lithiumBulk(comm):
    d = 3.5
    li = Atoms('Li2',
                  positions=[(0, 0, 0),
                             (0.5*d, 0.5*d, 0.5*d)
                            ],
                  cell = np.eye(3) * d)
    from ase.dft import kpoints
    kpointList = kpoints.monkhorst_pack((4,4,4))
    kpointList += np.asarray([.125, .125, .125])
    w = 1.0/len(kpointList)
    if "--gpu" in sys.argv:
        li.calc = ElectronicMinimizeGPU(atoms = li, log = True,
                                 kpts = [(k,w) for k in kpointList],
                                 comm = comm)
    else:
        li.calc = ElectronicMinimize(atoms = li, log = True,
                                 kpts = [(k,w) for k in kpointList],
                                 comm = comm)
    li.calc.dragWavefunctions = False

    print(li.get_potential_energy() / Hartree)

def die():
    d = 0.9575
    t = np.pi / 180 * 104.51
    water = Atoms('H2O',
                  positions=[(d, 0, 0),
                             (d * np.cos(t), d * np.sin(t), 0),
                             (0, 0, 0)])
    water.calc = ElectronicMinimize(atoms = water, log = True)
    dyn = LBFGS(water)
    dyn.run(fmax=0.05)

if __name__=="__main__":
    comm = MPI.COMM_WORLD
    lithiumBulk(comm)
