import copy
import numpy as np

from ase.calculators.calculator import Calculator, all_changes
from ase import Atoms

from ase.units import Bohr, Hartree
from JDFTCalculator import JDFTCalculator

import sys

class ElectronicMinimize(Calculator,JDFTCalculator):
    implemented_properties = ['energy', 'forces']

    @staticmethod
    def _changeOrder(x, indexList):
        out = copy.copy(x)
        if isinstance(x, np.ndarray):
            for i, ind in enumerate(indexList):
                out[i] = x[ind]
            return out
        elif isinstance(x, Atoms):
            for i, ind in enumerate(indexList):
                out[i] = copy.copy(x[ind])
            return out
        else:
            raise TypeError("Can change the order of np.ndarray or ase.Atoms")

    @staticmethod
    def _createIndexLists(atoms):
        #JDFT has atoms ordered by their symbols so we need conversion tables of indices:
        symbols = {} #count number of occurances
        species_order = []
        for atom in atoms:
            try:
                symbols[atom.symbol] += 1
            except KeyError:
                species_order.append(atom.symbol)
                symbols[atom.symbol] = 0
        i = 0
        for sp in species_order:
            number_of_sp = symbols[sp] + 1
            symbols[sp] = i
            i += number_of_sp
        toJDFTOrderIndexList = [0] * len(atoms)
        fromJDFTOrderIndexList = [0] * len(atoms)
        for ind, atom in enumerate(atoms):
            toJDFTOrderIndexList[ind] = symbols[atom.symbol]
            fromJDFTOrderIndexList[symbols[atom.symbol]] = ind
            symbols[atom.symbol] += 1
        return (toJDFTOrderIndexList, fromJDFTOrderIndexList)

    def _toJDFTOrder(self, x):
        return self._changeOrder(x, self._toJDFTOrderIndexList)

    def _fromJDFTOrder(self, x):
        return self._changeOrder(x, self._fromJDFTOrderIndexList)

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label="JDFT", atoms=None, log = True, **kwargs):
        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            label, atoms, **kwargs)
        JDFTCalculator.__init__(self)
        if log == False:
            JDFTCalculator.disableLog(self)
        if atoms == None:
            return
        elif not isinstance(atoms, Atoms):
            raise TypeError("atoms should be ase.Atoms type.")

        self._toJDFTOrderIndexList, self._fromJDFTOrderIndexList = self._createIndexLists(atoms)
        self.R = atoms.cell
        for atom in atoms:
            self.add_ion(atom)
        self.setup()

    def updateAtomicPositions(self):
        dpos = self.atoms.positions - self._fromJDFTOrder(self.readIonicPositions()*Bohr)
        super().updateIonicPositions(self._toJDFTOrder(dpos/Bohr))

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        print("""Run one electronic minimize loop""")
        super().calculate(atoms, properties, system_changes)
        if 'positions' in system_changes:
            self.updateAtomicPositions()
        else:
            print(system_changes)
        self.runElecMin()
        energy = self.readTotalEnergy() * Hartree
        forces = np.asarray(self.readForces(), dtype = np.double)
        forces.resize((len(atoms),3))
        forces = forces * Hartree / Bohr
        self.results = {'energy': energy,
                       'forces': forces,
                       'stress': np.zeros(6),
                       'dipole': np.zeros(3),
                       'charges': np.zeros(len(atoms)),
                       'magmom': 0.0,
                       'magmoms': np.zeros(len(atoms))}
