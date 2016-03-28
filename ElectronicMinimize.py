#!/usr/bin/python3
# Author: Yalcin Ozhabes
# email: yalcinozhabes@gmail.com

import copy
import time
import numpy as np
from mpi4py import MPI

from ase.calculators.calculator import Calculator, all_changes
from ase import Atoms

from ase.units import Bohr, Hartree
from JDFTCalculatorCPU import JDFTCalculatorCPU
from JDFTCalculatorGPU import JDFTCalculatorGPU

class ElectronicMinimize(JDFTCalculatorCPU, Calculator):
    """
    A calculator derived from JDFTCalculatorCPU.
    """
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
        """JDFT has atoms ordered by their symbols so we need conversion tables
        of indices:"""
        symbols = {}  # count number of occurances
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
                 atoms=None, log=True, comm=None, **kwargs):
        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            "JDFT", atoms, **kwargs)
        super(ElectronicMinimize, self).__init__(comm=comm, log = log)

        if 'kpts' in kwargs:
            self.kpts = kwargs['kpts']
        if 'settings' in kwargs:
            self.settings = kwargs['settings']

        if atoms is None:
            return
        elif not isinstance(atoms, Atoms):
            raise TypeError("atoms should be ase.Atoms type.")

        self._toJDFTOrderIndexList, self._fromJDFTOrderIndexList = self._createIndexLists(atoms)
        self.R = atoms.cell

        if 'pseudopotential' in atoms.info:
            pspots = [atoms.info[pseudopotential]] * len(atoms)
        elif 'pseudopotentials' in atoms.info:
            pspots = atoms.info['pseudopotentials']
            assert len(pspots) == len(atoms)
        else:
            pspots = None
        for i, atom in enumerate(atoms):
            if pspots:
                atom.data['pseudopotential'] = pspots[i]
            self.add_ion(atom)
        t0 = time.time()
        c0 = time.clock()
        self.setup()
        print("Wall Time for e.setup()", time.time()-t0, "seconds")
        print("Process Time for e.setup()", time.clock()-c0, "seconds")

    def updateAtomicPositions(self):
        """"""
        dpos = self.atoms.positions - self._fromJDFTOrder(self.getIonicPositions() * Bohr)
        super(ElectronicMinimize, self).updateIonicPositions(self._toJDFTOrder(dpos / Bohr))

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        """Run one electronic minimize loop"""
        super(ElectronicMinimize, self).calculate(atoms, properties, system_changes)
        if 'positions' in system_changes:
            self.updateAtomicPositions()
        else:
            print(system_changes)

        t0 = time.time()
        c0 = time.clock()
        self.runElecMin()
        print("Wall Time for self.runElecMin()", time.time()-t0, "seconds")
        print("Process Time for self.runElecMin()", time.clock()-c0, "seconds")

        energy = self.readTotalEnergy() * Hartree
        forces = np.asarray(self.readForces(), dtype=np.double)
        forces.resize((len(atoms), 3))
        forces = forces * Hartree / Bohr
        self.results = {'energy': energy,
                        'forces': forces,
                        'stress': np.zeros(6),
                        'dipole': np.zeros(3),
                        'charges': np.zeros(len(atoms)),
                        'magmom': 0.0,
                        'magmoms': np.zeros(len(atoms))}

class ElectronicMinimizeGPU(JDFTCalculatorGPU, Calculator):
    """
    A calculator derived from JDFTCalculatorGPU.
    """
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
        """JDFT has atoms ordered by their symbols so we need conversion tables
        of indices:"""
        symbols = {}  # count number of occurances
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
                 atoms=None, log=True, comm=None, **kwargs):
        print "ElecMinGPU init running"
        Calculator.__init__(self, restart, ignore_bad_restart_file,
                            "JDFT", atoms, **kwargs)
        super(ElectronicMinimizeGPU, self).__init__(comm=comm, log = log)

        if 'kpts' in kwargs:
            self.kpts = kwargs["kpts"]

        if atoms is None:
            return
        elif not isinstance(atoms, Atoms):
            raise TypeError("atoms should be ase.Atoms type.")

        self._toJDFTOrderIndexList, self._fromJDFTOrderIndexList = self._createIndexLists(atoms)
        self.R = atoms.cell
        for atom in atoms:
            self.add_ion(atom)
        t0 = time.time()
        c0 = time.clock()
        self.setup()
        print("Wall Time for e.setup()", time.time()-t0, "seconds")
        print("Process Time for e.setup()", time.clock()-c0, "seconds")

    def updateAtomicPositions(self):
        """"""
        dpos = self.atoms.positions - self._fromJDFTOrder(self.getIonicPositions() * Bohr)
        super(ElectronicMinimizeGPU, self).updateIonicPositions(self._toJDFTOrder(dpos / Bohr))

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        """Run one electronic minimize loop"""
        super(ElectronicMinimizeGPU, self).calculate(atoms, properties, system_changes)
        if 'positions' in system_changes:
            self.updateAtomicPositions()
        else:
            print(system_changes)

        t0 = time.time()
        c0 = time.clock()
        self.runElecMin()
        print("Wall Time for self.runElecMin()", time.time()-t0, "seconds")
        print("Process Time for self.runElecMin()", time.clock()-c0, "seconds")

        energy = self.readTotalEnergy() * Hartree
        forces = np.asarray(self.readForces(), dtype=np.double)
        forces.resize((len(atoms), 3))
        forces = forces * Hartree / Bohr
        self.results = {'energy': energy,
                        'forces': forces,
                        'stress': np.zeros(6),
                        'dipole': np.zeros(3),
                        'charges': np.zeros(len(atoms)),
                        'magmom': 0.0,
                        'magmoms': np.zeros(len(atoms))}
