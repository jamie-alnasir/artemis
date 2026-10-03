#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#//==============================================================================
#// Artemis - Python molecular library for Brookhaven PDB files
#// By Jamie Alnasir, 07/2014
#// Royal Holloway University of London
#// CSSB - Centre for Systems and Synthetic Biology, dept. Computer Science
#// Copyright (c) 2014 Jamie J. Alnasir, All Rights Reserved
#//==============================================================================
#// Version: Python 2.7 edition - corrected October 2026
#//==============================================================================

# Requires Python 2.7; uses only the standard library.
# Existing public class/function names and the in-main implementation area remain.
# Corrections: fixed-width parsing; chain/insertion/model/TER separation; consistent
# alternate conformers; repeatable bonds with shared atoms; guarded geometry; and
# in-memory, field-preserving save with existing serials retained.
# HETATM records are preserved, but are not included in residue/bond analysis.
# Templates infer standard intra-residue bonds only; they are not a general
# bond-perception, protonation-state or chemical-valence validation engine.
# A retained serial-linked record that would become invalid causes saving to fail.
# Optional MASTER count records are omitted if atom/residue additions or removals
# invalidate them. Nonstandard modified residues remain RESTYPE.Unknown.

import sys;
import math;
import os;
import tempfile;


# The default example is to produce a report of the atoms in the PDB model by residue
# and list the bonds computed/found within those residues.

# Implementation to go in Implementation commented block towards the end of this
# source code file.

# Activities to perform (functional implementation provided):
_PRINT_RESIDUES_  = True;    # Print list of residue.
_PRINT_RES_BONDS_ = True;     # Print list of residue bonds (No effect if _PRINT_RESIDUES_=False).
_PRINT_DIHEDRALS_ = True;    # Print list of torsional/dihedral angles, by residue


def main():
    # Define our own Enum (Python 3.4 has it's own enum type)
    class Enum(object):
        def __init__(self, lstTuple):
            self.lstTuple = lstTuple

        def __getattr__(self, name):
            try:
                return self.lstTuple.index(name)
            except ValueError:
                raise AttributeError(name)

    # Enums for object attributes
    RESTYPE = Enum(('Unknown','AA', 'DNA'));
    HYBRIDISATION = Enum(('HYBRID_Unknown', 'HYBRID_SP1', 'HYBRID_SP2', 'HYBRID_SP3'));
    BOND_TYPE = Enum(('BOND_Ionic','BOND_Single','BOND_Double','BOND_Triple','BOND_HBond','BOND_VanDerWaals','BOND_Hydrogen','BOND_MeasureDist'));
    ATOM_SHELL_TYPE = Enum(('SHELL_S','SHELL_P','SHELL_D','SHELL_F'));

    # Consts
    KSTR_DNA_CODES = 'ATGCUI';

    # Tuples (NB: Tuples are immutable)
    from collections import namedtuple
    tplPDBAtom = namedtuple('tplPDBAtom', "serial,name,alt_loc,chain_id,res_name,res_seq,icode,x,y,z,occ,temp,element,charge");

    class TRealPoint:
        "3d vertex class"
        x=0;
        y=0;
        z=0;
        def __init__(self, a_x, a_y, a_z):
            self.x = a_x;
            self.y = a_y;
            self.z = a_z;

    def Distance3d(x1, y1, z1, x2, y2, z2):
        dx = (x2 - x1);
        dy = (y2 - y1);
        dz = (z2 - z1);
        tmp = dx * dx + dy * dy + dz * dz;
        return math.sqrt(tmp);

    def VectorSubtract(v1_x,v1_y,v1_z, v2_x,v2_y,v2_z):
    # Subtract V2 from V1
        return [v1_x - v2_x, v1_y - v2_y, v1_z - v2_z];

    def CrossProduct(v1_x, v1_y, v1_z,  v2_x, v2_y, v2_z):
            return [v1_y * v2_z - v1_z * v2_y, v1_z * v2_x - v1_x * v2_z, v1_x * v2_y - v1_y * v2_x];

    def DotProduct3d(v1_x, v1_y, v1_z,  v2_x, v2_y, v2_z):
    # Return Dot product of two vectors
        return v1_x * v2_x  +  v1_y * v2_y  +  v1_z * v2_z;

    def CosAfromDotProduct(aDP, aLine1Len, aLine2Len):
        if aLine1Len <= 0 or aLine2Len <= 0:
            raise ValueError("Angle undefined for a zero-length vector")
        result = float(aDP) / (aLine1Len * aLine2Len)
        if math.isnan(result) or math.isinf(result):
            raise ValueError("Non-finite angle input")
        return max(-1.0, min(1.0, result))

    def RadToDeg(Radians):
        return Radians * (180 / math.pi);


    def Angle(a0_x, a0_y, a0_z, a2_x, a2_y, a2_z, a3_x, a3_y, a3_z):
    # Three Atoms: a1 is Center atom (picked 1st), a2 and a3 are subsequently picked atoms
    # Lines are A and B which join at 0

        # dx, dy, dz;   # Difference in position (to get vector values)
        # lenA, lenB;   # Lengths
        # DProd, CosA;

        lenA = Distance3d(a0_x, a0_y, a0_z, a2_x, a2_y, a2_z);
        lenB = Distance3d(a0_x, a0_y, a0_z, a3_x, a3_y, a3_z);

        dx = a2_x - a0_x;
        dy = a2_y - a0_y;
        dz = a2_z - a0_z;
        vA_x = dx;
        vA_y = dy;
        vA_z = dz;

        dx = a3_x - a0_x;
        dy = a3_y - a0_y;
        dz = a3_z - a0_z;
        vB_x = dx;
        vB_y = dy;
        vB_z = dz;

        DProd = DotProduct3d(vA_x, vA_y, vA_z, vB_x, vB_y, vB_z);
        CosA  = CosAfromDotProduct(DProd, lenA, lenB);
        return RadToDeg(math.acos(CosA));


    # ---------------------------------------------------------------------------------------
    # Jamie Al-Nasir, Notes for computing torsional angles:
    # In order to calculate torsional/dihedral angles within a peptide we need to inspect
    # the x,y,z cartessian coordinates for various main chain atoms in the amino acids of
    # the first two residues (n, n+1) for Phi and Psi angles and the second and third residues
    # (n+1, n+2) in the case of Ohmega. Then we iterate n accordingly and repeat the same
    # procedure for the other residues in the chain. NB The first Phi and last Psi angles
    # cannot be computed for the first and last residues in the chain.
    # ---------------------------------------------------------------------------------------

    def DihedralAngle(vA_x, vA_y, vA_z, vB_x, vB_y, vB_z,
                      vC_x, vC_y, vC_z, vD_x, vD_y, vD_z):
        """Signed torsion in degrees, in [-180, 180]; reject degenerate planes."""
        b0 = VectorSubtract(vA_x, vA_y, vA_z, vB_x, vB_y, vB_z)
        b1 = VectorSubtract(vC_x, vC_y, vC_z, vB_x, vB_y, vB_z)
        b2 = VectorSubtract(vD_x, vD_y, vD_z, vC_x, vC_y, vC_z)
        length = math.sqrt(sum(t * t for t in b1))
        if length <= 1e-12:
            raise ValueError("Dihedral undefined: coincident central atoms")
        unit = [t / length for t in b1]
        d0 = sum(a * b for a, b in zip(b0, unit))
        d2 = sum(a * b for a, b in zip(b2, unit))
        v = [a - d0 * b for a, b in zip(b0, unit)]
        w = [a - d2 * b for a, b in zip(b2, unit)]
        if (sum(t * t for t in v) <= 1e-24 or
                sum(t * t for t in w) <= 1e-24):
            raise ValueError("Dihedral undefined: collinear atoms")
        cross = CrossProduct(*(unit + v))
        x = sum(a * b for a, b in zip(v, w))
        y = sum(a * b for a, b in zip(cross, w))
        if any(math.isnan(t) or math.isinf(t) for t in (x, y)):
            raise ValueError("Non-finite dihedral input")
        return RadToDeg(math.atan2(y, x))

    def DihedralAngleHandler(vA, vB, vC, vD):
    # Handle extraction of coordinates from vector lists
    # for call ti DihedralAngle function
        vA_x = vA[0]; vA_y = vA[1]; vA_z = vA[2];
        vB_x = vB[0]; vB_y = vB[1]; vB_z = vB[2];
        vC_x = vC[0]; vC_y = vC[1]; vC_z = vC[2];
        vD_x = vD[0]; vD_y = vD[1]; vD_z = vD[2];
        return DihedralAngle(vA_x, vA_y, vA_z, vB_x, vB_y, vB_z, vC_x, vC_y, vC_z, vD_x, vD_y, vD_z);

    def ListTupleByAttrVal(aList, aAttrIndex, aVal):
    # Searches aList of Tuples for attribute (tuple attr index) for matching aVal
        tmp = [i for i in aList if i[aAttrIndex] == aVal];
        if len(tmp) < 1:
            return None;
        else:
            return tmp[0];

    def CoordsListFromAtomTuple(aTuple):
        if aTuple is not None:
            return [float(aTuple.x), float(aTuple.y), float(aTuple.z)];
        else:
            raise ValueError("Missing atom has no coordinates")

    def CreateAtomFromTuple(aAtomTuple):
        if not aAtomTuple is None:
            anAtom = TAtom(aAtomTuple);
            return anAtom;
        else:
            return None;

    AA_NAMES = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())

    def NucleotideCode(name):
        name = name.strip().upper()
        if name in ("DA", "DT", "DG", "DC", "DU", "DI"):
            name = name[1:]
        return name if len(name) == 1 and name in KSTR_DNA_CODES else None

    def IsDNARes(aResName):
        # Historical API name: includes recognised RNA residues as well as DNA.
        return NucleotideCode(aResName) is not None

    def GetPDBResType(aResName):
        if IsDNARes(aResName):
            return RESTYPE.DNA
        return RESTYPE.AA if aResName.strip().upper() in AA_NAMES else RESTYPE.Unknown

    class TAtom(object):
        "Atom object class"
        __slots__ = ('xyz', 'bBackBone', 'bHydrogen', 'bComputed', 'Hybridisation', 'PDBData');
        def __init__(self, PDBAtomTuple):
            self.xyz = TRealPoint(0,0,0);
            self.bBackBone = 0;
            self.bHydrogen = 0;
            self.bComputed = 0;
            self.Hybridisation = HYBRIDISATION.HYBRID_Unknown;
            self.PDBData = PDBAtomTuple;
            if PDBAtomTuple is not None:
                self.xyz.x = float(PDBAtomTuple.x);
                self.xyz.y = float(PDBAtomTuple.y);
                self.xyz.z = float(PDBAtomTuple.z);

    class TBond(object):
        "Bond object class"
        def __init__(self, a1=None, a2=None):
            self.Atom1 = a1;
            self.Atom2 = a2;
            self.SciAtomBondType   = BOND_TYPE.BOND_Single;
            self.bBackBone = 0;
            self.bComputed = 0;

        # Keep the historical name as an alias, never a second independent value.
        @property
        def BondType(self):
            return self.SciAtomBondType

        @BondType.setter
        def BondType(self, value):
            self.SciAtomBondType = value

        def getBondName(self):
            r = "";
            if self.Atom1 is None or self.Atom2 is None:
                return r;
            if self.Atom1.PDBData is None or self.Atom2.PDBData is None:
                return "Error: Atom.PDBData is null" + str(type(self.Atom1.PDBData));
            r = self.Atom1.PDBData.name + " - " + self.Atom2.PDBData.name;
            return r;

    class TMolecule(object):
        "Molecule object base class"
        def __init__(self):
            self.lstBonds = [];
        def addBond(self, aBond):
            if not (aBond.Atom1 is None or aBond.Atom2 is None):
                key = frozenset((aBond.Atom1.PDBData, aBond.Atom2.PDBData))
                if not any(frozenset((b.Atom1.PDBData, b.Atom2.PDBData)) == key
                           for b in self.lstBonds):
                    self.lstBonds.append(aBond);

    class TResidue(TMolecule):
        "Residue object class"
        def __init__(self):
            TMolecule.__init__(self);
            self.res_name = "";
            self.res_seq = "";
            self.lstAtoms = [];
            self.chain_id = "";
            self.res_type = RESTYPE.Unknown;
            self.icode = ""
            self.model_index = 0
            self.model_id = "1"
            self.segment_id = 0
            self._atom_cache = {}


        def addAtomTuple(self, PDBAtomTuple):
            self.lstAtoms.append(PDBAtomTuple);

        def SelectedAtoms(self):
            """Use shared blank-altLoc atoms plus one consistent residue conformer.

            Choose the label with the highest mean occupancy; ties prefer A,
            then lexical order. Preserve every conformer in lstAtoms and on save.
            """
            scores = {}
            for atom in self.lstAtoms:
                if atom.alt_loc:
                    scores.setdefault(atom.alt_loc, []).append(float(atom.occ or 0))
            label = None
            if scores:
                label = sorted(scores, key=lambda k:
                               (-sum(scores[k]) / len(scores[k]), k != "A", k))[0]
            selected = {}
            for atom in self.lstAtoms:
                if atom.alt_loc not in ("", label):
                    continue
                name = atom.name.replace("*", "'")
                name = {"O1P": "OP1", "O2P": "OP2", "O3P": "OP3"}.get(name, name)
                old = selected.get(name)
                if old is None or (old.alt_loc and not atom.alt_loc):
                    selected[name] = atom
            return selected

        def AtomByName(self, name):
            name = name.replace("*", "'")
            name = {"O1P": "OP1", "O2P": "OP2", "O3P": "OP3"}.get(name, name)
            atom = self._selected_atoms.get(name)
            if atom is None:
                return None
            if name not in self._atom_cache:
                self._atom_cache[name] = TAtom(atom)
            return self._atom_cache[name]

        def ComputeResBonds(self):
            self.lstBonds = []
            self._atom_cache = {}
            self._selected_atoms = self.SelectedAtoms()
            res_name = NucleotideCode(self.res_name) or self.res_name.upper()
            if GetPDBResType(self.res_name) == RESTYPE.Unknown:
                return

            # Adapted from the Zeus Molecular visualisation framework developed by Jamie Al-Nasir

            # Amino Acids ===============================================================

            # Check a-Amino backbone
            # Add bonds for a-Amino backbone, N-C-C=O
            # Covers GLY

            # [N-C]-C=O
            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("N");
            Bond.Atom2 = self.AtomByName("CA");
            Bond.bBackBone = 1;
            self.addBond(Bond);

            # N-[C-C]=0
            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("CA");
            Bond.Atom2 = self.AtomByName("C");
            Bond.bBackBone = 1;
            self.addBond(Bond);

            # N-C-[C=0]
            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Double;
            Bond.Atom1 = self.AtomByName("C");
            Bond.Atom2 = self.AtomByName("O");

            # C= is normally SP2 Hybridised
            if not Bond.Atom1 is None:
                Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # C
            Bond.bBackBone = 1;
            self.addBond(Bond);

            # == ERROR a-backbone

            # Check a-Amino alkyl side-chain
            # Add bonds for alkyl side-chain, i.e. alpha-beta C, beta-gamma c, gamma-delta C
            # N-C-C=O
            #[|]
            #[C]
            # Caters for ALA, ILE
            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("CA");
            Bond.Atom2 = self.AtomByName("CB");
            if (Bond.Atom1 is not None and Bond.Atom2 is not None):
                Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CA is SP3
                Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CB is SP3
                self.addBond(Bond);
            else:
                del Bond;

            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("CB");
            Bond.Atom2 = self.AtomByName("CG");
            if (Bond.Atom1 is not None and Bond.Atom2 is not None):
                Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CA is SP3
                Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CG is SP3
                self.addBond(Bond);
            else:
                del Bond;

            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("CG");
            Bond.Atom2 = self.AtomByName("CD");
            if (Bond.Atom1 is not None and Bond.Atom2 is not None):
                Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CG is SP3
                Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CD is SP3
                self.addBond(Bond);
            else:
                del Bond;
            # End Alkyl backbone

            if (res_name == "ASP"):
                # Add COO- to D-carbon of Alkyl side-chain
                # which one is the C=O double bond?? OD1 or OD2? use OD1 for now
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("OD1");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("OD2");
                self.addBond(Bond);

            if (res_name == "GLU"):
            # Add COO- to E-carbon of Alkyl side-chain
                # which one is the C=O double bond?? OD1 or OD2? use OD1 for now
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CD");
                Bond.Atom2 = self.AtomByName("OE1");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CD");
                Bond.Atom2 = self.AtomByName("OE2");
                self.addBond(Bond);

            if (res_name == "ILE"):
            # Add B-Carbon to 2x G-Carbon bonds and one G-Carbon-D-Carbon bond
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CB");
                Bond.Atom2 = self.AtomByName("CG1");
                if not Bond.Atom1 is None:
                    Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CB SP3 hybridised
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CG1 SP3 hybridised
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CB");
                Bond.Atom2 = self.AtomByName("CG2");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CG2 SP3 hybridised
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CG1");
                Bond.Atom2 = self.AtomByName("CD1");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CD1 SP3 hybridised
                self.addBond(Bond);

            if (res_name == "LEU"):
            # Add 2x G-Carbon to D-Carbon bonds
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("CD1");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CD1 SP3 hybridised
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("CD2");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CD2 SP3 hybridised
                self.addBond(Bond);


            if (res_name == "MET"):
            # Add G-Carbon to D-Sulphur and D-Sulphur to E-Carbon bonds
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("SD");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("SD");
                Bond.Atom2 = self.AtomByName("CE");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CE SP3 Hybridised
                self.addBond(Bond);

            if (res_name == "VAL"):
            # Add 2x B-Carbon to G-Carbon bonds
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CB");
                Bond.Atom2 = self.AtomByName("CG1");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CB");
                Bond.Atom2 = self.AtomByName("CG2");
                self.addBond(Bond);

            if (res_name == "LYS"):
            # Add D-Carbon to E-Carbon bond and E-Carbon to Z-Nitrogen (NH3) bond
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CD");
                Bond.Atom2 = self.AtomByName("CE");
                if not Bond.Atom1 is None:
                    Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;
                self.addBond(Bond);
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CE");
                Bond.Atom2 = self.AtomByName("NZ");
                self.addBond(Bond);


            if (res_name == "SER"):
            # Add B-Carbon to G-Oxygen (OH) bond
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CB");
                Bond.Atom2 = self.AtomByName("OG");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # OG SP3 Hybridised (4 SP3 orbitals, 2 lone pairs, so bent)
                self.addBond(Bond);

            if (res_name == "THR"):
            # Add B-Carbon to G-Oxygen (OH) and B-Carbon to G-Carbon bonds
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CB");
                Bond.Atom2 = self.AtomByName("OG1");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CB");
                Bond.Atom2 = self.AtomByName("CG2");
                self.addBond(Bond);

            if (res_name == "ASN"):
            # Add G-Carbon to D-Oxygen (C=O) bond
            # Add G-Carbon to D-Nitrogen (NH3) bond
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("OD1");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("ND2");
                self.addBond(Bond);

            if (res_name == "GLN"):
            # Add D-Carbon to E-Oxygen (C=O) bond
            # Add D-Carbon to E-Nitrogen (NH3) bond
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CD");
                Bond.Atom2 = self.AtomByName("OE1");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CD");
                Bond.Atom2 = self.AtomByName("NE2");
                self.addBond(Bond);

            if (res_name == "ARG"):
            # Add D-Carbon to E-Nitrogen bond
            # Add E-Nitrogen to Z-Carbon bond
            # Add Z-Carbon double bond to H-N (NH2+)
            # Add Z-Carbon single bond to H-N (NH2)
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CD");
                Bond.Atom2 = self.AtomByName("NE");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("NE");
                Bond.Atom2 = self.AtomByName("CZ");
                self.addBond(Bond);

                # which one is the C=N double bond?? NH1 or NH2? use NH1 for now
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CZ");
                Bond.Atom2 = self.AtomByName("NH1");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CZ");
                Bond.Atom2 = self.AtomByName("NH2");
                self.addBond(Bond);

            if (res_name == "CYS"):
            # Add B-Carbon to G-Sulphur bond
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CB");
                Bond.Atom2 = self.AtomByName("SG");
                self.addBond(Bond);

            if (res_name == "PRO"):
            # Add D-Carbon to Nitrogen bond which closes the ring in Proline
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CD");
                Bond.Atom2 = self.AtomByName("N");
                self.addBond(Bond);

            if (res_name == "PHE"):
            # Create Phenyl Ring from two parallel alkyl-chains which meet up at Z-Carbon
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("CD1");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("CD2");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CD1");
                Bond.Atom2 = self.AtomByName("CE1");
                if not Bond.Atom1 is None:
                    Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CD2");
                Bond.Atom2 = self.AtomByName("CE2");
                if not Bond.Atom1 is None:
                    Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CE1");
                Bond.Atom2 = self.AtomByName("CZ");
                if not Bond.Atom1 is None:
                    Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CE2");
                Bond.Atom2 = self.AtomByName("CZ");
                self.addBond(Bond);

            if (res_name == "TYR"):
                # Create Phenyl Ring from two parallel alkyl-chains which meet up at Z-Carbon
                # Add bond to Z-Carbon to H-O (Phenolic-OH)
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("CD1");
                if not Bond.Atom1 is None:
                    Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CD1 SP2 hybridised
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CD1 SP2 hybridised
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("CD2");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CD2 SP2 hybridised
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CD1");
                Bond.Atom2 = self.AtomByName("CE1");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CE1 SP2 hybridised
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CD2");
                Bond.Atom2 = self.AtomByName("CE2");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CE2 SP2 hybridised
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CE1");
                Bond.Atom2 = self.AtomByName("CZ");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CE2");
                Bond.Atom2 = self.AtomByName("CZ");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CZ");
                Bond.Atom2 = self.AtomByName("OH");
                self.addBond(Bond);

            if (res_name == "HIS"):
                # Create Hetero Ring from two parallel alkyl-chains which meet up at Z-Carbon
                # Add bond to Z-Carbon to H-O (Phenolic-OH)
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("CD2");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CD2 SP2 hybridised
                self.addBond(Bond);


                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("ND1");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CD2 SP2 hybridised/protonated, + charge
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("ND1");
                Bond.Atom2 = self.AtomByName("CE1");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CE1 SP2 hybridised
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CD2");
                Bond.Atom2 = self.AtomByName("NE2");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("NE2");
                Bond.Atom2 = self.AtomByName("CE1");
                self.addBond(Bond);

            if (res_name == "TRP"):
                # Create Phenyl Ring from two parallel alkyl-chains which meet up at Z-Carbon
                # Add bond to Z-Carbon to H-O (Phenolic-OH)

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("CD1");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CD1 SP2 hybridised
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CG");
                Bond.Atom2 = self.AtomByName("CD2");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CD1");
                Bond.Atom2 = self.AtomByName("NE1");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("NE1");
                Bond.Atom2 = self.AtomByName("CE2");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CE2");
                Bond.Atom2 = self.AtomByName("CD2");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CE2");
                Bond.Atom2 = self.AtomByName("CZ2");
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CZ2");
                Bond.Atom2 = self.AtomByName("CH2");
                if not Bond.Atom1 is None:
                    Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CZ2 SP2 hybridised
                self.addBond(Bond);

                #
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CH2");
                Bond.Atom2 = self.AtomByName("CZ3");
                if not Bond.Atom1 is None:
                    Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CH2 SP2 hybridised
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("CZ3");
                Bond.Atom2 = self.AtomByName("CE3");
                if not Bond.Atom1 is None:
                    Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CZ3 SP2 hybridised
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CE3 SP2 hybridised
                self.addBond(Bond);

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("CE3");
                Bond.Atom2 = self.AtomByName("CD2");
                self.addBond(Bond);

            # Nucleic Acids =============================================================

            if (res_name == "C"):
            # Cytosine
                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C6");
                Bond.Atom2 = self.AtomByName("C5");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C5");
                Bond.Atom2 = self.AtomByName("C4");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C4");
                Bond.Atom2 = self.AtomByName("N3");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C4");
                Bond.Atom2 = self.AtomByName("N4");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N3");
                Bond.Atom2 = self.AtomByName("C2");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("O2");
                Bond.Atom2 = self.AtomByName("C2");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C2");
                Bond.Atom2 = self.AtomByName("N1");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N1");
                Bond.Atom2 = self.AtomByName("C6");
                Bond.bBackBone = False;
                self.addBond(Bond);

            if (res_name in ("T", "U")):
            # Thymine

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N1");
                Bond.Atom2 = self.AtomByName("C2");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C2");
                Bond.Atom2 = self.AtomByName("O2");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C2");
                Bond.Atom2 = self.AtomByName("N3");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N3");
                Bond.Atom2 = self.AtomByName("C4");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C4");
                Bond.Atom2 = self.AtomByName("C5");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C4");
                Bond.Atom2 = self.AtomByName("O4");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C5");
                Bond.Atom2 = self.AtomByName("C6");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N1");
                Bond.Atom2 = self.AtomByName("C6");
                Bond.bBackBone = False;
                self.addBond(Bond);

            if (res_name == "A"):
            # Adenine

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N9");
                Bond.Atom2 = self.AtomByName("C8");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C8");
                Bond.Atom2 = self.AtomByName("N7");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N7");
                Bond.Atom2 = self.AtomByName("C5");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C5");
                Bond.Atom2 = self.AtomByName("C6");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C6");
                Bond.Atom2 = self.AtomByName("N6");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C6");
                Bond.Atom2 = self.AtomByName("N1");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N1");
                Bond.Atom2 = self.AtomByName("C2");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C2");
                Bond.Atom2 = self.AtomByName("N3");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C5");
                Bond.Atom2 = self.AtomByName("C4");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C4");
                Bond.Atom2 = self.AtomByName("N9");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C4");
                Bond.Atom2 = self.AtomByName("N3");
                Bond.bBackBone = False;
                self.addBond(Bond);

            if (res_name == "G"):
            # Guanine

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N9");
                Bond.Atom2 = self.AtomByName("C8");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C8");
                Bond.Atom2 = self.AtomByName("N7");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N7");
                Bond.Atom2 = self.AtomByName("C5");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C5");
                Bond.Atom2 = self.AtomByName("C6");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C6");
                Bond.Atom2 = self.AtomByName("O6");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C6");
                Bond.Atom2 = self.AtomByName("N1");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N1");
                Bond.Atom2 = self.AtomByName("C2");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C2");
                Bond.Atom2 = self.AtomByName("N3");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C5");
                Bond.Atom2 = self.AtomByName("C4");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C4");
                Bond.Atom2 = self.AtomByName("N9");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C4");
                Bond.Atom2 = self.AtomByName("N3");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C2");
                Bond.Atom2 = self.AtomByName("N2");
                Bond.bBackBone = False;
                self.addBond(Bond);

            if (res_name == "I"):
            # Inosine

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N9");
                Bond.Atom2 = self.AtomByName("C8");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C8");
                Bond.Atom2 = self.AtomByName("N7");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N7");
                Bond.Atom2 = self.AtomByName("C5");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C5");
                Bond.Atom2 = self.AtomByName("C6");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C6");
                Bond.Atom2 = self.AtomByName("O6");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C6");
                Bond.Atom2 = self.AtomByName("N1");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N1");
                Bond.Atom2 = self.AtomByName("C2");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C2");
                Bond.Atom2 = self.AtomByName("N3");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Double;
                Bond.Atom1 = self.AtomByName("C5");
                Bond.Atom2 = self.AtomByName("C4");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C4");
                Bond.Atom2 = self.AtomByName("N9");
                Bond.bBackBone = False;
                self.addBond(Bond);

                # [N-C]-C=O
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("C4");
                Bond.Atom2 = self.AtomByName("N3");
                Bond.bBackBone = False;
                self.addBond(Bond);

            if (res_name == "A" or res_name == "G" or res_name == "I"):
            # Attach Purine to Sugar Phosphate

                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N9");
                Bond.Atom2 = self.AtomByName("C1*");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # C1* SP3 Hybridised
                Bond.bBackBone = False;
                self.addBond(Bond);

            if (res_name in ("C", "T", "U")):
            # Attach Pyramidine to Sugar Phosphate
                Bond = TBond();
                Bond.BondType = BOND_TYPE.BOND_Single;
                Bond.Atom1 = self.AtomByName("N1");
                Bond.Atom2 = self.AtomByName("C1*");
                if not Bond.Atom2 is None:
                    Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # C1* SP3 Hybridised
                Bond.bBackBone = False;
                self.addBond(Bond);

            # Phosphate Backbone
            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("O1P");
            Bond.Atom2 = self.AtomByName("P");
            Bond.bBackBone = 1;
            self.addBond(Bond);

            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Double;
            Bond.Atom1 = self.AtomByName("O2P");
            Bond.Atom2 = self.AtomByName("P");
            Bond.bBackBone = 1;
            self.addBond(Bond);


            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("C1*");
            Bond.Atom2 = self.AtomByName("C2*");
            if not Bond.Atom1 is None:
                Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # C1* SP3 Hybridised
            if not Bond.Atom2 is None:
                Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # C2* SP3 Hybridised
            Bond.bBackBone = 1;
            self.addBond(Bond);

            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("C2*");
            Bond.Atom2 = self.AtomByName("C3*");
            if not Bond.Atom2 is None:
                Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # C1* SP3 Hybridised
            Bond.bBackBone = 1;
            self.addBond(Bond);

            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("C3*");
            Bond.Atom2 = self.AtomByName("C4*");
            Bond.bBackBone = 1;
            self.addBond(Bond);

            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("C4*");
            Bond.Atom2 = self.AtomByName("O4*");
            Bond.bBackBone = 1;
            self.addBond(Bond);

            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("O4*");
            Bond.Atom2 = self.AtomByName("C1*");
            Bond.bBackBone = 1;
            self.addBond(Bond);

            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("C4*");
            Bond.Atom2 = self.AtomByName("C5*");
            Bond.bBackBone = 1;
            self.addBond(Bond);

            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("C3*");
            Bond.Atom2 = self.AtomByName("O3*");
            Bond.bBackBone = 1;
            self.addBond(Bond);

            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("C5*");
            Bond.Atom2 = self.AtomByName("O5*");
            if not Bond.Atom1 is None:
                Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # C5* SP3 Hybridised
            Bond.bBackBone = 1;
            self.addBond(Bond);

            Bond = TBond();
            Bond.BondType = BOND_TYPE.BOND_Single;
            Bond.Atom1 = self.AtomByName("P");
            Bond.Atom2 = self.AtomByName("O5*");
            Bond.bBackBone = 1;
            self.addBond(Bond);



            # Additional standard terminal/sugar substituents.
            for first, second in (("C", "OXT"), ("C5", "C7"), ("C2*", "O2*"), ("P", "O3P")):
                if first == "C" and res_name not in AA_NAMES:
                    continue
                if first == "C5" and res_name != "T":
                    continue
                if first in ("C2*", "P") and not IsDNARes(self.res_name):
                    continue
                self.addBond(TBond(self.AtomByName(first), self.AtomByName(second)))

            # Hybridisation belongs to the shared atom, not an independent copy
            # in each bond. Resolve it after all bond orders are known.
            for atom in self._atom_cache.values():
                atom.Hybridisation = HYBRIDISATION.HYBRID_Unknown
            for bond in self.lstBonds:
                for atom in (bond.Atom1, bond.Atom2):
                    if bond.BondType == BOND_TYPE.BOND_Double:
                        atom.Hybridisation = HYBRIDISATION.HYBRID_SP2
            if res_name in AA_NAMES:
                planar = set(["C", "O", "OXT"])
                planar.update({
                    "ARG": ["NE", "CZ", "NH1", "NH2"],
                    "ASN": ["CG", "OD1", "ND2"],
                    "ASP": ["CG", "OD1", "OD2"],
                    "GLN": ["CD", "OE1", "NE2"],
                    "GLU": ["CD", "OE1", "OE2"],
                    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
                    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
                    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
                    "TRP": ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]
                }.get(res_name, []))
                for name, atom in self._atom_cache.items():
                    if name in planar:
                        atom.Hybridisation = HYBRIDISATION.HYBRID_SP2
                    elif name.startswith("C"):
                        atom.Hybridisation = HYBRIDISATION.HYBRID_SP3
            else:
                for name, atom in self._atom_cache.items():
                    if "'" in name and name.startswith("C"):
                        atom.Hybridisation = HYBRIDISATION.HYBRID_SP3
                    elif name in ("N1", "N2", "N3", "N4", "N6", "N7", "N9", "C2", "C4", "C5", "C6", "C7", "C8", "O2", "O4", "O6"):
                        atom.Hybridisation = (HYBRIDISATION.HYBRID_SP3 if name == "C7"
                                              else HYBRIDISATION.HYBRID_SP2)

    # Brookhaven PDB (Protein Databank file) format:
    # PDB file fixed-width columns
    #1 -  6       Record name      "ATOM    "
    #7 - 11       Integer          serial     Atom serial number.
    #13 - 16      Atom             name       Atom name.
    #17           Character        altLoc     Alternate location indicator.
    #18 - 20      Residue name     resName    Residue name.
    #22           Character        chainID    Chain identifier.
    #23 - 26      Integer          resSeq     Residue sequence number.
    #27           AChar            iCode      Code for insertion of residues.
    #31 - 38      Real(8.3)        x          Orthogonal coordinates for X in
    #                                         Angstroms
    #39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in
    #                                         Angstroms
    #47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in
    #                                         Angstroms
    #55 - 60      Real(6.2)        occupancy  Occupancy.
    #61 - 66      Real(6.2)        tempFactor Temperature factor.
    #77 - 78      LString(2)       element    Element symbol, right-justified.
    #79 - 80      LString(2)       charge     Charge on the atom.

    class TPDBModel(object):
        """ATOM-based model; HETATM and other records are preserved on saving.

        All MODEL blocks are kept separate. Bond templates cover recognised
        amino acids and nucleotides, not arbitrary ligands or inter-residue bonds.
        Edit the immutable atom tuples with _replace() in a residue's lstAtoms.
        Existing serials are preserved; use serial='auto' for newly added atoms.
        """
        def __init__(self):
            self._data = []
            self._pdb_file = ""
            self.lstAtomTuples = []
            self.lstChains = []
            self.lstRes = []
            self._records = []
            self._explicit_models = False

        def _load(self):
            self.lstAtomTuples = []
            self.lstChains = []
            self.lstRes = []
            self._records = []
            self._explicit_models = any(line[:6].strip() == 'MODEL' for line in self._data)
            model_index = -1 if self._explicit_models else 0
            model_id = "1"
            in_model = not self._explicit_models
            segment = 0
            previous_key = None
            residue = None
            for line_number, source_line in enumerate(self._data, 1):
                line = source_line.rstrip('\r\n')
                kind = line[:6].strip()
                if kind == 'MODEL':
                    if in_model:
                        raise ValueError('Nested MODEL at line %d' % line_number)
                    model_index += 1
                    model_id = line[10:14].strip() or str(model_index + 1)
                    segment = 0
                    previous_key = None
                    in_model = True
                elif kind == 'ENDMDL':
                    if not self._explicit_models or not in_model:
                        raise ValueError('Unmatched ENDMDL at line %d' % line_number)
                    in_model = False
                    previous_key = None
                elif kind == 'TER':
                    segment += 1
                    previous_key = None
                if kind != 'ATOM':
                    self._records.append(('raw', line, model_index))
                    continue
                if not in_model:
                    raise ValueError('ATOM outside MODEL at line %d' % line_number)
                if len(line) < 54:
                    raise ValueError('Short ATOM record at line %d' % line_number)
                atom = tplPDBAtom(
                    line[6:11].strip(), line[12:16].strip(), line[16:17].strip(),
                    line[21:22].strip(), line[17:20].strip(), line[22:26].strip(),
                    line[26:27].strip(), line[30:38].strip(), line[38:46].strip(),
                    line[46:54].strip(), line[54:60].strip(), line[60:66].strip(),
                    line[76:78].strip(), line[78:80].strip())
                try:
                    int(atom.serial)
                    int(atom.res_seq)
                    for value in (atom.x, atom.y, atom.z):
                        if math.isnan(float(value)) or math.isinf(float(value)):
                            raise ValueError('non-finite coordinate')
                    for value in (atom.occ, atom.temp):
                        if value and (math.isnan(float(value)) or math.isinf(float(value))):
                            raise ValueError('non-finite optional field')
                    if not atom.name or not atom.res_name:
                        raise ValueError('missing atom or residue name')
                except ValueError as error:
                    raise ValueError('Invalid ATOM at line %d: %s' % (line_number, error))
                key = (model_index, segment, atom.chain_id, atom.res_seq,
                       atom.icode, atom.res_name)
                if key != previous_key:
                    residue = TResidue()
                    residue.res_name = atom.res_name
                    residue.res_seq = atom.res_seq
                    residue.icode = atom.icode
                    residue.chain_id = atom.chain_id
                    residue.model_index = model_index
                    residue.model_id = model_id
                    residue.segment_id = segment
                    residue.res_type = GetPDBResType(atom.res_name)
                    self.lstRes.append(residue)
                    previous_key = key
                residue.addAtomTuple(atom)
                self.lstAtomTuples.append(atom)
                if atom.chain_id not in self.lstChains:
                    self.lstChains.append(atom.chain_id)
                self._records.append(('atom', residue, atom, line, model_index))
            if self._explicit_models and in_model:
                raise ValueError('MODEL block has no ENDMDL')

        def getResByChainID(self, chain_id, model_index=None, segment_id=None):
            return [res for res in self.lstRes if res.chain_id == chain_id
                    and (model_index is None or res.model_index == model_index)
                    and (segment_id is None or res.segment_id == segment_id)]

        def LoadFromFile(self, pdb_file):
            with open(pdb_file, 'r') as handle:
                data = handle.readlines()
            self._pdb_file = pdb_file
            self._data = data
            self._load()

        def LoadFromStream(self):
            self._pdb_file = ""
            self._data = list(sys.stdin)
            self._load()

        def _format_atom(self, atom, serial, original=''):
            # Retain segment IDs, original atom-name alignment and extra columns.
            columns = list(original.ljust(80))
            def put(start, end, value, right=False):
                value = str(value)
                width = end - start
                if len(value) > width or '\n' in value or '\r' in value:
                    raise ValueError('PDB field exceeds columns %d-%d: %r' %
                                     (start + 1, end, value))
                columns[start:end] = list(value.rjust(width) if right else value.ljust(width))
            def real(value, width, decimals, optional=False):
                if optional and str(value).strip() == '':
                    return ''
                number = float(value)
                if math.isnan(number) or math.isinf(number):
                    raise ValueError('Non-finite PDB numeric field')
                result = ('%*.*f' % (width, decimals, number))
                if len(result) > width:
                    raise ValueError('Numeric value does not fit PDB field: %s' % value)
                return result
            put(0, 6, 'ATOM')
            put(6, 11, str(serial), True)
            put(11, 12, '')
            if original and original[12:16].strip() == atom.name:
                name_field = original[12:16]
            elif len(atom.name) < 4 and not atom.name[:1].isdigit() and len(atom.element.strip()) != 2:
                name_field = ' ' + atom.name
            else:
                name_field = atom.name
            put(12, 16, name_field)
            put(16, 17, atom.alt_loc)
            put(17, 20, atom.res_name, True)
            put(20, 21, '')
            put(21, 22, atom.chain_id)
            put(22, 26, str(int(atom.res_seq)), True)
            put(26, 27, atom.icode)
            put(27, 30, '')
            for start, value in ((30, atom.x), (38, atom.y), (46, atom.z)):
                put(start, start + 8, real(value, 8, 3), True)
            put(54, 60, real(atom.occ, 6, 2, True), True)
            put(60, 66, real(atom.temp, 6, 2, True), True)
            put(76, 78, atom.element, True)
            put(78, 80, atom.charge, True)
            return ''.join(columns)

        def RebuildPDB(self):
            """Rebuild in memory without reopening the source or renumbering atoms.

            Existing non-ATOM records retain their positions. A change that would
            invalidate a retained serial-linked record is rejected before writing.
            New residues are inserted before ENDMDL or the final connectivity/tail
            records. Set model_index explicitly for additions to multi-model files.
            """
            original_residues = set(record[1] for record in self._records if record[0] == 'atom')
            current_residues = set(self.lstRes)
            if len(current_residues) != len(self.lstRes):
                raise ValueError('The same residue object occurs more than once')
            known_models = set(record[2] for record in self._records
                               if record[0] == 'raw' and record[1][:6].strip() == 'MODEL')
            if not self._explicit_models:
                known_models = set([0])
            reserved = {}
            original_by_serial = {}
            for record in self._records:
                if record[0] == 'atom':
                    model = record[4]
                    serial = int(record[2].serial)
                    reserved.setdefault(model, set()).add(serial)
                    key = (model, serial)
                    if key in original_by_serial:
                        raise ValueError('Duplicate atom serial within a model: %s' % serial)
                    original_by_serial[key] = record[2]
                else:
                    line, model = record[1:]
                    if line[:6].strip() in ('HETATM', 'TER') and line[6:11].strip():
                        reserved.setdefault(model, set()).add(int(line[6:11]))
            # Reserve all explicit serials first so auto-assigned values cannot clash.
            seen = {}
            for record in self._records:
                if record[0] == 'raw' and record[1][:6].strip() in ('HETATM', 'TER'):
                    value = record[1][6:11].strip()
                    if value:
                        seen.setdefault(record[2], set()).add(int(value))
            for res in self.lstRes:
                if res.model_index not in known_models:
                    raise ValueError('New residue refers to an unknown model')
                for atom in res.lstAtoms:
                    if str(atom.serial).strip().lower() not in ('', 'auto'):
                        serial = int(atom.serial)
                        if not 1 <= serial <= 99999:
                            raise ValueError('Atom serial must be between 1 and 99999')
                        if serial in seen.setdefault(res.model_index, set()):
                            raise ValueError('Duplicate serial in model: %s' % serial)
                        seen[res.model_index].add(serial)
                        reserved.setdefault(res.model_index, set()).add(serial)
            serials = {}
            current_by_serial = {}
            for res in self.lstRes:
                for index, atom in enumerate(res.lstAtoms):
                    if str(atom.serial).strip().lower() in ('', 'auto'):
                        used = reserved.setdefault(res.model_index, set())
                        serial = max(used or set([0])) + 1
                        if serial > 99999:
                            raise ValueError('No serial space remaining in this PDB model')
                        used.add(serial)
                    else:
                        serial = int(atom.serial)
                    serials[(res, index)] = serial
                    current_by_serial[(res.model_index, serial)] = atom
            # Do not silently leave CONECT/ANISOU/etc pointing to renamed/deleted atoms.
            def identity(atom):
                return (atom.name, atom.alt_loc, atom.chain_id, atom.res_name,
                        atom.res_seq, atom.icode)
            for record in self._records:
                if record[0] != 'raw':
                    continue
                line, model = record[1:]
                kind = line[:6].strip()
                if kind not in ('CONECT', 'ANISOU', 'SIGATM', 'SIGUIJ'):
                    continue
                fields = ([line[i:i+5] for i in range(6, len(line), 5)]
                          if kind == 'CONECT' else [line[6:11]])
                for field in fields:
                    if not field.strip():
                        continue
                    serial = int(field)
                    candidates = ([key for key in original_by_serial if key[1] == serial]
                                  if kind == 'CONECT' else [(model, serial)])
                    for key in candidates:
                        old = original_by_serial.get(key)
                        new = current_by_serial.get(key)
                        if old is not None and (new is None or identity(old) != identity(new)):
                            raise ValueError('Edit invalidates retained %s record for atom %d; '
                                             'update/remove the dependent record explicitly' % (kind, serial))
            last_slot = {}
            for index, record in enumerate(self._records):
                if record[0] == 'atom':
                    last_slot[record[1]] = index
            positions = dict((res, 0) for res in self.lstRes)
            new_residues = [res for res in self.lstRes if res not in original_residues]
            emitted_new = set()
            output = []
            def emit_atom(res, index, original=''):
                output.append(self._format_atom(res.lstAtoms[index], serials[(res, index)], original))
            def emit_new(model):
                pending = [res for res in new_residues if res.model_index == model and res not in emitted_new]
                if not pending:
                    return
                # Explicitly separate additions from the preceding polymer segment.
                if any(line[:6].strip() == 'ATOM' for line in output):
                    output.append('TER')
                previous = None
                for res in pending:
                    group = (res.chain_id, res.segment_id)
                    if previous is not None and group != previous:
                        output.append('TER')
                    for index in range(len(res.lstAtoms)):
                        emit_atom(res, index)
                    previous = group
                    emitted_new.add(res)
            atom_count_changed = (sum(len(res.lstAtoms) for res in self.lstRes) !=
                                  sum(record[0] == 'atom' for record in self._records))
            topology_changed = bool(new_residues or original_residues - current_residues or atom_count_changed)
            for index, record in enumerate(self._records):
                if record[0] == 'raw':
                    line, model = record[1:]
                    kind = line[:6].strip()
                    if kind == 'ENDMDL' or (not self._explicit_models and kind in ('CONECT', 'MASTER', 'END')):
                        emit_new(model if self._explicit_models else 0)
                    if kind == 'MASTER' and topology_changed:
                        # Omit the optional count summary only when edits invalidate it.
                        continue
                    if kind == 'TER' and len(line) >= 27:
                        preceding = next((item for item in reversed(output)
                                          if item[:6].strip() in ('ATOM', 'MODEL', 'ENDMDL', 'TER')), None)
                        if preceding and preceding[:6].strip() == 'ATOM':
                            line = line[:17] + preceding[17:27] + line[27:]
                    output.append(line)
                    continue
                res = record[1]
                if res not in current_residues:
                    continue
                slot = positions[res]
                if slot < len(res.lstAtoms):
                    emit_atom(res, slot, record[3])
                    positions[res] += 1
                if index == last_slot[res]:
                    while positions[res] < len(res.lstAtoms):
                        emit_atom(res, positions[res])
                        positions[res] += 1
            if not self._explicit_models:
                emit_new(0)
            if any(res not in emitted_new for res in new_residues):
                raise ValueError('Could not place a new residue in its model')
            self._data = output

        def SaveToFile(self, pdb_file):
            self.RebuildPDB()
            # Write completely before replacing the destination (including source=destination).
            destination = os.path.abspath(pdb_file)
            descriptor, temporary = tempfile.mkstemp(prefix='.artemis-', dir=os.path.dirname(destination))
            try:
                with os.fdopen(descriptor, 'w') as handle:
                    for line in self._data:
                        handle.write(line.rstrip('\r\n') + '\n')
                # POSIX rename replaces atomically; Windows Python 2 cannot replace
                # an existing target atomically, so fail rather than deleting it.
                os.rename(temporary, destination)
            finally:
                if os.path.exists(temporary):
                    os.remove(temporary)


    def PeptideLinked(left, right, left_atoms, right_atoms):
        if (left is None or right is None or
                left.model_index != right.model_index or
                left.segment_id != right.segment_id or left.chain_id != right.chain_id or
                left.res_type != RESTYPE.AA or right.res_type != RESTYPE.AA):
            return False
        carbon = left_atoms.get('C')
        nitrogen = right_atoms.get('N')
        if carbon is None or nitrogen is None:
            return False
        if carbon.alt_loc and nitrogen.alt_loc and carbon.alt_loc != nitrogen.alt_loc:
            return False
        distance = Distance3d(*(CoordsListFromAtomTuple(carbon) + CoordsListFromAtomTuple(nitrogen)))
        return 0.5 < distance <= 2.0


    def ResidueDihedrals(residues, index):
        """Compute each available torsion independently, including terminal residues.

        A 2-Angstrom C--N cutoff is a conservative connectivity heuristic.
        Undefined/missing/broken-chain angles are returned as None.
        """
        residue = residues[index]
        atoms = residue.SelectedAtoms()
        result = dict(Phi=None, Psi=None, Ohmega=None)
        if residue.res_type != RESTYPE.AA:
            return result
        def torsion(atom_list):
            if any(atom is None for atom in atom_list):
                return None
            try:
                return DihedralAngleHandler(*[CoordsListFromAtomTuple(atom) for atom in atom_list])
            except ValueError:
                return None
        if index > 0:
            previous = residues[index - 1]
            before = previous.SelectedAtoms()
            if PeptideLinked(previous, residue, before, atoms):
                result['Phi'] = torsion([before.get('C'), atoms.get('N'), atoms.get('CA'), atoms.get('C')])
        if index + 1 < len(residues):
            following = residues[index + 1]
            after = following.SelectedAtoms()
            if PeptideLinked(residue, following, atoms, after):
                result['Psi'] = torsion([atoms.get('N'), atoms.get('CA'), atoms.get('C'), after.get('N')])
                result['Ohmega'] = torsion([atoms.get('CA'), atoms.get('C'), after.get('N'), after.get('CA')])
        return result

    #// Main ==========================================================================
    p = TRealPoint(1,2,3);

    b_stdin = not sys.stdin.isatty();

    if ( (len(sys.argv) == 1 and (b_stdin == True)) or (len(sys.argv) == 2) ):

        if (len(sys.argv) == 1):
            PDBModel = TPDBModel();
            PDBModel.LoadFromStream();

        if (len(sys.argv) == 2):
            PDBModel = TPDBModel();
            #PDBModel.load("./default.pdb");
            PDBFile = sys.argv[1];
            PDBModel.LoadFromFile(PDBFile);

        print "Artemis - Python molecular library for Brookhaven PDB files";
        print "Copyright (c) 2014 Jamie J. Al-Nasir, All Rights Reserved";
        print "";

        if (len(sys.argv) == 2):
            print "Loaded: " + PDBFile;
        else:
            print "Loaded from <stdin> ";

        print "";

        #// Implementation goes here ==================================================

        # PDBModel.lstChains contains a list of chain_id
        # PDBModel.getResByChainID returns a list of residues for given chain_id
        # PDBModel.lstRes contains a list of all residues
        # PDBModel.lstAtomTuples contains a list of all atom tuples (NB tuples are immutable)
        # aRes[n].lstBonds contains a list of all TBond *Objects* in the residue
        #   a TBond contains two Atom *Objects* (not tuples)
        #
        # NB molecular geometry calculations can be used for a variety of computation tasks
        # such as calculating bond lengths and torsional/dihedral angles.




        # -----------------------------------------------------------------------------
        # Saving Example: Create a new Residue and assign it to a new chain D then save

        #aNewRes = TResidue();

        # Create the atoms for the new residue
        # aNewAtomTuple =             "serial,name,alt_loc,chain_id,res_name,res_seq,icode,   x,       y,        z,    occ,temp,element,charge");
        # NB: Existing serials are preserved; use "auto" for new atoms.
        #aNewAtomTupleN   = tplPDBAtom("auto",   "N",   "",  "D",      "TST",   "155",  "",  "-1.001", "+2.002", "-3.003", "", "",  "N",   "");
        #aNewAtomTupleCA  = tplPDBAtom("auto",   "CA",  "",  "D",      "TST",   "155",  "",  "-3.003", "+2.001", "-1.002", "", "",  "CA",  "");
        #aNewAtomTupleC   = tplPDBAtom("auto",   "C",   "",  "D",      "TST",   "155",  "",  "-4.003", "+2.001", "-3.002", "", "",  "C",   "");

        #aNewRes.addAtomTuple(aNewAtomTupleN);
        #aNewRes.addAtomTuple(aNewAtomTupleCA);
        #aNewRes.addAtomTuple(aNewAtomTupleC);
        #aNewRes.chain_id =  "D";
        #PDBModel.lstChains.append("D");
        #PDBModel.lstRes.append(aNewRes);

        # If necessary call PDBModel.Rebuild (SaveToFile also performs rebuild prior to saving)

        # Save the modified PDB structure
        #PDBModel.SaveToFile('test_out.pdb');
        #------------------------------------------------------------------------------


        # Report contiguous model/chain/segment runs independently.
        groups = []
        for residue in PDBModel.lstRes:
            key = (residue.model_index, residue.chain_id, residue.segment_id)
            if not groups or groups[-1][0] != key:
                groups.append((key, []))
            groups[-1][1].append(residue)
        for key, residues in groups:
            print "Model: %s; Chain: %s; Segment: %d" % (residues[0].model_id, key[1] or '(blank)', key[2] + 1)
            if _PRINT_RESIDUES_:
                for aRes in residues:
                    print "\nResidue: " + aRes.res_name + aRes.res_seq + aRes.icode
                    print "chain \t name \t residue \t x \t y \t z \t altLoc"
                    for atom in aRes.lstAtoms:
                        print "\t".join((atom.chain_id or '(blank)', atom.name,
                                         atom.res_name + atom.res_seq + atom.icode,
                                         atom.x, atom.y, atom.z, atom.alt_loc))
                    if _PRINT_RES_BONDS_:
                        aRes.ComputeResBonds()
                        print "\n%d Covalent bond(s) computed within this residue:" % len(aRes.lstBonds)
                        for bond in aRes.lstBonds:
                            print bond.getBondName()
            if _PRINT_DIHEDRALS_:
                print "\nDihedral/Torsional angles (N/A = undefined, missing atoms or chain break):"
                for index, residue in enumerate(residues):
                    if residue.res_type != RESTYPE.AA:
                        continue
                    angles = ResidueDihedrals(residues, index)
                    values = ["N/A" if angles[name] is None else "%.3f" % angles[name]
                              for name in ('Phi', 'Psi', 'Ohmega')]
                    print "Residue %s%s%s: Phi=%s, Psi=%s, Ohmega=%s" % (
                        residue.res_name, residue.res_seq, residue.icode,
                        values[0], values[1], values[2])

        # Final spacer line
        print "\n";

        #// EndImplementation =========================================================


    else:
        print "Artemis - Python molecular library for Brookhaven PDB files";
        print "Copyright (c) 2014 Jamie J. Al-Nasir, All Rights Reserved";
        print "";
        print "Build on the code, or provide a PDB file to process:";
        print "";
        print "Usage Artemis.py <file.pdb>";
        print "";


# Main program execution in try-except block to catch un-caught exceptions
# from elsewhere

if __name__ == '__main__':
    try:
        main()
    except Exception as ErrMsg:
        sys.stderr.write('An error occurred: %s\n' % ErrMsg);
        sys.exit(1);
