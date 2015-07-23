#!/usr/bin/python
#//==============================================================================
#// Artemis - Python molecular library for Brookhaven PDB files
#// By Jamie Al-Nasir, 07/2014 
#// Royal Holloway University of London
#// CSSB - Centre for Systems and Synthetic Biology, dept. Computer Science
#// Copyright (c) 2014 Jamie J. Al-Nasir, All Rights Reserved
#//==============================================================================
#// Version: Python edition
#//==============================================================================

import sys;
import math;

def main():
	# Define our own Enum (Python 3.4 has it's own enum type)
	class Enum(object): 
		def __init__(self, lstTuple):
		        self.lstTuple = lstTuple

		def __getattr__(self, name):
		        return self.lstTuple.index(name)

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
		# Return CosA of Angle between lines aLine1 and aLine2 originating from 0
		# of Lengths aLine1Len and aLine2Len respectively
		try:
			return aDP / (aLine1Len * aLine2Len);
		except ZeroDivisionError as errormsg:
			print 'An error occured: ', errormsg;
			return 0;


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

	def DihedralAngle(vA_x, vA_y, vA_z, vB_x, vB_y, vB_z, vC_x, vC_y, vC_z, vD_x, vD_y, vD_z):
	#
	#   A
	#	\
	#	 B----C
	#		   \
	#			D

		arr_dist21 = VectorSubtract(vC_x, vC_y, vC_z, vB_x, vB_y, vB_z);
		arr_dist01 = VectorSubtract(vA_x, vA_y, vA_z, vB_x, vB_y, vB_z);
		arr_dist32 = VectorSubtract(vD_x, vD_y, vD_z, vC_x, vC_y, vC_z);		
	
		dist21_x = arr_dist21[0]; dist21_y = arr_dist21[1]; dist21_z = arr_dist21[2];
		dist01_x = arr_dist01[0]; dist01_y = arr_dist01[1]; dist01_z = arr_dist01[2];
		dist32_x = arr_dist32[0]; dist32_y = arr_dist32[1]; dist32_z = arr_dist32[2];
	
		arr_dd1 = CrossProduct(dist21_x, dist21_y, dist21_z, dist01_x, dist01_y, dist01_z);
		arr_dd3 = CrossProduct(dist21_x, dist21_y, dist21_z, dist32_x, dist32_y, dist32_z);

		dd1_x = arr_dd1[0]; dd1_y = arr_dd1[1]; dd1_z = arr_dd1[2];
		dd3_x = arr_dd3[0]; dd3_y = arr_dd3[1]; dd3_z = arr_dd3[2];
	
		r = Angle(0,0,0, dd1_x, dd1_y, dd1_z, dd3_x, dd3_y, dd3_z);
		arr_pos_d  = CrossProduct(dist21_x, dist21_y, dist21_z, dd1_x, dd1_y, dd1_z);

		pos_d_x = arr_pos_d[0]; pos_d_y = arr_pos_d[1]; pos_d_z = arr_pos_d[2]; 
		if (DotProduct3d(dd3_x, dd3_y, dd3_z, pos_d_x, pos_d_y, pos_d_z) < 0):
			r = -r;
		return r;


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
			return [0,0,0]
	
	def CreateAtomFromTuple(aAtomTuple):
		if not aAtomTuple is None:
			anAtom = TAtom(aAtomTuple);
			return anAtom;
		else:
			return None;

	def IsDNARes(aResName):
		if aResName[0] in KSTR_DNA_CODES:
			return True;
		else:
			return False;

	def GetPDBResType(aResName):
		r = RESTYPE.AA;
		tmpstr = aResName.strip();
		if (len(tmpstr) == 1) and IsDNARes(tmpstr[0]):
			r = RESTYPE.DNA;
		else:
			if (len(tmpstr) == 2):  
				if IsDNARes(tmpstr[1]):
					r = RESTYPE.DNA; # case of DA1, DG1 etc..
				else:
					r = RESTYPE.AA;
		return r;

		if IsDNARes(aResName):
			return RESTYPE.DNA;
		else:
			return RESTYPE.AA;
		
	class TAtom:
		"Atom object class"
		__slots__ = ('bBackbone', 'bHydrogen', 'bComputed', 'Hybridisation', 'PDBData');
		def __init__(self, PDBAtomTuple):
			self.xyz = TRealPoint(0,0,0);
			self.bBackBone = 0;
			self.bHydrogen = 0;
			self.bComputed = 0;
			self.Hybridisation = HYBRIDISATION.HYBRID_Unknown;
			self.PDBData = PDBAtomTuple;
			if PDBAtomTuple is not None:
				self.xyz.x = PDBAtomTuple.x;
				self.xyz.y = PDBAtomTuple.y;
				self.xyz.z = PDBAtomTuple.z;
		
	class TBond:
		"Bond object class"
		def __init__(self, a1=None, a2=None):
			self.Atom1 = a1;
			self.Atom2 = a2;
			if not self.Atom1 is None:
				self.Atom1.Hybridisation = HYBRIDISATION.HYBRID_Unknown;
			if not self.Atom2 is None:
				self.Atom1.Hybridisation = HYBRIDISATION.HYBRID_Unknown;
			self.SciAtomBondType   = BOND_TYPE.BOND_Single;
			self.bBackBone = 0;
			self.bComputed = 0;				
	
		def getBondName(self):
			r = "";
			if self.Atom1 is None or self.Atom2 is None:
				return r;			
			if self.Atom1.PDBData is None or self.Atom2.PDBData is None:
				return "Error: Atom.PDBData is null" + str(type(self.Atom1.PDBData));
			r = self.Atom1.PDBData.name + " - " + self.Atom2.PDBData.name;
			return r;
		
	class TMolecule:
		"Molecule object base class"
		def __init__(self):
			self.lstBonds = [];		
		def addBond(self, aBond):
			if not (aBond.Atom1 is None or aBond.Atom2 is None):
				if aBond not in self.lstBonds:
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
		
		def addAtomTuple(self, PDBAtomTuple):
			self.lstAtoms.append(PDBAtomTuple);
		
		def ComputeResBonds(self):
			# Adapted from the Zeus Molecular visualisation framework developed by Jamie Al-Nasir
	
			# Amino Acids ===============================================================

			# Check a-Amino backbone
			# Add bonds for a-Amino backbone, N-C-C=O
			# Covers GLY			

			# [N-C]-C=O
			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CA"));
			Bond.bBackBone = 1;
			self.addBond(Bond);

			# N-[C-C]=0
			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CA"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C"));
			Bond.bBackBone = 1;
			self.addBond(Bond);

			# N-C-[C=0]
			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Double;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "O"));

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
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CA"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CB"));			
			if (Bond.Atom1 is not None and Bond.Atom2 is not None):
				Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CA is SP3
				Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CB is SP3
				self.addBond(Bond);
			else:
				del Bond;

			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CB"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));			
			if (Bond.Atom1 is not None and Bond.Atom2 is not None):
				Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CA is SP3
				Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CG is SP3
				self.addBond(Bond);
			else:
				del Bond;

			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD"));			
			if (Bond.Atom1 is not None and Bond.Atom2 is not None):
				Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CG is SP3
				Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CD is SP3
				self.addBond(Bond);
			else:
				del Bond;
			# End Alkyl backbone

			if (self.res_name == "ASP"):
				# Add COO- to D-carbon of Alkyl side-chain						
				# which one is the C=O double bond?? OD1 or OD2? use OD1 for now
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "OD1"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "OD2"));
				self.addBond(Bond);

			if (self.res_name == "GLU"):
			# Add COO- to E-carbon of Alkyl side-chain					
				# which one is the C=O double bond?? OD1 or OD2? use OD1 for now
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "OE1"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "OE2"));
				self.addBond(Bond);

			if (self.res_name == "ILE"):
			# Add B-Carbon to 2x G-Carbon bonds and one G-Carbon-D-Carbon bond					
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CB"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG1"));
				if not Bond.Atom1 is None:  
					Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CB SP3 hybridised
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CG1 SP3 hybridised
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CB"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG2"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CG2 SP3 hybridised
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD1"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CD1 SP3 hybridised
				self.addBond(Bond);

			if (self.res_name == "LEU"):
			# Add 2x G-Carbon to D-Carbon bonds
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD1"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CD1 SP3 hybridised
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD2"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CD2 SP3 hybridised
				self.addBond(Bond);


			if (self.res_name == "MET"):
			# Add G-Carbon to D-Sulphur and D-Sulphur to E-Carbon bonds					
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "SD"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "SD"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # CE SP3 Hybridised
				self.addBond(Bond);		

			if (self.res_name == "VAL"):
			# Add 2x B-Carbon to G-Carbon bonds
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CB"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG1"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CB"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG2"));
				self.addBond(Bond);

			if (self.res_name == "LYS"):
			# Add D-Carbon to E-Carbon bond and E-Carbon to Z-Nitrogen (NH3) bond
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE"));
				if not Bond.Atom1 is None:  
					Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;
				self.addBond(Bond);
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "NZ"));
				self.addBond(Bond);


			if (self.res_name == "SER"):
			# Add B-Carbon to G-Oxygen (OH) bond
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CB"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "OG"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # OG SP3 Hybridised (4 SP3 orbitals, 2 lone pairs, so bent)
				self.addBond(Bond);

			if (self.res_name == "THR"):
			# Add B-Carbon to G-Oxygen (OH) and B-Carbon to G-Carbon bonds
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CB"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "OG1"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CB"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG2"));
				self.addBond(Bond);

			if (self.res_name == "ASN"):
			# Add G-Carbon to D-Oxygen (C=O) bond
			# Add G-Carbon to D-Nitrogen (NH3) bond
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "OD1"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "ND2"));
				self.addBond(Bond);

			if (self.res_name == "GLN"):
			# Add D-Carbon to E-Oxygen (C=O) bond
			# Add D-Carbon to E-Nitrogen (NH3) bond
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "OE1"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "NE2"));
				self.addBond(Bond);

			if (self.res_name == "ARG"):
			# Add D-Carbon to E-Nitrogen bond
			# Add E-Nitrogen to Z-Carbon bond
			# Add Z-Carbon double bond to H-N (NH2+)
			# Add Z-Carbon single bond to H-N (NH2)
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "NE"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "NE"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CZ"));
				self.addBond(Bond);

				# which one is the C=N double bond?? NH1 or NH2? use NH1 for now
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CZ"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "NH1"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CZ"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "NH2"));
				self.addBond(Bond);

			if (self.res_name == "CYS"):
			# Add B-Carbon to G-Sulphur bond
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CB"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "SG"));
				self.addBond(Bond);

			if (self.res_name == "PRO"):
			# Add D-Carbon to Nitrogen bond which closes the ring in Proline
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N"));
				self.addBond(Bond);

			if (self.res_name == "PHE"):
			# Create Phenyl Ring from two parallel alkyl-chains which meet up at Z-Carbon
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD1"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD2"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE1"));
				if not Bond.Atom1 is None:  
					Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE2"));
				if not Bond.Atom1 is None:  
					Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CZ"));
				if not Bond.Atom1 is None:  
					Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2; 
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2; 
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CZ"));
				self.addBond(Bond);
		
			if (self.res_name == "TYR"):
				# Create Phenyl Ring from two parallel alkyl-chains which meet up at Z-Carbon
				# Add bond to Z-Carbon to H-O (Phenolic-OH)			
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD1"));
				if not Bond.Atom1 is None:  
					Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CD1 SP2 hybridised
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CD1 SP2 hybridised
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD2"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CD2 SP2 hybridised
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE1"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CE1 SP2 hybridised
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE2"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CE2 SP2 hybridised
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CZ"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CZ"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CZ"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "OH"));
				self.addBond(Bond);

			if (self.res_name == "HIS"):
				# Create Hetero Ring from two parallel alkyl-chains which meet up at Z-Carbon
				# Add bond to Z-Carbon to H-O (Phenolic-OH)
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD2"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CD2 SP2 hybridised
				self.addBond(Bond);


				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "ND1"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CD2 SP2 hybridised/protonated, + charge
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "ND1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE1"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CE1 SP2 hybridised
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "NE2"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "NE2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE1"));
				self.addBond(Bond);
				
			if (self.res_name == "TRP"):
				# Create Phenyl Ring from two parallel alkyl-chains which meet up at Z-Carbon
				# Add bond to Z-Carbon to H-O (Phenolic-OH)			
			
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD1"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CD1 SP2 hybridised
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CG"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD2"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "NE1"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "NE1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE2"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD2"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CZ2"));
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CZ2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CH2"));
				if not Bond.Atom1 is None:  
					Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CZ2 SP2 hybridised
				self.addBond(Bond);

				#
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CH2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CZ3"));
				if not Bond.Atom1 is None:  
					Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CH2 SP2 hybridised
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CZ3"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE3"));
				if not Bond.Atom1 is None:  
					Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CZ3 SP2 hybridised
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP2;  # CE3 SP2 hybridised
				self.addBond(Bond);

				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CE3"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "CD2"));
				self.addBond(Bond);

			# Nucleic Acids =============================================================

			if (self.res_name == "C"):
			# Cytosine
				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C6"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N3"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N4"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N3"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "O2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N1"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C6"));
				Bond.bBackBone = False;
				self.addBond(Bond);

			if (self.res_name == "T"):
			# Thymine
		
				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "O2"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N3"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N3"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "O4"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C6"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C6"));
				Bond.bBackBone = False;
				self.addBond(Bond);

			if (self.res_name == "A"):
			# Adenine

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N9"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C8"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C8"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N7"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N7"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C6"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C6"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N6"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C6"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N1"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N3"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N9"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N3"));
				Bond.bBackBone = False;
				self.addBond(Bond);

			if (self.res_name == "G"):
			# Guanine

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N9"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C8"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C8"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N7"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N7"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C6"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C6"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "O6"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C6"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N1"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N3"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N9"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N3"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N2"));
				Bond.bBackBone = False;
				self.addBond(Bond);

			if (self.res_name == "I"):
			# Inosine

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N9"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C8"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C8"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N7"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N7"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C6"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C6"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "O6"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C6"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N1"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N3"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Double;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N9"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N3"));
				Bond.bBackBone = False;
				self.addBond(Bond);

				# [N-C]-C=O
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N2"));
				Bond.bBackBone = False;
				self.addBond(Bond);

			if (self.res_name == "C" or self.res_name == "G" or self.res_name ==  "I"):
			# Attach Purine to Sugar Phosphate
		
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N9"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C1*"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # C1* SP3 Hybridised
				Bond.bBackBone = False;
				self.addBond(Bond);

			if (self.res_name == "C" or self.res_name == "t"):
			# Attach Pyramidine to Sugar Phosphate
				Bond = TBond();
				Bond.BondType = BOND_TYPE.BOND_Single;
				Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "N1"));
				Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C1*"));
				if not Bond.Atom2 is None:  
					Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # C1* SP3 Hybridised
				Bond.bBackBone = False;
				self.addBond(Bond);

			# Phosphate Backbone
			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "O1P"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "P"));
			Bond.bBackBone = 1;
			self.addBond(Bond);

			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Double;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "O2P"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "P"));
			Bond.bBackBone = 1;
			self.addBond(Bond);


			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C1*"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2*"));
			if not Bond.Atom1 is None:  
				Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # C1* SP3 Hybridised
			if not Bond.Atom2 is None:  
				Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # C2* SP3 Hybridised
			Bond.bBackBone = 1;
			self.addBond(Bond);

			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C2*"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C3*"));
			if not Bond.Atom2 is None:  
				Bond.Atom2.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # C1* SP3 Hybridised
			Bond.bBackBone = 1;
			self.addBond(Bond);

			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C3*"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4*"));
			Bond.bBackBone = 1;
			self.addBond(Bond);

			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4*"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "O4*"));
			Bond.bBackBone = 1;
			self.addBond(Bond);

			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "O4*"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C1*"));
			Bond.bBackBone = 1;
			self.addBond(Bond);

			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C4*"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5*"));
			Bond.bBackBone = 1;
			self.addBond(Bond);

			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C3*"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "O3*"));
			Bond.bBackBone = 1;
			self.addBond(Bond);

			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "C5*"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "O5*"));
			if not Bond.Atom1 is None:  
				Bond.Atom1.Hybridisation = HYBRIDISATION.HYBRID_SP3;  # C5* SP3 Hybridised
			Bond.bBackBone = 1;
			self.addBond(Bond);

			Bond = TBond();
			Bond.BondType = BOND_TYPE.BOND_Single;
			Bond.Atom1 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "P"));
			Bond.Atom2 = CreateAtomFromTuple(ListTupleByAttrVal(self.lstAtoms, 1, "O5*"));
			Bond.bBackBone = 1;
			self.addBond(Bond);


		
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

	class TPDBModel:    
		def __init__(self):
			self._data = [];
			self._pdb_file = "";
			self.lstAtomTuples = [];   # array to hold immutable atom tuples
			self.lstChains = [];	   # array to hold list of chain ids strings
			self.lstRes = []; 		   # array to hold residue objects

		def _load(self):	
			cur_res = "";
			cur_chain = "";
			aRes = TResidue();
			for line in self._data:
				record_type = line[:6].strip();
				atom_serial = line[7:11].strip();
				atom_name   = line[13:16].strip();
				alt_loc     = ""#line[17].strip();
				res_name    = line[17:20].strip();
				chain_id    = line[20:22].strip();
				res_seq     = line[23:26].strip();
				icode       = ""#line[27].strip();
				x           = line[31:38].strip();
				y           = line[39:46].strip();
				z           = line[47:54].strip();
				occupancy   = line[55:60].strip();
				temp_factor = line[61:66].strip();
				element     = line[77:78].strip();
				charge      = line[79:80].strip();
				#if (res_name <> res):

				if (record_type == "ATOM"):
					anAtomTuple = tplPDBAtom(atom_serial, atom_name, alt_loc, chain_id, res_name, res_seq, icode, x, y, z, occupancy, temp_factor, element, charge);
										
					#print record_type, atom_serial, chain_id, res_name, res_seq, x, y, z;
					self.lstAtomTuples.append(anAtomTuple);	
					if (res_name + res_seq <> cur_res):
						if not aRes is None:
							self.lstRes.append(aRes);
						del aRes;
						aRes = TResidue();					
						aRes.res_name = res_name;
						aRes.res_seq = res_seq;
						cur_res = res_name + res_seq;
						aRes.res_type = GetPDBResType(res_name);
						#print id(aRes);
					aRes.addAtomTuple(anAtomTuple);
					aRes.chain_id =  anAtomTuple.chain_id;
					#print anAtomTuple;
					if (chain_id <> cur_chain):
						if chain_id not in self.lstChains:
							self.lstChains.append(chain_id);
						cur_chain = chain_id;
			# last residue
			if not aRes is None:
				self.lstRes.append(aRes);
			
		def getResByChainID(self, chain_id):
			result = [];
			for aRes in PDBModel.lstRes:
				if (aRes.chain_id == chain_id):
					result.append(aRes);
			return result;

	
		def LoadFromFile(self, pdb_file):
			self._pdb_file = pdb_file;
			f = open(self._pdb_file, 'r');
			self._data = f.readlines();
			self._load();

		def LoadFromStream(self):
			self._data = sys.stdin;
			self._load();

		def RebuildPDB(self):
			print "Rebuilding PDB structure...";
			self._data = [];
			serial = 0;
			
			# Copy non ATOM records/lines from before/after the block of ATOM record lines
			# so as to preserve user/program encoded data from the original PDB
			f = open(self._pdb_file, 'r');
			_tmpdata = f.readlines();
			
			lstHead = [];
			lstTail = [];
			bHead = True;
			
			for aHeaderLine in _tmpdata:
				if aHeaderLine.startswith('ATOM  '):
					bHead = False;
					continue;
				else:
					if bHead:
						lstHead.append(aHeaderLine);
					else:
						lstTail.append(aHeaderLine);
						
			
			for aChain in PDBModel.lstChains:
                        	for aRes in PDBModel.getResByChainID(aChain):
                                	for anAtom in aRes.lstAtoms:
											serial = serial + 1;
                                        	#print aChain + "\t" + anAtom.name + "\t" + anAtom.res_name + anAtom.res_seq + "\t" + anAtom.x + "\t" + anAtom.y + "\t" + anAtom.z + "\t";
											pdb_line = '{:6s}{:>5s}  {:4s}{:>3s} {:1s}{:>4s}    {:>8s}{:>8s}{:>8s}'.format("ATOM",str(serial),anAtom.name,anAtom.res_name,anAtom.chain_id,anAtom.res_seq, anAtom.x, anAtom.y, anAtom.z);
											self._data.append(pdb_line);
											
			# Let _data now contain Header records + _data (Atom records) + Tail records 
			self._data = lstHead + self._data + lstTail;
			

		def SaveToFile(self, pdb_file):
			self.RebuildPDB();
			print "Saving " + pdb_file;
			f = open(pdb_file, 'w');
			for line in self._data:
			  f.write("%s\n" % line);
		   

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
		# PDBModel.lstAtoms contains a list of all atom tuples (NB tuples are immutable)
		# aRes[n].lstBonds contains a list of all TBond *Objects* in the residue
		# 	a TBond contains two Atom *Objects* (not tuples)
		#
		# NB molecular geometry calculations can be used for a variety of computation tasks
		# such as calculating bond lengths and torsional/dihedral angles.
	
		# Default example is to produce a report of the atoms in the PDB model by residue
		# and list the bonds computed/found within those residues.

		# Activities to perform:
		bBonds     = True;
		bDihedrals = True;
		
		
		# -----------------------------------------------------------------------------
		# Saving Example: Create a new Residue and assign it to a new chain D then save
		
		#aNewRes = TResidue();
		
		# Create the atoms for the new residue
		# aNewAtomTuple =             "serial,name,alt_loc,chain_id,res_name,res_seq,icode,   x,       y,        z,    occ,temp,element,charge");
		# NB: ATOM serial is auto generated during rebuild/save
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
		

		for aChain in PDBModel.lstChains:
			print "Chain: " + aChain;
			for aRes in PDBModel.getResByChainID(aChain):			
				print "\n";
				print "Residue: " + aRes.res_name + aRes.res_seq;
				aRes.ComputeResBonds();
				print "chain \t name \t residue x \t y \t z";
				for anAtom in aRes.lstAtoms:
					print aChain + "\t" + anAtom.name + "\t" + anAtom.res_name + anAtom.res_seq + "\t" + anAtom.x + "\t" + anAtom.y + "\t" + anAtom.z + "\t";
				if bBonds:
					print "\n";
					print str(len(aRes.lstBonds)) + " Covalent bond(s) computed within this Residue:";
					for aBond in aRes.lstBonds:
						print aBond.getBondName();

			if bDihedrals:
				print "\n";
				print "Dihedral/Torsional angles for residues in chain " + aChain + ":";
				aRes = PDBModel.getResByChainID(aChain);			
				x = 1;
				while (x < len(aRes) - 1):
					x = x + 1;
				
					# We don't do DNA!
					if (aRes[x - 2].res_type == RESTYPE.DNA) or (aRes[x - 1].res_type == RESTYPE.DNA) or (aRes[x].res_type == RESTYPE.DNA):
						continue;
				
					#print aRes[x - 2].res_name + " " + aRes[x - 2].res_seq + str(CoordsListFromAtomTuple( ListTupleByAttrVal(aRes[x - 2].lstAtoms, 1, "C") ));
								

					res_atoms = aRes[x - 2].lstAtoms;				
					tplN = ListTupleByAttrVal(res_atoms, 1, "N");
					tplCA = ListTupleByAttrVal(res_atoms, 1, "CA");
					tplC = ListTupleByAttrVal(res_atoms, 1, "C");				
					# We need a C from residue 1
					if (tplC is None):
						continue;				
					r1_N = CoordsListFromAtomTuple( tplN );
					r1_CA = CoordsListFromAtomTuple( tplCA );
					r1_C = CoordsListFromAtomTuple( tplC );
					del res_atoms, tplN, tplCA, tplC;

					res_atoms = aRes[x - 1].lstAtoms;
					res_id = aRes[x - 1].res_name + aRes[x - 1].res_seq; # used for labelling
					tplN = ListTupleByAttrVal(res_atoms, 1, "N");
					tplCA = ListTupleByAttrVal(res_atoms, 1, "CA");
					tplC = ListTupleByAttrVal(res_atoms, 1, "C");
					# We need N, CA, C from residue 2
					if ((tplN is None) or (tplCA is None) or (tplC is None)):
						continue;
					r2_N = CoordsListFromAtomTuple( tplN );
					r2_CA = CoordsListFromAtomTuple( tplCA );
					r2_C = CoordsListFromAtomTuple( tplC );
					del res_atoms, tplN, tplCA, tplC;

					res_atoms = aRes[x].lstAtoms;
					tplN = ListTupleByAttrVal(res_atoms, 1, "N");
					tplCA = ListTupleByAttrVal(res_atoms, 1, "CA");
					tplC = ListTupleByAttrVal(res_atoms, 1, "C");
					# We need N, CA from residue 3
					if ((tplN is None) or (tplCA is None)):
						continue;				
					r3_N = CoordsListFromAtomTuple( tplN );
					r3_CA = CoordsListFromAtomTuple( tplCA );
					r3_C = CoordsListFromAtomTuple( tplC );
					del res_atoms, tplN, tplCA, tplC;

					# Phi angle
					# x,y,z of First residue: C
					# x,y,z of Second residue: N
					# x,y,z of Second residue: Ca 
					# x,y,z of Second residue: C
					# COMMENTED OUT UNTIL CAN EXPAND r1_N etc... to pass in x,y,z coordinates
					# AS DihedralAngle takes four sets of x,y,z coordinates
					Phi = DihedralAngleHandler(r1_C, r2_N, r2_CA, r2_C);

	
					# Psi
					# x,y,z of Second residue: N
					# x,y,z of Second residue: CA
					# x,y,z of Second residue: C
					# x,y,z of Third residue: N
					Psi = DihedralAngleHandler(r2_N, r2_CA, r2_C, r3_N);

					# Ohmega
					# x,y,z of Second residue: CA
					# x,y,z of Second residue: C
					# x,y,z of Third residue: N
					# x,y,z of Third residue: CA
					Ohmega = DihedralAngleHandler(r2_CA, r2_C, r3_N, r3_CA);
			
					#print "\t".join(["angles", res_id, str(Phi), str(Psi), str(Ohmega)]);
					print "Residue %s: Phi=%.3f, Psi=%.3f, Ohmega=%.3f" % (res_id, Phi, Psi, Ohmega);


		# Final spacer line
		print "\n";
				
		#// EndImplementation =========================================================
				
							
	else:
		print "Artemis - Python molecular library for Brookhaven PDB files";	
		print "Copyright (c) 2014 Jamie J. Al-Nasir, All Rights Reserved";
		print "";
		print "Usage Artmis.py <file.pdb>";
		print "";


# Main program execution in try-except block to catch un-caught exceptions
# from elsewhere

if __name__ == '__main__':
    try:
        main()
    except Exception as ErrMsg:
        print 'An error occured: ', ErrMsg;
        sys.exit(0);

