# Artemis
Artemis is a Python library and standalone executable for computing Dihedral (Torsional) angles and reading/writing PDB-files. It allows users to incorporate calculation of dihedral angles in their own code and applications. Artemis also builds data models of protein and nucleic moeities in PDB macromolecules and computes bonds between atoms of these residues. The models and their constituent records can be easily accessed and iterated -- this library facilitates both reading and writing of Brookhaven PDB files. 

Artemis has been used in benchmarking and demonstrating the author's PDB-Hadoop framework at the 3D-Sig conference (part of ISMB -- Institute of Molecular Biology) in Dublin 2015 and has been used by MSc Datascience students in their final year project. Artemis is freely available via the Github repository.


# Dihedral Angles

Dihedral angles (also known as Torsional angles) in peptides are measured between select atoms in neighbouring residues and yield important information about the structural conformation of the residues in the peptide (termed secondary structure). Three such dihedral angles that can be measured: φ, ψ, and ω. The φ (Phi) angle is measured from the C of one residue to the C of the next residue. The ψ (Psi) angle is measured from the N of one residue to the N of the next. The ω (Ohmega) angle is measured from the C) of one residue to the C) of the next.


