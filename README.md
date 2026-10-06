# Artemis
Artemis is a Python library and standalone executable for computing Dihedral (Torsional) angles and reading/writing PDB-files. It allows users to incorporate calculation of dihedral angles in their own code and applications. Artemis also builds data models of protein and nucleic moeities in PDB macromolecules and computes bonds between atoms of these residues. The models and their constituent records can be easily accessed and iterated -- this library facilitates both reading and writing of Brookhaven PDB files. 

Artemis has been used in benchmarking and demonstrating the author's PDB-Hadoop framework at the 3D-Sig conference (part of ISMB -- Institute of Molecular Biology) in Dublin 2015 and has been used by MSc Datascience students in their final year project. Artemis is freely available via the Github repository.


# Python 3 version and command-line options

`Artemis3.py` provides a Python 3 version of Artemis with command-line options for controlling reports, selecting chains and models, and saving rebuilt PDB files. It requires **Python 3.8 or later** and uses only the Python standard library, with no external dependencies. The Python 2.7 version remains available as `Artemis.py`.

By default, Artemis reports atoms by residue, computes intra-residue bonds, and calculates available backbone dihedral angles:

```bash
python3 Artemis3.py structure.pdb
```

## Available options

| Option | Description |
| --- | --- |
| `--no-residues` | Suppress the atom/residue listing; bond reports remain independently available. |
| `--no-bonds` | Skip bond calculations and bond reports. |
| `--no-dihedrals` | Skip dihedral calculations and reports. |
| `--chain ID` | Report only the specified chain. Chain identifiers are case-sensitive; use `--chain ''` for a blank identifier. |
| `--model N` | Report only the specified PDB `MODEL` identifier. Files without `MODEL` records use model `1`. |
| `-o PDB`, `--output PDB` | Save the full rebuilt structure to a PDB file. Use `-o -` to write PDB data to standard output and reports to standard error. |
| `-q`, `--quiet` | Suppress all reports and skip their calculations; errors remain visible. |
| `-h`, `--help` | Display usage information and examples. |

## Examples

Report dihedral angles only for chain A:

```bash
python3 Artemis3.py structure.pdb --chain A --no-residues --no-bonds
```

Report atoms and bonds for model 2, without calculating dihedral angles:

```bash
python3 Artemis3.py structure.pdb --model 2 --no-dihedrals
```

Read a PDB file from standard input:

```bash
cat structure.pdb | python3 Artemis3.py - --no-bonds
```

Save a rebuilt PDB without printing reports:

```bash
python3 Artemis3.py structure.pdb -q -o rebuilt.pdb
```

Display the complete command-line help:

```bash
python3 Artemis3.py --help
```

**Chain and model filters apply to reports only.** Saved PDB files retain the full structure, including all alternate conformers and preserved non-`ATOM` records. `HETATM` records are preserved when saving but are excluded from analysis. Undefined dihedral angles, including those affected by missing atoms or chain breaks, are displayed as `N/A`. Bond calculations use standard residue templates rather than general chemical bond perception.

The script returns exit status `0` on success, `1` for processing or I/O errors, and `2` for invalid command-line arguments.


# Dihedral Angles

Dihedral angles (also known as Torsional angles) in peptides are measured between select atoms in neighbouring residues and yield important information about the structural conformation of the residues in the peptide (termed secondary structure). Three such dihedral angles that can be measured: φ, ψ, and ω. The φ (Phi) angle is measured from the C of one residue to the C of the next residue. The ψ (Psi) angle is measured from the N of one residue to the N of the next. The ω (Ohmega) angle is measured from the C) of one residue to the C) of the next.
