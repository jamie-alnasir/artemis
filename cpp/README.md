# Artemis C++

A portable C++11 port of `Artemis3.py`, preserving its command-line options,
PDB parsing, residue templates, alternate-conformer selection and signed
backbone dihedral calculations.

The executable uses solely the C++ standard library. There is no dependency
on Python, Boost, Eigen, a package manager, or an external chemistry library.

## Compile

With GCC or Clang on Linux/macOS, from this folder:

```sh
c++ -std=c++11 -O2 -Wall -Wextra -pedantic artemis.cpp -o artemis
./artemis --help
```

Alternatively, with GNU Make:

```sh
make
make test
```

To choose the compiler explicitly:

```sh
make CXX=clang++
```

Run `make clean` before changing compilers or target architectures. The build
uses no processor-specific optimisation flags. Avoid `-ffast-math`: geometry
validation depends on normal finite/NaN floating-point behaviour.

### CMake (optional)

CMake is a convenience, not a dependency of the program. With CMake 3.10 or later:

```sh
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
ctest -C Release --output-on-failure
```

Use `-DBUILD_TESTING=OFF` to build only the executable. For cross-compilation,
pass your target's CMake toolchain file when configuring.

### Other architectures

Build natively with a C++11-capable compiler on the target machine, or select a
cross-compiler and its matching target sysroot. Examples, assuming the named
toolchain is installed:

| Target | Example command |
| --- | --- |
| Linux ARM64 | `aarch64-linux-gnu-g++ -std=c++11 -O2 artemis.cpp -o artemis` |
| Linux ARM hard-float | `arm-linux-gnueabihf-g++ -std=c++11 -O2 artemis.cpp -o artemis` |
| Linux RISC-V 64-bit | `riscv64-linux-gnu-g++ -std=c++11 -O2 artemis.cpp -o artemis` |
| Linux little-endian MIPS | `mipsel-linux-gnu-g++ -std=c++11 -O2 artemis.cpp -o artemis` |
| Linux x86 32-bit | `g++ -m32 -std=c++11 -O2 artemis.cpp -o artemis` |

The 32-bit x86 command needs multilib headers and libraries. Cross-compilers
must match the destination ABI, operating system and C++ runtime; the executable
is not automatically compatible with every device of the same CPU family.
There is no pointer-width, byte-order, SIMD or operating-system-specific code.
Static linking is optional and depends on the target toolchain.

On macOS, use the native compiler for the current architecture, or, with an
appropriate Apple SDK, build a universal binary:

```sh
clang++ -std=c++11 -O2 -arch x86_64 -arch arm64 artemis.cpp -o artemis
```

On Windows, from a Visual Studio Developer Command Prompt:

```bat
cl /EHsc /O2 /W4 /std:c++14 artemis.cpp /Fe:artemis.exe
artemis.exe --help
```

MSVC uses `/std:c++14` because it has no `/std:c++11` mode. MinGW-w64 can use the
same `g++ -std=c++11` command as Linux, with `-o artemis.exe`. The CMake build is
also suitable for Windows.

## Usage

```sh
# All reports
./artemis structure.pdb

# Dihedrals only for chain A
./artemis structure.pdb --chain A --no-residues --no-bonds

# Atoms and bonds for PDB MODEL 2
./artemis structure.pdb --model 2 --no-dihedrals

# Read standard input
cat structure.pdb | ./artemis - --no-bonds

# Save the full rebuilt structure without reports
./artemis structure.pdb -q -o rebuilt.pdb

# Write only PDB data to standard output
./artemis structure.pdb -q -o - > rebuilt.pdb

# Blank chain identifier (POSIX shells)
./artemis structure.pdb --chain ''
```

| Option | Behaviour |
| --- | --- |
| `--no-residues` | Suppress the atom/residue listing. Bond reports are independent. |
| `--no-bonds` | Skip bond calculations and reports. |
| `--no-dihedrals` | Skip dihedral calculations and reports. |
| `--chain ID` | Report only this case-sensitive chain; an empty value selects the blank chain. |
| `--model N` | Report only this positive PDB `MODEL` identifier, not its ordinal position. Files without `MODEL` records use `1`. |
| `-o PDB`, `--output PDB` | Save the full rebuilt structure. `-` writes PDB data to stdout and reports to stderr. |
| `-q`, `--quiet` | Suppress all reports and skip their calculations. Errors still go to stderr. |
| `-h`, `--help` | Display usage and examples. |
| `--` | End option parsing, allowing an input filename beginning with `-`. |

Options can precede or follow the input filename. Long value options also accept
`--chain=A`, `--model=2` and `--output=rebuilt.pdb`. Use full option names and
separate short options; Python argparse's abbreviations and combined short-option
forms are not implemented.

Omitting the input filename reads stdin until EOF, including in an interactive
terminal. Use `--help` to display help without reading input.

Exit codes: `0` success, `1` processing/I/O error, `2` invalid arguments or no
matching residues for the requested chain/model.

## Analysis and saving behaviour

- Chain/model selection filters **reports only**. Saved files retain the full
  structure, all alternate conformers and preserved non-`ATOM` records.
- `HETATM` records are preserved on saving but are not analysed.
- Residues are separated by model, chain, sequence number, insertion code and
  `TER` boundaries. Blank chain identifiers are supported.
- Analysis selects one residue conformer by highest mean occupancy, preferring
  `A` then lexical order on ties, together with shared blank-altLoc atoms.
- Bond templates match Artemis3's recognised amino acids and nucleotides,
  including legacy star/primes and phosphate atom-name aliases. They describe
  intra-residue bonds, not arbitrary ligands or inter-residue connectivity.
- Each available phi, psi and omega angle is calculated independently. Missing
  atoms, degenerate geometry and chain breaks produce `N/A`. Peptide continuity
  uses Artemis3's C--N distance heuristic of greater than 0.5 and at most 2.0 A.
- Atom serials and non-atom record order are preserved. Rebuilt numeric fields
  are formatted to PDB widths. Overflow and non-finite values are rejected.
- Library edits that invalidate retained serial-linked metadata such as `CONECT`
  or `ANISOU` are rejected on saving. Optional `MASTER` count records are omitted
  when structural edits invalidate their counts.
- Saving validates and builds the complete output before opening the destination.
  Unlike the Python version's atomic replacement, the standard C++11 file writer
  then writes directly to the destination. A write failure can leave a partial
  output; use a separate destination when preserving the source matters.

The report layout follows Artemis3, apart from a C++ version banner. The original
`Ohmega` output label is retained for compatibility.

## Use in another C++ program

Define `ARTEMIS_NO_MAIN` and include the source in one translation unit:

```cpp
#define ARTEMIS_NO_MAIN
#include "artemis.cpp"

int main() {
    artemis::Model model;
    model.load_file("structure.pdb");
    for (const auto& residue : model.residues) {
        residue->compute_bonds();
        std::cout << residue->res_name << residue->res_seq
                  << ": " << residue->bonds.size() << " bonds\n";
    }
    model.save("rebuilt.pdb");
}
```

Compile this example alone; do not separately link another copy of
`artemis.cpp`. The interface uses C++ containers and is not a drop-in Python API.

Useful types and functions in namespace `artemis`:

- `Model`: `load(std::istream&)`, `load_file(path)`, `residues_by_chain(...)`,
  `atom_count()`, `rebuild()` and `save(path)`.
- `Residue`: `atoms`, `bonds`, `selected_atoms()` and `compute_bonds()`.
- `Atom`: parsed PDB fields and `coords()`.
- `Bond`: indices `atom1` and `atom2` into the owning residue's atom vector,
  plus `order` and `backbone`.
- `Vec3`, `distance`, `angle`, `dihedral`, `residue_dihedrals` and `Torsions`.
  Geometry functions return degrees; undefined residue torsions use NaN.
- Familiar aliases: `TPDBModel`, `TResidue`, `TAtom` and `TBond`.

Models own residues through `std::shared_ptr`. Recompute bonds after editing an
atom vector. Atom field values are stored as strings to preserve PDB field text;
geometry methods convert them to doubles. Add atoms with `serial = "auto"` for
serial allocation during rebuilding. New residues must identify an existing
`model_index`; existing record slots determine the saved order. Preserved
records are exposed as `Model::records` for deliberate metadata editing.

## Tests and validation

Run the self-contained C++ regression suite:

```sh
make test
```

Optionally compare the CLI against the Python implementation at the repository
root (Python is needed only for this comparison):

```sh
python3 tests/compare_python.py --binary ./artemis --python ../Artemis3.py
```

Validation performed for this port:

- GCC 13.3 on Linux x86-64, C++11, with `-Wall -Wextra -Wpedantic -Wconversion
  -Wshadow -Werror`.
- Successful Makefile builds and C++ regression checks covering parsing,
  conformers, bond identity/order, geometry, boundaries, reloads and rebuilding.
- AddressSanitizer and UndefinedBehaviorSanitizer regression run passed. Leak
  checking was disabled because the execution environment blocks its process
  inspection.
- 160 CLI comparisons with Artemis3 passed, including exact successful report
  and rebuilt-PDB comparisons after normalising the version banner, plus error
  exit-code comparisons.

Other architectures, Windows/macOS compilers and the optional CMake build were
not executed in this environment. Their commands are build instructions, not
claims of cross-platform validation.

Original work: Jamie J. Alnasir, Royal Holloway University of London, CSSB, 2014.
The original copyright notice is preserved in the source.
