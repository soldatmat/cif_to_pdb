# original script: 2023 Jean-Marie Bourhis
#  adapted script: 2026 Matouš Soldát

import argparse
import os
import sys

from tqdm import tqdm

from Bio.PDB import PDBIO, Select
from Bio.PDB.MMCIFParser import MMCIFParser


class _ResidueSelect(Select):
    """Select which residues to write.

    - standard (protein/nucleic) residues -> ATOM records (always kept)
    - ligands / ions / cofactors          -> HETATM records (kept when keep_hetero)
    - water                               -> always dropped

    Biopython's hetero flag is residue.id[0]: " " for standard residues, "W" for
    water, and "H_<resname>" for other hetero groups (ligands, ions, cofactors).
    """

    def __init__(self, keep_hetero: bool = True):
        self.keep_hetero = keep_hetero

    def accept_residue(self, residue):
        hetflag = residue.id[0]
        if hetflag == " ":
            return True            # protein / nucleic -> ATOM
        if hetflag == "W":
            return False           # water -> drop
        return self.keep_hetero    # ligand / ion / cofactor -> HETATM


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert an mmCIF to PDB (Biopython), preserving the B-factor "
        "(AlphaFold pLDDT). Ligands/ions/cofactors are kept as HETATM by default."
    )
    parser.add_argument("--input_cif", required=True, type=str)
    parser.add_argument("--output_pdb", required=True, type=str)
    parser.add_argument(
        "--no_hetero", dest="keep_hetero", action="store_false",
        help="Write protein/nucleic atoms only; drop ligands/ions/cofactors.",
    )
    parser.set_defaults(keep_hetero=True)
    return parser.parse_args()


def main(args):
    # Convert CIF -> PDB with Biopython (NOT Open Babel). Open Babel did not carry
    # the mmCIF B_iso_or_equiv column into the PDB B-factor (wrote 0.00), discarding
    # AlphaFold pLDDT, and it emitted ligands as (mislabeled) ATOM records. Biopython
    # MMCIFParser reads B_iso into atom.bfactor and PDBIO writes it back as the
    # B-factor (pLDDT 0-100), and writes ligands/ions correctly as HETATM. By default
    # the full complex is kept (protein ATOM + ligand/ion/cofactor HETATM); callers
    # that want a protein-only structure pass keep_hetero=False (CLI: --no_hetero).
    # Water is always dropped.
    keep_hetero = getattr(args, "keep_hetero", True)
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("s", args.input_cif)
    io = PDBIO()
    io.set_structure(structure)
    io.save(args.output_pdb, _ResidueSelect(keep_hetero=keep_hetero))


def convert_multiple():
    # Define the input and output directories
    output_pdb_directory = "output_pdbs/"

    # Check if the command line arguments are provided correctly
    if len(sys.argv) < 2:
        print("Please provide the necessary arguments: script_name.py Cifs_files-Folder/")
        sys.exit(1)

    # Path to the sequence
    input_cif_directory = str(sys.argv[1])

    # Create output and temporary directories if they don't exist
    os.makedirs(output_pdb_directory, exist_ok=True)

    # Get a list of CIF files in the input directory and sort them
    cif_files = [f for f in os.listdir(input_cif_directory) if f.endswith(".cif")]
    cif_files.sort()  # Sort alphabetically

    # Wrap the loop with tqdm to add a progress bar
    for cif_file in tqdm(cif_files, desc="Converting files"):
        cif_path = os.path.join(input_cif_directory, cif_file)
        pdb_file = cif_file.replace(".cif", ".pdb")
        pdb_path = os.path.join(output_pdb_directory, pdb_file)
        main(argparse.Namespace(input_cif=cif_path, output_pdb=pdb_path))

    print("Conversion complete.")


if __name__ == "__main__":
    main(parse_args())
