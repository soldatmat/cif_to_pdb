# original script: 2023 Jean-Marie Bourhis
#  adapted script: 2026 Matouš Soldát

import argparse
import os
import sys
from tqdm import tqdm

from Bio.PDB import PDBIO, Select
from Bio.PDB.MMCIFParser import MMCIFParser


class _StandardResidueSelect(Select):
    """Keep only standard (non-hetero, non-water) residues, i.e. ATOM records.
    Matches the original obabel+awk behaviour that filtered to `^ATOM` lines and
    dropped ligands/ions/water."""

    def accept_residue(self, residue):
        return residue.id[0] == " "


def main(args):
    # Convert CIF -> PDB with Biopython (NOT Open Babel). Open Babel did not carry
    # the mmCIF B_iso_or_equiv column into the PDB B-factor (wrote 0.00), which
    # discarded AlphaFold pLDDT. MMCIFParser reads B_iso into atom.bfactor and PDBIO
    # writes it back to the B-factor column, so the extracted PDB keeps per-residue
    # pLDDT. Only standard residues are written (protein/nucleic ATOM records),
    # preserving the previous ATOM-only output shape.
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("s", args.input_cif)
    io = PDBIO()
    io.set_structure(structure)
    io.save(args.output_pdb, _StandardResidueSelect())


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
    args = parse_args()
    main(args)
