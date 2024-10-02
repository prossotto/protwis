from build.management.commands.base_build import Command as BaseBuild
from protein.models import Protein
from residue.models import Residue
from structure.functions import PdbChainSelector
from structure.models import Structure, StructureModel
from Bio.PDB import PDBIO, PDBParser, Select, Structure as BioPDBStructure
# from Bio.PDB import PDBParser, Structure as BioPDBStructure
from Bio.PDB.Model import Model  # Import Model here
from collections import OrderedDict
import pandas as pd
from structure.assign_generic_numbers_gpcr import GenericNumbering as as_gn
from io import StringIO
import os
import shutil

class ResidueBFactorSelect(Select):
    """
    A selection class for filtering residues based on the B-factor of their CA atoms.

    Attributes:
    -----------
    bfactor_range : tuple
        A tuple specifying the inclusive range of B-factors to select residues.
    """
    def __init__(self, bfactor_range=(1.0, 7.0)):
        """
        Initializes the selection object with the specified B-factor range.

        Parameters:
        -----------
        bfactor_range : tuple, optional
            The inclusive range of B-factors to select residues. Default is (1.0, 7.0).
        """
        self.bfactor_range = bfactor_range

    def accept_residue(self, residue):
        """
        Checks if a residue should be accepted based on the B-factor of its CA atom.

        Parameters:
        -----------
        residue : Bio.PDB.Residue
            The residue to check.

        Returns:
        --------
        bool
            True if the residue should be accepted (CA atom's B-factor is within the range), False otherwise.
        """
        for atom in residue:
            if atom.get_id() == 'CA':
                bfactor = atom.get_bfactor()
                if self.bfactor_range[0] <= bfactor <= self.bfactor_range[1]:
                    return True
        return False


def extract_preferred_chain(pdb_data, preferred_chain, pdb_code):
    """
    Extracts the preferred chain from the PDB data and returns a model containing only that chain.

    Parameters:
    -----------
    pdb_data : str
        The PDB data as a string.
    preferred_chain : str
        The identifier of the preferred chain to extract.
    pdb_code : str
        The PDB code for identification purposes.

    Returns:
    --------
    Bio.PDB.Model.Model or None
        A model containing only the preferred chain, or None if the chain is not found.
    """

    parser = PDBParser(QUIET=True)
    pdb_io = StringIO(pdb_data)
    structure = parser.get_structure(pdb_code, pdb_io)

    # Get the first model
    model = structure[0]  # Assuming only one model

    # Create a new model with the same id
    new_model = Model(model.id)

    # Extract the preferred chain
    if preferred_chain in model:
        chain = model[preferred_chain]
        new_model.add(chain)
        return new_model
    else:
        print(f"Chain {preferred_chain} not found in structure {pdb_code}")
        return None

def create_db_dir(structure_type, output_dir):

    db_dir = f'{structure_type}_foldseek_db_trim'
    output_db_dir = os.path.join(output_dir, db_dir)
    if os.path.exists(output_db_dir):
        shutil.rmtree(output_db_dir)
    os.mkdir(output_db_dir)
    return output_db_dir

class Command(BaseBuild):
    help = "Assigns generic numbers to the preferred chain of PDB structures and extracts residues with generic numbers from 1 to 7."

    def handle(self, *args, **options):
        # Define the output directory
        output_dir = '/Users/vtk842/Peptide_complexes/foldseek/dbs'


        # Define the structure types to process
        structure_types = ['raw', 'ref', 'af']

        for structure_type in structure_types:
            if structure_type == 'raw':
                print(structure_type)
                # Query experimental structures from the database
                exp_structures = Structure.objects.filter(
                    structure_type__slug__in=['x-ray-diffraction', 'electron-microscopy', 'electron-crystallography']
                ).values_list(
                    'pdb_code__index',
                    'pdb_data__pdb',
                    'preferred_chain'
                )

                # Create directory
                output_raw_dir = create_db_dir(structure_type, output_dir)

                for struct in exp_structures:
                    pdb_code = struct[0]
                    pdb_data = struct[1]  # PDB data as a string
                    preferred_chain = struct[2]

                    print(f"Processing structure {pdb_code} with preferred chain {preferred_chain}")

                    # Extract the preferred chain
                    preferred_chain_structure = extract_preferred_chain(pdb_data, preferred_chain, pdb_code)
                    if preferred_chain_structure is None:
                        continue  # Skip to next structure if chain not found

                    # Initialize GenericNumbering with the preferred chain structure
                    gn = as_gn(structure=preferred_chain_structure, pdb_code=pdb_code)
                    gn.assign_generic_numbers()
                    annotated_structure = gn.get_annotated_structure()

                    # Use PDBIO to write out the selected residues
                    io = PDBIO()
                    io.set_structure(annotated_structure)
                    output_filename = f"{output_raw_dir}/{pdb_code}_raw_info.pdb"
                    io.save(output_filename, ResidueBFactorSelect())
                    print(f"Saved selected residues to {output_filename}")


            elif structure_type == 'af':

                print('af')

                af_structures = StructureModel.objects.filter(main_template_id__isnull=True
                ).values_list(
                    'protein__entry_name',
                    'pdb_data__pdb'
                )

                # Create directory
                output_af_dir = create_db_dir(structure_type, output_dir)

                for struct in af_structures:
                    pdb_code = struct[0]
                    pdb_data = struct[1]  # PDB data as a string

                    print(f"Processing structure {pdb_code} of  alphafold multistate structures")

                    parser = PDBParser(QUIET=True)
                    pdb_io = StringIO(pdb_data)
                    structure = parser.get_structure(pdb_code, pdb_io)

                    # Get the first model
                    model = structure[0]

                    # Initialize GenericNumbering with the preferred chain structure
                    gn = as_gn(structure=model, pdb_code=pdb_code)
                    gn.assign_generic_numbers()
                    annotated_structure = gn.get_annotated_structure()

                    # Use PDBIO to write out the selected residues
                    io = PDBIO()
                    io.set_structure(annotated_structure)
                    output_filename = f"{output_af_dir}/{pdb_code}_af_info.pdb"
                    io.save(output_filename, ResidueBFactorSelect())
                    print(f"Saved selected residues to {output_filename}")


            else:

                inactive_structures = StructureModel.objects.filter(main_template_id__isnull=False
                ).values_list(
                    'protein__entry_name',
                    'pdb_data__pdb'
                )

                active_structures = Structure.objects.filter(structure_type__slug__in=['af-signprot-refined-cem', 'af-signprot-refined-xray']
                ).values_list(
                    'pdb_code__index',
                    'pdb_data__pdb',
                    'preferred_chain'
                )

                # Create directory
                output_af_dir = create_db_dir(structure_type, output_dir)

                # As Structure Objects
                for struct in active_structures:
                    pdb_code = struct[0]
                    pdb_data = struct[1]  # PDB data as a string
                    preferred_chain = struct[2]

                    print(f"Processing structure {pdb_code} with preferred chain {preferred_chain}")

                    # Extract the preferred chain
                    preferred_chain_structure = extract_preferred_chain(pdb_data, preferred_chain, pdb_code)
                    if preferred_chain_structure is None:
                        continue  # Skip to next structure if chain not found

                    # Initialize GenericNumbering with the preferred chain structure
                    gn = as_gn(structure=preferred_chain_structure, pdb_code=pdb_code)
                    gn.assign_generic_numbers()
                    annotated_structure = gn.get_annotated_structure()

                    # Use PDBIO to write out the selected residues
                    io = PDBIO()
                    io.set_structure(annotated_structure)
                    output_filename = f"{output_af_dir}/{pdb_code}_ref_info.pdb"
                    io.save(output_filename, ResidueBFactorSelect())
                    print(f"Saved selected residues to {output_filename}")

                # As StructureModel Objects
                for struct in inactive_structures:
                    pdb_code = struct[0]
                    pdb_data = struct[1]  # PDB data as a string

                    print(f"Processing structure {pdb_code} of  alphafold multistate structures")

                    parser = PDBParser(QUIET=True)
                    pdb_io = StringIO(pdb_data)
                    structure = parser.get_structure(pdb_code, pdb_io)

                    # Get the first model
                    model = structure[0]

                    # Initialize GenericNumbering with the preferred chain structure
                    gn = as_gn(structure=model, pdb_code=pdb_code)
                    gn.assign_generic_numbers()
                    annotated_structure = gn.get_annotated_structure()

                    # Use PDBIO to write out the selected residues
                    io = PDBIO()
                    io.set_structure(annotated_structure)
                    output_filename = f"{output_af_dir}/{pdb_code}_ref_info.pdb"
                    io.save(output_filename, ResidueBFactorSelect())
                    print(f"Saved selected residues to {output_filename}")











