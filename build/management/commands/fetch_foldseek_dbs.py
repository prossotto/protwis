from protein.models import Protein, ProteinSegment, Residue, ProteinConformation, ProteinCouplings
from structure.models import Structure, StructureModel
from structure.assign_generic_numbers_gpcr import *
import csv
from django.core.management.base import *
import pandas as pd 
from Bio.PDB import PDBParser, PDBIO, Select
from io import StringIO

import re
import time
from django.db.models import Q, OuterRef, Prefetch
import shutil
from pathlib import Path

# class ResidueBFactorSelect(Select):
#     """
#     A selection class for filtering residues based on the B-factor of their CA atoms.
    
#     Attributes:
#     -----------
#     bfactor_range : tuple
#         A tuple specifying the inclusive range of B-factors to select residues.
#     """
#     def __init__(self, bfactor_range=(1.0, 9.0)):
#         """
#         Initializes the selection object with the specified B-factor range.
        
#         Parameters:
#         -----------
#         bfactor_range : tuple, optional
#             The inclusive range of B-factors to select residues. Default is (7.0, 8.9).
#         """
#         self.bfactor_range = bfactor_range

#     def accept_residue(self, residue):
#         """
#         Checks if a residue should be accepted based on the B-factor of its CA atom.
        
#         Parameters:
#         -----------
#         residue : Bio.PDB.Residue
#             The residue to check.
        
#         Returns:
#         --------
#         bool
#             True if the residue should be accepted (CA atom's B-factor is within the range), False otherwise.
#         """
#         for atom in residue:
#             if atom.get_id() == 'CA':
#                 bfactor = atom.get_bfactor()
#                 if self.bfactor_range[0] <= bfactor < self.bfactor_range[1]:
#                     return True
#         return False

# class Command(BaseCommand):


#     def handle(self, *args, **options):

#         def process_pdb(input_pdb, output_pdb, bfactor_range=(1.0, 9.0)):
#             # Parse the PDB file
#             parser = PDBParser()
#             structure = parser.get_structure('structure', input_pdb)

#             # Create an instance of ResidueBFactorSelect with the specified B-factor range
#             selector = ResidueBFactorSelect(bfactor_range=bfactor_range)

#             # Write the selected residues to a new PDB file
#             io = PDBIO()
#             io.set_structure(structure)
#             io.save(output_pdb, selector)

#             # # Example usage
#             # input_pdb = 'input.pdb'  # Path to your input PDB file
#             # output_pdb = 'output.pdb'  # Path to save the filtered PDB file
#             # bfactor_range = (1.0, 9.0)  # Set the desired B-factor range

#             # process_pdb(input_pdb, output_pdb, bfactor_range)


#         print('HERE we ARE ##################')
#         start_time = time.time()

#         output_dir = Path('/Users/vtk842/Peptide_complexes/tm7_h8_db')

#         # deleting db dirs
#         for dir in output_dir.iterdir():
#             if dir.is_dir():
#                 shutil.rmtree(dir)

#         db_id = 'foldseek_db'
#         dirs = [f'raw_{db_id}', f'raw_{db_id}_trim', f'ref_{db_id}', f'ref_{db_id}_trim', f'af_{db_id}', f'af_{db_id}_trim']

#         # creating db dirs
#         for dir in dirs:
#             create_dir = output_dir / dir
#             create_dir.mkdir(exist_ok=True, parents=True)

#         #############################################
#         ############ Generate 7TM H8 Data base ######
#         #############################################

#         raw_types = ['1', '2', '3']
#         ref_types = ['4', '5']
#         # af_types = ['7']

#         pdb_files = Structure.objects.filter(
#                 structure_type__in=raw_types + ref_types
#             ).select_related(
#                 "pdb_data",       # Fetch the related pdb_data object
#                 "pdb_code"        # Fetch the related pdb_code object
#             )
        
#         af_multistate = StructureModel.objects.filter(
#             main_template__isnull=True
#         ).select_related(
#             'pdb_data',
#             'state',
#             'protein'
#         ).prefetch_related(
#             Prefetch(
#                 'protein__proteinconformation_set',
#                 queryset=ProteinConformation.objects.filter(protein__isnull=False),
#                 to_attr='protein_conformation'
#             )
#         )

#         for structure in pdb_files:
#             pdb_data = structure.pdb_data.pdb       # Accessing pdb_data
#             pdb_id = structure.id                   # Accessing ID
#             pdb_code = structure.pdb_code.index     # Accessing pdb_code
#             s_type = str(structure.structure_type.id)    # Accessing structure_type

#             print(f' S TYPE: {s_type}')

#             # Set structure_type based on the query value
#             if s_type in raw_types:
#                 structure_type = "raw"
#             elif s_type in ref_types:
#                 structure_type = "ref"
#             elif s_type in af_types:
#                 structure_type = "af"
#             else:
#                 structure_type = 'empty'

#             print(pdb_id,  pdb_code)  # This will print the ID

#             # Assuming GenericNumberingFromDB expects a Structure object and pdb_data
#             assign_grn = GenericNumberingFromDB(structure, pdb_data)
#             structure_grn = assign_grn.assign_generic_numbers()
#             io = PDBIO()
#             io.set_structure(structure_grn)
#             print('STRUCTURE #####################')
#             output_pdb = f'{output_dir}/{structure_type}_{db_id}/{pdb_code}_{structure_type}_info.pdb'
#             io.save(output_pdb)

#             out_put_7tm_h8 = f'{output_dir}/{structure_type}_{db_id}_trim/{pdb_code}_{structure_type}_info.pdb'
#             process_pdb(output_pdb, out_put_7tm_h8)
        
#         for structure in af_multistate:
#             structure_type = 'af'
#             pdb_data = structure.pdb_data.pdb       # Accessing pdb_data
#             pdb_id = structure.id                   # Accessing ID
#             pdb_code = structure.protein.entry_name     # Accessing pdb_code
#             state = structure.state.slug

#             print('Protein_conformation')
#             print(structure.protein.protein_conformation)

#             print(pdb_id,  pdb_code)  # This will print the ID

#             assign_grn = GenericNumberingFromDB1(structure, pdb_data)
#             structure_grn = assign_grn.assign_generic_numbers()
#             io = PDBIO()
#             io.set_structure(structure_grn)
#             print('MultiState #####################')
#             output_pdb = f'{output_dir}/{structure_type}_{db_id}/{pdb_code}_{state}_{structure_type}_info.pdb'
#             io.save(output_pdb)

#             out_put_7tm_h8 = f'{output_dir}/{structure_type}_{db_id}_trim/{pdb_code}_{state}_{structure_type}_info.pdb'
#             process_pdb(output_pdb, out_put_7tm_h8)

class ResidueBFactorSelect(Select):
    """A selection class for filtering residues based on the B-factor of their CA atoms."""

    def __init__(self, bfactor_range=(1.0, 9.0)):
        self.bfactor_range = bfactor_range

    def accept_residue(self, residue):
        for atom in residue:
            if atom.get_id() == 'CA' and self.bfactor_range[0] <= atom.get_bfactor() < self.bfactor_range[1]:
                return True
        return False


class Command(BaseCommand):
    def handle(self, *args, **options):
        def process_pdb(input_pdb, output_pdb, bfactor_range=(1.0, 9.0)):
            parser = PDBParser()
            structure = parser.get_structure('structure', input_pdb)
            selector = ResidueBFactorSelect(bfactor_range)
            io = PDBIO()
            io.set_structure(structure)
            io.save(output_pdb, selector)

        start_time = time.time()
        output_dir = Path('/Users/vtk842/Peptide_complexes/tm7_h8_db')
        db_id = 'foldseek_db'
        raw_types = ['1', '2', '3']
        ref_types = ['4', '5']
        af_types = ['7']

        # Clear and recreate output directories
        shutil.rmtree(output_dir, ignore_errors=True)
        dirs = [f'{prefix}_{db_id}{suffix}' for prefix in ['raw', 'ref', 'af'] for suffix in ['', '_trim']]
        for dir in dirs:
            (output_dir / dir).mkdir(parents=True, exist_ok=True)

        # Process raw and reference structures
        structures = Structure.objects.filter(structure_type__in=raw_types + ref_types).select_related("pdb_data", "pdb_code")
        for structure in structures:
            process_structure(structure, raw_types, ref_types, af_types, db_id, output_dir)

        # Process AF multistate structures
        af_multistate = StructureModel.objects.filter(main_template__isnull=True).select_related('pdb_data', 'state', 'protein').prefetch_related(
            Prefetch('protein__proteinconformation_set', queryset=ProteinConformation.objects.filter(protein__isnull=False), to_attr='protein_conformation')
        )
        for structure in af_multistate:
            process_af_structure(structure, db_id, output_dir)

def process_structure(structure, raw_types, ref_types, af_types, db_id, output_dir):
    pdb_data = structure.pdb_data.pdb
    pdb_code = structure.pdb_code.index
    s_type = str(structure.structure_type.id)
    structure_type = "raw" if s_type in raw_types else "ref" if s_type in ref_types else "af" if s_type in af_types else 'empty'
    assign_grn = GenericNumberingFromDB(structure, pdb_data)
    save_structure(assign_grn.assign_generic_numbers(), structure_type, db_id, pdb_code, output_dir)
    

def process_af_structure(structure, db_id, output_dir):
    pdb_data = structure.pdb_data.pdb
    pdb_code = structure.protein.entry_name
    state = structure.state.slug
    assign_grn = GenericNumberingFromDB1(structure, pdb_data)
    save_structure(assign_grn.assign_generic_numbers(), 'af', db_id, f'{pdb_code}_{state}', output_dir)


def save_structure(structure, structure_type, db_id, pdb_code, output_dir):
    io = PDBIO()
    io.set_structure(structure)
    output_pdb = f'{output_dir}/{structure_type}_{db_id}/{pdb_code}_{structure_type}_info.pdb'
    io.save(output_pdb)
    trimmed_output_pdb = f'{output_dir}/{structure_type}_{db_id}_trim/{pdb_code}_{structure_type}_info.pdb'
    process_pdb(output_pdb, trimmed_output_pdb)


def process_pdb(input_pdb, output_pdb, bfactor_range=(1.0, 9.0)):
    parser = PDBParser()
    structure = parser.get_structure('structure', input_pdb)
    selector = ResidueBFactorSelect(bfactor_range)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, selector)


