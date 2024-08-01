from protein.models import Protein, ProteinSegment, Residue, ProteinConformation, ProteinCouplings
from structure.models import Structure
from structure.assign_generic_numbers_gpcr import *
import csv
from django.core.management.base import *
import pandas as pd 
from Bio.PDB import PDBParser, PDBIO, Select
from io import StringIO

import re
import time
from django.db.models import Q
import shutil
from pathlib import Path

class ResidueBFactorSelect(Select):
    """
    A selection class for filtering residues based on the B-factor of their CA atoms.
    
    Attributes:
    -----------
    bfactor_range : tuple
        A tuple specifying the inclusive range of B-factors to select residues.
    """
    def __init__(self, bfactor_range=(7.0, 9.0)):
        """
        Initializes the selection object with the specified B-factor range.
        
        Parameters:
        -----------
        bfactor_range : tuple, optional
            The inclusive range of B-factors to select residues. Default is (7.0, 8.9).
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
                if self.bfactor_range[0] <= bfactor < self.bfactor_range[1]:
                    return True
        return False

class Command(BaseCommand):
    

    def handle(self, *args, **options):
        print('HERE we ARE ##################')
        start_time = time.time()

        output_dir = Path('/Users/vtk842/Peptide_complexes/tm7_h8_db')

        # deleting db dirs
        for dir in output_dir.iterdir():
            if dir.is_dir():
                shutil.rmtree(dir)

        dirs = ['raw', 'raw_1', 'ref', 'ref_1', 'af', 'af_1']

        # creating db dirs
        for dir in dirs:
            create_dir = output_dir / dir
            create_dir.mkdir(exist_ok=True, parents=True)

        # Testing for a single PDB
        #---------------------------

        # input_pdb = '/Users/vtk842/checking/protwis/test_foldseek/1F88_raw_info.pdb'
        # output = '/Users/vtk842/checking/protwis/test_foldseek/tm7_test/1F88_test_raw_info.pdb'
        input_pdb = Path('/Users/vtk842/checking/protwis/test_foldseek/1F88_raw_T.pdb')
        output = f'/Users/vtk842/checking/protwis/test_foldseek/{input_pdb.stem}_annotated_filtered.pdb'
        output_normal = f'/Users/vtk842/checking/protwis/test_foldseek/{input_pdb.stem}_annotated.pdb'

        annotator = GenericNumbering(pdb_file=input_pdb)
        annotated = annotator.assign_generic_numbers()
        # annotated = annotator.assign_generic_numbers_with_sequence_parser()

        print(annotated)
        print(type(annotated))

        io = PDBIO()
        io.set_structure(annotated)
        io.save(str(output), ResidueBFactorSelect())
        io.save(output_normal)


        #############################################
        ############ Generate TM7 H8 Data base ######
        #############################################

        # raw_types = ['1', '2', '3']
        # ref_types = ['4', '5']
        # af_types = ['7']

        # pdb_files = Structure.objects.filter(
        #         structure_type__in=raw_types + ref_types + af_types
        #     ).select_related(
        #         "pdb_data",       # Fetch the related pdb_data object
        #         "pdb_code"        # Fetch the related pdb_code object
        #     )

        # for structure in pdb_files:
        #     pdb_data = structure.pdb_data.pdb       # Accessing pdb_data
        #     pdb_id = structure.id                   # Accessing ID
        #     pdb_code = structure.pdb_code.index     # Accessing pdb_code
        #     s_type = str(structure.structure_type.id)    # Accessing structure_type

        #     print(f' S TYPE: {s_type}')

        #     # Set structure_type based on the query value
        #     if s_type in raw_types:
        #         structure_type = "raw"
        #     elif s_type in ref_types:
        #         structure_type = "ref"
        #     elif s_type in af_types:
        #         structure_type = "af"
        #     else:
        #         structure_type = 'empty'

        #     # print(pdb_id,  pdb_code)  # This will print the ID

        #     # Assuming GenericNumberingFromDB expects a Structure object and pdb_data
        #     assign_grn = GenericNumberingFromDB(structure, pdb_data)
        #     structure_grn = assign_grn.assign_generic_numbers()
        #     io = PDBIO()
        #     io.set_structure(structure_grn)
        #     print('STRUCTURE #####################')
        #     io.save(f'{output_dir}/{structure_type}/{pdb_code}_{structure_type}_info.pdb')
