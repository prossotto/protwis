from protein.models import Protein, ProteinSegment, Residue, ProteinConformation, ProteinCouplings
from structure.models import Structure
import csv
from django.core.management.base import BaseCommand
import pandas as pd 

import re
import time
from django.db.models import Q

class Command(BaseCommand):

    patterns = {
    #              P    x    P    P
    'PxPP': r'((?:S|T)[^ST](?:S|T)(?:S|T))',
    #       P    x    P    x     x  P/E/D
    'PxPxxPED':r'((?:S|T)[^ST](?:S|T)[^ST][^ST](?:S|T|E|D))',
    #       P    x    x    P    x    x    P/E/D
    'PxxPxxPED':r'((?:S|T)[^ST][^ST](?:S|T)[^ST][^ST](?:S|T|E|D))',
    #
    'longest':r'((?:(?:S|T|E|D)[^STED]{1,2})+(?:S|T|E|D))',

    }

    def analyze_sequences(self, data, column, patterns):

        data[f'{column}_seq'] = data[column]
        data[f'{column}_ST'] = data[column].count('S') + data[column].count('T') 
        data[f'{column}_DE'] = data[column].count('D') + data[column].count('E') 
        # print(f' column: {column}\n{data[f"{column}"]}\n{data[f"{column}_ST"]}')
        data[f'{column}_STDE'] = data[f'{column}_ST'] + data[f'{column}_DE']


        for key, val in patterns.items():
            if key != 'longest':
                # data[f'{column}_{key}'] = data[column].count(val)
                data[f'{column}_{key}'] = [len(re.findall(val, seq)) for seq in data[column]]

            else:
                continue
        keys_to_sum = [f'{column}_{key}' for key in patterns.keys() if key != 'longest']
        data[f'{column}_All'] = [sum(data[key] for key in keys_to_sum)]
        data[f'{column}_sites'] = re.findall(patterns['longest'], data[column])
        data[f'{column}_longest'] = [max(len(x) for x in site_list) if site_list else 0 for site_list in data[f'{column}_sites']]
        data[f'{column}_long_STDE'] = [self.count_specific_letters_(row) for row in data[f'{column}_sites']]


    def amino_s(self, iterator='', gpcr=''):
        sequence = ''
        sequence_pos = ''
        res = Residue.objects.filter(protein_conformation__protein__entry_name=gpcr, protein_segment__slug=iterator).all()
        for residue in res:
            sequence += str(residue.amino_acid)
            sequence_pos += f'{str(residue.sequence_number)}'
        return sequence

    def handle(self, *args, **options):

        state_pref = 'inactive'

        maper = {'active':'2', 'inactive':'1'}

        output = f'/Users/vtk842/Peptide_complexes/multistate/classification/exp_{state_pref}'

        structure_info = Structure.objects.filter(
            state=f'{maper[state_pref]}', structure_type__in=['1', '2', '3'],
        ).values_list(
            'pdb_data__pdb',                                     #0
            'state__slug',                                       #1
            'protein_conformation__protein__entry_name',         #2
        )

        for struc in structure_info:

            name = f'{struc[2]}_{state_pref}.pdb'

            with open(f'{output}/{name}', 'w') as f:
                f.write(struc[0])

    