from protein.models import Protein, ProteinSegment, Residue, ProteinConformation, ProteinCouplings
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
        start_time = time.time()
        # output = '/Users/vtk842/human_g_protein/gpcr_fasta_files/pre_phosphorylation.csv'
        output = '/Users/vtk842/human_g_protein/ArrestinDB_2023/pre_phosphorylation.csv'

        headers = ['REC','sequence', 'i2', 'i3', 'Ct']
        with open(output, 'w',newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(headers)

        start = time.time()
        prots = Protein.objects.filter(
                    accession__isnull=False, 
                    species__latin_name='Homo sapiens'
                ).exclude(
                    Q(family__slug__startswith='100') |
                    Q(family__slug__startswith='200')
                ).prefetch_related(
                    'proteinconformation_set__residue_set__protein_segment'
                )

        protein_data = []

        for protein in prots:
            sequences = {}
            sequences['gpcr'] =(protein.entry_name.split('_')[0].upper())
            sequences['gtodb'] = (protein.name.replace('receptor', ''))
            sequences['family'] = (protein.family.parent.parent.parent.name.split(' ')[1])
            sequences['class'] = (protein.family.parent.parent.name.split(' ')[0])
            segment_sequences = {'ICL2': '', 'ICL3': '', 'C-term': ''}

            # Iterate through each protein conformation
            for conf in protein.proteinconformation_set.all():
                # Build sequences for each segment
                for seg in conf.residue_set.all():
                    slug = seg.protein_segment.slug
                    if slug in segment_sequences:
                        segment_sequences[slug] += seg.amino_acid

            # sequences = {segment.lower(): seq for segment, seq in segment_sequences.items()}

            for segment, seq in segment_sequences.items():
                sequences[segment.lower()] = seq

            protein_data.append(sequences)

        print(time.time() - start)
        print(protein_data[0])
        


##################### analyse sequence



            



