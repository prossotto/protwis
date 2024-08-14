from protein.models import Protein, ProteinSegment, Residue, ProteinConformation, ProteinCouplings
import csv
from django.core.management.base import BaseCommand
import pandas as pd 

import time
from django.db.models import Q

class Command(BaseCommand):

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

        protein_segments = Protein.objects.filter(
            species__common_name='Human', 
            entry_name__endswith='human',
            family__parent__slug__startswith='00'
            ).prefetch_related(
                "proteinconformation_set__residue_set",
                "family"
            )

        with open('/Users/vtk842/human_g_protein/gpcr_fasta_files/CN-term_seq/TMs.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['gpcr', 'class', 'sequence', 
                             'TM1', 'tm1_n_start', 'tm1_n_end',
                             'TM2', 'tm2_n_start', 'tm2_n_end',
                             'TM3', 'tm3_n_start', 'tm3_n_end',
                             'TM4', 'tm4_n_start', 'tm4_n_end',
                             'TM5', 'tm5_n_start', 'tm5_n_end',
                             'TM6', 'tm6_n_start', 'tm6_n_end',
                             'TM7', 'tm7_n_start', 'tm7_n_end',
                             ])

        for protein in protein_segments:
            print(protein.entry_name)
            for conformation in protein.proteinconformation_set.all():
                TM1_list = []
                tm1_n_list =[]

                TM2_list = []
                tm2_n_list = []

                TM3_list = []
                tm3_n_list = []

                TM4_list = []
                tm4_n_list = []

                TM5_list = []
                tm5_n_list = []

                TM6_list = []
                tm6_n_list = []

                TM7_list = []
                tm7_n_list = []

                for residue in conformation.residue_set.all():
                    if residue.protein_segment.slug == 'TM1':
                        TM1_list.append(residue.amino_acid)
                        tm1_n_list.append(residue.sequence_number)
                    if residue.protein_segment.slug == 'TM2':
                        TM2_list.append(residue.amino_acid)
                        tm2_n_list.append(residue.sequence_number)
                    if residue.protein_segment.slug == 'TM3':
                        TM3_list.append(residue.amino_acid)
                        tm3_n_list.append(residue.sequence_number)
                    if residue.protein_segment.slug == 'TM4':
                        TM4_list.append(residue.amino_acid)
                        tm4_n_list.append(residue.sequence_number)
                    if residue.protein_segment.slug == 'TM5':
                        TM5_list.append(residue.amino_acid)
                        tm5_n_list.append(residue.sequence_number)
                    if residue.protein_segment.slug == 'TM6':
                        TM6_list.append(residue.amino_acid)
                        tm6_n_list.append(residue.sequence_number)
                    if residue.protein_segment.slug == 'TM7':
                        TM7_list.append(residue.amino_acid)
                        tm7_n_list.append(residue.sequence_number)

                try:
                    TM1 = ''.join(TM1_list)
                    tm1_start = tm1_n_list[0]
                    tm1_end = tm1_n_list[-1]

                    TM2 = ''.join(TM2_list)
                    tm2_start = tm2_n_list[0]
                    tm2_end = tm2_n_list[-1]

                    TM3 = ''.join(TM3_list)
                    tm3_start = tm3_n_list[0]
                    tm3_end = tm3_n_list[-1]

                    TM4 = ''.join(TM4_list)
                    tm4_start = tm4_n_list[0]
                    tm4_end = tm4_n_list[-1]

                    TM5 = ''.join(TM5_list)
                    tm5_start = tm5_n_list[0]
                    tm5_end = tm5_n_list[-1]

                    TM6 = ''.join(TM6_list)
                    tm6_start = tm6_n_list[0]
                    tm6_end = tm6_n_list[-1]

                    TM7 = ''.join(TM7_list)
                    tm7_start = tm7_n_list[0]
                    tm7_end = tm7_n_list[-1]


                    with open('/Users/vtk842/human_g_protein/gpcr_fasta_files/CN-term_seq/TMs.csv', 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([protein.entry_name, protein.family.parent.parent.parent.name, protein.sequence, 
                                            TM1, tm1_start, tm1_end, 
                                            TM2, tm2_start, tm2_end, 
                                            TM3, tm3_start, tm3_end, 
                                            TM4, tm4_start, tm4_end, 
                                            TM5, tm5_start, tm5_end, 
                                            TM6, tm6_start, tm6_end, 
                                            TM7, tm7_start, tm7_end,
                                            ])
                
                except Exception as e:
                    print(e)
                    print(TM1)
                    print(tm1_n_list)




