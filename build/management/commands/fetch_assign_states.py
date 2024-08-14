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
        output = '/Users/vtk842/human_g_protein/gpcr_fasta_files/CN-term_seq/TM_gns.csv'

        headers = ['REC','sequence', 'TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7']

        protein_segments = Protein.objects.filter(
            species__common_name='Human', 
            entry_name__endswith='human',
            family__parent__slug__startswith='00'
            ).prefetch_related(
                "proteinconformation_set__residue_set",
                "family"
            )

        with open(output, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['gpcr', 'class', 'sequence', 
                             'TM1', 'tm1_n', 'tm1_gn',
                             'TM2', 'tm2_n', 'tm2_gn',
                             'TM3', 'tm3_n', 'tm3_gn',
                             'TM4', 'tm4_n', 'tm4_gn',
                             'TM5', 'tm5_n', 'tm5_gn',
                             'TM6', 'tm6_n', 'tm6_gn',
                             'TM7', 'tm7_n', 'tm7_gn',
                             ])

        for protein in protein_segments:
            print(protein.entry_name)
            for conformation in protein.proteinconformation_set.all():
                TM1_list = []
                tm1_n_list =[]
                tm1_gn_list = []
                
                TM2_list = []
                tm2_n_list = []
                tm2_gn_list = []

                TM3_list = []
                tm3_n_list = []
                tm3_gn_list = []

                TM4_list = []
                tm4_n_list = []
                tm4_gn_list = []

                TM5_list = []
                tm5_n_list = []
                tm5_gn_list = []

                TM6_list = []
                tm6_n_list = []
                tm6_gn_list = []

                TM7_list = []
                tm7_n_list = []
                tm7_gn_list = []

                for residue in conformation.residue_set.all():
                    if residue.protein_segment.slug == 'TM1':
                        TM1_list.append(residue.amino_acid)
                        tm1_n_list.append(residue.sequence_number)
                        tm1_gn_list.append(residue.generic_number.label)
                    if residue.protein_segment.slug == 'TM2':
                        TM2_list.append(residue.amino_acid)
                        tm2_n_list.append(residue.sequence_number)
                        tm2_gn_list.append(residue.generic_number.label)
                    if residue.protein_segment.slug == 'TM3':
                        TM3_list.append(residue.amino_acid)
                        tm3_n_list.append(residue.sequence_number)
                        tm3_gn_list.append(residue.generic_number.label)
                    if residue.protein_segment.slug == 'TM4':
                        TM4_list.append(residue.amino_acid)
                        tm4_n_list.append(residue.sequence_number)
                        tm4_gn_list.append(residue.generic_number.label)
                    if residue.protein_segment.slug == 'TM5':
                        TM5_list.append(residue.amino_acid)
                        tm5_n_list.append(residue.sequence_number)
                        tm5_gn_list.append(residue.generic_number.label)
                    if residue.protein_segment.slug == 'TM6':
                        TM6_list.append(residue.amino_acid)
                        tm6_n_list.append(residue.sequence_number)
                        tm6_gn_list.append(residue.generic_number.label)
                    if residue.protein_segment.slug == 'TM7':
                        TM7_list.append(residue.amino_acid)
                        tm7_n_list.append(residue.sequence_number)
                        tm7_gn_list.append(residue.generic_number.label)

                try:
                    TM1 = ''.join(TM1_list)
                    tm1_n = tm1_n_list
                    tm1_gn = tm1_gn_list

                    TM2 = ''.join(TM2_list)
                    tm2_n = tm2_n_list
                    tm2_gn = tm2_gn_list

                    TM3 = ''.join(TM3_list)
                    tm3_n = tm3_n_list
                    tm3_gn = tm3_gn_list

                    TM4 = ''.join(TM4_list)
                    tm4_n = tm4_n_list
                    tm4_gn = tm4_gn_list

                    TM5 = ''.join(TM5_list)
                    tm5_n = tm5_n_list
                    tm5_gn = tm5_gn_list

                    TM6 = ''.join(TM6_list)
                    tm6_n = tm6_n_list
                    tm6_gn = tm6_gn_list

                    TM7 = ''.join(TM7_list)
                    tm7_n = tm7_n_list
                    tm7_gn = tm7_gn_list


                    with open(output, 'a') as f:
                        writer = csv.writer(f)
                        writer.writerow([protein.entry_name, protein.family.parent.parent.parent.name, protein.sequence, 
                                            TM1, tm1_n, tm1_gn, 
                                            TM2, tm2_n, tm2_gn, 
                                            TM3, tm3_n, tm3_gn, 
                                            TM4, tm4_n, tm4_gn, 
                                            TM5, tm5_n, tm5_gn, 
                                            TM6, tm6_n, tm6_gn, 
                                            TM7, tm7_n, tm7_gn, 
                                            ])
                
                except Exception as e:
                    print(e)
                    print(TM1)
                    print(tm1_n_list)

