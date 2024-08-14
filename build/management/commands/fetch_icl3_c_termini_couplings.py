from protein.models import Protein, ProteinSegment, Residue, ProteinConformation, ProteinCouplings
import csv
from django.core.management.base import BaseCommand
from ligand.models import Endogenous_GTP
import pandas as pd 
import time
from django.db.models import Subquery, OuterRef, Max, FloatField, Window, F
from collections import defaultdict


class Command(BaseCommand):

    def handle(self, *args, **options):
        start_time = time.time()
        # output = '/Users/vtk842/human_g_protein/gpcr_fasta_files/pre_phosphorylation.csv'
        output = '//Users/vtk842/Peptide_complexes/end_human_ligands.csv'

        headers = ['ligand', 'ligand_name', 'ligand_type', 'ligand_seq', 'receptor', 'receptor family', 'receptor class', 'receptor_sequence']
        with open(output, 'w',newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(headers)

        # data = {
        #     'gpcr':[],
        #     'sequence':[],
        #     ('icl2', 'sequence'):[],
        #     ('icl3', 'sequence'):[],
        #     ('c-term', 'sequence'):[],
        # }

        coupling_data = ProteinCouplings.objects.filter(
                protein__family__slug__startswith='00',
                protein__entry_name__endswith='human',
                g_protein_subunit__family__slug__startswith='100'
            ).values_list(
                "protein__entry_name",                          #0 GPCR name
                "g_protein_subunit__entry_name",                #1 g protein
                "g_protein_subunit__sequence",                  #2 g protein sequence
                "logemaxec50",                                  #3 log max ec50
                "protein__family__parent__parent__parent__name"       #4 class
            )
        
        # Dictionary to store coupling data
        detailed_coupling_data = defaultdict(lambda: defaultdict(list))

        # Open CSV file for writing
        with open('/Users/vtk842/Peptide_complexes/g_protein_best_mean_couplings.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['gpcr', 'class', 'logmax', 'g_protein', 'g_protein_seq'])
            
            for data in coupling_data:
                if data[3] is not None and data[3] > 0:  # Check if logemaxec50 is greater than 0
                    gpcr = data[0]
                    current_logmax = data[3]
                    current_g_protein = data[1]
                    current_g_protein_seq = data[2]
                    gpcr_class = data[4]
                    
                    # Store the logmax value
                    detailed_coupling_data[(gpcr, current_g_protein)]['logmax'].append(current_logmax)
                    detailed_coupling_data[(gpcr, current_g_protein)]['g_protein_seq'] = current_g_protein_seq
                    detailed_coupling_data[(gpcr, current_g_protein)]['class'] = gpcr_class

            # Calculate the mean logmax and write to CSV
            for (gpcr, g_protein), values in detailed_coupling_data.items():
                mean_logmax = sum(values['logmax']) / len(values['logmax'])
                writer.writerow([gpcr, values['class'], mean_logmax, g_protein, values['g_protein_seq']])

        # detailed_coupling_data = {}
        # with open('/Users/vtk842/Peptide_complexes/g_protein_best_couplings.csv', 'w') as f:
        #     writer = csv.writer(f)
        #     writer.writerow(['gpcr','class' , 'logmax', 'g_protein', 'g_protein_seq'])
        # for data in coupling_data:
        #     if data[3] is not None and data[3] > 0:  # Check if logemaxec50 is greater than 0
        #         gpcr = data[0]
        #         current_logmax = data[3]
        #         current_g_protein = data[1]
        #         current_g_protein_seq = data[2]
        #         gpcr_class = data[4]

        #         # Check if the GPCR is already in the dictionary
        #         if gpcr in detailed_coupling_data:
        #             # Update if the current logmax is greater than the one stored
        #             if current_logmax > detailed_coupling_data[gpcr]['logmax']:
        #                 detailed_coupling_data[gpcr] = {'logmax': current_logmax, 'g_protein': current_g_protein, 'g_protein_seq': current_g_protein_seq, 'class': gpcr_class}
        #         else:
        #             # Add the new GPCR with its logmax and details
        #             detailed_coupling_data[gpcr] = {'logmax': current_logmax, 'g_protein': current_g_protein, 'g_protein_seq': current_g_protein_seq, 'class': gpcr_class}

        for key, value in  detailed_coupling_data.items():
            with open('/Users/vtk842/Peptide_complexes/g_protein_best_couplings.csv', 'a') as f:
                writer = csv.writer(f)
                writer.writerow([key, value['class'], value['logmax'], value['g_protein'], value['g_protein_seq']])

        # for data in endogenous_data:
        #     ligand = data[6]
        #     ligand_name = data[5]
        #     ligand_type = data[9]
        #     ligand_seq = data[19]
        #     receptor = data[2]
        #     receptor_family = data[1]
        #     receptor_class = data[0]
        #     receptor_sequence = data[20]

        #     print(receptor)

        #     with open(output, 'a') as f:
        #         writer = csv.writer(f)
        #         writer.writerow([ligand, ligand_name, ligand_type, ligand_seq, receptor, receptor_family, receptor_class, receptor_sequence])


