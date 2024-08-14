from protein.models import Protein, ProteinSegment, Residue, ProteinConformation
import csv
from django.core.management.base import BaseCommand
from ligand.models import Endogenous_GTP
import pandas as pd 
import time

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
        output = '/Users/vtk842/Peptide_complexes/end_Xhuman_ligands.csv'

        headers = ['ligand', 'ligand_name', 'ligand_type', 'ligand_status', 'ligand_uniprot', 'ligand_seq', 'receptor', 'receptor family', 'receptor class', 'receptor_sequence']
        # headers = ['ligand', 'ligand_name', 'ligand_type', 'ligand_seq', 'receptor', 'receptor family', 'receptor class', 'receptor_sequence']
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

        endogenous_data = Endogenous_GTP.objects.filter(
            ligand__ligand_type__name__in=['peptide', 'protein'],
                        receptor__species__common_name='Human',
            ).values_list(
                            "receptor__family__parent__parent__parent__name", #0 Class
                            "receptor__family__parent__name",                 #1 Receptor Family
                            "receptor__entry_name",                           #2 UniProt
                            "receptor__name",                                 #3 IUPHAR
                            "receptor__species__common_name",                 #4 Species
                            "ligand__name",                                   #5 Ligand
                            "ligand",                                         #6 Ligand ID
                            "endogenous_status",                              #7 Principal/Secondary
                            "potency_ranking",                                #8 Potency Ranking
                            "ligand__ligand_type__name",                      #9 Type
                            "pec50",                                          #10 pEC50 - min - med - max
                            "pKi",                                            #11 pKi - min - med - max
                            "publication__authors",                           #12 Pub Authors
                            "publication__year",                              #13 Pub Year
                            "publication__title",                             #14 Pub Title
                            "publication__journal__name",                     #15 Pub Journal
                            "publication__reference",                         #16 Pub Reference
                            "publication__web_link__index",                   #17 DOI/PMID
                            "receptor",                                       #18 Receptor ID
                            "ligand__sequence",                               #19 sequence
                            'receptor__sequence',                             #20 receptor sequence
                            'endogenous_status',                              #21 status
                            'ligand__uniprot',                                #22 ligand_uniprot
                            "receptor__accession").distinct()                 #23 Accession (UniProt link)
        
        # endogenous_data = Endogenous_GTP.objects.all().values_list(
        #             "receptor__family__parent__parent__parent__name", #0 Class
        #             "receptor__family__parent__name",                 #1 Receptor Family
        #             "receptor__entry_name",                           #2 UniProt
        #             "receptor__name",                                 #3 IUPHAR
        #             "receptor__species__common_name",                 #4 Species
        #             "ligand__name",                                   #5 Ligand
        #             "ligand",                                         #6 Ligand ID
        #             "endogenous_status",                              #7 Principal/Secondary
        #             "potency_ranking",                                #8 Potency Ranking
        #             "ligand__ligand_type__name",                      #9 Type
        #             "pec50",                                          #10 pEC50 - min - med - max
        #             "pKi",                                            #11 pKi - min - med - max
        #             "publication__authors",                           #12 Pub Authors
        #             "publication__year",                              #13 Pub Year
        #             "publication__title",                             #14 Pub Title
        #             "publication__journal__name",                     #15 Pub Journal
        #             "publication__reference",                         #16 Pub Reference
        #             "publication__web_link__index",                   #17 DOI/PMID
        #             "receptor",                                       #18 Receptor ID
        #             "receptor__sequence"                              #19 sequence
        #             "receptor__accession").distinct()                 #19 Accession (UniProt link)

        for data in endogenous_data:
            ligand = data[6]
            ligand_name = data[5]
            ligand_type = data[9]
            ligand_seq = data[19]
            receptor = data[2]
            receptor_family = data[1]
            receptor_class = data[0]
            receptor_sequence = data[20]
            ligand_status = data[21]
            ligand_uniprot = data[22]
        # i = 0
        # for data in endogenous_data:
        #     i += 1
        #     ligand = data[6]
        #     ligand_name = data[5]
        #     ligand_type = data[9]
        #     ligand_seq = data[19]
        #     receptor = data[2]
        #     receptor_family = data[1]
        #     receptor_class = data[0]

            print(receptor, ligand_name)

            with open(output, 'a') as f:
                writer = csv.writer(f)
                writer.writerow(
                    [
                        ligand, ligand_name, ligand_type, ligand_status, ligand_uniprot, ligand_seq, 
                        receptor, receptor_family, receptor_class, receptor_sequence]
                    )
            
            # with open(output, 'a') as f:
            #     writer = csv.writer(f)
            #     writer.writerow(
            #         [
            #             ligand, ligand_name, ligand_type, ligand_seq, 
            #             receptor, receptor_family, receptor_class]
            #         )


