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
        output = '/Users/vtk842/Peptide_complexes/RFAA/set/end_Xhuman_ligands_sm.csv'

        headers = ['ligand', 'ligand_name', 'ligand_type', 'ligand_status', 'pec50', 'ligand_inchikey', 'receptor', 'receptor family', 'receptor class', 'receptor_sequence']
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
            ligand__ligand_type__name='small-molecule',
            receptor__species__common_name='Human'
            ).values_list(
                            "receptor__family__parent__parent__parent__name", #0 Class
                            "receptor__family__parent__name",                 #1 Receptor Family
                            "receptor__entry_name",                           #2 UniProt
                            "receptor__name",                                 #3 IUPHAR
                            "receptor__species__common_name",                 #4 Species
                            "ligand__name",                                   #5 Ligand
                            "ligand",                                         #6 Ligand ID
                            "ligand_action__slug",                            #7 action (agonist, antagonist, etc)
                            "endogenous_status",                              #8 Principal/Secondary
                            "potency_ranking",                                #9 Potency Ranking
                            "ligand__ligand_type__name",                      #10 Type
                            "pec50",                                          #11 pEC50 - min - med - max
                            "pKi",                                            #12 pKi - min - med - max
                            "publication__authors",                           #13 Pub Authors
                            "publication__year",                              #14 Pub Year
                            "publication__title",                             #15 Pub Title
                            "publication__journal__name",                     #16 Pub Journal
                            "publication__reference",                         #17 Pub Reference
                            "publication__web_link__index",                   #18 DOI/PMID
                            "receptor",                                       #19 Receptor ID
                            "ligand__inchikey",                               #20 sequence
                            'receptor__sequence',                             #21 receptor sequence
                            'endogenous_status',                              #22 status
                            'pec50',                                          #23 pec50
                            "receptor__accession").distinct()                 #24 Accession (UniProt link)

        for data in endogenous_data:
            ligand = data[6]
            ligand_name = data[5]
            ligand_type = data[10]
            ligand_inchikey = data[20]
            receptor = data[2]
            receptor_family = data[1]
            receptor_class = data[0]
            receptor_sequence = data[21]
            ligand_status = data[22]
            pec50 = data[23]

            print(receptor)

            with open(output, 'a') as f:
                writer = csv.writer(f)
                writer.writerow(
                    [
                        ligand, ligand_name, ligand_type, ligand_status, pec50, ligand_inchikey, 
                        receptor, receptor_family, receptor_class, receptor_sequence]
                    )