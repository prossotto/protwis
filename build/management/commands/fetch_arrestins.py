from protein.models import Protein, ProteinSegment, Residue, ProteinConformation
import csv
from django.core.management.base import BaseCommand

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
        output = '/Users/vtk842/human_g_protein/gpcr_fasta_files/pp_arrestins.csv'

        headers = 'arrestin,sequence'.split(',')
        with open(output, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(headers)

        proteins = Protein.objects.filter(accession__isnull=False, species__common_name='Human', family__slug__startswith='200')


        for protein in proteins:
            prot_name = protein.entry_name
            prot_seq = protein.sequence
            prot_clas = protein.family.parent.parent.parent.name

            with open(output, 'a') as f:
                writer = csv.writer(f)
                writer.writerow([prot_name, prot_seq,prot_clas])

            print(f'{prot_name}\n{prot_clas}\n')
            # print(prot_name)
