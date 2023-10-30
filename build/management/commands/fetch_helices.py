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
        output = '/Users/vtk842/human_g_protein/gpcr_fasta_files/pp_TMs.csv'

        headers = 'gpcr,sequence,TM1,TM2,TM3,TM4,TM5,TM6,TM7,clas'.split(',')
        with open(output, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(headers)

        proteins = Protein.objects.filter(accession__isnull=False, species__common_name='Human').exclude(family__slug__startswith='100').exclude(family__slug__startswith='200')


        for protein in proteins:
            prot_name = protein.entry_name
            prot_seq = protein.sequence
            prot_clas = protein.family.parent.parent.parent.name.split(' ')[1]
            TM1 = self.amino_s(iterator='TM1', gpcr=prot_name)
            TM2 = self.amino_s(iterator='TM2', gpcr=prot_name)
            TM3 = self.amino_s(iterator='TM3', gpcr=prot_name)
            TM4 = self.amino_s(iterator='TM4', gpcr=prot_name)
            TM5 = self.amino_s(iterator='TM5', gpcr=prot_name)
            TM6 = self.amino_s(iterator='TM6', gpcr=prot_name)
            TM7 = self.amino_s(iterator='TM7', gpcr=prot_name)

            with open(output, 'a') as f:
                writer = csv.writer(f)
                writer.writerow([prot_name, prot_seq, TM1, TM2, TM3, TM4, TM5, TM6, TM7, prot_clas])
            print(prot_clas)
            # print(prot_name)
