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

        with open('/Users/vtk842/human_g_protein/gpcr_fasta_files/CN-term_seq/new_CNc.csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['gpcr', 'sequence', 'c-term', 'n-term', 'icl3', 'clas'])

        proteins = Protein.objects.filter(accession__isnull=False, species__common_name='Human').exclude(family__slug__startswith='100').exclude(family__slug__startswith='200')

        for protein in proteins:
            prot_name = protein.entry_name
            prot_seq = protein.sequence
            prot_clas = protein.family.parent.parent.parent.name.split(' ')[1]
            c_term = self.amino_s(iterator='C-term', gpcr=prot_name)
            n_term = self.amino_s(iterator='N-term', gpcr=prot_name)
            icl3 = self.amino_s(iterator='ICL3', gpcr=prot_name)
            

            with open('/Users/vtk842/human_g_protein/gpcr_fasta_files/CN-term_seq/new_CNc.csv', 'a') as f:
                writer = csv.writer(f)
                writer.writerow([prot_name, prot_seq, c_term, n_term, icl3, prot_clas])
            print(prot_clas)
            # print(prot_name)

            
        
            
