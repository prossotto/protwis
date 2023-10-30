from protein.models import Protein, ProteinSegment, Residue, ProteinConformation
import csv
from django.core.management.base import BaseCommand

class Command(BaseCommand):

    def amino_s(self, iterator='', gpcr=''):
        sequence = ''
        sequence_pos = ''
        res = Residue.objects.filter(protein_conformation__protein=gpcr, protein_segment__slug=iterator).all()
        for residue in res:
            sequence += str(residue.amino_acid)
            sequence_pos += f'{str(residue.sequence_number)} '
        return sequence

    def handle(self, *args, **options):

        with open('/Users/vtk842/human_g_protein/new_jobs/changed_sequence/new_seqs.csv', 'w') as f:
                writer = csv.writer(f)
                writer.writerow(['gpcr', 'sequence', 'c-term', 'n-term', 'icl3'])

        names = '5ht2c_human, calcr_human, calrl_human, casr_human, drd3_human, drd4_human, gpr42_human, lpar2_human, lt4r2_human, opsb_human, opsr_human, t2r42_human'.split(', ')
        
        for name in names:
            gpcr = Protein.objects.filter(entry_name=name).first()
            gpcr_seq = gpcr.sequence
            c_term = self.amino_s(iterator='C-term', gpcr=gpcr)
            n_term = self.amino_s(iterator='N-term', gpcr=gpcr)
            icl3 = self.amino_s(iterator='ICL3', gpcr=gpcr)
                
            with open('/Users/vtk842/human_g_protein/new_jobs/changed_sequence/new_seqs.csv', 'a') as f:
                writer = csv.writer(f)
                writer.writerow([name, gpcr_seq, c_term, n_term, icl3])
            print(gpcr_seq)
            print(name)
        
