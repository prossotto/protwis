from protein.models import Protein, ProteinCouplings

from django.core.management.base import BaseCommand

class Command(BaseCommand):

    def handle(self, *args, **options):
        all_bc = ''
        proteins = Protein.objects.filter(accession__isnull=False, species__common_name='Human').exclude(family__slug__startswith='100').exclude(family__slug__startswith='200')
        for p in proteins:
            couplings = ProteinCouplings.objects.filter(protein=p, source__in=['Bouvier'], g_protein_subunit__family__slug__startswith='100').order_by('g_protein_subunit')
            
            for c in couplings:
                print(f'{p.entry_name}-{c.g_protein_subunit}')
                all_bc += f'{p.entry_name}_{str(c.g_protein_subunit).upper()}\n'
                # print(f"Protein: {p.entry_name}, G Protein Subunit: {c.g_protein_subunit}, LogMaxEC50: {c.logmaxec50}")

        with open('/Users/vtk842/human_g_protein/new_jobs/bouvier_c.txt', 'w') as f:
            f.write(all_bc)
