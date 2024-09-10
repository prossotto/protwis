from django.core.management.base import BaseCommand
from django.conf import settings

from protein.models import Protein, Gene
from residue.models import Residue
from structure.models import Structure
from structure.functions import ParseStructureCSV
from tools.management.commands.build_structure_angles import NonHetSelect
from contactnetwork.interaction import InteractingPair
from django.db.models import Count, Q, Prefetch, TextField, Avg, Case, When, CharField, Subquery, OuterRef



class Command(BaseCommand):
    help = '''Adds segment end annotations for new structures from structure related csv files into the yaml files. 
              Run it on a complete db after parse_excel_annotations.py was run and the csv files got updated with the 
              updated GPCRdb_structure_info.xlsx'''


    def handle(self, *args, **options):

        # StructureModel
        filter_structures = []
        filter_structures.extend(['x-ray-diffraction', 'electron-microscopy', 'electron-crystallography'])

        # Subquery definition
        gene_subquery_models = Gene.objects.filter(proteins=OuterRef('pk')).values('name')[:1]
        # gene_subquery_structs = Gene.objects.filter(proteins=OuterRef('protein_conformation__protein__pk')).values('name')[:1]
        gene_subquery_structs = Gene.objects.filter(proteins=OuterRef('protein_conformation__protein__pk')).values('name')[:1]


        # Annotating the Protein queryset
        prots = Protein.objects.filter(id=181).annotate(gene_name=Subquery(gene_subquery_models))
        structs = Structure.objects.all(structure_type__slug__in=filter_structures).annotate(gene_name=Subquery(gene_subquery_structs))



        proti = Protein.objects.get(id=181)
        gene = Gene.objects.get(proteins=proti)

        for struct in structs:
            # if struct.gene_name == None:
            print(struct.protein_conformation)
            print(struct.gene_name)

        # for prot in prots:
        #     print('From Subquery of Protein objects')
        #     print(prot.gene_name)

        # print('From GET')
        # print(gene)