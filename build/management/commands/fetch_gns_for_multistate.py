


class PdbStateIdentifier():
    def __init__(self, structure, tm2_gn='2x41', tm6_gn='6x38', tm3_gn='3x44', tm7_gn='7x52', inactive_cutoff=2, intermediate_cutoff=7.15):
        self.structure_type = None

        try:
            if structure.protein_conformation.protein.parent is None:
                raise Exception
            self.structure = structure
            self.structure_type = 'structure'
            family = structure.protein_conformation.protein.family
        except:
            try:
                structure.protein_conformation.protein
                self.structure = structure
                self.structure_type = 'refined'
                family = structure.protein_conformation.protein.family
            except:
                try:
                    structure.protein
                    self.structure = structure
                    self.structure_type = 'hommod'
                    family = structure.protein.family
                except:
                    structure.receptor_protein
                    self.structure = structure
                    self.structure_type = 'complex'
                    family = structure.receptor_protein.family
        if tm2_gn=='2x41' and tm6_gn=='6x38' and tm3_gn=='3x44' and tm7_gn=='7x52' and inactive_cutoff==2 and intermediate_cutoff==7.15:
            if family.slug.startswith('002') or family.slug.startswith('003'):
                tm6_gn, tm7_gn = '6x33', '7x51'
                inactive_cutoff, intermediate_cutoff = 2.5, 5.5
            elif family.slug.startswith('004'):
                inactive_cutoff, intermediate_cutoff = 5, 7.15
        self.tm2_gn, self.tm6_gn, self.tm3_gn, self.tm7_gn = tm2_gn, tm6_gn, tm3_gn, tm7_gn
        self.inactive_cutoff = inactive_cutoff
        self.intermediate_cutoff = intermediate_cutoff
        self.state = None
        self.activation_value = None
        self.line = False

    def run(self):
        if self.structure_type=='structure':
            self.parent_prot_conf = ProteinConformation.objects.get(protein=self.structure.protein_conformation.protein.parent)
            ssno = StructureSeqNumOverwrite(self.structure)
            ssno.seq_num_overwrite('pdb')
        elif self.structure_type=='refined':
            self.parent_prot_conf = ProteinConformation.objects.get(protein=self.structure.protein_conformation.protein)
        elif self.structure_type=='hommod':
            self.parent_prot_conf = ProteinConformation.objects.get(protein=self.structure.protein)
        elif self.structure_type=='complex':
            self.parent_prot_conf = ProteinConformation.objects.get(protein=self.structure.receptor_protein)
        # class A and T
        if self.parent_prot_conf.protein.family.slug.startswith('001') or self.parent_prot_conf.protein.family.slug.startswith('007'):
            tm6 = self.get_residue_distance(self.tm2_gn, self.tm6_gn)
            tm7 = self.get_residue_distance(self.tm3_gn, self.tm7_gn)
            print(tm6, tm7, tm6-tm7)
            if tm6 is not False and tm7 is not False:
                self.activation_value = tm6-tm7
                if self.activation_value<self.inactive_cutoff:
                    self.state = ProteinState.objects.get(slug='inactive')
                elif self.inactive_cutoff<=self.activation_value<=self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='intermediate')
                elif self.activation_value>self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='active')
        # class B
        elif self.parent_prot_conf.protein.family.slug.startswith('002') or self.parent_prot_conf.protein.family.slug.startswith('003'):
            tm2_gn_b = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm2_gn, scheme__short_name='GPCRdb(B)').label
            tm6_gn_b = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm6_gn, scheme__short_name='GPCRdb(B)').label
            tm3_gn_b = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm3_gn, scheme__short_name='GPCRdb(B)').label
            tm7_gn_b = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm7_gn, scheme__short_name='GPCRdb(B)').label

            tm6 = self.get_residue_distance(tm2_gn_b, tm6_gn_b)
            tm7 = self.get_residue_distance(tm3_gn_b, tm7_gn_b)
            if tm6 is not False and tm7 is not False:
                self.activation_value = tm6-tm7
                if self.activation_value<self.inactive_cutoff:
                    self.state = ProteinState.objects.get(slug='inactive')
                elif self.inactive_cutoff<=self.activation_value<=self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='intermediate')
                elif self.activation_value>self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='active')
        # class C
        elif self.parent_prot_conf.protein.family.slug.startswith('004'):
            tm2_gn_c = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm2_gn, scheme__short_name='GPCRdb(C)').label
            tm6_gn_c = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm6_gn, scheme__short_name='GPCRdb(C)').label
            tm3_gn_c = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm3_gn, scheme__short_name='GPCRdb(C)').label
            tm7_gn_c = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm7_gn, scheme__short_name='GPCRdb(C)').label

            tm6 = self.get_residue_distance(tm2_gn_c, tm6_gn_c)
            tm7 = self.get_residue_distance(tm3_gn_c, tm7_gn_c)
            if tm6 is not False and tm7 is not False:
                self.activation_value = tm6-tm7
                if self.activation_value<self.inactive_cutoff:
                    self.state = ProteinState.objects.get(slug='inactive')
                elif self.inactive_cutoff<=self.activation_value<=self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='intermediate')
                elif self.activation_value>self.intermediate_cutoff:
                    self.state = ProteinState.objects.get(slug='active')
        # class D
        #########
        # class F
        elif self.parent_prot_conf.protein.family.slug.startswith('006'):
            tm2_gn_f = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm2_gn, scheme__short_name='GPCRdb(F)').label
            tm6_gn_f = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm6_gn, scheme__short_name='GPCRdb(F)').label
            tm3_gn_f = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm3_gn, scheme__short_name='GPCRdb(F)').label
            tm7_gn_f = ResidueGenericNumberEquivalent.objects.get(default_generic_number__label=self.tm7_gn, scheme__short_name='GPCRdb(F)').label

            tm6 = self.get_residue_distance(tm2_gn_f, tm6_gn_f)
            tm7 = self.get_residue_distance(tm3_gn_f, tm7_gn_f)
            if tm6 is not False and tm7 is not False:
                self.activation_value = tm6-tm7
                if self.activation_value<0:
                    self.state = ProteinState.objects.get(slug='inactive')
                elif 0<=self.activation_value<=2:
                    self.state = ProteinState.objects.get(slug='intermediate')
                elif self.activation_value>2:
                    self.state = ProteinState.objects.get(slug='active')
        else:
            print('{} is not class A,B,C,F'.format(self.structure))
        if self.structure_type=='structure':
            ssno.seq_num_overwrite('pdb')

    def get_residue_distance(self, residue1, residue2):
        try:
            res1 = Residue.objects.get(protein_conformation__protein=self.structure.protein_conformation.protein.parent, display_generic_number__label=dgn(residue1, self.parent_prot_conf))
            res2 = Residue.objects.get(protein_conformation__protein=self.structure.protein_conformation.protein.parent, display_generic_number__label=dgn(residue2, self.parent_prot_conf))
            print(res1, res1.id, res2, res2.id)
            try:
                rota1 = Rotamer.objects.filter(structure=self.structure, residue__sequence_number=res1.sequence_number)
                if len(rota1)==0:
                    raise Exception
            except:
                rota1 = Rotamer.objects.filter(structure=self.structure, residue__display_generic_number__label=dgn(residue1, self.structure.protein_conformation))
            rota1 = right_rotamer_select(rota1, self.structure.preferred_chain[0])
            try:
                rota2 = Rotamer.objects.filter(structure=self.structure, residue__sequence_number=res2.sequence_number)
                if len(rota2)==0:
                    raise Exception
            except:
                rota2 = Rotamer.objects.filter(structure=self.structure, residue__display_generic_number__label=dgn(residue2, self.structure.protein_conformation))
            rota2 = right_rotamer_select(rota2, self.structure.preferred_chain[0])
            rotas = [rota1, rota2]
            io1 = StringIO(rotas[0].pdbdata.pdb)
            rota_struct1 = PDB.PDBParser(QUIET=True).get_structure('structure', io1)[0]
            io2 = StringIO(rotas[1].pdbdata.pdb)
            rota_struct2 = PDB.PDBParser(QUIET=True).get_structure('structure', io2)[0]

            for chain1, chain2 in zip(rota_struct1, rota_struct2):
                for r1, r2 in zip(chain1, chain2):
                    # print(self.structure, r1.get_id()[1], r2.get_id()[1], self.calculate_CA_distance(r1, r2), self.structure.state.name)
                    line = '{},{},{},{},{}\n'.format(self.structure, self.structure.state.name, round(self.calculate_CA_distance(r1, r2), 2), r1.get_id()[1], r2.get_id()[1])
                    self.line = line
                    return self.calculate_CA_distance(r1, r2)
        except:
            try:
                res1 = Residue.objects.get(protein_conformation=self.parent_prot_conf, display_generic_number__label=dgn(residue1, self.parent_prot_conf))
                res2 = Residue.objects.get(protein_conformation=self.parent_prot_conf, display_generic_number__label=dgn(residue2, self.parent_prot_conf))
                if self.structure_type=='refined':
                    pdb_data = self.structure.pdb_data.pdb
                elif self.structure_type=='hommod':
                    pdb_data = self.structure.pdb_data.pdb
                io = StringIO(pdb_data)
                struct = PDB.PDBParser(QUIET=True).get_structure('structure', io)[0]
                for chain in struct:
                    r1 = chain[res1.sequence_number]
                    r2 = chain[res2.sequence_number]
                    print(self.structure, r1.get_id()[1], r2.get_id()[1], self.calculate_CA_distance(r1, r2), self.structure.state.name)
                    line = '{},{},{},{},{}\n'.format(self.structure, self.structure.state.name, round(self.calculate_CA_distance(r1, r2), 2), r1.get_id()[1], r2.get_id()[1])
                    self.line = line
                    return self.calculate_CA_distance(r1, r2)

            except:
                print('Error: {} no matching rotamers ({}, {})'.format(self.structure.pdb_code.index, residue1, residue2))
                return False

    def calculate_CA_distance(self, residue1, residue2):
        diff_vector = residue1['CA'].get_coord()-residue2['CA'].get_coord()
        return numpy.sqrt(numpy.sum(diff_vector * diff_vector))


