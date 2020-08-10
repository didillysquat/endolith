#!/usr/bin/env python3

import pandas as pd
import os

class SplitSeqsByTaxa:
    def __init__(self):
        self.tax_dict = self._get_tax_dict()
        self.seq_seq_to_seq_represenative_name_dict = self._get_seq_dict()
        self.sample_fasta_dict, self.fasta_name = self._get_sample_fasta_dict()
        self.names_line_dict, self.names_name = self._get_name_objects()
        self.c_fasta_name = self.fasta_name.replace('.fasta', '.c.fasta')
        self.c_names_name = self.names_name.replace('.names', '.c.names')
        self.c_fasta = []
        self.c_names = []
        self.s_fasta_name = self.fasta_name.replace('.fasta', '.s.fasta')
        self.s_names_name = self.names_name.replace('.names', '.s.names')
        self.s_fasta = []
        self.s_names = []
        self.o_fasta_name = self.fasta_name.replace('.fasta', '.o.fasta')
        self.o_names_name = self.names_name.replace('.names', '.o.names')
        self.o_fasta = []
        self.o_names = []
        self._split()
        self._write()

    def _split(self):
        """
        Go through the fasta file and create split into three sets of new names and fastas
        For each sequnce, look up in the tax designation and assign accordingly
        """
        for seq_name, seq_seq in self.sample_fasta_dict.items():
            representative_seq = self.seq_seq_to_seq_represenative_name_dict[seq_seq]
            try:
                tax_designation = self.tax_dict[representative_seq]
            except KeyError:
                tax_designation = 'o'
            if tax_designation == 'c':
                self.c_fasta.append(f'>{seq_name}')
                self.c_fasta.append(seq_seq)
                self.c_names.append(self.names_line_dict[seq_name])
            elif tax_designation == 's':
                self.s_fasta.append(f'>{seq_name}')
                self.s_fasta.append(seq_seq)
                self.s_names.append(self.names_line_dict[seq_name])
            else:
                self.o_fasta.append(f'>{seq_name}')
                self.o_fasta.append(seq_seq)
                self.o_names.append(self.names_line_dict[seq_name])

    def _write(self):
        with open(self.c_fasta_name, 'w') as f:
            for line in self.c_fasta:
                f.write(f'{line}\n')
        with open(self.c_names_name, 'w') as f:
            for line in self.c_names:
                f.write(f'{line}\n')

        with open(self.s_fasta_name, 'w') as f:
            for line in self.s_fasta:
                f.write(f'{line}\n')
        with open(self.s_names_name, 'w') as f:
            for line in self.s_names:
                f.write(f'{line}\n')

        with open(self.o_fasta_name, 'w') as f:
            for line in self.o_fasta:
                f.write(f'{line}\n')
        with open(self.o_names_name, 'w') as f:
            for line in self.o_names:
                f.write(f'{line}\n')

    def _get_tax_dict(self):
        """
        We will process the taxonomyResult.tsv into a dict that related seq_name to either Symbiodiniaceae, Scleractinian or other

        """
        df = pd.read_table('taxonomyResult.tsv', index_col=0, header=0, names=['tax_id', 'tax_level', 'tax_identity', 'tax_levels'])
        tax_designations_dict = {}
        for seq_name, tax_level in df['tax_levels'].items():
            try:
                if 'o_Scleractinia' in tax_level:
                    tax_designations_dict[seq_name] = 'c'
                elif 'f_Symbiodiniaceae' in tax_level:
                    tax_designations_dict[seq_name] = 's'
                else:
                    tax_designations_dict[seq_name] = 'o'
            except TypeError:
                continue
        return tax_designations_dict

    def _get_seq_dict(self):
        with open('tax_query.fasta', 'r') as f:
            fasta_as_list = [_.rstrip() for _ in f]
        fasta_dict = {}
        for i in range(0, len(fasta_as_list), 2):
            fasta_dict[fasta_as_list[i+1]] = fasta_as_list[i][1:]
        return fasta_dict

    def _get_sample_fasta_dict(self):
        fasta_files = [_ for _ in os.listdir('.') if _.endswith('filtered.fasta')]
        assert(len(fasta_files) == 1)
        sample_fasta_name = fasta_files[0]
        with open(sample_fasta_name, 'r') as f:
            sample_fasta_as_list = [_.rstrip() for _ in f]
        fasta_dict = {}
        for i in range(0, len(sample_fasta_as_list), 2):
            fasta_dict[sample_fasta_as_list[i][1:]] = sample_fasta_as_list[i+1]
        return fasta_dict, sample_fasta_name

    def _get_name_objects(self):
        names_files = [_ for _ in os.listdir('.') if _.endswith('filtered.names')]
        assert(len(names_files) == 1)
        sample_name = names_files[0]
        with open(sample_name, 'r') as f:
            name_line_dict = {line.rstrip().split('\t')[0]: line.rstrip() for line in f}
        return name_line_dict, sample_name

SplitSeqsByTaxa()