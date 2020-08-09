#!/usr/bin/env python3

import os
from collections import defaultdict
import numpy as np

class MakeMmseqsInFasta:
    """
    Go through all of the .fasta files and create a count dict of the sequences where we count the
    number of samples a given sequence was found in
    Go back through all of the .fasta .name files and remove any sequences that are not found below a 1% abundance
    and not found in more than 1 sample.
    Also perform a size crop.
    Output from this should be a query fasta that we will use as input to the mmseqs of blast 
    and a new .fasta .names pair for each of the extant pairs
    If a pair no longer exists then simply don't output anything
    """
    def __init__(self):
        self.count_dict = defaultdict(int)
        # set of the sequences that we should make a fasta from to do the blast or mmseq from
        self.query_set = set()
        self.upper, self.lower = self._first_count()
        self._screening()
        self._write_out_query_fasta()

    def _first_count(self):
        for fasta_file in [_ for _ in os.listdir('.') if _.endswith('.fasta')]:
            print(f'Counting {fasta_file}')
            with open(fasta_file, 'r') as f:
                fasta_file_as_list = [_.rstrip() for _ in f]
            
            for line in fasta_file_as_list:
                if not line.startswith('>'):
                    self.count_dict[line] += 1
        
        # At this point we have the count dict populated
        # Now let's work out a size cutoff
        # Let's work this out as whatever is larger, 30 bp +- the average
        # or the 5% and 95% percentile
        sizes = np.array([len(_) for _ in self.count_dict.keys()])
        fith = int(np.percentile(sizes, 5))
        average = int(np.mean(sizes))
        ninety_fith = int(np.percentile(sizes, 95))
        lower = min((average-30), fith)
        upper = max((average+30), ninety_fith)
        return upper, lower
    
    def _screening(self):
        # Now when we do the second pass of the .fasta files we can screen by size and sample count
        
        for fasta_file in [_ for _ in os.listdir('.') if _.endswith('.fasta')]:
            new_fasta = []
            new_names = []
            # we will make an abundance dict out of the names file
            names_file = fasta_file.replace('.unique.fasta', '.names')
            with open(names_file, 'r') as f:
                name_file_list = [_.rstrip() for _ in f]
            abs_abund_dict = {line.split('\t')[0]: len(line.split('\t')[1].split(',')) for line in name_file_list}
            name_line_dict = {line.split('\t')[0]: line for line in name_file_list}
            with open(fasta_file, 'r') as f:
                fasta_file_list = [_.rstrip() for _ in f]
            
            fasta_dict = {}
            for i in range(0,len(fasta_file_list),2):
                fasta_dict[fasta_file_list[i][1:].split('\t')[0]] = fasta_file_list[i+1]

            tot = sum(abs_abund_dict.values())
            for seq_name, seq_sequence in fasta_dict.items():
                # If any of the below are true, the seq will not make it into the new fasta and name files
                if len(seq_sequence) < self.lower:
                    continue
                elif len(seq_sequence) > self.upper:
                    continue
                elif (abs_abund_dict[seq_name]/tot) < 0.01 and self.count_dict[seq_sequence] == 1:
                    continue
                new_fasta.append(f'>{seq_name}')
                new_fasta.append(seq_sequence)
                new_names.append(name_line_dict[seq_name])
                # this is one of the sequences that should be in the query fasta
                self.query_set.add(seq_sequence)
            
            # At this point we have a new fasta and names file populated
            # write them out
            new_fasta_name = fasta_file.replace('.fasta', '.filtered.fasta')
            new_names_name = names_file.replace('.names', '.filtered.names')
            if new_fasta:
                with open(new_fasta_name, 'w') as f:
                    for line in new_fasta:
                        f.write(f'{line}\n')
                with open(new_names_name, 'w') as f:
                    for line in new_names:
                        f.write(f'{line}\n')

    def _write_out_query_fasta(self):
        """
        write out the self.query set as a fasta
        """
        with open('tax_query.fasta', 'w') as f:
            for i, seq in enumerate(self.query_set):
                f.write(f'>seq_{i}\n')
                f.write(f'{seq}\n')
    
mmif = MakeMmseqsInFasta()
            


            
            



        
        