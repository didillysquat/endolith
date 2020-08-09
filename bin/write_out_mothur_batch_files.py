#!/usr/bin/env python3

import sys

class MothurBatchOutput:
    """
    Create and write out a batch file that can be used by mothur to run the QC of one sample.
    Command line input will give us the directory where the batch files should be written out to
    Comand line will also give us the fwd read name.
    """
    def __init__(self):
        self.read_one_name = sys.argv[1]
        self.read_two_name = self.read_one_name.replace('R1', 'R2')
        self._write_out_oligo_seqs()
        self.mothur_batch_file = [
            f'make.contigs(ffastq={self.read_one_name}, rfastq={self.read_two_name}, processors=1)',

            f'screen.seqs(fasta={self.read_one_name.replace(".gz", ".trim.contigs.fasta")}, maxambig=0, processors=1)',

            f'unique.seqs(fasta={self.read_one_name.replace(".gz", ".trim.contigs.good.fasta")})',

            f'pcr.seqs(fasta={self.read_one_name.replace(".gz", ".trim.contigs.good.unique.fasta")}, '
            f'name={self.read_one_name.replace(".gz", ".trim.contigs.good.names")}, '
            f'oligos={self.read_one_name.replace("R1.fastq.gz", "primers.oligos")}, pdiffs=2, rdiffs=2, processors=1)',

            f'unique.seqs(fasta={self.read_one_name.replace(".gz", ".trim.contigs.good.unique.pcr.fasta")}, '
            f'name={self.read_one_name.replace(".gz", ".trim.contigs.good.pcr.names")})'
        ]

    def write_out_batch(self):
        with open(f'{self.read_one_name.replace("R1.fastq.gz", "mothur_batch")}', 'w') as f:
            for line in self.mothur_batch_file:
                f.write(f'{line}\n')

    def _write_out_oligo_seqs(self):
        """
        Write out the primers.oligos file
        """
        oligo_file = ['forward\tTTGTACACACCGCCC', 'reverse\tCCTTCYGCAGGTTCACCTAC']
        with open(f'{self.read_one_name.replace("R1.fastq.gz", "primers.oligos")}', 'w') as f:
            for line in oligo_file:
                f.write(f'{line}\n')

mbo = MothurBatchOutput()
mbo.write_out_batch()