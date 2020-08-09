#!/usr/bin/env nextflow
params.input_data_path = "${params.base_dir}/input_data"
params.path_to_seqtaxdb = "/home/humebc/nt/nt.fnaDB"
Channel.fromFilePairs("${params.input_data_path}/*{1,2}.fastq.gz").map{it[1]}.into{ch_input_file_pairs_batch; ch_input_file_pairs_run}

// Write out the mothur batch files that will be used to:
// make.contigs
// screen.seqs maxambig=0
// pcr.seqs pdiffs=2, rdiffs=2
// We wil ignore the errors as the negative samples fail and this is fine. 
// Presumably they will simply be missing from the output channel
process write_out_mothur_batch_fies{
    cache 'lenient'
    tag "${file_one}"
    conda "envs/nf_18S.yml"
    storeDir "${params.base_dir}/mothur_qc/mothur_out"
    errorStrategy 'ignore'

    input:
    tuple file(file_one), file(file_two) from ch_input_file_pairs_batch

    output:
    tuple file("${file_one.getName().replaceAll('.gz', '.trim.contigs.good.unique.pcr.unique.fasta')}"), file("${file_one.getName().replaceAll('.gz', '.trim.contigs.good.unique.pcr.names')}") into ch_mothur_out, ch_tax_screen_in
    // tuple file("${file_one.getName().replaceAll('R1.fastq.gz', 'mothur_batch')}"), file("${file_one.getName().replaceAll('R1.fastq.gz', 'primers.oligos')}") into ch_mothur_batch

    script:
    mothur_batch_file_name = file_one.getName().replaceAll('R1.fastq.gz', 'mothur_batch')
    """
    python3 ${params.base_dir}/bin/write_out_mothur_batch_files.py $file_one
    mothur $mothur_batch_file_name
    """
}

// After the mothur QC we will have pairs of .name .fasta files
// We can come back to doing a size screening at a later date.
// For the time being we will proceed straight onto the taxonomic screening
// Gather all distinct sequences into a single and fasta and run this against the nt database using mmseqs
// Do tax screening of the results to find out which sequences are the keepers.

// Pseudo-code for taxonomy. We end with .names and .fasta files
// Python script that makes a distinct sequence fasta
process make_distinct_query_fasta_run_mmseqs{
    cache 'lenient'
    tag "run_taxonomy"
    conda "envs/nf_18S.yml"
    // storeDir "${params.base_dir}/taxonomy_QC/mmseqs_input_fasta"
    
    input:
    file(pairs) from ch_tax_screen_in.collect()

    output:
    file "*.filtered.{fasta,names}" into ch_filtered_files
    //"mmseq.out" into ch_mmseq_out

    script:
    """
    python3 ${params.base_dir}/bin/make_mmseqs_in_fasta.py
    """
}

// for the time being, let's try to get output channel into shape

// At this point we'll want to work with the filtered files
// It would be good if we could go back to paired files from the collect
ch_filtered_files.toSortedList().flatten().groupBy { String str -> str[0..3] }.println()
// fasta is input to process that runs against an mmseqs database
// output is passed into python script that splits the .fasta and .names files into coral symbiodiniaceae and other
// process mmseq_fasta_against_nt{

// }

// mmseqs createdb tax_query.fasta queryDB
// mmseqs taxonomy queryDB ${params.path_to_seqtaxdb} taxonomyResult tmp
// mmseqs createtsv queryDB taxonomyResult taxonomyResult.tsv
// mmseqs taxonomyreport seqTaxDB taxonomyResult report.html --report-mode 1




