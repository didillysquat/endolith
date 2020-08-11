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
    tag "make_tax_query_db"
    conda "envs/nf_18S.yml"
    
    input:
    file(pairs) from ch_tax_screen_in.collect()

    output:
    file "*.filtered.{fasta,names}" into ch_filtered_files
    file "tax_query.fasta" into ch_mmseq_out

    script:
    """
    python3 ${params.base_dir}/bin/make_mmseqs_in_fasta.py
    """
}

// Make the query db from the fasta that was ouput from prev process
// Use the queryDB made in prev process to do an actual taxonomy query
// NB it was a lot of effort to make the mmseqs taxonomy database
// We started from from the nt database that had already been downloaded using
// the ususal update_blastdb.pl
// We then followed the docs here:
// https://github.com/soedinglab/MMseqs2/wiki#create-a-seqtaxdb-from-an-existing-blast-database
// We then made an index of the resulting nt.fnaDB
// This index was absolutely massive though - about 4T so it may be worth not doing the index in future
// In general the massive file sizes could be a down side to using mmseqs in general and it could be that
// blast will be a more realistic way forwards. But for the time being will try to work with
// mmseqs as we've put the effort in so far.
process make_mmseqs_querydb{
    cache 'lenient'
    tag 'make_mmseqs_querydb'
    conda "envs/nf_18S.yml"
    publishDir "${params.base_dir}/tax_results"

    input:
    file in_fasta from ch_mmseq_out

    output:
    tuple file('taxonomyResult.tsv'), file('tax_query.fasta') into ch_tax_results

    script:
    """
    mmseqs createdb $in_fasta queryDB
    mmseqs taxonomy queryDB ${params.path_to_seqtaxdb} taxonomyResult tmp --tax-lineage 1
    mmseqs createtsv queryDB taxonomyResult taxonomyResult.tsv
    """
}


// Once we have the taxonomyResult file we can then go on a pair by pair basis and split the .fasta and .name files
// into scleractinia, symbiodiniaceae and other 
// The channel manipulation was relatively complex. Becuase we had done a collect() for the make_distinct_query_fasta_run_mmseqs
// process, we had to work out how to pairs this back up i.e. one tuple containing the .fasta and .names pair.
// To do the pairing we used .groupBy, and used the first 4 characters as the key. Importantly,
// this produced a single list hash map object that we then extracted the keys from, converted to list and then importantly
// called with flatMap, so that the single list was returned as individual tuples.
// Finally, we wanted to have a copy of the taxonomyResults.tsv and the tax_query.fasta for each of the fasta/names pairs.
// To do this we were able to use combine. This worked well.
process split_seqs_by_taxa{
    cache 'lenient'
    tag "${fasta_file}"
    conda "envs/nf_18S.yml"
    publishDir "${params.base_dir}/post_qc_seqs"

    input:
    tuple file(names_file), file(fasta_file), file(tax_tsv), file(query_fasta) from ch_filtered_files.toSortedList().flatten().groupBy { file_obj -> file_obj.getName()[0..3] }.flatMap{it.values().toList()}.combine(ch_tax_results)

    output:
    tuple file('*.c.fasta'), file('*.c.names'), file('*.s.fasta'), file('*.s.names'), file('*.o.fasta'), file('*.o.names') into ch_plotting

    script:
    """
    python3 ${params.base_dir}/bin/split_seqs_by_taxa.py
    """
}

// At this point we have the filtered reads split by taxa. Now would be a good time to start analysing the read data and plotting
// Let's work towards doing plots that are similar to those that we did for the TARA 18S. I.e. we should start by doing a plot
// that is staked bar of the sequences. We should complete this roughly and quickly so that we can judge how to move forwards.
process plot_stacked_bars_all_seqs{
    cache false
    tag "plot_stacked_bars_all_seqs"
    conda "envs/nf_18S.yml"
    publishDir "${params.base_dir}/plotting"

    input:
    tuple file(c_fasta), file(c_names), file(s_fasta), file(s_names), file(o_fasta), file(o_names) from ch_plotting

    output:
    tuple file('*.svg'), file('*.png') into ch_plotting_out

    script:
    """
    python3 ${params.base_dir}/bin/plotting.py
    """
}