#!/usr/bin/env nextflow

include { plink } from '../../processes/plink.nf'
include { split_plink } from '../../processes/split_plink.nf'
include { overlap } from '../../processes/overlap.nf'
include { combine_overlap } from '../../processes/overlap.nf'

workflow ROH {

	take:
	merged_vcf

	main:

	chromosomes= file( "${params.chromsfile}" )
	bed_file = file("${params.bedfile}")
	split_results = file("${projectDir}/scripts/split_plink.py")
	overlap_script = file("${projectDir}/scripts/overlap_roh.py")

	plink(merged_vcf)

	split_plink(split_results, plink.out.roh_file, chromosomes)
		.flatten()
		.set { individual_tsvs }

	overlap(overlap_script, individual_tsvs, bed_file)
	| collect
	| set { all_tsvs }

	combine_overlap(all_tsvs)
}
