#!/usr/bin/env nextflow

include { overlap } from '../../processes/overlap.nf'
include { combine_overlap } from '../../processes/overlap.nf'
include { plot_results } from '../../processes/plots.nf'
include { splitROH } from '../../processes/splitROH.nf'
include { fixROH } from '../../processes/fixROH.nf'

workflow ROH {

	take:
	collected_tsvs
	state

	main:

	chromosomes= file( "${params.settings.chromsfile}" )
	bed_file = file("${params.settings.bedfile}")
	overlap_script = file("${projectDir}/scripts/overlap_roh.py")
	plot_script = file("${projectDir}/scripts/plots.R")
	split_script = file("${projectDir}/scripts/splitROH.R")
	fix_script = file("${projectDir}/scripts/fixROH.R")

	collected_tsvs
		| flatten()
		| set { individual_tsvs }

	overlap(overlap_script, individual_tsvs, bed_file)
		| collect
		| set { all_tsvs }

	combine_overlap(all_tsvs, chromosomes)

	splitROH(split_script, combine_overlap.out.ROH_file)
		| flatten()
		| set { individual_ROH }

	fixROH(fix_script, individual_ROH)
		| collect
		| set { fixed_tsv }

	plot_results(plot_script, fixed_tsv, state)

	ROH_results = plot_results.out.plots

	emit:
	ROH_results
}
