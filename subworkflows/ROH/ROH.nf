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

	overlap_map = overlap(overlap_script, individual_tsvs, bed_file)

	overlap_version = overlap_map.OL_version.collect().map{ it[0] }

	overlap_map.overlap_bed
		| collect
		| set { all_tsvs }

	combine_overlap_map = combine_overlap(all_tsvs, chromosomes)

	splitROH_map = splitROH(split_script, combine_overlap_map.ROH_file)

	splitROH_map.split_tsv
		| flatten()
		| set { individual_ROH }

	fixROH_map = fixROH(fix_script, individual_ROH)

	fixROH_version = fixROH_map.fixROH_version.collect().map{ it[0] }

	fixROH_map.fixed_tsv
		| collect
		| set { fixed_tsv }

	plot_results_map = plot_results(plot_script, fixed_tsv, state)

	ROH_results = plot_results_map.plots

	all_versions = overlap_version
		| combine(combine_overlap_map.combineOL_version)
		| combine(splitROH_map.splitROH_version)
		| combine(fixROH_version)
		| combine(plot_results_map.plots_version)

	emit:
	ROH_results
	all_versions
}
