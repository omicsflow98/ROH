#!/usr/bin/env nextflow

include {runfastqc} from '../../processes/fastqc.nf'
include {trim_galore} from '../../processes/trim_galore.nf'
include {bwa} from '../../processes/bwa.nf'
include {markduplicates} from '../../processes/markdup.nf'
include {BQSR} from '../../processes/bqsr.nf'

workflow alignment {

main:
first_ch = Channel.fromPath(params.data_csv, checkIfExists: true)
	| splitCsv(header: true, sep: '\t')
	| map { row -> tuple(row.SampName,
			file(row.File1),
			file(row.File2),
			row.Adapter1,
			row.Adapter2,
			row.LibName,
			row.Barcode,
			row.Platform) }

bwa_index = file( "${params.reference}/index/bwa/genome.fa.{,amb,ann,bwt,pac,sa}" )

runfastqc(first_ch)

if (params.trim) {
	trimmed_reads = trim_galore(first_ch)
} else {
	trimmed_reads = Channel.fromPath(params.data_csv, checkIfExists: true)
	| splitCsv(header: true, sep: '\t')
	| map { row -> tuple(row.SampName,
			row.LibName,
			row.Barcode,
			row.Platform,
			[file(row.File1), file(row.File2)]) 
			}
}

bwa(bwa_index, trimmed_reads)

markduplicates(params.temp_dir, bwa.out.bam_files)

if (params.recalibrate) {
	BQSR(params.reference, params.temp_dir, markduplicates.out.mark_dup)
	bam_files = BQSR.out.bam_bqsr
} else {
	bam_files = markduplicates.out.mark_dup
}

emit:
bam_files

}