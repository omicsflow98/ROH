#!/usr/bin/env nextflow

include {gvcf} from '../../processes/VariantCalling.nf'
include {mergedgvcf} from '../../processes/VariantCalling.nf'

workflow variant_calling {
	take:
	bam_files

	main:

	gvcf(params.reference, params.temp_dir, bam_files)

	gvcf.out.vcf_files
		| collect
		| set {collected_vcf}

	mergedgvcf(params.reference, params.temp_dir, collected_vcf)

	emit:
	vcf_file = mergedgvcf.out.merged_vcf

}