#!/usr/bin/env nextflow

include { dnaseq } from './subworkflows/dnaseq/dnaseq.nf'
include { admixture } from './subworkflows/admixture/admixture.nf'

workflow {
	
    if (!params.vcf_ready) {
        merged_vcf = dnaseq()
    } else {
        merged_vcf = file(params.vcf_file)
    }

    admixture(merged_vcf)
}
