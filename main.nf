#!/usr/bin/env nextflow

include { dnaseq } from './subworkflows/dnaseq/dnaseq.nf'
include { ROH } from './subworkflows/ROH/ROH.nf'

workflow {
	
    if (!params.vcf_ready) {
        merged_vcf = dnaseq()
    } else {
        merged_vcf = file(params.vcf_file)
    }

    ROH(merged_vcf)
}
