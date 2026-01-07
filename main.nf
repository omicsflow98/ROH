#!/usr/bin/env nextflow

include { dnaseq } from './subworkflows/dnaseq/dnaseq.nf'
include { ROH as ROHOM } from './subworkflows/ROH/ROH.nf'
include { ROH as ROHET } from './subworkflows/ROH/ROH.nf'
include { multiqc } from './processes/multiqc.nf'
include { plink as plinkhet } from './processes/plink.nf'
include { plink as plinkhom } from './processes/plink.nf'
include { detectruns } from './processes/detectruns.nf'
include { split_plink } from './processes/split_plink.nf'

workflow {
	
    mqc_config = file("${workflow.projectDir}/multiqc_config.yaml")
    dr_script = file("${projectDir}/scripts/detectRUNS.r")
    chromosomes= file( "${params.settings.chromsfile}" )
    split_results = file("${projectDir}/scripts/split_plink.py")
    ROHOM_out = channel.empty()
    ROHET_out = channel.empty()

    if (!params.main.vcf_ready) {
        merged_vcf = dnaseq()
    } else {
        merged_vcf = file(params.vcf.vcf_file)
    }
    
    if (params.settings.homozygosity) {
        homstate = "HOM"
        plinkhom(merged_vcf, params.plink.parameters, "")

        split_plink(split_results, plinkhom.out.roh_file, chromosomes)

        ROHOM_out = ROHOM(split_plink.out.separated_plink, homstate)
    }

    if (params.settings.heterozygosity) {
        
        chroms = channel.from(chromosomes.readLines())
        hetstate = "HET"
        plinkhet(merged_vcf, params.plink, chroms)
        detectruns(dr_script, plinkhet.out.dr_files, params.detectruns.parameters, plinkhet.out.chromosome)
            | collect
            | set {all_runs}

        ROHET_out = ROHET(all_runs, hetstate)
    }

    ROH_all = ROHOM_out.mix(ROHET_out).collect()

    multiqc(ROH_all, mqc_config)
}
