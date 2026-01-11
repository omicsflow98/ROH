#!/usr/bin/env nextflow

include { dnaseq } from './subworkflows/dnaseq/dnaseq.nf'
include { ROH as ROHOM } from './subworkflows/ROH/ROH.nf'
include { ROH as ROHET } from './subworkflows/ROH/ROH.nf'
include { multiqc } from './processes/multiqc.nf'
include { plink as plinkhet } from './processes/plink.nf'
include { plink as plinkhom } from './processes/plink.nf'
include { detectruns } from './processes/detectruns.nf'
include { split_plink } from './processes/split_plink.nf'
include { all_versions } from './processes/MergeVersions.nf'

workflow {
	
    mqc_config = file("${workflow.projectDir}/multiqc_config.yaml")
    dr_script = file("${projectDir}/scripts/detectRUNS.r")
    chromosomes= file( "${params.settings.chromsfile}" )
    split_results = file("${projectDir}/scripts/split_plink.py")
    version_script = file("${projectDir}/scripts/version_merge.py")
    pipeline_version = "1.1.0"
    ROHOM_out = channel.empty()
    ROHET_out = channel.empty()

    if (!params.main.vcf_ready) {
        merged_vcf = dnaseq()
    } else {
        merged_vcf = file(params.vcf.vcf_file)
    }
    
    if (params.settings.homozygosity) {
        homstate = "HOM"

        plinkhom_map = plinkhom(merged_vcf, params.plink.parameters, "")

        split_plink_map = split_plink(split_results, plinkhom_map.roh_file, chromosomes)

        ROHOM_map = ROHOM(split_plink_map.separated_plink, homstate)
        ROHOM_out = ROHOM_map.ROH_results
        ROHOM_version = ROHOM_map.all_versions
    }

    if (params.settings.heterozygosity) {
        
        chroms = channel.from(chromosomes.readLines())
        hetstate = "HET"
        plinkhet_map = plinkhet(merged_vcf, params.plink, chroms)
        plinkhet_version = plinkhet_map.plink_version.collect().map{ it[0] }

        detectruns_map = detectruns(dr_script, plinkhet_map.dr_files, params.detectruns.parameters, plinkhet_map.chromosome)

        detectruns_map.dr_tsv
            | collect
            | set {all_runs}

        detectruns_version = detectruns_map.dr_version.collect().map{ it[0] }

        ROHET_map = ROHET(all_runs, hetstate)
        ROHET_out = ROHET_map.ROH_results
        ROHET_version = ROHET_map.all_versions
    }

    ROH_all = ROHOM_out.mix(ROHET_out).collect()

    multiqc_map = multiqc(ROH_all, mqc_config)

    all_versions_ch = multiqc_map.multiqc_version

    if (params.settings.homozygosity) {
        all_versions_ch = all_versions_ch
        | combine(ROHOM_version)
        | combine(plinkhom_map.plink_version)
        | combine(split_plink_map.splitplink_version)
    }
    if (params.settings.heterozygosity) {
        all_versions_ch = all_versions_ch
        | combine(ROHET_version)
        | combine(plinkhet_version)
        | combine(detectruns_version)
    }

    all_versions(version_script, all_versions_ch, pipeline_version)
}
