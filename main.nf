#!/usr/bin/env nextflow

include { ROH as ROHOM } from './subworkflows/ROH/ROH.nf'
include { ROH as ROHET } from './subworkflows/ROH/ROH.nf'
include { quarto } from './processes/quarto.nf'
include { plinkhet } from './processes/plink.nf'
include { plink as plinkhom } from './processes/plink.nf'
include { plink as plinkhetbim } from './processes/plink.nf'
include { detectruns } from './processes/detectruns.nf'
include { gtftobed } from './processes/bedops.nf'
include { split_plink } from './processes/split_plink.nf'
include { all_versions } from './processes/MergeVersions.nf'
include { vcftools } from './processes/vcftools.nf'

workflow {
	
    dr_script = file("${projectDir}/scripts/detectRUNS.r")
    chromosomes= file( "${params.settings.chromsfile}" )
    split_results = file("${projectDir}/scripts/split_plink.py")
    version_script = file("${projectDir}/scripts/version_merge.py")
    pipeline_version = "1.1.0"
    ROHOM_out = channel.empty()
    ROHET_out = channel.empty()
    ROHOM_indiv = channel.empty()
    ROHET_indiv = channel.empty()
    ROHOM_islands = channel.empty()
    ROHET_islands = channel.empty()
    dr_params = channel.value(params.detectruns.parameters)
    plink_params = channel.value(params.plink.parameters)
    quarto_files = channel.of(
    tuple(
        file("${projectDir}/scripts/Quarto/demo.qmd"),
        file("${projectDir}/scripts/Quarto/test.R"),
        file("${projectDir}/scripts/Quarto/subpopulations.R"),
        file("${projectDir}/scripts/Quarto/style.scss")
        )
    )
    island_info = channel.value(
        tuple (
            params.island.type,
            params.island.minSNP,
            params.island.indiv_count
        )
    )
    subpop_tsv = params.settings.subpopulations? channel.fromPath(params.subpopulation.subpop_file):channel.of([])

    if (params.settings.genefile == "gtf") {
        def genefile = file("${params.genepath.gtf}")
        convertgtf = gtftobed(genefile)
        bed_file = convertgtf.bedfile
    } else {
        bed_file = file("${params.genepath.bed}")
    }

    if (params.inputfiles.start_vcf) {
        infile = channel.value(
            tuple(
                'VCF',
                file(params.inputfiles.vcf_file),
                [],
                [],
                []
            )
        )
        vcf_info = vcftools(file(params.inputfiles.vcf_file))
    } else {
        prefix = params.inputfiles.plink_file
        def bed = file("${prefix}.bed")
        def bim = file("${prefix}.bim")
        def fam = file("${prefix}.fam")
            
        infile = channel.value(
            tuple(
                'plink',
                [],
                bed,
                bim,
                fam
            )
        )
        vcf_info = [ vcf_info: file("${prefix}.fam") ]
    }

    if (params.settings.homozygosity) {
        homstate = "HOM"

        plinkhom_map = plinkhom(infile, plink_params, "")

        split_plink_map = split_plink(split_results, plinkhom_map.roh_file, chromosomes)

        ROHOM_map = ROHOM(split_plink_map.separated_plink, plinkhom_map.bim_file, bed_file, homstate)
        ROHOM_out = ROHOM_map.ROH_results
        ROHOM_version = ROHOM_map.all_versions
        ROHOM_indiv = ROHOM_map.indiv_info
        ROHOM_islands = ROHOM_map.finalislands
    }

    if (params.settings.heterozygosity) {
        
        chroms = channel.fromList(chromosomes.readLines())
        hetstate = "HET"
        plinkhet_map = plinkhet(infile, chroms)
        plinkhet_version = plinkhet_map.plink_version.collect().map{ it[0] }

        if (params.inputfiles.start_vcf) {
            make_bim = plinkhetbim(infile, plink_params, "")
        } else {
            make_bim = [ bim_file: file("${prefix}.bim") ]
        }

        detectruns_map = detectruns(dr_script, plinkhet_map.dr_files, dr_params, plinkhet_map.chromosome)

        detectruns_map.dr_tsv
            | collect
            | set {all_runs}

        detectruns_version = detectruns_map.dr_version.collect().map{ it[0] }

        ROHET_map = ROHET(all_runs, make_bim.bim_file, bed_file, hetstate)
        ROHET_out = ROHET_map.ROH_results
        ROHET_version = ROHET_map.all_versions
        ROHET_indiv = ROHET_map.indiv_info
        ROHET_islands = ROHET_map.finalislands
    }

    ROH_all = ROHOM_out.mix(ROHET_out).collect()
    ROH_island = ROHOM_islands.mix(ROHET_islands).collect()
    indiv_all = ROHOM_indiv.mix(ROHET_indiv).collect()

    quarto_map = quarto(ROH_all,
    ROH_island, 
    indiv_all, 
    vcf_info.vcf_info, 
    quarto_files, 
    params.settings.homozygosity, 
    params.settings.heterozygosity, 
    params.settings.subpopulations, 
    params.inputfiles.start_vcf,
    subpop_tsv, 
    bed_file,
    island_info )

    all_versions_ch = quarto_map.quarto_version
//    | combine(vcf_info.vcft_version)

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
