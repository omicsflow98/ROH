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
	
    //Pipeline version
    pipeline_version = "1.1.0"

    // Read in the scripts and files needed for the workflow
    dr_script = file("${projectDir}/scripts/detectRUNS.r")
    chromosomes= file( "${params.settings.chromsfile}" )
    split_results = file("${projectDir}/scripts/split_plink.py")
    version_script = file("${projectDir}/scripts/version_merge.py")
    subpop_tsv = params.settings.subpopulations? channel.fromPath(params.subpopulation.subpop_file):channel.of([])

    // Initialize empty channels for outputs
    ROHOM_out = channel.empty()
    ROHET_out = channel.empty()
    ROHOM_indiv = channel.empty()
    ROHET_indiv = channel.empty()
    ROHOM_islands = channel.empty()
    ROHET_islands = channel.empty()
    all_versions_ch = channel.empty()

    // Read in parameters from the config file
    dr_params = channel.value(params.detectruns.parameters)
    plink_params = channel.value(params.plink.parameters)
    island_info = channel.value(
        tuple (
            params.island.type,
            params.island.minSNP,
            params.island.indiv_count
        )
    )


    // Read in files needed to generate final report
    quarto_files = channel.of(
    tuple(
        file("${projectDir}/scripts/Quarto/demo.qmd"),
        file("${projectDir}/scripts/Quarto/test.R"),
        file("${projectDir}/scripts/Quarto/subpopulations.R"),
        file("${projectDir}/scripts/Quarto/style.scss")
        )
    )

    // Read in gene file if user wants to include gene information in the final report
    if (params.settings.genes) {
        if (params.genes.genefile == "gtf") {
            def genefile = file("${params.genes.gtf}")
            convertgtf = gtftobed(genefile)
            bed_file = convertgtf.bedfile
        } else {
            bed_file = channel.value(file(params.genes.bed))
        }        
    } else {
        bed_file = channel.value([])
    }

    // Read in input files and determine if the user is starting with a VCF or PLINK file. Extract VCF info if starting with a VCF file
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

    //Perform homozygosity analysis if user has selected this option in the config file
    if (params.settings.homozygosity) {
        homstate = "HOM"

        plinkhom_map = plinkhom(infile, plink_params, "")

        split_plink_map = split_plink(split_results, plinkhom_map.roh_file, chromosomes)

        ROHOM_map = ROHOM(split_plink_map.separated_plink, plinkhom_map.bim_file, bed_file, homstate, params.settings.genes)
        ROHOM_out = ROHOM_map.ROH_results
        ROHOM_version = ROHOM_map.all_versions
        ROHOM_indiv = ROHOM_map.indiv_info
        ROHOM_islands = ROHOM_map.finalislands
    }

    //Perform homozygosity analysis if user has selected this option in the config file
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

        ROHET_map = ROHET(all_runs, make_bim.bim_file, bed_file, hetstate, params.settings.genes)
        ROHET_out = ROHET_map.ROH_results
        ROHET_version = ROHET_map.all_versions
        ROHET_indiv = ROHET_map.indiv_info
        ROHET_islands = ROHET_map.finalislands
    }

    //Collect all output information from subworkflows
    ROH_all = ROHOM_out.mix(ROHET_out).collect()
    ROH_island = ROHOM_islands.mix(ROHET_islands).collect()
    indiv_all = ROHOM_indiv.mix(ROHET_indiv).collect()

    //Use all information to generate final report using QUARTO
    quarto_map = quarto(ROH_all,
    ROH_island, 
    indiv_all, 
    vcf_info.vcf_info, 
    quarto_files, 
    params.settings.homozygosity, 
    params.settings.heterozygosity, 
    params.settings.subpopulations, 
    params.inputfiles.start_vcf,
    params.settings.genes,
    subpop_tsv, 
    bed_file,
    island_info )

    //populate version information channel
    all_versions_ch = quarto_map.quarto_version

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

    //output version information to a text file
    all_versions(version_script, all_versions_ch, pipeline_version)
}
