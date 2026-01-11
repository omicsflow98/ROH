process split_plink {
        ext version: '1.0.1'

	      label 'split_plink'

        publishDir "${params.main.outdir}/output/split_plink"

        input:
        path(split_script)
        path(hom_file)
        path(chrom_file)

        output:
        path("*.tsv"), emit: separated_plink
        path("${task.process}_versions.yml"), emit: splitplink_version

        script:

        """
        python3 ${split_script} \
        --hom ${hom_file} \
        --chroms ${chrom_file} \
        --version ${task.ext.version}

        mv versions.yml ${task.process}_versions.yml
        
        """
}
