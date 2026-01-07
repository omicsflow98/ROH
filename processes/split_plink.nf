process split_plink {

	label 'split_plink'

        publishDir "${params.main.outdir}/output/split_plink"

        input:
        path(split_script)
        path(hom_file)
        path(chrom_file)

        output:
        path("*.tsv"), emit: separated_plink

        script:
        """
        python3 ${split_script} --hom ${hom_file} --chroms ${chrom_file}
        
        """
}
