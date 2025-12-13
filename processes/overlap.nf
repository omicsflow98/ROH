process overlap {

	label 'overlap'

        publishDir "${params.outdir}/output/overlap"
	container "${params.apptainer}/python.sif"

        input:
        path(overlap_script)
        path(tsv_file)
        path(bed_file)

        output:
        path("OL*"), emit: overlap_bed

        script:
        """
        python3 ${overlap_script} --tsv ${tsv_file} --bed ${bed_file}
        
        """
}

process combine_overlap {
        
        label 'overlap'

        publishDir "${params.outdir}/output/overlap"

        input:
        path(all_tsvs)

        output:
        path("merged.tsv")

        script:
        """
        awk 'FNR==1 && NR!=1 { next } { print }' *.tsv > merged.tsv

        """
}
