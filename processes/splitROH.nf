process splitROH {

	label 'splitROH'

        input:
        path split_script
        path merged_tsv

        output:
        path("*.tsv"), emit: split_tsv

        script:
        """
        Rscript ${split_script} -T ${merged_tsv}
        
        """
}
