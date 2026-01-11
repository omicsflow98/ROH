process splitROH {
        ext version: '1.0.6'
	label 'splitROH'

        input:
        path split_script
        path merged_tsv

        output:
        path("*.tsv"), emit: split_tsv
        path("${task.process}_versions.yml"), emit: splitROH_version

        script:

        """
        Rscript ${split_script} \
        -T ${merged_tsv} \
        -V ${task.ext.version}

        mv versions.yml ${task.process}_versions.yml
        
        """
}
