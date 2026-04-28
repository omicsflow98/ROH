process splitROH {
        ext version: '1.0.6'
	label 'splitROH'

        input:
        path split_script
        path merged_tsv
        val state
        val genome

        output:
        path("*.tsv"), emit: split_tsv
        path("${task.process}_versions.yml"), emit: splitROH_version
        path("*.indiv"), emit: indiv_file

        script:

        """
        Rscript ${split_script} \
        -T ${merged_tsv} \
        -V ${task.ext.version} \
        -L ${genome} \
        -S ${state}

        mv versions.yml ${task.process}_versions.yml
        mv indiv.tsv ${state}.indiv
        
        """
}
