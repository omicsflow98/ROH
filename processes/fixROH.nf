process fixROH {

        ext version: '1.0.7'
	label 'fixROH'

        input:
        path fix_script
        path tsv_file

        output:
        path("*_fixed.tsv"), emit: fixed_tsv
        path("${task.process}_versions.yml"), emit: fixROH_version

        script:

        def cleanname = task.process.split(':')[-1]
        
        """
        Rscript ${fix_script} -T ${tsv_file} -V ${task.ext.version}

        mv versions.yml ${task.process}_versions.yml
        
        """
}
