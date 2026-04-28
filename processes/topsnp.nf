process topsnp {
        ext version: '1.0.1'

	label 'topsnp'

        input:
        val(state)
        path(topsnp_script)
        path(countfile)

        output:
        path("${state}_topsnps.tsv"), emit: topsnps
        path("${task.process}_versions.yml"), emit: topsnps_version

        script:

        """
        python3 ${topsnp_script} \
        --counts ${countfile} \
        --version ${task.ext.version}

        mv versions.yml ${task.process}_versions.yml
        mv topsnps.tsv ${state}_topsnps.tsv
        
        """
}
