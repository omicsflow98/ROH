process multiqc {
        ext version: '1.0.3'

        label 'multiqc'

        publishDir "${params.main.outdir}/output/multiqc"

        input:
        path(ROHall)
        file(mqc_config)

        output:
	path("*.html")
        path("*.zip")
        path("${task.process}_versions.yml"), emit: multiqc_version

        script:

        def cleanname = task.process.split(':')[-1]
               
        """
        multiqc --config ${mqc_config} .

	zip -r test_multiqc_report_data.zip test_multiqc_report_data

        cat <<-EOF > ${task.process}_versions.yml
        Process version:
          ${cleanname}: ${task.ext.version}

        Tool version:
          multiqc: \$(multiqc --version | awk '{print \$3}')
	EOF

        """
}
