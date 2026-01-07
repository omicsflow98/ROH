process multiqc {

        label 'multiqc'

        publishDir "${params.main.outdir}/output/multiqc"

        input:
        path(ROHall)
        file(mqc_config)

        output:
	path("*.html")
        path("*.zip")

        script:
               
        """
        multiqc --config ${mqc_config} .

	zip -r test_multiqc_report_data.zip test_multiqc_report_data

        """
}
