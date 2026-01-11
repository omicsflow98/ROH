process all_versions {

        label 'versions'

        publishDir "${params.main.outdir}/output", mode: 'move'

        input:
        path(version_script)
        path(all_versions)
        val(pipeline_version)

        output:
        path("version_info.yml")

        script:
               
        """
        python3 ${version_script} --files ${all_versions}

        {
        echo "Nextflow: ${nextflow.version}"
        echo "Pipeline: ${pipeline_version}"
        echo
        cat version_info.yml
        } > tmp && mv tmp version_info.yml

        """
}
