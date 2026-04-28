process plot_results {

        ext version: '1.0.9'
	label 'plots'

        publishDir "${params.main.outdir}/output/final_plots"

        input:
        path plot_script
        path tsv_tuple
        val state

        output:
        path("*.tsv"), emit: tsv_files
        path("*.png"), emit: plots
        path("${state}_versions.yml"), emit: plots_version

        script:
        """
        awk 'FNR==1 && NR!=1 {next} {print}' ${tsv_tuple} > temp.tsv

        Rscript ${plot_script} -T temp.tsv -V ${task.ext.version}
        
        rm *_fixed.tsv

        mv temp.tsv temp_file

        for f in * ; do mv "\$f" ${state}_"\$f" ; done
        
        """
}
