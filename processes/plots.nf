process plot_results {

	label 'plots'

        publishDir "${params.main.outdir}/output/final_plots"

        input:
        path plot_script
        path tsv_tuple
        val state

        output:
        path("*.tsv"), emit: tsv_files
        path("*.png"), emit: plots

        script:
        """
        awk 'FNR==1 && NR!=1 {next} {print}' ${tsv_tuple} > temp.tsv

        Rscript ${plot_script} -T temp.tsv
        
        rm temp.tsv *_fixed.tsv

        for f in * ; do mv "\$f" ${state}_"\$f" ; done
        
        """
}
