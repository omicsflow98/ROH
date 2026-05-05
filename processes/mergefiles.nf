process merge_files {

        ext version: '1.0.9'
	label 'plots'

        input:
        path tsv_tuple
        val state

        output:
        path("${state}_temp.tsv"), emit: tsv_file

        script:
        """
        awk 'FNR==1 && NR!=1 {next} {print}' ${tsv_tuple} > ${state}_temp.tsv
        
        """
}
