process fixROH {

	label 'fixROH'

        input:
        path fix_script
        path tsv_file

        output:
        path("*_fixed.tsv"), emit: fixed_tsv

        script:
        """
        Rscript ${fix_script} -T ${tsv_file}
        
        """
}
