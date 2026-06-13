process quarto {
        ext version: '1.0.3'

        label 'quarto'

        publishDir "${params.main.outdir}/output/quarto"

        input:
        path(ROHall)
        path(ROHisland)
        path(indiv_all)
        path(vcf_info)
        tuple path(quartofile), path(plotsfile), path(subpopfile), path(stylesfile)
        val(homstate)
        val(hetstate)
        val(substate)
        val(vcfstate)
        path(subpop)
        path(bedfile)
        tuple val(islandtype), val(minsnp), val(indiv_count)

        output:
        path("*.zip")
        path("*combined.tsv")
        path("gene_*")
        path("${task.process}_versions.yml"), emit: quarto_version

        script:

        def cleanname = task.process.split(':')[-1]
               
        """
        export HOME=\${PWD}

        mv ${bedfile} tempfilename
        mv tempfilename genes.bed
        ${ subpop ? "mv ${subpop} subfile.tsv" : "" }

        quarto render ${quartofile} \\
        -P hom:${homstate} \\
        -P het:${hetstate} \\
        -P substate:${substate} \\
        -P vcf:${vcfstate} \\
        -P type:${islandtype} \\
        -P minsnp:${minsnp} \\
        -P indiv:${indiv_count}

        zip -r demo_files.zip demo.html demo_files

        cat <<-EOF > ${task.process}_versions.yml
        Process version:
          ${cleanname}: ${task.ext.version}

        Tool version:
          Quarto: \$(quarto --version)
	EOF

        """
}
