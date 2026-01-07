process overlap {

	label 'overlap'

        publishDir "${params.main.outdir}/output/overlap"

        input:
        path(overlap_script)
        path(tsv_file)
        path(bed_file)

        output:
        path("OL*"), emit: overlap_bed

        script:
        """
        python3 ${overlap_script} --tsv ${tsv_file} --bed ${bed_file}
        
        """
}

process combine_overlap {
        
        label 'combineoverlap'

        publishDir "${params.main.outdir}/output/overlap"

        input:
        path(all_tsvs)
        path(chroms)

        output:
        path("merged.tsv"), emit: ROH_file

        script:
        """
        awk -F'\\t' '
        NR==FNR {
            order[++n] = \$1
            next
        }
        FNR==1 {
            header = \$0
            next
        }
        {
            rows[\$4] = rows[\$4] \$0 ORS
        }
        END {
            print header
            for (i=1; i<=n; i++) {
                if (order[i] in rows)
                    printf "%s", rows[order[i]]
            }
        }
    ' ${chroms} ${all_tsvs} > merged.tsv

        """
}
