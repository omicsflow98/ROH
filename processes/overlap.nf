process overlap {
    ext version: '1.0.4'
	label 'overlap'

    publishDir "${params.main.outdir}/output/overlap"

    input:
    path(overlap_script)
    path(tsv_file)
    path(bed_file)

    output:
    path("OL*"), emit: overlap_bed
    path("${task.process}_versions.yml"), emit: OL_version

    script:

    """
    python3 ${overlap_script} \
    --tsv ${tsv_file} \
    --bed ${bed_file} \
    --version ${task.ext.version}

    mv versions.yml ${task.process}_versions.yml
    """
}

process combine_overlap {
    ext version: '1.0.5'    
    label 'combineoverlap'

    publishDir "${params.main.outdir}/output/overlap"

    input:
    path(all_tsvs)
    path(chroms)

    output:
    path("merged.tsv"), emit: ROH_file
    path("${task.process}_versions.yml"), emit: combineOL_version

    script:

    def cleanname = task.process.split(':')[-1]

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

    cat <<-EOF > ${task.process}_versions.yml
    Process version:
      ${cleanname}: ${task.ext.version}
    EOF

        """
}
