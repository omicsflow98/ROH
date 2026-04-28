process snpislands {

        ext version: '1.0.7'
	label 'snpislands'

        input:
        path ROHfile
        path bim_file

        output:
        path("count.txt"), emit: snpcount
        path("${task.process}_versions.yml"), emit: snpislands_version

        script:

        def cleanname = task.process.split(':')[-1]
        
        """
        mkdir -p tmp
        export TMPDIR="\$PWD/tmp"

        awk 'NR>1 {print \$4"\t"\$7"\t"\$8"\t"\$1}' ${ROHfile} > roh_regions.bed

        less ${bim_file} | \
        awk 'BEGIN{i=1} {print \$1"\t"\$4"\t"\$4"\tSNP"i; i++}' > snpisland.bed

        bedtools intersect \
        -a snpisland.bed \
        -b roh_regions.bed \
        -wa | sort | uniq -c | awk '{print \$5"\t"\$2"\t"\$3"\t"\$4"\t"\$1}'| sort -k 5nr > count.txt


        cat <<-EOF > ${task.process}_versions.yml
                Process version:
                  ${cleanname}: ${task.ext.version}

                Tool version:
                  bcftools: \$(bcftools --version | head -n1 | cut -f2 -d' ')
                  bedtools: \$(bedtools --version | cut -f2 -d' ')
        
        """
}
