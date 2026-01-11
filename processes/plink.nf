process plink {
        ext version: '1.0.0'

	label 'plink'

        publishDir "${params.main.outdir}/output/plink"

        input:
        path(vcf_file)
        val plink_map
        val chrom

        output:
        path("plink.*"), emit: plink_files
        path("plink.hom"), emit: roh_file
        val(chrom), emit: chromosome
        tuple path("plink.map"), path("plink.ped"), emit: dr_files
        path("${task.process}_versions.yml"), emit: plink_version

        script:

        def cleanname = task.process.split(':')[-1]

        """
        plink \
        --vcf ${vcf_file} \
        --homozyg \
        ${ plink_map.snp ? "--homozyg-snp ${plink_map.snp}" : "" } \
        ${ plink_map.kb ? "--homozyg-kb ${plink_map.kb}" : "" } \
        ${ plink_map.density ? "--homozyg-density ${plink_map.density}" : "" } \
        ${ plink_map.gap ? "--homozyg-gap ${plink_map.gap}" : "" } \
        ${ plink_map.windowsnp ? "--homozyg-window-snp ${plink_map.windowsnp}" : "" } \
        ${ plink_map.WindowHet ? "--homozyg-window-het ${plink_map.WindowHet}" : "" } \
        ${ plink_map.WindowMissing ? "--homozyg-window-missing ${plink_map.WindowMissing}" : "" } \
        ${ plink_map.WindowThreshold ? "--homozyg-window-threshold ${plink_map.WindowThreshold}" : "" } \
        ${ plink_map.HomozygHet ? "--homozyg-het ${plink_map.HomozygHet}" : "" } \
        --allow-extra-chr \
        ${ chrom ? "--chr ${chrom}" : "" } \
        --recode \
        --threads 2 \
        --out plink

        cat <<-EOF > ${task.process}_versions.yml
        Process version:
          ${cleanname}: ${task.ext.version}

        Tool version:
          plink: \$(plink --version | awk '{print \$2}')
        EOF
        
        """
}
