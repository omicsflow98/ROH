process detectruns {
        ext version: '1.0.2'

	label 'detectruns'

        input:
        path (dr_script)
        tuple path(map_file), path(ped_file)
        val dr_map 
        val chrom

        output:
        path("*.tsv"), emit: dr_tsv
        path("${task.process}_versions.yml"), emit: dr_version

        script:

        def cleanname = task.process.split(':')[-1]

        """
        Rscript ${dr_script} \
        --ped ${ped_file} \
        --map ${map_file} \
        ${ dr_map.window ? "--window ${dr_map.window}" : "" } \
        ${ dr_map.threshold ? "--threshold ${dr_map.threshold}" : "" } \
        ${ dr_map.minSNP ? "--minSNP ${dr_map.minSNP}" : "" } \
        ${ dr_map.maxhom ? "--maxhom ${dr_map.maxhom}" : "" } \
        ${ dr_map.maxmissing ? "--maxmissing ${dr_map.maxmissing}" : "" } \
        ${ dr_map.maxgap ? "--maxgap ${dr_map.maxgap}" : "" } \
        ${ dr_map.minlength ? "--minlength ${dr_map.minlength}" : "" } \
        ${ dr_map.mindensity ? "--mindensity ${dr_map.mindensity}" : "" } \
        ${ dr_map.maxhomrun ? "--maxhomrun ${dr_map.maxhomrun}" : "" } \
        ${ chrom ? "--chromosome ${chrom}" : "" } \
        --version ${task.ext.version}

        mv versions.yml ${task.process}_versions.yml 

        """
}
