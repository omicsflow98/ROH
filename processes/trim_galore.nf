process trim_galore {
        
        label 'trim_galore'

        publishDir "${params.outdir}/output/trimmed_fastq"
        container "${params.apptainer}/trim_galore.sif"

        input:
        tuple val(sample_id), path(File1), path(File2), val(adapter1), val(adapter2), val(LibName), val(Barcode), val(Platform)

        output:
        tuple val(sample_id), val(LibName), val(Barcode), val(Platform), path("*.gz"), emit: trimmed_fastq
	path("*.txt"), emit: fastqc

        script:

        """
        trim_galore \
        --paired \
        --phred33 \
        --gzip \
	--cores 2 \
        --adapter ${adapter1} \
        --adapter2 ${adapter2} \
        --basename ${sample_id} \
        ${File1} \
        ${File2}

        """ 

}
