process runfastqc {
        
        label 'fastqc'
       
        publishDir "${params.outdir}/output/fastqc", mode: 'move'
        container "${params.apptainer}/fastqc.sif"

        input:
        tuple val(sample_id), path(File1), path(File2), val(adapter1), val(adapter2), val(LibName), val(Barcode), val(Platform)

        output:
        path("*.{html,zip}")
	val("process_complete"), emit: control_1

        script:

        """
        fastqc -o . \
        ${File1} \
        ${File2}
        
        """

}
