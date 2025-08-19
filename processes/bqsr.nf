process BQSR {

	label 'bqsr'

        publishDir "${params.outdir}/output/bqsr/"
	container "${params.apptainer}/gatk.sif"

        input:
        path reference_dir
        path temp_dir
        tuple val(name), path(bam)

        output:
        tuple val(name), path("*.recal.bam"), emit: bam_bqsr

        script:

        """
	gatk \
	--java-options "-Djava.io.tmpdir=${temp_dir}" \
	BaseRecalibratorSpark \
        --input ${bam} \
        --reference ${reference_dir}/genome.fa \
        --known-sites ${reference_dir}/genome.vcf \
	--spark-master local[*] \
        --output ${name}.table

        gatk ApplyBQSRSpark \
        --input ${bam} \
        --reference ${reference_dir}/genome.fa \
        --bqsr-recal-file ${name}.table \
        --tmp-dir ${temp_dir} \
	--spark-master local[*] \
        --output ${name}.recal.bam

        """
}
