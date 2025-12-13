process plink {

	label 'plink'

        publishDir "${params.outdir}/output/plink"
	container "${params.apptainer}/plink.sif"

        input:
        path(vcf_file)

        output:
        path("plink.*"), emit: plink_files
        path("plink.hom"), emit: roh_file

        script:
        """
        plink \
        --vcf ${vcf_file} \
        --homozyg \
        --homozyg-snp 100 \
        --homozyg-kb 100 \
        --homozyg-density 10 \
        --homozyg-gap 1000 \
        --homozyg-window-snp 100 \
        --homozyg-window-het 1 \
        --homozyg-window-missing 5 \
        --homozyg-window-threshold 0.05 \
        --homozyg-het 3 \
        --allow-extra-chr \
        --threads 7 \
        --out plink
        
        """
}
