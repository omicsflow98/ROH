process vcftools {
    ext version: '1.0.4'
	label 'vcftools'

    input:
    path(vcf_file)

    output:
    path("vcf_info.txt"), emit: vcf_info
    path("${task.process}_versions.yml"), emit: vcft_version

    script:

    def cleanname = task.process.split(':')[-1]

    """
    if [[ "${vcf_file}" == *.gz ]]; then
        vcftools --gzvcf ${vcf_file} --site-mean-depth
        vcftools --gzvcf ${vcf_file} --site-quality
        nsnp=\$(zcat ${vcf_file} | grep -v "^#" | wc -l)
    else 
        vcftools --vcf ${vcf_file} --site-mean-depth
        vcftools --vcf ${vcf_file} --site-quality
        nsnp=\$(grep -v "^#" ${vcf_file} | wc -l
)
    fi

    vcfname=\$(basename ${vcf_file})
    avgqual=\$(awk -F'\t' '{sum+=\$3; count++} END {if(count>0) print sum/count}' out.lqual)
    avgdepth=\$(awk -F'\t' '{sum+=\$3; count++} END {if(count>0) print sum/count}' out.ldepth.mean)

    printf "%s\t%s\n" "Type" "Value" > vcf_info.txt
    printf "%s\t%s\n" "VCF file name" "\$vcfname" >> vcf_info.txt
    printf "%s\t%s\n" "Number of SNP" "\$nsnp" >> vcf_info.txt
    printf "%s\t%s\n" "Average SNP Quality" "\$avgqual" >> vcf_info.txt
    printf "%s\t%s" "Average Read depth" "\$avgdepth" >> vcf_info.txt

    cat <<-EOF > ${task.process}_versions.yml
        Process version:
          ${cleanname}: ${task.ext.version}

        Tool version:
          vcftools: \$(vcftools --version | cut -f2 -d' ' | tr -d "()")
    """
}
