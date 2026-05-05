process gtftobed {

        ext version: '1.0.9'
	label 'bedops'

        input:
        path gtffile

        output:
        path("genes.bed"), emit: bedfile

        script:
        """
        convert2bed -i gtf \
        < ${gtffile} \
        > temp.bed

        awk '\$8 ~ /gene/' temp.bed > genes.bed
        
        """
}
