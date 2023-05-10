
#!/bin/bash

Genome_index=$PATH_Genome_index
fastq_path=$PATH_fastq 
Gene_annotation=$PATH_Genomeannotation_file
NAME="3Dvs2D"
N_R1 ="${PATH_fastq}/${NAME}*R1*.fastq"
N_R2 ="${PATH_fastq}/${NAME}*R2*.fastq"

STAR --genomeDir $Genome_index --readFilesIn $N_R1 $N_R2 --runThreadN 10 --outFileNamePrefix ${NAME}.sam


featureCounts -M -p -T 10 -a $PATH_Genomeannotation_file
-o count_matrix.txt  ${NAME}.sam > ${NAME}.featureCounts.interface
