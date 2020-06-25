#Script for generating coverage file for a fasta file. Uses bowtie2 nad samtools
NAME="fasta_file"
READS1=/path/to/forward/reads
READS2=/path/to/reverse/reads
CPUS=12

bowtie2-build "$NAME".fa bam/index
bowtie2 -x bam/index -1 "$READS1" -2 "$READS2" -S "$NAME".sam -p "$CPUS" 2>bowtie.out 1>bowtie.err
samtools view -b -S "$NAME".sam > "$NAME".bam
samtools sort "$NAME".bam > "$NAME".bam.sorted.bam
samtools depth "$NAME".bam.sorted.bam > "$NAME".coverage
#see gc-coverage_plot.R for further analysis steps






