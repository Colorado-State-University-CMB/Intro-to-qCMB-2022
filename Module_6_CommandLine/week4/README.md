[See the alignment/variant calling lesson here.](week4/alignment.md)

## Putting the steps together into a script

Using jupyterhub.rc.colorado.edu, 
1. Use the launcher to open a terminal
2. `cd` to your `/projects/$USER/dc_workshop` directory.
3. Sync your file browser to the current directory
   1. Do `pwd` and copy the full path with your mouse.
   2. Paste the copied text into _File -> Open From Path_
3. Use the launcher (+) icon to create a text file. Rename it `pipeline.bash`.
4. Paste in the following text:

```bash
#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=0:15:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --job-name=pipeline

bwa mem -t $SLURM_NTASKS data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq_small/SRR2584866_1.trim.sub.fastq data/trimmed_fastq_small/SRR2584866_2.trim.sub.fastq > results/sam/SRR2584866.aligned.sam

samtools view --threads $SLURM_NTASKS -S -b results/sam/SRR2584866.aligned.sam > results/bam/SRR2584866.aligned.bam

samtools sort --threads $SLURM_NTASKS -o results/bam/SRR2584866.aligned.sorted.bam results/bam/SRR2584866.aligned.bam 

bcftools mpileup --threads $SLURM_NTASKS -O b -o results/bcf/SRR2584866_raw.bcf \
-f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.aligned.sorted.bam

bcftools call --threads $SLURM_NTASKS --ploidy 1 -m -v -o results/vcf/SRR2584866_variants.vcf results/bcf/SRR2584866_raw.bcf 

vcfutils.pl varFilter results/vcf/SRR2584866_variants.vcf > results/vcf/SRR2584866_final_variants.vcf
```
