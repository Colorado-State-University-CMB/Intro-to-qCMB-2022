[See the alignment/variant-calling lesson here.](week4/alignment.md)

## Putting the steps together into a script

Why repeat the steps when we have already processed everything?

**Reproducible Research**

We have done all of the variant-calling steps individually, but we would like to put all of the steps in a script. 

**Scalability**

We have only processed the subsampled data. The original data files are 6-12 times larger than the "sub" files.  

We will take full advantage of cluster computing power in order to handle the larger volume.

**Generalizability**

Once you have a pipeline, you can reuse it on new data, or alternatively, copy it and modify it for a new application.


## Starting a pipeline script

We will begin by using code that we already know works, and modify it to run on all of the samples.

### If you are unable to launch a jupyterhub session

```
ssh eid@colostate.edu@login.rc.colorado.edu
Password:
```

1. Type "password,push" (the typing will be invisible) and approve the Duo notice.
2. `cd` to your `/projects/$USER/dc_workshop` directory.
3. nano `pipeline.bash`
4. Paste in the script below.

### Using jupyterhub.rc.colorado.edu, 
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
#SBATCH --time=0:05:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --job-name=pipeline

module purge
source /curc/sw/anaconda3/latest
conda activate variant-calling

bwa mem -t $SLURM_NTASKS data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq_small/SRR2584866_1.trim.sub.fastq data/trimmed_fastq_small/SRR2584866_2.trim.sub.fastq > results/sam/SRR2584866.aligned.sam

samtools view --threads $SLURM_NTASKS -S -b results/sam/SRR2584866.aligned.sam > results/bam/SRR2584866.aligned.bam

samtools sort --threads $SLURM_NTASKS -o results/bam/SRR2584866.aligned.sorted.bam results/bam/SRR2584866.aligned.bam 

bcftools mpileup --threads $SLURM_NTASKS -O b -o results/bcf/SRR2584866_raw.bcf \
-f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.aligned.sorted.bam

bcftools call --threads $SLURM_NTASKS --ploidy 1 -m -v -o results/vcf/SRR2584866_variants.vcf results/bcf/SRR2584866_raw.bcf 

vcfutils.pl varFilter results/vcf/SRR2584866_variants.vcf > results/vcf/SRR2584866_final_variants.vcf
```

### Replace the explicit "SRR2584866" to ${ACCESSION}

We will set a variable `ACCESSION` to a prefix in our data. This variable can be changed before running the script to run the entire pipeline on the forward/reverse trimmed reads for a given sample.

Add `ACCESSION=SRR2584866` before the comment about the bwa index (starts with #).

Replace all of the _remaining_ occurances of `SRR2584866` with `${ACCESSION}`. !! Try using CTRL-F/CMD-F to search and replace !!

```bash
#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=0:05:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --job-name=pipeline

module purge
source /curc/sw/anaconda3/latest
conda activate variant-calling

ACCESSION=SRR2584866

# bwa index was run to create the alignment index in data/ref_genome/

bwa mem -t $SLURM_NTASKS data/ref_genome/ecoli_rel606.fasta \
data/trimmed_fastq_small/${ACCESSION}_1.trim.sub.fastq \
data/trimmed_fastq_small/${ACCESSION}_2.trim.sub.fastq > results/sam/${ACCESSION}.aligned.sam

samtools view --threads $SLURM_NTASKS -S -b results/sam/${ACCESSION}.aligned.sam > results/bam/${ACCESSION}.aligned.bam
samtools sort --threads $SLURM_NTASKS -o results/bam/${ACCESSION}.aligned.sorted.bam results/bam/${ACCESSION}.aligned.bam 

bcftools mpileup --threads $SLURM_NTASKS -O b -o results/bcf/${ACCESSION}_raw.bcf \
-f data/ref_genome/ecoli_rel606.fasta results/bam/${ACCESSION}.aligned.sorted.bam

bcftools call --threads $SLURM_NTASKS --ploidy 1 -m -v -o results/vcf/${ACCESSION}_variants.vcf results/bcf/${ACCESSION}_raw.bcf 

vcfutils.pl varFilter results/vcf/${ACCESSION}_variants.vcf > results/vcf/${ACCESSION}_final_variants.vcf
```

Run the script with `sbatch pipeline.bash`

Check the pending/running status with `squeue -u $USER`

Check `sacct -X` for another view.  

#### !TIP! Set an alias for sacct to return more fields

`sacct -X --format JobID,JobName,AllocCPUS,State,ExitCode,Elapsed,TimeLimit,Submit,Start,End`

I set an alias, and then only have to type `sa` to get my version of the command:

```bash
alias sa='sacct -X --format JobID,JobName,AllocCPUS,State,ExitCode,Elapsed,TimeLimit,Submit,Start,End'
sa
```

Example output:

```
       JobID    JobName  AllocCPUS      State ExitCode    Elapsed  Timelimit              Submit               Start                 End 
------------ ---------- ---------- ---------- -------- ---------- ---------- ------------------- ------------------- ------------------- 
9800221      spawner-j+          1 CANCELLED+      0:0   02:09:17   12:00:00 2022-03-29T13:55:18 2022-03-29T13:55:38 2022-03-29T16:04:55 
9800576        pipeline          2     FAILED    127:0   00:00:01   00:15:00 2022-03-29T14:09:03 2022-03-29T14:09:18 2022-03-29T14:09:19 
9800585        pipeline          2  COMPLETED      0:0   00:01:47   00:15:00 2022-03-29T14:10:22 2022-03-29T14:10:49 2022-03-29T14:12:36 
9800918        pipeline          2  COMPLETED      0:0   00:01:40   00:15:00 2022-03-29T15:51:26 2022-03-29T15:51:54 2022-03-29T15:53:34 
9800973        pipeline          5  COMPLETED      0:0   00:04:52   00:15:00 2022-03-29T16:03:34 2022-03-29T16:04:59 2022-03-29T16:09:51
```


### How long did it take?




