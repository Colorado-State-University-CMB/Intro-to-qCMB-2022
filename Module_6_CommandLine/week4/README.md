[See Monday's alignment/variant-calling lesson here.](alignment.md)

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


### Check output, try other parameters

- Where is the output file? What's in it?
- How long did it take?
- Replace ACCESSION with SRR2584863 and rerun
- Increase `--ntasks` to 5 and rerun
- Did it go faster?

## Modify to run on all samples 

The following modification to the script makes use of `basename`, which was used in the trimmomatic script. 

We will use it to extract the accession number from the files matched by `data/trimmed_fastq_small/*_1.trim.sub.fastq`

### Example basename usage

```bash
$ ls data/trimmed_fastq_small/*_1.trim.sub.fastq
data/trimmed_fastq_small/SRR2584863_1.trim.sub.fastq  data/trimmed_fastq_small/SRR2584866_1.trim.sub.fastq  data/trimmed_fastq_small/SRR2589044_1.trim.sub.fastq
$ basename data/trimmed_fastq_small/SRR2584863_1.trim.sub.fastq 
SRR2584863_1.trim.sub.fastq
```

Adding a second argument, a file extension, will clip that extension from the filename. For our example, the returned string is the accession number of the sample.

```bash
$ basename data/trimmed_fastq_small/SRR2584863_1.trim.sub.fastq  _1.trim.sub.fastq
SRR2584863
```

To save the value to a variable, wrap the entire command in `$(  )`

```bash
acc=$(basename data/trimmed_fastq_small/SRR2584863_1.trim.sub.fastq  _1.trim.sub.fastq)
echo $acc
SRR2584863
```

### Script with loop

Using the for-loop introduced in `trimmomatic.bash`, we'll run the `basename` command on each file in turn.

```bash
#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --time=0:15:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --job-name=pipeline-step2

module purge
source /curc/sw/anaconda3/latest
conda activate variant-calling

for fname in data/trimmed_fastq_small/*_1.trim.sub.fastq
do
    # It is convention to indent the loop contents to make it visually stand out

    # this part extracts the accession number from the string stored in $fname
    ACCESSION=$(basename $fname _1.trim.sub.fastq)
    echo "Working on ACCESSION: $ACCESSION"

    # The rest of the body is the same, but notice the done statement at the end!

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

done  ### <-- DON'T FORGET TO ADD THIS 'done' statement !
```

Submit the job again using `sbatch`.

## Array jobs: Modify to run on all samples in parallel

Sometimes the pipeline has too much data to process to complete all samples in the available time (24 hours for normal qos, shas). A solution is to submit a different job individually for each sample. `Array jobs` provide a built-in way to do this in a single command.

We will take advantage of the "array" feature of sbatch to submit 3 instances of the script (3 separate jobs), automatically.

Syntax:

`sbatch --array=1-3 pipeline.bash`

The above example submits pipeline.bash 3 times, but setting the variable `SLURM_ARRAY_TASK_ID` to 1, 2, and 3 for each respective instance.

We can change our script to take advantage of that variable by counting our passes through the loop, and only executing the contents when the loop iteration matches `SLURM_ARRAY_TASK_ID`

```bash
i=0 # ADD A COUNTER
for fname in data/trimmed_fastq_small/*_1.trim.sub.fastq
do
    i=$((i+1)) # INCREMENT COUNTER

    # SKIP INPUT FILE WHEN THE ARRAY INDEX DOESN'T MATCH
    if [ $i -ne $SLURM_ARRAY_TASK_ID ]
    then
        continue
    fi  
    
    ACCESSION=$(basename $fname _1.trim.sub.fastq)
    echo "Working on ACCESSION: $ACCESSION"

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

done  ### <-- DON'T FORGET TO ADD THIS 'done' statement !
```

Try it!

`sbatch --array=1-3 pipeline.bash`

What do you see with `squeue -u $USER` ?

## Script on full data

There are two major changes in this script:

1. Changed all "sub" input files, which are in decompressed fastq, to the original, trimmed samples, which are gzipped (fastq.gz).
2. Replaced the loop by setting `data/trimmed_fastq/*_1.trim.fastq.gz` to an array, and then select the array entry using $SLURM_ARRAY_TASK_ID.

This array script is run with `sbatch --array=0-2 pipeline.bash`, because bash arrays are 0-based (R is 1-based).

```bash
#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --time=0:20:00
#SBATCH --qos=normal
#SBATCH --partition=shas
#SBATCH --job-name=pipeline
args=( data/trimmed_fastq/*_1.trim.fastq.gz )
fname=${args[$SLURM_ARRAY_TASK_ID]}

ACCESSION=$(basename $fname _1.trim.fastq.gz)

echo "Working on ACCESSION: $ACCESSION"

# bwa index was run to create the alignment index in data/ref_genome/

bwa mem -t $SLURM_NTASKS data/ref_genome/ecoli_rel606.fasta \
data/trimmed_fastq/${ACCESSION}_1.trim.fastq.gz \
data/trimmed_fastq/${ACCESSION}_2.trim.fastq.gz > results/sam/${ACCESSION}.aligned.sam

samtools view --threads $SLURM_NTASKS -S -b results/sam/${ACCESSION}.aligned.sam > results/bam/${ACCESSION}.aligned.bam
samtools sort --threads $SLURM_NTASKS -o results/bam/${ACCESSION}.aligned.sorted.bam results/bam/${ACCESSION}.aligned.bam 

bcftools mpileup --threads $SLURM_NTASKS -O b -o results/bcf/${ACCESSION}_raw.bcf \
-f data/ref_genome/ecoli_rel606.fasta results/bam/${ACCESSION}.aligned.sorted.bam

bcftools call --threads $SLURM_NTASKS --ploidy 1 -m -v -o results/vcf/${ACCESSION}_variants.vcf results/bcf/${ACCESSION}_raw.bcf 

vcfutils.pl varFilter results/vcf/${ACCESSION}_variants.vcf > results/vcf/${ACCESSION}_final_variants.vcf
```

Run with `sbatch --array=0-2 pipeline.bash`

The full data takes much longer, even with 6 CPUs (ntasks):

```
       JobID    JobName  AllocCPUS      State ExitCode    Elapsed  Timelimit              Submit               Start                 End 
------------ ---------- ---------- ---------- -------- ---------- ---------- ------------------- ------------------- -------------------
9801084_0    pipeline-+          6  COMPLETED      0:0   00:09:37   00:20:00 2022-03-29T17:15:52 2022-03-29T17:15:52 2022-03-29T17:25:29 
9801084_1    pipeline-+          6  COMPLETED      0:0   00:15:58   00:20:00 2022-03-29T17:15:52 2022-03-29T17:15:52 2022-03-29T17:31:50 
9801084_2    pipeline-+          6  COMPLETED      0:0   00:06:47   00:20:00 2022-03-29T17:15:52 2022-03-29T17:25:41 2022-03-29T17:32:28 
```

## Next steps

The fastqc and trimmomatic steps can be added to the beginning of this script, creating a full pipeline that only requires a few steps to set up from scratch:

1. Download data and organize directory structure
2. Create bwa index
3. Run pipeline script

What would you have to do to integrate all of the previous commands of fastqc and trimmomatic?
