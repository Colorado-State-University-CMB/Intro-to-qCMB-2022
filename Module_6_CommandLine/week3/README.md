### From last time

From your home directory (use `cd` without arguments to get there).

```
$ ls dc_workshop/results/fastqc_untrimmed_reads
SRR2584863_1_fastqc.html SRR2584863_2_fastqc.html SRR2584866_1_fastqc.html SRR2584866_2_fastqc.html SRR2589044_1_fastqc.html SRR2589044_2_fastqc.html
SRR2584863_1_fastqc.zip  SRR2584863_2_fastqc.zip  SRR2584866_1_fastqc.zip  SRR2584866_2_fastqc.zip  SRR2589044_1_fastqc.zip  SRR2589044_2_fastqc.zip
```

Mine also expanded directories in `dc_workshop/data/untrimmed_fastq` (they end in \_fastqc).

```
$ ls dc_workshop/results/fastqc_untrimmed_reads
SRR2584863_2_fastqc          SRR2584866_2_fastqc          SRR2584863_1.fastq           
SRR2584866_1.fastq           SRR2589044_1.fastq.gz        SRR2589044_2.fastq.gz
SRR2584863_1_fastqc          SRR2584866_1_fastqc          SRR2584863_2.fastq.gz        
SRR2584866_2.fastq.gz        SRR2589044_1_fastqc          SRR2589044_2_fastqc
```

Let's move those to the results directory:

```
$ mv dc_workshop/data/untrimmed_fastq/*_fastqc dc_workshop/results/fastqc_untrimmed_reads
$ ls dc_workshop/results/fastqc_untrimmed_reads                                         
SRR2584863_1_fastqc      SRR2584863_2_fastqc      SRR2584866_1_fastqc      SRR2584866_2_fastqc      SRR2589044_1_fastqc      SRR2589044_2_fastqc
SRR2584863_1_fastqc.html SRR2584863_2_fastqc.html SRR2584866_1_fastqc.html SRR2584866_2_fastqc.html SRR2589044_1_fastqc.html SRR2589044_2_fastqc.html
SRR2584863_1_fastqc.zip  SRR2584863_2_fastqc.zip  SRR2584866_1_fastqc.zip  SRR2584866_2_fastqc.zip  SRR2589044_1_fastqc.zip  SRR2589044_2_fastqc.zip
```


```
$ cd dc_workshop/data/untrimmed_fastq 
$ ls
SRR2584863_1.fastq    SRR2584863_2.fastq.gz SRR2584866_1.fastq    
SRR2584866_2.fastq.gz SRR2589044_1.fastq.gz SRR2589044_2.fastq.gz
```

It should look something like the above, although I unzipped more than one fastq file.


### pwd check

Check to make sure you are in `dc_workshop/data/untrimmed_fastq` before proceeding.

```
$ pwd
/Users/david/dc_workshop/data/untrimmed_fastq
```

My home directory `/Users/david/` will differ from yours.

# Trimming and Filtering

## Cleaning reads
In the previous episode, we took a high-level look at the quality of each of our samples using FastQC. We visualized per-base quality graphs showing the distribution of read quality at each base across all reads in a sample and extracted information about which samples fail which quality checks. Some of our samples failed quite a few quality metrics used by FastQC. This does not mean, though, that our samples should be thrown out! It is very common to have some quality metrics fail, and this may or may not be a problem for your downstream application. For our variant calling workflow, we will be removing some of the low quality sequences to reduce our false positive rate due to sequencing error.

We will use a program called [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) to filter poor quality reads and trim poor quality bases from our samples.

## Trimmomatic options
Trimmomatic has a variety of options to trim your reads. If we run the following command, we can see some of our options.

```bash
$ trimmomatic
```

```bash
Usage: 
       PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   or: 
       SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...
   or: 
       -version

```
This output shows us that we must first specify whether we have paired end (`PE`) or single end (`SE`) reads. Next, we specify what flag we would like to run. For example, you can specify threads to indicate the number of processors on your computer that you want Trimmomatic to use. In most cases using multiple `threads` (processors) can help to run the trimming faster. These flags are not necessary, but they can give you more control over the command. The flags are followed by positional arguments, meaning the order in which you specify them is important. In paired end mode, Trimmomatic expects the two input files, and then the names of the output files. These files are described below. While, in single end mode, Trimmomatic will expect 1 file as input, after which you can enter the optional settings and lastly the name of the output file.

| option |	meaning |
| --- | --- |
| `<inputFile1>`	| Input reads to be trimmed. Typically the file name will contain an `_1` or `_R1` in the name. |
| `<inputFile2>`	| Input reads to be trimmed. Typically the file name will contain an `_2` or `_R2` in the name. |
| `<outputFile1P>`	| Output file that contains surviving pairs from the `_1` file.  |
| `<outputFile1U>`	| Output file that contains orphaned reads from the `_1` file.  |
| `<outputFile2P>`	| Output file that contains surviving pairs from the `_2` file. |
| `<outputFile2U>`	| Output file that contains orphaned reads from the `_2` file.

The last thing trimmomatic expects to see is the trimming parameters:

|step	| meaning |
| --- | --- |
| `ILLUMINACLIP`	| Perform adapter removal. | 
| `SLIDINGWINDOW` | 	Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. |
| `LEADING`	| Cut bases off the start of a read, if below a threshold quality. |
| `TRAILING`	| Cut bases off the end of a read, if below a threshold quality. |
| `CROP`	| Cut the read to a specified length. |
| `HEADCROP`	| Cut the specified number of bases from the start of the read. |
| `MINLEN`	| Drop an entire read if it is below a specified length. |
| `TOPHRED33`	| Convert quality scores to Phred-33. |
| `TOPHRED64`	|Convert quality scores to Phred-64. |

We will use only a few of these options and trimming steps in our analysis. It is important to understand the steps you are using to clean your data. For more information about the Trimmomatic arguments and options, see the [Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf).

However, a complete command for Trimmomatic will look something like the command below. **This command is an example and will not work, as we do not have the files it refers to:**

```bash
$ trimmomatic PE -threads 4 SRR_1056_1.fastq SRR_1056_2.fastq  \
              SRR_1056_1.trimmed.fastq SRR_1056_1un.trimmed.fastq \
              SRR_1056_2.trimmed.fastq SRR_1056_2un.trimmed.fastq \
              ILLUMINACLIP:SRR_adapters.fa SLIDINGWINDOW:4:20
```

In this example, we have told Trimmomatic:

| code	| meaning |
| ---  | --- |
| `PE`	| that it will be taking a paired end file as input |
| `-threads 4`	| to use four computing threads to run (this will speed up our run) |
| `SRR_1056_1.fastq`	| the first input file name |
| `SRR_1056_2.fastq`	| the second input file name |
| `SRR_1056_1.trimmed.fastq` |	the output file for surviving pairs from the `_1` file |
| `SRR_1056_1un.trimmed.fastq`	| the output file for orphaned reads from the `_1` file |
| `SRR_1056_2.trimmed.fastq`	| the output file for surviving pairs from the `_2` file |
| `SRR_1056_2un.trimmed.fastq` |	the output file for orphaned reads from the `_2` file |
| `ILLUMINACLIP:SRR_adapters.fa` | to clip the Illumina adapters from the input file using the adapter sequences listed in `SRR_adapters.fa` |
| `SLIDINGWINDOW:4:20`	 | to use a sliding window of size 4 that will remove bases if their phred score is below 20 |

---

**Some of the commands we ran in this lesson are long!** When typing a long command into your terminal, you can use the `\` character to separate code chunks onto separate lines. This can make your code more readable.

---

## Running Trimmomatic

Now we will run Trimmomatic on our data. To begin, navigate to your untrimmed_fastq data directory:

```bash
$ cd ~/dc_workshop/data/untrimmed_fastq
```

We are going to run Trimmomatic on one of our paired-end samples. While using FastQC we saw that Nextera adapters were present in our samples. The adapter sequences came with the installation of trimmomatic, so we will first copy these sequences into our current directory.

```bash
$ cp ~/.miniconda3/pkgs/trimmomatic-0.38-0/share/trimmomatic-0.38-0/adapters/NexteraPE-PE.fa .
```

### Locating the installed adapter sequences

!!! We have to search for our installation. !!!

**Do not use cd in the following commands**, stay in the `~/dc_workshop/data/untrimmed_fastq directory`.

#### 1. find your conda install location

```bash
$ conda env list
# conda environments:
#
base                     /Users/david/opt/miniconda3
derptools                /Users/david/opt/miniconda3/envs/derptools
elt-2-rev                /Users/david/opt/miniconda3/envs/elt-2-rev
variant-calling       *  /Users/david/opt/miniconda3/envs/variant-calling

```

#### 2. Add '/share' to that directory to find your version of trimmomatic

```
$ ls /Users/david/opt/miniconda3/envs/variant-calling/share
aclocal            et                 fontconfig         info               man                tabset             trimmomatic        xml
doc                examples           gettext            locale             nghttp2            terminfo           trimmomatic-0.39-2 zoneinfo
```

#### 3. Look in the directory with the version attached to find the adapter files

```
ls /Users/david/opt/miniconda3/envs/variant-calling/share/trimmomatic-0.39-2 
LICENSE                   build_env_setup.sh        metadata_conda_debug.yaml trimmomatic.jar
adapters                  conda_build.sh            trimmomatic
```

There is an `adapters` directory.

```
$ ls /Users/david/opt/miniconda3/envs/variant-calling/share/trimmomatic-0.39-2/adapters
NexteraPE-PE.fa TruSeq2-PE.fa   TruSeq2-SE.fa   TruSeq3-PE-2.fa TruSeq3-PE.fa   TruSeq3-SE.fa
```

#### 4. Copy the adapters to your current directory (should be ~/dc_workshop/untrimmed_fastq)

The adapter sequences are in `NexteraPE-PE.fa`. Copy *yours* to your current directory. 

My copy command is this:

```
$ cp /Users/david/opt/miniconda3/envs/variant-calling/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa .
```

---

### Finishing the command

We will also use a sliding window of size 4 that will remove bases if their phred score is below 20 (like in our example above). We will also discard any reads that do not have at least 25 bases remaining after this trimming step. Three additional pieces of code are also added to the end of the ILLUMINACLIP step. These three additional numbers (2:40:15) tell Trimmimatic how to handle sequence matches to the Nextera adapters. A detailed explanation of how they work is advanced for this particular lesson. For now we will use these numbers as a default and recognize they are needed to for Trimmomatic to run properly. This command will take a few minutes to run.

CODE:

```
$ trimmomatic PE SRR2589044_1.fastq.gz SRR2589044_2.fastq.gz \
                SRR2589044_1.trim.fastq.gz SRR2589044_1un.trim.fastq.gz \
                SRR2589044_2.trim.fastq.gz SRR2589044_2un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
```

My OUTPUT:

```
TrimmomaticPE: Started with arguments:
 SRR2589044_1.fastq.gz SRR2589044_2.fastq.gz SRR2589044_1.trim.fastq.gz SRR2589044_1un.trim.fastq.gz SRR2589044_2.trim.fastq.gz SRR2589044_2un.trim.fastq.gz SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
Using PrefixPair: 'AGATGTGTATAAGAGACAG' and 'AGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'
ILLUMINACLIP: Using 1 prefix pairs, 4 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Quality encoding detected as phred33
Input Read Pairs: 1107090 Both Surviving: 885220 (79.96%) Forward Only Surviving: 216472 (19.55%) Reverse Only Surviving: 2850 (0.26%) Dropped: 2548 (0.23%)
TrimmomaticPE: Completed successfully
```

Mine took a little over a minute and a half.

### Exercise

Use the output from your Trimmomatic command to answer the following questions.

1) What percent of reads did we discard from our sample? 
2) What percent of reads did we keep both pairs?

<!--
1) 0.23% 2) 79.96%
-->

You may have noticed that Trimmomatic automatically detected the quality encoding of our sample. It is always a good idea to double-check this or to enter the quality encoding manually.

We can confirm that we have our output files:

```
$ ls SRR2589044*
SRR2589044_1.fastq.gz        SRR2589044_1un.trim.fastq.gz SRR2589044_2.trim.fastq.gz
SRR2589044_1.trim.fastq.gz   SRR2589044_2.fastq.gz        SRR2589044_2un.trim.fastq.gz
```
The output files are also FASTQ files. It should be smaller than our input file, because we have removed reads. We can confirm this:

```
ls -lh SRR2589044*
ls -lh SRR2589044*
-rw-r--r--  1 david  staff   123M Mar  9 08:51 SRR2589044_1.fastq.gz
-rw-r--r--  1 david  staff    93M Mar 20 15:40 SRR2589044_1.trim.fastq.gz
-rw-r--r--  1 david  staff    17M Mar 20 15:40 SRR2589044_1un.trim.fastq.gz
-rw-r--r--  1 david  staff   128M Mar  9 08:52 SRR2589044_2.fastq.gz
-rw-r--r--  1 david  staff    91M Mar 20 15:40 SRR2589044_2.trim.fastq.gz
-rw-r--r--  1 david  staff   271K Mar 20 15:40 SRR2589044_2un.trim.fastq.gz
```
