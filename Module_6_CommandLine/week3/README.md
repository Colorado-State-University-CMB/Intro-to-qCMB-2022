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
