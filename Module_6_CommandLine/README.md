# Module 6 - Command Line, Summit and Variant-Calling Workflow

[Week1 - basics](week1)

[Week2 - Summit/Variant Calling](week2)

[Week3 - Variant Calling Part 2 - Read Trimming](week3)

[Week4 - Variant Calling Part 3 - Alignment](week4)


## Resources and tips

### CURC Documentation

https://curc.readthedocs.io/en/latest/

Conda instructions: https://curc.readthedocs.io/en/latest/software/python.html#configuring-conda-with-condarc


#### Special setup for CSU users 
... involves using a directory that omits the `@` sign. It is not always necessary, but some PERL code fails.

Below is how I would set up my configuration (username dcking@colostate.edu).

(The warning below is because I already have these settings)
```
$ conda activate /curc/sw/anaconda3/latest
$ conda config --add pkgs_dirs /projects/.colostate.edu/dcking/.conda_pkgs
Warning: '/projects/.colostate.edu/dcking/.conda_pkgs' already in 'pkgs_dirs' list, moving to the top
$ conda config --add envs_dirs '/projects/.colostate.edu/dcking/software/anaconda/envs'
Warning: '/projects/.colostate.edu/dcking/software/anaconda/envs' already in 'envs_dirs' list, moving to the top
```
Files will still be installed in `/projects/$USER`, but through a link which doesn't use an `@` sign.

#### An example R environment

\*Make sure bioconda and conda-forge are added to your list of channels.

You can install R packages directly (search for them with `conda search` to get the precise name).

```
$ conda create --name deseq2 r-base=4 bioconductor-deseq2
...
a LOT
...
Proceed ([y]/n)? 
```

Hit enter to proceed [the default], or `y`. Hit `n` to abort.

 * Specifying r-base=4 is just to ensure the major version of R is 4. It actually installed r-base-4.1.3 for me.
 * Many packages used by deseq2 were also installed (these are dependencies). Some examples: `r-matrix-1.4_1`, `r-nlme-3.1_157`, `r-magrittr-2.0.3`, `r-rcolorbrewer-1.1_3`
 * You will not have to install these dependencies using `install.packages()`

Once created, enter the environment and start R.

```
$ conda activate deseq2
$ R

R version 4.1.3 (2022-03-10) -- "One Push-Up"
Copyright (C) 2022 The R Foundation for Statistical Computing
Platform: x86_64-conda-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> library(DESeq2)
Loading required package: S4Vectors
Loading required package: stats4
Loading required package: BiocGenerics

Attaching package: ‘BiocGenerics’
... more messages ...
>
```

You can now use R like you would on your own computer. This is not RStudio however. It is equivalent to the console pane in RStudio.

You can install new packages from within R: e.g.: `install.packages('dplyr')` It might ask to set a mirror, I choose Kansas or Iowa.

You can also add packages through conda instead, as in `conda install r-dplyr`.

---


#### Slides on Partitions and QoS
https://github.com/ResearchComputing/Partition_and_QoS_Fall_2021

### Stopping your jupyterhub server

It is good practice to stop your server once you are done working.

#### 1. Go to Hub Control Panel

![This is an image](img/HubControlPanel.png)

#### 2. Click Stop My Server

![This is an image](img/StopMyServer.png)

#### 3. Close remaining jupyterhub tabs

Otherwise they will constantly prompt you to restart.
