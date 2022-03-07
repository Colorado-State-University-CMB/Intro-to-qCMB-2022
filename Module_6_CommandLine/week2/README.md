# Week 2 - command line

## 1. Detective work

Say you inherit someone's project and they have downloaded a bunch of CIF (protein structure files) that aren't named correctly. Luckily, they are text format, and you can read them easily.

**Create a directory structure in your working directory.**

```bash
proj01/
	pdb/
		cif/
	loci/
		bed/
		gtf/
		fasta/

```

**Download and expand CIFs.tgz**

Go to CIFs.tgz in this github directory and download them. Move the file into `proj01/pdb/cif`.

`tar -zxvf CIFs.tgz` expands the file, creating a new directory with cif files in them. See [Command reference](https://github.com/Colorado-State-University-CMB/Intro-to-qCMB-2022/edit/main/Module_6_CommandLine/week2/README.md#command-reference) for more on [tar](https://github.com/Colorado-State-University-CMB/Intro-to-qCMB-2022/edit/main/Module_6_CommandLine/week2/README.md#tar).

**Figure out a field to *grep* on**

Find a way to identify what locus the file belongs to by inspecting it (*head*, *less*).

**Wildcards**

`cd` to your *cif* directory (`pwd` should end with `proj01/pdb/cif`).

You can "autofill" the command line with pattern matching. Example:

```
wc -l *.cif
    4075 cif_0.cif
   25956 cif_1.cif
   23050 cif_2.cif
    7354 cif_3.cif
   60435 total
```

### Command reference

#### tar 

**T**ape **ar**chive

Use `man tar` to figure out the flags:
<dl>
  <dt>-z</dt>
  <dd>Compress/uncompress</dd>
  <dt>-f</dt>
  <dd>Specify file.</dd>
  <dt>-x</dt>
  <dd>...</dd>
  <dt>-v</dt>
  <dd>...</dd>
<dl>
