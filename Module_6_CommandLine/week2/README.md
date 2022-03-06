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
