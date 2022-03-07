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

In fact, you can match _all_ files (except dot-files), with the asterisk alone.

```
wc -l *
    4075 cif_0.cif
   25956 cif_1.cif
   23050 cif_2.cif
    7354 cif_3.cif
   60435 total
```

**Using *less*, choose some fields to grep at in all files**

Once you choose something, you grep it with:

`grep the_pattern *.cif`

**Task**
  1. Rename each **cif** file to a meaninful name, (hint: an accession number). Use `mv orig_name new_name` to rename a file.
  2. Find a field that mentions the species source (only 3 have one).
  3. Use the single grep command, plus the informative field, to save the result to `species.txt`.
  4. Go to https://www.rcsb.org/ to search your accession numbers and validate your investigative work.

**Get sequences for genes**

Try searching the PDB ID's in ncbi. 

*WebAPI for downloading a sequence from NCBI*

#### 1

![google pdb](/Module_6_CommandLine/week2/images/google.png)

#### 2

![google pdb](/Module_6_CommandLine/week2/images/ncbi_PDB.png)

```curl 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=NP_001347554&rettype=fasta&retmode=text'```

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
