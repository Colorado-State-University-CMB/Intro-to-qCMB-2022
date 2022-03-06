# Week 2 - command line

## 1. Detective work

Say you inherit someone's project and they have downloaded a bunch of CIF (protein structure files) that aren't named correctly. Luckily, they are text format, and you can read them easily.

**Download and expand CIFs.tgz**

`tar -zxvf CIFs.tgz` expands the file, creating a new directory with cif files in them. See [Command reference](https://github.com/Colorado-State-University-CMB/Intro-to-qCMB-2022/edit/main/Module_6_CommandLine/week2/README.md#command-reference)

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
