# Installation notebook

Start a document. It can be a Word document. This will be your notebook for all the installation steps for the different requirements of the class.

First, record information about your OS, followed by the steps you take to install software. This notebook is critical if complications arise, but it's also handy to just have around.

### Example Notebook

[This is my notebook](/Module_0_Installation/R_install_notes.txt) for installing R and RStudio fresh on an MacOS Monterey laptop. It follows the steps below to about line 60.

You can read further to see how I troubleshoot other issues that we won't deal with right away.

## 1 Accounts 

1. At the top of your notebook, state your computer username and whether or not you have administrative permissions, if you don't know- leave this blank until step 3.

2. Then get a github account and log your github username. 

_See the following sections for details._

**Keep your passwords in a separate file.** in case others need to see this notebook.

### Administrative access on your computer

Make sure you know the password to your computer, and have administrative access to install software to it. If you don't know this, you'll find out when you try to install anything (step 3). **If you can't install anything, you need to figure out your admin password.** You can proceed until step 3 while you figure it out.

Just write "admin OK" if you can install. Again, keep write down your password somewhere safe (but not here).

### A Github Account

Get an account at github.com. It's free. Remember your username and password. Put your github username in the notebook.

## 2 Your computer - Operating system version

Note **which version of Mac or Windows you have**. This will be useful when you need to choose a certain program version for your system.

Also, everyone is probably on a 64 bit machine, but Windows people should check.

* Windows: 
  * [Which version of Windows operating system am I running?](https://support.microsoft.com/en-us/windows/which-version-of-windows-operating-system-am-i-running-628bec99-476a-2c13-5296-9dd081cdd808)
  * [How can I tell if my computer is running a 32-bit or a 64-bit version of Windows?](https://support.microsoft.com/en-us/windows/32-bit-and-64-bit-windows-frequently-asked-questions-c6ca9541-8dce-4d48-0415-94a3faa2e13d)
* Mac: 
  * [Find out which macOS your Mac is using](https://support.apple.com/en-us/HT201260)
  * Should be 64 bit (unless your Mac is pre-2007).

## 3 R, RStudio and Swirl

Note where you downloaded from, and the version. 

### R
[Download R](https://cran.rstudio.com/)

### RStudio
[Quicklink to download RStudio](https://www.rstudio.com/products/rstudio/download/#download).
* Windows installer ends with *.exe*
* Mac installer ends with *.dmg*

### Swirl

Find your RStudio console to type in the following commands. [How to find the console.](https://www.google.com/search?q=rstudio+where+is+the+console&oq=rstudio+where+is+the+console)

Type the following into the console and hit return:
```r
install.packages("swirl")
```
In your notes, **record the output of this command**, or note if there isn't any.

Check the installation by typing the following into the console and hitting return.
```r
library("swirl")
```

Again, **record any messages**, or that none were given.
(These steps are based on [the swirl installation guide](https://swirlstats.com/students.html))

## 4. Tex software for creating PDF documents.

### Mac

Download and install MacText from https://tug.org/mactex/

This was very large when I downloaded it: 7GB!


### Windows

Download and install MiKTeX from http://miktex.org

### Create the RMarkdown to PDF

**Go back to RStudio**

Record in your notes where you downloaded the Tex from, and any errors or whether it was successful.

Now try to create a test document from RStudio. **Go to File -> New File -> Rmarkdown...**

You (probably) will be prompted to install a package. Do so. Record in your notebook the package name and any output from that step. 

Now, create a new Rmarkdown file. Leave it "Untitled", but **change the default output format to PDF.** 

It will open a file for editing, but j_ust click on the ball of yarn with the word Knit next to it_. Save it as knit_test. This should open a new window with the document in PDF format. **Record whether you see this document, or any errors that might have popped up in the console.**

If it didn't work, try restarting RStudio and repeat the "Create the RMarkdown to PDF" steps. Record in your notebook if this is necessary.

## 5. Install more R packages.

Now you will install R libaries that we'll use in class. There can be a lot of messages here, so they don't need to go into your notebook. **the thing to record is anything that says "ERROR".** 

### Common messages
 * Warnings. "Warning" is OK and fairly common.
 * `The package blahblah is not available for your version of R.` This doesn't mean what it implies- that there are other versions that _do_ have that package. It most likely just means you mistyped the package name. eyeroll.
 * Out of date packages. `Do you want to update packages that have a newer version? All / Some / None?` **Choose None.**
 * `Do you want to compile from source packages that have a newer version.` **No.**
 * `Choose a local mirror` for download. If you get this, it might list ~100 locations to download from. Look for the one in Kansas (it's slightly closer than Iowa). 


### Find your console again, and do the following:

Record in your notes each packages name and whether it was sucessful, but there can be verbose output to sift through. Just be on the look out for ERRORS (tend to be red).

```r
install.packages("tidyverse")
library(tidyverse)
```

```r
install.packages("BiocManager")
library(BiocManager)
```

To install from Bioconductor, we use a slightly different install command.

```r
BiocManager::install("biomaRt")
library(biomaRt)
```



## 6. Compilers for R source packages

Not now. We'll do #6 as needed.
