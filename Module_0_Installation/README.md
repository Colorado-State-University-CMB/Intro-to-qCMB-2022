# Installation notebook

Start a document. It can be a Word document. This will be your notebook for all the installation steps for the different requirements of the class.

First, record information about your OS, followed by the steps you take to install software. This notebook is critical if complications arise, but it's also handy to just have around.

## 1. Operating system version

Note **which version of Mac or Windows you have**. This will be useful when you need to choose a certain program version for your system.

Also, everyone is probably on a 64 bit machine, but Windows people should check.

* Windows: 
  * [Which version of Windows operating system am I running?](https://support.microsoft.com/en-us/windows/which-version-of-windows-operating-system-am-i-running-628bec99-476a-2c13-5296-9dd081cdd808)
  * [How can I tell if my computer is running a 32-bit or a 64-bit version of Windows?](https://support.microsoft.com/en-us/windows/32-bit-and-64-bit-windows-frequently-asked-questions-c6ca9541-8dce-4d48-0415-94a3faa2e13d)
* Mac: 
  * [Find out which macOS your Mac is using](https://support.apple.com/en-us/HT201260)
  * Should be 64 bit (unless your Mac is pre-2007).

## 2. R, RStudio and Swirl

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

## 3. Tex software for creating PDF documents.

### Mac

Download and install MacText from https://tug.org/mactex/


### Windows

Download and install MiKTeX from http://miktex.org

### Back to RStudio

Record in your notes where you downloaded the Tex from, and any errors or whether it was successful.

Now try to create a test document from RStudio. **Go to File -> New File -> Rmarkdown...**

You (probably) will be prompted to install a package. Do so. Record in your notebook the package name and any output from that step. 

Now, create a new Rmarkdown file. Leave it "Untitled", but **change the default output format to PDF.** 

It will open a file for editing, but j_ust click on the ball of yarn with the word Knit next to it_. Save it as knit_test. This should open a new window with the document in PDF format. **Record whether you see this document, or any errors that might have popped up in the console.**

## 4. Compilers for R source packages

