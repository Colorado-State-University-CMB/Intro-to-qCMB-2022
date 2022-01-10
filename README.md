# Introduction to Quantitative Cell and Molecular Biology, 2022 (CM580A3)

![Bioinformatics Image](/dna-g4efa38871_1920.jpeg)

## Course material will be available here.

## For Instructors

<details><summary>Knitting RMarkdown documents</summary>


### TeX for knitting to PDF
* MiKTeX on Windows
* MaCTeX 2013+ on Mac.

### Too much output (max.print)
Many students last year had an RStudio with a default `max.print` of 10000... which led to documents being turned in that were hundreds of pages. It can be easily handled in the setup chunk.

````r
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

options(max.print=100)
```
````
</details>
