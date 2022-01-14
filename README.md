# Introduction to Quantitative Cell and Molecular Biology, 2022 (CM580A3)

![Bioinformatics Image](/images/dna-g4efa38871_1920.jpeg)

## Course material is available here.

1. [Module 1 - R and RStudio](Module_1_RStudio/README.md)
1. [Module 2 - RMarkdown](Module_2_RMarkdown/README.md)
1.
1.
1.


## For Instructors

<details><summary>Knitting RMarkdown documents</summary>

### Too much output (max.print)
Many students last year had an RStudio with a default `max.print` of 10000... which led to documents being turned in that were hundreds of pages. It can be easily handled in the setup chunk.

````r
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

options(max.print=100)
```
````
</details>
