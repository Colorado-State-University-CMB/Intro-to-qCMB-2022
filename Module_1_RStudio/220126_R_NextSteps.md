# R Next Steps

January 26, 2022

## Lessons for today

  * [Functions - the verbs of the R language](#functions---the-verbs-of-the-r-language)
  * Import and Export of data
  * Packages expand R
  * Plotting - R is beautiful

## Useful references


----

# Functions - the verbs of the R language

By now, you are already experienced in calling R functions. Here are some examples:

```r

date()

dim(model_systems)

mean(chromosomes)

```

Functions typically involve parentheses. Sometimes these parentheses are empty and sometimes we type something in them. Think of this similar to how some verbs take direct objects (I threw the ball) and some don't (I run). 

**ARGUMENTS**: In the case of R, this type of required input is called an **argument**. 

**OPTIONS**: R functions can also take optional inputs called **options**.

The `help()` function gives us information about each function. Place the function name within the parentheses:

➡️ Try the help function

```r
help(dim)
```

  * The `x` in the example tells you that this function takes an argument. 

➡️ Let's look at the help menu for `mean` and see how options are specified.

```r
help(mean)
```

The help menu species the **defalt** values that are set for each option. As a default, trim is set to 0. In other words, all the values are used to calculate the mean. However, this value can be changed to 0.2, in which case, the most extreme 20 % of datapoints will be removed before the mean is calculated. 

➡️ Give it a try:

```r
mean(chromosomes, trim = 0.2)
```

⚠️ **SUMMARY** 


## Packages

## Import & Export of Data

## Plotting




