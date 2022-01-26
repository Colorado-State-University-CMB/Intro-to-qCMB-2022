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

Two different types of information can go in the parentheses: **arguments** and **options**.

  * **ARGUMENTS**: These are object names that are added inside the parentheses. They are what the function will operate on. They are analogous to direct objects in English. 
  * **OPTIONS**: Think of these as adverbs. These are optional content you can add that changes **how** the function will operate.

The `help()` function gives us information about how a function operates. The help page will tell us
   * What a function does
   * Whether it requires an argument
   * what object classes are allowed as arguments
   * The potential list of options
   * The default values associated with each option

➡️ Try the help function

```r
help(dim)
```

  * The `x` in the example tells you that this function takes an argument. If we read under **ARGUMENTS** we can learn whith object classes are allowed.

➡️ Let's look at the help menu for `mean`.

```r
help(mean)
```

This help menu also species the **options** that the mean function takes and their **default** values. As a default, trim is set to 0. In other words, all the values are used to calculate the mean. However, this value can be changed to 0.2, in which case, the most extreme 20 % of all datapoints will be removed before the mean is calculated. 

➡️ Give it a try:

```r
mean(chromosomes, trim = 0.2)
```

⚠️ **GRAPHICAL SUMMARY** 

<img src="webContent/WebContent_Powerpoint_functionGrammar.jpg" width="600">


## Import and Export of Data

So far, we've created objects by assignment expressions that directly specify their values. Next, we'll learn how to **import** data into R through an special assignment expression.

First, let's download a dataset to import. 

To do this, we will first need to learn a little bit about how to navigate file structures within R.

To determine where R "thinks it is" on your computer, use the command `getwd()` for **get working directory**.

```r
getwd()
```

For today, we're going to make things easy. We will set our working directory to the directory that contains our dataset.

➡️ Go to the **Files** Panel of RStudio.

➡️ Navigate to the location containing the downloaded dataset.

➡️ Change the working directory by going to the **Files Menu Banner**, selecting **More**, and selecting **Set As Working Directory**



## Plotting

## Packages



