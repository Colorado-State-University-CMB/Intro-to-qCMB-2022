# Assignment 1

**Due:** January 31, 2022

**Instructions:** 
  * Please turn in the answers to this assignment as a .txt document. To create a .txt document in R, go to **New File**, then select **Text file**. You can use any other text editor if you like. Please do not use Mac's Text Edit application, though.
  * DO NOT include the questions in the document you turn in. Answers only!
  * TURN in your assignment on canvas

-----

## QUESTION 1 (5 pts)

We learned that vectors come in different classes depending on the data type they house. 

```r
users <- c("alvin", "viet", "leila")
logins <- c(12, 5, 34)
```

A. What are the classes of each of these vectors? 

```r
super_vector <- c(users, logins)
```

B. If we merge these vectors together into super_vector, how do their classes change? 

```r
users+1
logins+1
super_vector+1
mean(users)
mean(logins)
mean(super_vector)
```

C. If we try to add 1 to each vector, what happens? If we try to take the mean of each vector, what happens? Explain the results you obtain and why you think this is happening.


-----

## QUESTION 2 (5 pts)

Heather has written some code to create a data frame. Each line of her code has a bug, or error. Explain each error. 

```
languages <- ("English", "Spanish", "Japanese", "French")
_greetings_ <- c("hello", "hola", "ohio", "bonjour")
partings < c("bye", "adios", "mata", "salut")
dictionary <- DataFrame(languages, _greetings_, partings)
len(dictionary)

```

-----

## QUESTION 3 (5 pts)

The following code uses a base R plotting function called **barplot**. 

```r
 barplot(scottish_towns$Population[1:10], names.arg=row.names(scottish_towns)[1:10], las = 2)
 ```

A. Using the help pages for barplot, learn about the option `horiz`. Hack the line of code above to set horiz to "TRUE". What is the new line of code?

B. What does the option `horiz` control?

C. Change the values for the option `las` to `0,1, 2, or 3`. What do they control?

D. What line of code would add a main title to your plot that says "Ten Most Populous Scottish Towns"? Write the full expression.

E. Save your plot as a .pdf and turn it in.




-----

## QUESTION 4 (5 pts)

Let's practice importing some data. Here is a real supplementary dataset that my lab recently published for a manuscript. 

[Table_S4_Signal_to_noise_quantification_table](https://drive.google.com/file/d/1bJy_ELikr5F264xRe-ASNI4iXBVYuxIP/view?usp=sharing)

  * Download the file to your computer.
  * Ensure your working directory is set properly
  * Import the dataset into R using either `read.table()` or `read.csv()` and save it as an object called `signal_to_noise`
  * What is the output of `str(signal_to_noise)`? Copy and paste it here as the answer to this question.

-----

## QUESTION 5 (5 pts)

A. Go to the [R Graph Gallery](https://www.r-graph-gallery.com/index.html). Choose a category of R plots that you would like to learn more about. Using the R Graph Gallery pages, wikipedia, and other internet resources, learn how these plots generate their data. In 2 - 3 sentences, describe the plots in the category you have chosen. Be sure to include what type of data work best in this plot, how the plot is generated (what is on which axis, how measurements are calculated, etc), and what the benefits are of this type of data visualization (for the reader/viewer). 

B. Next, read through some of the R code in the Gallery associated with each plot. You may not understand the R code itself, but just try. You will probably better understand the commenting codes. Answer the following: What functions are being used in this code? Are they part of the base R package or an add on package? If an add-on, what is the suggested add on package? You may need to google search some of this information.

