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


-----

## QUESTION 4 (5 pts)


-----

## QUESTION 5 (5 pts)

