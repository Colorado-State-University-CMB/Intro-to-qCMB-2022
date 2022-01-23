
<img src="webContent/RBanner1.jpeg" width="600">

# R Basics

Let's get started!

January 24, 2022

## Lessons for today

 - What is R?
 - What is RStudio
 - Interacting with R within RStudio
 - Meet the different R objects
   - Values & Numbers
   - Vectors
   - Dataframes
   - Matrices

-----


# What is R?



## R is a programming language
 * designed for statistics, plotting, & data analysis
 * comprised of a base package that can be functionality expanded through add-on packages
 * open-source
 * an interpreted language

## R is an ecosystem
 * the base package - R
 * the extended environment - packages, RStudio, documentation
 * the repository - CRAN (Comprehensive R Archive Network)
 * the governance, support, & organization - R Core Development Team, R Foundation, R Consortium, RStudio
 * the contribution community - developers
 * the user community - us



## How did R come to be?


<img src="webContent/Founders.jpeg" width="600">

  * 1976 - Bell Labs develops S, an open-source statistical programming language. The team was led by John Chambers.
  * 1991 - Ross Ihaka and Robert Gentleman at the University of Aukland in New Zealand embarks on a research venture to adapt S into R. 
    * Motivated to create a better user experience for academics
  * 1995 - first official beta release
  * 1997 - CRAN (Comprehensive R Archive Network) distributes R packages
  * 1997 - R Core Team, a governing body, is founded
  * 2000 - R 1.0.0 is released 
  * 2011 - R versions start getting funny nicknames. R 2.14.0 Great Pumpkin
  * 2010 - CRAN exceeds 10,000 published packages!
  * 2021 - R's latest version is 4.1.2 Bird Hippie



>  “R changed my opinion of humanity to some extent, to see how people are really willing to freely give of themselves and produce something larger than themselves without any thought of personal glory. There’s a lot of work with no recognition.” - Ross Ihaka

## Why do we use R in computational biology?

Why do we use any programming language in biology? 
  * To increase our speed & efficiency
  * To increase our scale
  * To promote reproducibility
  * To save money
  * To allow for transparency

Why use R in particular? What are its benefits for the life sciences?
  * Statistics is the heart of R
  * Packages extend R for biology-specific tasks
  * We can benefit from a wide, supportive, and dedicated community of users and developers who have created stereotyped documentation
  * R is beautiful

-----

# What is RStudio?

<img src="webContent/Screen Shot RStudio.png" width="700">

RStudio is an **Integrated Development Environment** or **IDE**. An IDE is a software application that allows programmerss to develop software within an organized user space. Think of it as a dedicated application in which we can interact with a programming language to generate organized software projects. IDE's combine a lot of useful tools all into one application with an organized layout.

:heavy_exclamation_mark: **EXERCISE: Exploring RStudio**

  * :arrow_right: Open R Studio. What does it look like?
  * It has a menu bar on the top and a series of windows.
  * **Console panel** - this is where R is loaded and ready to go. It shows you the version of R you are running. A prompt, which is denoted by the ">" character shows you that R is ready for your input. This is where we will **interact** with R. More on that later.
  * **Environment panel** - this is a list that will populate with **objects** as we create them in R. More on objects later. 
  * **Plots panel** - this is where plots and graphics will appear when we generate them.
  * _But wait! That's not all. We can toggle through other panels behind the Console, Environment, and Plot panels that are shown as defaults._
  * **Files panel** - This gives you the ability to browse and interact with the file structure on your own local computer.
  * **Help panel** - This will pop up with useful help information when you execute the help command. More on that later.
  * _Finally, there is one panel not shown upon opening RStudio, and that is the Text Editor_
  * **Text editor** - The text editor is a panel that allows you to write a script of R commands into a document that can be saved for later, executed one line at a time, or executed all at once.
  * :arrow_right: Let's start a new script so we can see the text editor.
    * :arrow_right: Go to **File** on the top menu bar
    * :arrow_right: Go to **New File**
    * :arrow_right: Select **R Script**

**Graphical re-cap:**

<img src="webContent/WebContent_Powerpoint_interface.jpg" width="1200">

:heavy_exclamation_mark: **BEST PRACTICES: Closing RStudio**

When RStudio closes, by default, it will prompt you whether you want to "Save workspace image to ~/.Rdata". I recommend selecting **Don't Save**. Otherwise, the next time R starts up, your old objects will already be there waiting for you. It's better form to start from scratch with each new session. However, when working on very long projects, I will admit, this is a useful feature. 

  * What to do? 
  * :arrow_right: Either remember to select **Don't Save** each time. 
  * :arrow_right: OR, you can change this setting by going to **RStudio** on the top menu, then **Preferences ...**, then changing **"Save workshopspace to .RData on exit"** from **Ask** to **Never**. It'll look like so:

<img src="webContent/WebContent_Powerpoint_RData.jpg" width="900">


-----

# Interacting with R within RStudio


