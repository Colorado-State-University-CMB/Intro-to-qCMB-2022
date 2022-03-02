# Review

What does `ls` do?

What is the command to count words, lines, etc? `____`

What flag to the previous command prints *only lines?*

What does the command `file` do?

What output will it give for an executable?

The operator for a pipe in tidyr is `%>%`, what is it in the shell?

The operator to save variable output to the right hand side in R is `->`, what shell operator saves output to a file?

How do you get the *first **10** lines* of a file?

How do you get the *first **5** lines* of a file?

How do you get the *last **3** lines* of a file?

How do you get just the third line?

# New commands

## echo

The command `echo` can print to the screen (or be redirected to a file).

Try:

`echo hello`

`echo hello > greeting.txt`

## history

You can access your command history for your current terminal session.

Try:

`history | tail`

`history | less` How many commands are there?

You can rerun a numbered command from the prompt. (Make sure you have quit out of `less`)

``` bash
david@Cumbernauld % !1011
ls
Assignment.html Assignment.md README.md
```

## grep

Stands for global regular expression and print. (See a hardcore explanation of its origin with [this video](https://www.youtube.com/watch?v=NTfOnGZUZDk).)

Kind of synonymous with tidyr `filter` accept that it doesn't specify a column, it matches lines by string matching.

``` bash
david@Cumbernauld % grep elt promoters.hilo.classC.bed

chrX 10480240 10481439 1200 + 1200 WBGene00001250 elt-2 13939.048604 29186.796875 13.7854357493514 14.8371145726847 NA NA NA NA NA NA NA 6074.28154993418 6.7355897659042 0.569852085237524 1.07416954319249e-33 1.85175905488319e-31 classC

chrX 10474527 10475726 1200 + 1200 WBGene00001252 elt-4 3684.935242 6905.955078 11.9165225519833 12.7708169132894 NA NA NA NA NA NA NA 40.6441481292203 4.16556268383076 0.868313775372244 3.13136982534824e-07 3.47508592401713e-06 classC
```

A more specific search with `elt-2` instead of `elt`:

`david@Cumbernauld % grep elt-2 promoters.hilo.classC.bed`

What if we change the pattern to `clec`?

## cat

Concatenate (stitch together) files.

`david@Cumbernauld % cat promoters.hilo.classC.bed promoters.hilo.classC.bed > double.bed`

## diff

`david@Cumbernauld % diff file1 file2`

# In-class Exercise/Homework

**How to fill this out**

You will execute the commands in this document and then save your work with `history > assignment.txt`

*You will turn in assignment.txt.*

If there is a specific question:

Q: What is blah blah blah?

A: type `echo your answer`

Otherwise, the question might ask you to come up with, and execute a command.

You still will have an A: for answer, but just execute the command (the `history > assignment.txt`) will capture the work.

**Turn in assignment.txt.**

## 1: Basics

### 1.1: Running a terminal

#### Exercise: 1 {#exercise-2}

Q: By examining the menu items for your terminal program, figure out how to create a new tab. *What is the keyboard shortcut for creating a new tab.*

A:

------------------------------------------------------------------------

## 1.3: Man pages

### Exercise 1 {#exercise-1}

Q: According to the man page, what is the official description of **echo** on your system?

A:

### Exercise 2 {#exercise-2-1}

Q: By reading the man page for echo, determine the command needed to print out "hello" without the trailing newline, and verify using your terminal that it works as expected.

A:

------------------------------------------------------------------------

## 1.2: More practice

### Exercise: 1 {#exercise-1-1}

Q: Use the command "echo" to print the string "hello, world".

A:

### Exercise 2 {#exercise-2-2}

Q: Type the command echo 'hello (with a mismatched single quote), and then get out of trouble using a ctrl- key combination. What is it?

A:

### Exercise 3 {#exercise-1-2}

Q: Using the up arrow, print to the screen the strings "fee", "fie", "foe", and "fum" without retyping echo each time.

A:

------------------------------------------------------------------------

## 1.4: Practice bailing with ctrl-c

### Exercise 1 {#exercise-2-4}

Q: By running man sleep, figure out how to make the terminal "sleep" for 5 seconds, and execute the command to do so. After waiting the requisite 5 seconds, execute the command to sleep for 5000 seconds, realize that's well over an hour, use `ctrl-c` to get out of trouble.

A:

------------------------------------------------------------------------

# 2: Manipulating files

## 2.1: Redirecting and appending

### Exercise 1 {#exercise-1-4}

Q: Using echo and \>, make files called line_1.txt and line_2.txt containing the first and second lines of Sonnet 1, respectively.

`echo "From fairest creatures we desire increase," > line_1.txt` , `echo "That thereby beauty's rose might never die," > line_2.txt`

A:

### Exercise 2 {#exercise-2-5}

Q: Create sonnet_1.txt (containing the first two lines of the sonnet) by first redirecting the contents of line_1.txt and then appending the contents of line_2.txt. `cat line_2.txt >> sonnet_1.txt.`

`echo "From fairest creatures we desire increase," > sonnet_1_copy.txt` , `echo "That thereby beauty's rose might never die," >> sonnet_1_copy.txt` , `diff sonnet_1.txt sonnet_1_copy.txt`

A:

### Exercise 3

Q: Use cat to combine the contents of line_1.txt and line_2.txt in reverse order using a single command, yielding the file sonnet_1\_reversed.txt.

A: `cat line_2.txt line_1.txt > sonnet_1_reversed.txt`

## 2.2: Listing

### Exercise 1 {#exercise-1-5}

Q: What's the command to list all the files (and directories) that start with the letter "s"?

A: `ls s*`

### Exercise 2 {#exercise-2-6}

Q: What is the command to list all the files that contain the string "onnet", long-form by reverse modification time?

A: `ls -rtl *onnet*`

### Exercise 3 {#exercise-3-1}

Q: Using the man page, find the command to list all files (including hidden ones) by reverse modification time, in long form.

A:

------------------------------------------------------------------------

## 2.3: Renaming, copying, deleting

### Exercise 1 {#exercise-1-6}

Q: Use the echo command and the redirect operator \> to make a file called foo.txt containing the text "hello, world". Then, using the cp command, make a copy of foo.txt called bar.txt. Using the diff command, confirm that the contents of both files are the same.

A:

### Exercise 2 {#exercise-2-7}

Q: By combining the cat command and the redirect operator \>, create a copy of foo.txt called baz.txt without using the cp command.

A:

### Exercise 3 {#exercise-3-2}

Q: Create a file called quux.txt containing the contents of foo.txt followed by the contents of bar.txt.

A:

### Exercise 4 {#exercise-4}

Q: How do rm nonexistent and rm -f nonexistent differ for a nonexistent file?

A:

## 2.4: Summary

### Exercise 1 {#exercise-1-7}

Q: Use echo to make a file called sonnets.txt containing the full (original) text of Shakespeare's first sonnet.

`echo "FRom faireſt creatures we deſire increaſe, That thereby beauties Roſe might neuer die, But as the riper ſhould by time deceaſe, His tender heire might beare his memory: But thou contracted to thine owne bright eyes, Feed’ſt thy lights flame with ſelfe ſubſtantiall fewell, Making a famine where aboundance lies, Thy ſelfe thy foe,to thy ſweet ſelfe too cruell: Thou that art now the worlds freſh ornament, And only herauld to the gaudy ſpring, Within thine owne bud burieſt thy content, And tender chorle makſt waſt in   Pitty the world,or elſe this glutton be,    To eate the worlds due,by the graue and thee." > sonnets.txt`

`cat sonnets.txt`

`file sonnets.txt`

A:

### Exercise 2 {#exercise-2-8}

Q: Type the sequence of commands needed to create an empty file called foo, rename it to bar, and copy it to baz.

A:

### Exercise 3 {#exercise-3-3}

Q: What is the command to list only the files starting with the letter "b"?

A:

### Exercise 4 {#exercise-4-1}

Q: Remove both bar and baz using a single call to rm.

A:

------------------------------------------------------------------------

# 3: Inspecting files

## 3.1: Downloading a file

### Exercise 1 {#exercise-1-8}

Q: Use the command curl -I www.learnenough.com to fetch the HTTP header for the Learn Enough website. What is the HTTP status code of the address? How does this differ from the status code of learnenough.com?

A:

### Exercise 2 {#exercise-2-9}

Q: Using ls, confirm that sonnets.txt exists on your system. How big is it in bytes?

A:

### Exercise 3 {#exercise-3-4}

Q: Using the -h ("human-readable") option to ls, list the long form of the sonnets file with a human-readable byte count.

A:

### Exercise 4 {#exercise-4-2}

Q: Suppose you wanted to list the files and directories using human-readable byte counts, all, by reverse time-sorted long-form. Find and execute this command:

A:

## 3.2: Making heads and tails of it

### Exercise 1 {#exercise-1-9}

Q: By piping the results of tail sonnets.txt through wc, confirm that (like head) the tail command outputs 10 lines by default.

A:

### Exercise 2 {#exercise-2-10}

Q: By experimenting with different values of n, find a head command to print out just enough lines to display the first sonnet in its entirety

A:

### Exercise 3 {#exercise-3-5}

Q: Pipe the results of the previous exercise through tail (with the appropriate options) to print out only the 14 lines composing Sonnet 1.

A:

### Exercise 4 {#exercise-4-3}

Q: To simulate the creation of a log file, run ping learnenough.com \> learnenough.log in one terminal tab. (The ping command "pings" a server to see if it's working.) In a second tab, type the command to tail the log file. (At this point, both tabs will be stuck, so once you've gotten the gist of tail -f you should use the technique from Box 4 to get out of trouble.)

`ping learnenough.com > learnenough.log` , `tail -f learnenough.log` , `ctrl-c`

------------------------------------------------------------------------

## 3.3: Less is more

### Exercise 1 {#exercise-1-10}

Q: Run less on sonnets.txt. Go down three pages and then back up three pages. Go to the end of the file, then to the beginning, then quit. *Practice only, no answer/output.*

A:

### Exercise 2 {#exercise-2-11}

Q: Search for the string "All" (case-sensitive). Go forward a few occurrences, then back a few occurrences. Then go to the beginning of the file and count the occurrences by searching forward until you hit the end. Compare your count to the result of running grep All sonnets.txt \| wc.

A:

### Exercise 3 {#exercise-3-6}

Q: Using less and / ("slash"), find the sonnet that begins with the line "Let me not". Are there any other occurrences of this string in the Sonnets?

A: `No output, just practice.`

### Exercise 4 {#exercise-4-4}

Q: By searching for the string "sort" in the man page for ls, discover the option to sort files by size. What is the command to display the long form of files sorted so the largest files appear at the bottom?

A:

## 3.4: Grepping

### Exercise 1

Q: By searching man grep for "line number", construct a command to find the line numbers in sonnets.txt where the string "rose" appears.

A:

### Exercise 2

Q: You should find that the last occurrences of "rose" is (via "roses") on line 2203. Figure out how to go directly to this line when running less sonnets.txt.

A:

### Exercise 3

Q: By piping the output of grep to head, print out the first (and only the first) line in sonnets.txt containing "rose".

A:

### Exercise 4

Q: In Listing 16, we saw two additional lines that case-insensitively matched "rose". Execute a command confirming that both of the lines contain the string "Rose" (and not, e.g., "rOSe").

A:

### Exercise 5 {#exercise-5}

Q: Write a command confirming that the number of lines matching "Rose" but not matching "rose" is equal to the expected 2.

A:

------------------------------------------------------------------------

## 3.5: Summary

### Exercise 1

Q: Pipe history to less to examine your command history. What was your 17th command?

A:

### Exercise 2

Q: By piping the output of history to wc, count how many commands you've executed so far.

A:

### Exercise 3

Q: By piping the output of history to grep, determine the number for the last occurrence of curl.

A:

### Exercise 4

Q: Use the result from the previous exercise to re-run the last occurrence of curl.

A:

### Exercise 5

# 4: Directories

## 4.1: Structure

### Exercise 1

Q: Write in words how you might speak the directory \~/foo/bar.

A:

### Exercise 2

Q: In /Users/bill/sonnets, what is the home directory? What is the username? Which directory is deepest in the hierarchy?

A:

### Exercise 3

Q: For a user with username bill, how do /Users/bill/sonnets and \~/sonnets differ (if at all)?

A:

------------------------------------------------------------------------

## 4.2: Making directories

### Exercise 1

Q: What is the option for making intermediate directories as required, so that you can create, e.g., \~/foo and \~/foo/bar with a single command?

A:

### Exercise 2

Q: Use the option from the previous exercise to make the directory foo and, within it, the directory bar (i.e., \~/foo/bar) with a single command.

A:

### Exercise 3

Q: By piping the output of ls to grep, list everything in the home directory that contains the letter "o".

A:

## 4.3: Navigating directories

### Exercise 1

Q: How do the effects of cd and cd \~ differ (or do they)?

A:

### Exercise 2

Q: Change to text_directory, then change to second_directory using the "one directory up" double dot operator:

A:

### Exercise 3

Q: From wherever you are, create an empty file called nil in text_directory using whatever method you wish.

A:

### Exercise 4

Q: Remove nil from the previous exercises using a different path from the one you used before. (In other words, if you used the path \~/text_directory before, use something like ../text_directory or /Users/<username>/text_directory.)

A:

------------------------------------------------------------------------

## 4.4: Renaming, copying, deleting directories

### Exercise 1

Q: Make a directory foo with a subdirectory bar, then rename the subdirectory to baz.

A:

### Exercise 2

Q: Copy all the files from `text_files` into foo. Use `man cp` to find a flag that copies them all at once, i.e.: `cp _flag_ ../text_files _dest_`

### Exercise 3

Q: Using the same flag as above, add "Copy all the files in text_files, without directory, into bar. Add

A: `cp -r ../text_files/ bar`

### Exercise 4

Q: Remove foo and everything in it using a single command.

A: `rm -rf foo`

## 
