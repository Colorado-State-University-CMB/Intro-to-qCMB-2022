# Introduction to the command line

## What is the command line?

User interfaces to the operating system (kernel) are considered a (shell). These include the windowing systems that we're used to, called graphical user interfaces
(GUIs).  The command line interface is a text-based shell where you launch programs, handle files, and manage the computer using typed commands.

## What does it look like?

The command line consists of 1) the prompt and 2) blank space for you to type your command. The prompt is usually terminated by a single character.

```
% command1 command2
$ command1 command2
```

By default, most setups make it so the prompt shows useful information, such as your present working directory, username, computer name, etc.

```bash
david@Cumbernauld /Users/David % cd /
david@Cumbernauld / % whoami
david
david@Cumbernauld / % 
```

## The command line is closely tied to the filesystem

Using the command line makes you learn your file organization and navigation together.

A Mac user directory:

```
/Users/
       yournamehere/
              Desktop/
              Documents/
              Downloads/
```

Windows is similar, but the _path seperator is a backslash (**\\**)_.
```
C:\users\
      yournamehere\
```

### Navigation

`ls` List your directory contents. For Ubuntu users, this may be blank, as you have just started using it.

_however_,

`ls -a` should reveal more files. These files all start with a '.' and are called _dot files._ They are hidden by convention in order to declutter the output of `ls`. They usually contain configuration issues.

`cd destination` takes you to a place in the file hierarchy. From the top level, you can go "down" in the hierarchy by supplying just the directory/folder name you wish to travel to.

`pwd` - present working directory. In addition to being displayed by the prompt, you can call up your present working directory explicitly.

```
david@Cumbernauld /Users/David % pwd
/Users/David
david@Cumbernauld /Users/David % cd /
david@Cumbernauld / % pwd
/
```


Starting at root (slash), I will navigate back to my home directory, using `cd`, `pwd` and `ls`
```
david@Cumbernauld / % pwd
/
david@Cumbernauld / % ls      
Applications System       Volumes      cores        etc          opt          sbin         usr
Library      Users        bin          dev          home         private      tmp          var
david@Cumbernauld / % cd Users
david@Cumbernauld /Users % ls
Shared david
david@Cumbernauld /Users % cd david
david@Cumbernauld /Users/david % pwd
/Users/david
```

#### Relative path

You can go "up" using the `..` shortcut.

```
david@Cumbernauld /Users/david % pwd
/Users/david
david@Cumbernauld /Users/david % cd ..
david@Cumbernauld /Users % pwd
/Users
david@Cumbernauld /Users % cd ..
david@Cumbernauld / % pwd
/
```

## Terminus

Play the in-browser game [Terminus](https://web.mit.edu/mprat/Public/web/Terminus/Web/main.html) to get use to file hiearchy navigation. 
It is especially good at getting you into the habit of ls'ing after cd'ing.

## Command reference

### Navigation
<dl>
       <dt>ls</dt>
       <dd>list directory contents. Use `-a` to list "all", including the dot-files.</dd>
       <dt>cd</dt>
       <dd>Change directories. Use `cd ..` to navigate upward in the directory tree. Use `cd` by itself to go to your home directory.</dd>
       <dt>pwd</dt>
       <dd>Print Working Directory.</dd>
</dl>

### Inspection
<dl>
       <dt>head</dt>
       <dd>Like `head` in R, except you supply the number of lines as a separate word like `head -5`, `head -n 5`, or `head -n5`</dd>
       <dt>tail</dt>
       <dd>Like head, but tail.</dd>
       <dt>less</dt>
       <dd>Page through a document</dd>
       <dt>file</dt>
       <dd>Read the file type of a file.  For example, `file textfile` might give `ASCII text`, whereas `file /bin/cp` might give `Mach-O 64-bit executable x86_64`
       </dd>
</dl>
              
           
