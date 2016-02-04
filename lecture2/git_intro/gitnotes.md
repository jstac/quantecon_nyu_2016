# Version Control with Git

![Comic](images/phd101212s.gif)

I suspect that the most of us have participated in some sort of "home-rolled" version of version control. For example, below are three forms of "collaboration" that we might have engaged in

* Tediously "careful" approach: This is when a person (or people) maintain a type of version control by placing the date edited at the end of the file name and then maintaining a full folder of these types of files.
* Stake your claim approach: This appraoch takes the form of calling "dibs" on the use of a specific file for the time being by shouting across the office or sending emails to your collaborators and informing others that the file is "currently under your control" and that others should stay away!
* Where's that email approach: Perhaps worst of all is the approach in which collaborators each maintain their own file and email it back and forth with changes as they make them.

All of these are dangerous workflows and will result in mistakes and unwanted material showing up in your files (or wanted material being erased). Luckily for us, a variety of version control products have evolved in order to help us eliminate the need for these types of behaviors. Dropbox, Google Drives, SVN, Mercurial, and git are all examples of version control systems.

Version control allows us to keep track of what changes have been made over time. Careful maintenance of code and data is vital to reproducibility which "is the hallmark of good science." While Dropbox and Google Drives may be useful for sharing certain types of files (I am not an advocate of their complete abandonment), they are not suitable for the fine tuned version control that we need to maintain good science (or software).

## What is `git`?

Git is a _distributed_ version control system. _Distributed version control_ means that the entire history of every file is kept on your computer. It was originally written by Linus Torvalds (creator of Linux) to help maintain the Linux project (Fun Fact: Git was originally written because Torvalds found all other alternatives of version control to be too slow to manage a project as large as Linux, so he decided to write his own version which he began on 3 April 2005 and started being used on 7 April 2005. Read more about the [history](https://en.wikipedia.org/wiki/Git_(software))).

## Using `git`

## Resources and References

Below is a sequence of references to things that I have found useful. I suggest reading a few of them and at least skimming the majority of them. In particular, the software carpentry lectures on git are very useful for learning the basics -- The Pro Git book is the biblical reference for git and can likely answer any question you will have for at least a few years to come. I have organized both sections by how relevant/convincing I found the documents.

### Git Technical References
[Software Carpentry: Git Lectures](http://swcarpentry.github.io/git-novice/02-setup.html)
[Pro Git](http://git-scm.com/)
[Git Tower](https://www.git-tower.com/learn)

### Why Git References
[Difference between git and Dropbox](https://gist.github.com/magicseth/1951984)
[Version control for scientific research](http://blogs.biomedcentral.com/bmcblog/2013/02/28/version-control-for-scientific-research/)
[Why Physicists Should Use Git Or Why Everyone Should Try Git](http://openmetric.org/assets/slides/whygit/#/)
[Git can facilitate greater reproducibility and increased transparency in science](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3639880/)

