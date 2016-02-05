# Version Control with Git

<div style="max-width: 600px; margin:0 auto;">
  <img src="images/phd101212s.gif" style="max-width:100%;"/>
</div>

I suspect that the most of us have participated in some sort of "home-rolled" version of version control. For example, below are three forms of "collaboration" that we might have engaged in

* Tediously "careful" approach: This is when a person (or people) maintain a type of version control by placing the date edited at the end of the file name and then maintaining a full folder of these types of files.
* Stake your claim approach: This appraoch takes the form of calling "dibs" on the use of a specific file for the time being by shouting across the office or sending emails to your collaborators and informing others that the file is "currently under your control" and that others should stay away!
* Where's that email approach: Perhaps worst of all is the approach in which collaborators each maintain their own file and email it back and forth with changes as they make them.

All of these are dangerous workflows and will result in mistakes and unwanted material showing up in your files (or wanted material being erased). Luckily for us, a variety of version control products have evolved in order to help us eliminate the need for these types of behaviors. Dropbox, Google Drives, SVN, Mercurial, and git are all examples of version control systems.

Version control allows us to keep track of what changes have been made over time. Careful maintenance of code and data is vital to reproducibility which "is the hallmark of good science." While Dropbox and Google Drives may be useful for sharing certain types of files (I am not an advocate of their complete abandonment), they are not suitable for the fine tuned version control that we need to maintain good science (or software).

## What is `git`?

Git is a _distributed_ version control system. _Distributed version control_ means that the entire history of every file is kept on your computer. It was originally written by Linus Torvalds (creator of Linux) to help maintain the Linux project (Fun Fact: Git was originally written because Torvalds found all other alternatives of version control to be too slow to manage a project as large as Linux, so he decided to write his own version which he began on 3 April 2005 and started being used on 7 April 2005. Read more about the [history](https://en.wikipedia.org/wiki/Git_(software))).

## Using `git`

Many people associate `git` with the cloud based repository service Github, but `git` can be run independently either just on your own computer or on a self-hosted server. In this short tutorial, we are going to create a `git` repository hosted on our computer. We will then talk about some of the day to day commads that will be used in `git`.

WARNING: Until you know what you're doing and exactly how they work, NEVER use the `-f` or `--force` flags no matter what the internet tells you.

We are now going to walk through some of the basic settings you should set and an example of some commands.

### Configuration

Here we deal with configuration details such as our default editor, user name, email, colors, etc...

* Set name: `git config --global user.name FirstName LastName`
* Set email: `git config --global user.email email@email.com`
* Sets git colors: `git config --global color.ui "auto"`
* Sets editor to vim: `git config --global core.editor "vim"`

### Creating a folder

We will create a folder called `<MyFirstGitRepo>` using `mkdir <MyFirstGitRepo>`

### Initializing a `git` repository

Now enter that directory using `cd <MyFirstGitRepo>`. To initialize this directory as a `git` repository (which in the background will create a series of directories and files) we will use the command `git init`.

We can see the things that were created by typing `ls .git`, but don't worry too much about what is being kept inside yet.

### Four Stages of Files

Files in a git repository can be in four different stages: untracked, unstaged, staged, and committed.

* Untracked: This is a new file that your repository hasn't seen before.
* Unstaged: A file that has previously been tracked and saved, but has changed since its last version.
* Staged: A file that has been changed and is prepared to be committed. Nothing is set in stone yet though and new changes can be made.
* Committed: When a file is committed it becomes a piece of the history of the repository. This moment of time in the file will be able to be referenced or referred to in the future.

The picture below illustrates this "life cycle"

<div style="max-width: 600px; margin:0 auto;">
<img src="images/gitlifecycle.png" style="max-width:100%;"/>
</div>

Let's illustrate this through an example and to expose ourselves to the commands that will be helpful: `git add`, `git commit`, `git diff`, and `git status`.

First let's check what is new in our repository. Type `git status`. What does it say? It should say that there is nothing to commit because we haven't done anything yet.

Now let's create a file called `README.md`. Open this file and type your name in it (could also use the command `echo "FirstName LastName" >> README.md`). If you type `git status` now, what do you see? It should list `README.md` as an untracked file.

We can move this file from untracked to staged by using `git add README.md`. Type this command and then once again type `git status`. Our `README.md` file now shows up in green as a change to be committed. This means it is staged.

We can take our "snapshot" of the file by using the command `git commit -m "Type a meaningful msg here"`. Commit the file and then once again type `git status`. What do you see now? The repository should be clean again.

Add something new to the `README.md` file. Once again, we can check the status of our git repository by typing `git status`. It will tell us that `README.md` is unstaged because we have made changes to it.

Imagine we wanted to double check the things that we changed. We could type `git diff README.md` which will show us the changes that have been made to that file since its previous commit.

These commands will be the core of your `git` workflow so I suggest familiarizing yourself with what they do.

### Commit History

Remember how we can leave ourselves commit messages. It is useful to leave meaningful commit messages because they are left as a guide for yourself. We can see our history of commits by typing `git log`. We can do this in several formats: Try `git log --pretty=oneline`, `git log --stat`, `git log --since=2.weeks`, etc... See the Git Pro book for more options.

If we want to return to a previous commit in our history we can use the information from `git log` to get back. There is a sequence of numbers and letters, which we will call commit hash, at the top of each entry in your history. If you copy this and type `git reset <commit hash>` then it will return us to that moment of our history. All of our more recent changes will be there, but they will be as if they had been just staged. Can use `git reset <commit hash> --hard` to reset to a previous point in time and delete all changes, but be very careful with that command as it can erase your history (in fact, I suggest not using it until you really know enough that you know it is what you want).

### Branching

A branch is essentially an specific version of your folder. You start with a main branch which is called master. When you create another branch, it is an exact replica of the branch it is being created from (typically master) and includes the full history of the repository. After creating a new branch you can make changes and this new branch will develop its own history (without changing the history of the original branch -- such as master). If you decide you like the changes that you made then you can bring them into the original branch.

The command `git checkout` is used for a few different things, but we will mostly use it for switching between branches. When used with the `-b` flag, it creates a new branch and switches to it.

Let's create a new branch in our repository. `git checkout -b test`. What has changed? Let's look at the output of `git branch`. Now let's make some changes in our branch.

```
> echo "New line" >> README.md
> git add README.md
> git commit -m "Branch update"
> git log
```

We can see the history has changed and we have a new commit. Let's return to the master branch by typing `git checkout master`. Now check `git log`, notice that we no longer have the changes from our test branch.

Imagine that we decided to bring those changes into our master branch. We could use the `git merge` command by typing `git merge test` which will merge the branch `test` into the current branch.

### Ignore Files

Sometimes there are files that get created via an intermediate process (such as `.aux`, `.synctex`, etc... in latex). We often don't care about keeping track of these files. Another wonderful thing about `git` is that it allows us to ignore files we don't care about!

We do this by creating a file called `.gitignore` at the initial directory of our git repository.

For example, `echo *.garbage >> .gitignore`.

Then create a file called `foo.garbage` with `touch foo.garbage`.

Now type `git status`. Notice that this file doesn't show up! It is because it recognizes that anything with an ending of `.garbage` should be ignored (this is because `*` is treated as a wildcard. Read more about regular expressions or [wild cards](http://tldp.org/LDP/GNU-Linux-Tools-Summary/html/x11655.htm) for more information on how to use them).

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

