---
title:  'Intro to GitHub'
author: "Spencer Lyon"
date : "2016-02-05"
---

# Git `remote`s

## What's a `remote`?

- A git repo can be on your hard drive: called local
- A `remote` is a copy of the repository on someone else's hard drive or server
- You `git push` commits from local to remote
- `git pull` commits from remote to local

## Collaboration

- Working with remotes enables many [workflows](https://git-scm.com/book/en/v2/Distributed-Git-Distributed-Workflows)
- One common workflow (image taken from [here](https://git-scm.com/book/en/v2/Distributed-Git-Distributed-Workflows))

![Centralized workflow](images/git_centralized_workflow.png)

# [GitHub](https://github.com)

## Facts

- GitHub is a common place to have a `remote`
    - [Over 2.2 million](http://githut.info) active repositories
    - Unlimited free public repositories
    - 5 free private repositories for academics

## GitHub Collaboration

- Permission management: only certain users can `push`
- Forking: anyone can create a copy of any (visible) repository under their account
- Pull requests: system for project maintainers to review proposed changes before accepting them
- Automated testing (CI) and coverage

# Example

## quantecon_nyu_2016

- Add this lecture to course repo
- Steps:
    1. Login to github (online)
    2. Clone repository
    ```sh
    git clone https://github.com/jstac/quantecon_nyu_2016.git
    ```
    3. Fork repository (online)
    4. Add `remote` that points to fork
    ```sh
    git remote add fork https://github.com/spencerlyon2/quantecon_nyu_2016.git
    ```

## Steps (continued)

- More steps
    5. `add` + `commit` files
    ```
    git add files
    git commit -m "Adding github intro slides"
    ```
    6. `push` to fork
    ```
    git push fork master
    ```
    7. Open pull request (online)

# Exercise

## Pushing homework

- We will let you practice hw submission process
- Steps
    - Create account on github
    - Clone [homework repository](https://github.com/jstac/quantecon_nyu_2016_homework)
    - Fork homework repo
    - Add your fork as remote
    - Create file `firstname_lastname` in folder `hw_github_intro`
    - Leave short message in file
    - Add and commit the file
    - Push to fork
    - Submit pull request


# Extras

## SSH keys

- Two modes of authenticating to github: https, SSH
- https will require you to enter password on every push
- Adding ssh-key removes this requirement
- Follow these steps:
    - `ssh-keygen`: follow prompts
    - Copy the output of `cat ~/.ssh/id_rsa.pub` (use shift-c)
        - NOTE: path `~/.ssh/id_rsa.pub` might be different depending on your answers to previous step
        - Important part is to get the `.pub` version
    - Go to github account [settings](https://github.com/settings/profile)
    - Click "ssh keys" link on left, then "New ssh key"
    - Give it a title (name) and paste clipboard contents in key
    - Click big green "add SSH key button"
