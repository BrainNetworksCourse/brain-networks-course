## A brief git tutorial

This tutorial assumes that you have already installed git and have a github account.

### A mental model for version control

### Creating a repository

```git init```

### Create a file

```echo 'print("hello world")' > test.py```

### check the status

```git status```

### Add the file to the repository

```git add test.py```

### commit the file

```git commit -m"initial commit"```

### modify the file

```echo 'print("hello world!")' >| test.py```

### check the status

```git status```

### look at the difference

```git diff```

### add/commit the changed file

```git add test.py```

```git commit -m"add enthusiasm"```

### look at the log

```git log```

### revert to the initial commit

```git revert <initial commit hash>```

```git log```

## Using github

This assumes that you already have an account on github, which you should log into.

### Create a repository

go to github.com, log into your account, and create a new repository (e.g. called "git_test").
