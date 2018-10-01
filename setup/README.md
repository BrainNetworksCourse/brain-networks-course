### Environment setup for Brain Networks course

For this course we will use a virtual machine that will allow everyone to have
an identical environment regardless of their operating system. You are of course
free to install all of the dependencies on your own system, but they are not
guaranteed to work (and some will certainly not work if you are on Windows).

To use the VM, first set up the following software on your computer:

- Install [VirtualBox](https://www.virtualbox.org/)

- Install [Vagrant](https://www.vagrantup.com/)
    - Windows users will likely need to [turn off Hyper-V](https://ugetfix.com/ask/how-to-disable-hyper-v-in-windows-10/)

- Install [git](https://git-scm.com/) if you don't already have it.

### Creating your own fork of the class repository.

Go to github.com and log into your account.  Then go to the class repository:

https://github.com/poldrack/brain-networks-course

Click the Fork button to create a fork of the repository in your github account.  This will allow you to make changes and submit them for inclusion in the main repository, and is the standard approach used for collaborative software development on github.


Clone the repository onto your computer:

```git clone https://github.com/<your github username>/brain-networks-course.git```

inserting your github username in the appropriate place in the URL.

### Setting up vagrant

Go into the repo directory and open the file called "vagrant_setup.rb" in your favorite text editor.  Change the following line:

```GITHUB_USERNAME = "poldrack"```

by replacing poldrack with your own github username.  Be sure to save the file and then commit it to your fork:

```
git add vagrant_setup.rb
git commit -m"changing github username"
git push origin master
```

### Provision the virtual machine

Run the following command to provision the
virtual machine (this will take a little while the first time you do it, and involves lots of downloading so be sure to do it from a good network):

  ```vagrant up```

Once the installation is done, then restart the virtual machine; to do this, go to the VirtualBox window for the VM (which should just show a login window), close the VM window, and choose "Power off the machine".  Then run ```vagrant up``` again to restart it.  At this point, a GUI window should appear for the VM, which we will use for our class exercises.

## Updating the VM

On occasion we will update the VM with new packages.  To pull the latest changes into your fork, you will first need to add the main repository as a remote to your repo:

```
git remote add upstream https://github.com/poldrack/brain-networks-course.git
```

After doing this once, you can then use the following command within the repo directory to pull the latest changes:

```
git pull upstream master
```

You may need to commit any changes you have made to other files in your repo, or use ```git checkout <filename>``` to pull a clean version if you don't need to keep the changes.

Once you have done this it should pull down the new Vagrantfile, and you can then run ```vagrant provision``` in the main repo directory which will rerun the provisioning and install any new software.  

### Some tips for using the VM

- I find the default screensaver to use up a lot of processing power, so I generally change it to a blank screen (under Preferences->Screensaver).
- You can add an ssh key for your VM so that you won't have to enter your github credentials each time you push or pull.  See [here](https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/) for info on how to do this.  


## Alternative: Local installation

Some people have had trouble installing the virtual machine on their systems, in which case it may be preferable to install the software locally.  Here are the steps for installing on a Mac.  For Windows users, I would suggest sticking with the VM or installing the Windows Linux Subsystem (which requires Windows 10 Professional):

1.  Download and install the [Anaconda](https://www.anaconda.com/download/) Python distribution.
2. Go to your local clone of this reposistory, and pull the latest version: ```git pull upstream master```
3. Execute the setup script: ```sh setup/local_setup_mac.sh```

This will install a new python virtual environment called py36 - whenever you want to use this in a new session you will first need to activate the virtual environment using ```source activate py36```

All downloaded software will be installed in ~/brain-networks-software
