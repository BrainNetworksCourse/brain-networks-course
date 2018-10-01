# this was adapted from nipype Vagrantfile

VAGRANTFILE_API_VERSION = "2"

# load github username from file
require_relative 'vagrant_setup.rb'
include MyVars

$script = <<SCRIPT



if [ ! -d $HOME/anaconda3 ]
then
 # install anaconda
wget https://repo.anaconda.com/archive/Anaconda3-5.2.0-Linux-x86_64.sh -O anaconda.sh
 chmod +x anaconda.sh
 ./anaconda.sh -b
 echo "export PATH=$HOME/anaconda3/bin:\\$PATH" >> .bashrc
 echo "source activate $HOME/anaconda3/envs/py3" >> .bashrc
 echo "export PATH=$HOME/anaonda3/bin:\\$PATH" >> .env
 $HOME/anaconda3/bin/conda create -n py3 python=3.6 anaconda
fi

source $HOME/anaconda3/bin/activate $HOME/anaconda3/envs/py3

sudo add-apt-repository -y ppa:webupd8team/atom # install Atom

sudo apt-get update -y

# workaround http://forums.debian.net/viewtopic.php?f=10&t=101659
sudo /usr/share/debconf/fix_db.pl
sudo apt-get install -y dictionaries-common
sudo /usr/share/debconf/fix_db.pl

sudo apt-get install -y --force-yes build-essential
sudo apt-get install -y --force-yes  python3-dev libxml2-dev libxslt1-dev zlib1g-dev
sudo apt-get install -y --force-yes git


conda install --yes -c anaconda biopython
conda install --yes -c conda-forge python-igraph
conda install --yes -c conda-forge nibabel
conda install --yes -c conda-forge nilearn
conda install --yes -c conda-forge nipy
conda install --yes -c conda-forge xlrd

pip install --upgrade https://github.com/nipy/nipype/archive/master.zip
conda install --yes -c conda-forge dipy
conda install --yes vtk


pip install python-louvain
pip install SimpleITK

wget -O- http://neuro.debian.net/lists/trusty.us-nh.full | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
sudo apt-key adv --recv-keys --keyserver hkp://pgp.mit.edu:80 0xA5D32F012649A5A9
sudo apt-get update -y

sudo apt-get install -y --force-yes r-base \
fsl-core \
fsl-atlases \
lxde \
lightdm \
chromium-browser \
wget \
tar \
unzip \
default-jre \
eog \
geany \
zip

sudo apt-get install -y atom

if [ ! -d $HOME/mcr ]
then
  echo "destinationFolder=$HOME/mcr" > mcr_options.txt
  echo "agreeToLicense=yes" >> mcr_options.txt
  echo "outputFile=/tmp/matlabinstall_log" >> mcr_options.txt
  echo "mode=silent" >> mcr_options.txt
  mkdir -p matlab_installer
  echo "downloading matlab runtime..."
  wget -nc --quiet http://www.mathworks.com/supportfiles/downloads/R2015a/deployment_files/R2015a/installers/glnxa64/MCR_R2015a_glnxa64_installer.zip -O $HOME/matlab_installer/installer.zip
  unzip $HOME/matlab_installer/installer.zip -d $HOME/matlab_installer/
  ./matlab_installer/install -inputFile mcr_options.txt
  echo "export PATH=$HOME/mcr/v85/bin/:\\$PATH" >> .bashrc
  echo "export PATH=$HOME/mcr/v85/bin/:\\$PATH" >> .env
  rm -rf matlab_installer mcr_options.txt
fi

if [ ! -d $HOME/bctpy ]
then
  git clone https://github.com/aestrivex/bctpy.git
  cd bctpy
  python setup.py install
  cd
fi

if [ ! -d $HOME/spm12 ]
then
  echo "downloading SPM..."
  wget --quiet http://www.fil.ion.ucl.ac.uk/spm/download/restricted/utopia/dev/spm12_latest_BI_Linux_R2018z.zip -O spm12.zip
  unzip spm12.zip
  echo 'alias spm="$HOME/spm12/run_spm12.sh $HOME/mcr/v85/"' >> .bashrc
  echo 'alias spm="$HOME/spm12/run_spm12.sh $HOME/mcr/v85/"' >> .env
  rm -rf spm12.zip
  cd
fi

sudo sh -c 'echo "[SeatDefaults]
user-session=LXDE
autologin-user=vagrant
autologin-user-timeout=0" >> /etc/lightdm/lightdm.conf'

echo 'export FSLDIR=/usr/share/fsl/5.0' >> .bashrc
echo ". /usr/share/fsl/5.0/etc/fslconf/fsl.sh"  >> .bashrc
echo "export FMRIDATADIR=$HOME/data" >> .bashrc

if [ ! -d $HOME/R_libs ]
then
  mkdir $HOME/R_libs
  echo "export R_LIBS_USER=$HOME/R_libs" >> .bashrc
  echo "export R_LIBS_USER=$HOME/R_libs" >> .env
fi

# get example data
if [ ! -d $HOME/data ]
then
  mkdir $HOME/data
fi

# get connectome workbench
if [ ! -d $HOME/workbench ]
then
  wget https://ftp.humanconnectome.org/workbench/workbench-linux64-v1.3.2.zip
  unzip workbench-linux64-v1.3.2.zip
fi


if [ ! -d $HOME/tetrad ]
then
 echo "setting up tetrad"
 mkdir $HOME/tetrad
 wget --quiet  http://www.phil.cmu.edu/projects/tetrad_download/download/tetrad-4.3.10-7.jar -O $HOME/tetrad/tetrad-4.3.10-7.jar
fi

# get this repo
if [ ! -d $HOME/brain-networks-course ]
then
	git clone https://github.com/#{GITHUB_USERNAME}/brain-networks-course
  cd brain-networks-course
  git remote add upstream https://github.com/BrainNetworksCourse/brain-networks-course.git
  cd
fi

# add wxpython dependencies for FSLeyes
sudo apt-get install -y --force-yes freeglut3 libsdl1.2debian
sudo apt-get install -y --force-yes libgtk2.0-dev libgtk-3-dev libwebkitgtk-dev libwebkitgtk-3.0-dev
sudo apt-get install -y --force-yes libjpeg-turbo8-dev libtiff5-dev libsdl1.2-dev libgstreamer1.0-dev libgstreamer-plugins-base1.0-dev libnotify-dev freeglut3-dev

# install FSLeyes
pip install -f https://extras.wxpython.org/wxPython4/extras/linux/gtk2/ubuntu-16.04 wxpython
pip install fsleyes


# clean up

sudo apt-get clean -y
sudo apt-get autoclean -y
sudo apt-get autoremove -y

sudo VBoxClient --display -d
sudo VBoxClient --clipboard -d
sudo VBoxClient --draganddrop -d
sudo VBoxClient --checkhostversion -d
sudo VBoxClient --seamless -d

sudo sed -i 's/allowed_users=.*$/allowed_users=anybody/' /etc/X11/Xwrapper.config

SCRIPT

Vagrant.configure(VAGRANTFILE_API_VERSION) do |config|

  config.ssh.forward_x11 = true

  config.vm.box = "ubuntu/trusty64"
  config.vm.box_version = "20180814.0.0"

  config.vm.provider :virtualbox do |vb|
      vb.customize ["modifyvm", :id, "--ioapic", "on"]
      vb.customize ["modifyvm", :id, "--memory", "5120"]
      vb.customize ["modifyvm", :id, "--cpus", "2"]
      vb.customize ["setextradata", :id, "GUI/MaxGuestResolution", "any"]
      vb.customize ["modifyvm", :id, "--natdnshostresolver1", "on"]
      vb.customize ["modifyvm", :id, "--natdnsproxy1", "on"]
      vb.customize ["modifyvm", :id, "--nictype1", "virtio"]
      vb.customize ["modifyvm", :id, "--vram", "64"]
      vb.gui = true
      vb.name = "brain-networks-course"

  end
    # uncomment following line to allow syncing to local machine
    config.vm.synced_folder ".", "/vagrant"
    config.vm.provision "shell", :privileged => false, inline: $script

end
