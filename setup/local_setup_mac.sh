# this setup file adapted from the Vagrantfile for local installation

BASEDIR=$HOME/brain-networks-software

if [ ! -d $BASEDIR ]
then
    mkdir $BASEDIR
fi

conda update --yes conda
conda create --yes -n py36 python=3.6 anaconda

source activate py36

conda install --yes -c anaconda biopython
conda install --yes -c conda-forge python-igraph
conda install --yes -c conda-forge nibabel
conda install --yes -c conda-forge nilearn
conda install --yes -c conda-forge nipy
conda install --yes -c conda-forge xlrd

pip install --upgrade https://github.com/nipy/nipype/archive/master.zip
conda install --yes -c conda-forge dipy
conda install --yes vtk

pip install SimpleITK
pip install --upgrade brainnetworks

cd $BASEDIR
git clone https://github.com/aestrivex/bctpy.git
cd bctpy
python setup.py install
cd


echo "downloading SPM..."
cd $BASEDIR
wget https://www.fil.ion.ucl.ac.uk/spm/download/restricted/utopia/spm12_r7219.zip -O spm12.zip
unzip spm12.zip
rm -rf spm12.zip
cd

# get connectome workbench
echo "downloading connectome workbench"
cd $BASEDIR
wget https://ftp.humanconnectome.org/workbench/workbench-mac64-v1.3.2.zip
unzip workbench-mac64-v1.3.2.zip


echo "downloading tetrad"
cd $BASEDIR
wget --quiet  http://www.phil.cmu.edu/projects/tetrad_download/download/tetrad-4.3.10-7.jar


# install FSLeyes
cd $BASEDIR
wget https://fsl.fmrib.ox.ac.uk/fsldownloads/fsleyes/FSLeyes-latest-macos.zip
unzip FSLeyes-latest-macos.zip

echo "downloading MATLAB runtime"
echo "Installer will open automatically"
wget -nc --quiet http://ssd.mathworks.com/supportfiles/downloads/R2018b/deployment_files/R2018b/installers/maci64/MCR_R2018b_maci64_installer.dmg.zip -O $HOME/Downloads/MCR_R2018b_maci64_installer.dmg.zip
unzip $HOME/Downloads/MCR_R2018b_maci64_installer.dmg.zip -d $HOME/Downloads
open $HOME/Downloads/MCR_R2018b_maci64_installer.dmg
