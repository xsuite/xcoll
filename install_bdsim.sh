echo "This script will download, build and install geant4 and bdsim,"
echo " which are necessary for using the geant4 scattering module."
echo
echo "The recommended way is to set up a new environment using mamba or conda."
echo

function setupEnvironment(){
    read -p "Do you prefer mamba or conda? [mamba/conda]" envexe
    if command -v $envexe &> /dev/null
    then
        echo $envexe" could not be found, checking if already installed..."
        if [ -d "~/miniconda3" ]; then
            echo "found it, initializing..."
            eval "$(~/miniconda3/bin/"$envexe" shell.bash hook)"
        else
            mkdir miniconda
            cd miniconda
            wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
            bash Miniconda3-latest-Linux-x86_64.sh
            eval "$(~/miniconda3/bin/"$envexe" shell.bash hook)"
        fi
    fi
    read -p "Input name of your new environment: " envname
    $envexe create --name $envname python=3.9
    $envexe activate $envname
    pip install numpy scipy matplotlib pandas ipython pytest pyyaml numexpr schema tables

    read -p "Do you want to download and install the other xsuite packages?" yn
    case $yn in
        [Yy]* ) echo "downloading and installing in the parent directory..."; cd ../ ; git clone https://github.com/xsuite/xobjects ; git clone https://github.com/xsuite/xdeps ; git clone https://github.com/xsuite/xpart ; git clone https://github.com/xsuite/xtrack ; git clone https://github.com/xsuite/xfields ; pip install -e xobjects xdeps xpart xtrack xfields ; cd xcoll ;;
        [Nn]* ) echo "please pip install them manually in this environment."; ;;
        * ) echo "Please answer y or n";;
    esac

    $envexe config --add channels conda-forge
    $envexe config --set channel_priority strict
    $envexe install compilers cmake make
    $envexe install root=6.28.12  # to get root version compiled with C++17
    $envexe install clhep=2.4.5.3
    $envexe install xerces-c
    $envexe install boost
    $envexe install qt-main
    $envexe install mesa-libgl-devel-cos7-x86_64 mesa-dri-drivers-cos7-x86_64 libselinux-cos7-x86_64 libxdamage-cos7-x86_64 libxxf86vm-cos7-x86_64 libxext-cos7-x86_64
    $envexe install xorg-libxmu
    $envexe install flex bison   
}

xcollpath=$(pwd)
cd ../
rootpath=$(pwd)
cd $xcollpath

read -p "Set up conda or mamba environment? [y/n]" yn
case $yn in
    [Yy]* ) echo "setting up environment..."; setupEnvironment ;;
    [Nn]* ) echo "will not set up environment or install any packages, if geant4 or bdsim cannot be built, please do it manually"; ;;
    * ) echo "Please answer y or n";;
esac

cd $rootpath
echo "downloading geant4..."
wget https://gitlab.cern.ch/geant4/geant4/-/archive/v10.4.3/geant4-v10.4.3.tar.gz
tar -xvf geant4-v10.4.3.tar.gz
mkdir -p geant4-v10.4.3/build
cd geant4-v10.4.3/build
cmake $rootpath/geant4-v10.4.3 \
    -DCMAKE_INSTALL_PREFIX=$rootpath/geant4-v10.4.3 \
    -DGEANT4_INSTALL_DATA=ON \
    -DGEANT4_INSTALL_DATADIR=$rootpath/geant4-v10.4.3 \
    -DUSE_GDML=ON \
    -DUSE_QT=ON \
    -DUSE_OPENGL_X11=ON \
    -DUSE_RAYTRACER_X11=ON \
    -DUSE_SYSTEM_CLHEP=ON \
    -DUSE_SYSTEM_EXPAT=ON
make -j $(nproc) 
make install
cd $rootpath

echo
echo "downloading bdsim..."
git clone --recursive https://bitbucket.org/jairhul/bdsim

# edit bdsim/CMakeLists.txt to change CMAKE_CXX_STANDARD to 17
sed -i 's/set(CMAKE_CXX_STANDARD [0-9]\+)/set(CMAKE_CXX_STANDARD 17)/' bdsim/CMakeLists.txt
mkdir -p bdsim/build
source $rootpath/geant4-v10.4.3/bin/geant4.sh
cd bdsim/build
cmake $rootpath/bdsim \
    -DCMAKE_INSTALL_PREFIX=$rootpath/bdsim \
    -DUSE_SIXTRACKLINK=ON \
    -DBDSIM_BUILD_STATIC_LIBS=ON \
    -DGEANT4_CONFIG=$rootpath/geant4-v10.4.3/bin/geant4-config \
make -j $(nproc) 
make install
unset LD_LIBRARY_PATH

echo
echo "compiling collimasim..."
cd $xcollpath
# check if collimasim is already downloaded
if [ ! -d "xcoll/scattering_routines/geant4/collimasim" ];
then
    echo "collimasim not found, downloading..."
    pushd xcoll/scattering_routines/geant4
    git clone --recurse-submodules https://gitlab.cern.ch/anabramo/collimasim.git
    if [ $? -ne 0 ]; then
        echo "ERROR: failued to clone collimasim repository"
        popd
        exit 1
    fi
    popd
fi     
./compile_collimasim.sh
