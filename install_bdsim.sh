#!/bin/bash

function GOTO() {
    cmd=$(sed -n "/:$1:/{:a;n;p;ba};" $0 | grep -v ':$')
    eval "$cmd"
}     

GOTO 'colors'

:start:
echo -e "${LIMEGREEN}This script will download, build and install geant4 and bdsim,"
echo -e " which are necessary for using the geant4 scattering module.\n"
echo -e "The recommended way is to set up a new environment using mamba or conda.\n${RESET}"

GOTO 'fundefs'

:main:
xcollpath=$(pwd)
cd ../
rootpath=$(pwd)
cd $xcollpath

while true; do
    read -p "$(echo -e ${CYAN}Set up conda or mamba environment? [y/n] ${RESET}\n${RED}WARNING: this will use $HOME/miniconda3 as dir for your conda setup. If you have it in a different directory, please set up the environment manually before building.${RESET})" yn
    case $yn in
        [Yy]* ) echo -e "${YELLOW}setting up environment...${RESET}"; setupEnvironment ;
                break ;;
        [Nn]* ) echo -e "${YELLOW}will not set up environment or install any packages, if geant4 or bdsim cannot be built, please do it manually${RESET}";
                break ;;
        * ) echo -e "${YELLOW}Please answer y or n ${RESET}";;
    esac
done

cd $rootpath
echo -e "${YELLOW}downloading geant4...${RESET}"
wget https://gitlab.cern.ch/geant4/geant4/-/archive/v10.4.3/geant4-v10.4.3.tar.gz
tar -xvf geant4-v10.4.3.tar.gz
mkdir -p geant4-v10.4.3/build
cd geant4-v10.4.3/build
cmake $rootpath/geant4-v10.4.3 \
    -DCMAKE_INSTALL_PREFIX=$rootpath/geant4-v10.4.3 \
    -DGEANT4_INSTALL_DATA=ON \
    -DGEANT4_INSTALL_DATADIR=$rootpath/geant4-v10.4.3 \
    -DGEANT4_USE_GDML=ON \
    -DGEANT4_USE_QT=ON \
    -DGEANT4_USE_OPENGL_X11=ON \
    -DGEANT4_USE_RAYTRACER_X11=ON \
    -DGEANT4_USE_SYSTEM_CLHEP=ON \
    -DGEANT4_USE_SYSTEM_EXPAT=ON
make -j $(nproc) 
make install
cd $rootpath

echo
echo -e "${YELLOW}downloading bdsim...${RESET}"
git clone --recursive https://github.com/bdsim-collaboration/bdsim.git

# edit bdsim/CMakeLists.txt to change CMAKE_CXX_STANDARD to 17
sed -i 's/set(CMAKE_CXX_STANDARD [0-9]\+)/set(CMAKE_CXX_STANDARD 17)/' bdsim/CMakeLists.txt
mkdir -p bdsim/build
cd $rootpath/geant4-v10.4.3/bin
source geant4.sh
cd $rootpath/bdsim/build
cmake $rootpath/bdsim \
    -DCMAKE_INSTALL_PREFIX=$rootpath/bdsim \
    -DUSE_SIXTRACKLINK=ON \
    -DBDSIM_BUILD_STATIC_LIBS=ON \
    -DGEANT4_CONFIG=$rootpath/geant4-v10.4.3/bin/geant4-config
make -j $(nproc) 
make install
unset LD_LIBRARY_PATH

echo
echo -e "${YELLOW}compiling collimasim...${RESET}"
cd $xcollpath
# check if collimasim is already downloaded
if [ ! -d "xcoll/scattering_routines/geant4/collimasim" ];
then
    echo -e "${YELLOW}collimasim not found, downloading...${RESET}"
    pushd xcoll/scattering_routines/geant4
    git clone --recurse-submodules https://gitlab.cern.ch/anabramo/collimasim.git
    if [ $? -ne 0 ]; then
        echo -e "${RED}ERROR: failed to clone collimasim repository${RESET}"
        popd
        exit 1
    fi
    popd
fi
source $rootpath/bdsim/bin/bdsim.sh     
./compile_collimasim.sh

echo -e "${LIMEGREEN}Installation complete!${RESET}"
echo
echo -e "${CYAN}Do not forget to pip install all the xsuite packages, including xcoll. To use xcoll-geant4, you will need to run the following commands in the terminal where you will be simulating:\n"
echo -e "eval \"\$($HOME/miniconda3/bin/$envexe shell.bash hook)\""
echo -e "source $HOME/miniconda3/bin/activate $envname"
echo -e "source $rootpath/geant4-v10.4.3/bin/geant4.sh"
echo -e "unset LD_LIBRARY_PATH"
echo -e "source $rootpath/bdsim/bin/bdsim.sh${RESET}"
echo
echo -e "${BOLD}${MAGENTA}HAVE FUN!${RESET}"


exit 0

:colors:
# Color codes
RED='\033[31m'
GREEN='\033[32m'
YELLOW='\033[33m'
BLUE='\033[34m'
MAGENTA='\033[35m'
CYAN='\033[36m'
WHITE='\033[37m'
BOLD='\033[1m'
RESET='\033[0m'
LIMEGREEN='\033[38;5;10m'  

GOTO 'start'      

:fundefs:
# Global error handler
handle_error() {
    echo "${RED}ERROR on line $1"
    echo "exiting...${RESET}"
    exit 1
}
trap 'handle_error $LINENO' ERR    

function setupEnvironment(){
    while true; do
        read -p "$(echo -e ${CYAN}Do you prefer mamba or conda? [mamba/conda] ${RESET})" envexe
        if [[ "$envexe" == "mamba" || "$envexe" == "conda" ]]; then
            break
        else
            echo -e "${YELLOW}Please enter 'mamba' or 'conda'${RESET}"
        fi
    done    
    ### checking if the executable is accessible
    if ! command -v $envexe &> /dev/null; then
        echo -e "${YELLOW}$envexe could not be found, checking if already installed...${RESET}"
        if [ -d "$HOME/miniconda3" ]; then
            echo -e "${LIMEGREEN}found it, initializing...${RESET}"
            eval "$($HOME/miniconda3/bin/$envexe shell.bash hook)"
        else
            mkdir -p miniconda
            cd miniconda
            wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
            bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
            eval "$($HOME/miniconda3/bin/$envexe shell.bash hook)"
            cd ..
        fi
    fi

    read -p "$(echo -e ${CYAN}Input name of your new environment: ${RESET})" envname
    $envexe create --name $envname python=3.9
    source $HOME/miniconda3/bin/activate $envname
    if [ $? -ne 0 ]; then
        echo -e "${RED}Error: Failed to activate environment $envname.${RESET}"
        exit 1
    fi
    pip install numpy scipy matplotlib pandas ipython pytest pyyaml numexpr schema tables

    while true; do
        read -p "$(echo -e ${CYAN}Do you want to download and install the other xsuite packages? ${RESET})" yn
        case $yn in
            [Yy]* ) echo -e "${YELLOW}downloading and installing in the parent directory...${RESET}"; cd ../ ; 
                    git clone https://github.com/xsuite/xobjects ; 
                    git clone https://github.com/xsuite/xdeps ; 
                    git clone https://github.com/xsuite/xpart ; 
                    git clone https://github.com/xsuite/xtrack ; 
                    git clone https://github.com/xsuite/xfields ; 
                    pip install -e xobjects xdeps xpart xtrack xfields ; 
                    cd xcoll ;
                    break ;;
            [Nn]* ) echo -e "${YELLOW}please pip install them manually in this environment.${RESET}";
                    break ;;
            * ) echo -e "${YELLOW}Please answer y or n${RESET}";;
        esac
    done

    $envexe config --add channels conda-forge
    $envexe config --set channel_priority strict
    $envexe install compilers cmake make -y
    $envexe install root=6.28.12 -y # to get root version compiled with C++17
    $envexe install clhep=2.4.5.3 -y
    $envexe install xerces-c -y
    $envexe install boost -y
    $envexe install qt-main -y
    $envexe install mesa-libgl-devel-cos7-x86_64 mesa-dri-drivers-cos7-x86_64 libselinux-cos7-x86_64 libxdamage-cos7-x86_64 libxxf86vm-cos7-x86_64 libxext-cos7-x86_64 -y
    $envexe install xorg-libxmu -y
    $envexe install flex bison -y
    $envexe install xorg-renderproto xorg-xextproto xorg-xproto -y
}

GOTO 'main'        
