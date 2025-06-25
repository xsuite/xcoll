#!/bin/bash

# python_ver=3.9
python_ver=3.13  # does not work
root_ver=6.34.4 #"6.34.4=py313h7eb37b6_0"
# clhep_ver=2.4.5.3
clhep_ver=2.4.7.1   # does not work
# geant4_ver=10.4.3
geant4_ver=11.3.2


function GOTO() {
    cmd=$(sed -n "/:$1:/{:a;n;p;ba};" $0 | grep -v ':$')
    eval "$cmd"
}

GOTO 'colors'

:start:
echo -e "${LIMEGREEN}This script will download, build and install geant4 and bdsim,"
echo -e "which are necessary for using the geant4 scattering module.\n"
echo -e "The recommended way is to set up a new environment using mamba or conda.\n${RESET}"

GOTO 'fundefs'

:main:
xcoll_path=$(pwd)
if [ ! -f $xcoll_path/install_bdsim.sh ]
then
    echo -e "${RED}This script should be run from the xcoll directory.${RESET}"
    exit 1
fi
geant4_path=${xcoll_path}/xcoll/lib/
mkdir $geant4_path
cd ${xcoll_path}/..

while true
do
    read -p "$(echo -e ${CYAN}Set up conda or mamba environment? [y/n] ${RESET})" yn
    case $yn in
        [Yy]* ) echo -e "${YELLOW}setting up environment...${RESET}"; setupEnvironment ;
                break ;;
        [Nn]* ) echo -e "${YELLOW}will not set up environment or install any packages, if geant4 or bdsim cannot be built, please do it manually${RESET}";
                break ;;
        * ) echo -e "${YELLOW}Please answer y or n ${RESET}";;
    esac
done

cd $geant4_path
echo -e "${YELLOW}downloading geant4 v${geant4_ver}...${RESET}"
wget https://gitlab.cern.ch/geant4/geant4/-/archive/v${geant4_ver}/geant4-v${geant4_ver}.tar.gz
tar -xvf geant4-v${geant4_ver}.tar.gz
rm geant4-v${geant4_ver}.tar.gz
mkdir -p geant4-v${geant4_ver}/build
cd geant4-v${geant4_ver}/build
cmake ${geant4_path}/geant4-v${geant4_ver} \
    -DCMAKE_INSTALL_PREFIX=${geant4_path}/geant4-v${geant4_ver} \
    -DGEANT4_BUILD_MULTITHREADED=OFF \
    -DGEANT4_INSTALL_DATA=ON \
    -DGEANT4_INSTALL_DATADIR=${geant4_path}/geant4-v${geant4_ver} \
    -DGEANT4_USE_GDML=ON \
    -DGEANT4_USE_OPENGL_X11=ON \
    -DGEANT4_USE_QT=ON \
    -DGEANT4_USE_SYSTEM_CLHEP=ON \
    -DGEANT4_USE_SYSTEM_ZLIB=OFF \
    -DGEANT4_USE_RAYTRACER_X11=ON \
    -DGEANT4_USE_SYSTEM_EXPAT=ON \
    -DGEANT4_USE_XM=ON
make -j $(nproc)
make install
cd $geant4_path

echo
echo -e "${YELLOW}downloading latest version of bdsim...${RESET}"
git clone --recursive https://github.com/bdsim-collaboration/bdsim.git

# edit bdsim/CMakeLists.txt to change CMAKE_CXX_STANDARD to 17
# sed -i 's/set(CMAKE_CXX_STANDARD [0-9]\+)/set(CMAKE_CXX_STANDARD 17)/' bdsim/CMakeLists.txt
mkdir -p bdsim/build
cd ${geant4_path}/geant4-v${geant4_ver}/bin
source geant4.sh
cd ${geant4_path}/bdsim/build
cmake ${geant4_path}/bdsim \
    -DCMAKE_INSTALL_PREFIX=${geant4_path}/bdsim \
    -DUSE_SIXTRACKLINK=ON \
    -DBDSIM_BUILD_STATIC_LIBS=ON \
    -DGEANT4_CONFIG=$(which geant4-config)
    # -DGEANT4_CONFIG=${geant4_path}/geant4-v${geant4_ver}/bin/geant4-config
make -j $(nproc) 
make install
unset LD_LIBRARY_PATH

echo
echo -e "${YELLOW}compiling collimasim...${RESET}"
cd $xcoll_path
source $geant4_path/bdsim/bin/bdsim.sh
./compile_collimasim.sh

echo -e "${LIMEGREEN}Installation complete!${RESET}"
echo
echo -e "${CYAN}Do not forget to pip install all the xsuite packages, including xcoll. To use xcoll-geant4, you will need to run the following commands in the terminal where you will be simulating:\n"
echo -e 'eval "$(conda shell.bash hook)"'
echo -e "conda activate $envname"
echo -e "cd ${geant4_path}/geant4-v${geant4_ver}/bin"
echo -e "source geant4.sh"
echo -e "cd -"
echo -e "unset LD_LIBRARY_PATH"
echo -e "source ${geant4_path}/bdsim/bin/bdsim.sh${RESET}"
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
    ### checking if the executable is accessible
    if ! command -v conda &> /dev/null
    then
        echo -e "${RED}conda could not be found! Please install it first.${RESET}"
        exit 1
    fi
    ### activate conda scripts (cannot init mamba this way unfortunately)
    eval "$(conda shell.bash hook)"

    while true; do
        read -p "$(echo -e ${CYAN}Create new environment? [y/n] ${RESET})" yn
        case $yn in
            [Yy]* ) read -p "$(echo -e ${CYAN}Input name of your new environment: ${RESET})" envname
                    conda create --name $envname python=${python_ver} -y
                    conda activate $envname
                    if [ $? -ne 0 ]; then
                        echo -e "${RED}Error: Failed to activate environment $envname.${RESET}"
                        exit 1
                    fi
                    echo -e "${YELLOW}Installing dependencies...${RESET}";
                    break ;;
            [Nn]* ) echo -e "${YELLOW}Installing dependencies in the current environment...${RESET}";
                    break ;;
            * ) echo -e "${YELLOW}Please answer y or n ${RESET}";;
        esac
    done

    python -m pip install numpy scipy matplotlib pandas ruamel.yaml

    while true; do
        read -p "$(echo -e ${CYAN}Do you want to download and install the other xsuite packages? ${RESET})" yn
        case $yn in
            [Yy]* ) echo -e "${YELLOW}downloading and installing in the parent directory...${RESET}"; cd ../ ; 
                    git clone https://github.com/xsuite/xobjects ;
                    git clone https://github.com/xsuite/xdeps ;
                    git clone https://github.com/xsuite/xpart ;
                    git clone https://github.com/xsuite/xtrack ;
                    git clone https://github.com/xsuite/xfields ;
                    python -m pip install -e ./xobjects -e ./xdeps -e ./xtrack -e ./xpart -e ./xfields ;
                    cd xcoll ;
                    break ;;
            [Nn]* ) echo -e "${YELLOW}please pip install them manually in this environment.${RESET}";
                    break ;;
            * ) echo -e "${YELLOW}Please answer y or n${RESET}";;
        esac
    done

    conda config --add channels conda-forge
    conda config --set channel_priority strict
    conda install compilers cmake make -y
    conda install root=${root_ver} -y # to get root version compiled with C++17
    conda install clhep=${clhep_ver} -y
    conda install xerces-c -y
    conda install boost -y
    conda install qt-main -y
    conda install mesa-libgl-devel-cos7-x86_64 mesa-dri-drivers-cos7-x86_64 libselinux-cos7-x86_64 libxdamage-cos7-x86_64 libxxf86vm-cos7-x86_64 libxext-cos7-x86_64 -y
    conda install xorg-libxmu -y
    conda install flex bison -y
    conda install xorg-renderproto xorg-xextproto xorg-xproto -y
}

GOTO 'main'
