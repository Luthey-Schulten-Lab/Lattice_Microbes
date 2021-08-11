#! /bin/bash -login
###############################################################################
# install.sh <option> <directory> <optional: cores>
###############################################################################
# Author:         Joseph R. Peterson
# Date Created:   09/17/2012
# Date Modified:  02/05/2012
###############################################################################



#############
# Functions #
#############

# USAGE Function #
function usage {
  echo "Usage:"
  echo "------"
  echo " install.sh <architecture> <directory> <optional: device> <optional: cores>"
  echo
  echo " Attempts to install Lattice Microbes in the user's home directory"
  echo " using the required 'directory' option as a workspace."
  echo " The 'architecture' argument may be either an installation architecture."
  echo " The 'device' option allows you to specify what builds should be made."
  echo " Default for device is 'CPUCUDA' which builds one version of Lattice "
  echo " Microbes for the CPU and one for the GPU.  Other options are 'CPU' and"
  echo " 'CUDA' which build only for the CPU and GPU respectively."
  echo " The 'cores' argument is optional and can be used to compile with more"
  echo " than one processor.  This option defaults to 1 core."
  echo " "
  echo " The '-h' or '--help' option displays this message."
  echo 
  echo "Valid architectures:"
  echo "--------------------"
  echo " osx"
  echo " keeneland"
  echo " titan"
  echo 

  exit
}

# ERROR Function #
# Takes three parameters:
#   1) Exit status
#   2) error string
#   3) log file
function checkError {
  if [[ $1 != 0 ]] ; then
    echo "ERROR: Failed with error: '$2'.  Bailing out." | tee -a $3 
    exit
  fi
}

# Download Function #
# Takes three parameters:
#   1) URL
#   2) log file
#   3) architecture
function webDownload {
  if [[ $3 == "osx" ]] ; then 
    curl -O -L $1 >> $2 2>&1
    checkError $? "Error downloading: $1" $2
  else 
    wget $1 -a $2
    checkError $? "Error downloading: $1" $2
  fi
}

# Print Functions #
# Takes three parameters:
#   1) operation
#   2) name
#   3) log file
function printOperation {
  str="$1 $2"
  dots="-"
  len=`expr ${#1} + ${#2}`
  for (( i=1; i<=$len; i++ ))
  do
    dots="$dots-"
  done 
  
  echo $str | tee -a $3
  echo $dots >> $3
}




####################################
# Parsing input and Error Checking #
####################################
# Make sure at least one argument was specified
if [[ $# -lt 1  ]] ; then
  usage
fi

# User asked for help prompt
if [[ $1 == "-h" ]] ; then
  usage
fi
if [[ $1 == "--help" ]] ; then
  usage
fi

# If not the help prompt, make sure the user specified a directory as well
if [[ $# -lt 2 ]] ; then
  usage
fi

# check which machine
architecture="unknown"
if [[ $1 == "titan" ]] ; then
  architecture="Titan"
elif [[ $1 == "keeneland" ]] ; then
  architecture="Keeneland"
elif [[ $1 == "osx" ]] ; then
  architecture="osx"
  # check for correct version and os
  os=`uname`
  if [[ $os != "Darwin"  ]] ; then
    echo "WARNING: Expected operating system 'Darwin' for architecture 'osx', recieved: $os. Trying anyway."
  fi
  osVersion=`uname -r`
  if [[ $osVersion  != "11.4.2" ]] ; then
    echo "WARNING: The 'osx' option has only been tested on version 11.4.0, recieved: $osVersion. Trying anyway."
  fi
else
  usage
fi

# create the installation directory
mkdir -p $2
dir=$2


forCPU=1
forGPU=1
if [[ $# -gt 2 ]] ; then
  if [[ $3 == "CPU" ]] ; then
    forCPU=1
    forGPU=0
  elif [[ $3 == "CUDA" ]] ; then
    forGPU=1
    forCPU=0
  elif [[ $3 == "CPUCUDA" ]] ; then
    forGPU=1
    forCPU=1
  fi
fi

cores=1
if [[ $# -gt 3 ]] ; then
  cores=$4
fi

# get the root of lm-x.y
lmdir=`pwd`

##############
# Installing #
##############
# Print to the user to let them know what is going on
echo
echo "-------------------------------------------------------------------"
echo "Installing LM with working $dir for the '$architecture' computer system"
echo "-------------------------------------------------------------------"
echo


cd $dir
export dirr=`pwd`
# set up log file
if [[ -f $lmdir/LMinstall.log ]] ; then
  rm $lmdir/LMinstall.log
fi
touch $lmdir/LMinstall.log
export log="$lmdir/LMinstall.log"


# build some working directories
mkdir -p work
mkdir -p build
export WORK=`pwd`/work
export BUILD=`pwd`/build

counter=1

# ProtoBuffers #
echo "$counter) Installing Google Protocol Buffers" | tee -a $log
echo "-------------------------------------" | tee -a $log
cd $WORK
if [[ ! -f protobuf-2.4.1.tar.gz ]] ; then
  printOperation "Downloading" "Protocol Buffers" $log
  webDownload http://protobuf.googlecode.com/files/protobuf-2.4.1.tar.gz $log $architecture
fi
echo >> $log
printOperation "Extracting" "Protocol Buffers" $log
tar -xvf protobuf-2.4.1.tar.gz >> $log 2>&1
checkError $? "Failed to extract protobuf-2.4.1.tar.gz" $log
cd protobuf-2.4.1 >> $log
mkdir -p $BUILD/protobuf 
echo >> $log
printOperation "Configuring" "Protocol Buffers" $log
./configure --prefix=$BUILD/protobuf >> $log 2>&1
checkError $? "Failed to configure protobuf-2.4.1" $log
echo >> $log
printOperation "Making" "Protocol Buffers" $log
make -j$cores >> $log 2>&1
checkError $? "Failed to make protobuf-2.4.1" $log
make install >> $log 2>&1
checkError $? "Failed to install protobuf-2.4.1" $log
counter=`expr $counter + 1`


# SBML #
echo | tee -a $log
echo "$counter) Installing SBML" | tee -a $log
echo "------------------" | tee -a $log
cd $BUILD
mkdir -p sbml
echo >> $log 
cd sbml
if [[ $architecture == "osx" ]] ; then
  # check if the right autotools is installed
  autoconfVer=`autoconf --version | head -1 | awk '{print $NF}'`
  if [ $(echo "$autoconfVer <  2.62" | bc) -ne 0 ] ; then
    printOperation "Installing" "m4 1.4.13" $log
    if [[ ! -f m4-1.4.13.tar.gz ]] ; then
      webDownload http://mirrors.kernel.org/gnu/m4/m4-1.4.13.tar.gz $log $architecture 
      checkError $? "Failed to download m4-1.4.13.tar.gz" $log
    fi
    tar -xvf m4-1.4.13.tar.gz >> $log 2>&1
    checkError $? "Failed to extract m4-1.4.13.tar.gz" $log
    mkdir -p m4
    cd m4
    M4dir=`pwd`
    cd ../m4-1.4.13
    ./configure --prefix=$M4dir >> $log 2>&1
    checkError $? "Faled to configure m4" $log
    make >> $log 2>&1
    checkError $? "Failed to make m4" $log
    make install >> $log 2>&1
    checkError $? "Failed to install m4" $log
    cd ../

    printOperation "Installing" "autoconf 2.65" $log
    if [[ ! -f autoconf-2.65.tar.gz ]] ; then
      webDownload http://mirrors.kernel.org/gnu/autoconf/autoconf-2.65.tar.gz $log $architecture 
      checkError $? "Failed to download autoconf-2.65.tar.gz" $log
    fi
    tar -xvf autoconf-2.65.tar.gz >> $log 2>&1
    checkError $? "Failed to extract autoconf-2.65.tar.gz" $log
    mkdir -p autoconf
    cd autoconf
    ACdir=`pwd`
    cd ../autoconf-2.65
    ./configure --prefix=$ACdir >> $log 2>&1
    checkError $? "Faled to configure autoconf" $log
    make  >> $log 2>&1
    checkError $? "Failed to make autoconf" $log
    make install >> $log 2>&1
    checkError $? "Failed to install autoconf" $log
    cd ../
    PATH=$ACdir/bin:$M4dir/bin:$PATH

    printOperation "Installing" "automake 1.11" $log
    if [[ ! -f automake-1.11.tar.gz ]] ; then
      webDownload http://mirrors.kernel.org/gnu/automake/automake-1.11.tar.gz $log $architecture
      checkError $? "Failed to download automake-1.11.tar.gz" $log
    fi
    tar -xvf automake-1.11.tar.gz >> $log 2>&1
    checkError $? "Failed to extract automake-1.11.tar.gz" $log
    mkdir -p automake
    cd automake
    AMdir=`pwd`
    cd ../automake-1.11
    ./configure --prefix=$AMdir >> $log 2>&1
    checkError $? "Faled to configure automake" $log
    make >> $log 2>&1
    checkError $? "Failed to make automake" $log
    make install >> $log 2>&1
    checkError $? "Failed to install automake" $log
    cd ../
    PATH=$AMdir/bin:$PATH

    printOperation "Installing" "libtool 2.2.6b" $log
    if [[ ! -f libtool-2.2.6b.tar.gz ]] ; then
      webDownload http://mirrors.kernel.org/gnu/libtool/libtool-2.2.6b.tar.gz $log $architecture
      checkError $? "Failed to download libtool-2.2.6b.tar.gz" $log
    fi
    tar -xvf libtool-2.2.6b.tar.gz >> $log 2>&1
    checkError $? "Failed to extract libtool-2.2.6b.tar.gz" $log
    mkdir -p libtool
    cd libtool
    LTdir=`pwd`
    cd ../libtool-2.2.6b
    ./configure --prefix=$LTdir >> $log 2>&1
    checkError $? "Faled to configure libtool" $log
    make >> $log 2>&1
    checkError $? "Failed to make libtool" $log
    make install >> $log 2>&1
    checkError $? "Failed to install libtool" $log
    cd ../
    PATH=$LTdir/bin:$PATH


    if [[ ! -f libsbml-5.6.0-src.tar.gz ]] ; then
      printOperation "Downloading" "SBML" $log
      webDownload http://downloads.sourceforge.net/project/sbml/libsbml/5.6.0/stable/libsbml-5.6.0-src.tar.gz $log $architecture
      checkError $? "Failed to download libsbml-5.6.0-src.tar.gz" $log
    fi
    tar -xvf libsbml-5.6.0-src.tar.gz >> $log 2>&1

    printOperation "Installing" "SBML" $log
    SBML=`pwd`
    cd libsbml-5.6.0
    printOperation "Configuring" "SBML" $log
    ./configure --prefix=$SBML >> $log 2>&1
    checkError $? "Faled to configure libsbml" $log
    printOperation "Making" "SBML" $log
    make >> $log 2>&1
    checkError $? "Failed to make libsbml" $log
    printOperation "Installing" "SBML" $log
    make install >> $log 2>&1
    checkError $? "Failed to install libsbml" $log

  else 
    if [[ ! -f libsbml-5.6.0-src.tar.gz ]] ; then
      printOperation "Downloading" "SBML" $log
      webDownload http://downloads.sourceforge.net/project/sbml/libsbml/5.6.0/stable/libsbml-5.6.0-src.tar.gz $log $architecture
      checkError $? "Failed to download libsbml-5.6.0-src.tar.gz" $log
    fi
    tar -xvf libsbml-5.6.0-src.tar.gz >> $log 2>&1

    printOperation "Installing" "SBML" $log
    SBML=`pwd`
    cd libsbml-5.6.0
    ./configure --prefix=$SBML >> $log 2>&1
    checkError $? "Faled to configure libsbml" $log
    make >> $log 2>&1
    checkError $? "Failed to make libsbml" $log
    make install >> $log 2>&1
    checkError $? "Failed to install libsbml" $log
  fi
else 
  if [[ ! -f libSBML-5.6.0-Linux-x64.rpm ]] ; then
    printOperation "Downloading" "SBML" $log
    webDownload http://downloads.sourceforge.net/project/sbml/libsbml/5.6.0/stable/Linux/64-bit/libSBML-5.6.0-Linux-x64.rpm $log $architecture
    checkError $? "Failed to download libSBML-5.6.0-Linux-x64.rpm" $log
  fi
  echo >> $log 
  printOperation "Installing" "SBML" $log
  rpm2cpio libSBML-5.6.0-Linux-x64.rpm | cpio -idmv >>  $log 2>&1
  checkError $? "Failed to install libSBML-5.6.0-Linux-x64.rpm" $log
fi
counter=`expr $counter + 1`


# SZip #
echo | tee -a $log
echo "$counter) Installing SZip" | tee -a $log
echo "------------------" | tee -a $log
cd $WORK
if [[ ! -f szip-2.1.tar.gz ]] ; then
  printOperation "Downloading" "SZip" $log
  webDownload http://www.hdfgroup.org/ftp/lib-external/szip/2.1/src/szip-2.1.tar.gz $log $architecture
  checkError $? "Failed to download szip-2.1.tar.gz" $log
fi
echo >> $log
printOperation "Extracting" "SZip" $log
tar -xvf szip-2.1.tar.gz >> $log 2>&1
checkError $? "Failed to extract szip-2.1.tar.gz" $log
cd szip-2.1
mkdir -p $BUILD/szip
echo >> $log
printOperation "Configuring" "SZip" $log
./configure --prefix=$BUILD/szip >> $log  2>&1
checkError $? "Failed to configure szip-2.1" $log
echo >> $log
printOperation "Making" "SZip" $log
make -j$cores >> $log 2>&1
checkError $? "Failed to make szip-2.1" $log
printOperation "Installing" "SZip" $log
make install >> $log 2>&1
checkError $? "Failed to install szip-2.1" $log
counter=`expr $counter + 1`


# HDF5 #
echo | tee -a $log
echo "$counter) Installing HDF5" | tee -a $log
echo "------------------" | tee -a $log
cd $WORK
if [[ ! -f hdf5-1.8.9.tar.gz ]] ; then
  printOperation "Downloading" "HDF5" $log
  webDownload http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8.9/src/hdf5-1.8.9.tar.gz $log $architecture
  checkError $? "Failed to download hdf5-1.8.9.tar.gz" $log
fi
echo >> $log
printOperation "Extracting" "HDF5" $log
tar -xvf hdf5-1.8.9.tar.gz >> $log 2>&1
checkError $? "Failed to extract hdf5-1.8.9.tar.gz" $log
cd hdf5-1.8.9
mkdir -p $BUILD/hdf5
echo >> $log
printOperation "Configuring" "HDF5" $log
./configure  --prefix=$BUILD/hdf5 --with-szlib=$BUILD/szip >> $log
checkError $? "Failed to configure hdf5" $log
echo >> $log
printOperation "Making" "HDF5" $log
make -j$cores >> $log 2>&1
checkError $? "Failed to make hdf5" $log
printOperation "Installing" "HDF5" $log
make install >> $log 2>&1
checkError $? "Failed to install hdf5" $log
counter=`expr $counter + 1`


# VMD #
if [[ ! $architecture == "osx" ]] ; then
  echo | tee -a $log
  echo "$counter) Installing VMD" | tee -a $log
  echo "-----------------" | tee -a $log
  cd $BUILD
  if [[ ! -f vmd-1.9.1.bin.LINUXAMD64.opengl.tar.gz ]] ; then
    printOperation "Downloading" "VMD" $log
    webDownload http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.1/files/final/vmd-1.9.1.bin.LINUXAMD64.opengl.tar.gz $log $architecture
    checkError $? "Failed to download vmd-1.9.1.bin.LINUXAMD64.opengl.tar.gz" $log
  fi
  echo >> $log
  printOperation "Extracting" "VMD" $log
  tar -xvf vmd-1.9.1.bin.LINUXAMD64.opengl.tar.gz >> $log 2>&1
  checkError $? "Failed to extract vmd-1.9.1.bin.LINUXAMD64.opengl.tar.gz" $log
  counter=`expr $counter + 1`
fi

# Swig #
if [[ $architecture == "osx" ]] || [[ $architecture == "Titan" ]] ; then
  echo | tee -a $log
  echo "$counter) Installing SWIG" | tee -a $log
  echo "------------------" | tee -a $log
  cd $BUILD
  printOperation "Installing" "PCRE" $log
  rm -rf pcre-8.21
  if [[ ! -f pcre-8.21.tar.gz ]] ; then
    webDownload ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.21.tar.gz $log $architecture
    checkError $? "Failed to download pcre-8.21.tar.gz" $log
  fi
  tar -xvf pcre-8.21.tar.gz >> $log 2>&1
  checkError $? "Failed to extract pcre-8.21.tar.gz" $log
  mkdir -p pcre
  cd pcre
  PCdir=`pwd`
  cd ../pcre-8.21
  make clean >> $log 2>&1
  ./configure --prefix=$PCdir >> $log 2>&1
  checkError $? "Failed to configure PCRE" $log
  make >> $log 2>&1
  checkError $? "Failed to make PCRE" $log
  make install >> $log 2>&1
  checkError $? "Failed to install PCRE" $log
  cd ../

  if [[ ! -f swig-2.0.9.tar.gz ]] ; then
    printOperation "Downloading" "SWIG" $log
    webDownload http://downloads.sourceforge.net/project/swig/swig/swig-2.0.9/swig-2.0.9.tar.gz $log $architecture
    checkError $? "Failed to download swig-2.0.9.tar.gz" $log
  fi
  printOperation "Extracting" "SWIG" $log
  tar -xvf swig-2.0.9.tar.gz >> $log 2>&1
  mkdir -p swig
  cd swig 
  SWDir=`pwd`
  cd ../swig-2.0.9
  printOperation "Configuring" "SWIG" $log
  ./configure --prefix=$SWDir --with-python=/usr/bin/python2.6 --with-pcre-prefix=$PCdir>> $log 2>&1
  checkError $? "Failed to configure swig" $log
  printOperation "Making" "SWIG" $log
  make -j$cores >> $log 2>&1
  checkError $? "Failed to make swig" $log
  printOperation "Installing" "SWIG" $log
  make install >> $log 2>&1
  checkError $? "Failed to install swig" $log
  counter=`expr $counter + 1`
fi


# OpenMPI #
cd $BUILD
if [ "a" == "b" ] ; then
# NOTE This titan build should be replaced later with the on machine openMPI when it becomes available
if [ $architecture == "osx" -o $architecture == "Titan" ]; then
  echo | tee -a $log
  echo "$counter) Installing OpenMPI" | tee -a $log
  echo "---------------------" | tee -a $log
  if [[ ! -f openmpi-1.6.3.tar.gz ]] ; then
    printOperation "Downloading" "OpenMPI" $log
    webDownload http://www.open-mpi.org/software/ompi/v1.6/downloads/openmpi-1.6.3.tar.gz $log $architecture
    checkError $? "Failed to download openmpi-1.6.3.tar.gz" $log
  fi 
  printOperation "Extracting" "OpenMPI" $log
  tar -xvf openmpi-1.6.3.tar.gz >> $log 2>&1
  checkError $? "Failed to extract openmpi-1.6.3.tar.gz" $log
  mkdir -p openmpi
  cd openmpi
  MPIDir=`pwd`
  cd ../openmpi-1.6.3
  printOperation "Configuring" "OpenMPI" $log
  ./configure --prefix=$MPIDir >> $log 2>&1
  checkError $? "Failed to configure openmpi-1.6.3" $log
  printOperation "Making" "OpenMPI" $log
  make >> $log 2>&1
  checkError $? "Failed to make openmpi-1.6.3" $log
  printOperation "Installing" "OpenMPI" $log
  make install >> $log 2>&1
  checkError $? "Failed to install openmpi-1.6.3" $log
  counter=`expr $counter + 1`
  export PATH=$PATH:$MPIDir/bin
fi 
fi 

# Lattice Microbes #
echo | tee -a $log
echo "$counter) Installing Lattice Microbes" | tee -a $log
echo "------------------------------" | tee -a $log
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$BUILD/protobuf/lib:$BUILD/sbml/usr/lib64:$BUILD/szip/lib:$BUILD/hdf5/lib:/sw/keeneland/openmpi/1.5.1/centos5.5_intel11.1.073/bin/mpicc:$BUILD/openmpi/lib"
mkdir -p $BUILD/LM_CPU
mkdir -p $BUILD/LM_CUDA



if [[ $forCPU -eq 1 ]] ; then
  cd $lmdir
  echo >> $log
  printOperation "Building" "LM with CPU support" $log
  if [[ $architecture == "Keeneland" ]] ; then
    cp docs/config/local.mk.keeneland-CPU-dynamic local.mk
    checkError $? "Failed to find 'local.mk.keeneland-CPU-dynamic'" $log
  elif [[ $architecture == "Titan" ]] ; then
    cp docs/config/local.mk.titan-CPU-dynamic local.mk
    checkError $? "Failed to find 'local.mk.titan-CPU-dynamic'" $log

    #load required modules
#    module load cray-hdf5-parallel/1.8.9
    module load gcc/4.7.2
    module load szip/2.1
    sed "s|/home/erobert3/share/bin/protoc|$BUILD/protobuf/bin/protoc|" local.mk > mk2
    sed "s|PROTOBUF_INCLUDE_DIR := -I/home/erobert3/share/include|PROTOBUF_INCLUDE_DIR := -I$BUILD/protobuf/include|" mk2 > mk1
    sed "s|PROTOBUF_LIB_DIR := -L/home/erobert3/share/lib|PROTOBUF_LIB_DIR := -L$BUILD/protobuf/lib|" mk1 > mk2
    sed "s|/share/apps/hdf5/include|$BUILD/hdf5/include|" mk2 > mk1
    sed "s|/usr/bin/swig|$BUILD/swig/bin/swig|" mk1 > mk2
    sed "s|SBML_INCLUDE_DIR := -I/home/erobert3/share/include|SBML_INCLUDE_DIR := -I$BUILD/sbml/usr/include|" mk2 > mk1
    sed "s|SBML_LIB_DIR := -L/home/erobert3/share/lib|SBML_LIB_DIR := -L$BUILD/sbml/usr/lib64|" mk1 > mk2
    sed "s|/share/apps/hdf5/lib|$BUILD/hdf5/lib|" mk2 > mk1
    sed "s|INSTALL_PREFIX := /home/erobert3/share|INSTALL_PREFIX := $BUILD/LM_CPU|" mk1 > mk2
    mv mk2 local.mk

  elif [[ $architecture == "osx" ]] ; then
    cp docs/config/local.mk.osx-x86_64_cuda-2.x_vmd-1.9.1_static local.mk
    sed 's|Build\-osx\-x86\_64\_cuda\-2\.x\_vmd\-1\.9\.1\_static|LM\_CPU|' local.mk > mk2
    sed "s|/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/bin/protoc|$BUILD/protobuf/bin/protoc|" mk2 > mk1
    sed "s|PROTOBUF_INCLUDE_DIR := -I/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/include| PROTOBUF_INCLUDE_DIR := -I$BUILD/protobuf/include|" mk1 > mk2
    sed "s|HDF5_INCLUDE_DIR := -I/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/include|HDF5_INCLUDE_DIR := -I$BUILD/hdf5/include|" mk2 > mk1
    sed "s|/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/lib/libprotobuf.a|$BUILD/protobuf/lib/libprotobuf.a|" mk1 > mk2
    sed "s|/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/lib/libhdf5.a|$BUILD/hdf5/lib/libhdf5.a|" mk2 > mk1 
    sed "s|/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/lib/libhdf5_hl.a|$BUILD/hdf5/lib/libhdf5_hl.a $BUILD/szip/lib/libsz.a|" mk1 > mk2 
#  sed "s|USE_MPI := 0|USE_MPI := 1|" mk1 > mk2
    sed "s|/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/include|$BUILD/sbml/include|" mk2 > mk1 
    sed "s|/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/lib/libsbml.a|$BUILD/sbml/lib/libsbml.a|" mk1 > mk2 
    sed "s|/usr/bin/swig|$BUILD/swig/bin/swig|" mk2 > mk1
    sed "s|USE_CUDA := 1|USE_CUDA := 0|" mk1 > mk2
    sed "s|USE_VMD := 1|USE_VMD := 0|" mk2 > mk1
    sed "s|INSTALL_PREFIX := |INSTALL_PREFIX := $BUILD/LM_CPU|" mk2 > mk1
    if [[ -d /Applications/VMD1.9.1.app ]] ; then
      sed "s|VMD 1.9|VMD1.9.1|" mk1 > mk2
      mv mk2 mk1
    fi  
    if [[ -d /Applications/VMD\ 1.9.1.app ]] ; then
      sed "s|VMD 1.9|VMD 1.9.1|" mk1 > mk2
      mv mk2 mk1
    fi  
    cp mk1 local.mk
    rm mk1 mk2
  fi
  make -j$cores >> $log 2>&1
  checkError $? "Failed to make Lattice Microbes for CPU" $log
  make install >> $log 2>&1
  checkError $? "Failed to install Lattice Microbes for CPU" $log

  if [[ $architecture == "Titan" ]] ; then
    # Unload modules to be safe
#    module unload cray-hdf5-parallel/1.8.9
    module unload gcc/4.7.2
    module unload szip/2.1
  fi

fi

if [[ $forGPU -eq 1 ]] ; then
  cd $lmdir
  make clean >> $log 2>&1
  rm install
  echo >> $log
  printOperation "Building" "LM with GPU support" $log
  if [[ $architecture == "Keeneland" ]] ; then
    cp docs/config/local.mk.keeneland-CUDA-dynamic local.mk
    checkError $? "Failed to find 'local.mk.keeneland-CUDA-dynamic'" $log
  elif [[ $architecture == "Titan" ]] ; then
    cp docs/config/local.mk.titan-CUDA-dynamic local.mk
    checkError $? "Failed to find 'local.mk.titan-CUDA-dynamic'" $log

    #load required modules
#    module load cray-hdf5-parallel/1.8.9
    module load gcc/4.7.2
    module load szip/2.1
    sed "s|/home/erobert3/share/bin/protoc|$BUILD/protobuf/bin/protoc|" local.mk > mk2
    sed "s|PROTOBUF_INCLUDE_DIR := -I/home/erobert3/share/include|PROTOBUF_INCLUDE_DIR := -I$BUILD/protobuf/include|" mk2 > mk1
    sed "s|PROTOBUF_LIB_DIR := -L/home/erobert3/share/lib|PROTOBUF_LIB_DIR := -L$BUILD/protobuf/lib|" mk1 > mk2
    sed "s|/share/apps/hdf5/include|$BUILD/hdf5/include|" mk2 > mk1
    sed "s|/usr/bin/swig|$BUILD/swig/bin/swig|" mk1 > mk2
    sed "s|SBML_INCLUDE_DIR := -I/home/erobert3/share/include|SBML_INCLUDE_DIR := -I$BUILD/sbml/usr/include|" mk2 > mk1
    sed "s|SBML_LIB_DIR := -L/home/erobert3/share/lib|SBML_LIB_DIR := -L$BUILD/sbml/usr/lib64|" mk1 > mk2
    sed "s|/share/apps/hdf5/lib|$BUILD/hdf5/lib|" mk2 > mk1
    sed "s|INSTALL_PREFIX := /home/erobert3/share|INSTALL_PREFIX := $BUILD/LM_CUDA|" mk1 > mk2
    sed "s|/usr/local/cuda/include|/opt/nvidia/cudatoolkit/5.0.35.102/include|" mk2 > mk1
    sed "s|/usr/local/cuda/lib64|/opt/nvidia/cudatoolkit/5.0.35.102/lib64 -L/opt/cray/nvidia/default/lib64|" mk1 > mk2
    sed "s|/usr/local/cuda/bin/nvcc|/opt/nvidia/cudatoolkit/5.0.35.102/bin/nvcc|" mk2 > mk1
    mv mk1 local.mk
  elif [[ $architecture == "osx" ]] ; then
    cp docs/config/local.mk.osx-x86_64_cuda-2.x_vmd-1.9.1_static local.mk
    sed 's|Build\-osx\-x86\_64\_cuda\-2\.x\_vmd\-1\.9\.1\_static|LM\_CPU|' local.mk > mk2
    sed "s|/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/bin/protoc|$BUILD/protobuf/bin/protoc|" mk2 > mk1
    sed "s|PROTOBUF_INCLUDE_DIR := -I/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/include| PROTOBUF_INCLUDE_DIR := -I$BUILD/protobuf/include|" mk1 > mk2
    sed "s|HDF5_INCLUDE_DIR := -I/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/include|HDF5_INCLUDE_DIR := -I$BUILD/hdf5/include|" mk2 > mk1
    sed "s|/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/lib/libprotobuf.a|$BUILD/protobuf/lib/libprotobuf.a|" mk1 > mk2
    sed "s|/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/lib/libhdf5.a|$BUILD/hdf5/lib/libhdf5.a|" mk2 > mk1
    sed "s|/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/lib/libhdf5_hl.a|$BUILD/hdf5/lib/libhdf5_hl.a $BUILD/szip/lib/libsz.a|" mk1 > mk2
#  sed "s|USE_MPI := 0|USE_MPI := 1|" mk1 > mk2
    sed "s|/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/include|$BUILD/sbml/include|" mk2 > mk1
    sed "s|/Network/Servers/sol.scs.uiuc.edu/Volumes/HomeRAID2/Homes/erobert3/usr/Darwin-i386/lib/libsbml.a|$BUILD/sbml/lib/libsbml.a|" mk1 > mk2
    sed "s|/usr/bin/swig|$BUILD/swig/bin/swig|" mk2 > mk1
    if [[ -d /Applications/VMD1.9.1.app ]] ; then
      sed "s|VMD 1.9|VMD1.9.1|" mk1 > mk2
      sed "s|VMD_INSTALL_DIR :=|VMD_INSTALL_DIR := /Applications/VMD1.9.1.app/Contents/vmd/plugins/MACOSXX86/molfile/|" mk2 > mk1
      mv mk1 mk2
    fi  
    if [[ -d "/Applications/VMD\ 1.9.1.app" ]] ; then
      sed "s|VMD 1.9|VMD 1.9.1|" mk1 > mk2
      sed "s|VMD_INSTALL_DIR :=|VMD_INSTALL_DIR := /Applications/VMD\ 1.9.1.app/Contents/vmd/plugins/MACOSXX86/molfile/|" mk2 > mk1
      mv mk1 mk2
    fi  
    sed "s|INSTALL_PREFIX := |INSTALL_PREFIX := $BUILD/LM_CUDA|" mk2 > mk1
    cp mk2 local.mk
    rm mk2 mk1
  fi
  make -j$cores >> $log 2>&1
  checkError $? "Failed to make Lattice Microbes for CPU" $log
  if [[ $architecture == "osx" ]] ; then
    echo
    echo "In order to isntall Lattice Microbes with support for VMD, administrative access is needed."
    echo "Please type your password now to complete installation:"
    echo
    sudo make install >> $log 2>&1
    checkError $? "Failed to install Lattice Microbes for GPU" $log
  else 
    make install >> $log 2>&1
    checkError $? "Failed to install Lattice Microbes for GPU" $log
  fi

  if [[ $architecture == "Titan" ]] ; then
    # Unload modules to be safe
#    module unload cray-hdf5-parallel/1.8.9
    module unload gcc/4.7.2
    module unload szip/2.1
  fi
fi

echo
echo
echo "The LD_LIBRARY_PATH variable should now be updated.  To do this run and add to your ~/.bash_profile:"
echo
echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
echo

if [[ $architecture == "Titan" ]] ; then
  echo "To run on titan you must use GNU compiler version 4.7.* and several other modules."
  echo "For these add the following module lines to your login script: "
  echo 
  echo "module load cray-hdf5-parallel/1.8.9"
  echo "module load gcc/4.7.2"
  echo "module load szip/2.1"
  echo
  echo "If this is unacceptable, you should build Lattice Microbes from scratch."
  echo 
fi

rm install

