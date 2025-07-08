# Portable conda build
## Prerequisites
Install Anaconda 3 if not already installed on your computer, download and follow the instructions at
https://www.anaconda.com/

If you want to do spatial RDME simulations which require GPUs, install CUDA, download and instructions can be found at
https://developer.nvidia.com/cuda-downloads

Add the following line at the end of /.bashrc so LM can find CUDA.  If path to CUDA is different,

```bash
export path to where CUDA is installed:
export PATH="/usr/local/cuda/bin/:$PATH"
```

## LM installation via conda
The conda environment can be created from the yml file included here.  Go into the file and change "lmpretest" to whatever you want your environment name to be, then run:

replace `envname` with your choice of name for the environment.

```bash
    conda env create -n envname -f lm_precomp.yml
    
    conda activate envname
```

**NOT RECOMMENDED(deprecated)**: You can also make your Anaconda environment with your choice of name "envname" using Python 3 through the following procedure:

```bash
    conda create -n envname -c conda-forge python=3.7.3
    conda activate envname
    conda install -c conda-forge -y git==2.33.1 gcc_linux-64==7.3.0 gxx_linux-64==7.3.0 binutils_linux-64==2.31.1 libstdcxx-ng==11.2.0 libgcc-ng==11.2.0
    conda install -c anaconda hdf5==1.10.4 h5py==2.10.0 protobuf==3.13.0.1
    conda install -c conda-forge -y numpy==1.19.2 cython==0.29.4 cmake==3.14.0 xlrd==1.2.0 swig==4.0.2 jupyter==1.0.0 matplotlib==3.4.3 ipywidgets==7.6.5 tqdm==4.62.3 pillow==8.3.2 jinja2==3.0.2 scipy==1.5.2 pybind11==2.8.1 pandas==1.3.4 pytables==3.5.2 biopython==1.79
```

Once you have created and activated your environment, go to a directory where you want to build Lattice Microbes and do the build running the following commands:
(You do not need to run simulations where the build is located, lm will be installed in your environment as an executable)

```bash
cd .. (make sure you are in the directory which has Lattice-Microbes folder)
```

create a build directory:

```bash
mkdir build
cd build
cmake ../Lattice-Microbes/src/ -D MPD_GLOBAL_T_MATRIX=True -D MPD_GLOBAL_R_MATRIX=True
# The -D flags can be removed if you build with smaller transition matrix and reaction matrix size limits in the cmake options.
make && make install
# Run make install as root or sudo if possible
```

manually set CUDA_ARCHITECTURES if needed:

```bash
cmake ../Lattice-Microbes/src/ -D MPD_GLOBAL_T_MATRIX=True -D MPD_GLOBAL_R_MATRIX=True -D CUDA_ARCHITECTURES="60;70;75;80;86"
make && make install
```

## Build for docker container
It is important to make sure the CUDA_ARCHITECTURES argument matches the architecture of the GPU you are using. If you are not sure what architecture your GPU is, you can run:
```bash
nvidia-smi
```
And check compute capability of your GPU. You can find the information in the last two digits of the GPU model. For example, RTX 3060 Ti has compute capability 86.
```bash
docker build --build-arg CUDA_ARCHITECTURES="70;75;80;86" -t lm2.5dev:tag .
```
## Build for apptainer
For x86_64 architecture:
```bash
apptainer build your-container.sif apptainer.def
```
For aarch64 architecture:
```bash
apptainer build your-container.sif apptainer-dtai.def
```

at run time, specify the architecture:
```bash
   CUDA_ARCHITECTURES="75;80;86" apptainer run your-container.sif
```
## Testing the installation

Check to make sure the installation worked, run:

```bash
lm --help
# Should return help commands for lm
```

There are example executable simulation files in example_rdme_files you can try running.

```bash
cd example_rdme_files
python run_simulation.py
```