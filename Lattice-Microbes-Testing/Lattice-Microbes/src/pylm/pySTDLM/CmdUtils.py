
#!/usr/bin/env python

import os
import argparse

def estimateGPUMemory():
    # Get arguments
    parser = argparse.ArgumentParser(description="Estimate the GPU memory required for a RDME simulation.")
    parser.add_argument('-d', '--dims', type=int, nargs=3, help="The dimensions of the lattice; e.g. --dims x y z")
    parser.add_argument('-p', '--pps',  default=8, type=int, required=False, help="Particles per site.")
    parser.add_argument('-g', '--gpus', default=1, type=int, help="Numer of GPUs to use (to correctly compute overlap.")
    args = parser.parse_args()

    # Constants
    bytesPerParticle = 1
    bytesPerSite     = 1

    # Get lattice dimensions
    x = args.dims[0]
    y = args.dims[1]
    z = args.dims[2]
    if args.gpus > 1:
        z = round(z/args.gpus) + 2

    # Compute memory
    latticeSites = x*y*z
    memoryDiffusionLattice = 2*bytesPerParticle*args.pps*latticeSites
    memorySiteLattice = bytesPerSite*latticeSites
    totalMemory = float(memoryDiffusionLattice + memorySiteLattice)

    unit = "B"
    if totalMemory > 2**30:
        unit = 'GB'
        totalMemory /= 2**30
    elif totalMemory > 2**20:
        unit = 'MB'
        totalMemory /= 2**20
    elif totalMmory > 2**10:
        unit = 'KB'
        totalMemory /= 2**10

    # Print information
    print("Lattice Size:\t\t(%d,%d,%d)"%(x,y,z))
    print("Lattice Per GPU:\t%d"%(latticeSites))
    print("Memory Required:\t%0.2f%s"%(totalMemory,unit))

def peekFile():
    parser = argparse.ArgumentParser(description="Print info on simulation file")
    parser.add_argument('file', metavar='LM_HDF5', type=str, help='HDF5 simulation file')
    args = parser.parse_args()

    with  h5py.File(args.file) as hdf:
        tf = float(hdf['Parameters'].attrs['maxTime'])
        for r,sim in hdf['Simulations'].items():
            t = sim['LatticeTimes'][-1]
            n = len(sim['Lattice'])
            print("%s:%s: %d frames at %.2f sec.; %.2f complete" % (args.file, r, n, t, 100*t/tf) )
