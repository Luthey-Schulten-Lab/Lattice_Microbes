/*
 * University of Illinois Open Source License
 * Copyright 2011-2018 Luthey-Schulten Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the Software), to deal with 
 * the Software without restriction, including without limitation the rights to 
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
 * of the Software, and to permit persons to whom the Software is furnished to 
 * do so, subject to the following conditions:
 * 
 * - Redistributions of source code must retain the above copyright notice, 
 * this list of conditions and the following disclaimers.
 * 
 * - Redistributions in binary form must reproduce the above copyright notice, 
 * this list of conditions and the following disclaimers in the documentation 
 * and/or other materials provided with the distribution.
 * 
 * - Neither the names of the Luthey-Schulten Group, University of Illinois at
 * Urbana-Champaign, nor the names of its contributors may be used to endorse or
 * promote products derived from this Software without specific prior written
 * permission.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Elijah Roberts
 */

#ifndef LM_RDME_RDMESOLVER_H_
#define LM_RDME_RDMESOLVER_H_

#include "cme/CMESolver.h"
#include "DiffusionModel.pb.h"
#include "rdme/Lattice.h"

using lm::cme::CMESolver;
using lm::io::DiffusionModel;
using lm::rdme::Lattice;

namespace lm {
namespace rdme {

class RDMESolver : public CMESolver
{
public:
    RDMESolver(RandomGenerator::Distributions neededDists);
    virtual ~RDMESolver();
    virtual void setDiffusionModel(DiffusionModel * dm, const uint8_t * lattice, size_t latticeSize, const uint8_t * latticeSites, size_t latticeSitesSize);
    virtual void buildDiffusionModel(const uint numberSiteTypesA, const double * DFA, const uint * RLA, lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing, const uint8_t * latticeData, const uint8_t * latticeSitesData, bool rowMajorData=true);

protected:
    virtual void allocateDiffusionModel(uint numberSiteTypesA, lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, unsigned int bytes_per_particle, si_dist_t latticeSpacing);
    virtual void allocateLattice(lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, unsigned int bytes_per_particle, si_dist_t latticeSpacing);
	virtual void setLatticeData(const uint8_t* latticeData);
	virtual void setLatticeSitesData(const uint8_t* latticeSitesData);
    virtual void destroyDiffusionModel();

protected:
    uint numberSiteTypes;
    double * DF;                            // diffusion matrix: numberSpecies x numberSiteTypes x numberSiteTypes
    uint * RL;								// reaction location matrix: numberReactions x numberSiteTypes
    Lattice * lattice;
};

}
}

#endif
