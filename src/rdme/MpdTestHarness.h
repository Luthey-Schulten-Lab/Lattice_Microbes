/*
 * University of Illinois Open Source License
 * Copyright 2015-2018 Luthey-Schulten Group,
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
 * Urbana-Champaign, the Roberts Group, Johns Hopkins University, nor the names
 * of its contributors may be used to endorse or promote products derived from
 * this Software without specific prior written permission.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Mike Hallock
 */

#ifndef LM_RDME_MPDALTHARNESS
#define LM_RDME_MPDALTHARNESS

#include "cuda/lm_cuda.h"
#include "rdme/RDMESolver.h"
#include "rdme/MpdRdmeSolver.h"
#include "rdme/CudaByteLattice.h"

using lm::rdme::RDMESolver;
using lm::rdme::Lattice;

namespace lm {

namespace io{
class Lattice;
class SpeciesCounts;
}
namespace rdme {

class MpdTestHarness : public MpdRdmeSolver 
{
public:
    MpdTestHarness();
    virtual ~MpdTestHarness();
    virtual void generateTrajectory();

protected:
    virtual void runTimestep(CudaByteLattice * lattice, uint32_t timestep);

	cudaEvent_t original_start, original_end, jit_start, jit_end;
	float total_jit, total_orig;

};

}
}

#endif
