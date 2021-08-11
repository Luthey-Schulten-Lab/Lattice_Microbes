/*
 * University of Illinois Open Source License
 * Copyright 2010-2018 Luthey-Schulten Group,
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

#ifndef LM_CME_SELFREGULATINGGENESWITCH_H_
#define LM_CME_SELFREGULATINGGENESWITCH_H_

#include <cmath>
#include <map>
#include <string>
#include "cme/NextReactionSolver.h"
#include "ParameterValues.pb.h"

using std::string;
using std::map;
using lm::io::ParameterValues;

namespace lm {
namespace cme {

class SelfRegulatingGeneSwitch: public NextReactionSolver
{
protected:
    struct OUKHillPropensityArgs
    {
        OUKHillPropensityArgs(uint xi, uint oui, double kmin, double kmax, double x50, double noiseVariance, double noiseTau, RandomGenerator * rng):xi(xi),oui(oui),kmin(kmin),dk(kmax-kmin),x50sq(x50*x50),noiseVariance(noiseVariance),noiseTau(noiseTau),noise(0.0),lastOUJumpNumber(0),previousTime(0.0),rng(rng),rngNext(TUNE_LOCAL_RNG_CACHE_SIZE) {}
        ~OUKHillPropensityArgs() {}
        uint xi;
        uint oui;
        double kmin;
        double dk;
        double x50sq;
        double noiseVariance;
        double noiseTau;
        double noise;
        uint lastOUJumpNumber;
        double previousTime;
        RandomGenerator * rng;
        double normRngValues[TUNE_LOCAL_RNG_CACHE_SIZE];
        int rngNext;
    };

public:
    SelfRegulatingGeneSwitch();
    virtual ~SelfRegulatingGeneSwitch();
    virtual bool needsReactionModel() {return false;}
    virtual bool needsDiffusionModel()  {return false;}
    virtual void generateTrajectory();
	
protected:
    static double ouKHillPropensity(double time, uint * speciesCounts, void * pargs);
};

}
}

#endif
