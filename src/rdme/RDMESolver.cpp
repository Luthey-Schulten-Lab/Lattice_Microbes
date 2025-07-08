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

#include "config.h"
#include "core/Exceptions.h"
#include "core/Print.h"
#include "cme/CMESolver.h"
#include "DiffusionModel.pb.h"
#include "rdme/Lattice.h"
#include "rdme/ByteLattice.h"
#include "rdme/RDMESolver.h"
#include "rng/RandomGenerator.h"
#include "lptf/Profile.h"

using lm::io::DiffusionModel;
using lm::rdme::Lattice;
using lm::rng::RandomGenerator;

namespace lm {
namespace rdme {

RDMESolver::RDMESolver(RandomGenerator::Distributions neededDists)
:CMESolver(neededDists),numberSiteTypes(0),DF(NULL),RL(NULL),lattice(NULL)
{
}

RDMESolver::~RDMESolver()
{
    // Free any model memory.
    destroyDiffusionModel();
}

void RDMESolver::allocateDiffusionModel(uint numberSiteTypesA, lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, unsigned int bytes_per_particle, si_dist_t latticeSpacing)
{
    // Set the number of site types.
    numberSiteTypes = numberSiteTypesA;

    // Allocate the resources.
    DF = new double[numberSpecies*numberSiteTypes*numberSiteTypes];
    RL = new uint[numberReactions*numberSiteTypes];
    allocateLattice(latticeXSize, latticeYSize, latticeZSize, particlesPerSite, bytes_per_particle, latticeSpacing);
}

void RDMESolver::allocateLattice(lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, unsigned int bytes_per_particle, si_dist_t latticeSpacing)
{
	assert(bytes_per_particle == 1);
    lattice = new ByteLattice(latticeXSize, latticeYSize, latticeZSize, latticeSpacing, particlesPerSite);
}

void RDMESolver::destroyDiffusionModel()
{
    // Free any resources.
    if (DF != NULL) {delete[] DF; DF = NULL;}
    if (RL != NULL) {delete[] RL; RL = NULL;}
    if (lattice != NULL) {delete lattice; lattice = NULL;}

    // Reset the number of sites types.
    numberSiteTypes = 0;
}

void RDMESolver::setDiffusionModel(DiffusionModel * dm, const uint8_t * lattice, size_t latticeSize, const uint8_t * latticeSites, size_t latticeSitesSize)
{
    // Validate the model.
    if (dm->number_species() != numberSpecies) throw InvalidArgException("dm.number_species", "number of species in the diffusion model does not agree with the number in the reaction model");
    if (dm->number_reactions() != numberReactions) throw InvalidArgException("dm.number_reactions", "number of reactions in the diffusion model does not agree with the number in the reaction model");
    if ((uint)dm->diffusion_matrix_size() != dm->number_species()*dm->number_site_types()*dm->number_site_types()) throw InvalidArgException("dm", "diffusion matrix size does not agree with the number of species and site types");
    if (size_t(dm->lattice_x_size())*size_t(dm->lattice_y_size())*size_t(dm->lattice_z_size())*size_t(dm->particles_per_site())*size_t(dm->bytes_per_particle()) != latticeSize) throw InvalidArgException("latticeSize", "the lattice data size does not agree with the lattice dimensions", dm->lattice_x_size()*dm->lattice_y_size()*dm->lattice_z_size()*dm->particles_per_site(), latticeSize);
    if (dm->lattice_x_size()*dm->lattice_y_size()*dm->lattice_z_size() != latticeSitesSize) throw InvalidArgException("latticeSitesSize", "the lattice data size does not agree with the lattice dimensions");

    // Build the model.
    buildDiffusionModel(dm->number_site_types(), dm->diffusion_matrix().data(), dm->reaction_location_matrix().data(), dm->lattice_x_size(), dm->lattice_y_size(), dm->lattice_z_size(), dm->particles_per_site(), dm->bytes_per_particle(), dm->lattice_spacing(), lattice, latticeSites);
}

void RDMESolver::buildDiffusionModel(const uint numberSiteTypesA, const double * DFA, const uint * RLA, lattice_size_t latticeXSize, lattice_size_t latticeYSize, lattice_size_t latticeZSize, site_size_t particlesPerSite, const unsigned int bytes_per_particle, si_dist_t latticeSpacing, const uint8_t * latticeData, const uint8_t * latticeSitesData, bool rowMajorData)
{
    // Destroy the previous model, if we have one.
    destroyDiffusionModel();

    // Allocate space for the new model.
    allocateDiffusionModel(numberSiteTypesA, latticeXSize, latticeYSize, latticeZSize, particlesPerSite, bytes_per_particle, latticeSpacing);

    // Set the diffusion matrix.
    for (uint i=0; i<numberSpecies*numberSiteTypes*numberSiteTypes; i++)
    {
        DF[i] = DFA[i];
    }

    // Set the reaction location matrix.
    for (uint i=0; i<numberReactions*numberSiteTypes; i++)
    {
        RL[i] = RLA[i];
    }

    // Set the lattice.
    if (rowMajorData)
    {
		setLatticeData(latticeData);
		setLatticeSitesData(latticeSitesData);
	}
    else
    {
        throw lm::InvalidArgException("rowMajorData","columns major lattice data is not currently supported");
    }

    Print::printf(Print::DEBUG, "Set diffusion model.");
}

void RDMESolver::setLatticeData(const uint8_t* latticeData)
{
    lattice_coord_t s = lattice->getSize();
    site_size_t p = lattice->getMaxOccupancy();

    // Set the lattice data.
    for (uint i=0, index=0; i<s.x; i++)
    {
        for (uint j=0; j<s.y; j++)
        {
            for (uint k=0; k<s.z; k++)
            {
                for (uint l=0; l<p; l++, index++)
                {
                    if (latticeData[index] != 0)
                    {
                        if (latticeData[index] > numberSpecies) throw InvalidArgException("latticeData", "an invalid species was found",latticeData[index]);
                        lattice->addParticle(i,j,k,latticeData[index]);
                    }
                }
            }
        }
    }
}

void RDMESolver::setLatticeSitesData(const uint8_t* latticeSitesData)
{
    lattice_coord_t s = lattice->getSize();
    for (uint i=0, index=0; i<s.x; i++)
    {
        for (uint j=0; j<s.y; j++)
        {
            for (uint k=0; k<s.z; k++, index++)
            {
                if (latticeSitesData[index] != 0)
                {
                    if (latticeSitesData[index] >= numberSiteTypes)
                    {
                        throw InvalidArgException("latticeSitesData", "an invalid species was found",latticeSitesData[index]);
                    }
                    lattice->setSiteType(i,j,k,latticeSitesData[index]);
                }
            }
        }
    }
}

}
}
