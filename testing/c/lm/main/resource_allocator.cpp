/*
 * University of Illinois Open Source License
 * Copyright 2011 Luthey-Schulten Group,
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

#include <vector>
#include "TestingHelper.h"
#include "lm/Exceptions.h"
#include "lm/Types.h"
#include "lm/io/SpeciesCounts.pb.h"
#include "lm/main/ResourceAllocator.h"
#define BOOST_TEST_MODULE LatticeMicrobe
#include <boost/test/unit_test.hpp>

using std::vector;
using lm::main::ResourceAllocator;

class ResourceAllocatorTester : public ResourceAllocator
{
public:
    ResourceAllocatorTester(int numberCpuCoresA, float cpuCoresPerReplicateA, vector<int> cudaDevicesA, float cudaDevicesPerReplicateA)
    :ResourceAllocator(numberCpuCoresA,cpuCoresPerReplicateA,cudaDevicesA,cudaDevicesPerReplicateA),numberCpuCoresP(&numberCpuCores),cudaDevicesP(&cudaDevices),cpuSlotsPerCoreP(&cpuSlotsPerCore),cpuSlotsPerReplicateP(&cpuSlotsPerReplicate),cudaSlotsPerDeviceP(&cudaSlotsPerDevice),cudaSlotsPerReplicateP(&cudaSlotsPerReplicate),maxSimultaneousReplicatesP(&maxSimultaneousReplicates),cpuSlotsP(&cpuSlots),cudaSlotsP(&cudaSlots)
    {}

    int * numberCpuCoresP;
    vector<int> * cudaDevicesP;
    int * cpuSlotsPerCoreP;
    int * cpuSlotsPerReplicateP;
    int * cudaSlotsPerDeviceP;
    int * cudaSlotsPerReplicateP;
    int * maxSimultaneousReplicatesP;
    int *** cpuSlotsP;
    int *** cudaSlotsP;
};

BOOST_AUTO_TEST_SUITE(ResourceAllocatorTest)

BOOST_AUTO_TEST_CASE(Constructor)
{
    ResourceAllocatorTester * t = NULL;
    vector<int> cudaDevices;

    // Test a simple constructor with 1:1 mapping.
    cudaDevices.clear(); cudaDevices.push_back(0);
    BOOST_REQUIRE_NO_THROW(t=new ResourceAllocatorTester(1,1.0,cudaDevices,1.0));
    BOOST_CHECK_EQUAL(*(t->numberCpuCoresP), 1);
    BOOST_CHECK_EQUAL(t->cudaDevicesP->size(), 1);
    BOOST_CHECK_EQUAL((*(t->cudaDevicesP))[0], 0);
    BOOST_CHECK_EQUAL(*(t->cpuSlotsPerCoreP), 1);
    BOOST_CHECK_EQUAL(*(t->cpuSlotsPerReplicateP), 1);
    BOOST_CHECK_EQUAL(*(t->cudaSlotsPerDeviceP), 1);
    BOOST_CHECK_EQUAL(*(t->cudaSlotsPerReplicateP), 1);
    BOOST_CHECK_EQUAL(*(t->maxSimultaneousReplicatesP), 1);
    BOOST_CHECK_EQUAL(t->getMaxSimultaneousReplicates(), *(t->maxSimultaneousReplicatesP));
    for (int i=0; i<*(t->numberCpuCoresP); i++)
    {
        for (int j=0; j<*(t->cpuSlotsPerCoreP); j++)
        {
            BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[i][j], -1);
        }

    }
    for (int i=0; i<(int)t->cudaDevicesP->size(); i++)
    {
        for (int j=0; j<*(t->cudaSlotsPerDeviceP); j++)
        {
            BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[i][j], -1);
        }

    }
    if (t != NULL) delete t; t = NULL;

    // Test a constructor without any resources.
    cudaDevices.clear();
    BOOST_REQUIRE_THROW(t=new ResourceAllocatorTester(0,0.0,cudaDevices,0.0), lm::Exception);

    // Test a constructor without any resources.
    cudaDevices.clear();
    BOOST_REQUIRE_THROW(t=new ResourceAllocatorTester(10,0.0,cudaDevices,0.0), lm::Exception);

    // Test a constructor without any resources.
    cudaDevices.clear(); cudaDevices.push_back(0);
    BOOST_REQUIRE_THROW(t=new ResourceAllocatorTester(0,0.0,cudaDevices,0.0), lm::Exception);

    // Test a constructor without any resources.
    cudaDevices.clear(); cudaDevices.push_back(0);
    BOOST_REQUIRE_THROW(t=new ResourceAllocatorTester(10,0.0,cudaDevices,0.0), lm::Exception);

    // Test a constructor with only multiple cpus.
     cudaDevices.clear(); cudaDevices.push_back(0);
     BOOST_REQUIRE_NO_THROW(t=new ResourceAllocatorTester(23,1.0,cudaDevices,0.0));
     BOOST_CHECK_EQUAL(*(t->numberCpuCoresP), 23);
     BOOST_CHECK_EQUAL(t->cudaDevicesP->size(), 1);
     BOOST_CHECK_EQUAL(*(t->cpuSlotsPerCoreP), 1);
     BOOST_CHECK_EQUAL(*(t->cpuSlotsPerReplicateP), 1);
     BOOST_CHECK_EQUAL(*(t->cudaSlotsPerDeviceP), 0);
     BOOST_CHECK_EQUAL(*(t->cudaSlotsPerReplicateP), 0);
     BOOST_CHECK_EQUAL(*(t->maxSimultaneousReplicatesP), 23);
     BOOST_CHECK_EQUAL(t->getMaxSimultaneousReplicates(), *(t->maxSimultaneousReplicatesP));
     for (int i=0; i<*(t->numberCpuCoresP); i++)
     {
         for (int j=0; j<*(t->cpuSlotsPerCoreP); j++)
         {
             BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[i][j], -1);
         }

     }
     for (int i=0; i<(int)t->cudaDevicesP->size(); i++)
     {
         for (int j=0; j<*(t->cudaSlotsPerDeviceP); j++)
         {
             BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[i][j], -1);
         }

     }
     if (t != NULL) delete t; t = NULL;

     // Test a constructor with only multiple cpus and multiple cpus per slot.
      cudaDevices.clear();
      BOOST_REQUIRE_NO_THROW(t=new ResourceAllocatorTester(23,4.0,cudaDevices,0.0));
      BOOST_CHECK_EQUAL(*(t->numberCpuCoresP), 23);
      BOOST_CHECK_EQUAL(t->cudaDevicesP->size(), 0);
      BOOST_CHECK_EQUAL(*(t->cpuSlotsPerCoreP), 1);
      BOOST_CHECK_EQUAL(*(t->cpuSlotsPerReplicateP), 4);
      BOOST_CHECK_EQUAL(*(t->cudaSlotsPerDeviceP), 0);
      BOOST_CHECK_EQUAL(*(t->cudaSlotsPerReplicateP), 0);
      BOOST_CHECK_EQUAL(*(t->maxSimultaneousReplicatesP), 5);
      BOOST_CHECK_EQUAL(t->getMaxSimultaneousReplicates(), *(t->maxSimultaneousReplicatesP));
      for (int i=0; i<*(t->numberCpuCoresP); i++)
      {
          for (int j=0; j<*(t->cpuSlotsPerCoreP); j++)
          {
              BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[i][j], -1);
          }

      }
      for (int i=0; i<(int)t->cudaDevicesP->size(); i++)
      {
          for (int j=0; j<*(t->cudaSlotsPerDeviceP); j++)
          {
              BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[i][j], -1);
          }

      }
      if (t != NULL) delete t; t = NULL;

      // Test a constructor with only multiple cpus and multiple slots per cpu.
       cudaDevices.clear();
       BOOST_REQUIRE_NO_THROW(t=new ResourceAllocatorTester(23,0.25,cudaDevices,0.0));
       BOOST_CHECK_EQUAL(*(t->numberCpuCoresP), 23);
       BOOST_CHECK_EQUAL(t->cudaDevicesP->size(), 0);
       BOOST_CHECK_EQUAL(*(t->cpuSlotsPerCoreP), 4);
       BOOST_CHECK_EQUAL(*(t->cpuSlotsPerReplicateP), 1);
       BOOST_CHECK_EQUAL(*(t->cudaSlotsPerDeviceP), 0);
       BOOST_CHECK_EQUAL(*(t->cudaSlotsPerReplicateP), 0);
       BOOST_CHECK_EQUAL(*(t->maxSimultaneousReplicatesP), 92);
       BOOST_CHECK_EQUAL(t->getMaxSimultaneousReplicates(), *(t->maxSimultaneousReplicatesP));
       for (int i=0; i<*(t->numberCpuCoresP); i++)
       {
           for (int j=0; j<*(t->cpuSlotsPerCoreP); j++)
           {
               BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[i][j], -1);
           }

       }
       for (int i=0; i<(int)t->cudaDevicesP->size(); i++)
       {
           for (int j=0; j<*(t->cudaSlotsPerDeviceP); j++)
           {
               BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[i][j], -1);
           }

       }
       if (t != NULL) delete t; t = NULL;

       // Test a constructor with only multiple gpus.
        cudaDevices.clear(); cudaDevices.push_back(0); cudaDevices.push_back(1); cudaDevices.push_back(2);
        BOOST_REQUIRE_NO_THROW(t=new ResourceAllocatorTester(11,0.0,cudaDevices,1.0));
        BOOST_CHECK_EQUAL(*(t->numberCpuCoresP), 11);
        BOOST_CHECK_EQUAL(t->cudaDevicesP->size(), 3);
        BOOST_CHECK_EQUAL(*(t->cpuSlotsPerCoreP), 0);
        BOOST_CHECK_EQUAL(*(t->cpuSlotsPerReplicateP), 0);
        BOOST_CHECK_EQUAL(*(t->cudaSlotsPerDeviceP), 1);
        BOOST_CHECK_EQUAL(*(t->cudaSlotsPerReplicateP), 1);
        BOOST_CHECK_EQUAL(*(t->maxSimultaneousReplicatesP), 3);
        BOOST_CHECK_EQUAL(t->getMaxSimultaneousReplicates(), *(t->maxSimultaneousReplicatesP));
        for (int i=0; i<*(t->numberCpuCoresP); i++)
        {
            for (int j=0; j<*(t->cpuSlotsPerCoreP); j++)
            {
                BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[i][j], -1);
            }

        }
        for (int i=0; i<(int)t->cudaDevicesP->size(); i++)
        {
            for (int j=0; j<*(t->cudaSlotsPerDeviceP); j++)
            {
                BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[i][j], -1);
            }

        }
        if (t != NULL) delete t; t = NULL;

        // Test a constructor with only multiple gpus and multiple gpus per slot.
         cudaDevices.clear(); cudaDevices.push_back(0); cudaDevices.push_back(1); cudaDevices.push_back(2);
         BOOST_REQUIRE_NO_THROW(t=new ResourceAllocatorTester(0,0.0,cudaDevices,2.0));
         BOOST_CHECK_EQUAL(*(t->numberCpuCoresP), 0);
         BOOST_CHECK_EQUAL(t->cudaDevicesP->size(), 3);
         BOOST_CHECK_EQUAL(*(t->cpuSlotsPerCoreP), 0);
         BOOST_CHECK_EQUAL(*(t->cpuSlotsPerReplicateP), 0);
         BOOST_CHECK_EQUAL(*(t->cudaSlotsPerDeviceP), 1);
         BOOST_CHECK_EQUAL(*(t->cudaSlotsPerReplicateP), 2);
         BOOST_CHECK_EQUAL(*(t->maxSimultaneousReplicatesP), 1);
         BOOST_CHECK_EQUAL(t->getMaxSimultaneousReplicates(), *(t->maxSimultaneousReplicatesP));
         for (int i=0; i<*(t->numberCpuCoresP); i++)
         {
             for (int j=0; j<*(t->cpuSlotsPerCoreP); j++)
             {
                 BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[i][j], -1);
             }

         }
         for (int i=0; i<(int)t->cudaDevicesP->size(); i++)
         {
             for (int j=0; j<*(t->cudaSlotsPerDeviceP); j++)
             {
                 BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[i][j], -1);
             }

         }
         if (t != NULL) delete t; t = NULL;

         // Test a constructor with only multiple gpus and multiple gpus per slot.
          cudaDevices.clear(); cudaDevices.push_back(0); cudaDevices.push_back(1); cudaDevices.push_back(2);
          BOOST_REQUIRE_NO_THROW(t=new ResourceAllocatorTester(0,0.0,cudaDevices,0.5));
          BOOST_CHECK_EQUAL(*(t->numberCpuCoresP), 0);
          BOOST_CHECK_EQUAL(t->cudaDevicesP->size(), 3);
          BOOST_CHECK_EQUAL(*(t->cpuSlotsPerCoreP), 0);
          BOOST_CHECK_EQUAL(*(t->cpuSlotsPerReplicateP), 0);
          BOOST_CHECK_EQUAL(*(t->cudaSlotsPerDeviceP), 2);
          BOOST_CHECK_EQUAL(*(t->cudaSlotsPerReplicateP), 1);
          BOOST_CHECK_EQUAL(*(t->maxSimultaneousReplicatesP), 6);
          BOOST_CHECK_EQUAL(t->getMaxSimultaneousReplicates(), *(t->maxSimultaneousReplicatesP));
          for (int i=0; i<*(t->numberCpuCoresP); i++)
          {
              for (int j=0; j<*(t->cpuSlotsPerCoreP); j++)
              {
                  BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[i][j], -1);
              }

          }
          for (int i=0; i<(int)t->cudaDevicesP->size(); i++)
          {
              for (int j=0; j<*(t->cudaSlotsPerDeviceP); j++)
              {
                  BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[i][j], -1);
              }

          }
          if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_CASE(AssignReplicate)
{
    ResourceAllocatorTester * t = NULL;
    vector<int> cudaDevices;
    ResourceAllocator::ComputeResources resources;

    // Test a simple assignment policy.
    cudaDevices.clear(); cudaDevices.push_back(2); cudaDevices.push_back(3); cudaDevices.push_back(4); cudaDevices.push_back(5);
    BOOST_REQUIRE_NO_THROW(t=new ResourceAllocatorTester(10,1.0,cudaDevices,1.0));
    BOOST_CHECK_EQUAL(t->getMaxSimultaneousReplicates(), 4);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(0));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 1);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 0);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 2);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][0], 0);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][0], 0);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(100));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 1);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 3);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][0], 100);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][0], 100);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(205));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 1);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 2);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 4);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[2][0], 205);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[2][0], 205);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(189));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 1);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 3);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 5);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[3][0], 189);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[3][0], 189);

    BOOST_REQUIRE_THROW(resources=t->assignReplicate(190), lm::Exception);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][0], 0);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][0], 100);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[2][0], 205);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[3][0], 189);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[4][0], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[5][0], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[6][0], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[7][0], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[8][0], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[9][0], -1);

    if (t != NULL) delete t; t = NULL;

    // Test a complex assignment policy.
    cudaDevices.clear(); cudaDevices.push_back(2); cudaDevices.push_back(3);
    BOOST_REQUIRE_NO_THROW(t=new ResourceAllocatorTester(10,2.0,cudaDevices,0.25));
    BOOST_CHECK_EQUAL(t->getMaxSimultaneousReplicates(), 5);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(0));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 2);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 0);
    BOOST_CHECK_EQUAL(resources.cpuCores[1], 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 2);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(1));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 2);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 2);
    BOOST_CHECK_EQUAL(resources.cpuCores[1], 3);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 3);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(2));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 2);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 4);
    BOOST_CHECK_EQUAL(resources.cpuCores[1], 5);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 2);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(3));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 2);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 6);
    BOOST_CHECK_EQUAL(resources.cpuCores[1], 7);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 3);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(4));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 2);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 8);
    BOOST_CHECK_EQUAL(resources.cpuCores[1], 9);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 2);

    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][0], 0);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][0], 0);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[2][0], 1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[3][0], 1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[4][0], 2);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[5][0], 2);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[6][0], 3);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[7][0], 3);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[8][0], 4);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[9][0], 4);

    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][0], 0);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][1], 2);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][2], 4);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][3], -1);

    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][0], 1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][1], 3);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][2], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][3], -1);

    if (t != NULL) delete t; t = NULL;

    // Test a complex assignment policy.
    cudaDevices.clear(); cudaDevices.push_back(2); cudaDevices.push_back(3);
    BOOST_REQUIRE_NO_THROW(t=new ResourceAllocatorTester(3,1.0/3.0,cudaDevices,0.25));
    BOOST_CHECK_EQUAL(t->getMaxSimultaneousReplicates(), 8);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(0));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 1);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 0);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 2);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(1));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 1);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 3);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(2));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 1);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 2);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 2);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(3));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 1);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 0);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 3);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(4));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 1);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 2);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(5));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 1);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 2);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 3);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(6));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 1);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 0);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 2);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(7));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 1);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 3);

    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][0], 0);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][1], 3);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][2], 6);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][0], 1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][1], 4);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][2], 7);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[2][0], 2);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[2][1], 5);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[2][2], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][0], 0);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][1], 2);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][2], 4);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][3], 6);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][0], 1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][1], 3);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][2], 5);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][3], 7);

    if (t != NULL) delete t; t = NULL;

    // Test a complex assignment policy.
    cudaDevices.clear(); cudaDevices.push_back(2); cudaDevices.push_back(3); cudaDevices.push_back(4); cudaDevices.push_back(5);
    BOOST_REQUIRE_NO_THROW(t=new ResourceAllocatorTester(13,5.0,cudaDevices,2.0));
    BOOST_CHECK_EQUAL(t->getMaxSimultaneousReplicates(), 2);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(0));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 5);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 0);
    BOOST_CHECK_EQUAL(resources.cpuCores[1], 1);
    BOOST_CHECK_EQUAL(resources.cpuCores[2], 2);
    BOOST_CHECK_EQUAL(resources.cpuCores[3], 3);
    BOOST_CHECK_EQUAL(resources.cpuCores[4], 4);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 2);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 2);
    BOOST_CHECK_EQUAL(resources.cudaDevices[1], 3);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(1));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 5);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 5);
    BOOST_CHECK_EQUAL(resources.cpuCores[1], 6);
    BOOST_CHECK_EQUAL(resources.cpuCores[2], 7);
    BOOST_CHECK_EQUAL(resources.cpuCores[3], 8);
    BOOST_CHECK_EQUAL(resources.cpuCores[4], 9);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 2);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 4);
    BOOST_CHECK_EQUAL(resources.cudaDevices[1], 5);

    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_CASE(RemoveReplicate)
{
    ResourceAllocatorTester * t = NULL;
    vector<int> cudaDevices;
    ResourceAllocator::ComputeResources resources;

    // Test a complex assignment policy.
    cudaDevices.clear(); cudaDevices.push_back(2); cudaDevices.push_back(3);
    BOOST_REQUIRE_NO_THROW(t=new ResourceAllocatorTester(10,2.0,cudaDevices,0.25));
    BOOST_CHECK_EQUAL(t->getMaxSimultaneousReplicates(), 5);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(0));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 2);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 0);
    BOOST_CHECK_EQUAL(resources.cpuCores[1], 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 2);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(1));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 2);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 2);
    BOOST_CHECK_EQUAL(resources.cpuCores[1], 3);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 3);

    BOOST_REQUIRE_NO_THROW(t->removeReplicate(0));

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(2));
    BOOST_CHECK_EQUAL(resources.cpuCores.size(), 2);
    BOOST_CHECK_EQUAL(resources.cpuCores[0], 0);
    BOOST_CHECK_EQUAL(resources.cpuCores[1], 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices.size(), 1);
    BOOST_CHECK_EQUAL(resources.cudaDevices[0], 2);

    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][0], 2);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][0], 2);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[2][0], 1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[3][0], 1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[4][0], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[5][0], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[6][0], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[7][0], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[8][0], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[9][0], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][0], 2);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][1], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][2], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][3], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][0], 1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][1], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][2], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][3], -1);

    if (t != NULL) delete t; t = NULL;

    // Test a complex assignment policy.
    cudaDevices.clear(); cudaDevices.push_back(2); cudaDevices.push_back(3);
    BOOST_REQUIRE_NO_THROW(t=new ResourceAllocatorTester(2,0.25,cudaDevices,0.25));
    BOOST_CHECK_EQUAL(t->getMaxSimultaneousReplicates(), 8);

    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(0));
    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(1));
    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(2));
    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(3));
    BOOST_REQUIRE_NO_THROW(t->removeReplicate(0));
    BOOST_REQUIRE_NO_THROW(t->removeReplicate(1));
    BOOST_REQUIRE_NO_THROW(t->removeReplicate(2));
    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(10));
    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(11));
    BOOST_REQUIRE_NO_THROW(resources=t->assignReplicate(12));
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][0], 10);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][1], 12);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][2], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][3], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][0], 11);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][1], 3);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][2], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][3], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][0], 10);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][1], 12);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][2], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][3], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][0], 11);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][1], 3);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][2], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][3], -1);

    BOOST_REQUIRE_NO_THROW(t->removeReplicate(10));
    BOOST_REQUIRE_NO_THROW(t->removeReplicate(11));
    BOOST_REQUIRE_NO_THROW(t->removeReplicate(12));
    BOOST_REQUIRE_NO_THROW(t->removeReplicate(3));
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][0], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][1], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][2], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[0][3], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][0], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][1], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][2], -1);
    BOOST_CHECK_EQUAL((*(t->cpuSlotsP))[1][3], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][0], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][1], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][2], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[0][3], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][0], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][1], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][2], -1);
    BOOST_CHECK_EQUAL((*(t->cudaSlotsP))[1][3], -1);

    if (t != NULL) delete t; t = NULL;
}

BOOST_AUTO_TEST_SUITE_END()
