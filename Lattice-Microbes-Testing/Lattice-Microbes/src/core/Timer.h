/*
 * University of Illinois Open Source License
 * Copyright 2016-2018 Luthey-Schulten Group,
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
 * Author(s): Tyler M. Earnest
 */

#ifndef TIMER_HH_
#define TIMER_HH_

#ifdef LINUX
#include <ctime>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

class Timer {
    timespec now,then;
    void get_time(timespec* t) { 
#ifdef __MACH__
      clock_serv_t cclock;
      mach_timespec_t mts;
      host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
      clock_get_time(cclock, &mts);
      mach_port_deallocate(mach_task_self(), cclock);
      t->tv_sec  = mts.tv_sec;
      t->tv_nsec = mts.tv_nsec;
#else
      clock_gettime(CLOCK_MONOTONIC, t); 
#endif
    }

public:

    double
    tock()
    {
      long long int sec,nsec;

      get_time(&now);

      if ((now.tv_nsec-then.tv_nsec)<0) 
        {
          sec = now.tv_sec-then.tv_sec-1;
          nsec = 1000000000L+now.tv_nsec-then.tv_nsec;
        } 
      else 
        {
          sec = now.tv_sec-then.tv_sec;
          nsec = now.tv_nsec-then.tv_nsec;
        }
      return sec + nsec*1e-9;
    }

    void tick() { get_time(&then); }

    Timer() { tick(); }
};
#endif /* LINUX */

#ifdef MACOSX
/* https://developer.apple.com/library/mac/qa/qa1398/_index.html */
#include <mach/mach_time.h>

class Timer {
	uint64_t now,then;
	mach_timebase_info_data_t base;
    void get_time(uint64_t* t) { *t=mach_absolute_time(); } 

public:

	Timer()
	{
		mach_timebase_info(&base);
		tick();
	}

    double
    tock()
    {
      long long int nsec;

      get_time(&now);
		nsec = (now - then) * base.numer / base.denom;
      return nsec*1e-9;
    }

    void tick() { get_time(&then); }
};
#endif /* MACOSX */

#endif /* TIMER_HH_ */
