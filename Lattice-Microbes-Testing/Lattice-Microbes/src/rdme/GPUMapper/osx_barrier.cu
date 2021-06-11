/*
 * University of Illinois Open Source License
 * Copyright 2012-2018 Luthey-Schulten Group,
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
 * Author(s): Mike Hallock
 */

#if defined(MACOSX)

#include <pthread.h>
#include "config.h"
#include "osx_barrier.h"

int pthread_barrier_init(pthread_barrier_t *b, void *attrs, int count)
{
	b->participants=count;
	b->current=0;
	pthread_mutex_init(&b->mutex,NULL);
	pthread_cond_init(&b->condition,NULL);
	return 0;
}

int pthread_barrier_wait(pthread_barrier_t *b)
{
	pthread_mutex_lock(&b->mutex);
	int rv=0;
	b->current++;
	if(b->current == b->participants)
	{
		b->current=0;
		pthread_cond_broadcast(&b->condition);
		rv=PTHREAD_BARRIER_SERIAL_THREAD;
	}	
	else
	{
		pthread_cond_wait(&b->condition, &b->mutex);
	}

	pthread_mutex_unlock(&b->mutex);
	return rv;
}

#endif
