/*
 * University of Illinois Open Source License
 * Copyright 2010-2018 Luthey-Schulten Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 *               University of Illinois at Urbana-Champaign
 *               http://www.scs.uiuc.edu/~schulten
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
 * Author(s): Elijah Roberts, Mike Hallock
 */
#ifndef LPTF_PROF_H_
#define LPTF_PROF_H_




#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifndef PROF_ENABLE

#include "ProfileCodes.h"
// No-op definitions.
#define PROF_ALLOC
#define PROF_INIT
#define PROF_FORK_INIT
#define PROF_SET_THREAD(thread)
#define PROF_TEVENT(thread,event)
#define PROF_EVENT(event)
#define PROF_TBEGIN(thread,event)
#define PROF_BEGIN(event)
#define PROF_TEND(thread,event)
#define PROF_END(event)
#define PROF_WRITE
#define PROF_FORK_WRITE

#ifdef CUDA_VERSION
#define PROF_CUDA_TSTART(thread,stream)
#define PROF_CUDA_START(stream)
#define PROF_CUDA_TEVENT(thread,event,stream)
#define PROF_CUDA_EVENT(event,stream)
#define PROF_CUDA_TBEGIN(thread,event,stream)
#define PROF_CUDA_BEGIN(event,stream)
#define PROF_CUDA_TEND(thread,event,stream)
#define PROF_CUDA_END(event,stream)
#define PROF_CUDA_TFINISH(thread,stream)
#define PROF_CUDA_FINISH(stream)
#endif

#elif defined PROF_USE_NVTX
	#include "nvToolsExt.h"
	#include "Profile_NVTX.h"

	#define PROF_ALLOC
	#define PROF_INIT
	#define PROF_FORK_INIT
	#define PROF_SET_THREAD(thread)
	#define PROF_TEVENT(thread,event)
	#define PROF_EVENT(event)
	#define PROF_TBEGIN(thread,event)
	#define PROF_BEGIN(event) {nvtxRangePushA(profile_description[event]);}
	#define PROF_TEND(thread,event)
	#define PROF_END(event) {nvtxRangePop();}
	#define PROF_WRITE
	#define PROF_FORK_WRITE

	#define PROF_CUDA_TSTART(thread,stream)
	#define PROF_CUDA_START(stream)
	#define PROF_CUDA_TEVENT(thread,event,stream)
	#define PROF_CUDA_EVENT(event,stream)
	#define PROF_CUDA_TBEGIN(thread,event,stream)
	#define PROF_CUDA_BEGIN(event,stream)
	#define PROF_CUDA_TEND(thread,event,stream)
	#define PROF_CUDA_END(event,stream)
	#define PROF_CUDA_TFINISH(thread,stream)
	#define PROF_CUDA_FINISH(stream)

#else // LPTF

#include "ProfileCodes.h"

#ifndef PROF_MAX_THREADS
#error Must specify PROF_MAX_THREADS.
#endif

#ifndef PROF_MAX_EVENTS
#error Must specify PROF_MAX_EVENTS.
#endif

#ifndef PROF_OUT_FILE
#error Must specify PROF_OUT_FILE.
#endif

#define PROF_STRINGIFY(a) #a
#define PROF_MAKE_STR(a) PROF_STRINGIFY(a)


#if !defined(MACOSX) && !defined(LINUX)
#error Unknown profile architecture.
#endif


// Include the headers.
#include <stdio.h>
#if defined(MACOSX)
#include <mach/mach_time.h>
#include <machine/endian.h>
#include <pthread.h>
#elif defined(LINUX)
#include <byteswap.h>
#include <arpa/inet.h>
#include <time.h>
#endif


// Define some types.
#if defined(MACOSX)
#define _prof_time_to_net(x)        __DARWIN_OSSwapInt64(x)
typedef uint64_t                    _prof_time_t;
extern double                       _prof_time_mult;
extern pthread_key_t                _prof_thread_key;
#elif defined(LINUX)
#if __BYTE_ORDER == __LITTLE_ENDIAN
#define _prof_time_to_net(x)        bswap_64(x)
#elif __BYTE_ORDER == __BIG_ENDIAN
#define _prof_time_to_net(x)        x
#endif
typedef unsigned long long          _prof_time_t;
extern __thread unsigned int        _prof_thread_id;
extern __thread struct timespec    _prof_timespec;
#endif
#define _prof_event_to_net(x)       htonl(x)
typedef unsigned int                _prof_event_t;
extern _prof_time_t                 _prof_start_time;
extern unsigned int                 _prof_next_event[PROF_MAX_THREADS];
extern _prof_event_t                _prof_event_types[PROF_MAX_THREADS][PROF_MAX_EVENTS];
extern _prof_time_t                 _prof_event_times[PROF_MAX_THREADS][PROF_MAX_EVENTS];
extern unsigned int                 _prof_cuda_next_event[PROF_MAX_THREADS];
extern _prof_time_t                 _prof_cuda_start_time[PROF_MAX_THREADS];\

struct _eventList
{
    _prof_event_t evt;
    void *begin;
    void *end;
    struct _eventList *next;
};
typedef struct _eventList _eventList_s;
extern _eventList* _prof_cuda_event_list[PROF_MAX_THREADS];
extern void* _prof_cuda_ref_event[PROF_MAX_THREADS];



// Allocate storage.
#define PROF_ALLOC \
    double                      _prof_time_mult;\
    _prof_time_t                _prof_start_time;\
    unsigned int                _prof_next_event[PROF_MAX_THREADS];\
    _prof_event_t               _prof_event_types[PROF_MAX_THREADS][PROF_MAX_EVENTS];\
    _prof_time_t                _prof_event_times[PROF_MAX_THREADS][PROF_MAX_EVENTS];\
    unsigned int                _prof_cuda_next_event[PROF_MAX_THREADS];\
    _prof_time_t                _prof_cuda_start_time[PROF_MAX_THREADS];\
    _eventList_s*               _prof_cuda_event_list[PROF_MAX_THREADS];\
    void *                      _prof_cuda_ref_event[PROF_MAX_THREADS]; \
    PROF_ALLOC_ARCH
#if defined(MACOSX)
#define PROF_ALLOC_ARCH \
    pthread_key_t               _prof_thread_key;
#elif defined(LINUX)
#define PROF_ALLOC_ARCH \
    __thread unsigned int       _prof_thread_id;\
    __thread struct timespec    _prof_timespec;
#endif


// Perform any initialization.
#define PROF_INIT \
    {\
    PROF_INIT_ARCH \
    memset(_prof_next_event,0,sizeof(unsigned int)*PROF_MAX_THREADS);\
    memset(_prof_event_types,0,sizeof(_prof_event_t)*PROF_MAX_THREADS*PROF_MAX_EVENTS);\
    memset(_prof_event_times,0,sizeof(_prof_time_t)*PROF_MAX_THREADS*PROF_MAX_EVENTS);\
    memset(_prof_cuda_next_event,0,sizeof(unsigned int)*PROF_MAX_THREADS);\
    memset(_prof_cuda_start_time,0,sizeof(_prof_time_t)*PROF_MAX_THREADS);\
    memset(_prof_cuda_event_list,0,sizeof(_eventList_s*)*PROF_MAX_THREADS);\
    memset(_prof_cuda_ref_event,0,sizeof(void *)*PROF_MAX_THREADS); \
    }
#if defined(MACOSX)
#define PROF_INIT_ARCH \
    mach_timebase_info_data_t _prof_info;\
    mach_timebase_info(&_prof_info);\
    _prof_time_mult=((double)_prof_info.numer)/((double)_prof_info.denom);\
    _prof_start_time=mach_absolute_time();\
    if (pthread_key_create(&_prof_thread_key, NULL)) _prof_thread_key=0;
#elif defined(LINUX)
#define PROF_INIT_ARCH \
    _prof_start_time=PROF_GET_TIME;
#endif


// Perform any initialization for a new FORK.
#define PROF_FORK_INIT \
    {\
    memset(_prof_next_event,0,sizeof(unsigned int)*PROF_MAX_THREADS);\
    memset(_prof_event_types,0,sizeof(_prof_event_t)*PROF_MAX_THREADS*PROF_MAX_EVENTS);\
    memset(_prof_event_times,0,sizeof(_prof_time_t)*PROF_MAX_THREADS*PROF_MAX_EVENTS);\
    }


// Get the current time.
#if defined(MACOSX)
#define PROF_GET_TIME ((_prof_time_t)(((double)(mach_absolute_time()-_prof_start_time))*_prof_time_mult))
#elif defined(LINUX)
#define PROF_GET_TIME ((clock_gettime(CLOCK_MONOTONIC, &_prof_timespec)==0)?((_prof_time_t)((((_prof_time_t)_prof_timespec.tv_sec)*1000000000ULL)+_prof_timespec.tv_nsec-_prof_start_time)):((_prof_time_t)0ULL))
#endif


// Set the index of this thread.
#if defined(MACOSX)
#define PROF_SET_THREAD(thread) \
    {\
    if (_prof_thread_key!=0){\
        void* _prof_thread_id=NULL;\
        if ((_prof_thread_id=pthread_getspecific(_prof_thread_key)) == NULL){\
            if((_prof_thread_id=malloc(sizeof(unsigned int)))!=NULL){\
                pthread_setspecific(_prof_thread_key, _prof_thread_id);\
            }\
        }\
        if(_prof_thread_id!=NULL){\
            *((unsigned int *)_prof_thread_id)=thread;\
        }\
    }}
#elif defined(LINUX)
#define PROF_SET_THREAD(thread) _prof_thread_id=thread;
#endif


// Set the index of this thread.
#if defined(MACOSX)
#define PROF_GET_THREAD ((_prof_thread_key != 0 && pthread_getspecific(_prof_thread_key) != NULL)?(*((unsigned int *)pthread_getspecific(_prof_thread_key))):(0))
#elif defined(LINUX)
#define PROF_GET_THREAD (_prof_thread_id)
#endif


// Record an event.
#define PROF_TEVENT_AT_TIME(thread,event,event_time) \
    {if (thread < PROF_MAX_THREADS && _prof_next_event[thread] < PROF_MAX_EVENTS){\
        _prof_event_types[thread][_prof_next_event[thread]]=_prof_event_to_net(event);\
        _prof_event_times[thread][_prof_next_event[thread]]=_prof_time_to_net(event_time);\
        _prof_next_event[thread]++;\
    }}
#define PROF_TEVENT(thread, event) PROF_TEVENT_AT_TIME(thread,event,PROF_GET_TIME)
#define PROF_EVENT(event) PROF_TEVENT_AT_TIME(PROF_GET_THREAD,event,PROF_GET_TIME)


// Record an interval event start.
#define PROF_TBEGIN(thread,event) PROF_TEVENT(thread,event+10000)
#define PROF_BEGIN(event) PROF_EVENT(event+10000)


// Record an interval event end.
#define PROF_TEND(thread,event) PROF_TEVENT(thread,event+20000)
#define PROF_END(event) PROF_EVENT(event+20000)

#ifdef CUDA_VERSION

#if !defined(PROF_CUDA_ENABLE) || PROF_CUDA_ENABLE == 0

#define PROF_CUDA_TSTART(thread,stream)
#define PROF_CUDA_START(stream)
#define PROF_CUDA_TEVENT(thread,event,stream)
#define PROF_CUDA_EVENT(event,stream)
#define PROF_CUDA_TBEGIN(thread,event,stream)
#define PROF_CUDA_BEGIN(event,stream)
#define PROF_CUDA_TEND(thread,event,stream)
#define PROF_CUDA_END(event,stream)
#define PROF_CUDA_TFINISH(thread,stream)
#define PROF_CUDA_FINISH(stream)

#else

#define PROF_CUDA_TINIT(thread, event) \
    {if(thread < PROF_MAX_THREADS) { \
        _eventList_s *e=new _eventList_s; e->evt=event; cudaEventCreate((cudaEvent_t*)&e->begin); cudaEventCreate((cudaEvent_t*)&e->end); \
        e->next=_prof_cuda_event_list[thread]; _prof_cuda_event_list[thread]=e; \
    }}
#define PROF_CUDA_INIT(event) PROF_CUDA_TINIT(PROF_GET_THREAD,event)

#define PROF_CUDA_TSTART(thread,stream) \
    {if (thread < PROF_MAX_THREADS) \
        { \
            if(!_prof_cuda_ref_event[thread]) \
                cudaEventCreate((cudaEvent_t*)&_prof_cuda_ref_event[thread]); \
            if (cudaEventRecord((cudaEvent_t)_prof_cuda_ref_event[thread], stream)==cudaSuccess){\
                _prof_cuda_start_time[thread]=PROF_GET_TIME;\
            }\
        }\
    }
#define PROF_CUDA_START(stream) PROF_CUDA_TSTART(PROF_GET_THREAD,stream)

#define PROF_CUDA_TBEGIN(thread, event, stream) \
    { if(thread < PROF_MAX_THREADS) { \
            _eventList_s *e=_prof_cuda_event_list[thread]; \
            while(e!=NULL) \
            { \
                if(e->evt==event){cudaEventRecord((cudaEvent_t)e->begin, stream); break;} \
                e=e->next; \
            } \
            if(e==NULL) \
            { \
                PROF_CUDA_TINIT(thread, event); \
                e=_prof_cuda_event_list[thread]; \
                cudaEventRecord((cudaEvent_t)e->begin, stream); \
            } \
        } \
    }

#define PROF_CUDA_BEGIN(event,stream) PROF_CUDA_TBEGIN(PROF_GET_THREAD,event,stream)

#define PROF_CUDA_TEND(thread, event, stream) \
    { \
        if(thread < PROF_MAX_THREADS) \
        { \
            _eventList_s *e=_prof_cuda_event_list[thread]; \
            while(e!=NULL) \
            { \
                if(e->evt==event){cudaEventRecord((cudaEvent_t)e->end, stream); break;} \
                e=e->next; \
            } \
        } \
    }
#define PROF_CUDA_END(event,stream) PROF_CUDA_TEND(PROF_GET_THREAD,event,stream)

#define PROF_CUDA_TFINISH(thread,stream) \
    {if (thread < PROF_MAX_THREADS){\
        cudaEvent_t _prof_cuda_start_handle=(cudaEvent_t)_prof_cuda_ref_event[thread]; \
        _eventList_s *e=_prof_cuda_event_list[thread]; \
        while(e!=NULL) \
        { \
            float _prof_cuda_elaped_ms=0.0f;\
            if (cudaEventElapsedTime(&_prof_cuda_elaped_ms,(cudaEvent_t)_prof_cuda_start_handle,(cudaEvent_t)e->begin)==cudaSuccess){\
                PROF_TEVENT_AT_TIME(thread, (e->evt+10000), (_prof_cuda_start_time[thread]+(_prof_time_t)(1000000.0*((double)_prof_cuda_elaped_ms)))); \
            } \
            if (cudaEventElapsedTime(&_prof_cuda_elaped_ms,(cudaEvent_t)_prof_cuda_start_handle,(cudaEvent_t)e->end)==cudaSuccess){\
                PROF_TEVENT_AT_TIME(thread, (e->evt+20000), (_prof_cuda_start_time[thread]+(_prof_time_t)(1000000.0*((double)_prof_cuda_elaped_ms)))); \
            } \
            e=e->next; \
        } \
    }}

#define PROF_CUDA_FINISH(stream) PROF_CUDA_TFINISH(PROF_GET_THREAD,stream)

#endif //PROF_CUDA_ENABLE
#endif //CUDA_VERSION

// Write out the profile data.
#define PROF_WRITE \
    {\
    int _prof_buf_max=strlen(PROF_MAKE_STR(PROF_OUT_FILE))+256;\
    char* _prof_buf = (char*)malloc(_prof_buf_max+1);\
    if (_prof_buf!=NULL){\
        snprintf(_prof_buf,_prof_buf_max,"%s.%d",PROF_MAKE_STR(PROF_OUT_FILE),getpid());\
        FILE * _prof_fp=fopen(_prof_buf,"w");\
        if (_prof_fp!=NULL){\
            char _prof_magic[]="SBPT";\
            fwrite(_prof_magic, 1, strlen(_prof_magic), _prof_fp);\
            unsigned int _prof_tmp=htonl(1);\
            fwrite(&_prof_tmp, sizeof(_prof_tmp), 1, _prof_fp);\
            _prof_tmp=htonl(PROF_MAX_THREADS);\
            fwrite(&_prof_tmp, sizeof(_prof_tmp), 1, _prof_fp);\
            _prof_tmp=htonl(PROF_MAX_EVENTS);\
            fwrite(&_prof_tmp, sizeof(_prof_tmp), 1, _prof_fp);\
            fwrite(_prof_event_types, sizeof(_prof_event_t), PROF_MAX_THREADS*PROF_MAX_EVENTS, _prof_fp);\
            fwrite(_prof_event_times, sizeof(_prof_time_t), PROF_MAX_THREADS*PROF_MAX_EVENTS, _prof_fp);\
            fclose(_prof_fp);\
        }\
        free(_prof_buf);\
    }}

#define PROF_FORK_WRITE PROF_WRITE

#endif
#endif
