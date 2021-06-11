/*
 * University of Illinois Open Source License
 * Copyright 2008-2018 Luthey-Schulten Group,
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

#ifndef LM_MAIN_DATAOUTPUTQUEUE
#define LM_MAIN_DATAOUTPUTQUEUE

#include <queue>
#include <list>
#include <google/protobuf/message.h>
#include "core/Exceptions.h"
#include "core/Types.h"
#include "core/Print.h"
#include "thread/Thread.h"
#include "io/ArbitraryH5.h"

using std::queue;
using lm::thread::PthreadException;



namespace lm {
namespace main {

/// @class DataOutputQueue
/// @brief A queue class that writes out data in the order it is recieved
class DataOutputQueue
{
protected:
    /// @class DataSet
    /// @brief A set of data in the form of bytes to be written to disk
    class DataSet
    {
    public:
        /// @brief Create an empty DataSet object
        DataSet() : size(0), data(NULL) {}
        DataSet(size_t size):size(size) {data=new byte[size];}
        /// @brief Create a DataSet object with a set of data
        /// @param origData A pointer to the data
        /// @param size The size of the data in bytes
        DataSet(void * origData, size_t size):size(size) {data=new byte[size]; memcpy(data, origData, size);}
        virtual ~DataSet() {delete [] data; data=NULL;}
        
        // Internal variables
        byte * data;
        size_t size;
    };

public:
    static const uint SPECIES_COUNTS        =   10;
    static const uint FIRST_PASSAGE_TIMES   =   11;
    static const uint PARAMETER_VALUES      =   12;
    static const uint BYTE_LATTICE          =   20;
    static const uint INT_LATTICE           =   21;
    static const uint SITE_LATTICE          =   22;
    static const uint ARBITRARY_H5          =   23;
    static const uint ARBITRARY_H5_READ     =   24;

public:
    /// @brief Sets the current DataOutputQueue that is active
    /// @param instance The new "current" queue
    static void setInstance(DataOutputQueue * instance);
    /// @brief Get the current DataOuputQueue tht is active
    /// @return Active output queue
    static DataOutputQueue * getInstance();

private:
    static DataOutputQueue * instance;  // Singleton like "active" output queue

public:
    /// @brief Create a DataOutputQueue
    DataOutputQueue();
    virtual ~DataOutputQueue();

    /// @brief Put some data on the queue to be output
    /// @param type Type of data that is written to the data as a metaheader
    /// @param message A message to put to protocol buffers
    /// @param payload The data to be output
    /// @param payloadSize The size of the payload in bytes
    /// @param payloadSerializer A function that converts the object/data to a stream of bytes (If this is NULL the object is written out as pure bytes)
    virtual void pushDataSet(uint type, uint replicate, ::google::protobuf::Message * message, void * payload=NULL, size_t payloadSize=0, void (*payloadSerializer)(void *, void *, size_t)=NULL);
    /// @brief Put some data on the queue 
    /// @param data The data to be output
    /// @param dataSize Size of the data in bytes
    virtual void pushDataSet(void * data, size_t dataSize);
    /// @brief Put a DataSet on the queue
    /// @param dataSet The data set to be put on the queue for output
    virtual void pushDataSet(DataSet * dataSet);

    virtual void pushDataSet(const H5MetaData md, uint replicate, void * payload);
protected:
    pthread_mutex_t dataMutex;  // A lock for multithreaded applications
    queue<DataSet *> dataQueue; // Queue for data to be added to

    pthread_mutex_t returnMutex;
    std::list<DataSet *> returnList;
    pthread_cond_t returnAvailable;

public:

    template <typename T>
    void
    queryH5(const H5Lookup::Mode mode, H5MetaData &hdr, std::vector<T> &payload, const uint replicate,
            const std::string path, const std::string attr="")
    {
        Print::printf(Print::INFO, "in queryH5");
        H5Lookup req = make_H5_lookup(mode, path, attr);
        size_t buffer_sz = sizeof(req)+req.payloadSize;
        req.h5type = get_h5_type_id<T>();
        DataSet* dataSet = new DataSet(buffer_sz);
        H5Lookup *md_ptr = reinterpret_cast<H5Lookup *>(dataSet->data);
        *md_ptr = req;
        md_ptr->type = ARBITRARY_H5_READ;
        md_ptr->replicate = replicate;
        pushDataSet(dataSet);
        payload.clear();

        while (payload.size() == 0) {
        //// BEGIN CRITICAL SECTION: returnMutex
            Print::printf(Print::INFO, "pre returnMutex");
            PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&returnMutex));
            Print::printf(Print::INFO, "pre cond_wait");
            PTHREAD_EXCEPTION_CHECK(pthread_cond_wait(&returnAvailable, &returnMutex));
            Print::printf(Print::INFO, "post cond_wait");
            std::list<DataSet *>::iterator it=returnList.begin();
            while (it != returnList.end()) {
                hdr = *reinterpret_cast<H5MetaData *>((*it)->data);
                if (hdr.replicate == replicate) {
                    T *ds = reinterpret_cast<T*>((*it)->data + sizeof(H5MetaData));
                    payload.resize(hdr.payloadSize/sizeof(T));
                    std::copy(ds, ds+hdr.payloadSize/sizeof(T), payload.begin());
                    returnList.erase(it);
                    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&returnMutex));
                    return;
                }
                it++;
            }
            PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&returnMutex));
        //// END CRITICAL SECTION: returnMutex
        }
    }

    template <typename T>
    T
    queryH5attr(const uint replicate, const std::string path, const std::string attr)
    {
        std::vector<T> payload;
        H5MetaData hdr;
        queryH5(H5Lookup::ATTR, hdr, payload, replicate, path, attr);
        if (hdr.ndim == 1 && hdr.shape[0] == 1) {
            return payload[0];
        }
        else {
            throw Exception("Got non-scalar from HDF5 attribute lookup");
        }
    }
};

}
}


#endif
