/*
 * University of Illinois Open Source License
 * Copyright 2018-2018 Luthey-Schulten Group,
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

#ifndef _ARBITRARYH5_H
#define _ARBITRARYH5_H

#include <hdf5.h>
#include <algorithm>
#include <vector>
#include <cstddef>
#include <string>
#include <functional>
#include <numeric>
#include "core/Exceptions.h"
#include "core/Print.h"

#define MAX_H5_NAME  64
#define MAX_H5L_PATH  64

struct H5Lookup {
    enum Mode {
        DATASET,
        ATTR
    };
    unsigned int type;           // expected by LocalDataOutputWorker::run()
    unsigned int replicate;      // expected by LocalDataOutputWorker::run()
    size_t messageSize;          // expected by LocalDataOutputWorker::run()
    Mode mode;                   // "message" data
    char path[MAX_H5L_PATH];     // "message" data
    char attr[MAX_H5L_PATH];     // "message" data
    hid_t h5type;                // "message" data
    size_t payloadSize;          // expected by LocalDataOutputWorker::run()
};


struct H5MetaData {
    enum Mode {
        APPEND_TO_GROUP,
        APPEND_TO_DATASET,
        NEW_DATASET
    };
    unsigned int type;           // expected by LocalDataOutputWorker::run()
    unsigned int replicate;      // expected by LocalDataOutputWorker::run()
    size_t messageSize;          // expected by LocalDataOutputWorker::run()
    Mode mode;                   // "message" data
    hid_t h5type;                // "message" data
    size_t ndim;                 // "message" data
    size_t shape[H5S_MAX_RANK];  // "message" data
    char name[MAX_H5_NAME];      // "message" data
    size_t payloadSize;          // expected by LocalDataOutputWorker::run()
};

template <typename T> hid_t get_h5_type_id();

template <> inline hid_t get_h5_type_id<char>() { return H5T_NATIVE_CHAR; }
template <> inline hid_t get_h5_type_id<unsigned char>() { return H5T_NATIVE_UCHAR; }
template <> inline hid_t get_h5_type_id<short>() { return H5T_NATIVE_SHORT; }
template <> inline hid_t get_h5_type_id<unsigned short>() { return H5T_NATIVE_USHORT; }
template <> inline hid_t get_h5_type_id<int>() { return H5T_NATIVE_INT; }
template <> inline hid_t get_h5_type_id<unsigned int>() { return H5T_NATIVE_UINT; }
template <> inline hid_t get_h5_type_id<long>() { return H5T_NATIVE_LONG; }
template <> inline hid_t get_h5_type_id<unsigned long>() { return H5T_NATIVE_ULONG; }
template <> inline hid_t get_h5_type_id<long long>() { return H5T_NATIVE_LLONG; }
template <> inline hid_t get_h5_type_id<unsigned long long>() { return H5T_NATIVE_ULLONG; }
template <> inline hid_t get_h5_type_id<float>() { return H5T_NATIVE_FLOAT; }
template <> inline hid_t get_h5_type_id<double>() { return H5T_NATIVE_DOUBLE; }
template <> inline hid_t get_h5_type_id<long double>() { return H5T_NATIVE_LDOUBLE; }

template <typename T>
H5MetaData
make_H5_meta(H5MetaData::Mode mode, const std::vector<size_t> shape, const std::string name)
{
    if (shape.size() > H5S_MAX_RANK) {
        throw lm::Exception("Dimension exceeds HDF5 limits");
    }

    if (name.size() > MAX_H5_NAME-1) {
        throw lm::Exception("Name length too long");
    }

    H5MetaData md;
    md.h5type = get_h5_type_id<T>();
    md.ndim = shape.size();
    md.payloadSize = sizeof(T)*std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<size_t>());
    md.messageSize = offsetof(H5MetaData, payloadSize) - offsetof(H5MetaData, mode);

    md.mode = mode;
    std::copy(shape.begin(), shape.end(), md.shape);
    std::copy(name.begin(), name.end(), md.name);
    md.name[name.size()] = 0;
    return md;
}


inline H5Lookup
make_H5_lookup(H5Lookup::Mode mode, const std::string path, const std::string attr="")
{
    if (path.size() > MAX_H5L_PATH-1) {
        throw lm::Exception("Path length too long");
    }
    if (attr.size() > MAX_H5L_PATH-1) {
        throw lm::Exception("Attr length too long");
    }

    H5Lookup md;
    md.payloadSize = 0;
    md.messageSize = offsetof(H5Lookup, payloadSize) - offsetof(H5Lookup, mode);

    md.mode = mode;

    std::copy(attr.begin(), attr.end(), md.attr);
    md.attr[attr.size()] = 0;

    std::copy(path.begin(), path.end(), md.path);
    md.path[path.size()] = 0;
    return md;
}



#undef MAX_H5_NAME
#undef MAX_H5L_PATH

#endif /* _ARBITRARYH5_H */
