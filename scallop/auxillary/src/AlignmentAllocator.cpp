/*	This file AlignmentAllocator.cpp is part of scallop.
 *
 *  scallop is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  scallop is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with scallop.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Created on: Nov 8, 2016
 *      Author: A. Linscheid
 */

#include "scallop/auxillary/AlignmentAllocator.h"
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <new>
#include <errno.h>

namespace scallop
{
namespace auxillary
{

void*
detail::allocate_aligned_memory(size_t align, size_t size)
{
    assert(align >= sizeof(void*));
    //is align a power of two?
    assert( (align & (align - 1)) == 0 );

    if (size == 0) {
        return nullptr;
    }

#ifdef __GNUC__
    void* ptr = nullptr;
    int rc = posix_memalign(&ptr, align, size);

    if (rc != 0)
    {
        if ( rc == ENOMEM )
          throw std::bad_alloc{};
        return nullptr;
    }
#endif
    return ptr;
}


void
detail::deallocate_aligned_memory(void *ptr) noexcept
{
#ifdef __GNUC__
    return free(ptr);
#endif
}

} /*namespace auxillary */
} /*namespace scallop */
