/*	This file AlignmentAllocator.h is part of scallop.
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
 *  Created on: Nov 7, 2016
 *      Author: Copied from http://stackoverflow.com/questions/8456236/how-is-a-vectors-data-aligned
 *      		A. Linscheid
 */

#ifndef SCALLOP_AUXILLARY_ALIGNMENTALLOCATOR_H_
#define SCALLOP_AUXILLARY_ALIGNMENTALLOCATOR_H_

#include <cstdlib>
#include <malloc.h>
#include <type_traits>
#include <new>
#include <memory>

namespace scallop
{
namespace auxillary
{

namespace detail {
    void* allocate_aligned_memory(size_t align, size_t size);
    void deallocate_aligned_memory(void* ptr) noexcept;
}

template <typename T, size_t Align = 32>
class AlignmentAllocator;


template <size_t Align>
class AlignmentAllocator<void, Align>
{
public:
    typedef void*             pointer;
    typedef const void*       const_pointer;
    typedef void              value_type;

    template <class U> struct rebind { typedef AlignmentAllocator<U, Align> other; };
};


template <typename T, size_t Align>
class AlignmentAllocator
{
public:
    typedef T         value_type;
    typedef T*        pointer;
    typedef const T*  const_pointer;
    typedef T&        reference;
    typedef const T&  const_reference;
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef std::true_type propagate_on_container_move_assignment;

    template <class U>
    struct rebind { typedef AlignmentAllocator<U, Align> other; };

public:
    AlignmentAllocator() noexcept
    {}

    template <class U>
    AlignmentAllocator(const AlignmentAllocator<U, Align>&) noexcept
    {}

    size_type
    max_size() const noexcept
    { return (size_type(~0) - size_type(Align)) / sizeof(T); }

    pointer
    address(reference x) const noexcept
    { return std::addressof(x); }

    const_pointer
    address(const_reference x) const noexcept
    { return std::addressof(x); }

    pointer
    allocate(size_type n, typename AlignmentAllocator<void, Align>::const_pointer = 0)
    {
        const size_type alignment = static_cast<size_type>( Align );
        void* ptr = detail::allocate_aligned_memory(alignment , n * sizeof(T));
        if (ptr == nullptr) {
            throw std::bad_alloc();
        }

        return reinterpret_cast<pointer>(ptr);
    }

    void
    deallocate(pointer p, size_type) noexcept
    { detail::deallocate_aligned_memory(p); }

    template <class U, class ...Args>
    void
    construct(U* p, Args&&... args)
    { ::new(reinterpret_cast<void*>(p)) U(std::forward<Args>(args)...); }

    void
    destroy(pointer p)
    { p->~T(); }
};


template <typename T, size_t Align>
class AlignmentAllocator<const T, Align>
{
public:
    typedef T         value_type;
    typedef const T*  pointer;
    typedef const T*  const_pointer;
    typedef const T&  reference;
    typedef const T&  const_reference;
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef std::true_type propagate_on_container_move_assignment;

    template <class U>
    struct rebind { typedef AlignmentAllocator<U, Align> other; };

public:
    AlignmentAllocator() noexcept
    {}

    template <class U>
    AlignmentAllocator(const AlignmentAllocator<U, Align>&) noexcept
    {}

    size_type
    max_size() const noexcept
    { return (size_type(~0) - size_type(Align)) / sizeof(T); }

    const_pointer
    address(const_reference x) const noexcept
    { return std::addressof(x); }

    pointer
    allocate(size_type n, typename AlignmentAllocator<void, Align>::const_pointer = 0)
    {
        const size_type alignment = static_cast<size_type>( Align );
        void* ptr = detail::allocate_aligned_memory(alignment , n * sizeof(T));
        if (ptr == nullptr) {
            throw std::bad_alloc();
        }

        return reinterpret_cast<pointer>(ptr);
    }

    void
    deallocate(pointer p, size_type) noexcept
    { detail::deallocate_aligned_memory(p); }

    template <class U, class ...Args>
    void
    construct(U* p, Args&&... args)
    { ::new(reinterpret_cast<void*>(p)) U(std::forward<Args>(args)...); }

    void
    destroy(pointer p)
    { p->~T(); }
};

template <typename T, size_t TAlign, typename U, size_t UAlign>
inline
bool
operator== (const AlignmentAllocator<T,TAlign>&, const AlignmentAllocator<U, UAlign>&) noexcept
{ return TAlign == UAlign; }

template <typename T, size_t TAlign, typename U, size_t UAlign>
inline
bool
operator!= (const AlignmentAllocator<T,TAlign>&, const AlignmentAllocator<U, UAlign>&) noexcept
{ return TAlign != UAlign; }

} /*namespace auxillary */
} /*namespace scallop */

#endif /* SCALLOP_AUXILLARY_ALIGNMENTALLOCATOR_H_ */
