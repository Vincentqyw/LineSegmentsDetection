/*
 * Copyright (C) 2008 Leandro A. F. Fernandes and Manuel M. Oliveira
 *
 * author   : Fernandes, Leandro A. F.
 * e-mail   : laffernandes@gmail.com
 * home page: http://www.inf.ufrgs.br/~laffernandes
 *
 *
 * The complete description of the implemented techinique can be found at
 *
 *      Leandro A. F. Fernandes, Manuel M. Oliveira
 *      Real-time line detection through an improved Hough transform voting scheme
 *      Pattern Recognition (PR), Elsevier, 41:1, 2008, 299-314.
 *      DOI: http://dx.doi.org/10.1016/j.patcog.2007.04.003
 *      Project Page: http://www.inf.ufrgs.br/~laffernandes/kht.html
 *
 * If you use this implementation, please reference the above paper.
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 */

#ifndef _TYPES_
#define _TYPES_

#include <cmath>
#include <cstdlib>
#include <memory.h>
#include "buffer_2d.h"

// A simple list implementation (use it only with aggregate types).
template<typename item_type, size_t capacity_inc>
class list
{
private:

	// Specifies the size of allocated storage for the container.
	size_t m_capacity;

	// Specifies the list of items.
	item_type *m_items;

	// Counts the number of elements.
	size_t m_size;

public:

	// Erases the elements of the list.
	inline
	void clear()
	{
		m_size = 0;
	}

	// Tests if the list is empty.
	inline
	bool empty() const
	{
		return (m_size == 0);
	}

	// Returns a pointer to the list of items.
	inline
	item_type* items()
	{
		return m_items;
	}
	
	// Returns a pointer to the list of items.
	inline
	const item_type* items() const
	{
		return m_items;
	}
	
	// Class constructor.
	list() :
		m_capacity(0),
		m_items(0),
		m_size(0)
	{
	}

	// Class destructor.
	~list()
	{
		free( m_items );
	};

	// Deletes the element at the end of the list.
	inline
	void pop_back()
	{
		m_size--;
	}
	
	// Adds a new last element and returns a reference to it.
	inline
	item_type& push_back()
	{
		if (m_capacity == m_size)
		{
			m_items = (item_type*)realloc( m_items, (m_capacity += capacity_inc) * sizeof( item_type ) );
			memset( &m_items[m_size], 0, capacity_inc * sizeof( item_type ) );
		}
		return m_items[m_size++];
	}

	// Specifies a new capacity for a list.
	inline
	void reserve(const size_t capacity)
	{
		if (m_capacity < capacity)
		{
			size_t first = m_capacity;
			m_items = (item_type*)realloc( m_items, (m_capacity = capacity) * sizeof( item_type ) );
			memset( &m_items[first], 0, (capacity - first) * sizeof( item_type ) );
		}

		if (m_size > capacity)
		{
			m_size = capacity;
		}
	}

	// Specifies a new size for a list.
	inline
	void resize(const size_t size)
	{
		if (m_capacity < size)
		{
			size_t first = m_capacity;
			m_items = (item_type*)realloc( m_items, (m_capacity = size) * sizeof( item_type ) );
			memset( &m_items[first], 0, (size - first) * sizeof( item_type ) );
		}
		m_size = size;
	}

	// Returns the number of elements.
	inline
	size_t size() const
	{
		return m_size;
	}

	// Returns a reference to the list element at a specified position.
	inline
	item_type& operator [] (const size_t index)
	{
		return m_items[index];
	}

	// Returns a reference to the list element at a specified position.
	inline
	const item_type& operator [] (const size_t index) const
	{
		return m_items[index];
	}
};

// A feature pixel.
struct pixel_t
{
	float x;
	float y;
};

// A cluster of approximately collinear feature pixels.
struct cluster_t
{
	const pixel_t *pixels;
	size_t size;
};

// Specifies a list of approximately collinear feature pixels.
typedef list<cluster_t,1000> clusters_list_t;

// Specifies a string of adjacent feature pixels.
typedef list<pixel_t,1000> string_t;

// Specifies a list of string of feature pixels.
typedef list<string_t,1000> strings_list_t;

typedef list<size_t,1000> size_list_t;

struct line_t
{
	int ID;
	string_t points;
	float xs,ys;
	float xe,ye;
	int dir; //0 for k<1, 1 for k >=1
	float k,b;
	float dev;
};
// Specifies a list of line
typedef list<line_t,10000> lines_list_t;

#endif // !_TYPES_
