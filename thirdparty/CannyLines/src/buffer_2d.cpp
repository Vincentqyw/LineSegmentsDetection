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

#include <cstdlib>
#include <memory.h>
#include "buffer_2d.h"

// Allocates 2D memory blocks.
void*
malloc_2d(const size_t size1, const size_t size2, const size_t data_size)
{
	const size_t pointers_size = size1 * sizeof( void* );
	const size_t items_size = size1 * size2 * data_size;

	void *buffer = malloc( pointers_size + items_size );

	void **pointers = static_cast<void**>( buffer );
	char *items = &(static_cast<char*>( buffer ))[pointers_size];

	for (size_t i=0, j=0, j_inc=size2*data_size; i!=size1; ++i, j+=j_inc)
	{
		pointers[i] = &items[j];
	}

	return buffer;
}

// Sets 2D buffers to a specified character.
void*
memset_2d(void *dest, const int c, const size_t size1, const size_t size2, const size_t data_size)
{
	if (dest)
	{
		const size_t pointers_size = size1 * sizeof( void* );
		const size_t items_size = size1 * size2 * data_size;

		char *buffer = static_cast<char*>( dest );
		
		memset( &buffer[pointers_size], c, items_size );
	}

	return dest;
}

// Reallocate 2D memory blocks.
void*
realloc_2d(void *memblock, const size_t size1, const size_t size2, const size_t data_size)
{
	const size_t pointers_size = size1 * sizeof( void* );
	const size_t items_size = size1 * size2 * data_size;

	memblock = realloc( memblock, pointers_size + items_size );

	void **pointers = static_cast<void**>( memblock );
	char *items = &(static_cast<char*>( memblock ))[pointers_size];

	for (size_t i=0, j=0, j_inc=size2*data_size; i!=size1; ++i, j+=j_inc)
	{
		pointers[i] = &items[j];
	}

	return memblock;
}
