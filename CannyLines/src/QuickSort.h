/* --- --- ---
 * Copyright (C) 2008--2010 Idiap Research Institute (.....@idiap.ch)
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
// QuickSort.h: interface for the CQuickSort class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_QUICK_SORT_H_)
#define _QUICK_SORT_H_

#include "stdio.h"

template <class TD, class TI>		/* class TD - the type of sorted data, class TI - the type of indices */
class QuickSort
{
public:
	/* sort everything inbetween `low' <-> `high' */
	static void Sort(TD *pData, long low, long high, bool bAscent=true, TI* pIdxes=NULL)
	{
		long i = low;
		long j = high;
		TD y = 0;
		TI idx;

		/* compare value */
		TD z = pData[(low + high) / 2];

		/* partition */
		do {
			if ( bAscent ) {
				/* find member above ... */
				while(pData[i] < z) i++;

				/* find element below ... */
				while(pData[j] > z) j--;
			}
			else {
				/* find member below ... */
				while(pData[i] > z) i++;

				/* find element above ... */
				while(pData[j] < z) j--;
			}

			if(i <= j) {
				/* swap two elements */
				y = pData[i];
				pData[i] = pData[j];
				pData[j] = y;

				if ( pIdxes ) {
					idx = (pIdxes)[i];
					(pIdxes)[i] = (pIdxes)[j];
					(pIdxes)[j] = idx;
				}

				i++;
				j--;
			}
		} while(i <= j);

		/* recurse */
		if(low < j)
			Sort(pData, low, j, bAscent, pIdxes);

		if(i < high)
			Sort(pData, i, high, bAscent, pIdxes);
	};

	/* sort everything inbetween `low' <-> `high' */
	static void SortAscent(TD *pData, long low, long high)
	{
		long i = low;
		long j = high;
		TD y = 0;

		/* compare value */
		TD z = pData[(low + high) / 2];

		/* partition */
		do {
			/* find member above ... */
			while(pData[i] < z) i++;
			
			/* find element below ... */
			while(pData[j] > z) j--;

			if(i <= j) {
				/* swap two elements */
				y = pData[i];
				pData[i] = pData[j];
				pData[j] = y;

				i++;
				j--;
			}
		} while(i <= j);

		/* recurse */
		if(low < j)
			SortAscent(pData, low, j);

		if(i < high)
			SortAscent(pData, i, high);
	};

	/* sort everything inbetween `low' <-> `high' */
	static void SortDescent(TD *pData, long low, long high)
	{
		long i = low;
		long j = high;
		TD y = 0;

		/* compare value */
		TD z = pData[(low + high) / 2];

		/* partition */
		do {
			/* find member below ... */
			while(pData[i] > z) i++;

			/* find element above ... */
			while(pData[j] < z) j--;

			if(i <= j) {
				/* swap two elements */
				y = pData[i];
				pData[i] = pData[j];
				pData[j] = y;

				i++;
				j--;
			}
		} while(i <= j);

		/* recurse */
		if(low < j)
			SortDescent(pData, low, j);

		if(i < high)
			SortDescent(pData, i, high);
	};

	/* sort everything inbetween `low' <-> `high' */
	static void SortAscent(TD *pData, long low, long high, TI* pIdxes)
	{
		long i = low;
		long j = high;
		TD y = 0;
		TI idx;

		/* compare value */
		TD z = pData[(low + high) / 2];

		/* partition */
		do {
			/* find member above ... */
			while(pData[i] < z) i++;
			
			/* find element below ... */
			while(pData[j] > z) j--;

			if(i <= j) {
				/* swap two elements */
				y = pData[i];
				pData[i] = pData[j];
				pData[j] = y;

				idx = (pIdxes)[i];
				(pIdxes)[i] = (pIdxes)[j];
				(pIdxes)[j] = idx;

				i++;
				j--;
			}
		} while(i <= j);

		/* recurse */
		if(low < j)
			SortAscent(pData, low, j, pIdxes);

		if(i < high)
			SortAscent(pData, i, high, pIdxes);
	};

	/* sort everything inbetween `low' <-> `high' */
	static void SortDescent(TD *pData, long low, long high, TI* pIdxes)
	{
		long i = low;
		long j = high;
		TD y = 0;
		TI idx;

		/* compare value */
		TD z = pData[(low + high) / 2];

		/* partition */
		do {
			/* find member below ... */
			while(pData[i] > z) i++;

			/* find element above ... */
			while(pData[j] < z) j--;

			if(i <= j) {
				/* swap two elements */
				y = pData[i];
				pData[i] = pData[j];
				pData[j] = y;

				idx = (pIdxes)[i];
				(pIdxes)[i] = (pIdxes)[j];
				(pIdxes)[j] = idx;

				i++;
				j--;
			}
		} while(i <= j);

		/* recurse */
		if(low < j)
			SortDescent(pData, low, j, pIdxes);

		if(i < high)
			SortDescent(pData, i, high, pIdxes);
	};

	/* sort everything inbetween `low' <-> `high' */
	static void SortAscent(TD *pData, long low, long high, TI* pIdxes,TI* pID)
	{
		long i = low;
		long j = high;
		TD y = 0;
		TI idx;

		/* compare value */
		TD z = pData[(low + high) / 2];

		/* partition */
		do {
			/* find member above ... */
			while(pData[i] < z) i++;

			/* find element below ... */
			while(pData[j] > z) j--;

			if(i <= j) {
				/* swap two elements */
				y = pData[i];
				pData[i] = pData[j];
				pData[j] = y;

				idx = (pIdxes)[i];
				(pIdxes)[i] = (pIdxes)[j];
				(pIdxes)[j] = idx;

				pID[(pIdxes)[i]]=i;
				pID[(pIdxes)[j]]=j;

				i++;
				j--;
			}
		} while(i <= j);

		/* recurse */
		if(low < j)
			SortAscent(pData, low, j, pIdxes,pID);

		if(i < high)
			SortAscent(pData, i, high, pIdxes,pID);
	};

	/* sort everything inbetween `low' <-> `high' */
	static void SortDescent(TD *pData, long low, long high, TI* pIdxes,TI* pID)
	{
		long i = low;
		long j = high;
		TD y = 0;
		TI idx;

		/* compare value */
		TD z = pData[(low + high) / 2];

		/* partition */
		do {
			/* find member below ... */
			while(pData[i] > z) i++;

			/* find element above ... */
			while(pData[j] < z) j--;

			if(i <= j) {
				/* swap two elements */
				y = pData[i];
				pData[i] = pData[j];
				pData[j] = y;

				idx = (pIdxes)[i];
				(pIdxes)[i] = (pIdxes)[j];
				(pIdxes)[j] = idx;

				pID[(pIdxes)[i]]=i;
				pID[(pIdxes)[j]]=j;

				i++;
				j--;
			}
		} while(i <= j);

		/* recurse */
		if(low < j)
			SortDescent(pData, low, j, pIdxes,pID);

		if(i < high)
			SortDescent(pData, i, high, pIdxes,pID);
	};
};

#endif // !defined(_QUICK_SORT_H_)
