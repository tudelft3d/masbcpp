/***********************************************************************
MatrixHelperFunctions - Helper functions dealing with n x m matrices.
Copyright (c) 2004-2005 Oliver Kreylos

This file is part of the Templatized Geometry Library (TGL).

The Templatized Geometry Library is free software; you can redistribute
it and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of the
License, or (at your option) any later version.

The Templatized Geometry Library is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with the Templatized Geometry Library; if not, write to the Free
Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
02111-1307 USA
***********************************************************************/

#ifndef GEOMETRY_MATRIXHELPERFUNCTIONS_INCLUDED
#define GEOMETRY_MATRIXHELPERFUNCTIONS_INCLUDED

#include <Misc/Utility.h>

namespace Geometry {

/*******************************************
Calculate subdeterminants of a 3 x m matrix:
*******************************************/

template <class ScalarParam,int numColumnsParam>
inline
void
calcSubdeterminants(
	const Matrix<ScalarParam,3,numColumnsParam>& m,
	double subdets[3][3])
	{
	/* Calculate all subdeterminants: */
	subdets[0][0]=double(m(1,1))*double(m(2,2))-double(m(2,1))*double(m(1,2));
	subdets[0][1]=double(m(1,2))*double(m(2,0))-double(m(2,2))*double(m(1,0));
	subdets[0][2]=double(m(1,0))*double(m(2,1))-double(m(2,0))*double(m(1,1));
	subdets[1][0]=double(m(2,1))*double(m(0,2))-double(m(0,1))*double(m(2,2));
	subdets[1][1]=double(m(2,2))*double(m(0,0))-double(m(0,2))*double(m(2,0));
	subdets[1][2]=double(m(2,0))*double(m(0,1))-double(m(0,0))*double(m(2,1));
	subdets[2][0]=double(m(0,1))*double(m(1,2))-double(m(1,1))*double(m(0,2));
	subdets[2][1]=double(m(0,2))*double(m(1,0))-double(m(1,2))*double(m(0,0));
	subdets[2][2]=double(m(0,0))*double(m(1,1))-double(m(1,0))*double(m(0,1));
	}

/********************************************************************
Perform Gaussian elimination with column pivoting on an n x m matrix:
********************************************************************/

template <int numRowsParam,int numColumnsParam>
inline
void
gaussElimination(
	double matrix[numRowsParam][numColumnsParam])
	{
	/* Perform Gaussian elimination with column pivoting on the extended matrix: */
	for(int step=0;step<numRowsParam-1;++step)
		{
		/* Find the column pivot: */
		double pivot=Math::abs(matrix[step][step]);
		int pivotRow=step;
		for(int i=step+1;i<numRowsParam;++i)
			{
			double val=Math::abs(matrix[i][step]);
			if(pivot<val)
				{
				pivot=val;
				pivotRow=i;
				}
			}
		
		/* Swap current and pivot rows if necessary: */
		if(pivotRow!=step)
			{
			/* Swap rows step and pivotRow: */
			for(int j=step;j<numColumnsParam;++j)
				Misc::swap(matrix[step][j],matrix[pivotRow][j]);
			}
		
		/* Combine all rows with the current row: */
		for(int i=step+1;i<numRowsParam;++i)
			{
			/* Combine rows i and step: */
			double factor=-matrix[i][step]/matrix[step][step];
			for(int j=step+1;j<numColumnsParam;++j)
				matrix[i][j]+=matrix[step][j]*factor;
			}
		}
	}

}

#endif
