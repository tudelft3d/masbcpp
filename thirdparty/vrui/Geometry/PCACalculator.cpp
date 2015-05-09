/***********************************************************************
PCACalculator - Helper class to calculate the principal component
analysis matrix of a set of 3D points by a single traversal over the set
of points.
Copyright (c) 2009 Oliver Kreylos

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

#include <Geometry/PCACalculator.h>

#include <Misc/Utility.h>
#include <Math/Math.h>
#include <Math/Constants.h>

namespace Geometry {

namespace {

template <int dimensionParam>
inline
Geometry::Vector<double,dimensionParam>
calcEigenvector(
	const Geometry::Matrix<double,dimensionParam,dimensionParam>& cov,
	double eigenvalue)
	{
	/* Create the modified covariance matrix: */
	Geometry::Matrix<double,dimensionParam,dimensionParam> c=cov;
	for(int i=0;i<dimensionParam;++i)
		c(i,i)-=eigenvalue;
	
	/* Find the null space of the modified covariance matrix: */
	int rowIndices[dimensionParam];
	for(int i=0;i<dimensionParam;++i)
		rowIndices[i]=i;
	for(int step=0;step<dimensionParam-1;++step)
		{
		/* Find the full pivot: */
		double pivot=Math::abs(c(step,step));
		int pivotRow=step;
		int pivotCol=step;
		for(int i=step;i<dimensionParam;++i)
			for(int j=step;j<dimensionParam;++j)
				{
				double val=Math::abs(c(i,j));
				if(pivot<val)
					{
					pivot=val;
					pivotRow=i;
					pivotCol=j;
					}
				}
		
		/* Swap current and pivot rows if necessary: */
		if(pivotRow!=step)
			{
			/* Swap rows step and pivotRow: */
			for(int j=0;j<dimensionParam;++j)
				Misc::swap(c(step,j),c(pivotRow,j));
			}
		
		/* Swap current and pivot columns if necessary: */
		if(pivotCol!=step)
			{
			/* Swap columns step and pivotCol: */
			for(int i=0;i<dimensionParam;++i)
				Misc::swap(c(i,step),c(i,pivotCol));
			Misc::swap(rowIndices[step],rowIndices[pivotCol]);
			}
		
		/* Combine all rows with the current row: */
		for(int i=step+1;i<dimensionParam;++i)
			{
			/* Combine rows i and step: */
			double factor=-c(i,step)/c(step,step);
			for(int j=step+1;j<dimensionParam;++j)
				c(i,j)+=c(step,j)*factor;
			}
		}
	
	/* Calculate the swizzled result using backsubstitution: */
	double x[3];
	x[dimensionParam-1]=1.0;
	for(int i=dimensionParam-2;i>=0;--i)
		{
		x[i]=0.0;
		for(int j=i+1;j<dimensionParam;++j)
			x[i]-=c(i,j)*x[j];
		x[i]/=c(i,i);
		}
	
	/* Unswizzle and normalize the result: */
	Geometry::Vector<double,dimensionParam> result;
	for(int i=0;i<dimensionParam;++i)
		result[rowIndices[i]]=x[i];
	result.normalize();
	return result;
	}

}

/*********************************
Methods of class PCACalculator<2>:
*********************************/

//template <>
unsigned int
PCACalculator<2>::calcEigenvalues(
	double eigenvalues[2]) const
	{
	/* Calculate the coefficients of the covariance matrix' characteristic polynomial: */
	double mph=0.5*(cov(0,0)+cov(1,1));
	double q=cov(0,0)*cov(1,1)-cov(0,1)*cov(1,0);
	double det=Math::sqr(mph)-q;
	if(det>0.0)
		{
		det=Math::sqrt(det);
		eigenvalues[0]=mph-det;
		eigenvalues[1]=mph+det;
		if(Math::abs(eigenvalues[0])<Math::abs(eigenvalues[1]))
			Misc::swap(eigenvalues[0],eigenvalues[1]);
		return 2;
		}
	else if(det==0.0)
		{
		eigenvalues[0]=eigenvalues[1]=mph;
		return 1;
		}
	else
		return 0;
	}

//template <>
PCACalculator<2>::Vector
PCACalculator<2>::calcEigenvector(
	double eigenvalue) const
	{
	return Geometry::calcEigenvector(cov,eigenvalue);
	}

/*********************************
Methods of class PCACalculator<3>:
*********************************/

//template <>
unsigned int
PCACalculator<3>::calcEigenvalues(
	double eigenvalues[3]) const
	{
	/* Calculate the coefficients of the covariance matrix' characteristic polynomial: */
	double cp[3];
	cp[0]=-cov(0,0)-cov(1,1)-cov(2,2);
	cp[1]=cov(0,0)*cov(1,1)+cov(0,0)*cov(2,2)+cov(1,1)*cov(2,2)-cov(0,1)*cov(1,0)-cov(0,2)*cov(2,0)-cov(1,2)*cov(2,1);
	cp[2]=-cov(0,0)*(cov(1,1)*cov(2,2)-cov(1,2)*cov(2,1))+cov(0,1)*(cov(1,0)*cov(2,2)-cov(1,2)*cov(2,0))-cov(0,2)*(cov(1,0)*cov(2,1)-cov(1,1)*cov(2,0));
	
	/* Find all roots of the characteristic polynomial: */
	double q=(Math::sqr(cp[0])-3.0*cp[1])/9.0;
	double q3=Math::sqr(q)*q;
	double r=((2.0*Math::sqr(cp[0])-9.0*cp[1])*cp[0]+27.0*cp[2])/54.0;
	if(Math::sqr(r)<q3)
		{
		/* There are three real roots: */
		double theta=Math::acos(r/Math::sqrt(q3));
		eigenvalues[0]=-2.0*Math::sqrt(q)*Math::cos(theta/3.0)-cp[0]/3.0;
		eigenvalues[1]=-2.0*Math::sqrt(q)*Math::cos((theta+2.0*Math::Constants<double>::pi)/3.0)-cp[0]/3.0;
		eigenvalues[2]=-2.0*Math::sqrt(q)*Math::cos((theta-2.0*Math::Constants<double>::pi)/3.0)-cp[0]/3.0;
		
		/* Use Newton iteration to clean up the roots: */
		for(int i=0;i<3;++i)
			for(int j=0;j<5;++j)
				{
				double f=((eigenvalues[i]+cp[0])*eigenvalues[i]+cp[1])*eigenvalues[i]+cp[2];
				double fp=(3.0*eigenvalues[i]+2.0*cp[0])*eigenvalues[i]+cp[1];
				double s=f/fp;
				eigenvalues[i]-=s;
				}
		
		/* Sort the roots by descending absolute value: */
		if(Math::abs(eigenvalues[0])<Math::abs(eigenvalues[1]))
			Misc::swap(eigenvalues[0],eigenvalues[1]);
		if(Math::abs(eigenvalues[1])<Math::abs(eigenvalues[2]))
			Misc::swap(eigenvalues[1],eigenvalues[2]);
		if(Math::abs(eigenvalues[0])<Math::abs(eigenvalues[1]))
			Misc::swap(eigenvalues[0],eigenvalues[1]);
		
		return 3;
		}
	else
		{
		/* There is only one real root: */
		double a=Math::pow(Math::abs(r)+Math::sqrt(Math::sqr(r)-q3),1.0/3.0);
		if(r>0.0)
			a=-a;
		double b=a==0.0?0.0:q/a;
		eigenvalues[0]=a+b-cp[0]/3.0;
		
		/* Use Newton iteration to clean up the root: */
		for(int j=0;j<5;++j)
			{
			double f=((eigenvalues[0]+cp[0])*eigenvalues[0]+cp[1])*eigenvalues[0]+cp[2];
			double fp=(3.0*eigenvalues[0]+2.0*cp[0])*eigenvalues[0]+cp[1];
			double s=f/fp;
			eigenvalues[0]-=s;
			}
		
		/* Copy the eigenvalue twice: */
		eigenvalues[1]=eigenvalues[2]=eigenvalues[0];
		
		return 1;
		}
	}

//template <>
PCACalculator<3>::Vector
PCACalculator<3>::calcEigenvector(
	double eigenvalue) const
	{
	return Geometry::calcEigenvector(cov,eigenvalue);
	}

/***********************************************************************
Force instantiation of all standard PCACalculator classes and functions:
***********************************************************************/

template class PCACalculator<2>;
template class PCACalculator<3>;

}
