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

#ifndef GEOMETRY_PCACALCULATOR_INCLUDED
#define GEOMETRY_PCACALCULATOR_INCLUDED

#include <Geometry/Point.h>
#include <Geometry/Vector.h>
#include <Geometry/Matrix.h>

namespace Geometry {

template <int dimensionParam>
class PCACalculator // Generic class for n-dimensional PCA
	{
	/* Embedded classes: */
	public:
	typedef double Scalar; // Scalar type for points, vectors, and matrices
	static const int dimension=dimensionParam; // Dimension of point space
	typedef Geometry::Point<double,dimensionParam> Point; // Point type
	typedef Geometry::Vector<double,dimensionParam> Vector; // Vector type
	typedef Geometry::Matrix<double,dimensionParam,dimensionParam> Matrix; // Matrix type
	};

template <>
class PCACalculator<2> // Class for two-dimensional PCA
	{
	/* Embedded classes: */
	public:
	typedef double Scalar; // Scalar type for points, vectors, and matrices
	static const int dimension=2; // Dimension of point space
	typedef Geometry::Point<double,2> Point; // Point type
	typedef Geometry::Vector<double,2> Vector; // Vector type
	typedef Geometry::Matrix<double,2,2> Matrix; // Matrix type
	
	/* Elements: */
	private:
	double pxpxs,pxpys,pypys,pxs,pys; // Accumulated components of covariance matrix and centroid
	size_t numPoints; // Number of accumulated points
	Matrix cov; // The covariance matrix of all accumulated points
	
	/* Constructors and destructors: */
	public:
	PCACalculator(void)
		:pxpxs(0.0),pxpys(0.0),pypys(0.0),pxs(0.0),pys(0.0),
		 numPoints(0)
		{
		}
	
	/* Methods: */
	template <class PointParam>
	void accumulatePoint(const PointParam& point) // Accumulates the given point into the covariance matrix
		{
		/* Accumulate the point: */
		pxpxs+=double(point[0])*double(point[0]);
		pxpys+=double(point[0])*double(point[1]);
		pypys+=double(point[1])*double(point[1]);
		pxs+=double(point[0]);
		pys+=double(point[1]);
		++numPoints;
		}
	size_t getNumPoints(void) const // Returns the number of accumulated points
		{
		return numPoints;
		}
	void merge(const PCACalculator& other) // Merges the accumulated covariance matrix of another PCA calculator
		{
		/* Add the other's accumulated covariance matrix: */
		pxpxs+=other.pxpxs;
		pxpys+=other.pxpys;
		pypys+=other.pypys;
		pxs+=other.pxs;
		pys+=other.pys;
		numPoints+=other.numPoints;
		}
	Point calcCentroid(void) const // Returns the centroid of all accumulated points
		{
		return Point(pxs/double(numPoints),pys/double(numPoints));
		}
	void calcCovariance(void) // Calculates the covariance matrix of all accumulated points
		{
		double np=double(numPoints);
		cov(0,0)=(pxpxs-pxs*pxs/np)/np;
		cov(0,1)=(pxpys-pxs*pys/np)/np;
		cov(1,0)=cov(0,1);
		cov(1,1)=(pypys-pys*pys/np)/np;
		}
	const Matrix& getCovariance(void) const // Returns the already computed covariance matrix
		{
		return cov;
		}
	unsigned int calcEigenvalues(double eigenvalues[2]) const; // Calculates the eigenvalues of the covariance matrix in order of decreasing absolute value; returns the number of distinct real roots
	Vector calcEigenvector(double eigenvalue) const; // Returns the eigenvector of the covariance matrix for the given eigenvalue
	};

template <>
class PCACalculator<3> // Class for three-dimensional PCA
	{
	/* Embedded classes: */
	public:
	typedef double Scalar; // Scalar type for points, vectors, and matrices
	static const int dimension=3; // Dimension of point space
	typedef Geometry::Point<double,3> Point; // Point type
	typedef Geometry::Vector<double,3> Vector; // Vector type
	typedef Geometry::Matrix<double,3,3> Matrix; // Matrix type
	
	/* Elements: */
	private:
	double pxpxs,pxpys,pxpzs,pypys,pypzs,pzpzs,pxs,pys,pzs; // Accumulated components of covariance matrix and centroid
	size_t numPoints; // Number of accumulated points
	Matrix cov; // The covariance matrix of all accumulated points
	
	/* Constructors and destructors: */
	public:
	PCACalculator(void)
		:pxpxs(0.0),pxpys(0.0),pxpzs(0.0),pypys(0.0),pypzs(0.0),pzpzs(0.0),pxs(0.0),pys(0.0),pzs(0.0),
		 numPoints(0)
		{
		}
	
	/* Methods: */
	template <class PointParam>
	void accumulatePoint(const PointParam& point) // Accumulates the given point into the covariance matrix
		{
		/* Accumulate the point: */
		pxpxs+=double(point[0])*double(point[0]);
		pxpys+=double(point[0])*double(point[1]);
		pxpzs+=double(point[0])*double(point[2]);
		pypys+=double(point[1])*double(point[1]);
		pypzs+=double(point[1])*double(point[2]);
		pzpzs+=double(point[2])*double(point[2]);
		pxs+=double(point[0]);
		pys+=double(point[1]);
		pzs+=double(point[2]);
		++numPoints;
		}
	size_t getNumPoints(void) const // Returns the number of accumulated points
		{
		return numPoints;
		};
	void merge(const PCACalculator& other) // Merges the accumulated covariance matrix of another PCA calculator
		{
		/* Add the other's accumulated covariance matrix: */
		pxpxs+=other.pxpxs;
		pxpys+=other.pxpys;
		pxpzs+=other.pxpzs;
		pypys+=other.pypys;
		pypzs+=other.pypzs;
		pzpzs+=other.pzpzs;
		pxs+=other.pxs;
		pys+=other.pys;
		pzs+=other.pzs;
		numPoints+=other.numPoints;
		}
	Point calcCentroid(void) const // Returns the centroid of all accumulated points
		{
		return Point(pxs/double(numPoints),pys/double(numPoints),pzs/double(numPoints));
		}
	void calcCovariance(void) // Calculates the covariance matrix of all accumulated points
		{
		double np=double(numPoints);
		cov(0,0)=(pxpxs-pxs*pxs/np)/np;
		cov(0,1)=(pxpys-pxs*pys/np)/np;
		cov(0,2)=(pxpzs-pxs*pzs/np)/np;
		cov(1,0)=cov(0,1);
		cov(1,1)=(pypys-pys*pys/np)/np;
		cov(1,2)=(pypzs-pys*pzs/np)/np;
		cov(2,0)=cov(0,2);
		cov(2,1)=cov(1,2);
		cov(2,2)=(pzpzs-pzs*pzs/np)/np;
		}
	const Matrix& getCovariance(void) const // Returns the already computed covariance matrix
		{
		return cov;
		}
	unsigned int calcEigenvalues(double eigenvalues[3]) const; // Calculates the eigenvalues of the covariance matrix in order of decreasing absolute value; returns the number of distinct real roots
	Vector calcEigenvector(double eigenvalue) const; // Returns the eigenvector of the covariance matrix for the given eigenvalue
	};

}

#endif
