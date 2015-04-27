/***************************************************************************
 *   Copyright (C) 2006-2008 by Florian Goth   *
 *   CaptainSifff@gmx.de   *
 *                                                                         *
 *   Permission is hereby granted, free of charge, to any person obtaining *
 *   a copy of this software and associated documentation files (the       *
 *   "Software"), to deal in the Software without restriction, including   *
 *   without limitation the rights to use, copy, modify, merge, publish,   *
 *   distribute, sublicense, and/or sell copies of the Software, and to    *
 *   permit persons to whom the Software is furnished to do so, subject to *
 *   the following conditions:                                             *
 *                                                                         *
 *   The above copyright notice and this permission notice shall be        *
 *   included in all copies or substantial portions of the Software.       *
 *                                                                         *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,       *
 *   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF    *
 *   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. *
 *   IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR     *
 *   OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, *
 *   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR *
 *   OTHER DEALINGS IN THE SOFTWARE.                                       *
 ***************************************************************************/
#ifndef MATRIX_H
#define MATRIX_H

#include "MTL_ICCL.h"
#include "BoundsChecker.h"
#include "MatrixGenerator.h"
#include "TypePromoter.h"
#include "Generators.h"
#include "Multiply.h"
#include "Addition.h"
#include "Determinant.h"
#include "Inversion.h"
#include "eigensystems.h"
#include <algorithm>
#include <cmath>

namespace MTLICCL
{

template <typename T> inline double scalar_abs(T& a)
{
    return std::abs(a);
}

template < class A, class B>
class Sizer
{
public:
    enum
    {
        Acolumns = A::Config::columns,
        Arows =  A::Config::rows,
        Bcolumns = B::Config::columns,
        Brows =  B::Config::rows,
    };
};

template < class A, class B >
struct TranspositionSizer
{
    enum
    {
        Rows = A::Config::columns,
        Columns = A::Config::rows
    };
};

template<class A>
inline Matrix< typename ConfigParser<A, A, TranspositionSizer >::RET > operator~ ( const A &MTL_RESTRICT lhs ) throw()
{
    //Transponierung
    typedef typename A::Elementtype ScalarType;
    typedef Matrix< typename ConfigParser<A, A, TranspositionSizer >::RET > ReturnType;
    MatrixCreator<ReturnType, None<ScalarType> > M;
    unsigned int NrOfRows = lhs.Rows();
    unsigned int NrOfColumns = lhs.Columns();
    ReturnType RetVal = M.Create ( NrOfColumns , NrOfRows , None<ScalarType> ( 0 ) );//Columns and Rows must be swapped
    for ( unsigned int i ( 0 )  ; i < NrOfColumns ; ++i )
        for ( unsigned int j ( 0 ) ; j < NrOfRows ; ++j )
            RetVal (i, j , lhs.GetElement ( j,i ) );
    return RetVal;
}

/**
This calculates the LU Decomposition of a matrix. See "Numerical Recipes"
@param m the Matrix to decompose
@param indx an array of numbers to store the permutations
@return  the sign of the permutation
*/
template <typename M>
inline int ludecompose ( MTLICCL::Matrix<M> &MTL_RESTRICT m, unsigned int *MTL_RESTRICT indx )
{
    typedef typename MTLICCL::Matrix<M>::Elementtype elementtype;
    int d = 1;//stores the sign of the permutation
    unsigned int size = m.Rows();
    double *MTL_RESTRICT vv = new double[size];
    double stemp;
    for ( unsigned int i = 0; i < size; ++i )
    {
        double big = 0.0;
        for ( unsigned int j = 0; j < size; ++j )
            if ( (stemp = scalar_abs( m( i,j ) )) > big ) big = stemp;
        if ( unlikely(big == 0.0) ) throw ( "Singular Matrix in LU Decomposition" );
        vv[i] = 1.0 / big; //save the scaling
    }
    for ( unsigned int j = 0; j < size; ++j )
    {
        for ( unsigned int i = 0; i < j; ++i )
        {
            elementtype sum = m ( i,j );
            for ( unsigned int k = 0; k < i; ++k )
                sum -= m.GetElement( i, k ) * m.GetElement( k, j );
            m ( i, j, sum );
        }
        elementtype big = 0.0;
	unsigned int imax;
        for ( unsigned int i = j; i < size; ++i )
        {
            elementtype sum = m ( i,j );
            for ( unsigned int k = 0; k < j; ++k )
                sum -= m ( i,k ) * m ( k,j );
            m ( i,j, sum );
            elementtype dum;
            if ( scalar_abs( dum = vv[i] * scalar_abs ( sum ) ) >= scalar_abs(big) ) //Is the figure of merit for the pivot better than the best so far
            {
                big = dum;
                imax = i;
            }
        }
        if ( j != imax ) //Do we need to interchange rows?
        {
            for ( unsigned int k = 0; k < size; ++k ) //yes, do it
            {
                //the actual swapping
                elementtype dum = m ( imax, k );
                m ( imax, k, m ( j,k ) );
                m ( j, k, dum );
            }
            d = -d;//change the sign of the permutation
            vv[imax] = vv[j];
        }
        indx[j] = imax;
//the next line is something NR specific
//                if(m(j, j) == 0.0) m(j,j) = 1E-06
        if ( j != (size-1) )////////////FIXME
        {
            elementtype dum = elementtype(1.0)/ m ( j, j );
            for ( unsigned int i = j + 1; i < size; ++i ) m ( i, j ) *= dum;
        }
    }
    delete [] vv;
    return d;
}

/**
A function for solving linear Systems with righthandside b via the LU Decomposition given by ludecompose.
Again: the vector b gets overwritten with the solution!
See "Numerical Recipes"
@param m the linear System
@param idx the vector with the Permutations  as given by ludecompose
@param b the rhs of the equation. it is overwritten with the solution. NOTE THAT THIS SHOULD BE A VECTOR! THUS A MATRIX WITH DIMENSION (N,1) !!!!!!!!!!!!!!!!
*/
template <typename M, typename T>
inline void lubacksubstitute( const MTLICCL::Matrix<M> &MTL_RESTRICT m, const unsigned int *const MTL_RESTRICT idx, T *const MTL_RESTRICT b) throw()
{
    const int size = m.Columns();
    int ii = -1;
    T sum;
    for(int i = 0; i < size; i++)
    {
	    int ip = idx[i];
	    sum = b[ip];
	    b[ip] = b[i];
	    if(ii >= 0)
		    for(int j = ii; j <=(i-1); j++) sum -= m(i,j) *b[j];
	    else if(sum != T(0)) ii = i;
	    b[i] = sum;
    }
    for(int i = (size -1); i >= 0; i--)
    {
	    sum = b[i];
	    for(int j = i+1; j < size; j++) sum -= m(i, j) * b[j];
	    b[i] = sum / m(i, i);
    }
	    return;
}

}

template <typename M, typename T>
inline void lubacksubstitute( const MTLICCL::Matrix<M>&MTL_RESTRICT m, const unsigned int *const MTL_RESTRICT idx, MTLICCL::Matrix<T> &MTL_RESTRICT b) throw()
{
	const unsigned int size = m.Columns();
	int ii = -1;
	float sum;
	for(int i = 0; i < static_cast<int>(size); i++)
	{
		int ip = idx[i];
		sum = b(ip,0);
		b(ip,0) = b(i,0);
		if(ii >= 0)
			for(int j = ii; j <=(i-1); j++) sum -= m(i,j) *b(j,0);
		else if(sum != 0) ii = i;
		b(i,0) = sum;
	}
	for(int i = (size -1); i >= 0; i--)
	{
		sum = b(i,0);
		for(int j = i+1; j < static_cast<int>(size); j++) sum -= m(i, j) * b(j,0);
		b(i,0) = sum / m(i, i);
	}
	return;
}

/**
A function for swapping the Rows of a Matrix
@param m the Matrix whose rows you want to swap
@param r1 the first row
@param r2 the second row
*/
template <typename M>
inline void swapRows(MTLICCL::Matrix<M> &MTL_RESTRICT m, const unsigned int r1, const unsigned int r2)
{
for(unsigned int k = 0; k < m.Columns(); ++k)
std::swap(m(r1, k), m(r2,k));
}

/**
A function for swapping the columns of a matrix
@param m the matrix whose columns you want to swap
@param c1 the first column
@param c2 the second column
*/
template <typename M>
inline void swapColumns(MTLICCL::Matrix<M> &MTL_RESTRICT m, const unsigned int c1, const unsigned int c2)
{
for(unsigned int k = 0; k < m.Columns(); ++k)
std::swap(m(k, c1), m(k,c2));
}
#endif
