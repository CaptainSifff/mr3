/***************************************************************************
 *   Copyright (C) 2007 - 2011 by Florian Goth   *
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
#ifndef MATRIX_HEADER_H
#define MATRIX_HEADER_H

#include "BoundsChecker.h"
#include "Generators.h"
#include <iostream>
#include <limits>

template <typename T>
std::ostream& operator<<(std::ostream& out, const MTLICCL::Matrix<T>& rhs )
{
    std::ios_base::fmtflags flags = out.flags();
    out.setf(std::ios::fixed,std::ios::floatfield);
    for ( unsigned int k = 0; k < rhs.Rows(); ++k)
    {
        out<<"( ";
        for ( unsigned int i = 0; i < rhs.Columns(); ++ i )
        {
            out<<rhs(k, i)<<" ";
        }
        out<<")"<<std::endl;
    }
    out.setf(flags);
    return out;
}

template <typename T>
void MathematicaPrint(std::ostream& out, const MTLICCL::Matrix<T>& rhs)
{
    out<<"{ ";
    if (rhs.Rows() > 1)
    {
        for ( int k = 0; k < (static_cast<int>(rhs.Rows()) - 1); ++ k)
        {
            out<<"{ ";
            if (rhs.Columns() > 1)
            {
                for ( int i = 0; i < (static_cast<int>(rhs.Columns()) - 1); ++ i )
                    out<<rhs(k, i)<<", ";
                out<<rhs(k, rhs.Columns() - 1);
            }
            else
            {
                out<<rhs(k, 0)<<std::endl;
            }
            out<<"},"<<std::endl;
        }
        out<<"{ ";
        if (rhs.Columns() > 1)
        {
            for ( int i = 0; i < (static_cast<int>(rhs.Columns()) - 1); ++ i )
                out<<rhs(rhs.Rows() - 1, i)<<", ";
            out<<rhs(rhs.Rows() - 1, rhs.Columns() - 1);
        }
        else
        {
            out<<rhs(rhs.Rows() - 1, 0)<<std::endl;
        }
        out<<"}"<<std::endl;
    }
    else
    {
        out<<"{ ";
        if (rhs.Columns() > 1)
        {
            for ( int i = 0; i < (static_cast<int>(rhs.Columns()) - 1); ++ i )
                out<<rhs(0, i)<<", ";
            out<<rhs(0, rhs.Columns() - 1);
        }
        else
        {
            out<<rhs(0, 0)<<std::endl;
        }
        out<<"}"<<std::endl;
    }
    out<<"}"<<std::endl;
}

namespace MTLICCL
{
template <class OptBoundsCheckedMatrix>
class Matrix : public OptBoundsCheckedMatrix
{
private:
public:
    friend std::ostream& operator<< <>(std::ostream& out, const Matrix<OptBoundsCheckedMatrix>& rhs);
    typedef typename OptBoundsCheckedMatrix::Config Config;
    typedef typename Config::ElementType Elementtype;
    typedef Elementtype value_type;//to mimic STL Container behaviour
    enum {IsMatrix = true};

    template < class Generator>
    inline Matrix(unsigned int i, unsigned int k, Generator G = None<Elementtype>(0.0)) : OptBoundsCheckedMatrix(i,k)
    {
        for (unsigned int x = 0; x < this->Rows() ; ++x)
            for (unsigned int y = 0; y < this->Columns() ; ++y )
                this->SetElement(x, y, G(x,y) );
        return;
    }
    inline Matrix(unsigned int i, unsigned int k) : OptBoundsCheckedMatrix(i,k)
    {
        for (unsigned int x = 0; x < this->Rows() ; ++x)
            for (unsigned int y = 0; y < this->Columns() ; ++y )
                this->SetElement(x, y, 0.0 );
        return;
    }

    inline Matrix() : OptBoundsCheckedMatrix()
    {}

    template < class Generator >
    inline Matrix(Generator G = None<Elementtype>(0.0)) : OptBoundsCheckedMatrix()
    {
        for (unsigned int k = 0 ; k < this->Rows() ; ++k)
            for (unsigned int i = 0 ; i < this->Columns() ; ++i )
                this->SetElement(k, i, G(k,i) );
        return;
    }

    inline const Elementtype& operator() (unsigned int i , unsigned int k) const
    {
        return OptBoundsCheckedMatrix::GetElement(i,k);
    }
    inline Elementtype& operator() (unsigned int i , unsigned int k)
    {
        return OptBoundsCheckedMatrix::operator() (i,k);
    }
    inline void operator() (unsigned int i , unsigned int k , Elementtype val)
    {
        return OptBoundsCheckedMatrix::SetElement(i,k,val);
    }
    inline const Matrix<OptBoundsCheckedMatrix>& operator=(const Matrix<OptBoundsCheckedMatrix>& rhs)
    {
        if ( this != &rhs)
        {
          OptBoundsCheckedMatrix::operator=(static_cast<const OptBoundsCheckedMatrix&>(rhs));
        }
        return *this;
    }
    template <typename T>
    inline const Matrix<OptBoundsCheckedMatrix>& operator=(const Matrix<T>& rhs)
    {
        //we don't need to check for self-assignment since these are distinct types
        OptBoundsCheckedMatrix::operator=(rhs);
        return *this;
    }
    inline Matrix(const Matrix& rhs) : OptBoundsCheckedMatrix(rhs)
    {
            if (( rhs.Rows() != this->Rows() ) && (rhs.Columns() != this->Columns() ) )
                throw("Sizes don't match!");
            for ( unsigned int k = 0 ; k < this->Rows() ; ++k )
                for ( unsigned int i = 0 ; i < this->Columns() ; ++i )
                    this->SetElement( k , i , rhs(k,i) );
//            std::cout<<"Copy Ctor of Matrix"<<endl;
    }
    inline void chop(Elementtype r = std::numeric_limits<double>::epsilon())
    {
      for(uint i = 0; i < this->Rows(); ++i)
	for(uint k = 0; k < this->Columns(); ++k)
	  if(std::abs(this->GetElement(i,k)) < r) this->SetElement(i,k, 0);
	  return;
    }
};
}
#endif
