#ifndef PRECISION_H
#define PRECISION_H
#include <limits>

/*
 * This is an exact relation transforming {a, b} into the tuple {x,y} the property x = a + b in floating point arithmetic holds
 */
template <typename T>
inline void twosum(T a, T b, T& x, T& y)
{
x = a + b;
const T z = x - a;
y = (a - (x - z)) + ( b - z);
}

template <typename T>
inline void split(T a, T& x, T& y)
{
  const T factor = ceil(static_cast<T>(std::numeric_limits<T>::digits)/ 2.0);
  T c = factor * a;
  x = c - (c - a);
  y = a - x;
}

/*
 * This function transforms a*b = x + y mit a*b = x in fp
 */
template <typename T>
inline void twoproduct(T a, T b, T& x, T& y)
{
  x = a * b;
  T a1, a2;
  split(a, a1, a2);
  T b1, b2;
  split(b, b1, b2);
  y = a2 * b2 - ((( x - a1*b1) - a2*b1) - a1*b2);
}

template <typename T>
inline void squaretransform(T a, T& x, T& y)
{
  x = a * a;
  T a1, a2;
  split(a, a1, a2);
  y = a2*a2 - ((( x-a1*a1) - a2 * a1) - a1*a2);
}

/**
 * A function to reliably compare two floating point numbers
 * @param a
 * @param b
 * @return returns true if a and are equal to zero or differ only by max(a,b)*eps
 * */
template <typename FPType> 
bool fpequal(const FPType& a, const FPType& b)
{
  bool retval = true;
    if ((a != FPType(0)) || (b != FPType(0)))//looks mean, but is necessary that the next line has sense.
      retval = (std::fabs(a - b) < 5.0*std::max(std::fabs(a), std::fabs(b))*std::numeric_limits<FPType>::epsilon());
  return retval;
}
#endif
