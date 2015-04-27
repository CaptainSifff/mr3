#ifndef EIGENSYSTEMS_H
#define EIGENSYSTEMS_H
#include <inttypes.h>
#include <cmath>
#include <cstring>
#include <vector>
#include <list>
#include <limits>
#include <algorithm>
#include <complex>
#include <iomanip>
#include "Matrix_Header.h"
#include "Matrix.h"
#include "precision.h"

constexpr uint CLS = 64;

namespace MTLICCL
{

template <typename T>
struct Interval
{
    T xmin;
    T xmax;
    uint EVsxmin;
    uint EVsxmax;
    std::vector<Interval<T> > substructure;
    Interval(T i, T a, uint l, uint u) : xmin(i), xmax(a), EVsxmin(l), EVsxmax(u) {}
    Interval() : xmin(0), xmax(0), EVsxmin(0), EVsxmax(0) {}
    inline const Interval<T>& operator=(const Interval& rhs)
    {
        if (this != & rhs)
        {
            this->xmin = rhs.xmin;
            this->xmax = rhs.xmax;
            this->EVsxmin = rhs.EVsxmin;
            this->EVsxmax = rhs.EVsxmax;
            this->substructure = rhs.substructure;
        }
        return *this;
    }
    void print() const
    {
        std::cout<<std::scientific<<"xmin = "<<xmin<<" , xmax = "<<xmax<<" , EVs = "<<EVsxmax - EVsxmin<<std::endl;
    }
    void shiftsubstructure(T shift)
    {
        for (typename std::vector<Interval<T> >::iterator it = substructure.begin(); it != substructure.end(); ++it)
        {
            it->xmin -= shift;
            it->xmax -= shift;
        }
    }
};

template <typename T>
class TriDiagMatrix
{
public:
    TriDiagMatrix() : diag(NULL), off(NULL), size(0) {}
    T* diag;
    T* off;
    uint size;
private:
};

template <typename T>
class TwistedMatrix
{
public:
    T* diag;
    T* off;
    uint twist;
    uint len;

    TwistedMatrix() : diag(NULL), off(NULL), twist(0) {}
    TwistedMatrix(const TwistedMatrix<T>& arg) : diag(arg.diag), off(arg.off), twist(arg.twist), len(arg.len) {}
    const TwistedMatrix<T> operator=(const TwistedMatrix<T>& rhs)
    {
        if (this != &rhs)
        {
            diag = rhs.diag;
            off = rhs.off;
            twist = rhs.twist;
            len = rhs.len;
        }
        return *this;
    }
    void print()
    {
        for (int k = 0; k < len; ++k)
        {
            std::cout<<diag[k]<<" "<<off[k]<<std::endl;
        }
        return;
    }
private:
};

template <typename T>
class SolutionPair
{
public:
    SolutionPair(T ev, T *const evec, uint l) : evalue(ev), len(l)
    {
#ifndef __USE_XOPEN2K
      evector = new T[l];
#else
      int ret = posix_memalign((void**)(&evector), 16, sizeof(T)*l);
      if(unlikely(ret != 0))
      {
        throw std::bad_alloc();
      }
#endif
        memcpy(evector, evec, len*sizeof(T));
    }
    SolutionPair(T ev, uint l) : evalue(ev), len(l)
    {
#ifndef __USE_XOPEN2K
      evector = new T[l];
#else
      int ret = posix_memalign((void**)(&evector), 16, sizeof(T)*l);
      if(unlikely(ret != 0))
      {
        throw std::bad_alloc();
      }
#endif
        memset(evector, 0, l*sizeof(T));
    }
    SolutionPair(uint l) : evalue(0.0), len(l)
    {
#ifndef __USE_XOPEN2K
      evector = new T[l];
#else
      int ret = posix_memalign((void**)(&evector), 16, sizeof(T)*l);
      if(unlikely(ret != 0))
      {
        throw std::bad_alloc();
      }
#endif
        memset(evector, 0, l*sizeof(T));
    }
    ~SolutionPair()
    {
#ifndef __USE_XOPEN2K
      delete [] evector;
#else
      if(evector != NULL)
	free(evector);
#endif
    }
    SolutionPair(const SolutionPair<T>& rhs) : evalue(rhs.evalue), len(rhs.len)
    {
      #ifndef __USE_XOPEN2K
      evector = new T[rhs.len];
#else
      int ret = posix_memalign((void**)(&evector), 16, sizeof(T)*rhs.len);
      if(unlikely(ret != 0))
      {
        throw std::bad_alloc();
      }
#endif
        memcpy(evector, rhs.evector, len*sizeof(T));
    }
#ifdef HAS_RVALUE_REFERENCES
    SolutionPair(SolutionPair<T>&& rhs) : evalue(rhs.evalue), evector(rhs.evector), len(rhs.len)
    {
        rhs.evector = NULL;
    }
#endif
    void print()
    {
//    std::cout.precision(14);
        std::cout<<"Eigenvalue: "<<std::scientific<<evalue<<std::endl;
        std::cout<<"Eigenvector: (";
        for (uint k = 0; k < (len-1); ++k)
            std::cout<<k<<" "<<evector[k]<<", "<<std::endl;
        std::cout<<evector[len-1]<<")"<<std::endl;
    }
    T evalue;
    T* evector;
private:
    uint len;
};

/** The problem set contains information about the minimum and  maximum required Eigenvector. MRRR would allow an arbitrary non-continuous set
 */
template <typename T>
struct ProblemSet
{//pahole says this structure has 4bytes of padding at the end(64bit system)
    Interval<T> iv;
    uint min_EV_idx;
    uint max_EV_idx;
    T shift;
    TwistedMatrix<T> matrix;
    uint level;
};

template <typename T>
class CTDT//ComplexTriDiagonalizeTrait
{
public:
    typedef T value_type;
    CTDT(value_type* d, value_type* o, uint len) : diagarr(d), offarr(o) {}
    static inline T realpart( const T& a) {
        return a;
    }
    static inline T scale_abs(const T& a) {
        return std::abs(a);
    }
    static inline T abssquared(const T& a) {
        return a*a;
    }
    static inline T conjugate(const T& a) {
        return a;
    }
    /** This function copies over the data to the final storage and updates if necessary the transformation matrices.
    * It is a No-Op for real matrix entries.
    */
    template <class Mat>
    void copyover(Mat& mat, bool calcEVs) {}
    inline T get_g(T& f, T& h)
    {
        return (f >= 0.0) ? -std::sqrt(h) : std::sqrt(h);
    }
    T& diag(uint k) {
        return diagarr[k];
    }
    T& off(uint k) {
        return offarr[k];
    }
    T& diag(uint k) const {
        return diagarr[k];
    }
    T& off(uint k) const {
        return offarr[k];
    }
private:
    value_type* diagarr;
    value_type* offarr;
};

template <typename T>
class CTDT<std::complex<T> >
{
public:
    typedef T value_type;
    CTDT(value_type* d, value_type* o, uint l) : diagarr(d), offarr(o), len(l)
    {
#ifndef __USE_XOPEN2K
      diagcplx = new std::complex<T>[2*l];
#else
      int ret = posix_memalign((void**)(&diagcplx), 16, 2*sizeof(std::complex<T>)*l);
      if(unlikely(ret != 0))
      {
        throw std::bad_alloc();
      }
#endif
        offcplx = diagcplx + len;
    }
    static inline T realpart( const std::complex<T>& a) {
        return real(a);
    }
    static inline T scale_abs(const std::complex<T>& a) {
        return std::abs(real(a)) + std::abs(imag(a));    //that's the way the scale is determined in htridi.f
    }
    static inline T abssquared(const std::complex<T>& a) {
        return norm(a);
    }
    static inline std::complex<T> conjugate(const std::complex<T>& a) {
        return conj(a);
    }
    ~CTDT()
    {
#ifndef __USE_XOPEN2K
        delete [] diagcplx;
#else
	if(diagcplx != NULL)
	  free(diagcplx);
#endif
//        delete [] offcplx;
    }
    /** This function copies over the data to the final storage and updates if necessary the transformation matrices.
    * The Transformation matrices get their phase and the real tridiagonal representation is determined.
    */
    template <class Mat>
    void copyover(Mat& mat, bool calcEVs)
    {
        if (calcEVs)
        {
            offarr[len - 1] = std::abs(off(len - 1));
            diagarr[len -1] = real(diagcplx[len - 1]);
            std::complex<T> phase = 1.0;
            for (int j = len - 2; j >= 0; --j)
            {
                offarr[j] = std::abs(off(j));
                diagarr[j] = real(diagcplx[j]);
                if (offarr[j+1] != 0.0)
                {
                    phase *= off(j+1) / offarr[j+1];
                    for (uint i = 0; i < len; ++i)
                        mat(i, j) *= phase;//this is the multplication with the matrix that makes the final tridiagonal matrix real
                }
            }
        }
        else
        {
            for (uint i = 0; i < len; ++i)
            {
                offarr[i] = std::abs(offcplx[i]);
                diagarr[i] = real(diagcplx[i]);
            }
        }
    }
    inline std::complex<T> get_g(std::complex<T>& f, T& h)
    {
        std::complex<T> sigma = (norm(f) > 0 ? myangle(f) : 1);
        sigma *= h;
        std::complex<T> g = std::sqrt(sigma);
        std::complex<T> c1 = f+g;
        std::complex<T> c2 = f-g;
        if (norm(c1) > norm(c2))
            g = -g;
        return g;
    }
    std::complex<T>& diag(uint k) {
        return diagcplx[k];
    }
    std::complex<T>& off(uint k) {
        return offcplx[k];
    }
    std::complex<T>& diag(uint k) const {
        return diagcplx[k];
    }
    std::complex<T>& off(uint k) const {
        return offcplx[k];
    }
private:
    value_type* diagarr;
    value_type* offarr;
    std::complex<T>* diagcplx;
    std::complex<T>* offcplx;
    uint len;
    std::complex<T> myangle(std::complex<T>& z)//calculates z/conj(z)
    {
        T a = z.real();
        T b = z.imag();
        if (std::abs(a) < std::abs(b))
        {
            T r = a * (a/b);
            T n = r + b;
            return std::complex<T>(r -b, 2.0*a)/n;
        }
        T r = b * (b/a);
        T n = r + a;
        return std::complex<T>(a - r, 2.0*b)/n;
    }
};

template <class Matrix>
class EigenSystemSolver
{
public:
    typedef typename Matrix::value_type value_type;
    typedef typename CTDT<value_type>::value_type datatype;
    EigenSystemSolver(const Matrix& m, int mini = -1, int maxi = -1, bool cEVs = false) : orig(m), calcEVs(cEVs), diag(NULL), lower(NULL), min_idx((mini>-1?mini : 0)), max_idx((maxi >0? maxi : (m.Rows() - 1)))
    {
        if (unlikely(mini > maxi)) throw("invalid Interval for Eigenvalues!!");
 #ifndef __USE_XOPEN2K
     diag = new datatype[2*orig.Rows()];
#else
     int ret = posix_memalign((void**)(&diag), 16, 2*sizeof(datatype)*orig.Rows());
     if(unlikely(ret != 0))
     {
       throw std::bad_alloc();
     }
#endif
     lower = diag + orig.Rows();
    }
    ~EigenSystemSolver()
    {
#ifndef __USE_XOPEN2K
      delete [] diag;
#else
      free(diag);
#endif
    }
    void calculateEVs(bool evs = true) {
        calcEVs = evs;
    }
    Matrix tridiagonalize();
    Matrix tridiagonalize2();
    std::vector<SolutionPair<datatype> > eigensystem();
    bool testeigensystem(std::vector<SolutionPair<datatype> >& esystem, Matrix& z);
private:
#if defined(__GXX_EXPERIMENTAL_CXX0X__) && (GCC_VERSION > GCC_VER(4,5,0))
    static constexpr double gaptol = 0.005;//compromise for double from Willems paper
#else
    static const double gaptol = 0.005;//compromise for double from Willems paper
#endif
    const Matrix& orig;
    datatype* diag;
    datatype* lower;
    uint min_idx;
    uint max_idx;
    datatype pivmin;
    Matrix* orthoTrans;///< place to store the orthogonal transformation that does the similarity transform
    bool calcEVs;
    void mrrr(uint min, uint diaglen, std::vector<ProblemSet<datatype> >& nodes, std::vector<SolutionPair<datatype> >& esystem);
    /**
    @param t the tridiagonal matrix
    @param shift an eventual shift
    @param d a preallocated array for the diagonal entries
    @param l a preallocated array for the off-diagonal entries
    */
    void generatelowerBidiagonalfactorization(datatype *const d, datatype *const l, uint min, uint max, datatype shift) const;
    void generateupperBidiagonalfactorization(datatype *const r, datatype *const u);
    Interval<datatype> findinitialInterval(int, int) const;
    void bisectLDL(std::vector<Interval<datatype> >& finalintervals, ProblemSet<datatype>&) const;
    uint sturmcountLDL(const TwistedMatrix<datatype>& mat, datatype x) const;
    uint sturmcount(datatype x, datatype shift, uint min, uint max) const;
    void refineSingletonInterval(Interval<datatype>& iv, const TwistedMatrix<datatype>& m) const;
    void bisect(std::vector<Interval<datatype> >& finalintervals, Interval<datatype>& interval, datatype shift);
    SolutionPair<datatype> rqi(const TwistedMatrix<datatype>& m, Interval<datatype>& iv, uint min) const;
    /**
     * @param diag_src diagonal entries of representation
     * @param off_src other entries
     * @param diag_out new diagonal entries
     * @param off_out new off_diagonal entries
     * @param shift the shift value
     * @param twist_k the twist index of the src representation
     * @param twist_t the twist index of the new representation
     */
    void dtwqds(const datatype *const diag_src, const datatype *const off_src, datatype* diag_out, datatype* off_out, datatype shift, uint twist_k, uint twist_t, uint srclen) const;
    TwistedMatrix<datatype> initialdtwqds(const datatype *const d, const datatype *const l, const datatype shift, uint twist, uint diaglen) const;
    template <class A>
    void dstqds_part(const datatype *const diag_src, const datatype *const off_src, datatype* diag_out, datatype* off_out, datatype shift, uint a, datatype* s) const;
    template <class A>
    void flipped_dstqds_part(const datatype *const diag_src, const datatype *const off_src, datatype* diag_out, datatype* off_out, datatype shift, int a, const int len, datatype* s) const;
    template <class A>
    void dqds_part(const datatype *const diag_src, const datatype *const off_src, datatype* diag_out, datatype* off_out, datatype shift, uint lower_limit, uint upper_limit, datatype* p, datatype off) const;
    template <class A>
    void flipped_dqds_part(const datatype *const diag_src, const datatype *const off_src, datatype* diag_out, datatype* off_out, datatype shift, uint lower_limit, uint upper_limit, datatype* p, datatype off) const;
    ProblemSet<datatype> createInitialProblem(uint, uint);
};

struct PivMin
{
    template <typename T>
    static void call(T& d, const T& pivmin)
    {
        if (std::abs(d) < pivmin) d = -pivmin;
    }
};

struct NoPivMin
{
    template <typename T>
    static void call(T&, const T&)
    {}
};

/**
 * check for a NaN at s[a]
 * dstqds takes an LDL^t rep and transform it to another LDL^t rep
 */
template <class Matrix>
template <class A>
void EigenSystemSolver<Matrix>::dstqds_part(const datatype *const diag_src, const datatype *const off_src, datatype* diag_out, datatype* off_out, datatype shift, uint a, datatype* s) const
{//Algorithm 2.7 without line 7
    s[0] = 0;//the auxilliary data
    for (uint i = 0; i < a; ++i)
    {
        datatype temp = s[i] - shift;
        diag_out[i] = diag_src[i] + temp;
        A::template call<datatype>(diag_out[i], pivmin);
        //This is what is referred to as the N-representation
        {
            off_out[i] = off_src[i]*(diag_src[i]/diag_out[i]);
            s[i+1] = off_out[i] * off_src[i]*temp;
        }
    }
}

/**
 * check for a NaN at s[srclen - 1]
 */
template <class Matrix>
template <class A>
void EigenSystemSolver<Matrix>::flipped_dstqds_part(const datatype *const diag_src, const datatype *const off_src, datatype* diag_out, datatype* off_out, datatype shift, int a, const int len, datatype* s) const
{//Algorithm 2.7 without line 7
    s[a] = 0;//the auxilliary data
    for (int i = len - 1; i > a; --i)
    {
        datatype temp = s[len + a - i - 1] - shift;
        diag_out[i] = diag_src[i] + temp;
        A::template call<datatype>(diag_out[i], pivmin);
        //This is what is referred to as the N-representation
        {
            off_out[i-1] = off_src[i-1]*(diag_src[i]/diag_out[i]);
            s[len + a - i] = off_out[i-1] * off_src[i-1]*temp;
        }
    }
}

/**
 * check for NaN in p[upper_limit]
 * dqds takes an URU^t rep and transforms it to a LDL^t rep
 */
template <class Matrix>
template <class A>
void EigenSystemSolver<Matrix>::dqds_part(const datatype *const diag_src, const datatype *const off_src, datatype* diag_out, datatype* off_out, datatype shift, uint lower_limit, uint upper_limit, datatype* p, datatype off) const
{//Algorithm 2.9 without line 8, again we use the N - representation
//Note the fixed off-Diagonal-Entry indexing in contrast to Paul Willems Thesis
    datatype delta = off + (diag_src[lower_limit] - shift);
    for (uint i = lower_limit; i < upper_limit; ++i)
    {
        diag_out[i] = off_src[i] * off_src[i]*diag_src[i+1] + delta;
//	std::cout<<off_src[i] * off_src[i]*diag_src[i+1]<< " "<<lower[i+1]*lower[i+1]/diag_src[i]<<" "<<off_src[i]*diag_src[i+1]<<std::endl;
        A::template call<datatype>(diag_out[i], pivmin);
        datatype temp = (diag_src[i+1]/diag_out[i]);
        off_out[i] = temp*off_src[i];

        p[i+1] = temp * delta;
        delta = p[i+1] - shift;
    }
}

/**
 * check for NaN in p[lower_limit]
 * flipped dqds takes an LDL^t rep and transform it to a URU^t
 */
template <class Matrix>
template <class A>
void EigenSystemSolver<Matrix>::flipped_dqds_part(const datatype *const diag_src, const datatype *const off_src, datatype* diag_out, datatype* off_out, datatype shift, uint lower_limit, uint upper_limit, datatype* p, datatype off) const
{//Algorithm 2.9 without line 8, again we use the N - representation
//Note the fixed off-Diagonal-Entry indexing in contrast to Paul Willems Thesis
    datatype delta = off + (diag_src[upper_limit] - shift);
    for (int i = upper_limit; i > static_cast<int>(lower_limit); --i)
    {
        diag_out[i] = off_src[i-1] * off_src[i-1]*diag_src[i-1] + delta;
//	std::cout<<off_src[i] * off_src[i]*diag_src[i+1]<< " "<<lower[i+1]*lower[i+1]/diag_src[i]<<" "<<off_src[i]*diag_src[i+1]<<std::endl;
        A::template call<datatype>(diag_out[i], pivmin);
        datatype temp = (diag_src[i-1]/diag_out[i]);
        off_out[i-1] = temp*off_src[i-1];

        p[i-1] = temp * delta;
        delta = p[i-1] - shift;
    }
}

/**
 * check for NaN in p[lower_limit]
 * flipped dqds takes an LDL^t rep and transform it to a URU^t
 */
template <class Matrix>
void EigenSystemSolver<Matrix>::dtwqds(const datatype *const diag_src, const datatype *const off_src, datatype* diag_out, datatype* off_out, datatype shift, uint twist_k, uint twist_t, uint srclen) const
{//Algorithm 2.10
//again validated to work on 13.02.12...
    uint a = std::min(twist_k, twist_t);
    uint b = std::max(twist_k, twist_t);
    datatype* s =
    //new datatype[srclen+1];
//    memset(s, 0, (srclen+1)*sizeof(datatype));
    (datatype*)calloc(srclen + 1, sizeof(datatype));
    dstqds_part<NoPivMin>(diag_src, off_src, diag_out, off_out, shift, a, s);
    if (!std::isfinite(s[a]))
        dstqds_part<PivMin>(diag_src, off_src, diag_out, off_out, shift, a, s);
    //save last auxilliary
    datatype templ = s[a];

    flipped_dstqds_part<NoPivMin>(diag_src, off_src, diag_out, off_out, shift, b, srclen, s);
    if (!std::isfinite(s[srclen - 1]))
        flipped_dstqds_part<PivMin>(diag_src, off_src, diag_out, off_out, shift, b, srclen, s);
    datatype tempu = s[srclen - 1];
    if (twist_k < twist_t)
    {
        dqds_part<NoPivMin>(diag_src, off_src, diag_out, off_out, shift, twist_k, twist_t, s, templ );
        if (!std::isfinite(s[twist_t]))
            dqds_part<PivMin>(diag_src, off_src, diag_out, off_out, shift, twist_k, twist_t, s, templ );
        diag_out[twist_t] = (s[twist_t] + tempu) - shift;//write out gamma_t
//    std::cout<<diag_out[twist_t]<<" = "<<s[twist_t]<<" + "<<tempu<<" - "<<shift<<" k < t"<<std::endl;
    }
    else if (twist_k == twist_t)
    {
        diag_out[twist_k] = diag_src[twist_t] + ((templ + s[srclen - 1]) - shift);//write gamma_k
//     std::cout<<diag_out[twist_k]<<" = "<<diag_src[twist_t]<<" + "<<templ<<" + "<<s[orig.Rows() - 1]<<" - "<<shift<<" k == t" <<std::endl;
    }
    else
    {
        flipped_dqds_part<NoPivMin>(diag_src, off_src, /*flipped_*/diag_out, /*flipped_*/off_out, shift, twist_t, twist_k, s, tempu);
        if (!std::isfinite(s[twist_t]))
            flipped_dqds_part<PivMin>(diag_src, off_src, /*flipped_*/diag_out, /*flipped_*/off_out, shift, twist_t, twist_k, s, tempu);
        diag_out[twist_t]=(templ + s[twist_t]) - shift;
    }
    free(s);
//    delete [] s;
    return;
}

template <class Matrix>
void EigenSystemSolver<Matrix>::generatelowerBidiagonalfactorization(datatype *const d, datatype *const l, uint min, uint max, datatype shift) const
{
    d[0] = diag[min]-shift;
    for (uint i = 0; i < max - 1 - min; ++i)
    {
        l[i] = lower[min + i+1]/d[i];
        d[i+1] = (diag[min + i+1] -shift)- l[i]*l[i]* d[i];
    }
    return;
}

template <class Matrix>
void EigenSystemSolver<Matrix>::generateupperBidiagonalfactorization(datatype *const r, datatype *const u)
{
    r[orig.Rows()-1] = diag[orig.Rows()-1];
    for (int i = orig.Rows() - 1; i > 0; --i)
    {
        u[i] = lower[i]/r[i];
        r[i - 1] = diag[i-1] - u[i] * u[i] * r[i];
    }
}

template <typename Matrix>
uint EigenSystemSolver<Matrix>::sturmcount( datatype x, datatype shift, uint min, uint max) const
{
    datatype qk = diag[min] - shift - x;
    uint a = (qk < 0? 1 : 0);
    for (uint k = min + 1; k < max; ++k)
    {
        datatype sub = (std::abs(qk) - 100.0*std::numeric_limits<datatype>::epsilon() ?
                        lower[k]*(lower[k]/qk) :
                        lower[k] / std::numeric_limits<datatype>::epsilon());
        qk = (diag[k] - shift - x) - sub;
        if (qk < 0) a++;
    }
    return a;
}

template <typename T>
T relgap(T a, T b)
{
    T retval = std::abs(a-b)/std::max(std::abs(a) , std::abs(b));
    return retval;
}

template <class Matrix>
Interval<typename EigenSystemSolver<Matrix>::datatype> EigenSystemSolver<Matrix>::findinitialInterval(int min, int max) const
{
//  std::cout<<min<<" "<<max<<std::endl;
    Interval<datatype> interval(diag[max - 1] - std::abs(lower[max - 1]), diag[max - 1] + std::abs(lower[max - 1]), 0, max - min);
    interval.print();
    for (int k = max-2; k >= min; --k)//determine an initial interval by the Gershgorin circles
    {
        datatype dis = std::abs(lower[k]) + std::abs(lower[k+1]);
        datatype testmax = diag[k] + dis;
        datatype testmin = diag[k] - dis;
//	std::cout<<lower[k]<<" "<<lower[k+1]<<" "<<diag[k]<<" "<<testmax<<" "<<testmin<<std::endl;
        if (testmin < interval.xmin) interval.xmin = testmin;
        if (testmax > interval.xmax) interval.xmax = testmax;
    }
    return interval;
}

template <typename T>
struct BisectInterval
{
    Interval<T> iv;
    uint refinement;
    BisectInterval(Interval<T>& a) : iv(a), refinement(0) {}
};

/**
 * As the name implies, this is supposed to only work for LDL^T factorizations
 */
template <typename Matrix>
uint EigenSystemSolver<Matrix>::sturmcountLDL(const TwistedMatrix<datatype>& mat, datatype x) const
{
    datatype qk = mat.diag[0] - x;
    uint a = (qk < 0? 1 : 0);
    for (uint k = 1; k < mat.len; ++k)
    {
        datatype ldp = mat.diag[k-1] * mat.off[k-1];//this lower[k]
//      std::cout<<k-1 <<" "<<ldp<<"  "<<lower[k]<<std::scientific<< std::abs(ldp-lower[k])<<std::endl;
        datatype sub =  ldp*(ldp/qk);
        qk = (mat.diag[k] + ldp*mat.off[k-1] - x) - sub;
        if (qk < 0.0 ) ++a;
    }
    return a;
}

template <class Matrix>
void EigenSystemSolver<Matrix>::refineSingletonInterval(Interval<datatype>& iv, const TwistedMatrix<datatype>& m) const
{
    datatype oldiv = iv.xmax - iv.xmin;
    const datatype fac = 2.0;
    do
    {
        datatype x0 = (iv.xmax + iv.xmin)*0.5;
        uint nrofzeroes = sturmcountLDL(m, x0);
        if (nrofzeroes == iv.EVsxmin)//lower bound refinement
            iv.xmin = x0;
        else//upper bound refinement
            iv.xmax = x0;
    } while ((oldiv <= fac*(iv.xmax - iv.xmin)) && ((iv.xmax - iv.xmin) > std::numeric_limits<datatype>::epsilon()));
    return;
}

template <class Matrix>
void EigenSystemSolver<Matrix>::bisectLDL(std::vector<Interval<datatype> >& finalintervals, ProblemSet<datatype>& ps) const
{
    uint n = ps.matrix.len;
    //indices of the EVs
    std::list<BisectInterval<datatype> > testintervals;
    testintervals.push_back(ps.iv);
    testintervals.front().iv.substructure.clear();
    datatype avdensity = (ps.iv.EVsxmax - ps.iv.EVsxmin)/(ps.iv.xmax - ps.iv.xmin);// a measure for the average density of the eigenvalues
    const datatype densecluster = 4 * (ps.level + 1);
    for (uint k = 0; (k < 14) && (testintervals.size() > 0); ++k)// let's try to improve the approximately first three digits( -> 14 bisections)
    {
        for (typename std::list<BisectInterval<datatype> >::iterator it = testintervals.begin(); it != testintervals.end(); ++it)
        {
            std::cout<<"Test-Interval"<<std::endl;
            it->iv.print();
            datatype x0 = (it->iv.xmax + it->iv.xmin)*0.5;
            uint nrofzeroes = sturmcountLDL(ps.matrix, x0);
            uint lowerintervalEVs = nrofzeroes - it->iv.EVsxmin;
            uint upperintervalEVs = it->iv.EVsxmax - nrofzeroes;
            if (lowerintervalEVs == 0)//shorten interval
            {
                it->iv.xmin = x0;
            }
            else if (upperintervalEVs == 0)//shorten interval
            {
                it->iv.xmax = x0;
            }
            else
            {//need to introduce new interval
                Interval<datatype> temp(it->iv.xmin, x0, it->iv.EVsxmin, nrofzeroes);
                testintervals.insert(it, BisectInterval<datatype>(temp));
                it->iv.xmin = x0;
                it->iv.EVsxmin = nrofzeroes;
            }
        }
        //now test for separation
        for (typename std::list<BisectInterval<datatype> >::iterator it = testintervals.begin(); it != testintervals.end();/*nothing*/)
        {
            uint evs = it->iv.EVsxmax - it->iv.EVsxmin;
            if (testintervals.size() >1)
            {
                if (it == testintervals.begin())//take special care for handling the beginning
                {
//                    std::cout<<"testing beginning"<<std::endl;
                    typename std::list<BisectInterval<datatype> >::iterator successor = it;
                    ++successor;
                    if (relgap(it->iv.xmax, successor->iv.xmin) > gaptol)
                    {
                        if (evs == 1)
                        {
//                            std::cout<<"removing singleton!"<<std::endl;
                            finalintervals.push_back(it->iv);
                            testintervals.pop_front();
                            it = testintervals.begin();
                        }
                        else
                        {
                            datatype density = evs/(it->iv.xmax - it->iv.xmin);
                            if ((density > densecluster* avdensity) || (relgap(it->iv.xmax, it->iv.xmin) <= gaptol))//high density of Eigenvalues, expect clustering
                            {
                                finalintervals.push_back(it->iv);
                                testintervals.pop_front();
//                                std::cout<<"removing cluster!"<<std::endl;
                                it = testintervals.begin();
                            }
                            else
                                ++it;
                        }
                    }
                    else
                    {
                        ++it;
                    }
                }
                else if (it == --(testintervals.end()))//should be the iterator to the last element
                {
//                    std::cout<<"testing end"<<std::endl;
                    typename std::list<BisectInterval<datatype> >::iterator predecessor = it;
                    --predecessor;
                    if (relgap(predecessor->iv.xmax, it->iv.xmin) > gaptol)
                    {
                        if (evs == 1)
                        {
//                            std::cout<<"removing singleton!"<<std::endl;
                            finalintervals.push_back(it->iv);
                            testintervals.pop_back();
                            it = testintervals.end();
                        }
                        else
                        {
                            datatype density = evs/(it->iv.xmax - it->iv.xmin);
                            if ((density > densecluster* avdensity) || (relgap(it->iv.xmax, it->iv.xmin) <= gaptol))
                            {
                                finalintervals.push_back(it->iv);
                                testintervals.pop_back();
                                it = testintervals.end();
                            }
                            else
                                ++it;
                        }
                    }
                    else
                    {
                        ++it;
                    }
                }
                else
                {//not at beginning or end
                    typename std::list<BisectInterval<datatype> >::iterator predecessor = it;
                    --predecessor;
                    typename std::list<BisectInterval<datatype> >::iterator successor = it;
                    ++successor;
//                    std::cout<<"testing generic"<<std::endl;
                    if ((relgap(it->iv.xmin, predecessor->iv.xmax) > gaptol) && (relgap(it->iv.xmax, successor->iv.xmin) > gaptol) )
                    {//at least it's well seperated...
                        if (evs == 1)//singleton!
                        {
                            finalintervals.push_back(it->iv);
                            it = testintervals.erase(it);
                        }
                        else
                        {
                            datatype density = evs/(it->iv.xmax - it->iv.xmin);
                            if ((density > densecluster *avdensity) || (relgap(it->iv.xmax, it->iv.xmin) <= gaptol))
                            {
                                finalintervals.push_back(it->iv);
                                it = testintervals.erase(it);
//                                std::cout<<"removing cluster!"<<std::endl;
                            }
                            else
                                ++it;
                        }
                    }
                    else
                        ++it;
                }
            }
            else
            {
                //
                if (evs == 1)
                {
//                    std::cout<<"removing singleton!"<<std::endl;
                    finalintervals.push_back(it->iv);
                    testintervals.pop_back();
                    it = testintervals.end();
                }
                else
                {
                    datatype density = evs/(it->iv.xmax - it->iv.xmin);
                    if ((density > densecluster* avdensity) || (relgap(it->iv.xmax, it->iv.xmin) <= gaptol))
                    {
                        finalintervals.push_back(it->iv);
                        testintervals.pop_back();
                    }
                    it = testintervals.end();
                }
            }
        }
    }
    //now follows clean-up work over possible remaining intervals
    while (testintervals.size() > 1)
    {
        typename std::list<BisectInterval<datatype> >::iterator it = testintervals.begin();
        ++it;
        if (relgap(it->iv.xmin, testintervals.front().iv.xmax) > gaptol)
        {
            finalintervals.push_back(testintervals.front().iv);
        }
        else
        {
            //fuse intervals
            Interval<datatype> newint = testintervals.front().iv;
            Interval<datatype> myint = it->iv;
            myint.substructure.clear();
            if (testintervals.front().iv.substructure.size() > 0)
            {
                it->iv.substructure = testintervals.front().iv.substructure;
                testintervals.front().iv.substructure.clear();
            }
            else
            {
                it->iv.substructure.push_back(testintervals.front().iv);
            }
            it->iv.substructure.push_back(myint);
            it->iv.xmin = testintervals.front().iv.xmin;
            it->iv.EVsxmin = testintervals.front().iv.EVsxmin;
        }
        testintervals.pop_front();
    }
    if (testintervals.size() == 1)
        finalintervals.push_back(testintervals.front().iv);
    std::cout<<"============Final Intervals after bisection==========="<<std::endl;
    for (uint k = 0; k < finalintervals.size(); ++k)
    {
        finalintervals[k].print();
        for (uint j = 0; j < finalintervals[k].substructure.size(); ++j)
        {
            std::cout<<"   ";
            finalintervals[k].substructure[j].print();
        }
    }
    std::cout<<"======================================================"<<std::endl;
}

template <class Matrix>
void EigenSystemSolver<Matrix>::bisect(std::vector<Interval<datatype> >& finalintervals, Interval<datatype>& interval, datatype shift)
{//unused
    uint n = orig.Rows();
    //indices of the EVs
    std::list<BisectInterval<datatype> > testintervals;
    testintervals.push_back(interval);
    while (testintervals.size() > 0)
    {
        //split interval
        BisectInterval<datatype>& test = testintervals.back();
        datatype x0 = test.iv.xmin + (test.iv.xmax - test.iv.xmin)/2.0;
        uint nrofzeroes = sturmcount(x0, shift);
        uint lowerintervalEVs = nrofzeroes - test.iv.EVsxmin;
        uint upperintervalEVs = test.iv.EVsxmax - nrofzeroes;
        if (lowerintervalEVs == 0)//shortening of interval
        {
            test.iv.xmin = x0;
            if (testintervals.size() > 1)
            {
                if (relgap((++(testintervals.rbegin()))->iv.xmax, x0 ) > gaptol)//found sth.!
                {
                    finalintervals.push_back(test.iv);
                    testintervals.pop_back();
                }
            }
            else
            {
                finalintervals.push_back(test.iv);
                testintervals.pop_back();
            }
        }
        else if (upperintervalEVs == 0)//shortening of interval
        {
            test.iv.xmax = x0;
            if (testintervals.size() == 1)
            {
                finalintervals.push_back(test.iv);
                testintervals.pop_back();
            }
        }
        else
        {
            datatype maximum = test.iv.xmax;
            uint exmax = test.iv.EVsxmax;
            test.iv.xmax = x0;
            test.iv.EVsxmax = nrofzeroes;
            Interval<datatype> temp(x0, maximum, nrofzeroes, exmax);
            testintervals.push_back(BisectInterval<datatype>(temp));
//	  std::cout<<"Adding interval"<<std::endl;
//	  testintervals.back().iv.print();
        }
    }
    for (uint k = 0; k < finalintervals.size(); ++k)
        finalintervals[k].print();
}

template <class Matrix>
bool EigenSystemSolver<Matrix>::testeigensystem(std::vector<SolutionPair<datatype> >& esystem, Matrix& z)
{
  uint n = esystem.size();
  Matrix rfac(n, n);
  Matrix lfac(rfac);
  for(uint k = 0; k < n; ++k)
    for(uint j = 0; j < n; ++j)
      rfac(j, k) = esystem[k].evector[j];
  Matrix rfac2(z * rfac);
  for(uint k = 0; k < n; ++k)
    for(uint j = 0; j < n; ++j)
      lfac(k, j) = conj(rfac2(j, k));
    Matrix res = lfac * orig * rfac2;
  std::ofstream resultfile("tests.txt", std::ios::app);
  resultfile<<std::setprecision(std::numeric_limits<datatype>::digits10);
  for(uint i = 0; i < n; ++i)
  {
    resultfile<<std::scientific<<std::showpos<<"EV "<<esystem[i].evalue<<" compares to "<<res(i, i)<<" like "<<std::fabs(esystem[i].evalue - res(i, i));
    datatype offdiag = 0;
    for(uint k = 0; k < i; ++k)
    {
      offdiag += std::fabs(res(i, k));
    }
    for(uint k = i+1; k < n; ++k)
    {
      offdiag += std::fabs(res(i, k));
    }
    resultfile<<"   off diagonal sum: "<<offdiag<<std::endl;
  }
  bool eigenvalues_are_ok = true;
  for(uint k = 0; k < n; ++k)
    eigenvalues_are_ok = eigenvalues_are_ok && (std::abs(res(k,k) - esystem[k].evalue) < 1E-6);
  bool offs_are_small = true;
  for(uint k = 0; k < n; ++k)
  {
    for(uint j = 0; j < k; ++j)
      offs_are_small = offs_are_small && (std::abs(res(k,j)) < 1E-6);
        for(uint j = k+1; j < n; ++j)
      offs_are_small = offs_are_small && (std::abs(res(k,j)) < 1E-6);
  }
  return offs_are_small && eigenvalues_are_ok;
}

template <typename Matrix>
SolutionPair<typename EigenSystemSolver<Matrix>::datatype> EigenSystemSolver<Matrix>::rqi(const TwistedMatrix<datatype>& m, Interval<datatype>& iv, uint min) const
{
    uint n = m.len;
    datatype* diag_out;
#ifndef __USE_XOPEN2K
    diag_out = new datatype[3*n];
#else
    
    int errret = posix_memalign((void**)&diag_out, 16, 3*n*sizeof(datatype));
    if(unlikely(errret != 0))
      throw std::bad_alloc();
#endif
    datatype* off_out = diag_out + n;
    datatype* diag1 = off_out +n;
    datatype mingammak;//Is initialized in just four lines
    uint mingammak_idx = 0;
    datatype evnorm;
    SolutionPair<datatype> retval(0.5*(iv.xmax + iv.xmin), orig.Rows());//create eigenvector with the full size of the matrix. We only write to the non zero entries
    std::cout<<"initial guess: "<<retval.evalue<<std::endl;
    uint iter = 0;
    bool stoprqi = false;
    bool notconverged = false;
    datatype rqicorrection = 0.0;
    int envmin = 0;
    int envmax = n-1;
    do
    {
        dtwqds(m.diag, m.off, diag_out, off_out, retval.evalue, m.twist, n-1, n);//create an LDL version of the matrix
        dtwqds(m.diag, m.off, diag1, off_out, retval.evalue, m.twist, 0, n);//create URU form of the matrix
        //from those two decompositions all twist indices can be recovered
        mingammak = (std::isfinite(diag1[envmin]) ? diag1[envmin] : std::numeric_limits<datatype>::max());//replace by a large number if the first is NaN
        mingammak_idx = 0;//necessary here //FIXME: NaN?
        for (uint k = envmin + 1; k < envmax; ++k)
        {
//	  std::cout<<diag_out[k]<<" "<<diag1[k]<<" "<<off_out[k]<<" "<<diag1[k+1]<<" "<<off_out[k]*off_out[k]*diag1[k+1]<<" "<< diag_out[k] - off_out[k]*off_out[k]*diag1[k+1]<<std::endl;
            datatype temp = diag_out[k] -  off_out[k]*off_out[k]*diag1[k+1];
	    //the alternative is:  diag1[k]-off_out2[k-1]*off_out2[k-1]*diag_out[k-1]    , where off_out2 denotes the off-diagonal entries of the LDL factorization
            /*datatype x1,y1,x2,y2,x3,y3;
            datatype temp = diag_out[k];
            twoproduct(-off_out[k], diag1[k+1], x1, y1);
            twoproduct(x1, off_out[k], x2, y2);
            twoproduct(y1, off_out[k], x3, y3);
            twosum(temp, x2, temp, x1);
            twosum(temp, x3, temp, x2);
            twosum(y2, y3, y1, y2);
            twosum(temp, y1, temp, x3);
            temp = temp + ((x3 + y2) + (x1 + x2));*/
//std::cout<<temp<<"    "<<diag_out[k] -  off_out[k]*off_out[k]*diag1[k+1]<<std::endl;
            if (std::abs(temp) < std::abs(mingammak))
            {
                mingammak = temp;
                mingammak_idx = k;
            }
        }
        if ((std::abs(diag_out[envmax]) < std::abs(mingammak)))//last element
        {
            mingammak = diag_out[envmax];
            mingammak_idx = envmax;
        }
        std::cout<<mingammak<<std::endl;
        dtwqds(m.diag, m.off, diag_out, off_out, retval.evalue, m.twist, mingammak_idx, n);
        mingammak = diag_out[mingammak_idx];//let's retrieve it here for consistency
        std::cout<<mingammak<<"       "<<mingammak_idx<<std::endl;
        retval.evector[min + mingammak_idx] = 1;
//        evnorm = 1.0;
        datatype p = 1.0;
        datatype s = 0.0;
//we use an algorithm to sum the norm in twice the initial precision
        for (int i = (mingammak_idx - 1); i >= envmin; --i)
        {
            datatype h, r, q;
            retval.evector[min + i] = - off_out[i] * retval.evector[min + i+1];
            squaretransform(retval.evector[min + i], h, r);
            twosum(p, h, p, q);
            s = s + (q + r);
//            evnorm += retval.evector[min + i]*retval.evector[min + i];
        }
        for (uint i = (mingammak_idx + 1); i <= envmax; ++i)
        {
            datatype h, r, q;
//    	    std::cout<<"off: "<<off_out[i-1]<<std::endl;
            retval.evector[min + i] = - off_out[i-1]*retval.evector[min + i - 1];
            squaretransform(retval.evector[min + i], h, r);
            twosum(p, h, p, q);
            s = s + (q + r);
//            evnorm += retval.evector[min + i]*retval.evector[min + i];
        }
        evnorm = p + s;
        uint envminguess = envmin;
	uint envmaxguess = envmax;
//	const datatype myeps = evnorm * std::numeric_limits<datatype>::epsilon()*std::numeric_limits<datatype>::epsilon();//the vector is not normalized yet
//	while(std::abs(retval.evector[envminguess]) < myeps) ++envminguess;
//	while(std::abs(retval.evector[envmaxguess]) < myeps) --envmaxguess;
//        std::cout<<evnorm<<std::endl;
        rqicorrection = mingammak/evnorm;
	std::cout<<"RQI - Correction: "<<rqicorrection<<std::endl;
        datatype correction;
	datatype newev = retval.evalue + rqicorrection; 
        if ((newev < iv.xmax) && ((newev > iv.xmin )))
        {
            std::cout<<"Using RQI correction"<<std::endl;
            if (rqicorrection > 0)
            {
                iv.xmin = retval.evalue;
		iv.xmax = retval.evalue + 2.5*rqicorrection;
//	    iv.xmax = std::min(retval.evalue + 2.0*rqicorrection, iv.xmax);
            }
            else
            {
	      iv.xmin = retval.evalue + 2.5*rqicorrection;
	      iv.xmax = retval.evalue;
//	    iv.xmin = std::max(retval.evalue + 2.0*rqicorrection, iv.xmin);
            }
            envmin = envminguess;
            envmax = envmaxguess;
            retval.evalue += rqicorrection;
            evnorm = std::sqrt(evnorm);
            correction = std::abs(mingammak)/evnorm;
        }
        else
        {
            std::cout<<"No improvement! Refining..."<<std::endl;
            refineSingletonInterval(iv, m);
            retval.evalue = 0.5*(iv.xmax + iv.xmin);
            correction = std::abs(iv.xmax - iv.xmin);
            std::cout<<"Better Guess: "<<retval.evalue<<std::endl;
	    evnorm = std::sqrt(evnorm);
        }
        std::cout<<"envelope: "<<envmin<<" -> "<<envmax<<std::endl;
        std::cout<<"Current Guess for EV: "<<retval.evalue<<std::endl;
        notconverged = (iter >= static_cast<uint>(std::numeric_limits<datatype>::digits));
        ++iter;
        stoprqi = (correction < n*std::abs(retval.evalue)* std::numeric_limits<datatype>::epsilon());
    } while (!stoprqi && !notconverged );
    //Estimate as given in Willems thesis(section 2.2). Since we encountered one case where a small eigenvalue led to a very small accuracy limit we give some headroom
    //the denominator on the lhs comes from the RQI in Algorithm 2.5
    if (notconverged)
        std::cout<<"WARNING! RQI FAILED TO CONVERGE!"<<std::endl;
    for (uint k = min; k < min + n; ++k)
        retval.evector[k] /= evnorm;
    std::cout<<"Number of Iterations: "<<iter<<std::endl;
#ifndef __USE_XOPEN2K
    delete [] diag_out;
#else
    free(diag_out);//it's defintely a valid pointer down here... else we would have an error above
#endif
    return retval;
}

template <class Matrix>
TwistedMatrix<typename EigenSystemSolver<Matrix>::datatype> EigenSystemSolver<Matrix>::initialdtwqds( const datatype *const d, const datatype *const l, const datatype shift, uint twist, uint diaglen) const
{//the initial dtwqds as outlined in Dhillons thesis. Generates from a bidiagonal factorization a twisted factorization at twist t
    datatype* s = new datatype[2*diaglen];
    datatype* p = s + diaglen;
    TwistedMatrix<datatype> retval;
    retval.twist = twist;
    retval.diag = new datatype[diaglen];
    retval.off = new datatype[diaglen];
    retval.len = diaglen;
    s[0] = -shift;
    for (uint i = 0; i < twist; ++i)
    {
        retval.diag[i] = s[i] + d[i];
        retval.off[i] = (d[i]*l[i])/retval.diag[i];
        s[i+1] = retval.off[i]*l[i]*s[i] - shift;
    }
    p[diaglen-1] = d[diaglen-1] - shift;
    for (int i = static_cast<int>(diaglen) - 2; i >= static_cast<int>(twist); i--)
    {
        retval.diag[i+1] = d[i]*l[i]*l[i] + p[i+1];
        datatype t = d[i]/retval.diag[i+1];
        retval.off[i] = l[i]*t;
        p[i] = p[i+1]*t - shift;
    }

    if (twist != 0)
        retval.diag[twist] = p[twist] + retval.off[twist - 1]* l[twist - 1]  * s[twist -1];
    else
        retval.diag[twist] = s[twist] + d[twist]/retval.diag[twist + 1]*p[twist + 1];
    delete [] s;
    return retval;
}

template <class Matrix>
ProblemSet<typename EigenSystemSolver<Matrix>::datatype> EigenSystemSolver<Matrix>::createInitialProblem(uint min, uint max)
{
    ProblemSet<datatype> problem;
    problem.level = 0;
    int diaglen = max - min;
    datatype* d = new datatype[diaglen];
    datatype* l = new datatype[diaglen];
    Interval<datatype> initialinterval(findinitialInterval(min, max));
//FIXME: since we now use the possible block structure in the tridiagonal matrix, the specification of desired Eigenvalues is broken!!!
    initialinterval.print();
    if ((min == 0) && ( max == orig.Rows() - 1))
    {
        initialinterval.EVsxmax = max_idx + 1;
        initialinterval.EVsxmin = min_idx;
    }
    else
    {//get consistent counts
        initialinterval.EVsxmax = diaglen;
        initialinterval.EVsxmin = 0;
    }
    initialinterval.print();
    pivmin = 10.0/(std::numeric_limits<datatype>::max()-100000)* std::max( initialinterval.xmax * initialinterval.xmax, 1.0);
    //first let's refine the upper index
//    std::cout<<max_idx<<" "<<min_idx<<std::endl;
    if ((min == 0) && ( max == orig.Rows() - 1))//If the matrix has no degeneracy we allow for specification of the desired indices
    {
        if (max_idx != (orig.Rows() -1))
        {
            datatype xmin = initialinterval.xmin;
            datatype xmax = initialinterval.xmax;
            datatype x0 = xmin + (initialinterval.xmax - initialinterval.xmin)* static_cast<datatype>(initialinterval.EVsxmax)/(orig.Rows());
            uint nrofzeroes = sturmcount(x0, 0.0, min, max);
            while (nrofzeroes != initialinterval.EVsxmax)
            {
                if (nrofzeroes < initialinterval.EVsxmax)
                {
                    xmin = x0;
                    x0 = (x0 + xmax)*0.5;
                }
                if (nrofzeroes > initialinterval.EVsxmax)
                {
                    xmax = x0;
                    x0 = (x0 + xmin)*0.5;
                }
                nrofzeroes = sturmcount(x0, 0.0, min, max);
            }
            initialinterval.xmax = x0;
        }
        if (min_idx != 0)
        {
            datatype xmin = initialinterval.xmin;
            datatype xmax = initialinterval.xmax;
            initialinterval.print();
            datatype x0 = xmin; //+ (initialinterval.xmax - initialinterval.xmin)* static_cast<float>((initialinterval.EVsxmax - initialinterval.EVsxmin)/(orig.Rows()));
            uint nrofzeroes = sturmcount(x0, 0.0, min, max);
            uint iter = 0;
            while (nrofzeroes != initialinterval.EVsxmin)
            {
                std::cout<<x0<<" "<<nrofzeroes<<std::endl;
                if (nrofzeroes < initialinterval.EVsxmin)
                {
                    xmin = x0;
                    x0 = (x0 + xmax)*0.5;
                }
                if (nrofzeroes > initialinterval.EVsxmin)
                {
                    xmax = x0;
                    x0 = (x0 + xmin)*0.5;
                }
                nrofzeroes = sturmcount(x0, 0.0, min, max);
                iter++;
            }
            initialinterval.xmin = x0;
        }
    }
//    std::cout<<initialinterval.xmin<<" "<<sturmcount(initialinterval.xmin, 0.0,min, max)<<" ? "<<initialinterval.EVsxmin<<std::endl;
//    std::cout<<initialinterval.xmax<<" "<<sturmcount(initialinterval.xmax, 0.0, min, max)<<" ? "<<initialinterval.EVsxmax<<std::endl;
    //The next few lines determine a shift of the initial interval to sth. a bit larger than zero.
    //That way the matrix is guaranteed to be positive definite.
    datatype initialshift = initialinterval.xmin;
    initialshift = (1.0 - copysign(0.1, initialshift) ) *initialshift;
    std::cout<<"SHIFT: "<<initialshift<<std::endl;
    initialinterval.xmin -= initialshift;
    initialinterval.xmax -= initialshift;
    initialinterval.print();
    std::cout<<std::endl;
    generatelowerBidiagonalfactorization(d, l, min, max, initialshift);
    problem.matrix = initialdtwqds(d, l,/* initialshift*/ 0.0, diaglen - 1, diaglen);
    problem.shift = initialshift;
    if ((min == 0) && ( max == orig.Rows() - 1))
    {
        problem.min_EV_idx = min_idx;
        problem.max_EV_idx = max_idx;
    }
    else
    {
        problem.min_EV_idx = 0;
        problem.max_EV_idx = diaglen - 1;
    }
    problem.iv = initialinterval;
    delete [] d;
    delete [] l;
    return problem;
}

template <typename T>
T elg(T* newdiag, T* olddiag, T shift, uint len)
{
    T elg = 0.0;
    for (uint i = 0; i < len; ++i)
    {
        elg = std::max(elg, std::sqrt( std::abs( newdiag[i]/(olddiag[i] - shift)  ) ));
    }
    return elg;
}

template <typename T>
struct Shiftelg
{
    T shift;
    T elg;
    uint idx;
};

template <typename T>
inline bool sortelg(const Shiftelg<T>& a, const Shiftelg<T>& b)
{
    return (a.elg - b.elg) < 0;
}

template <class Matrix>
void EigenSystemSolver<Matrix>::mrrr(uint min, uint diaglen, std::vector<ProblemSet<datatype> >& nodes, std::vector<SolutionPair<datatype> >& esystem)
{
    //An implementation of the MRRR algorithm.
    while (!nodes.empty())
    {
        std::vector<Interval<datatype> > intervals;
        ProblemSet<datatype> mynode = nodes.back();
        bisectLDL(intervals, mynode);
        nodes.pop_back();
        for (typename std::vector<Interval<datatype> >::iterator it = intervals.begin(); it != intervals.end(); ++it)
        {
            if ((it->EVsxmax - it->EVsxmin) == 1)//found a singleton!
            {
                //FIXME: Heuristic(Flo): we need to refine the Eigenvalues to the point that their Interval is smaller than their estimated value
                //Therefore a condition like (max - min) <= alpha * 0.5 * (max + min)
                //should hold.
                static const float alpha = 1.0/4.0;
                while (std::abs(it->xmax - it->xmin) > alpha*0.5*std::abs(it->xmax + it->xmin))
                {
                    std::cout<<"Refining singleton interval!"<<std::endl;
                    refineSingletonInterval(*it, mynode.matrix);
                }
                it->print();
                esystem.push_back(rqi(mynode.matrix, *it, min));
                esystem.back().evalue += mynode.shift;
            }
            else// it is a cluster
            {
                datatype alpha = 100*diaglen*std::numeric_limits<datatype>::epsilon();
                std::cout<<"Dealing with cluster"<<std::endl;
                ProblemSet<datatype> newnode;
                datatype new_shift = 0.95*it->xmin;
                newnode.iv = *it;
                newnode.level = mynode.level + 1;
                if (new_shift < 0) alpha = -alpha;
                newnode.iv.xmin -= (1.0 + alpha)*new_shift;
                newnode.iv.xmax -= (1.0 - alpha)*new_shift;
                newnode.iv.print();
                newnode.shift = mynode.shift + new_shift;
                newnode.min_EV_idx = it->EVsxmin;
                newnode.max_EV_idx = it->EVsxmax;
                newnode.matrix.diag = new datatype[diaglen];
                newnode.matrix.off = new datatype[diaglen];
                newnode.matrix.twist = diaglen - 1;
                newnode.matrix.len = diaglen;
                std::cout<<"New shift: "<<newnode.shift<<std::endl;
                dtwqds(mynode.matrix.diag, mynode.matrix.off, newnode.matrix.diag, newnode.matrix.off, new_shift, mynode.matrix.twist, newnode.matrix.twist, diaglen);
                datatype localelg = elg(newnode.matrix.diag, diag + min, newnode.shift, diaglen);
                std::cout<<"ELG: "<<localelg<<std::endl;
//		newnode.matrix.print();
                bool is_consistent = (sturmcountLDL(newnode.matrix, newnode.iv.xmin) == newnode.iv.EVsxmin) && (sturmcountLDL(newnode.matrix, newnode.iv.xmax) == newnode.iv.EVsxmax);
                std::cout<<sturmcountLDL(newnode.matrix, newnode.iv.xmin)<<" == "<<newnode.iv.EVsxmin<<" && "<<sturmcountLDL(newnode.matrix, newnode.iv.xmax)<<" == "<<newnode.iv.EVsxmax<<std::endl;
                if (likely((localelg < 200) && is_consistent))
                {
                    nodes.push_back(newnode);
                }
                else
                {
                    std::cout<<"WARNING! left representation not robust"<<std::endl;
                    if (it->substructure.size() == 0)
                    {
                        std::cout<<"No substructure!"<<std::endl;
                        ProblemSet<datatype> newnodes[9];
                        datatype elg2[9] = {0};
                        uint bestidx = 9;
                        datatype currelg = is_consistent? localelg : std::numeric_limits<datatype>::max();
                        for (int k = 0; k < 9; ++k)
                        {
                            newnodes[k].iv = *it;
                            newnodes[k].level = newnode.level;
                            datatype new_shift = it->xmin + (k-1)*(it->xmax -it->xmin)*0.2;
                            newnodes[k].iv.xmin -= new_shift;
                            newnodes[k].iv.xmax -= new_shift;
                            newnodes[k].shift = mynode.shift + new_shift;
                            newnodes[k].min_EV_idx = it->EVsxmin;
                            newnodes[k].max_EV_idx = it->EVsxmax;
                            newnodes[k].matrix.len = diaglen;
                            newnodes[k].matrix.diag = new datatype[diaglen];
                            newnodes[k].matrix.off = new datatype[diaglen];
                            newnodes[k].matrix.twist = diaglen - 1;
                            std::cout<<"New shift: "<<newnodes[k].shift<<std::endl;
                            dtwqds(mynode.matrix.diag, mynode.matrix.off, newnodes[k].matrix.diag, newnodes[k].matrix.off, new_shift, mynode.matrix.twist, newnodes[k].matrix.twist, diaglen);
                            elg2[k] = elg(newnodes[k].matrix.diag, diag + min, newnodes[k].shift, diaglen);
                            std::cout<<"ELG2["<<k<<"]  =  "<<elg2[k]<<std::endl;
                            if (elg2[k] < currelg)
                            {
                                currelg = elg2[k];
                                bestidx = k;
                            }
                        }
                        if (bestidx == 9)
                        {
                            if (likely(is_consistent))//well at least it's consistent... very likely we end up here in the next iteration again...
                            {
                                nodes.push_back(newnode);
                                for (uint k = 0; k < 9; ++k)
                                {
                                    delete [] newnodes[k].matrix.diag;
                                    delete [] newnodes[k].matrix.off;
                                }
                            }
                            else
                            {
                                std::cout<<"CRAP!"<<std::endl;
                                nodes.push_back(newnode);
                                for (uint k = 0; k < 9; ++k)
                                {
                                    delete [] newnodes[k].matrix.diag;
                                    delete [] newnodes[k].matrix.off;
                                }
                            }
                        }
                        else
                        {
                            std::cout<<"choosing: "<<bestidx<<std::endl;
                            nodes.push_back(newnodes[bestidx]);
                            //let's check for consistency
                            std::cout<<sturmcountLDL(mynode.matrix, mynode.iv.xmin)<<" == "<<mynode.iv.EVsxmin<<" && "<<sturmcountLDL(mynode.matrix, mynode.iv.xmax)<<" == "<<mynode.iv.EVsxmax<<std::endl;
                            std::cout<<sturmcountLDL(newnode.matrix, newnode.iv.xmin)<<" == "<<newnode.iv.EVsxmin<<" && "<<sturmcountLDL(newnode.matrix, newnode.iv.xmax)<<" == "<<newnode.iv.EVsxmax<<std::endl;
                            std::cout<<sturmcountLDL(newnodes[bestidx].matrix, newnodes[bestidx].iv.xmin)<<" == "<<newnodes[bestidx].iv.EVsxmin<<" && "<<sturmcountLDL(newnodes[bestidx].matrix, newnodes[bestidx].iv.xmax)<<" == "<<newnodes[bestidx].iv.EVsxmax<<std::endl;
                            if (!((sturmcountLDL(newnodes[bestidx].matrix, newnodes[bestidx].iv.xmin) == newnodes[bestidx].iv.EVsxmin) && (sturmcountLDL(newnodes[bestidx].matrix, newnodes[bestidx].iv.xmax) == newnodes[bestidx].iv.EVsxmax)))
                            {
                                std::cout<<"sturm counts not consistent!!!!!!!"<<std::endl;
                                exit(-1);
                            }
                            for (uint k = 0; k < bestidx; ++k)
                            {
                                delete [] newnodes[k].matrix.diag;
                                delete [] newnodes[k].matrix.off;
                            }
                            for (uint k = bestidx + 1; k < 9; ++k)
                            {
                                delete [] newnodes[k].matrix.diag;
                                delete [] newnodes[k].matrix.off;
                            }
                            delete [] newnode.matrix.diag;
                            delete [] newnode.matrix.off;
                        }
                    }
                    else
                    {//we have substructure!
                        std::cout<<"Using substructure!"<<std::endl;
                        uint nrshifts = static_cast<uint>(3 * it->substructure.size() - 1);
                        Shiftelg<datatype>* shifts = new Shiftelg<datatype>[nrshifts];
                        ProblemSet<datatype>* data = new ProblemSet<datatype>[nrshifts];
                        Shiftelg<datatype>* beg = shifts;
                        //let's first go straight for the gaps
                        for (uint i = 0; i < (it->substructure.size() - 1); ++i)
                            (shifts++)->shift = 0.5*(it->substructure[i].xmax + it->substructure[i+1].xmin);
                        //Now let's zoom to the intervals that we have.
                        for (uint i = 0; i < it->substructure.size(); ++i)
                        {
                            (shifts++)->shift = 0.94*it->substructure[i].xmin;
                            (shifts++)->shift = 1.06*it->substructure[i].xmax;
                        }
                        shifts = beg;
                        for (uint k = 0; k < nrshifts; ++k)
                        {
                            shifts[k].idx = k;
                            data[k].iv = *it;
                            data[k].level = newnode.level;
                            datatype new_shift = shifts[k].shift;
                            data[k].iv.xmin -= new_shift;
                            data[k].iv.xmax -= new_shift;
                            data[k].shift = mynode.shift + new_shift;
                            data[k].min_EV_idx = it->EVsxmin;
                            data[k].max_EV_idx = it->EVsxmax;
                            data[k].matrix.len = diaglen;
                            data[k].matrix.diag = new datatype[diaglen];
                            data[k].matrix.off = new datatype[diaglen];
                            data[k].matrix.twist = diaglen - 1;
                            std::cout<<"New shift: "<<data[k].shift<<std::endl;
                            dtwqds(mynode.matrix.diag, mynode.matrix.off, data[k].matrix.diag, data[k].matrix.off, new_shift, mynode.matrix.twist, data[k].matrix.twist, diaglen);
                            shifts[k].elg = elg(data[k].matrix.diag, diag, data[k].shift, diaglen);
                            std::cout<<"ELG2["<<k<<"]  =  "<<shifts[k].elg<<std::endl;
                        }
                        std::sort(shifts, shifts + nrshifts, sortelg<datatype>);
                        uint bestidx = 0;
                        if ((shifts[bestidx].elg > localelg) && is_consistent)//none of our new shifts was better than the old one-> let's take that
                        {
                            std::cout<<"Didn't find a new good shift. take our first guess..."<<std::endl;
                            nodes.push_back(newnode);
                            for (uint k = 0; k < nrshifts; ++k)
                            {
                                delete [] data[shifts[k].idx].matrix.diag;
                                delete [] data[shifts[k].idx].matrix.off;
                            }
                        }
                        else
                        {
                            uint bestdataidx = shifts[bestidx].idx;
                            datatype alpha = 100*diaglen*std::numeric_limits<datatype>::epsilon();
                            bool consistent = false;
                            while (!consistent && (bestidx < nrshifts))
                            {
                                data[bestdataidx].iv = *it;//reset the interval to the old value...
                                if (shifts[bestidx].shift < 0) alpha = -alpha;
                                //inflate boundaries
                                data[bestdataidx].iv.xmin -= (1.0 + alpha)* shifts[bestidx].shift;
                                data[bestdataidx].iv.xmax -= (1.0 - alpha)* shifts[bestidx].shift;
                                consistent = (sturmcountLDL(data[bestdataidx].matrix, data[bestdataidx].iv.xmin) == data[bestdataidx].iv.EVsxmin) && (sturmcountLDL(data[bestdataidx].matrix, data[bestdataidx].iv.xmax) == data[bestdataidx].iv.EVsxmax);
                                if (!consistent)
                                {
                                    std::cout<< bestidx <<" -> "<<bestdataidx<<" was not consistent!"<<std::endl;
                                    ++bestidx;
                                    bestdataidx = shifts[bestidx].idx;
                                }
                            }
                            std::cout<<shifts[bestidx].elg<<std::endl;
                            if (bestidx < nrshifts)
                            {
                                nodes.push_back(data[bestdataidx]);
                                //now the cleanup part
                                for (uint k = 0; k < bestidx; ++k)
                                {
                                    delete [] data[shifts[k].idx].matrix.diag;
                                    delete [] data[shifts[k].idx].matrix.off;
                                }
                                for (uint k = bestidx + 1; k < nrshifts; ++k)
                                {
                                    delete [] data[shifts[k].idx].matrix.diag;
                                    delete [] data[shifts[k].idx].matrix.off;
                                }
                                delete [] newnode.matrix.diag;
                                delete [] newnode.matrix.off;
                            }
                            else
                            {
                                std::cout<<"No consistent shift found!"<<std::endl;
                                if (likely(is_consistent))//well the old and first guess was at least consistent...
                                {
                                    nodes.push_back(newnode);
                                    for (uint k = 0; k < nrshifts; ++k)
                                    {
                                        delete [] data[shifts[k].idx].matrix.diag;
                                        delete [] data[shifts[k].idx].matrix.off;
                                    }
                                }
                                else
                                {
                                    std::cout<<"CRAP!!!!!!!!!!!"<<std::endl;
                                    exit(-1);
                                }
                            }
                        }
                        delete [] data;
                        delete [] shifts;
                    }
                    /*                Matrix N(diaglen, diaglen);
                                    Matrix D(diaglen, diaglen);
                                    for (uint k = 0; k < diaglen; ++k)
                                    {
                                        N(k,k) = 1.0;
                                        D(k,k) = newnode.matrix.diag[k];
                                    }
                                    for (uint k = 0; k < newnode.matrix.twist; ++k)
                                    {
                                        if (k+1< diaglen)
                                            N(k+1,k) = newnode.matrix.off[k];
                                    }
                                    for (uint k = newnode.matrix.twist + 1; k < diaglen; ++k)
                                    {
                                        if (k+1< diaglen)
                                            N(k,k+1) = newnode.matrix.off[k];
                                    }
                                    std::cout<<N<<std::endl;
                                    std::cout<<D<<std::endl;*/
                    /*    std::cout<<N*D*(~N)<<std::endl;
                    	exit(-1);*/
                }
            }
        }
        delete [] mynode.matrix.diag;
        delete [] mynode.matrix.off;
    }
    return;
}

template <typename T>
void findmaxmin(const T& a, const T&b, T& max, T& min)
{
    if (std::abs(a) > std::abs(b))
    {
        max = a;
        min = b;
    }
    else
    {
        max = b;
        min = a;
    }
    return;
}

template <typename T>
T mymax(const T& a, const T&b)
{
    return std::abs(a) > std::abs(b) ? a : b;
}

template <typename T>
T mymin(const T& a, const T&b)
{
    return std::abs(a) <= std::abs(b) ? a : b;
}

template <class Matrix>
std::vector<SolutionPair<typename EigenSystemSolver<Matrix>::datatype> > EigenSystemSolver<Matrix>::eigensystem()
{/*
  datatype* di = new datatype[orig.Rows()];
  datatype* of = new datatype[orig.Rows()];
  datatype* diout = new datatype[orig.Rows()];
  datatype* ofout = new datatype[orig.Rows()];
  datatype* diout2 = new datatype[orig.Rows()];
  datatype* ofout2 = new datatype[orig.Rows()];
  datatype* diout3 = new datatype[orig.Rows()];
  datatype* ofout3 = new datatype[orig.Rows()];  
  
  generatelowerBidiagonalfactorization(diout, ofout, 0, orig.Rows(), 0.0);
  dtwqds(diout, ofout, diout2, ofout2, 0.5, orig.Rows()-1, 0, orig.Rows());
  
  
  for(uint k = 0; k < orig.Rows(); ++k)
    std::cout<<diag[k]- 1<<" "<<lower[k]<<std::endl;
  dtwqds(diout2, ofout2, diout3, ofout3, 0.5, 0, 0, orig.Rows());
  for(uint k = 0; k < orig.Rows(); ++k)
  {
    std::cout<<diout2[k]<<" "<<ofout2[k]<<" "<<diout2[k+1]*ofout2[k]<<" "<<lower[k+1]<<" "<<diout[k] * ofout[k]<<std::endl;
  }
  exit(-1);*/
    lower[0] = 0.0;
    datatype frobeniusnorm = diag[0]*diag[0];
//  the frobeniusnorm fulfills for a matrix A of rank n
    // ||A||_2 <= ||A||_F <= ||A||_2 * sqrt(n)
    //hence we use it to estimate a lower bound for the 2 norm
    for (uint k = 1; k < orig.Rows(); ++k)
    {
        frobeniusnorm += 2.0*lower[k]*lower[k] + diag[k] * diag[k];
    }
    datatype firstestimate = std::sqrt(frobeniusnorm/orig.Rows()*std::numeric_limits<datatype>::epsilon());
    std::cout<<firstestimate<<std::endl;
    std::vector<SolutionPair<datatype> > esystem;
    uint min = 0;
    for (uint k = 1; k < orig.Rows(); ++k)
    {
//      std::cout<<k<<"  "<<std::scientific<<lower[k]<<std::endl;
//      std::cout<<std::abs(lower[k])* (std::abs(lower[k-1])+std::abs(lower[k+1])) * 2.0 <<" < "<<std::numeric_limits<datatype>::epsilon()*(std::abs(diag[k]) + std::abs(diag[k+1])*std::abs(diag[k+1] - diag[k]))<<std::endl;
        if ((std::abs(lower[k]) < firstestimate)//a first crude estimate that is usually sufficient
	  /*&& (std::abs(lower[k])* (std::abs(lower[k-1])+std::abs(lower[k+1])) * 2.0 < 
	  std::numeric_limits<datatype>::epsilon()*(std::abs(diag[k]) + std::abs(diag[k+1])*std::abs(diag[k+1] - diag[k]))
	  )//this estimate is from B.N. Parlett: The SEP ex. 7.11.4 He thinks it is good. Why argue with the master?
	*/)
        {
            lower[k] = 0.0;
            std::cout<<"splitting found at: "<<k<<std::endl;
            switch (k - min)
            {
            case 1:
                esystem.push_back(SolutionPair<datatype>(diag[k-1], orig.Rows()));
                esystem.back().evector[k-1] = 1.0;
                break;
            case 2:
            {
                datatype B = diag[k-2] + diag[k-1];
                datatype mydet = -lower[k-1] * lower[k-1] + diag[k-2] * diag[k-1];
                datatype discriminant;
                if (std::signbit(diag[k-2]) != std::signbit(diag[k-1]))
                {
                    discriminant = std::hypot(diag[k-1] - diag[k-2], 2.0 * lower[k-1]);
                }
                else
                {
                    discriminant = std::sqrt(B*B - 4*mydet);
                }
                datatype q = 0.5*(B + ( B<0.0 ?-1.0 : 1.0) * discriminant);
                datatype evalue1 = q;
                datatype evalue2 = mydet/q;
                esystem.push_back(SolutionPair<datatype>(evalue1, orig.Rows()));
                datatype ev1 = -lower[k-1]/(diag[k-2] - evalue1);
                datatype norm = 1.0/std::hypot(ev1, 1.0);
                ev1 *= norm;
                esystem.back().evector[k-1] = norm;
                esystem.back().evector[k-2] = ev1;
                esystem.push_back(SolutionPair<datatype>(evalue2, orig.Rows()));
                esystem.back().evector[k-1] = -ev1;
                esystem.back().evector[k-2] =  norm;;
            }
            break;
            default:
	    {
                std::cout<<"looking at "<<min<<" "<<k<<std::endl;
		std::vector<ProblemSet<datatype> > nodes;
                nodes.push_back(createInitialProblem(min, k));
                mrrr(min, k - min, nodes, esystem);
	    }
            }
            min = k;
        }
    }
    std::cout<<"looking at "<<min<<" "<<orig.Rows()<<std::endl;
    if ((orig.Rows() - min) == 1)
    {
        esystem.push_back(SolutionPair<datatype>(diag[orig.Rows()-1], orig.Rows()));
        esystem.back().evector[orig.Rows()-1] = 1.0;
    }
    else
    {
        std::vector<ProblemSet<datatype> > nodes;
	nodes.push_back(createInitialProblem(min, orig.Rows()));
        mrrr(min, orig.Rows() - min, nodes, esystem);//last matrix, or if no degeneracy is present, then the full matrix
    }
    /*std::cout<<"Final Eigenpairs  "<<esystem.size()<<std::endl;
    for(uint k = 0; k < esystem.size(); ++k)
    esystem[k].print();*/
    return esystem;
}

template <class T>
struct Tridiagonalization_opti_trait
{
    static void g_loop ( T& z, uint i )
    {
        for ( uint j = 0; j < i; ++j )
        {
            typename T::value_type g ( 0.0 );
            for ( uint k = 0; k < i; ++k ) //this loop is responsible for 15% of cycles
                g += z ( i,k ) * z ( k,j );
            for ( uint k = 0; k < i; ++k )
            {
                z ( k, j ) -= g * CTDT<typename T::value_type>::conjugate ( z ( k,i ) );
            }
        }
    }
};

template <typename Q>
struct Tridiagonalization_opti_trait<MTLICCL::Matrix< MTLICCL::MemArray1D<MTLICCL::Config1D<Q> > > >
{//yields a speedup of about 25%
    typedef MTLICCL::Matrix< MTLICCL::MemArray1D<MTLICCL::Config1D<Q> > > Type;
    static inline void g_loop ( Type& z, uint i )
    {
        constexpr uint elemsperCL = CLS/sizeof ( Q );
        const uint maxj = ( i / elemsperCL ) * elemsperCL;
        /*        for ( uint j = 0; j < i; ++j )
                {
                    typename Type::value_type g ( 0.0 );
                    asm volatile ( "nop;nop;nop;" );
                    for ( uint k = 0; k < i; ++k ) //this loop is responsible for 15% of cycles
                        g += z ( i,k ) * z ( k,j );
                    asm volatile ( "nop;nop;nop;" );
                    for ( uint k = 0; k < i; ++k )
                    {
                        z ( k, j ) -= g * CTDT<typename Type::value_type>::conjugate ( z ( k,i ) );
                    }
                }*/
        Q* zi_ptr = & z ( i,0 );
        for ( uint j = 0; j < maxj; j += elemsperCL )
        {
            typename Type::value_type g[elemsperCL];
            for ( int t = 0; t < elemsperCL; ++t )
                g[t] = 0.0;
//            asm volatile ( "nop;nop;nop;" );
            for ( uint k = 0; k < i; ++k ) //this loop is responsible for 15% of cycles
	    {
	      Q* zk_ptr = & z(k, j);
                for ( uint t = 0; t < elemsperCL; ++t )
                    g[t] += /*z ( i,k )*/ zi_ptr[k] * /*z ( k,j + t );*/ zk_ptr[t];
	    }
//            asm volatile ( "nop;nop;nop;" );
            for ( uint k = 0; k < i; ++k )
            {
	      Q temp = CTDT<typename Type::value_type>::conjugate ( z ( k,i ) );
	      Q* zk_ptr = & z(k, j);
                for ( uint t = 0; t < elemsperCL; ++t )
                    /*z ( k, j + t )*/ zk_ptr[t] -= g[t] * /*CTDT<typename Type::value_type>::conjugate ( z ( k,i ) )*/ temp;
            }
        }
        for ( uint j = maxj; j < i; ++j )
        {
            typename Type::value_type g ( 0.0 );
            for ( uint k = 0; k < i; ++k ) //this loop is responsible for 15% of cycles
                g += /*z ( i,k )*/ zi_ptr[k] * z ( k,j );
            for ( uint k = 0; k < i; ++k )
            {
                z ( k, j ) -= g * CTDT<typename Type::value_type>::conjugate ( z ( k,i ) );
            }
        }
    }
};

template <typename Q>
struct Tridiagonalization_opti_trait<MTLICCL::Matrix< MTLICCL::Dynamic<MTLICCL::Config_Dynamic<Q> > > >//it is a continuous piece of memory in the second index
{
    typedef MTLICCL::Matrix< MTLICCL::Dynamic<MTLICCL::Config_Dynamic<Q> > > Type;
    static inline void g_loop ( Type& z, uint i )//yields in this version a speedup of 3%
    {
        constexpr uint elemsperCL = CLS/sizeof ( Q );
        const uint maxj = ( i / elemsperCL ) * elemsperCL;
        /*        for ( uint j = 0; j < i; ++j )
                {
                    typename Type::value_type g ( 0.0 );
                    asm volatile ( "nop;nop;nop;" );
                    for ( uint k = 0; k < i; ++k ) //this loop is responsible for 15% of cycles
                        g += z ( i,k ) * z ( k,j );
                    asm volatile ( "nop;nop;nop;" );
                    for ( uint k = 0; k < i; ++k )
                    {
                        z ( k, j ) -= g * CTDT<typename Type::value_type>::conjugate ( z ( k,i ) );
                    }
                }*/
        Q* zi_ptr = & z ( i,0 );
        for ( uint j = 0; j < maxj; j += elemsperCL )
        {
            typename Type::value_type g[elemsperCL];
            for ( int t = 0; t < elemsperCL; ++t )
                g[t] = 0.0;
//            asm volatile ( "nop;nop;nop;" );
            for ( uint k = 0; k < i; ++k ) //this loop is responsible for 15% of cycles
	    {
	      Q* zk_ptr = & z(k, j);
                for ( uint t = 0; t < elemsperCL; ++t )
                    g[t] += /*z ( i,k )*/ zi_ptr[k] * /*z ( k,j + t );*/ zk_ptr[t];
	    }
//            asm volatile ( "nop;nop;nop;" );
            for ( uint k = 0; k < i; ++k )
            {
	      Q temp = CTDT<typename Type::value_type>::conjugate ( z ( k,i ) );
	      Q* zk_ptr = & z(k, j);
                for ( uint t = 0; t < elemsperCL; ++t )
                    /*z ( k, j + t )*/ zk_ptr[t] -= g[t] * /*CTDT<typename Type::value_type>::conjugate ( z ( k,i ) )*/ temp;
            }
        }
        for ( uint j = maxj; j < i; ++j )
        {
            typename Type::value_type g ( 0.0 );
            for ( uint k = 0; k < i; ++k ) //this loop is responsible for 15% of cycles
                g += /*z ( i,k )*/ zi_ptr[k] * z ( k,j );
            for ( uint k = 0; k < i; ++k )
            {
                z ( k, j ) -= g * CTDT<typename Type::value_type>::conjugate ( z ( k,i ) );
            }
        }
    }
};

template <class Matrix>
Matrix EigenSystemSolver<Matrix>::tridiagonalize()
{//implementation adapted from NR for hermitian matrices
    lower[0] = 0;
    CTDT<value_type> ds(diag, lower, orig.Rows());
    Matrix z(orig);
    std::cout<<"calculate Eigenvectors: "<< (calcEVs?"true":"false")<<std::endl;
    for (uint i = orig.Rows() - 1; i > 0; i--)
    {
        int l = i - 1;
        datatype scale = 0.0;
        typename CTDT<value_type>::value_type h = 0.0;
        if (likely(l > 0))
        {
            for (uint k = 0; k < i; ++k)
                scale += CTDT<value_type>::scale_abs(z(i,k));
            if (unlikely(scale == 0.0)) //that's really like that in NR!
                ds.off(i) = CTDT<value_type>::conjugate(z(i,l));
            else
            {
                datatype c = 0.0;
                for (uint k = 0; k < i; ++k)//Kahan summation really needed here! already for a size 64 problem
                {
//                    z(i,k) /= scale;
		    value_type temp = z(i, k);
		    temp/=scale;
		    z(i,k) = temp;
                    datatype y = CTDT<value_type>::abssquared(/*z(i,k)*/temp) - c;
                    datatype t = h + y;
                    c = (t - h) -y;
                    h = t;
//                    h += CTDT<value_type>::abssquared(z(i,k));//z(i,k) * z(i,k);//we form what is called sigma in the explanation text of NR
                }
                value_type f = z(i,l);
                value_type g = ds.get_g(f, h);
                ds.off(i) = scale *CTDT<value_type>::conjugate(g);
                h -= CTDT<value_type>::realpart(CTDT<value_type>::conjugate(g) * f);
                z(i, l) = f-g;
                f = 0.0;
                for (uint j = 0; j < i; ++j)
                {
                    if (calcEVs)
                        z(j, i) = z(i, j)/h;
                    value_type gl = 0.0;//we calculate one entry of A*u and store the result in g
                    value_type cl = 0.0;//gl is evaluated using Kahan's summation formula applied to complex numbers
                    for (uint k = 0; k < j + 1; k++)
                    {
                        value_type y = z(j, k) * CTDT<value_type>::conjugate(z(i, k)) - cl;
                        value_type t = gl + y;
                        cl = (t - gl) - y;
                        gl = t;
//                        gl += z(j, k) * CTDT<value_type>::conjugate(z(i, k));
                    }
                    for (uint k = j + 1; k < i; k++)
                    {
                        value_type y = CTDT<value_type>::conjugate(z(k, j)) * CTDT<value_type>::conjugate(z(i, k)) - cl;
                        value_type t = gl + y;
                        cl = (t - gl) - y;
                        gl = t;
//                        gl += CTDT<value_type>::conjugate(z(k, j)) * CTDT<value_type>::conjugate(z(i, k));
                    }
                    ds.off(j) = gl/h;//we store an element of the vector p to the currently unused entry lower[j]
                    f += ds.off(j) * z(i, j);//in f we accumulate the scalar product u^\dagger * p. u is still stored in z
                }

                //the variable hh denotes the quantity K from NR. K = (u^\dagger * p) / (2 * H)
                //currently p is stored in lowercplx
                //and only the sub diagonal elements of z(i,j) are changed

                typename CTDT<value_type>::value_type hh = CTDT<value_type>::realpart(f) / (h + h);
                for (uint j = 0; j < i; j++)
                {
                    value_type flocal = CTDT<value_type>::conjugate(z(i, j));
                    value_type glocal = ds.off(j) - hh * flocal;
                    ds.off(j) = glocal;
                    for (uint k = 0; k < j + 1; k++)
                        z(j, k) -= (flocal * CTDT<value_type>::conjugate(ds.off(k)) + glocal * z(i,k));//the final reduction A' = A - q*u^\dagger - u * q^\dagger
                }
            }
        }
        else
            ds.off(i) = CTDT<value_type>::conjugate(z(i, l));
        ds.diag(i) = h;
    }
    if (calcEVs)ds.diag(0) = 0.0;
    ds.off(0) = 0.0;
    for (uint i = 0; i < orig.Rows(); ++i)
    {
        if (calcEVs)
        {
            if (ds.diag(i) != 0.0)
            {
	      Tridiagonalization_opti_trait<Matrix>::g_loop(z, i);
             /*   for (uint j = 0; j < i; ++j)
                {
                    value_type g(0.0);
		    asm volatile("nop;nop;nop;");
                    for (uint k = 0; k < i; ++k)//this loop is responsible for 15% of cycles
                        g += z(i,k) * z(k,j);
		    asm volatile("nop;nop;nop;");
                    for (uint k = 0; k < i; ++k)
                    {
                        z(k, j) -= g * CTDT<value_type>::conjugate(z(k,i));
                    }
                }*/
            }
            ds.diag(i) = z(i, i);
            z(i,i) = 1.0;
            for (uint j = 0; j < i; ++j)
            {
                z(j,i) = 0.0;
                z(i,j) = 0.0;
            }
        }
        else
        {
            ds.diag(i) = z(i,i);
        }
    }
    ds.copyover(z, calcEVs);
    /*    Matrix temp(orig.Rows(), orig.Rows());
             for (uint k = 0; k < orig.Rows(); ++k)
            {
                temp(k,k) = diag[k];
                if ((k+1) < orig.Rows())
                {
                    temp(k, k + 1) = lower[k+1];
                    temp(k + 1, k) = lower[k+1];
                }
            }
            std::cout<<temp<<std::endl;/
    //    delete [] diagcplx;
    //    delete [] lowercplx;
            //    std::cout.precision(10);
            std::cout<<z<<std::endl;

                Matrix zhc = ~z;
                for(uint k = 0; k < zhc.Rows(); ++k)
                  for(uint i = 0; i < zhc.Columns(); ++i)
            	zhc(k,i) = conj(zhc(k,i));

                  std::cout<<z * temp * zhc<<std::endl;
                for(uint k = 0; k < orig.Rows(); ++k)
                  std::cout<<diag[k]<<"      "<<lower[k]<<std::endl;
                std::cout<<zhc * orig * z<<std::endl;
                exit(-1);*/
    return z;
}
};
#endif
