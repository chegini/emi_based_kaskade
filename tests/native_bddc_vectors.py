"""Minimal scalar-vector storage for portable production-codec correctness tests.

These test-only headers provide array/vector arithmetic, not a Dune grid, PDE,
solver, quantizer or encoder. All codec mathematics is the production template.
They must never be included in application builds or performance measurements.
"""
from pathlib import Path


def write_headers(root):
    root=Path(root)
    files={
        'dune/common/fvector.hh':r'''
#pragma once
#include <array>
#include <cmath>
#include <ostream>
namespace Dune {
template<class T,int N> struct FieldVector:std::array<T,N> {
 using field_type=T;
 FieldVector(){this->fill(T{});}
 FieldVector(T v){this->fill(v);}
 FieldVector& operator=(T v){this->fill(v);return *this;}
 operator T&() requires(N==1){return (*this)[0];}
 operator T const&() const requires(N==1){return (*this)[0];}
 FieldVector& operator+=(FieldVector const& b){for(int i=0;i<N;++i)(*this)[i]+=b[i];return *this;}
 FieldVector& operator+=(T b){for(auto& x:*this)x+=b;return *this;}
 FieldVector& operator-=(FieldVector const& b){for(int i=0;i<N;++i)(*this)[i]-=b[i];return *this;}
 FieldVector& operator-=(T b){for(auto& x:*this)x-=b;return *this;}
 FieldVector& operator*=(T a){for(auto& x:*this)x*=a;return *this;}
 FieldVector& operator/=(T a){for(auto& x:*this)x/=a;return *this;}
 void axpy(T a,FieldVector const& b){for(int i=0;i<N;++i)(*this)[i]+=a*b[i];}
 T two_norm2()const{T z=0;for(auto x:*this)z+=x*x;return z;}
 T two_norm()const{return std::sqrt(two_norm2());}
};
template<class T,int N>auto operator+(FieldVector<T,N> a,FieldVector<T,N> const& b){return a+=b;}
template<class T,int N>auto operator-(FieldVector<T,N> a,FieldVector<T,N> const& b){return a-=b;}
template<class T,int N>auto operator*(T a,FieldVector<T,N> b){return b*=a;}
template<class T,int N>auto operator*(FieldVector<T,N> b,T a){return b*=a;}
template<class T,int N>auto operator/(FieldVector<T,N> b,T a){return b/=a;}
template<class T,int N>std::ostream& operator<<(std::ostream& o,FieldVector<T,N> const& v){for(auto x:v)o<<x<<' ';return o;}
}
''',
        'dune/istl/bvector.hh':r'''
#pragma once
#include <vector>
#include <cmath>
namespace Dune {
template<class B>struct BlockVector:std::vector<B>{
 using std::vector<B>::vector;using std::vector<B>::operator=;
 using block_type=B;using field_type=typename B::field_type;
 BlockVector& operator=(field_type a){for(auto& b:*this)b=a;return *this;}
 BlockVector& operator+=(BlockVector const& v){for(std::size_t i=0;i<this->size();++i)(*this)[i]+=v[i];return *this;}
 BlockVector& operator-=(BlockVector const& v){for(std::size_t i=0;i<this->size();++i)(*this)[i]-=v[i];return *this;}
 BlockVector& operator*=(field_type a){for(auto& b:*this)b*=a;return *this;}
 void axpy(field_type a,BlockVector const& v){for(std::size_t i=0;i<this->size();++i)(*this)[i].axpy(a,v[i]);}
 field_type two_norm2()const{field_type a=0;for(auto const& b:*this)a+=b.two_norm2();return a;}
 field_type two_norm()const{return std::sqrt(two_norm2());}
};
}
''',
        # The production transfer includes this header but uses no Boost timers.
        'boost/timer/timer.hpp':'#pragma once\n',
    }
    for name,text in files.items():
        p=root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text)
    return root
