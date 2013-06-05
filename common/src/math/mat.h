// +-------------------------------------------------------------------------
// | mat.h
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2012
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the TopTop library.
// |
// |    TopTop is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    TopTop is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with TopTop.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#pragma once

#include "vec.h"

#include <algorithm>
#include <cmath>
#include <iostream>

// **************************************************************************
// *  Mat2 stores 2x2 matrices using whatever base number type you choose
// **************************************************************************
template<class N>
class Mat2 {
public: // data is stored in row major order
    N data[2*2];
private: // data access macro
    inline N& m(uint r, uint c) { return data[2*r+c]; }
// +---------------------------------
public: // constructors
    inline Mat2() {}
    inline Mat2(const N &m00, const N &m01,
                const N &m10, const N &m11) {
        m(0,0) = m00;  m(0,1) = m01;
        m(1,0) = m10;  m(1,1) = m11;
    }
    template<class T>
    inline Mat2(const Mat2<T> &cp) {
        for(int i=0; i<2*2; i++)
            data[i] = cp.data[i];
    }
    inline Mat2(const N *ar) {
        for(int i=0; i<2*2; i++)
            data[i] = ar[i];
    }
    inline Mat2(const Vec2<N> &r0, const Vec2<N> &r1) {
        for(int c=0; c<2; c++) {
            m(0,c) = r0[c];
            m(1,c) = r1[c];
        }
    }
// +---------------------------------
public: // constants (only implemented for float and double)
    static Mat2<N> zero();
    static Mat2<N> id();
// +---------------------------------
public: // indexing
    inline N& operator()(uint r, uint c)       { return data[2*r+c]; }
    inline N  operator()(uint r, uint c) const { return data[2*r+c]; }
// +---------------------------------
public: // column/row vector extraction
    inline Vec2<N> row(uint r) const
        { return Vec2<N>(data[4*r], data[4*r+1]); }
    inline Vec2<N> col(uint c) const
        { return Vec2<N>(data[c], data[4+c]); }
// +---------------------------------
public: // destructive arithmetic
    inline Mat2<N>& operator+=(const Mat2<N> &rhs);
    inline Mat2<N>& operator-=(const Mat2<N> &rhs);
    inline Mat2<N>& operator*=(const N &rhs);
    inline Mat2<N>& operator/=(const N &rhs);
};
// +---------------------------------
// comparison
template<class N> inline
bool operator==(const Mat2<N> &lhs, const Mat2<N> &rhs);
template<class N> inline
bool operator!=(const Mat2<N> &lhs, const Mat2<N> &rhs);
// +---------------------------------
// non-destructive arithmetic (matrix-as-vector operations)
template<class N> inline
Mat2<N> operator+(const Mat2<N> &lhs, const Mat2<N> &rhs);
template<class N> inline
Mat2<N> operator-(const Mat2<N> &lhs, const Mat2<N> &rhs);
template<class N> inline
Mat2<N> operator*(const Mat2<N> &lhs, const N &rhs);
template<class N> inline
Mat2<N> operator*(const N &lhs, const Mat2<N> &rhs);
template<class N> inline
Mat2<N> operator/(const Mat2<N> &lhs, const N &rhs);
template<class N> inline
Mat2<N> operator-(const Mat2<N> &m);
// +---------------------------------
// non-destructive arithmetic (matrix-specific operations)
template<class N> inline
Vec2<N> operator*(const Mat2<N> &lhs, const Vec2<N> &rhs);
template<class N> inline
Mat2<N> operator*(const Mat2<N> &lhs, const Mat2<N> &rhs);
// +---------------------------------
// named operations (destructive)
template<class N> inline
void transposeMat(Mat2<N> &m);
// +---------------------------------
// named operations (non-destructive)
template<class N> inline
N det(const Mat2<N> &m);
template<class N> inline
N trace(const Mat2<N> &m);
template<class N> inline
Mat2<N> transpose(const Mat2<N> &m);
template<class N> inline
Mat2<N> outer(const Vec2<N> &lhs, const Vec2<N> &rhs);
// +---------------------------------
// NEED SPECIAL EIGENANALYSIS/MATRIX DECOMP FUNCTIONS HERE
// +---------------------------------
template<class N>
std::ostream& operator<<(std::ostream &out, const Mat2<N> &mat);

// **************************************************************************
// *  Mat3 stores 3x3 matrices using whatever base number type you choose
// **************************************************************************
template<class N>
class Mat3 {
public: // data is stored in row major order
    N data[3*3];
private: // data access macro
    inline N& m(uint r, uint c) { return data[3*r+c]; }
// +---------------------------------
public: // indexing
    inline N& operator()(uint r, uint c)       { return data[3*r+c]; }
    inline N  operator()(uint r, uint c) const { return data[3*r+c]; }
// +---------------------------------
public: // column/row vector extraction
    inline Vec3<N> row(uint r) const
        { return Vec3<N>(data[4*r], data[4*r+1], data[4*r+2]); }
    inline Vec3<N> col(uint c) const
        { return Vec3<N>(data[c], data[4+c], data[8+c]); }
// +---------------------------------
public: // constructors
    inline Mat3() {}
    inline Mat3(const N &m00, const N &m01, const N &m02,
                const N &m10, const N &m11, const N &m12,
                const N &m20, const N &m21, const N &m22) {
        m(0,0) = m00;  m(0,1) = m01;  m(0,2) = m02;
        m(1,0) = m10;  m(1,1) = m11;  m(1,2) = m12;
        m(2,0) = m20;  m(2,1) = m21;  m(2,2) = m22;
    }
    template<class T>
    inline Mat3(const Mat3<T> &cp) {
        for(int i=0; i<3*3; i++)
            data[i] = cp.data[i];
    }
    inline Mat3(const N *ar) {
        for(int i=0; i<3*3; i++)
            data[i] = ar[i];
    }
    inline Mat3(const Vec3<N> &r0, const Vec3<N> &r1, const Vec3<N> &r2) {
        for(int c=0; c<3; c++) {
            m(0,c) = r0[c];
            m(1,c) = r1[c];
            m(2,c) = r2[c];
        }
    }
// +---------------------------------
public: // constants (only implemented for float and double)
    static Mat3<N> zero();
    static Mat3<N> id();
// +---------------------------------
public: // destructive arithmetic
    inline Mat3<N>& operator+=(const Mat3<N> &rhs);
    inline Mat3<N>& operator-=(const Mat3<N> &rhs);
    inline Mat3<N>& operator*=(const N &rhs);
    inline Mat3<N>& operator/=(const N &rhs);
};
// +---------------------------------
// comparison
template<class N> inline
bool operator==(const Mat3<N> &lhs, const Mat3<N> &rhs);
template<class N> inline
bool operator!=(const Mat3<N> &lhs, const Mat3<N> &rhs);
// +---------------------------------
// non-destructive arithmetic (matrix-as-vector operations)
template<class N> inline
Mat3<N> operator+(const Mat3<N> &lhs, const Mat3<N> &rhs);
template<class N> inline
Mat3<N> operator-(const Mat3<N> &lhs, const Mat3<N> &rhs);
template<class N> inline
Mat3<N> operator*(const Mat3<N> &lhs, const N &rhs);
template<class N> inline
Mat3<N> operator*(const N &lhs, const Mat3<N> &rhs);
template<class N> inline
Mat3<N> operator/(const Mat3<N> &lhs, const N &rhs);
template<class N> inline
Mat3<N> operator-(const Mat3<N> &m);
// +---------------------------------
// non-destructive arithmetic (matrix-specific operations)
template<class N> inline
Vec3<N> operator*(const Mat3<N> &lhs, const Vec3<N> &rhs);
template<class N> inline
Mat3<N> operator*(const Mat3<N> &lhs, const Mat3<N> &rhs);
// +---------------------------------
// named operations (destructive)
template<class N> inline
void transposeMat(Mat3<N> &m);
// +---------------------------------
// named operations (non-destructive)
template<class N> inline
N det(const Mat3<N> &m);
template<class N> inline
N trace(const Mat3<N> &m);
template<class N> inline
Mat3<N> transpose(const Mat3<N> &m);
template<class N> inline
Mat3<N> outer(const Vec3<N> &lhs, const Vec3<N> &rhs);
// +---------------------------------
// NEED SPECIAL EIGENANALYSIS/MATRIX DECOMP FUNCTIONS HERE
// +---------------------------------
template<class N>
std::ostream& operator<<(std::ostream &out, const Mat3<N> &mat);

// **************************************************************************
// *  Mat4 stores 4x4 matrices using whatever base number type you choose
// **************************************************************************
template<class N>
class Mat4 {
public: // data is stored in row major order
    N data[4*4];
private: // data access macro
    inline N& m(uint r, uint c) { return data[4*r+c]; }
// +---------------------------------
public: // constructors
    inline Mat4() {}
    inline Mat4(const N &m00, const N &m01, const N &m02, const N &m03,
                const N &m10, const N &m11, const N &m12, const N &m13,
                const N &m20, const N &m21, const N &m22, const N &m23,
                const N &m30, const N &m31, const N &m32, const N &m33) {
        m(0,0) = m00;   m(0,1) = m01;   m(0,2) = m02;   m(0,3) = m03;
        m(1,0) = m10;   m(1,1) = m11;   m(1,2) = m12;   m(1,3) = m13;
        m(2,0) = m20;   m(2,1) = m21;   m(2,2) = m22;   m(2,3) = m23;
        m(3,0) = m30;   m(3,1) = m31;   m(3,2) = m32;   m(3,3) = m33;
    }
    template<class T>
    inline Mat4(const Mat4<T> &cp) {
        for(int i=0; i<4*4; i++)
            data[i] = cp.data[i];
    }
    inline Mat4(const N *ar) {
        for(int i=0; i<4*4; i++)
            data[i] = ar[i];
    }
    inline Mat4(const Vec4<N> &r0, const Vec4<N> &r1,
                const Vec4<N> &r2, const Vec4<N> &r3) {
        for(int c=0; c<4; c++) {
            m(0,c) = r0[c];
            m(1,c) = r1[c];
            m(2,c) = r2[c];
            m(3,c) = r3[c];
        }
    }
// +---------------------------------
public: // constants (only implemented for float and double)
    static Mat4<N> zero();
    static Mat4<N> id();
// +---------------------------------
public: // indexing
    inline N& operator()(uint r, uint c)       { return data[4*r+c]; }
    inline N  operator()(uint r, uint c) const { return data[4*r+c]; }
// +---------------------------------
public: // column/row vector extraction
    inline Vec4<N> row(uint r) const
        { return Vec4<N>(data[4*r], data[4*r+1], data[4*r+2], data[4*r+3]); }
    inline Vec4<N> col(uint c) const
        { return Vec4<N>(data[c], data[4+c], data[8+c], data[12+c]); }
// +---------------------------------
public: // destructive arithmetic
    inline Mat4<N>& operator+=(const Mat4<N> &rhs);
    inline Mat4<N>& operator-=(const Mat4<N> &rhs);
    inline Mat4<N>& operator*=(const N &rhs);
    inline Mat4<N>& operator/=(const N &rhs);
};
// +---------------------------------
// comparison
template<class N> inline
bool operator==(const Mat4<N> &lhs, const Mat4<N> &rhs);
template<class N> inline
bool operator!=(const Mat4<N> &lhs, const Mat4<N> &rhs);
// +---------------------------------
// non-destructive arithmetic (matrix-as-vector operations)
template<class N> inline
Mat4<N> operator+(const Mat4<N> &lhs, const Mat4<N> &rhs);
template<class N> inline
Mat4<N> operator-(const Mat4<N> &lhs, const Mat4<N> &rhs);
template<class N> inline
Mat4<N> operator*(const Mat4<N> &lhs, const N &rhs);
template<class N> inline
Mat4<N> operator*(const N &lhs, const Mat4<N> &rhs);
template<class N> inline
Mat4<N> operator/(const Mat4<N> &lhs, const N &rhs);
template<class N> inline
Mat4<N> operator-(const Mat4<N> &m);
// +---------------------------------
// non-destructive arithmetic (matrix-specific operations)
template<class N> inline
Vec4<N> operator*(const Mat4<N> &lhs, const Vec4<N> &rhs);
template<class N> inline
Mat4<N> operator*(const Mat4<N> &lhs, const Mat4<N> &rhs);
// +---------------------------------
// named operations (destructive)
template<class N> inline
void transposeMat(Mat4<N> &m);
// +---------------------------------
// named operations (non-destructive)
template<class N> inline
N det(const Mat4<N> &m);
template<class N> inline
N trace(const Mat4<N> &m);
template<class N> inline
Mat4<N> transpose(const Mat4<N> &m);
template<class N> inline
Mat4<N> outer(const Vec4<N> &lhs, const Vec4<N> &rhs);
// +---------------------------------
// NEED SPECIAL EIGENANALYSIS/MATRIX DECOMP FUNCTIONS HERE
// +---------------------------------
template<class N>
std::ostream& operator<<(std::ostream &out, const Mat4<N> &mat);

// **************************************************************************
// *  Matrix embedding, extraction conversion for different transforms
// **************************************************************************
// NOTE: These are dependent on Mat::zero() and Mat::id() being implemented
// Linear -> Affine
/*template<class N> inline
Mat3<N> affineEmbed(const Mat2<N> &m);
template<class N> inline
Mat3<N> affineEmbed(const Mat2<N> &m, const Vec2<N> &v);
template<class N> inline
Mat4<N> affineEmbed(const Mat3<N> &m);
template<class N> inline
Mat4<N> affineEmbed(const Mat3<N> &m, const Vec3<N> &v);
// +---------------------------------
// Extract Linear portion of Affine Transformation
template<class N> inline
Mat2<N> linearComponent(const Mat3<N> &m);
template<class N> inline
Mat3<N> linearComponent(const Mat4<N> &m);
// +---------------------------------
// Extract Translation portion of Affine Transformation
template<class N> inline
Vec2<N> translationComponent(const Mat3<N> &m);
template<class N> inline
Vec3<N> translationComponent(const Mat4<N> &m);
// +---------------------------------
*/
// **************************************************************************
// *  Matrix generators
// **************************************************************************
// convert (vec x) as an operator on another vector to matrix form
template<class N> inline
Mat3<N> cross(const Vec3<N> &vec);
// rotations
template<class N> inline
Mat2<N> rotByDeg(const N &degree);
template<class N> inline
Mat3<N> rotByDeg(const N &degree, const Vec3<N> &axis);
// provide yaw & pitch in degrees
template<class N> inline
Mat3<N> rotByYawPitch(const N &yaw, const N &pitch);

// **************************************************************************
// *  Aliases for Common Specializations
// **************************************************************************
typedef Mat2<float> Mat2f;
typedef Mat3<float> Mat3f;
typedef Mat4<float> Mat4f;

typedef Mat2<double> Mat2d;
typedef Mat3<double> Mat3d;
typedef Mat4<double> Mat4d;



// **************************************************************************
// **************************************************************************
// ***************************  IMPLEMENTATION  *****************************
// **************************************************************************
// **************************************************************************

// **************************************************************************
// *  Mat2 IMPLEMENTATION
// **************************************************************************
// destructive arithmetic
template<class N> inline
Mat2<N>& Mat2<N>::operator+=(const Mat2<N> &rhs) {
    for(int i=0; i<2*2; i++)
        data[i] += rhs.data[i];
    return *this;
}
template<class N> inline
Mat2<N>& Mat2<N>::operator-=(const Mat2<N> &rhs) {
    for(int i=0; i<2*2; i++)
        data[i] -= rhs.data[i];
    return *this;
}
template<class N> inline
Mat2<N>& Mat2<N>::operator*=(const N &rhs) {
    for(int i=0; i<2*2; i++)
        data[i] *= rhs;
    return *this;
}
template<class N> inline
Mat2<N>& Mat2<N>::operator/=(const N &rhs) {
    for(int i=0; i<2*2; i++)
        data[i] /= rhs;
    return *this;
}
// +---------------------------------
// comparison
template<class N> inline
bool operator==(const Mat2<N> &lhs, const Mat2<N> &rhs) {
    for(int i=0; i<2*2; i++)
        if(lhs.data[i] != rhs.data[i])  return false;
    return true;
}
template<class N> inline
bool operator!=(const Mat2<N> &lhs, const Mat2<N> &rhs) {
    for(int i=0; i<2*2; i++)
        if(lhs.data[i] != rhs.data[i])  return true;
    return false;
}
// +---------------------------------
// non-destructive arithmetic (matrix-as-vector operations)
template<class N> inline
Mat2<N> operator+(const Mat2<N> &lhs, const Mat2<N> &rhs) {
    Mat2<N> res(lhs);
    return res += rhs;
}
template<class N> inline
Mat2<N> operator-(const Mat2<N> &lhs, const Mat2<N> &rhs) {
    Mat2<N> res(lhs);
    return res -= rhs;
}
template<class N> inline
Mat2<N> operator*(const Mat2<N> &lhs, const N &rhs) {
    Mat2<N> res(lhs);
    return res *= rhs;
}
template<class N> inline
Mat2<N> operator*(const N &lhs, const Mat2<N> &rhs) {
    Mat2<N> res(rhs);
    return res *= lhs;
}
template<class N> inline
Mat2<N> operator/(const Mat2<N> &lhs, const N &rhs) {
    Mat2<N> res(lhs);
    return res /= rhs;
}
template<class N> inline
Mat2<N> operator-(const Mat2<N> &m) {
    Mat2<N> res;
    for(int i=0; i<2*2; i++)
        res.data[i] = -m.data[i];
    return res;
}
// +---------------------------------
// non-destructive arithmetic (matrix-specific operations)
template<class N> inline
Vec2<N> operator*(const Mat2<N> &lhs, const Vec2<N> &rhs) {
    Vec2<N> res;
    for(int r=0; r<2; r++)
        res[r] = lhs(r,0)*rhs[0] + lhs(r,1)*rhs[1];
    return res;
}
template<class N> inline
Mat2<N> operator*(const Mat2<N> &lhs, const Mat2<N> &rhs) {
    Mat2<N> res;
    for(int r=0; r<2; r++)
        for(int c=0; c<2; c++)
            res(r,c) = lhs(r,0)*rhs(0,c) + lhs(r,1)*rhs(1,c);
    return res;
}
// +---------------------------------
// named operations (destructive)
template<class N> inline
void transposeMat(Mat2<N> &m) {
    std::swap(m(0,1), m(1,0));
}
// +---------------------------------
// named operations (non-destructive)
template<class N> inline
N det(const Mat2<N> &m) {
    return m(0,0)*m(1,1) - m(0,1)*m(1,0);
}
template<class N> inline
N trace(const Mat2<N> &m) {
    return m(0,0) + m(1,1);
}
template<class N> inline
Mat2<N> transpose(const Mat2<N> &m) {
    Mat2<N> res(m);
    transposeMat(res);
    return res;
}
template<class N> inline
Mat2<N> outer(const Vec2<N> &lhs, const Vec2<N> &rhs) {
    Mat2<N> res;
    for(int r=0; r<2; r++)
        for(int c=0; c<2; c++)
            res(r,c) = lhs[r]*rhs[c];
    return res;
}
// +---------------------------------
// NEED SPECIAL EIGENANALYSIS/MATRIX DECOMP FUNCTIONS HERE
// +---------------------------------
// constant values (implemented only for double and float)
template<> inline
Mat2<float> Mat2<float>::zero() {
    Mat2<float> m;
    for(int i=0; i<2*2; i++)
        m.data[i] = 0;
    return m;
}
template<> inline
Mat2<double> Mat2<double>::zero() {
    Mat2<double> m;
    for(int i=0; i<2*2; i++)
        m.data[i] = 0;
    return m;
}
template<> inline
Mat2<float> Mat2<float>::id() {
    Mat2<float> m;
    for(int r=0; r<2; r++)
        for(int c=0; c<2; c++)
            m.m(r,c) = (r==c)? 1 : 0;
    return m;
}
template<> inline
Mat2<double> Mat2<double>::id() {
    Mat2<double> m;
    for(int r=0; r<2; r++)
        for(int c=0; c<2; c++)
            m.m(r,c) = (r==c)? 1 : 0;
    return m;
}

template<class N>
std::ostream& operator<<(std::ostream &out, const Mat2<N> &mat) {
    return out << "[[" << mat(0,0) << ',' << mat(0,1) << "];["
                       << mat(1,0) << ',' << mat(1,1) << "]]";
}

// **************************************************************************
// *  Mat3 IMPLEMENTATION
// **************************************************************************
// destructive arithmetic
template<class N> inline
Mat3<N>& Mat3<N>::operator+=(const Mat3<N> &rhs) {
    for(int i=0; i<3*3; i++)
        data[i] += rhs.data[i];
    return *this;
}
template<class N> inline
Mat3<N>& Mat3<N>::operator-=(const Mat3<N> &rhs) {
    for(int i=0; i<3*3; i++)
        data[i] -= rhs.data[i];
    return *this;
}
template<class N> inline
Mat3<N>& Mat3<N>::operator*=(const N &rhs) {
    for(int i=0; i<3*3; i++)
        data[i] *= rhs;
    return *this;
}
template<class N> inline
Mat3<N>& Mat3<N>::operator/=(const N &rhs) {
    for(int i=0; i<3*3; i++)
        data[i] /= rhs;
    return *this;
}
// +---------------------------------
// comparison
template<class N> inline
bool operator==(const Mat3<N> &lhs, const Mat3<N> &rhs) {
    for(int i=0; i<3*3; i++)
        if(lhs.data[i] != rhs.data[i])  return false;
    return true;
}
template<class N> inline
bool operator!=(const Mat3<N> &lhs, const Mat3<N> &rhs) {
    for(int i=0; i<3*3; i++)
        if(lhs.data[i] != rhs.data[i])  return true;
    return false;
}
// +---------------------------------
// non-destructive arithmetic (matrix-as-vector operations)
template<class N> inline
Mat3<N> operator+(const Mat3<N> &lhs, const Mat3<N> &rhs) {
    Mat3<N> res(lhs);
    return res += rhs;
}
template<class N> inline
Mat3<N> operator-(const Mat3<N> &lhs, const Mat3<N> &rhs) {
    Mat3<N> res(lhs);
    return res -= rhs;
}
template<class N> inline
Mat3<N> operator*(const Mat3<N> &lhs, const N &rhs) {
    Mat3<N> res(lhs);
    return res *= rhs;
}
template<class N> inline
Mat3<N> operator*(const N &lhs, const Mat3<N> &rhs) {
    Mat3<N> res(rhs);
    return res *= lhs;
}
template<class N> inline
Mat3<N> operator/(const Mat3<N> &lhs, const N &rhs) {
    Mat3<N> res(lhs);
    return res /= rhs;
}
template<class N> inline
Mat3<N> operator-(const Mat3<N> &m) {
    Mat3<N> res;
    for(int i=0; i<3*3; i++)
        res.data[i] = -m.data[i];
    return res;
}
// +---------------------------------
// non-destructive arithmetic (matrix-specific operations)
template<class N> inline
Vec3<N> operator*(const Mat3<N> &lhs, const Vec3<N> &rhs) {
    Vec3<N> res;
    for(int r=0; r<3; r++)
        res[r] = lhs(r,0)*rhs[0] + lhs(r,1)*rhs[1] + lhs(r,2)*rhs[2];
    return res;
}
template<class N> inline
Mat3<N> operator*(const Mat3<N> &lhs, const Mat3<N> &rhs) {
    Mat3<N> res;
    for(int r=0; r<3; r++)
        for(int c=0; c<3; c++)
            res(r,c) = lhs(r,0)*rhs(0,c) + lhs(r,1)*rhs(1,c)
                     + lhs(r,2)*rhs(2,c);
    return res;
}
// +---------------------------------
// named operations (destructive)
template<class N> inline
void transposeMat(Mat3<N> &m) {
    std::swap(m(0,1), m(1,0));
    std::swap(m(0,2), m(2,0));
    std::swap(m(1,2), m(2,1));
}
// +---------------------------------
// named operations (non-destructive)
template<class N> inline
N det(const Mat3<N> &m) {
    N d01 = m(0,0)*m(1,1) - m(0,1)*m(1,0);
    N d02 = m(0,0)*m(1,2) - m(0,2)*m(1,0);
    N d12 = m(0,1)*m(1,2) - m(0,2)*m(1,1);
    return d01*m(2,2) - d02*m(2,1) + d12*m(2,0);
}
template<class N> inline
N trace(const Mat3<N> &m) {
    return m(0,0) + m(1,1) + m(2,2);
}
template<class N> inline
Mat3<N> transpose(const Mat3<N> &m) {
    Mat3<N> res(m);
    transposeMat(res);
    return res;
}
template<class N> inline
Mat3<N> outer(const Vec3<N> &lhs, const Vec3<N> &rhs) {
    Mat3<N> res;
    for(int r=0; r<3; r++)
        for(int c=0; c<3; c++)
            res(r,c) = lhs[r]*rhs[c];
    return res;
}
// +---------------------------------
// NEED SPECIAL EIGENANALYSIS/MATRIX DECOMP FUNCTIONS HERE
// +---------------------------------
// constant values (implemented only for double and float)
template<> inline
Mat3<float> Mat3<float>::zero() {
    Mat3<float> m;
    for(int i=0; i<3*3; i++)
        m.data[i] = 0;
    return m;
}
template<> inline
Mat3<double> Mat3<double>::zero() {
    Mat3<double> m;
    for(int i=0; i<3*3; i++)
        m.data[i] = 0;
    return m;
}
template<> inline
Mat3<float> Mat3<float>::id() {
    Mat3<float> m;
    for(int r=0; r<3; r++)
        for(int c=0; c<3; c++)
            m.m(r,c) = (r==c)? 1 : 0;
    return m;
}
template<> inline
Mat3<double> Mat3<double>::id() {
    Mat3<double> m;
    for(int r=0; r<3; r++)
        for(int c=0; c<3; c++)
            m.m(r,c) = (r==c)? 1 : 0;
    return m;
}

template<class N>
std::ostream& operator<<(std::ostream &out, const Mat3<N> &mat) {
    out << '[';
    for(uint r=0; r<3; r++) {
        out << '[';
        for(uint c=0; c<3; c++) {
            out << mat(r,c);
            if(c < 2) out << ',';
        }
        out << ']';
        if(r < 2) out << ';';
    }
    return out << ']';
}

// **************************************************************************
// *  Mat4 IMPLEMENTATION
// **************************************************************************
// destructive arithmetic
template<class N> inline
Mat4<N>& Mat4<N>::operator+=(const Mat4<N> &rhs) {
    for(int i=0; i<4*4; i++)
        data[i] += rhs.data[i];
    return *this;
}
template<class N> inline
Mat4<N>& Mat4<N>::operator-=(const Mat4<N> &rhs) {
    for(int i=0; i<4*4; i++)
        data[i] -= rhs.data[i];
    return *this;
}
template<class N> inline
Mat4<N>& Mat4<N>::operator*=(const N &rhs) {
    for(int i=0; i<4*4; i++)
        data[i] *= rhs;
    return *this;
}
template<class N> inline
Mat4<N>& Mat4<N>::operator/=(const N &rhs) {
    for(int i=0; i<4*4; i++)
        data[i] /= rhs;
    return *this;
}
// +---------------------------------
// comparison
template<class N> inline
bool operator==(const Mat4<N> &lhs, const Mat4<N> &rhs) {
    for(int i=0; i<4*4; i++)
        if(lhs.data[i] != rhs.data[i])  return false;
    return true;
}
template<class N> inline
bool operator!=(const Mat4<N> &lhs, const Mat4<N> &rhs) {
    for(int i=0; i<4*4; i++)
        if(lhs.data[i] != rhs.data[i])  return true;
    return false;
}
// +---------------------------------
// non-destructive arithmetic (matrix-as-vector operations)
template<class N> inline
Mat4<N> operator+(const Mat4<N> &lhs, const Mat4<N> &rhs) {
    Mat4<N> res(lhs);
    return res += rhs;
}
template<class N> inline
Mat4<N> operator-(const Mat4<N> &lhs, const Mat4<N> &rhs) {
    Mat4<N> res(lhs);
    return res -= rhs;
}
template<class N> inline
Mat4<N> operator*(const Mat4<N> &lhs, const N &rhs) {
    Mat4<N> res(lhs);
    return res *= rhs;
}
template<class N> inline
Mat4<N> operator*(const N &lhs, const Mat4<N> &rhs) {
    Mat4<N> res(rhs);
    return res *= lhs;
}
template<class N> inline
Mat4<N> operator/(const Mat4<N> &lhs, const N &rhs) {
    Mat4<N> res(lhs);
    return res /= rhs;
}
template<class N> inline
Mat4<N> operator-(const Mat4<N> &m) {
    Mat4<N> res;
    for(int i=0; i<4*4; i++)
        res.data[i] = -m.data[i];
    return res;
}
// +---------------------------------
// non-destructive arithmetic (matrix-specific operations)
template<class N> inline
Vec4<N> operator*(const Mat4<N> &lhs, const Vec4<N> &rhs) {
    Vec4<N> res;
    for(int r=0; r<4; r++)
        res[r] = lhs(r,0)*rhs[0] + lhs(r,1)*rhs[1]
               + lhs(r,2)*rhs[2] + lhs(r,3)*rhs[3];
    return res;
}
template<class N> inline
Mat4<N> operator*(const Mat4<N> &lhs, const Mat4<N> &rhs) {
    Mat4<N> res;
    for(int r=0; r<4; r++)
        for(int c=0; c<4; c++)
            res(r,c) = lhs(r,0)*rhs(0,c) + lhs(r,1)*rhs(1,c)
                     + lhs(r,2)*rhs(2,c) + lhs(r,3)*rhs(3,c);
    return res;
}
// +---------------------------------
// named operations (destructive)
template<class N> inline
void transposeMat(Mat4<N> &m) {
    for(int r=0; r<4; r++)
        for(int c=r+1; c<4; c++)
            std::swap(m(r,c),m(c,r));
}
// +---------------------------------
// named operations (non-destructive)
template<class N> inline
N det4helper(uint r0, uint r1, const Mat4<N> &m, uint c0, uint c1) {
    return m(r0,c0) * m(r1,c1) - m(r0,c1) * m(r1,c0);
}
template<class N> inline
N det(const Mat4<N> &m) {
    return det4helper(0,1,m,0,1) * det4helper(2,3,m,2,3)
         - det4helper(0,1,m,0,2) * det4helper(2,3,m,1,3)
         + det4helper(0,1,m,0,3) * det4helper(2,3,m,1,2)
         + det4helper(0,1,m,1,2) * det4helper(2,3,m,0,3)
         - det4helper(0,1,m,1,3) * det4helper(2,3,m,0,2)
         + det4helper(0,1,m,2,3) * det4helper(2,3,m,0,1);
}
template<class N> inline
N trace(const Mat4<N> &m) {
    return m(0,0) + m(1,1) + m(2,2) + m(3,3);
}
template<class N> inline
Mat4<N> transpose(const Mat4<N> &m) {
    Mat4<N> res(m);
    transposeMat(res);
    return res;
}
template<class N> inline
Mat4<N> outer(const Vec4<N> &lhs, const Vec4<N> &rhs) {
    Mat4<N> res;
    for(int r=0; r<4; r++)
        for(int c=0; c<4; c++)
            res(r,c) = lhs[r]*rhs[c];
    return res;
}
// +---------------------------------
// NEED SPECIAL EIGENANALYSIS/MATRIX DECOMP FUNCTIONS HERE
// +---------------------------------
// constant values (implemented only for double and float)
template<> inline
Mat4<float> Mat4<float>::zero() {
    Mat4<float> m;
    for(int i=0; i<4*4; i++)
        m.data[i] = 0;
    return m;
}
template<> inline
Mat4<double> Mat4<double>::zero() {
    Mat4<double> m;
    for(int i=0; i<4*4; i++)
        m.data[i] = 0;
    return m;
}
template<> inline
Mat4<float> Mat4<float>::id() {
    Mat4<float> m;
    for(int r=0; r<4; r++)
        for(int c=0; c<4; c++)
            m.m(r,c) = (r==c)? 1 : 0;
    return m;
}
template<> inline
Mat4<double> Mat4<double>::id() {
    Mat4<double> m;
    for(int r=0; r<4; r++)
        for(int c=0; c<4; c++)
            m.m(r,c) = (r==c)? 1 : 0;
    return m;
}

template<class N>
std::ostream& operator<<(std::ostream &out, const Mat4<N> &mat) {
    out << '[';
    for(uint r=0; r<4; r++) {
        out << '[';
        for(uint c=0; c<4; c++) {
            out << mat(r,c);
            if(c < 3) out << ',';
        }
        out << ']';
        if(r < 3) out << ';';
    }
    return out << ']';
}

// **************************************************************************
// *  Matrix embedding, extraction conversion IMPLEMENTAION
// **************************************************************************
// Linear -> Affine
/*template<class N> inline
Mat3<N> affineEmbed(const Mat2<N> &m) {
    Mat3<N> res = Mat3<N>::id();
    for(uint r=0; r<2; r++)
        for(uint c=0; c<2; c++)
            res(r,c) = m(r,c);
    return res;
}
template<class N> inline
Mat3<N> affineEmbed(const Mat2<N> &m, const Vec2<N> &v) {
    Mat3<N> res = Mat3<N>::id();
    for(uint r=0; r<2; r++)
        for(uint c=0; c<2; c++)
            res(r,c) = m(r,c);
    for(uint k=0; k<2; k++)
        res(k,2) = v[k];
    return res;
}
template<class N> inline
Mat4<N> affineEmbed(const Mat3<N> &m) {
    Mat4<N> res = Mat4<N>::id();
    for(uint r=0; r<3; r++)
        for(uint c=0; c<3; c++)
            res(r,c) = m(r,c);
    return res;
}
template<class N> inline
Mat4<N> affineEmbed(const Mat3<N> &m, const Vec3<N> &v) {
    Mat4<N> res = Mat4<N>::id();
    for(uint r=0; r<3; r++)
        for(uint c=0; c<3; c++)
            res(r,c) = m(r,c);
    for(uint k=0; k<3; k++)
        res(k,3) = v[k];
    return res;
}
// +---------------------------------
// Extract Linear portion of Affine Transformation
template<class N> inline
Mat2<N> linearComponent(const Mat3<N> &m) {
    Mat2<N> res;
    for(uint r=0; r<2; r++)
        for(uint c=0; c<2; c++)
            res(r,c) = m(r,c);
    return res;
}
template<class N> inline
Mat3<N> linearComponent(const Mat4<N> &m) {
    Mat3<N> res;
    for(uint r=0; r<3; r++)
        for(uint c=0; c<3; c++)
            res(r,c) = m(r,c);
    return res;
}
// +---------------------------------
// Extract Translation portion of Affine Transformation
template<class N> inline
Vec2<N> translationComponent(const Mat3<N> &m) {
    return Vec2<N>(m(2,0), m(2,1));
}
template<class N> inline
Vec3<N> translationCompooent(const Mat4<N> &m) {
    return Vec3<N>(m(2,0), m(2,1), m(2,2));
}
// +---------------------------------
*/
// **************************************************************************
// *  Matrix generator IMPLEMENTATION
// **************************************************************************
template<> inline
Mat3d cross(const Vec3d &vec)
{
    return Mat3d(     0, -vec.z,  vec.y,
                  vec.z,      0, -vec.x,
                 -vec.y,  vec.x,      0);
}
template<> inline
Mat3f cross(const Vec3f &vec)
{
    return Mat3f(     0, -vec.z,  vec.y,
                  vec.z,      0, -vec.x,
                 -vec.y,  vec.x,      0);
}
template<> inline
Mat2d rotByDeg(const double &degree)
{
    using std::cos;     using std::sin;
    double cosd = cos(deg2rad(degree));
    double sind = sin(deg2rad(degree));
    
    return Mat2d(cosd, -sind,
                 sind,  cosd);
}
template<> inline
Mat3d rotByDeg(const double &degree, const Vec3d &axis)
{
    using std::cos;     using std::sin;
    double cosd = cos(deg2rad(degree));
    double sind = sin(deg2rad(degree));
    
    // from wikipedia article on axis-angle representation
    Vec3d n = normalized(axis);
    return cosd * Mat3d::id()
         + sind * cross(n)
         + (1.0-cosd) * outer(n,n);
}
template<> inline
Mat3f rotByDeg(const float &degree, const Vec3f &axis)
{
    using std::cos;     using std::sin;
    float cosd = cos(deg2rad(degree));
    float sind = sin(deg2rad(degree));
    
    // from wikipedia article on axis-angle representation
    Vec3f n = normalized(axis);
    return cosd * Mat3f::id()
         + sind * cross(n)
         + (1.0f-cosd) * outer(n,n);
}
// provide yaw & pitch in degrees
template<> inline
Mat3d rotByYawPitch(const double &yaw, const double &pitch)
{
    using std::cos;     using std::sin;
    double cosy = cos(deg2rad(yaw));
    double siny = sin(deg2rad(yaw));
    double cosp = cos(deg2rad(pitch));
    double sinp = sin(deg2rad(pitch));
    /*  
     *  [  cosy 0 siny ] [ 1    0     0 ]   [  cosy siny*sinp siny*cosp ]
     *  [     0 1    0 ] [ 0 cosp -sinp ] = [     0      cosp     -sinp ]
     *  [ -siny 0 cosy ] [ 0 sinp  cosp ]   [ -siny cosy*sinp cosy*cosp ]
     */
    return Mat3d(  cosy,  siny * sinp,  siny * cosp,
                      0,         cosp,        -sinp,
                  -siny,  cosy * sinp,  cosy * cosp);
    /*  
     *  [ 1    0     0 ] [  cosy 0 siny ]   [       cosy    0       siny ]
     *  [ 0 cosp -sinp ] [     0 1    0 ] = [  siny*sinp cosp -sinp*cosy ]
     *  [ 0 sinp  cosp ] [ -siny 0 cosy ]   [ -siny*cosp sinp  cosy*cosp ]
     */
     /*return Mat3d(        cosy,     0,         siny,
                   siny * sinp,  cosp, -sinp * cosy,
                  -siny * cosp,  sinp,  cosy * cosp);*/
}
template<> inline
Mat3f rotByYawPitch(const float &yaw, const float &pitch)
{
    using std::cos;     using std::sin;
    float cosy = cos(deg2rad(yaw));
    float siny = sin(deg2rad(yaw));
    float cosp = cos(deg2rad(pitch));
    float sinp = sin(deg2rad(pitch));
    /*  
     *  [  cosy 0 siny ] [ 1    0     0 ]   [  cosy siny*sinp siny*cosp ]
     *  [     0 1    0 ] [ 0 cosp -sinp ] = [     0      cosp     -sinp ]
     *  [ -siny 0 cosy ] [ 0 sinp  cosp ]   [ -siny cosy*sinp cosy*cosp ]
     */
    return Mat3f(  cosy,  siny * sinp,  siny * cosp,
                      0,         cosp,        -sinp,
                  -siny,  cosy * sinp,  cosy * cosp);
    /*  
     *  [ 1    0     0 ] [  cosy 0 siny ]   [       cosy    0       siny ]
     *  [ 0 cosp -sinp ] [     0 1    0 ] = [  siny*sinp cosp -sinp*cosy ]
     *  [ 0 sinp  cosp ] [ -siny 0 cosy ]   [ -siny*cosp sinp  cosy*cosp ]
     */
     /*return Mat3d(        cosy,     0,         siny,
                   siny * sinp,  cosp, -sinp * cosy,
                  -siny * cosp,  sinp,  cosy * cosp);*/
}


