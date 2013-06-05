// +-------------------------------------------------------------------------
// | transform.h
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

#include "mat.h"
#include "ray.h"

// affine transforms
class Transform
{
public:
    Transform() : forward(Mat4f::id()), backward(Mat4f::id()) {}
    
    Mat4f mat() const { return forward; }
    Mat4f invMat() const { return backward; }
    
    Vec3f translate() const {
        return Vec3f(forward(0, 3), forward(1, 3), forward(2, 3));
    }
    Vec3f invTranslate() const {
        return Vec3f(forward(0, 3), forward(1, 3), forward(2, 3));
    }
    // note: rotate may contain information which is not strictly rotation
    Mat3f rot() const {
        return Mat3f(forward(0,0), forward(0,1), forward(0,2),
                     forward(1,0), forward(1,1), forward(1,2),
                     forward(2,0), forward(2,1), forward(2,2));
    }
    Mat3f invRot() const { // may also include scale etc.
        return Mat3f(backward(0,0), backward(0,1), backward(0,2),
                     backward(1,0), backward(1,1), backward(1,2),
                     backward(2,0), backward(2,1), backward(2,2));
    }
    
    Transform inverse() const {
        Transform ret;
        ret.forward = backward;
        ret.backward = forward;
        return ret;
    }
    
    static Transform translation(const Vec3f &vec) {
        Transform ret;
            ret.forward(0,3) = vec.x;
            ret.forward(1,3) = vec.y;
            ret.forward(2,3) = vec.z;
            ret.backward(0,3) = -vec.x;
            ret.backward(1,3) = -vec.y;
            ret.backward(2,3) = -vec.z;
        return ret;
    }
    
    static Transform rotation(const Mat3f &rot, const Mat3f &invRot) {
        Transform ret;
        for(uint r=0; r<3; r++) {
            for(uint c=0; c<3; c++) {
                ret.forward(r,c) = rot(r,c);
                ret.backward(r,c) = invRot(r,c);
            }
        }
        return ret;
    }
    
    Transform operator*(const Transform &rhs) const {
        Transform ret;
        ret.forward = forward * rhs.forward;
        ret.backward = rhs.backward * backward;
        return ret;
    }
    
    Vec3f operator()(const Vec3f &vec) const {
        return rot() * vec + translate();
    }
    
    Ray3f operator()(const Ray3f &ray) const {
        return Ray3f((*this)(ray.p), rot() * ray.r);
    }
    
private:
    Mat4f forward;
    Mat4f backward;
};










