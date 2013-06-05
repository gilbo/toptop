// +-------------------------------------------------------------------------
// | collide.cpp
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
#include "collide.h"

#include "cubicRoots.h"

// returns true if any collisions occurred
bool collide(VertTriCollisionInput *in, VertTriCollisionOutput *out)
{
    // First, we re-center the problem on the point in question,
    // and skew temporally so that its velocity is zero-ed out
    Vec3d   t[3];
    Vec3d   tv[3];
    for(uint k=0; k<3; k++) {
        t[k]    = in->t[k] - in->p;
        tv[k]   = in->tv[k] - in->pv;
    }
    
    // we could perform a bounding box check here,
    // though it should probably be performed elsewhere...
    
    // We could also run Tyson's CCD routine here to cull out
    // the majority of intersections.
    
    // now, we want to find times of collision (time)
    // and the barycentric coordinates of collision (bary)
    // To do this, we first solve for time.
    
    // We look at the tetrahedron based at the origin (i.e. point p)
    // and how it changes over time.  In particular,
    // every time the determinant/volume of the tetrahedron zeros out,
    // the origin becomes collineated with the plane of the triangle
    // If that occurs during the time interval 0 to 1,
    // then we have a candidate intersection point.
    // Now, the determinant can be written
    // F(time) = DET( [ t0 ; t1 ; t2 ] + time * [ tv0 ; tv1 ; tv2 ] )
    double  F0  = det(             t[0],             t[1],             t[2]);
    double  F1  = det(     t[0] + tv[0],     t[1] + tv[1],     t[2] + tv[2]);
    double  Fn1 = det(     t[0] - tv[0],     t[1] - tv[1],     t[2] - tv[2]);
    double  F2  = det( t[0] + 2.0*tv[0], t[1] + 2.0*tv[1], t[2] + 2.0*tv[2]);
    // Given these four samples, we can reconstruct an explicit cubic
    // polynomial to represent F(time)
    //      F(time) = A + B*(time) + C*(time)^2 + D*(time)^3
    // which gives us the identities
    //      F0  =  A
    //      F1  =  A +  B +  C +  D
    //      Fn1 =  A -  B +  C -  D
    //      F2  =  A + 2B + 4C + 8D
    // from these, we can get out the coefficients
    double  A   = F0;
    double  X   = (F1 + Fn1) / 2.0;     // = (A + C)
    double  Y   = (F1 - Fn1) / 2.0;     // = (B + D)
    double  C   = X - A;
    double  Z   = F2 - A - 4.0*C;       // = 2B + 8D
    double  D   = (Z - 2.0*Y) / 6.0;
    double  B   = (Z - 8.0*D) / 2.0;
    
    // then we can do an interval search for polynomial roots
    double  times[3];
    uint    nTimes;
    cubicRoots(times, &nTimes, A, B, C, D, 0.0, 1.0);
    ENSURE(nTimes < 3); // should be the case...?
    
    uint write = out->nCollisions = 0;
    for(uint i=0; i<nTimes; i++) {
        double  time    = out->time[write] = times[i];
        Vec3d   v[3];
        for(uint k=0; k<3; k++)
                v[k]    = t[k] + time * tv[k];
        
        // In order to get barycentric coordinates,
        // We measure areas via cross products.
        // In order to check for containment, we
        // test whether these cross products point in the same direction
        Vec3d areas;
        Vec3d crosses[3];
        for(uint k=0; k<3; k++) {
            crosses[k]  = cross(v[(k+1)%3], v[(k+2)%3]);
            areas[k]    = len(crosses[k]);
        }
        bool    contained = true;
        for(uint k=0; k<3; k++) {
            if(dot(crosses[k], crosses[(k+1)%3]) < 0.0)
                contained = false;
        }
        double sumArea = areas[0] + areas[1] + areas[2];
        
        if(contained && sumArea > 0) {
            Vec3d bary          = areas / sumArea;
            out->bary[write]    = bary;
            write++;
        }
    }
    out->nCollisions = write;
    
    return out->nCollisions > 0;
}














