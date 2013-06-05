// +-------------------------------------------------------------------------
// | meshTypes.h
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

#include "mesh.decl.h"

#include "drawMesh.h"

struct ModelerTriangle;

struct ModelerVertex :
    public MinimalVertexData,
    public RenderVertexData,
    public DistanceVertexData,
    public RemeshVertexData,
    public IsctVertexData,
    public CollisionVertexData
{
    void merge(const ModelerVertex &v0, const ModelerVertex &v1) {
        double                              a0 = 0.5;
        if(v0.manifold && !v1.manifold)     a0 = 0.0;
        if(!v0.manifold && v1.manifold)     a0 = 1.0;
        double a1 = 1.0 - a0;
        
        pos         = a0 * v0.pos       + a1 * v1.pos;
        distField   = a0 * v0.distField + a1 * v1.distField;
        velocity    = a0 * v0.velocity  + a1 * v1.velocity;
        dead = UNKNOWN;
        color       = float(a0) * v0.color
                    + float(a1) * v1.color;
        pastCollisionParity =
            (v0.pastCollisionParity && v1.pastCollisionParity)? 1 : 0;
    }
    void interpolate(const ModelerVertex &v0, const ModelerVertex &v1) {
        double a0 = 0.5;
        double a1 = 0.5;
        pos         = a0 * v0.pos       + a1 * v1.pos;
        distField   = a0 * v0.distField + a1 * v1.distField;
        velocity    = a0 * v0.velocity  + a1 * v1.velocity;
        dead = UNKNOWN;
        color       = float(a0) * v0.color
                    + float(a1) * v1.color;
        pastCollisionParity =
            (v0.pastCollisionParity && v1.pastCollisionParity)? 1 : 0;
    }
    
    
    void isct(IsctVertEdgeTriInput<ModelerVertex,ModelerTriangle> input)
    {
        Vec2d       a_e     = Vec2d(1,1)/2.0;
        Vec3d       a_t     = Vec3d(1,1,1)/3.0;
        a_e /= 2.0;
        a_t /= 2.0;
        color       = Vec4f(0,0,0,0);
        distField   = 0.0;
        velocity    = Vec3d(0,0,0);
        dead        = UNKNOWN;
        for(uint k=0; k<2; k++) {
            color       += float(a_e.v[k])  * input.e[k]->color;
            distField   += a_e.v[k]         * input.e[k]->distField;
            velocity    += a_e.v[k]         * input.e[k]->velocity;
        }
        for(uint k=0; k<3; k++) {
            color       += float(a_t.v[k])  * input.t[k]->color;
            distField   += a_t.v[k]         * input.t[k]->distField;
            velocity    += a_t.v[k]         * input.t[k]->velocity;
        }
    }
    void isct(IsctVertTriTriTriInput<ModelerVertex,ModelerTriangle> input)
    {
        Vec3d       a[3];
        for(uint k=0; k<3; k++) {
            a[k]    = Vec3d(1,1,1)/3.0;
            a[k] /= 3.0;
        }
        color       = Vec4f(0,0,0,0);
        distField   = 0.0;
        dead        = UNKNOWN;
        for(uint i=0; i<3; i++) {
          for(uint j=0; j<3; j++) {
            color       += float(a[i].v[j]) * input.t[i][j]->color;
            distField   += a[i].v[j]        * input.t[i][j]->distField;
            velocity    += a[i].v[j]        * input.t[i][j]->velocity;
        }}
    }
    void isctInterpolate(const ModelerVertex &v0, const ModelerVertex &v1) {
        double a0 = len(v1.pos - pos);
        double a1 = len(v0.pos - pos);
        if(a0 + a1 == 0.0) a0 = a1 = 0.5; // safety
        double sum = a0+a1;
        a0 /= sum;
        a1 /= sum;
        
        distField   = a0 * v0.distField + a1 * v1.distField;
        velocity    = a0 * v0.velocity  + a1 * v1.velocity;
        dead = UNKNOWN;
        color       = float(a0) * v0.color
                    + float(a1) * v1.color;
        pastCollisionParity =
            (v0.pastCollisionParity && v1.pastCollisionParity)? 1 : 0;
    }
};
struct ModelerTriangle :
    public MinimalTriangleData,
    public RenderTriangleData,
    public DistanceTriangleData,
    public RemeshTriangleData,
    public IsctTriangleData,
    public CollisionTriangleData
{
    void merge(const ModelerTriangle &, const ModelerTriangle &) {}
    static void split(ModelerTriangle &, ModelerTriangle &,
                      const ModelerTriangle &) {}
    void move(const ModelerTriangle &) {}
    void subdivide(SubdivideTriInput<ModelerVertex,ModelerTriangle>) {}
};
using RawModelerMesh = RawMesh<ModelerVertex, ModelerTriangle>;
using ModelerMesh = Mesh<ModelerVertex, ModelerTriangle>;





