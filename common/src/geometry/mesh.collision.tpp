// +-------------------------------------------------------------------------
// | mesh.collision.tpp
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

#include "collide.h"

#include "aabvh.h"


template<class VertData, class TriData>
class Mesh<VertData,TriData>::CollisionProblem 
{
public:
    CollisionProblem(Mesh *owner)
        :   mesh(owner),
            qcoords(owner->verts.size()),
            qvels(owner->verts.size())
    {
        // Callibrate the quantization unit...
        double maxMag = 0.0;
        for(VertData &v : mesh->verts) {
            maxMag = std::max(maxMag, max(abs(v.pos)));
        }
        Quantization::callibrate(maxMag);
        
        // then compute quantized coordinates
        for(uint i=0; i<mesh->verts.size(); i++) {
            Vec3d pos = mesh->verts[i].pos;
            Vec3d vel = mesh->verts[i].velocity;
            for(uint k=0; k<3; k++) {
                qcoords[i][k] = Quantization::quantize(pos[k]);
                qvels[i][k]   = Quantization::quantize(vel[k]);
            }
        }
    }
    
    void perturb()
    {
        const double EPSILON = 1.0e-5; // perturbation epsilon
        for(Vec3d &coord : qcoords) {
            Vec3d perturbation(
                Quantization::quantize(drand(-EPSILON, EPSILON)),
                Quantization::quantize(drand(-EPSILON, EPSILON)),
                Quantization::quantize(drand(-EPSILON, EPSILON))
            );
            coord += perturbation;
        }
        //for(Vec3d &vel : qvels) {
        //    Vec3d perturbation(
        //        Quantization::quantize(drand(-EPSILON, EPSILON)),
        //        Quantization::quantize(drand(-EPSILON, EPSILON)),
        //        Quantization::quantize(drand(-EPSILON, EPSILON))
        //    );
        //    vel += perturbation;
        //}
    }
    
    
    void findCollisions()
    {
        // Uses a REALLY DUMB strategy to explore
        // potential intersections right now!
        //for(uint vi=0; vi<mesh->verts.size(); vi++) {
        //    VertData &v = mesh->verts[vi];
        //    
        //    v.collisions.resize(0);
        //    
        //    for(uint tid=0; tid<mesh->tris.size(); tid++) {
        //        Tri &t = mesh->tris[tid];
        //        
        //        collisionKernel(vi, tid, v, t);
        //    }
        //}
        
        // build an acceleration structure
        std::vector< GeomBlob<uint> > geoms(mesh->verts.size());
        for(uint vi=0; vi<mesh->verts.size(); vi++) {
            VertData &v = mesh->verts[vi];
            v.collisions.resize(0);
            populateBlob(vi, geoms[vi]);
        }
        AABVH<uint> vertBVH(geoms);
        
        // use the acceleration structure
        for(uint tid=0; tid<mesh->tris.size(); tid++) {
            Tri &t = mesh->tris[tid];
            BBox3d bbox = tribox(t);
            vertBVH.for_each_in_box(bbox, [&](uint vi) {
                VertData &v = mesh->verts[vi];
                collisionKernel(vi, tid, v, t);
            });
        }
    }
    
private:
    inline void populateBlob(uint vi, GeomBlob<uint> &blob)
    {
        Vec3d p0 = qcoords[vi];
        Vec3d p1 = p0 + qvels[vi];
        blob.bbox = convex(BBox3d(p0,p0),BBox3d(p1,p1));
        blob.point = (blob.bbox.minp + blob.bbox.maxp) / 2.0;
        blob.id = vi;
    }
    
    inline BBox3d tribox(const Tri &t)
    {
        BBox3d result;
        for(uint k=0; k<3; k++) {
            uint vi = t.v[k];
            Vec3d p0 = qcoords[vi];
            Vec3d p1 = p0 + qvels[vi];
            BBox3d vertbox = convex(BBox3d(p0, p0), BBox3d(p1,p1));
            result = convex(result, vertbox);
        }
        return result;
    }
    
    inline void collisionKernel(
        uint vi, uint tid,
        VertData &v, Tri &t
    ) {
        // filter out combinatorially ridiculous collisions
        if(vi == t.v[0] ||
           vi == t.v[1] ||
           vi == t.v[2])
            return;
        
        // Call the collision routine
        // pack the structures
        VertTriCollisionInput in;
            in.p        =   qcoords[vi];
            in.pv       =   qvels[vi];
        for(uint k=0; k<3; k++) {
            in.t[k]     =   qcoords[t.v[k]];
            in.tv[k]    =   qvels[t.v[k]];
        }
        VertTriCollisionOutput out;
        // call the collision routine
        if(collide(&in, &out)) {
        // and store collisions if we find any!
            for(uint k=0; k<out.nCollisions; k++) {
                double time = out.time[k];
                v.collisions.push_back({time, tid});
            }
        }
    }
    
private:
    Mesh *mesh;
    
    std::vector<Vec3d> qcoords; // quantized coordinates
    std::vector<Vec3d> qvels; // quantized velocities...
};



template<class VertData, class TriData>
void Mesh<VertData,TriData>::findCollisions()
{
    CollisionProblem prob(this);
    
    prob.perturb();
    
    prob.findCollisions();
}
























