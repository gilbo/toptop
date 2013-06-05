// +-------------------------------------------------------------------------
// | mesh.tpp
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

#include <algorithm>
#include "unsafeRayTriIsct.h"
#include <cfloat>
#include <cmath>
#include <sstream>

#include "unionFind.h"

// constructors
template<class VertData, class TriData>
Mesh<VertData,TriData>::Mesh() {}
template<class VertData, class TriData>
Mesh<VertData,TriData>::Mesh(Mesh &&cp)
    : tris(cp.tris), verts(cp.verts), surfps(std::move(cp.surfps))
{}
template<class VertData, class TriData>
Mesh<VertData,TriData>::Mesh(const RawMesh<VertData,TriData> &raw) :
    tris(raw.triangles.size()), verts(raw.vertices)
{
    // fill out the triangles
    for(uint i=0; i<raw.triangles.size(); i++) {
        tris[i].data = raw.triangles[i];
        tris[i].a = raw.triangles[i].a;
        tris[i].b = raw.triangles[i].b;
        tris[i].c = raw.triangles[i].c;
    }
}
template<class VertData, class TriData>
Mesh<VertData,TriData>::~Mesh() {}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::operator=(Mesh &&src)
{
    tris = src.tris;
    verts = src.verts;
    surfps = std::move(src.surfps);
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::valid() const
{
    for(uint i=0; i<verts.size(); i++) {
        if(!std::isfinite(verts[i].pos.x) ||
           !std::isfinite(verts[i].pos.y) ||
           !std::isfinite(verts[i].pos.z)) {
            std::ostringstream message;
            message << "vertex #" << i << " has non-finite coordinates: "
                    << verts[i].pos;
            ERROR(message.str());
            return false;
        }
    }
    
    for(uint i=0; i<tris.size(); i++) {
        if(tris[i].a >= verts.size() ||
           tris[i].b >= verts.size() ||
           tris[i].c >= verts.size()) {
            std::ostringstream message;
            message << "triangle #" << i << " should have indices in "
                    << "the range 0 to " << (verts.size()-1)
                    << ", but it has invalid indices: "
                    << tris[i].a << ", " << tris[i].b << ", " << tris[i].c;
            ERROR(message.str());
            return false;
        }
    }
    
    return true;
}

template<class VertData, class TriData>
RawMesh<VertData,TriData> Mesh<VertData,TriData>::raw() const
{
    RawMesh<VertData,TriData> result;
    result.vertices = verts;
    result.triangles.resize(tris.size());
    for(uint i=0; i<tris.size(); i++) {
        result.triangles[i]   = tris[i].data;
        result.triangles[i].a = tris[i].a;
        result.triangles[i].b = tris[i].b;
        result.triangles[i].c = tris[i].c;
    }
    return result;
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::disjointUnion(const Mesh &cp)
{
    uint oldVsize = verts.size();
    uint oldTsize = tris.size();
    uint cpVsize  = cp.verts.size();
    uint cpTsize  = cp.tris.size();
    uint newVsize = oldVsize + cpVsize;
    uint newTsize = oldTsize + cpTsize;
    
    std::vector<int> v_remap(cpVsize); // oh this is obvious...
    verts.resize(newVsize);
    tris.resize(newTsize);
    
    for(uint i=0; i<cpVsize; i++)
        verts[oldVsize + i] = cp.verts[i];
    
    for(uint i=0; i<cpTsize; i++) {
        auto &tri = tris[oldTsize + i];
        tri = cp.tris[i];
        tri.a += oldVsize;
        tri.b += oldVsize;
        tri.c += oldVsize;
    }
}





// Picking.
// Dumb Implementation just passes over all triangles w/o any precomputed
// acceleration structure
template<class VertData, class TriData>
typename Mesh<VertData,TriData>::Isct
    Mesh<VertData,TriData>::pick(Ray3d ray)
{
    Isct result;
    result.ray = ray;
    result.exists = false;
    
    double mint = DBL_MAX;
    
    // pass all triangles over ray
    for(uint i=0; i<tris.size(); i++) {
        const Tri  &tri = tris[i];
        
        uint   a = tri.a;
        uint   b = tri.b;
        uint   c = tri.c;
        Vec3d va = verts[a].pos;
        Vec3d vb = verts[b].pos;
        Vec3d vc = verts[c].pos;
        // normalize vertex order (to prevent leaks)
        if(a > b) { std::swap(a, b); std::swap(va, vb); }
        if(b > c) { std::swap(b, c); std::swap(vb, vc); }
        if(a > b) { std::swap(a, b); std::swap(va, vb); }
        
        double t;
        Vec3d  bary;
        if(isct_ray_triangle(ray, va, vb, vc, &t, &bary)) {
            if(t > 0 && t < mint) {
                result.exists = true;
                mint = t;
                result.tri_id = i;
                result.isct = ray.p + t * ray.r;
                result.bary = bary;
            }
        }
    }
    
    return result;
}




static inline
bool contains(const ShortVec<uint, 8> &list, uint item)
{
    for(uint k : list)
        if(k == item)
            return true;
    return false;
}

template<class VertData, class TriData>
typename Mesh<VertData,TriData>::NeighborCache
    Mesh<VertData,TriData>::createNeighborCache()
{
    NeighborCache result;
    result.skeleton.resize(verts.size());
    
    for(uint tid = 0; tid < tris.size(); tid++) {
        const Tri &tri = tris[tid];
        
        result(tri.a, tri.b).tids.push_back(tid);
        result(tri.b, tri.a).tids.push_back(tid);
        
        result(tri.a, tri.c).tids.push_back(tid);
        result(tri.c, tri.a).tids.push_back(tid);
        
        result(tri.b, tri.c).tids.push_back(tid);
        result(tri.c, tri.b).tids.push_back(tid);
    }
    
    return result;
}


template<class VertData, class TriData>
std::vector<uint> Mesh<VertData,TriData>::getComponentIds()
{
    UnionFind uf(verts.size());
    for(const Tri &tri : tris) {
        uf.unionIds(tri.a, tri.b);
        uf.unionIds(tri.a, tri.c);
    }
    
    return uf.dump();
}








// Surface Point Stuff

template<class VertData, class TriData> inline
Vec3d Mesh<VertData,TriData>::sp_pos(SPptr sp)
{
    Tri &tri = tris[sp->tid];
    Vec3d pos(0,0,0);
    for(uint k=0; k<3; k++)
        pos += sp->bary[k] * verts[tri.v[k]].pos;
    return pos;
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::getSurfPoints(std::vector<Vec3d> *points)
{
    ENSURE(valid_surfps());
    
    points->resize(surfps.size());
    uint write = 0;
    surfps.for_each([&](SPptr sp) {
        // compute position
        Tri         &tri        = tris[sp->tid];
        Vec3d       bary        = sp->bary;
        Vec3d       pos(0,0,0);
        for(uint i=0; i<3; i++)
            pos += bary[i] * verts[tri.v[i]].pos;
        // pack position
        (*points)[write] = pos;
        write++;
    });
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::initPointDistribution()
{
    //randomSurfacePointSeed(tris.size()/2);
    randomSurfacePointSeed(6);
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::randomSurfacePointSeed(uint N)
{
    for(uint k=0; k<N; k++) {
        uint tid = randMod(tris.size());
        if(!tris[tid].sp) {
            SPptr   sp              = tris[tid].sp      = surfps.alloc();
                    sp->tid         = tid;
                    sp->bary        = Vec3d(1,1,1)/3.0;
        }
    }
}

template<class VertData, class TriData>
bool Mesh<VertData,TriData>::valid_surfps()
{
    // check that every surface point
    // references a valid triangle, and
    // that that triangle references back to the surface point
    std::vector<bool> reffed(tris.size(), false);
    uint count = 0;
    surfps.for_each([&](SPptr sp) {
        ENSURE(sp->tid < tris.size());
        Tri &tri = tris[sp->tid];
        ENSURE(sp == tri.sp);
        reffed[sp->tid] = true;
        
        // check barycentric coordinates for sanity
        ENSURE(sp->bary[0] >= 0.0 &&
               sp->bary[1] >= 0.0 &&
               sp->bary[2] >= 0.0);
        double sum = sp->bary[0] + sp->bary[1] + sp->bary[2];
        ENSURE(sum < 1.1);
        
        count++;
    });
    //std::cout << surfps.size() << std::endl;
    ENSURE(count == surfps.size());
    
    // check that the rest of the triangles have null sp
    for(uint i=0; i<tris.size(); i++) {
        if(!reffed[i])
            ENSURE(tris[i].sp == nullptr);
    }
    return true;
}







