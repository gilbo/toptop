// +-------------------------------------------------------------------------
// | mesh.distances.tpp
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

#include <queue>
using std::priority_queue;
using std::pair;

template<class VertData, class TriData>
void Mesh<VertData,TriData>::minPathDistances(
    const Mesh<VertData,TriData>::Isct &isct
) {
    NeighborCache cache = createNeighborCache();
    
    // do Dijkstra's algorithm
    priority_queue< pair<double, uint> > work; // use negative distance
    std::vector<bool> visited(verts.size(), false);
    
    // Seed with the three vertices of the intersected triangle
    Tri *tri = &(tris[isct.tri_id]);
    Vec3d base_pos = isct.isct;
    for(uint k=0; k<3; k++) {
        uint idx = tri->v[k];
        double distance = len(verts[idx].pos - base_pos);
        verts[idx].distField = distance; // current best guess
        work.push(pair<double,int>(-distance, idx)); // stuff...
    }
    
    while(!work.empty()) {
        double  base_dist   = -work.top().first;
        uint    idx         = work.top().second;
        work.pop();
        
        // ignore if we've already processed this vertex...
        if(visited[idx])    continue;
        
        // otherwise, commit this value
        verts[idx].distField = base_dist;
        visited[idx] = true;
        
        // and update the guesses for neighbors' distances
        for(NeighborEntry &entry : cache.skeleton[idx]) {
            uint nidx = entry.vid;
            double dist = len(verts[nidx].pos - verts[idx].pos);
            if(!visited[nidx])
                work.push(pair<double,int>(-(base_dist + dist), nidx));
        }
    }
    
    // mark all disconnected nodes as infinitely far away
    for(uint i=0; i<verts.size(); i++) {
        if(!visited[i]) {
            verts[i].distField = DBL_MAX;
        }
    }
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::existsPathDistances(
    const Mesh<VertData,TriData>::Isct &isct
) {
    NeighborCache cache = createNeighborCache();
    
    // do a flood-fill (could speed up with union-join
    priority_queue<uint> work;
    std::vector<bool> visited(verts.size(), false);
    
    // Seed with the three vertices of the intersected triangle
    Tri &tri = tris[isct.tri_id];
    for(uint k=0; k<3; k++) {
        uint idx = tri.v[k];
        visited[idx] = true;
        work.push(idx);
    }
    
    while(!work.empty()) {
        uint    idx         = work.top();
        work.pop();
        
        // try to add all neighbors to queue
        for(NeighborEntry &entry : cache.skeleton[idx]) {
            uint nidx = entry.vid;
            if(!visited[nidx]) {
                visited[nidx] = true;
                work.push(nidx);
            }
        }
    }
    
    // convert visited values into the distance field
    for(uint i=0; i<verts.size(); i++) {
        verts[i].distField = (visited[i])? 0.0f : 1.0f;
    }
}


// BROKEN FOR SOME MESHES WITH ZERO/NEGATIVE COTAN LAPLACIAN WEIGHTS
template<class VertData, class TriData>
void Mesh<VertData,TriData>::heatDistanceToSurfacePoint()
{
    // timestep is computed from average edge length
    double sum_len2 = 0.0;
    uint count = 0;
    for(Tri &tri : tris) {
        for(uint k=0; k<3; k++) {
            Vec3d e = verts[tri.v[(k+1)%3]].pos - verts[tri.v[k]].pos;
            sum_len2 += len2(e);
            count++;
        }
    }
    double avg_len2 = sum_len2 / double(count);
    
    double timestep = avg_len2;
    std::cout << "timestep: " << timestep << std::endl;
    Eigen::SparseMatrix<double> laplacian    = cotanLaplacian();
    //Eigen::SparseMatrix<double> laplacian    = combinatorialLaplacian();
    Eigen::SparseMatrix<double> id           = buildIdentity();
    Eigen::SparseMatrix<double> step_system  = id - timestep * laplacian;
    
    Eigen::VectorXd indicator(verts.size());
    for(uint i=0; i<verts.size(); i++)  indicator(i) = 0.0f;
    for(Tri &tri : tris) {
        if(tri.sp) {
            for(uint k=0; k<3; k++)
                indicator(tri.v[k]) =
                    std::max(indicator(tri.v[k]),
                             tri.sp->bary[k]);
        }
    }
    /*NeighborCache cache = createNeighborCache();
    for(uint i=0; i<verts.size(); i++) {
        if(indicator(i) == 0.0f)    continue;
        for(uint crawl=0; crawl<cache.skeleton[i].size(); crawl++) {
            uint j = cache.skeleton[i][crawl].vid;
            step_system.coeffRef(i,j) = 0.0f;
        }
        step_system.coeffRef(i,i) = 1.0f;
    }*/
        
    // need to solve the laplacian equation with the indicator
    // to do a forward implicit step...
    Eigen::SimplicialCholesky< Eigen::SparseMatrix<double> >
        implicit_heat_step(step_system);
    Eigen::VectorXd heat = implicit_heat_step.solve(indicator);
    
    // construct a normalized vector field out of the heat
    // and compute its divergence
    std::vector<double> divergence(verts.size(), 0.0f);
    for(uint i=0; i<tris.size(); i++) {
        // compute normalized gradient
        Vec3d   p0      = verts[tris[i].v[0]].pos;
        Vec3d   p1      = verts[tris[i].v[1]].pos;
        Vec3d   p2      = verts[tris[i].v[2]].pos;
        Vec3d   e01     = p1 - p0;
        Vec3d   e12     = p2 - p1;
        Vec3d   e20     = p0 - p2;
        Vec3d   normal  = cross(e01, e12);
        double  area2   = len2(normal);
        
        Vec3d   grad    = heat(tris[i].v[0]) * e12
                        + heat(tris[i].v[1]) * e20
                        + heat(tris[i].v[2]) * e01;
                grad    = cross(normal, grad);
        Vec3d   dir     = -normalized(grad);
        
        // now compute divergence of the field
        double  cotan0  = cotanContribution(&(tris[i]), 0, area2);
        double  cotan1  = cotanContribution(&(tris[i]), 1, area2);
        double  cotan2  = cotanContribution(&(tris[i]), 2, area2);
        
        double  div0    = cotan2 * dot(e01, dir) - cotan1 * dot(e20, dir);
        double  div1    = cotan0 * dot(e12, dir) - cotan2 * dot(e01, dir);
        double  div2    = cotan1 * dot(e20, dir) - cotan0 * dot(e12, dir);
        divergence[tris[i].v[0]] += div0;
        divergence[tris[i].v[1]] += div1;
        divergence[tris[i].v[2]] += div2;
    }
    Eigen::VectorXd div(verts.size());
    for(uint i=0; i<verts.size(); i++)
        div(i) = divergence[i];
    
    // now we recover the distances by solving a Poisson equation
    // where the vector is given by the divergences
    Eigen::SimplicialCholesky< Eigen::SparseMatrix<double> >
        poisson_problem(laplacian);
    Eigen::VectorXd dist = poisson_problem.solve(div);
    // do a component by component offset
    std::vector<uint> component_ids = getComponentIds();
    // first find the minimum in each component
    std::vector<double> minvals(verts.size(), DBL_MAX); // over-allocate space
    for(uint i=0; i<verts.size(); i++) {
        uint c_id = component_ids[i];
        minvals[c_id] = std::min(minvals[c_id], dist(i));
    }
    // then, offset to make this minimum equal to zero
    for(uint i=0; i<verts.size(); i++) {
        double minimum = minvals[component_ids[i]];
        dist(i) -= minimum;
    }
    
    // get the distances out
    processScalarField(&dist, [](VertData &v, double d) {
        //if(d == 0.0)    v.distField = 0.0;
        //else            v.distField = -std::log(d);
        v.distField = d;
    });
}


template<class VertData, class TriData>
void Mesh<VertData,TriData>::minPathDistancesToSurfacePoints()
{
    NeighborCache cache = createNeighborCache();
    
    // do Dijkstra's algorithm
    priority_queue< pair<double, uint> > work; // use negative distance
    
    for(VertData &v : verts)
        v.distField = FLT_MAX;
    
    // run through and initialize the distance field for every triangle
    // that has a surface point
    surfps.for_each([&](SPptr sp) {
        Tri &tri = tris[sp->tid];
        // compute surface point coordinates
        Vec3d barycenter = sp_pos(sp);
        for(uint k=0; k<3; k++) {
            VertData    &vert               = verts[tri.v[k]];
            double      dist                = len(vert.pos - barycenter);
            if(dist < vert.distField) {
                        vert.distField      = dist;
                        work.push(std::make_pair(-dist, tri.v[k]));
            }
        }
    });
    
    // propagate outwards
    while(!work.empty()) {
        double  stored_dist     = -work.top().first;
        uint    idx             = work.top().second;
        double  base_dist       = verts[idx].distField;
        work.pop();
        
        // if we already have a better answer, then we
        // don't need to do this work...
        if(base_dist < stored_dist)     continue;
        
        // now, see if we can come up with closer values for
        // our neighbors than they already have...
        for(NeighborEntry &entry : cache.skeleton[idx]) {
            uint    nidx        = entry.vid;
            double  hop_dist    = len(verts[nidx].pos - verts[idx].pos);
            double  dist        = base_dist + hop_dist;
            if(dist < verts[nidx].distField) {
                    verts[nidx].distField = dist;
                    work.push(std::make_pair(-dist, nidx));
            }
        }
    }
}


// TODO: Should remove redundancy between this and the initial
// creation of a distance field...
template<class VertData, class TriData>
void Mesh<VertData,TriData>::updatePathDistance(
    const NeighborCache &cache, SPptr sp
) {
    priority_queue< pair<double, uint> > work;
    
    // update the immediate vertices...
    Vec3d barycenter = sp_pos(sp);
    for(uint k=0; k<3; k++) {
        VertData    &vert               = verts[tris[sp->tid].v[k]];
        double      dist                = len(vert.pos - barycenter);
        if(dist < vert.distField) {
                    vert.distField      = dist;
                    work.push(std::make_pair(-dist, tris[sp->tid].v[k]));
        }
    }
    
    // propagate outwards
    while(!work.empty()) {
        double  stored_dist     = -work.top().first;
        uint    idx             = work.top().second;
        double  base_dist       = verts[idx].distField;
        work.pop();
        
        // if we already have a better answer, then we
        // don't need to do this work...
        if(base_dist < stored_dist)     continue;
        
        // now, see if we can come up with closer values for
        // our neighbors than they already have...
        for(const NeighborEntry &entry : cache.skeleton[idx]) {
            uint    nidx        = entry.vid;
            double  hop_dist    = len(verts[nidx].pos - verts[idx].pos);
            double  dist        = base_dist + hop_dist;
            if(dist < verts[nidx].distField) {
                    verts[nidx].distField = dist;
                    work.push(std::make_pair(-dist, nidx));
            }
        }
    }
}






