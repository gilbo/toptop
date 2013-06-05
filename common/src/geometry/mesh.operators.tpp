// +-------------------------------------------------------------------------
// | mesh.operators.tpp
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

template<class VertData, class TriData>
Eigen::SparseMatrix<double> Mesh<VertData, TriData>::combinatorialLaplacian()
{
    NeighborCache neighbors = createNeighborCache();
    
    std::vector< Eigen::Triplet<double> > triplets;
    
    for(uint i=0; i<verts.size(); i++) {
        double sum = 0.0;
        for(NeighborEntry &entry : neighbors.skeleton[i]) {
            uint j = entry.vid;
            triplets.push_back(Eigen::Triplet<double>(i,j,1.0));
            sum += 1.0;
        }
        triplets.push_back(Eigen::Triplet<double>(i,i,-sum));
    }
    
    Eigen::SparseMatrix<double> laplacian(verts.size(), verts.size());
    laplacian.setFromTriplets(triplets.begin(), triplets.end());
    return laplacian;
}

template<class VertData, class TriData>
Eigen::SparseMatrix<double> Mesh<VertData, TriData>::meanValueLaplacian()
{
    NeighborCache neighbors = createNeighborCache();
    
    std::vector< Eigen::Triplet<double> > triplets;
    
    for(uint i=0; i<verts.size(); i++) {
        uint nNonManifold = 0; // number of non-manifold edges incident
        for(NeighborEntry &entry : neighbors.skeleton[i]) {
            if(entry.tids.size() != 2)  nNonManifold++;
        }
        if(nNonManifold > 0) {
            double sum = nNonManifold;
            for(NeighborEntry &entry : neighbors.skeleton[i]) {
                if(entry.tids.size() == 2)  continue;
                uint j = entry.vid;
                triplets.push_back(Eigen::Triplet<double>(i,j,1.0/sum));
            }
        } else { // manifold point
            double sum = neighbors.skeleton[i].size();
            for(NeighborEntry &entry : neighbors.skeleton[i]) {
                uint j = entry.vid;
                triplets.push_back(Eigen::Triplet<double>(i,j,1.0/sum));
            }
        }
        triplets.push_back(Eigen::Triplet<double>(i,i,-1.0));
    }
    
    Eigen::SparseMatrix<double> laplacian(verts.size(), verts.size());
    laplacian.setFromTriplets(triplets.begin(), triplets.end());
    return laplacian;
}

template<class VertData, class TriData>
Eigen::SparseMatrix<double> Mesh<VertData, TriData>::buildIdentity()
{
    return buildDiagonal([](VertData &) {
        return 1.0;
    });
}

struct LaplacianCacheEntry
{
    uint vid;
    double value;
    uint count;
    LaplacianCacheEntry() {}
    LaplacianCacheEntry(uint vert) : vid(vert), value(0.0), count(1) {}
};

struct LaplacianCache
{
    std::vector< ShortVec<LaplacianCacheEntry, 8> > edges;
    
    inline double& operator()(uint i, uint j) {
        uint N = edges[i].size();
        for(uint k = 0; k < N; k++) {
            if(edges[i][k].vid == j) {
                edges[i][k].count++;
                return edges[i][k].value;
            }
        }
        edges[i].push_back(LaplacianCacheEntry(j));
        return edges[i][N].value;
    }
    
    LaplacianCache(uint N) : edges(N) {}
};


template<class VertData, class TriData>
Eigen::SparseMatrix<double> Mesh<VertData, TriData>::cotanLaplacian()
{
    LaplacianCache cache(verts.size());
    
    // accumulate edge weights
    for(Tri &tri : tris) {
        // compute (2*triangle_area)^2
        Vec3d  p0    = verts[tri.v[0]].pos;
        Vec3d  p1    = verts[tri.v[1]].pos;
        Vec3d  p2    = verts[tri.v[2]].pos;
        Vec3d  prod  = cross(p1-p0, p2-p1);
        double area2 = len2(prod);
        
        double cotan_01 = cotanContribution(&tri, 2, area2);
        double cotan_02 = cotanContribution(&tri, 1, area2);
        double cotan_12 = cotanContribution(&tri, 0, area2);
        
        cache(tri.v[0], tri.v[1]) += cotan_01;
        cache(tri.v[1], tri.v[0]) += cotan_01;
        
        cache(tri.v[0], tri.v[2]) += cotan_02;
        cache(tri.v[2], tri.v[0]) += cotan_02;
        
        cache(tri.v[1], tri.v[2]) += cotan_12;
        cache(tri.v[2], tri.v[1]) += cotan_12;
    }
    
    // convert out to triples for the matrix package
    // and compute the diagonal as we go
    std::vector< Eigen::Triplet<double> > triplets;
    for(uint i=0; i<verts.size(); i++) {
        uint nNonManifold = 0; // number of non-manifold edges incident
        for(LaplacianCacheEntry &entry : cache.edges[i]) {
            if(entry.count != 2)  nNonManifold++;
        }
        
        double sum = 0.0;
        
        for(LaplacianCacheEntry &entry : cache.edges[i]) {
            //if(nNonManifold > 0 && entry.count == 2)  continue;
            uint j = entry.vid;
            sum += entry.value;
            triplets.push_back(Eigen::Triplet<double>(i,j,entry.value));
        }
        
        triplets.push_back(Eigen::Triplet<double>(i,i,-sum));
    }
    
    Eigen::SparseMatrix<double> laplacian(verts.size(), verts.size());
    laplacian.setFromTriplets(triplets.begin(), triplets.end());
    return laplacian;
}


template<class VertData, class TriData> inline
double Mesh<VertData, TriData>::cotanContribution(
    Tri *tri, uint vidx, double area2
) {
    // edge dot product
    // over squared normalizations?
    // times area squared
    Vec3d       p0          = verts[tri->v[(vidx+0)%3]].pos;
    Vec3d       p1          = verts[tri->v[(vidx+1)%3]].pos;
    Vec3d       p2          = verts[tri->v[(vidx+2)%3]].pos;
    Vec3d       edge01      = p1 - p0;
    Vec3d       edge20      = p0 - p2;
    
    double      denom       = len2(edge01) * len2(edge20);
    double      result      = dot(edge01, edge20) * area2;
    if(denom > 0) result   /= denom;
    else        result      = 0.0;
    
    return -result;
}




template<class VertData, class TriData>
Eigen::SparseMatrix<double> Mesh<VertData, TriData>::areaWeightLaplacian()
{
    LaplacianCache cache(verts.size());
    
    // accumulate edge weights
    for(Tri &tri : tris) {
        // compute barycenter and midpoints
        Vec3d  p0    = verts[tri.v[0]].pos;
        Vec3d  p1    = verts[tri.v[1]].pos;
        Vec3d  p2    = verts[tri.v[2]].pos;
        Vec3d  bary  = (p0 + p1 + p2) / 3.0;
        Vec3d  mid01 = (p0 + p1) / 2.0;
        Vec3d  mid02 = (p0 + p2) / 2.0;
        Vec3d  mid12 = (p1 + p2) / 2.0;
        double len01 = len(bary - mid01);
        double len02 = len(bary - mid02);
        double len12 = len(bary - mid12);
        
        cache(tri.v[0], tri.v[1]) += len01;
        cache(tri.v[1], tri.v[0]) += len01;
        
        cache(tri.v[0], tri.v[2]) += len02;
        cache(tri.v[2], tri.v[0]) += len02;
        
        cache(tri.v[1], tri.v[2]) += len12;
        cache(tri.v[2], tri.v[1]) += len12;
    }
    
    // convert out to triples for the matrix package
    // and compute the diagonal as we go
    std::vector< Eigen::Triplet<double> > triplets;
    for(uint i=0; i<verts.size(); i++) {
        double sum = 0.0;
        for(LaplacianCacheEntry &entry : cache.edges[i]) {
            uint j = entry.vid;
            sum += entry.value;
            triplets.push_back(Eigen::Triplet<double>(i,j,entry.value));
        }
        triplets.push_back(Eigen::Triplet<double>(i,i,-sum));
    }
    
    Eigen::SparseMatrix<double> laplacian(verts.size(), verts.size());
    laplacian.setFromTriplets(triplets.begin(), triplets.end());
    return laplacian;
}





template<class VertData, class TriData>
Eigen::SparseMatrix<double> Mesh<VertData, TriData>::invDistLaplacian()
{
    NeighborCache neighbors = createNeighborCache();
    
    std::vector< Eigen::Triplet<double> > triplets;
    
    for(uint i=0; i<verts.size(); i++) {
        double sum = 0.0;
        for(NeighborEntry &entry : neighbors.skeleton[i]) {
            //if(nNonManifold > 0 && entry.tids.size() == 2)  continue;
            
            uint            j           = entry.vid;
            Vec3d           pi          = verts[i].pos;
            Vec3d           pj          = verts[j].pos;
            Vec3d           eij         = pj - pi;
            double           elen        = len(eij);
            double           invlen      = 1.0 / elen;
            if(elen == 0)   invlen      = 0.0;
                            sum        += invlen;
            triplets.push_back(Eigen::Triplet<double>(i,j,invlen));
        }
        triplets.push_back(Eigen::Triplet<double>(i,i,-sum));
    }
    
    Eigen::SparseMatrix<double> laplacian(verts.size(), verts.size());
    laplacian.setFromTriplets(triplets.begin(), triplets.end());
    return laplacian;
}







