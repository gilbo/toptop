// +-------------------------------------------------------------------------
// | mesh.decl.h
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

#include "rawMesh.h"

#include "prelude.h"
#include <vector>
#include <set>

#include "vec.h"
#include "ray.h"
#include "shortVec.h"

#include "iterPool.h"

// for use of Linear Operator module
#include <Eigen/Sparse>



// in order to use distance functions,
// VertData and TriData must support
struct DistanceVertexData {
    double distField;
};
struct DistanceTriangleData {
};

struct CollisionPacket {
    double  time;
    uint    tid;
};

enum DeadState {
    ALIVE = 0,
    DEAD  = 1,
    UNKNOWN = 2,
    ERROR = 3,
};
// lattice join
inline DeadState joinDead(DeadState a, DeadState b) {
    //ENSURE(a == b || a == UNKNOWN || b == UNKNOWN);
    return (a == UNKNOWN)?  b :
                ((b == UNKNOWN)? a :
                      ((a == b)? a : ERROR ));
}

// also designed to support TopoChange module
struct CollisionVertexData {
    Vec3d   velocity;
    ShortVec<CollisionPacket, 2> collisions;
    uint        pastCollisionParity;
    DeadState   dead;
    double      smoothDead; // debug data for smoothing...
};
struct CollisionTriangleData {
    DeadState   dead;
};


// only for internal use, please do not use as client

    struct TopoVert;
    struct TopoEdge;
    struct TopoTri;
    
    using Vptr = TopoVert*;
    using Eptr = TopoEdge*;
    using Tptr = TopoTri*;

template<class VertData, class TriData>
struct IsctVertEdgeTriInput
{
    VertData*   e[2];
    VertData*   t[3];
};

template<class VertData, class TriData>
struct IsctVertTriTriTriInput
{
    VertData*   t[3][3];
};

template<class VertData, class TriData>
struct SubdivideTriInput
{
    TriData*    pt;
    VertData*   pv[3];
    VertData*   v[3];
};

// in order to perform intersections, VertData and TriData must support
struct IsctVertexData {
    // specify how to compute new data for vertices formed by intersections
    /*
    // vertices on edge and triangle forming an intersection...
    void isct(IsctVertEdgeTriInput input);
    void isct(IsctVertTriTriTriInput input);
    void isctInterpolate(const VertData &v0, const VertData &v1);
    */
};
struct IsctTriangleData {
    // specify how to compute new data for a triangle in the event
    // that it is merged with another triangle (merge)
    // split into two triangles (split)
    // or that the triangle is moved (move)
    /*
    void subdivide(SubdivideTriInput input);
    */
};



// in order to perform remeshing, VertData and TriData must support
struct RemeshVertexData {
    bool manifold; // whether this point is manifold.
                   // useful for modifying interpolation behavior
    // specify how to compute new data for a vertex in the event of
    // either an edge collapse (via merge) or edge split (via interpolate)
    /*
    void merge(const VertData &v0, const VertData &v1);
    void interpolate(const VertData &v0, const VertData &v1);
    */
};
struct RemeshTriangleData {
    // specify how to compute new data for a triangle in the event
    // that it is merged with another triangle (merge)
    // split into two triangles (split)
    // or that the triangle is moved (move)
    /*
    void merge(const TriData &t0, const TriData &t1);
    static void split(TriData &t0, TriData &t1, const TriData &t_orig);
    void move(const TriData &t_old);
    */
};

struct RemeshOptions
{
    double maxEdgeLength;
    double minEdgeLength;
    double minAngle;
    double maxAngle;
    RemeshOptions() :
        maxEdgeLength(1.0),
        minEdgeLength(0.3),
        minAngle(5.0),
        maxAngle(170.0)
    {}
};


template<class VertData, class TriData>
class Mesh
{
public:
    Mesh();
    Mesh(Mesh &&src);
    Mesh(const RawMesh<VertData,TriData> &raw);
    virtual ~Mesh();
    
    void operator=(Mesh &&src);
    
    // validity check:
    //  - all numbers are well-defined and finite
    //  - all triangle vertex indices are in the right range
    bool valid() const;
    
    RawMesh<VertData,TriData> raw() const;
    
    inline int numVerts() const { return verts.size(); }
    inline int numTris() const { return tris.size(); }
    
    inline void for_verts(std::function<void(VertData &)> func);
    inline void for_tris(std::function<void(TriData &,
                          VertData &, VertData &, VertData &)> func);
    inline void for_edges(
        std::function<void(VertData &, VertData &)> start,
        std::function<void(TriData &t,
                           VertData &, VertData &, VertData &)> each_tri
    );
    
    // form the disjoint union of two meshes
    void disjointUnion(const Mesh &cp);
    
    struct Isct {
        Ray3d   ray;
        bool    exists;
        
        uint    tri_id;
        Vec3d   isct;
        Vec3d   bary;
    };
    Isct pick(Ray3d ray);
    inline void accessIsct(const Isct &isct,
        std::function<void(TriData &,
            VertData &, VertData &, VertData &)> func);
    
    
public: // REMESHING module
    // REQUIRES:
    //  - MinimalData
    //  - RemeshData
    void remesh();
    RemeshOptions remesh_options;
    
public: // DISTANCES module
    // REQUIRES:
    //  - MinimalData
    //  - DistanceData
    // computes the distance field emanating from isct,
    // and stores values on the elements
    void minPathDistances(const Isct &isct);
        // existsDistances: 0 if same component, 1 otherwise
    void existsPathDistances(const Isct &isct);
    // computes distance to the surface point set
    // using the heat method (see paper "Geodesics in Heat" on the ArXiV)
    void heatDistanceToSurfacePoint();
    // just uses a simple minimum path on the edge skeleton
    void minPathDistancesToSurfacePoints();
    
public: // OPERATORS module
    // all operators on vertex fields:
        // the graph laplacian treats non-manifold/boundary points specially...
    Eigen::SparseMatrix<double> meanValueLaplacian();
        // use meanValue for explicit-solve smoothing only!
    Eigen::SparseMatrix<double> invDistLaplacian();
    Eigen::SparseMatrix<double> cotanLaplacian();
    Eigen::SparseMatrix<double> combinatorialLaplacian();
    Eigen::SparseMatrix<double> areaWeightLaplacian();
    Eigen::SparseMatrix<double> buildIdentity();
    inline Eigen::SparseMatrix<double> buildDiagonal(
        std::function<double(VertData &)> func // determine entry
    );
    // work with vector valued fields
    inline Eigen::MatrixX3d    buildVecField(
        std::function<Vec3d(VertData &)> func // determine entry
    );
    inline void                processVecField(
        const Eigen::MatrixX3d *field,
        std::function<void(VertData &, const Vec3d&)> func
    );
    // scalar valued fields
    inline Eigen::VectorXd     buildScalarField(
        std::function<double(VertData &)> func // determine entry
    );
    inline void                processScalarField(
        const Eigen::VectorXd *field,
        std::function<void(VertData &, double)> func
    );

public: // ISCT (intersections) module
    void resolveIntersections(); // makes all intersections explicit
    // TESTING
    void testingComputeStaticIsctPoints(std::vector<Vec3d> *points);
    void testingComputeStaticIsct(std::vector<Vec3d> *points,
               std::vector< std::pair<Vec3d,Vec3d> > *edges);
    
public: // COLLISION module
    void findCollisions();
    
public: // TOPO CHANGE module
   // run collision detection before using any of these routines
    void previewComponentMajorityVote();
    void previewLinearSolveVote();
    //void previewMinCutVote();
    // given an updated parity, we can 
    void applyDeathField();
    // run to do entire topo change from parity update through
    // cutting and gluing the final mesh
    void componentMajorityVote();
    void linearSolveVote();
    //void minCutVote();
    
    // patch: should probably move elsewhere.
    void conformOrientations();
    
private:    // Internal Formats
    struct SurfacePoint {
        uint tid; // triangle this point is attached to
        Vec3d bary; // barycentric coordinates for the point
    };
    using SPptr = SurfacePoint*;
    
    void getSurfPoints(std::vector<Vec3d> *points);
    void initPointDistribution();
    
    inline Vec3d sp_pos(SPptr sp);
    
    struct Tri {
        TriData data;
        union {
            struct {
                uint a, b, c; // vertex ids
            };
            uint v[3];
        };
        SPptr   sp;
        
        inline Tri() : sp(nullptr) {}
    };
    
    inline void merge_tris(uint tid_result, uint tid0, uint tid1);
           void split_helper(uint t0ref, uint t1ref, uint t_orig_ref);
    inline void split_tris(uint t0ref, uint t1ref, uint t_orig_ref);
    inline void move_tri(Tri &t_new, Tri &t_old);
           void subdivide_helper(uint t_piece_ref, uint t_parent_ref);
    inline void subdivide_tri(uint t_piece_ref, uint t_parent_ref);
    
private:    // DATA
    std::vector<Tri>        tris;
    std::vector<VertData>   verts;
    IterPool<SurfacePoint>  surfps;
    
private:    // caches
    struct NeighborEntry {
        uint vid;
        ShortVec<uint, 2> tids;
        inline NeighborEntry() {}
        inline NeighborEntry(uint vid_) : vid(vid_) {}
    };
    struct NeighborCache {
        std::vector< ShortVec<NeighborEntry, 8> > skeleton;
        inline NeighborEntry& operator()(uint i, uint j) {
            uint N = skeleton[i].size();
            for(uint k = 0; k < N; k++) {
                if(skeleton[i][k].vid == j)
                    return skeleton[i][k];
            }
            skeleton[i].push_back(NeighborEntry(j));
            return skeleton[i][N];
        }
    };
    NeighborCache createNeighborCache();
    
    // parallel to vertex array
    std::vector<uint> getComponentIds();
    
private:    // Surface Point Distribution Support
    void updatePathDistance(const NeighborCache &cache, SPptr sp);
    void randomSurfacePointSeed(uint n); // attempts to add n points;
                                         // may add fewer in case of collisions
    bool valid_surfps();
private:    // TopoChange Support
    class TopoChangeProblem;
private:    // Collision Support
    class CollisionProblem;
private:    // Operators Support
    // vidx should be 0,1,2
    // area2 = (2*triangle_area)^2
    inline double cotanContribution(Tri *tri, uint vidx, double area2);
private:    // TopoCache Support
    struct TopoCache;
private:    // Isct Support
    class  IsctProblem; // implements intersection functionality
        class TriangleProblem; // support type for IsctProblem
        using Tprob = TriangleProblem*;
    
private:    // Remeshing Support
    struct RemeshScratchpad;
    
    Eptr allocateRemeshEdge(RemeshScratchpad &);
    void deallocateRemeshEdge(RemeshScratchpad &, Eptr);
    
    void edgeSplit(RemeshScratchpad &,
                   Eptr e_split);
    void edgeCollapse(RemeshScratchpad &,
                      Eptr e_collapse,
                      bool collapsing_tetrahedra_disappear);
    
    // Need edge scoring routines...
    void scoreAndEnqueue(
        std::set< std::pair<double, Eptr> > &queue, Eptr edge);
    void dequeue(
        std::set< std::pair<double, Eptr> > &queue, Eptr edge);
    double computeEdgeScore(Eptr edge);
    
    // support functions
    void populateTriFromTopoTri(Tptr t);
    // calls the first function once, then the second once for each triangle
    inline void edgeNeighborhood(
        Eptr edge,
        std::function<void(VertData &v0, VertData &v1)> once,
        std::function<void(VertData &v0, VertData &v1,
                           VertData &vopp, TriData &t)> each_tri
    );
};


// inline functions

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::for_verts(
    std::function<void(VertData &v)> func
) {
    for(auto &v : verts)
        func(v);
}

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::for_tris(
    std::function<void(TriData &, VertData &, VertData &, VertData &)> func
) {
    for(auto &tri : tris) {
        auto &a = verts[tri.a];
        auto &b = verts[tri.b];
        auto &c = verts[tri.c];
        func(tri.data, a, b, c);
    }
}

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::for_edges(
    std::function<void(VertData &, VertData &)> start,
    std::function<void(TriData &t,
                       VertData &, VertData &, VertData &)> each_tri
) {
    NeighborCache cache = createNeighborCache();
    for(uint i=0; i<cache.skeleton.size(); i++) {
        for(auto &entry : cache.skeleton[i]) {
            uint j = entry.vid;
            start(verts[i], verts[j]);
            for(uint tid : entry.tids) {
                Tri &tri = tris[tid];
                each_tri(tri.data, verts[tri.a], verts[tri.b], verts[tri.c]);
            }
        }
    }
}

template<class VertData, class TriData>
inline void Mesh<VertData,TriData>::accessIsct(
    const Isct &isct,
    std::function<void(TriData &, VertData &, VertData &, VertData &)> func
) {
    Tri &tri = tris[isct.tri_id];
    auto &a = verts[tri.a];
    auto &b = verts[tri.b];
    auto &c = verts[tri.c];
    func(tri.data, a, b, c);
}

template<class VertData, class TriData>
inline Eigen::SparseMatrix<double> Mesh<VertData,TriData>::buildDiagonal(
    std::function<double(VertData &)> func
) {
    uint N = verts.size();
    
    std::vector< Eigen::Triplet<double> > triplets(N);
    for(uint i=0; i<N; i++)
        triplets[i] = Eigen::Triplet<double>(i,i,
            func(verts[i])
        );
    
    Eigen::SparseMatrix<double> identity(N,N);
    identity.setFromTriplets(triplets.begin(), triplets.end());
    return identity;
}

template<class VertData, class TriData>
inline Eigen::MatrixX3d    Mesh<VertData,TriData>::buildVecField(
    std::function<Vec3d(VertData &)> func // determine entry
) {
    Eigen::MatrixX3d field(verts.size(), 3);
    for(uint i=0; i<verts.size(); i++) {
        Vec3d vec = func(verts[i]);
        for(uint j=0; j<3; j++)
            field(i,j) = vec.v[j];
    }
    return field;
}

template<class VertData, class TriData>
inline void                Mesh<VertData,TriData>::processVecField(
    const Eigen::MatrixX3d *field,
    std::function<void(VertData &, const Vec3d&)> func
) {
    for(uint i=0; i<verts.size(); i++) {
        Vec3d vec((*field)(i,0), (*field)(i,1), (*field)(i,2));
        func(verts[i], vec);
    }
}

template<class VertData, class TriData>
inline Eigen::VectorXd    Mesh<VertData,TriData>::buildScalarField(
    std::function<double(VertData &)> func // determine entry
) {
    Eigen::VectorXd field(verts.size());
    for(uint i=0; i<verts.size(); i++) {
        field(i) = func(verts[i]);
    }
    return field;
}

template<class VertData, class TriData>
inline void                Mesh<VertData,TriData>::processScalarField(
    const Eigen::VectorXd *field,
    std::function<void(VertData &, double)> func
) {
    for(uint i=0; i<verts.size(); i++) {
        func(verts[i], (*field)(i));
    }
}





