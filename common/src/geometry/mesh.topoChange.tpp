// +-------------------------------------------------------------------------
// | mesh.topoChange.tpp
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

#include "mesh.isct.tpp"
#include "mesh.operators.tpp"

#include "unionFind.h"

// need to grab Vladimir's maxflow library to do mincut
// his code is omitted here to avoid incorporating the license
//#include "graph.h"

struct TopoIsctGraphEdge
{
    uint                vid;
    ShortVec<Vec3d, 2>  iscts;
    double              weight;
    double              flow;
    GluePt              glue_cut_pt;
    void accumulateIsct(const Vec3d &v) {
        for(uint k=0; k<iscts.size(); k++)
            if(iscts[k] == v)
                return;
        iscts.push_back(v);
    }
};
using TIGedge = TopoIsctGraphEdge*;

struct RedirectEdge {
    uint    vid;
    TIGedge tig;
};
using ReEdge = RedirectEdge*;

template<class VertData, class TriData>
class Mesh<VertData,TriData>::TopoChangeProblem : public IsctProblem
{
public:
    TopoChangeProblem(Mesh *owner) : IsctProblem(owner)
    {}
    // can call findIntersections() to compute that data...
    virtual ~TopoChangeProblem() {}
    
    void buildEdgeGraph();
    void annotateEdgeGraph();
    
    std::vector<uint> getIsctComponentIds();
    // components formed of triangles connected
    // by non-intersection edges
    std::vector<uint> getTriComponentIds();
    
public: // METHODS for SETTING death field
    // PRECONDITIONS: run collision testing
    // ignore current value of death field
    void populateDeathFieldFromParity();
    // use current value of death field
    void updateDeathFieldFromParity();
    
public: // METHODS for SMOOTHING death field
    // PRECONDITIONS: find intersections and build/annotate edge graph
    void smoothByComponentMajority();
    void smoothByLinearSolve();
    //void smoothByMinCut();
    
public: // METHODS for CUTTING and GLUING using the death field
    // PRECONDITIONS: find intersections and build/annotate edge graph
    //      must have been run
    void cullUnneededIscts();
    // PRECONDITIONS: find intersections must have been run
    void cutAndGlue();

private:
    // create edge (can force error if non-existent using third arg)
    TIGedge getEdge(uint i, uint j, bool allowCreate = true)
    {
        //std::cout << "req " << i << ", " << j << std::endl;
        uint N = egraph[i].size();
        for(uint k=0; k<N; k++)
            if(egraph[i][k]->vid == j)
                return egraph[i][k];
        // could not find this edge, so create it
        ENSURE(allowCreate);
        egraph[i].push_back(egraph_pool.alloc());
        egraph[i][N]->vid = j;
        egraph[i][N]->glue_cut_pt = nullptr;
        
        ReEdge re_ij = regraph_pool.alloc();
        ReEdge re_ji = regraph_pool.alloc();
        re_ij->tig = re_ji->tig = egraph[i][N];
        re_ij->vid = j;
        re_ji->vid = i;
        regraph[i].push_back(re_ij);
        regraph[j].push_back(re_ji);
        
        return egraph[i][N];
    }
    
    GluePt getGluePt(TIGedge tig, Eptr e) {
        if(tig->glue_cut_pt)
            return tig->glue_cut_pt;
        GluePt  glue                = IsctProblem::newGluePt();
                glue->split_type    = true;
                glue->e             = e;
        tig->glue_cut_pt = glue;
        return glue;
    }
    
    Vec3d getSplitCoords(Eptr e) {
        Vec3d p0 = IsctProblem::vPos(e->verts[0]);
        Vec3d p1 = IsctProblem::vPos(e->verts[1]);
        return (p0+p1) / 2.0;
    }
    
    inline
    void for_tig_edges(
        std::function<void(uint,uint,TIGedge)> action
    ) {
        for(uint i=0; i<egraph.size(); i++) {
            for(TIGedge tig : egraph[i]) {
                uint j = tig->vid;
                action(i,j,tig);
            }
        }
    }
    
    inline
    void for_neighbors(
        uint vert,
        std::function<void(uint,TIGedge)> action
    ) {
        for(ReEdge re : regraph[vert])
            action(re->vid, re->tig);
    }
    
    void recordEdgeWeights(const Eigen::SparseMatrix<double> &matrix)
    {
        for_tig_edges([&](uint i, uint j, TIGedge tig) {
            tig->weight = -matrix.coeff(i,j);
            //std::cout << -matrix.coeff(i,j) << std::endl;
        });
    }
    
private:
    // only store edges as
    //      egraph[i][k].vid == j
    // with i < j
    std::vector< ShortVec<TIGedge, 8> >     egraph;
    IterPool<TopoIsctGraphEdge>             egraph_pool;
    std::vector< double >                   egraph_mat_diag;
    std::vector< ShortVec<ReEdge, 8> >      regraph;
    IterPool<RedirectEdge>                  regraph_pool;
};



template<class VertData, class TriData>
void Mesh<VertData,TriData>::componentMajorityVote()
{
    TopoChangeProblem tcprob(this);
    
    tcprob.populateDeathFieldFromParity();
    
    tcprob.findIntersections();
    tcprob.buildEdgeGraph();
    tcprob.annotateEdgeGraph(); // with intersections
    
    tcprob.smoothByComponentMajority();
    
    tcprob.cullUnneededIscts();
    
    tcprob.cutAndGlue();
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::linearSolveVote()
{
    TopoChangeProblem tcprob(this);
    
    tcprob.populateDeathFieldFromParity();
    
    tcprob.findIntersections();
    tcprob.buildEdgeGraph();
    tcprob.annotateEdgeGraph(); // with intersections
    
    tcprob.smoothByLinearSolve();
    
    tcprob.cullUnneededIscts();
    
    tcprob.cutAndGlue();
}

//template<class VertData, class TriData>
//void Mesh<VertData,TriData>::minCutVote()
//{
//    TopoChangeProblem tcprob(this);
//    
//    tcprob.populateDeathFieldFromParity();
//    
//    tcprob.findIntersections();
//    tcprob.buildEdgeGraph();
//    tcprob.annotateEdgeGraph(); // with intersections
//    
//    tcprob.smoothByMinCut();
//    
//    tcprob.cullUnneededIscts();
//    
//    tcprob.cutAndGlue();
//}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::previewComponentMajorityVote()
{
    TopoChangeProblem tcprob(this);
    
    tcprob.updateDeathFieldFromParity();
    
    tcprob.findIntersections();
    tcprob.buildEdgeGraph();
    tcprob.annotateEdgeGraph(); // with intersections
    
    tcprob.smoothByComponentMajority();
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::previewLinearSolveVote()
{
    TopoChangeProblem tcprob(this);
    
    tcprob.updateDeathFieldFromParity();
    
    tcprob.findIntersections();
    tcprob.buildEdgeGraph();
    tcprob.annotateEdgeGraph(); // with intersections
    
    tcprob.smoothByLinearSolve();
}

//template<class VertData, class TriData>
//void Mesh<VertData,TriData>::previewMinCutVote()
//{
//    TopoChangeProblem tcprob(this);
//    
//    tcprob.updateDeathFieldFromParity();
//    
//    tcprob.findIntersections();
//    tcprob.buildEdgeGraph();
//    tcprob.annotateEdgeGraph(); // with intersections
//    
//    tcprob.smoothByMinCut();
//}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::applyDeathField()
{
    TopoChangeProblem tcprob(this);
    
    tcprob.findIntersections();
    tcprob.buildEdgeGraph();
    tcprob.annotateEdgeGraph(); // with intersections
    
    tcprob.cullUnneededIscts();
    
    tcprob.cutAndGlue();
}



template<class VertData, class TriData>
void Mesh<VertData,TriData>::TopoChangeProblem::populateDeathFieldFromParity()
{
    TopoCache::mesh->for_verts([](VertData &v) {
        v.pastCollisionParity = v.collisions.size()%2;
        v.dead = (v.pastCollisionParity == 0)? ALIVE : DEAD;
    });
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::TopoChangeProblem::updateDeathFieldFromParity()
{
    TopoCache::mesh->for_verts([](VertData &v) {
        v.pastCollisionParity =
            (v.pastCollisionParity + v.collisions.size()) % 2;
        v.dead = (v.pastCollisionParity == 0)? ALIVE : DEAD;
    });
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::TopoChangeProblem::smoothByComponentMajority()
{
    std::vector<uint> cids = getIsctComponentIds();
    std::vector<int> score(TopoCache::mesh->verts.size(), 0);
    // gather
    for(uint i=0; i<TopoCache::mesh->verts.size(); i++) {
        uint cid = cids[i];
        VertData &v = TopoCache::mesh->verts[i];
        // cast ballot
        if(v.dead == DEAD)  score[cid]--;
        else                score[cid]++;
    }
    // scatter
    for(uint i=0; i<TopoCache::mesh->verts.size(); i++) {
        uint cid = cids[i];
        VertData &v = TopoCache::mesh->verts[i];
        if(score[cid] < 0)      v.dead = DEAD;
        else                    v.dead = ALIVE;
    }
}

#include <set>
uint countComponents(std::vector<uint> &ids)
{
    std::set<uint> uniques;
    for(uint id : ids)
        uniques.insert(id);
    return uniques.size();
}

std::map< uint, std::vector<uint> >
    groupComponents(const std::vector<uint> &ids)
{
    std::map< uint, std::vector<uint> > components;
    
    for(uint i=0; i<ids.size(); i++)
        components[ids[i]].push_back(i);
    
    return components;
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::TopoChangeProblem::smoothByLinearSolve()
{
    // set up a laplacian-like system
    Eigen::SparseMatrix<double> laplacian =
    //    TopoCache::mesh->combinatorialLaplacian();
    //    TopoCache::mesh->invDistLaplacian();
        TopoCache::mesh->areaWeightLaplacian(); // perimeter driven actually
    // flip the Laplacian so the diagonal is positive
    laplacian = -laplacian;
    recordEdgeWeights(laplacian);
    
    Eigen::VectorXd rawParity = TopoCache::mesh->buildScalarField(
        [](VertData &v) {
            return ((v.dead == DEAD)? -1.0 :
                        ((v.dead == ALIVE)? 1.0 :
                                            0.0 ));
        });
    
    // right hand side
    Eigen::VectorXd rhs(TopoCache::mesh->verts.size());
    for(uint i=0; i<rhs.size(); i++)
        rhs[i] = 0.0;
    
    // compute areas
    Eigen::VectorXd areas(TopoCache::mesh->verts.size());
    for(uint i=0; i<rhs.size(); i++)
        areas[i] = 0.0;
    // go through tris
    double totalArea = 0.0;
    for(Tri &tri : TopoCache::mesh->tris) {
        Vec3d  p0       = TopoCache::mesh->verts[tri.a].pos;
        Vec3d  p1       = TopoCache::mesh->verts[tri.b].pos;
        Vec3d  p2       = TopoCache::mesh->verts[tri.c].pos;
        double area     = len(cross(p1-p0, p2-p0))/2.0;
        
        areas[tri.a]   += area / 3.0;
        areas[tri.b]   += area / 3.0;
        areas[tri.c]   += area / 3.0;
        totalArea += area;
    }
    
    std::vector<uint> cids = getIsctComponentIds();
    
    //uint nSprings = 0;
    
    // do something special for cut edges
    for_tig_edges([&](uint i, uint j, TIGedge tig) {
        if(tig->iscts.size() > 0) {
            double value = laplacian.coeff(i,j);
            double parity_ij = rawParity[i] - rawParity[j];
            if(parity_ij == 0.0 ||
               cids[i] != cids[j]) { // then zero out edge!
                laplacian.coeffRef(i,j) = 0.0;
                laplacian.coeffRef(j,i) = 0.0;
                // pull out the balancing negative value from the diagonal
                laplacian.coeffRef(i,i) += value;
                laplacian.coeffRef(j,j) += value;
            } else { // put in a spring that pushes i and j apart
                rhs[i] -= value * parity_ij;
                rhs[j] += value * parity_ij;
                //nSprings++;
                // encourages equation x_i - x_j = parity_ij to hold
            }
        }
    });
    //std::cout << "nsprings is " << nSprings << std::endl;
    
    Eigen::SparseMatrix<double> id = TopoCache::mesh->buildIdentity();
    
    // in order to ensure the system is not under-determined,
    // we regularize very slightly using the current parity field
    // To accomplish this, we need to add in a very little mix
    // of the identity matrix and use the
    // diminished parity field as a right hand side
    double regularization = 1.0e-2 / totalArea;
    Eigen::SparseMatrix<double> mixed = laplacian + regularization * id;
    for(uint i=0; i<rhs.size(); i++) {
        //std::cout << "a " << areas[i] << std::endl;
        rhs[i] += regularization * areas[i] * rawParity[i];
    }
    
    Eigen::SimplicialCholesky< Eigen::SparseMatrix<double> >
        factorization(mixed);
    Eigen::VectorXd solvedParity = factorization.solve(rhs);
    
    const double area_coeff = 2.0 * M_PI;
    
    // figure out the cutoff value, component-by-component
    std::vector<bool> visited(solvedParity.size(), false);
    std::map<uint, std::vector<uint> > components = groupComponents(cids);
    std::map<uint, double> wrongPerimInits;
    for(auto &pair : wrongPerimInits)
        pair.second = 0.0;
    for_tig_edges([&](uint i, uint j, TIGedge tig) {
        if(cids[i] == cids[j] &&
           tig->iscts.size() > 0 &&
           rawParity[i] != rawParity[j])
            wrongPerimInits[cids[i]] += tig->weight;
    });
    //std::cout << "num components is " << components.size() << std::endl;
    for(auto &pair : components) {
        auto &ids = pair.second;
        uint cid = cids[ids[0]];
        double areaWrong = 0.0;
        double perimeterWrong = wrongPerimInits[cid];
        std::vector< std::pair<double, uint> > values(ids.size());
        for(uint k=0; k<ids.size(); k++) {
            uint id = ids[k];
            if(rawParity[id] < 0.0) areaWrong += -(rawParity[id] * areas[id]);
            values[k] = std::make_pair(solvedParity[id], id);
        }
        
        std::sort(values.begin(), values.end());
        // sweep cut
        double min_score = area_coeff * areaWrong
                         + perimeterWrong*perimeterWrong;
        double min_cut = -DBL_MAX/10.0;
        double min_perim = perimeterWrong;
        double min_area = areaWrong;
        for(auto pair : values) {
            uint id = pair.second;
            visited[id] = true;
            
            areaWrong += rawParity[id] * areas[id];
            for_neighbors(id,[&](uint j, TIGedge tig) {
                if(cid == cids[j]) {
                    double agree = (visited[j])? 1.0 : -1.0;
                    if(tig->iscts.size() == 0)
                        perimeterWrong -= agree * tig->weight;
                    else {
                        if(rawParity[id] != rawParity[j]) {
                            perimeterWrong += agree * tig->weight;
                        }
                    }
                }
            });
            
            double score = area_coeff * areaWrong
                         + perimeterWrong * perimeterWrong;
            if(score < min_score) {
                min_cut = pair.first;
                min_score = score;
                min_perim = perimeterWrong;
                min_area = areaWrong;
            }
        }
        
        //std::cout << "min score: " << min_score << std::endl;
        //std::cout << "min perim: " << min_perim << std::endl;
        //std::cout << "min area:  " << min_area << std::endl;
        
        for(uint id : ids) {
            double d = solvedParity[id] - min_cut;
            VertData &v = TopoCache::mesh->verts[id];
            v.smoothDead = clamp(d, -1.0, 1.0);
            v.dead = (d > 0.0)? ALIVE : DEAD;
        }
    }
}

/*
template<class VertData, class TriData>
void Mesh<VertData,TriData>::TopoChangeProblem::smoothByMinCut()
{
    Eigen::SparseMatrix<double> laplacian =
    //    TopoCache::mesh->combinatorialLaplacian();
    //    TopoCache::mesh->invDistLaplacian();
        TopoCache::mesh->areaWeightLaplacian(); // perimeter driven actually
    // flip the Laplacian so the diagonal is positive
    laplacian = -laplacian;
    recordEdgeWeights(laplacian);
    
    // compute areas
    Eigen::VectorXd areas(TopoCache::mesh->verts.size());
    for(uint i=0; i<TopoCache::mesh->verts.size(); i++)
        areas[i] = 0.0;
    // go through tris
    double totalArea = 0.0;
    for(Tri &tri : TopoCache::mesh->tris) {
        Vec3d  p0       = TopoCache::mesh->verts[tri.a].pos;
        Vec3d  p1       = TopoCache::mesh->verts[tri.b].pos;
        Vec3d  p2       = TopoCache::mesh->verts[tri.c].pos;
        double area     = len(cross(p1-p0, p2-p0))/2.0;
        
        areas[tri.a]   += area / 3.0;
        areas[tri.b]   += area / 3.0;
        areas[tri.c]   += area / 3.0;
        totalArea += area;
    }
    
    double gamma = 1.0;
    // set up a maxflow graph
    typedef Graph<double,double,double> GraphType;
    GraphType *g = new GraphType(TopoCache::mesh->verts.size(),
                                 TopoCache::edges.size()/2); 
    g->add_node(TopoCache::mesh->verts.size());
    // run through and add connections to source and sink
    for(uint i=0; i<TopoCache::mesh->verts.size(); i++) {
        VertData &v = TopoCache::mesh->verts[i];
        double area = areas[i] * gamma;
        if(v.dead == DEAD)
            g->add_tweights(i, 0.0, area);
        else
            g->add_tweights(i, area, 0.0);
    }
    // run through and add connections between nodes...
    for_tig_edges([&](uint i, uint j, TIGedge tig) {
        if(tig->iscts.size() == 0) // only add uncut edges
            g->add_edge(i,j, tig->weight, tig->weight);
    });
    
    // solve the min-cut/max-flow problem
    g->maxflow();
    
    // read out the results of the cut
    for(uint i=0; i<TopoCache::mesh->verts.size(); i++) {
        VertData &v = TopoCache::mesh->verts[i];
        if(g->what_segment(i) == GraphType::SOURCE) {
            v.smoothDead = 1.0;
            v.dead = ALIVE;
        } else {
            v.smoothDead = -1.0;
            v.dead = DEAD;
        }
    }
    
    delete g;
}
*/

template<class VertData, class TriData>
void Mesh<VertData,TriData>::TopoChangeProblem::cutAndGlue()
{
    std::vector<Tri>            &meshTris   = TopoCache::mesh->tris;
    std::vector<VertData>       &meshVerts  = TopoCache::mesh->verts;
    
    // First off, we introduce new fake "intersection" geometry
    // to split the mesh between the live and dead vertices.
    TopoCache::tris.for_each([&](Tptr t) {
        // get the vertex dead field values...
        DeadState deads[3];
        for(uint k=0; k<3; k++) {
            deads[k] = meshVerts[t->verts[k]->ref].dead;
            ENSURE(deads[k] == ALIVE || deads[k] == DEAD);
        }
        
        // ignore if this isn't a bounary triangle
        if(deads[0] == deads[1] && deads[0] == deads[2])
            return; // continue
        
        // Consider two cases:
        //  1)  The triangle has no intersections yet
        //  2)  The triangle has pre-existing intersections
        if(t->data == nullptr) {
            // create tprob (this sets up existing vertices and edges)
            Tprob tprob = IsctProblem::getTprob(t);
            
            // find two edges that need to be split in this triangle
            // note: there must be exactly two edges based on possible
            //       entries in the deads[] array.
            for(uint k=0; k<3; k++) {
                uint ki = (k+1)%3;
                uint kj = (k+2)%3;
                if(deads[ki] != deads[kj]) {
                    uint i = t->verts[ki]->ref;
                    uint j = t->verts[kj]->ref;
                    if( i > j ) std::swap(i,j);
                    TIGedge tig = getEdge(i,j,false);
                    
                    Eptr e = t->edges[k];
                    GluePt glue = getGluePt(tig, e);
                    Vec3d coord = getSplitCoords(e);
                    
                    tprob->addBoundaryEndpoint(this, nullptr, // hack
                                               e, coord, glue);
                }
            }
        } else // 2) The triangle has pre-existing intersections
        {
            // In this case, we do not introduce a new
            // intersection edge crossing the triangle.
            // Instead, we simply introduce cut points on any
            // of the edges which separate alive vertices from dead vertices,
            // AND which also aren't already split by an intersection point
            Tprob tprob = IsctProblem::getTprob(t);
            
            for(uint k=0; k<3; k++) {
                uint ki = (k+1)%3;
                uint kj = (k+2)%3;
                if(deads[ki] != deads[kj]) {
                    uint i = t->verts[ki]->ref;
                    uint j = t->verts[kj]->ref;
                    if( i > j ) std::swap(i,j);
                    TIGedge tig = getEdge(i,j,false);
                    if(tig->iscts.size() > 0)   continue;
                    
                    Eptr e = t->edges[k];
                    GluePt glue = getGluePt(tig, e);
                    Vec3d coord = getSplitCoords(e);
                    
                    tprob->addBoundaryPointAlone(this, e, coord, glue);
                }
            }
        }
    });
    
    // Having spoofed intersections to separate the live part of the
    // mesh from the dead part, we can now resolve all the intersections
    // to produce a nicely subdivided and glued mesh.
    // WARNING: Because we chose not to insert new intersection edges
    //          in triangles that already had intersections,
    //          we will have slight gaps where some edges will
    //          not have the (e->data) flag set correctly.
    //      This is ok, because we will use a smoothing problem to
    //          more robustly assign triangle alive/dead values...
    IsctProblem::resolveAllIntersections();
    TopoCache::commit();
    
    // First, we identify triangle components as helper information
    std::vector<uint> trics = TopoChangeProblem::getTriComponentIds();
    
    // identify the status of different components
    std::vector<DeadState> tags(meshTris.size(), UNKNOWN);
    for(uint i=0; i<meshTris.size(); i++) {
        DeadState   d0      = meshVerts[meshTris[i].v[0]].dead;
        DeadState   d1      = meshVerts[meshTris[i].v[1]].dead;
        DeadState   d2      = meshVerts[meshTris[i].v[2]].dead;
        DeadState consensus = joinDead(d0, joinDead(d1, d2));
        meshTris[i].data.dead = consensus;
        tags[trics[i]] = joinDead(consensus, tags[trics[i]]);
        
        //meshTris[i].data.deadsmooth = 0.0;
        
        //TriData     &trirep = meshTris[trics[i]].data;
        //trirep.dead = joinDead(consensus, trirep.dead);
    }
    // Then we handle each kind of triangle as follows
    //      ALIVE:  mark ALIVE
    //      DEAD:   mark DEAD
    //      UNKNOWN:
    //          (if in a component tagged UNKNOWN) mark ALIVE
    //          (if in a component tagged ERROR) smooth in a value here...
    std::vector<Tptr> toSmooth; // list of triangles to smooth...
    std::vector<double> smoothVals(meshTris.size(), 0.0);
    TopoCache::tris.for_each([&](Tptr t) {
        uint i = t->ref;
        switch(meshTris[i].data.dead) {
        case ALIVE:     smoothVals[i] = 1.0;
            break;
        case DEAD:      smoothVals[i] = -1.0;
            break;
        case UNKNOWN:
        case ERROR:
        default:
            if(tags[trics[i]] == UNKNOWN) {
                meshTris[i].data.dead = ALIVE;
                smoothVals[i] = 1.0;
            } else {
                meshTris[i].data.dead = ERROR;
                smoothVals[i] = 0.0;
                toSmooth.push_back(t);
            }
            break;
        }
    });
    // perform iterative smoothing here...
    for(uint iter = 0; iter < 50; iter++) {
        for(Tptr t : toSmooth) {
            double sum = 0.0;
            double count = 0.0;
            for(uint k=0; k<3; k++) {
                for(Tptr other : t->edges[k]->tris) {
                    if(other == t)  continue;
                    sum += smoothVals[other->ref];
                    count += 1.0;
                }
            }
            smoothVals[t->ref] = sum / count;
        }
    }
    // read smoothing values back out into triangles
    for(Tptr t : toSmooth) {
        if(smoothVals[t->ref] < 0.0)
            meshTris[t->ref].data.dead = DEAD;
        else
            meshTris[t->ref].data.dead = ALIVE;
    }
    
    // now we delete the dead triangles...
    std::vector<Tptr>   dead_tris;
    TopoCache::tris.for_each([&](Tptr t) {
        if(meshTris[t->ref].data.dead == DEAD)
            dead_tris.push_back(t);
    });
    for(Tptr t : dead_tris)
        TopoCache::deleteTri(t);
    TopoCache::commit();
}


// this function will prevent intersections that
// don't need to be resolved from being resolved.
template<class VertData, class TriData>
void Mesh<VertData,TriData>::TopoChangeProblem::cullUnneededIscts()
{
    // assign every intersection vertex an index
    std::vector<IVptr> iverts(IsctProblem::ivpool.size());
    uint idx = 0;
    IsctProblem::ivpool.for_each([&](IVptr iv) {
        iv->idx = idx;
        iverts[idx] = iv;
        idx++;
    });
    
    // now, we can use the intersection edges with the union-find
    // in order to compute connected components of intersection curves
    UnionFind iv_comps(iverts.size());
    IsctProblem::iepool.for_each([&](IEptr ie) {
        IVptr iv0 = dynamic_cast<IVptr>(ie->ends[0]);
        ENSURE(iv0); // first endpoint should always be isct
        IVptr iv1 = dynamic_cast<IVptr>(ie->ends[1]);
        if(iv1) {
            iv_comps.unionIds(iv0->idx, iv1->idx);
        }
        for(IVptr ivx : ie->interior)
            iv_comps.unionIds(iv0->idx, ivx->idx);
    });
    
    // next we go through and determine whether or not to keep each
    // of these identified components.  We record whether any
    // edge of the component crosses an edge of the mesh with an
    // endpoint that is dead.
    std::vector<bool> hasDead(iverts.size(), false);
    for(uint i=0; i<iverts.size(); i++) {
        IVptr iv = iverts[i];
        uint cid = iv_comps.find(i);
        
        if(iv->glue_marker->edge_tri_type) {
            Eptr e = iv->glue_marker->e;
            uint vi0 = e->verts[0]->ref;
            uint vi1 = e->verts[1]->ref;
            if(TopoCache::mesh->verts[vi0].dead == DEAD ||
               TopoCache::mesh->verts[vi1].dead == DEAD)
            {
                hasDead[cid] = true;
            }
        }
    }
    
    // now we just need to delete all of the isct vertices
    // in components marked as having a dead neighbor
    
    // First, we need to unhook the vertices from the triangle problems
    std::vector<Tprob> tprobsToDelete;
    IsctProblem::tprobs.for_each([&](Tprob tprob) {
        // mark all isct edges here as keep
        for(IEptr ie : tprob->iedges)
            ie->idx = 1;
        // go through and kill isct vertices
        ShortVec<IVptr, 4> ivcopies(tprob->iverts);
        for(IVptr iv : ivcopies) {
            uint cid = iv_comps.find(iv->idx);
            if(!hasDead[cid]) {
                // remove from tprob
                tprob->iverts.erase(iv);
                // mark touching edges for future removal...
                for(GEptr ge : iv->edges) {
                    IEptr ie = dynamic_cast<IEptr>(ge);
                    if(ie)
                        ie->idx = 0;
                }
                // now kill this vertex
                IsctProblem::killIsctVert(iv);
            }
        }
        // now go through and kill isct edges
        ShortVec<IEptr, 2> iecopies(tprob->iedges);
        for(IEptr ie : iecopies) {
            if(ie->idx == 0) { // delete it
                tprob->iedges.erase(ie);
                IsctProblem::killIsctEdge(ie);
            }
        }
        // if we've removed all intersections from the triangle problem
        // then we should probably get rid of the triangle problem
        if(tprob->iverts.size() == 0) {
            tprob->the_tri->data = nullptr; // remove tptr reference
            // kill ov, oe
            for(uint k=0; k<3; k++) {
                IsctProblem::killOrigVert(tprob->overts[k]);
                IsctProblem::killOrigEdge(tprob->oedges[k]);
            }
            // kill the f-ing triangle problem
            tprobsToDelete.push_back(tprob);
        }
    });
    // kill tprobs that we don't need anymore
    for(Tprob tprob : tprobsToDelete)
        IsctProblem::tprobs.free(tprob);
}




template<class VertData, class TriData>
void Mesh<VertData,TriData>::TopoChangeProblem::buildEdgeGraph()
{
    egraph.resize(TopoCache::mesh->verts.size());
    regraph.resize(TopoCache::mesh->verts.size());
    
    for(Tri &tri : TopoCache::mesh->tris) {
        uint vs[3];
        for(uint k=0; k<3; k++) vs[k] = tri.v[k];
        // put in increasing order
        if(vs[0] > vs[1])   std::swap(vs[0], vs[1]);
        if(vs[1] > vs[2])   std::swap(vs[1], vs[2]);
        if(vs[0] > vs[1])   std::swap(vs[0], vs[1]);
        
        getEdge(vs[0], vs[1]);
        getEdge(vs[0], vs[2]);
        getEdge(vs[1], vs[2]);
    }
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::TopoChangeProblem::annotateEdgeGraph()
{
    // run through the original edges to collect intersections
    IsctProblem::oepool.for_each([&](OEptr oe) {
        ENSURE(oe->ends[0]->concrete && oe->ends[1]->concrete);
        uint        i       = oe->ends[0]->concrete->ref;
        uint        j       = oe->ends[1]->concrete->ref;
        if(i > j)   std::swap(i,j);
        TIGedge     e       = getEdge(i,j);
        for(uint k=0; k<oe->interior.size(); k++)
            e->accumulateIsct(oe->interior[k]->coord);
    });
}

template<class VertData, class TriData>
std::vector<uint>
Mesh<VertData,TriData>::TopoChangeProblem::getIsctComponentIds()
{
    UnionFind uf(TopoCache::mesh->verts.size());
    for(uint i=0; i<egraph.size(); i++) {
        for(uint k=0; k<egraph[i].size(); k++) {
            uint j = egraph[i][k]->vid;
            if(egraph[i][k]->iscts.size() <= 0)
                uf.unionIds(i, j);
        }
    }
    return uf.dump();
}

template<class VertData, class TriData>
std::vector<uint>
Mesh<VertData,TriData>::TopoChangeProblem::getTriComponentIds()
{
    UnionFind uf(TopoCache::mesh->tris.size());
    // iterate over edges, using them to union triangles
    TopoCache::edges.for_each([&](Eptr e) {
        if(e->data)     return; // continue if this is an isct edge
        
        Tptr        t0      = e->tris[0];
        for(uint k = 1; k < e->tris.size(); k++) {
            Tptr    t       = e->tris[k];
                    uf.unionIds(t0->ref, t->ref);
        }
    });
    return uf.dump();
}







struct ConformGraphEntry
{
    int  n[3]; // neighbors
    bool f[3]; // whether orientation is same or opposite
    ConformGraphEntry() {
        for(uint k=0; k<3; k++)     n[k] = -1;
    }
    inline void add(uint tid, bool flipped) {
        for(uint k=0; k<3; k++) {
            if(n[k] < 0) {
                n[k] = tid;
                f[k] = flipped;
                break;
            }
        }
    }
};
struct ConformGraph
{
    ConformGraph(uint nTris) : data(nTris) {}
    
    std::vector<ConformGraphEntry> data;
    
    inline void add(uint t0, uint t1, bool flipped) {
        data[t0].add(t1, flipped);
        data[t1].add(t0, flipped);
    }
};

inline bool isFlipped(uint i, uint j, uint *t)
{
    for(uint k=0; k<3; k++) {
        if(t[k] == i) {
            return t[(k+1)%3] != j;
        }
    }
    // This should never happen...
    return true;
}
inline bool isFlipped(uint i, uint j, uint *t0, uint *t1)
{
    // the two triangles ought to have opposite orientations
    // relative to the edge between them!
    return isFlipped(i,j,t0) == isFlipped(i,j,t1);
}

template<class VertData, class TriData>
void Mesh<VertData,TriData>::conformOrientations()
{
    if(tris.size() <= 0)    return;
    
    // build graph from cache.
    ConformGraph cgraph(tris.size());
    NeighborCache ncache = createNeighborCache();
    for(uint i=0; i<ncache.skeleton.size(); i++) {
        for(auto &entry : ncache.skeleton[i]) {
            uint j = entry.vid;
            if(i > j)   continue; // only use each edge once
            
            if(entry.tids.size() != 2)  continue;
            
            uint *t0 = tris[entry.tids[0]].v;
            uint *t1 = tris[entry.tids[1]].v;
            cgraph.add(entry.tids[0], entry.tids[1],
                       isFlipped(i, j, t0, t1));
        }
    }
    // graph built.
    
    // Now, we perform flood fills,
    // flipping vertex orders as necessary
    int begin_at = 0;
    std::vector<bool> visited(tris.size(), false);
    std::vector<bool> flipped(tris.size(), false);
    do {
        std::queue<uint> todo;
        std::vector<uint> component;
        
        visited[begin_at] = true;
        todo.push(begin_at);
        component.push_back(begin_at);
        
        bool orientable = true; // detect unorientable components
        while(!todo.empty()) {
            uint tbase = todo.front();
            bool fbase = flipped[tbase];
            todo.pop();
            
            for(uint k=0; k<3; k++) {
                int     t       = cgraph.data[tbase].n[k];
                bool    f       = cgraph.data[tbase].f[k];
                if(t < 0)       break;
                if(visited[t]) {
                    if(flipped[t] != (f != fbase))
                        orientable = false;
                    continue;
                }
                
                visited[t]      = true;
                todo.push(t);
                component.push_back(t);
                flipped[t]      = f != fbase; // XOR
            }
        }
        
        // don't flip anything in unorientable components of the surface!
        if(!orientable) {
            //std::cout << "unorientable" << std::endl;
            for(uint i : component)
                flipped[i] = false;
        }
        
        // find next component
        begin_at = -1;
        for(uint i=0; i<tris.size(); i++) {
            if(!visited[i]) {
                begin_at = i;
                break;
            }
        }
    } while(begin_at >= 0);
    
    // finally flip all triangles that need flipping
    for(uint i=0; i<tris.size(); i++) {
        if(flipped[i]) {
            Tri &t = tris[i];
            std::swap(t.b, t.c);
            std::swap(t.data.na, t.data.nb);
        }
    }
}
















