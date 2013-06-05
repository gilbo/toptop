// +-------------------------------------------------------------------------
// | meshgen.cpp
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

#include <iostream>
using std::endl;
using std::cout;
using std::cerr;
#include <string>
using std::string;

#include <set>
using std::set;
#include <vector>
using std::vector;

#include "files.h"

#include <cmath>
using std::exp;

// remove any unnecessary vertices from the mesh
void trimVertices(Files::FileMesh *mesh)
{
    int numverts = mesh->vertices.size();
    // identify all of the vertices in use
    set<int> used_vertices;
    for(const auto &tri : mesh->triangles) {
        used_vertices.insert(tri.a);
        used_vertices.insert(tri.b);
        used_vertices.insert(tri.c);
    }
    
    // compact the vertex list and build a remapping function
    vector<int> vertRemap(mesh->vertices.size());
    int write = 0;
    for(int i=0; i<numverts; i++) {
        if(used_vertices.find(i) == used_vertices.end()) {
            vertRemap[i] = -1;
        } else {
            vertRemap[i] = write;
            mesh->vertices[write] = mesh->vertices[i];
            write++;
        }
    }
    mesh->vertices.resize(write);
    
    // remap vertex references in the triangles
    for(auto &tri : mesh->triangles) {
        tri.a = vertRemap[tri.a];
        tri.b = vertRemap[tri.b];
        tri.c = vertRemap[tri.c];
    }
}

void genRectMesh(Files::FileMesh &mesh,
                 int x, int y,
                 double xmin, double xmax,
                 double ymin, double ymax)
{
    int numverts = (x+1)*(y+1);
    int numtris  = x*y*2;
    mesh.vertices.resize(numverts);
    mesh.triangles.resize(numtris);
    
    for(int i=0; i<y+1; i++) {
        for(int j=0; j<x+1; j++) {
            double ys = (double)(i)/(double)(y);
            double xs = (double)(j)/(double)(x);
            mesh.vertices[i*(x+1) + j].pos =
                                Vec3d((1.0-xs)*xmin + xs*xmax,
                                      (1.0-ys)*ymin + ys*ymax,
                                      0.0);
        }
    }
    
    for(int i=0; i<y; i++) {
        for(int j=0; j<x; j++) {
            auto &t0 = mesh.triangles[2*(i*x + j) + 0];
            auto &t1 = mesh.triangles[2*(i*x + j) + 1];
            int lower_left  =  i   *(x+1) + j;
            int lower_right =  i   *(x+1) + j+1;
            int upper_left  = (i+1)*(x+1) + j;
            int upper_right = (i+1)*(x+1) + j+1;
            t0.a = t1.a = lower_left;
            t0.b = lower_right;
            t0.c = t1.b = upper_right;
            t1.c = upper_left;
        }
    }
}

// break all connectivity
void disconnect(Files::FileMesh &result,
                const Files::FileMesh &input)
{
    result.vertices.resize(input.triangles.size()*3);
    result.triangles.resize(input.triangles.size());
    
    for(uint k=0; k<input.triangles.size(); k++) {
        const auto &t = input.triangles[k];
        result.triangles[k] = t;
        uint a = result.triangles[k].a = 3*k + 0;
        uint b = result.triangles[k].b = 3*k + 1;
        uint c = result.triangles[k].c = 3*k + 2;
        result.vertices[a] = input.vertices[t.a];
        result.vertices[b] = input.vertices[t.b];
        result.vertices[c] = input.vertices[t.c];
    }
}

void genDebrisField(Files::FileMesh &mesh,
                    int x, int y,
                    double xmin, double xmax,
                    double ymin, double ymax,
                    double perturb
) {
    Files::FileMesh temp;
    genRectMesh(temp, x, y, xmin, xmax, ymin, ymax);
    
    disconnect(mesh, temp);
    
    // do perturbation
    for(auto &v : mesh.vertices) {
        Vec3d diff(drand(-perturb, perturb),
                   drand(-perturb, perturb),
                   drand(-perturb, perturb));
        v.pos += diff;
    }
}

// F(ree) F(orm) D(eformation)
void gaussianFFD(Files::FileMesh &mesh,
                 Vec3d center,
                 Vec3d dir, // modulate magnitude for amplitude
                 double width)
{
    for(auto &vert : mesh.vertices) {
        Vec3d centered = vert.pos - center;
        Vec3d p_dir = (dot(centered, dir)/dot(dir,dir))*dir;
        Vec3d p_off = centered-p_dir;
        double exponent = dot(p_off,p_off)/(width*width);
        double displacement = exp(-exponent)/width;
        //cout << displacement << endl;
        displacement *= displacement;
        Vec3d final = displacement*dir;
        vert.pos += final;
    }
}

void gaussianCylinderFFD(Files::FileMesh &mesh,
                 Vec3d center,
                 Vec3d axis, // axis of cylinder
                 Vec3d dir, // modulate magnitude for amplitude
                 double width)
{
    normalize(axis);
    for(auto &vert : mesh.vertices) {
        // first compute a virtual center by projecting
        // the point onto the closest point on the cylinder's axis
        Vec3d centered = vert.pos - center;
        Vec3d offset = centered - dot(centered, axis) * axis;
        Vec3d p_dir = (dot(offset, dir)/dot(dir,dir))*dir;
        Vec3d p_off = offset-p_dir;
        double exponent = dot(p_off,p_off)/(width*width);
        double displacement = exp(-exponent)/width;
        //cout << displacement << endl;
        displacement *= displacement;
        Vec3d final = displacement*dir;
        vert.pos += final;
    }
}

int main(int argc, char *argv[])
{
    cout << "Hello, World!" << endl;
    
    Files::FileMesh meshfile;
    
    /*if(readTriMesh(firstarg, &meshfile) != 0) {
        cerr << "Unable to load in " << firstarg;
        exit(1);
    }
    trimVertices(&meshfile);*/
    
    genRectMesh(meshfile, 20, 20, -1.1234, 2.1235, -1.0912, 2.3098);
    
    string writefilename = "models/sheet.off";
    if(Files::writeTriMesh(writefilename, &meshfile) != 0) {
        cerr << "error writing to " << writefilename << endl;
        exit(1);
    }
    
    gaussianFFD(meshfile,
                Vec3d(0.0,0.0,0.0),
                Vec3d(0.0,0.0,1.0) * 2.0,
                1.0);
    
    writefilename = "models/bent_sheet.off";
    if(Files::writeTriMesh(writefilename, &meshfile) != 0) {
        cerr << "error writing to " << writefilename << endl;
        exit(1);
    }
    
    genRectMesh(meshfile, 20, 20, -1.1234, 2.1235, -1.0912, 2.3098);
    gaussianFFD(meshfile,
                Vec3d(-1.1,-1.1,0.0),
                Vec3d(0.0,0.0,1.0)*0.5,
                1.0);
    
    writefilename = "models/cornerBump.off";
    if(Files::writeTriMesh(writefilename, &meshfile) != 0) {
        cerr << "error writing to " << writefilename << endl;
        exit(1);
    }
    
    
    genRectMesh(meshfile, 20, 20, -1.1234, 2.1235, -1.0912, 2.3098);
    gaussianCylinderFFD(meshfile,
                Vec3d(0.2,-1.1,0.0),
                Vec3d(0,1,0),
                Vec3d(0.0,0.0,1.0)*0.5,
                1.0);
    
    writefilename = "models/cylBump.off";
    if(Files::writeTriMesh(writefilename, &meshfile) != 0) {
        cerr << "error writing to " << writefilename << endl;
        exit(1);
    }
    
    genDebrisField(meshfile, 
                   20, 20, -1.1234, 2.1235, -1.0912, 2.3098,
                   0.02);
    
    writefilename = "models/debris.off";
    if(Files::writeTriMesh(writefilename, &meshfile) != 0) {
        cerr << "error writing to " << writefilename << endl;
        exit(1);
    }
    
    return 0;
}