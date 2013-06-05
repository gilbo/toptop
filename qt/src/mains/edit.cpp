// +-------------------------------------------------------------------------
// | edit.cpp
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

#include <QtGui>

#include "glcanvas.h"

#include "perspectiveViewport.h"
#include "drawMesh.h"
#include "files.h"
#include "meshTypes.h"
#include "testObjects.h"
#include "drawPoints.h"
#include "drawEdges.h"

#include "color.h"
#include <cfloat>

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <sstream>
using std::stringstream;
using std::string;

// mesh conversion
RawModelerMesh file2modeler(const Files::FileMesh &in) {
    return transduce<ModelerVertex,     ModelerTriangle,
                     Files::FileVertex, Files::FileTriangle>
    (in,
    [](ModelerVertex &out, const Files::FileVertex &in) {
        out.pos = in.pos;
        out.color = Vec4f(1,1,1,1);
    },
    [](ModelerTriangle &out, const Files::FileTriangle &in) {
        out.a = in.a;
        out.b = in.b;
        out.c = in.c;
    });
}
Files::FileMesh modeler2file(const RawModelerMesh &in) {
    return transduce<Files::FileVertex, Files::FileTriangle,
                     ModelerVertex,     ModelerTriangle>
    (in,
    [](Files::FileVertex &out, const ModelerVertex &in) {
        out.pos = in.pos;
    },
    [](Files::FileTriangle &out, const ModelerTriangle &in) {
        out.a = in.a;
        out.b = in.b;
        out.c = in.c;
    });
}


class EditCanvas : public GLCanvas
{
public: // variables
    PerspectiveViewport viewport;
    ScnNodePtr          scene;
    
    DrawMeshPtr         drawMesh;
    ModelerMesh         mesh;
    
    DrawPointsPtr       drawPoints;
    std::vector<Vec3d>  points;
    DrawEdgesPtr                            drawEdges;
    std::vector< std::pair<Vec3d,Vec3d> >   edges;
    
public:
    void resizeMe(int w, int h) {
        viewport.resize(w,h);
    }
    
public:
    EditCanvas(QWidget *parent = NULL) : GLCanvas(parent) {
        scene = SceneNode::create();
        
        // camera positioned at ... looking at ...
        viewport.camera = Camera::create(Vec3d(0,0,5), Vec3d(0,0,0));
        scene->attach(viewport.camera);
        
        DirLightPtr light0 = DirLight::create(30, -45);
        viewport.camera->attach(light0);
        
        drawMesh = DrawMesh::create();
        drawMesh->showWireframe = true;
        drawMesh->wireframeMaterial.diffuse = rgb(180,180,255);
        scene->attach(drawMesh);
        
        drawPoints = DrawPoints::create();
        drawPoints->color = Vec4f(1.0,0.0,0.0,1.0);
        scene->attach(drawPoints);
        
        drawEdges = DrawEdges::create();
        drawEdges->color = Vec4f(1.0,0.0,0.0,1.0);
        scene->attach(drawEdges);
        
        
        // install UI handlers
        // ===================
                
        // Camera
        // -------------------
        // Rotation Control
        MouseDrag tumbler;
        tumbler.start  = [](int,int){ return true; };
        tumbler.drag   = [this](int, int, int dx, int dy)->bool {
            viewport.camera->orbitLeft(dx * 0.5);
            viewport.camera->orbitUp(-dy * 0.5);
            postRedisplay();
            return true;
        };
        tumbler.finish = [](int,int){};
        mouseDrag(mouseGuard(UI::LEFT,UI::NONE), tumbler);
        
        // Zoom Control
        wheelSpin(modGuard(UI::NONE), [this](double deg) {
            double dist = this->viewport.camera->zoomDist();
            this->viewport.camera->zoomIn(deg * dist / 500.0);
            this->postRedisplay();
        });
        
        // Lateral dollying control
        MouseDrag dolly;
        dolly.start  = [](int,int){ return true; };
        dolly.drag   = [this](int, int, int dx, int dy)->bool {
            double dist = viewport.camera->zoomDist();
            double scale = getHeight();
            scale = dist / scale * 1.5;
            viewport.camera->dollyLeft(dx * scale);
            viewport.camera->dollyUp(-dy * scale);
            postRedisplay();
            return true;
        };
        dolly.finish = [](int,int){};
        mouseDrag(mouseGuard(UI::LEFT, UI::SHIFT), dolly);
    }
    virtual ~EditCanvas() {}
    
    void draw() {
        glClearColor(0,0,0,1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        // sets up the viewport and camera transforms
        viewport.establishPerspective();
        glMatrixMode(GL_MODELVIEW);
        
        wireLightsInScene(scene);
        drawScene(scene);
    }
    
    void postMeshUpdate(bool useTriangleColor = false) {
        if(useTriangleColor) {
            drawMesh->useVertexColor = false;
        } else {
            drawMesh->useVertexColor = true;
        }
        drawMesh->postReload(&mesh);
    }
    
    void collisionSetup(string nextfilename)
    {
        Files::FileMesh filemesh;
        if(Files::readTriMesh(nextfilename, &filemesh) > 0) {
            cerr << "Unable to load in " << nextfilename << endl;
            exit(1);
        }
        RawModelerMesh nextFrame(file2modeler(filemesh));
        RawModelerMesh raw = mesh.raw();
        // some brief sanity checks...
        ENSURE(nextFrame.vertices.size() == raw.vertices.size());
        ENSURE(nextFrame.triangles.size() == raw.triangles.size());
        // now parallel run through the raw meshes
        // and annotate velocity data...
        for(uint i=0; i<raw.vertices.size(); i++) {
            Vec3d diff = nextFrame.vertices[i].pos - raw.vertices[i].pos;
            raw.vertices[i].velocity = diff;
        }
        // re-install the mesh with velocity annotations
        mesh = ModelerMesh(raw);
    }
};




int main(int argc, char* argv[])
{
    QApplication app(argc, argv);
    initRand();
    
    if(argc < 2) {
        cout << "Please type 'edit -help' for instructions" << endl;
        exit(0);
    }
    
    int arg_i = 1;
    
    string firstarg = argv[arg_i];
    arg_i++;
    if(firstarg == "-help") // then print the help info and exit
    {
        cout <<
    "Welcome to edit, the command line mesh editor!  Usage:" << endl <<
    "  > edit mesh_file [-command arg0 arg1 ... argn]*" << endl <<
    "for example," << endl <<
    "  > edit bunny.off -scale 2.0 -write bunny.ifs" << endl <<
    "will read in the bunny.off file, scale all the coordinates up by" <<endl<<
    "a factor of 2, and then write out the result to bunny.ifs" << endl <<
    "Options:" << endl <<
    "-scale factor       multiply all model coordinates by factor" << endl <<
    "-write filename     dump the current model to filename" << endl <<
    "-translate x y z    translate all of the point coordinates by x y z" << endl <<
    "-concat filename    append the contents of filename onto the current model" << endl <<
    "-noshow             don't display the final model" << endl <<
    "-rotatex angle      rotate the model angle degrees about the x axis" << endl <<
    "-rotatey angle      rotate the model angle degrees about the y axis" << endl <<
    "-rotatez angle      rotate the model angle degrees about the z axis" << endl <<
    "-jitter epsilon     jitter all of the coordinates by + or - epsilon" << endl <<
    "-remesh             remesh the mesh as it would be if used in modeler" << endl <<
    "-showisct           show where the mesh intersects itself" << endl <<
    "-resolveisct        subdivide to eliminate intersections" << endl;
        cout << endl;
        exit(0);
    }
    
    
    
    EditCanvas canvas;
    canvas.resize(500,500);
    canvas.move(0,0);
    
    Files::FileMesh filemesh;
    
    if(Files::readTriMesh(firstarg, &filemesh) > 0) {
        cerr << "Unable to load in " << firstarg << endl;
        exit(1);
    }
    canvas.mesh = ModelerMesh(file2modeler(filemesh));
    computeFlatNormals(canvas.mesh);
    canvas.postMeshUpdate();
    
    bool noshow = false;
    
    while(arg_i<argc) {
        string command = argv[arg_i];
        arg_i++;
        if(command == "-scale") {
            if(arg_i >= argc) {
                cerr << "no scaling factor provided" << endl;
                exit(1);
            }
            stringstream buf(argv[arg_i]);
            arg_i++;
            double factor;
            buf >> factor;
            canvas.mesh.for_verts([factor](ModelerVertex &v) {
                v.pos *= factor;
            });
            canvas.postMeshUpdate();
        } else
        if(command == "-write") {
            if(arg_i >= argc) {
                cerr << "no output filename provided" << endl;
                exit(1);
            }
            string writefilename = argv[arg_i];
            arg_i++;
            filemesh = modeler2file(canvas.mesh.raw());
            if(Files::writeTriMesh(writefilename, &filemesh) > 0) {
                cerr << "error writing to " << writefilename << endl;
                exit(1);
            }
        } else
        if(command == "-translate") {
            if(arg_i+2 >= argc) {
                cerr << "no x y z translation offset provided" << endl;
                exit(1);
            }
            stringstream bufx(argv[arg_i]);
            stringstream bufy(argv[arg_i+1]);
            stringstream bufz(argv[arg_i+2]);
            arg_i += 3;
            double x,y,z;
            bufx >> x;
            bufy >> y;
            bufz >> z;
            Vec3d trans(x,y,z);
            canvas.mesh.for_verts([trans](ModelerVertex &v) {
                v.pos += trans;
            });
            canvas.postMeshUpdate();
        } else
        if(command == "-concat") {
            if(arg_i >= argc) {
                cerr << "no concatenation filename provided" << endl;
                exit(1);
            }
            string concatfilename = argv[arg_i];
            arg_i++;
            if(Files::readTriMesh(concatfilename, &filemesh) > 0) {
                cerr << "Unable to load in " << concatfilename << endl;
                exit(1);
            }
            ModelerMesh othermesh(file2modeler(filemesh));
            computeFlatNormals(othermesh);
            canvas.mesh.disjointUnion(othermesh);
            canvas.postMeshUpdate();
        } else
        if(command == "-noshow") {
            noshow = true;
        } else
        if(command == "-rotatex") {
            if(arg_i >= argc) {
                cerr << "no rotation degree provided" << endl;
                exit(1);
            }
            stringstream buf(argv[arg_i]);
            arg_i++;
            double degree;
            buf >> degree;
            Mat3d rot = rotByDeg(degree, Vec3d(1,0,0));
            canvas.mesh.for_verts([&rot](ModelerVertex &v) {
                v.pos = rot * v.pos;
            });
            canvas.postMeshUpdate();
        } else
        if(command == "-rotatey") {
            if(arg_i >= argc) {
                cerr << "no rotation degree provided" << endl;
                exit(1);
            }
            stringstream buf(argv[arg_i]);
            arg_i++;
            double degree;
            buf >> degree;
            Mat3d rot = rotByDeg(degree, Vec3d(0,1,0));
            canvas.mesh.for_verts([&rot](ModelerVertex &v) {
                v.pos = rot * v.pos;
            });
            canvas.postMeshUpdate();
        } else
        if(command == "-rotatez") {
            if(arg_i >= argc) {
                cerr << "no rotation degree provided" << endl;
                exit(1);
            }
            stringstream buf(argv[arg_i]);
            arg_i++;
            double degree;
            buf >> degree;
            Mat3d rot = rotByDeg(degree, Vec3d(1,0,0));
            canvas.mesh.for_verts([&rot](ModelerVertex &v) {
                v.pos = rot * v.pos;
            });
            canvas.postMeshUpdate();
        } else
        if(command == "-reflectx") {
            canvas.mesh.for_verts([](ModelerVertex &v) {
                v.pos.x = - v.pos.x;
            });
            canvas.postMeshUpdate();
        } else
        if(command == "-jitter") {
            if(arg_i >= argc) {
                cerr << "no jitter epsilon provided" << endl;
                exit(1);
            }
            cout << "jitter currently a no-op" << endl;
            stringstream buf(argv[arg_i]);
            arg_i++;
            double epsilon;
            buf >> epsilon;
            //main_mesh.jitterMesh(epsilon);
            canvas.postMeshUpdate();
        } else
        if(command == "-remesh") {
            canvas.mesh.remesh();
            computeFlatNormals(canvas.mesh);
            canvas.postMeshUpdate();
        } else
        if(command == "-showisct") {
            canvas.mesh.testingComputeStaticIsct(
                &(canvas.points),
                &(canvas.edges)
            );
            std::cout << "# of intersection points: "
                      << canvas.points.size() << std::endl;
            std::cout << "# of intersection edges:  "
                      << canvas.edges.size() << std::endl;
            canvas.drawPoints->postReload(&(canvas.points));
            canvas.drawEdges->postReload(&(canvas.edges));
        } else
        if(command == "-resolveisct") {
            canvas.mesh.resolveIntersections();
            computeFlatNormals(canvas.mesh);
            canvas.postMeshUpdate();
        } else
        if(command == "-collide") {
            if(arg_i >= argc) {
                cerr << "no next frame mesh filename provided" << endl;
                exit(1);
            }
            canvas.collisionSetup(argv[arg_i]); // installs velocities
            arg_i++;
            // Run collision detection
            canvas.mesh.findCollisions();
            // Then, extract those collisions into the point list
            canvas.points.clear();
            canvas.mesh.for_verts([&](ModelerVertex &v) {
                for(CollisionPacket &cp : v.collisions) {
                    Vec3d pos = v.pos + cp.time * v.velocity;
                    canvas.points.push_back(pos);
                }
            });
            cout << "number of collision points: "
                 << canvas.points.size() << endl;
            // drawing update
            canvas.drawPoints->postReload(&(canvas.points));
            canvas.postMeshUpdate();
        } else
        if(command == "-rawparity") {
            if(arg_i >= argc) {
                cerr << "no next frame mesh filename provided" << endl;
                exit(1);
            }
            canvas.collisionSetup(argv[arg_i]); // installs velocities
            arg_i++;
            // Run collision detection
            canvas.mesh.findCollisions();
            // Then, assign colors to vertices,
            // depending on whether or not they were marked for deletion
            canvas.mesh.for_verts([](ModelerVertex &v) {
                uint nIsct = v.collisions.size();
                v.color = (nIsct % 2 == 0)? rgb(44,123,182) : rgb(215,25,28);
                if(nIsct % 2 == 1)  v.color.w = 0.2;
                //if(nIsct > 1)
                //    v.color = Vec4f(0,1,0,1);
                v.pos += v.velocity;
            });
            
            canvas.drawMesh->useTransparency = true;
            canvas.postMeshUpdate();
        } else
        if(command == "-parity") {
            if(arg_i >= argc) {
                cerr << "no next frame mesh filename provided" << endl;
                exit(1);
            }
            canvas.collisionSetup(argv[arg_i]); // installs velocities
            arg_i++;
            // Run collision detection
            canvas.mesh.findCollisions();
            canvas.mesh.for_verts([](ModelerVertex &v) {
                v.pos += v.velocity;
                v.pastCollisionParity = 0;
            });
            //canvas.mesh.previewComponentMajorityVote();
            canvas.mesh.previewLinearSolveVote();
            //canvas.mesh.previewMinCutVote();
            // colors
            canvas.mesh.for_verts([](ModelerVertex &v) {
                v.color = (v.dead == ALIVE)?   rgb(44,123,182) :
                            ((v.dead == DEAD)? rgb(215,25,28) :
                                               rgb(0,255,0)   );
                if(v.dead == DEAD) v.color.w = 0.2;
            });
            
            canvas.drawMesh->useTransparency = true;
            canvas.postMeshUpdate();
        } else
        if(command == "-smoothparity") {
            if(arg_i >= argc) {
                cerr << "no next frame mesh filename provided" << endl;
                exit(1);
            }
            canvas.collisionSetup(argv[arg_i]); // installs velocities
            arg_i++;
            // Run collision detection
            canvas.mesh.findCollisions();
            canvas.mesh.for_verts([](ModelerVertex &v) {
                v.pos += v.velocity;
                v.pastCollisionParity = 0;
            });
            //canvas.mesh.previewComponentMajorityVote();
            canvas.mesh.previewLinearSolveVote();
            //canvas.mesh.previewMinCutVote();
            // find minimum and maximum values
            double maxv = 0.0;
            canvas.mesh.for_verts([&](ModelerVertex &v) {
                maxv = std::max(maxv, v.smoothDead);
            });
            // colors
            canvas.mesh.for_verts([=](ModelerVertex &v) {
                double val = v.smoothDead / maxv;
                //val = val * 2.0 - 0.4;
                v.color = brewerRedBlueDiverging(val);
            });
            
            canvas.drawMesh->useTransparency = true;
            canvas.postMeshUpdate();
        } else
        if(command == "-topo") {
            if(arg_i >= argc) {
                cerr << "no next frame mesh filename provided" << endl;
                exit(1);
            }
            canvas.collisionSetup(argv[arg_i]); // installs velocities
            arg_i++;
            // Run collision detection
            canvas.mesh.findCollisions();
            // advance mesh
            canvas.mesh.for_verts([](ModelerVertex &v) {
                v.pos += v.velocity;
            });
            // then execute a component majority vote...
            canvas.mesh.linearSolveVote();
            //canvas.mesh.minCutVote();
            // and color triangles based on vertices being
            // dead or not...
            /*canvas.mesh.for_tris([](
                ModelerTriangle &tri,
                ModelerVertex &, ModelerVertex &, ModelerVertex &
            ) {
                byte consensus = tri.dead;
                switch(consensus) {
                case 0:  tri.color = Vec4f(0,0,1,1); break;
                case 1:  tri.color = Vec4f(1,0,0,1); break;
                default: tri.color = Vec4f(0,1,0,1);
                    cout << "hey this can happen" << endl;
                    break;
                }
            });*/
            canvas.mesh.conformOrientations();
            
            computeFlatNormals(canvas.mesh);
            //canvas.postMeshUpdate(true); // load triangle colors
            canvas.postMeshUpdate();
        } else
        if(command == "-showboundaries") {
            // extract boundary edges
            canvas.edges.clear();
            std::vector<uint> tcount;
            int i=-1;
            canvas.mesh.for_edges([&](ModelerVertex &v0,
                                                  ModelerVertex &v1)
            {
                i++;
                canvas.edges.push_back(std::make_pair(v0.pos, v1.pos));
                tcount.push_back(0);
            }, [&](ModelerTriangle &,
                   ModelerVertex&, ModelerVertex&, ModelerVertex&)
            {
                tcount[i]++;
            });
            // cull out all non-boundary edges
            uint write=0;
            for(uint read=0; read<canvas.edges.size(); read++) {
                if(tcount[read] == 1) {
                    canvas.edges[write] = canvas.edges[read];
                    write++;
                }
            }
            cout << "found " << write << " boundary edges on mesh" << endl;
            canvas.edges.resize(write);
            canvas.drawEdges->postReload(&(canvas.edges));
        } else
        if(command == "-showmulti") {
            // extract boundary edges
            canvas.edges.clear();
            std::vector<uint> tcount;
            int i=-1;
            canvas.mesh.for_edges([&](ModelerVertex &v0,
                                                  ModelerVertex &v1)
            {
                i++;
                canvas.edges.push_back(std::make_pair(v0.pos, v1.pos));
                tcount.push_back(0);
            }, [&](ModelerTriangle &,
                   ModelerVertex&, ModelerVertex&, ModelerVertex&)
            {
                tcount[i]++;
            });
            // cull out all non-boundary edges
            uint write=0;
            for(uint read=0; read<canvas.edges.size(); read++) {
                if(tcount[read] > 2) {
                    canvas.edges[write] = canvas.edges[read];
                    write++;
                }
            }
            cout << "found " << write << " >2 valence edges on mesh" << endl;
            canvas.edges.resize(write);
            canvas.drawEdges->postReload(&(canvas.edges));
        } else
        if(command == "-debug") {
            cout << "No current debugging behavior implemented" << endl;
        } else {
            cerr << "did not recognize command " << command << endl;
            exit(1);
        }
    }
    
    
    
    if(noshow)
        return 0;
    
    canvas.show();
    return app.exec();
}
