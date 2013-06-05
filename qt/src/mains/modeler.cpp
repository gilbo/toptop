// +-------------------------------------------------------------------------
// | modeler.cpp
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
#include "customWidgets.h"

#include "perspectiveViewport.h"
#include "drawMesh.h"
#include "files.h"
#include "meshTypes.h"
#include "testObjects.h"

#include "drawPoints.h"

#include "color.h"

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <string>
using std::string;
#include <vector>
using std::vector;


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



class ModelerCanvas : public GLCanvas
{
public: // variables
    
    PerspectiveViewport viewport;
    ScnNodePtr          scene;
    
    DrawMeshPtr         drawMesh;
    ModelerMesh         mesh;
    DrawMeshPtr         marker;
    TestObjMesh         markerMesh;
    DrawPointsPtr       drawPoints;
    vector<Vec3d>       points;
    
    
    string              currentTool;
    Timer               toolTimer;
    double              toolRadius;
    double              toolAlpha;
    
    Vec3d               grabPoint; // last point on plane
    Vec4d               grabPlane;
    
    string              topoMode;
    
public:
    void resizeMe(int w, int h) {
        viewport.resize(w,h);
    }
    
private: // toolbrush support
    Ray3d xy2ray(int x, int y) {
        return rootToNode(drawMesh)(viewport.pixelToRay(Vec2f(x,y)));
    }
    ModelerMesh::Isct xy2isct(int x, int y) {
        return mesh.pick(xy2ray(x,y));
    }
    void setDistanceFalloff(ModelerMesh::Isct isct) {
        // compute surface distances
        mesh.minPathDistances(isct);
        // mix the surface distance field with the 3d distances
        // and adjust into an influence map
        Vec3d center = isct.isct;
        double RADIUS = toolRadius;
        mesh.for_verts([center, RADIUS](ModelerVertex &v) {
            double dist3d = len(v.pos - center) / RADIUS;
            double distSurface = v.distField / RADIUS;
            double dist = 1.0f;
            if(std::max(distSurface, dist3d) < 1.0f)
                dist = std::sqrt(dist3d * distSurface);
            v.distField = dist;
        });
    }
    bool frontFacingIsct(ModelerMesh::Isct isct) {
        bool front;
        Vec3d ray_dir = isct.ray.r;
        mesh.accessIsct(isct, [&front, ray_dir](ModelerTriangle &,
          ModelerVertex &a, ModelerVertex &b, ModelerVertex &c) {
            Vec3d normal = cross(b.pos - a.pos, c.pos - a.pos);
            front = (dot(normal, -ray_dir) > 0);
        });
        return front;
    }
    MouseDrag surfaceBrush(
        std::function<void(ModelerMesh::Isct, double)> func
    ) {
        MouseDrag behavior;
        behavior.start = [this](int,int)->bool {
            toolTimer.start();
            return true;
        };
        behavior.finish = [this](int,int) {
            toolTimer.stop();
        };
        behavior.drag = [this, func](int x, int y, int,int)->bool {
            double dt = toolTimer.lap() / 1000.0;
            ModelerMesh::Isct isct = xy2isct(x,y);
            if(!isct.exists) {
                mesh.for_verts([](ModelerVertex &v) {
                    v.color = Vec4f(1,1,1,1);
                });
            } else {
                // update touch-point marker visualization
                marker->transform = Transform::translation(isct.isct);
                // compute distance-derived tool mask
                setDistanceFalloff(isct);
                
                func(isct, dt);
            }
            postMeshUpdate();
            return true;
        };
        return behavior;
    }
    // find plane perpendicular to the line of sight,
    // with given point lying on the plane
    Vec4d viewPlane(Vec3d onPlane) {
        Ray3d line_of_sight = xy2ray(viewport.getWidth()/2,
                                     viewport.getHeight()/2);
        Vec3d normal = line_of_sight.r;
        Vec4d result = toHom(normal);
        result.w = -dot(normal, onPlane);
        return result;
    }
    Vec3d xy2pointOnPlane(int x, int y, Vec4d plane, bool &error) {
        Ray3d ray = xy2ray(x,y);
        // (ray.p + t*ray.r).(plane.normal) + plane.d = 0
        // t*(ray.r . plane.normal) = -(ray.p . plane.normal + plane.d)
        // t = - plane(ray.p,1) / (ray.r . plane.normal)
        Vec3d normal(plane.x, plane.y, plane.z);
        double denom = dot(ray.r, normal);
        double t = 0.0;
        error = (fabs(denom) < 1.0e-7);
        if(!error)    
            t = - dot(plane, toHom(ray.p)) / denom;
        return ray.p + t * ray.r;
    }
    MouseDrag surfaceGrab(
        std::function<void(ModelerMesh::Isct)> start,
        std::function<void(Vec3d, Vec3d)> update,
        std::function<void()> finish
    ) {
        MouseDrag behavior;
        behavior.start = [this, start](int x, int y)->bool {
            ModelerMesh::Isct isct = xy2isct(x,y);
            if(isct.exists) {
                marker->transform = Transform::translation(isct.isct);
                grabPlane = viewPlane(isct.isct);
                grabPoint = isct.isct;
                
                start(isct);
            }
            mesh.for_verts([](ModelerVertex &v) {
                v.color = Vec4f(1,1,1,1);
            });
            postMeshUpdate();
            return isct.exists;
        };
        behavior.finish = [this, finish](int,int) {
            finish();
            postMeshUpdate();
        };
        behavior.drag = [this, update](int x, int y, int,int)->bool {
            bool error;
            Vec3d target = xy2pointOnPlane(x,y, grabPlane, error);
            if(error)       return false; // abort
            
            marker->transform = Transform::translation(target);
            Vec3f displacement = target - grabPoint;
            grabPoint = target;
            
            update(target, displacement);
            
            postMeshUpdate();
            return true;
        };
        return behavior;
    }
    
    
    
    void dumpMeshMotion()
    {
        RawMesh<ModelerVertex,ModelerTriangle> raw = mesh.raw();
        
        Files::FileMesh filemesh = modeler2file(raw);
        if(Files::writeTriMesh("resources/dump/modeler0.off", &filemesh) > 0) {
            ERROR("Unable to write to resources/dump/modeler0.off");
        }
        
        for(uint i=0; i<raw.vertices.size(); i++)
            raw.vertices[i].pos += raw.vertices[i].velocity;
        filemesh = modeler2file(raw);
        if(Files::writeTriMesh("resources/dump/modeler1.off", &filemesh) > 0) {
            ERROR("Unable to write to resources/dump/modeler1.off");
        }
    }
    
public:
    ModelerCanvas(QWidget *parent = NULL) : GLCanvas(parent) {
        scene                   = SceneNode::create();
        viewport.camera         = Camera::create(Vec3f(0,0,5), Vec3f(0,0,0));
                                scene->attach(viewport.camera);
        
        // stick a light onto the camera rig
        DirLightPtr light0      = DirLight::create(30, -45);
                                viewport.camera->attach(light0);
        
        drawMesh                = DrawMesh::create();
            drawMesh->showWireframe     = true;
            drawMesh->wireframeMaterial.diffuse = rgb(180,180,255);
            drawMesh->useTransparency   = true;
        scene->attach(drawMesh);
        
        currentTool             = "none";
        topoMode                = "none";
        
        marker                  = DrawMesh::create();
            marker->material.diffuse    = Vec4f(0.8, 0.8, 0.0, 1.0); // yellow
            markerMesh                  = createCube(0.1);
            markerMesh.for_verts([](TestObjVertex &v) {
                v.color = Vec4f(1.0,1.0,1.0,1.0);
            });
            marker->postReload(&markerMesh);
        scene->attach(marker);
        
        //drawPoints              = DrawPoints::create();
        //    drawPoints->color = Vec4f(0.0,0.7,0.9,1.0); // cyan/turquoise?
        //scene->attach(drawPoints);
        
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
        mouseDrag(mouseGuard(UI::RIGHT, UI::NONE), tumbler);
        
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
        mouseDrag(mouseGuard(UI::RIGHT, UI::SHIFT), dolly);
        
        // Toolbrushes
        // --------------------
        // Inflate / Deflate
        auto in_de_flate = [this](bool reverse) ->
            std::function<void(ModelerMesh::Isct, double)>
        {
          return [this,reverse](ModelerMesh::Isct isct, double dt) {
            double flipped = frontFacingIsct(isct) ? 1.0 : -1.0;
            double lambda  = toolAlpha * flipped * 0.1;
            mesh.for_verts([flipped, dt, lambda, reverse]
              (ModelerVertex &v) {
                double dist = v.distField;
                v.color = brewerRedRamp(dist);
                if(dist >= 1.0f) {
                    v.velocity = Vec3d(0,0,0);
                } else {
                    double distPerSecond = (1.0f-dist) * lambda;
                    v.velocity = distPerSecond * dt * Vec3d(v.normal);
                    if(reverse) v.velocity = -v.velocity;
                    //v.color.w = dist;
                }
            });
            
            // Save out mesh every frame for debugging!
            //dumpMeshMotion();
            // compute collisions, and advance
            mesh.findCollisions();
            mesh.for_verts([](ModelerVertex &v) {
                v.pos += v.velocity;
            });
            // determine what dies
            //mesh.componentMajorityVote();
            mesh.linearSolveVote();
            //mesh.minCutVote();
            
            mesh.remesh();
            mesh.conformOrientations();
            computeAveragedNormals(mesh);
            
            std::cout << "v " << mesh.numVerts()
                     << " t " << mesh.numTris()
                     << " dt " << dt << std::endl;
          };
        };
        auto inflate = in_de_flate(false);
        auto deflate = in_de_flate(true);
        // Inflate
        mouseDrag([this](UI::Buttons button, UI::Modifiers mod)->bool {
            return (button == UI::LEFT && mod == UI::NONE &&
                    currentTool == "Inflate");
        }, surfaceBrush(inflate));
        // Deflate
        mouseDrag([this](UI::Buttons button, UI::Modifiers mod)->bool {
            return (button == UI::LEFT && mod == UI::SHIFT &&
                    currentTool == "Inflate");
        }, surfaceBrush(deflate));
        
        // Smooth
        mouseDrag([this](UI::Buttons button, UI::Modifiers mod)->bool {
            return (button == UI::LEFT && mod == UI::NONE &&
                    currentTool == "Smooth") ||
                   (button == UI::LEFT && mod == UI::CTRL);
        }, surfaceBrush([this]
          (ModelerMesh::Isct, double dt) {
            Eigen::SparseMatrix<double> laplacian =
                mesh.meanValueLaplacian();
            Eigen::SparseMatrix<double> id =
                mesh.buildIdentity();
            Eigen::SparseMatrix<double> mask =
                mesh.buildDiagonal([](ModelerVertex &v) -> double {
                    double dist = v.distField;
                    v.color = brewerBlueRamp(dist);
                    return clamp(1.0 - dist, 0.0, 1.0);
                });
            
            Eigen::MatrixX3d positions = mesh.buildVecField(
                [](ModelerVertex &v) {
                    return v.pos;
                });
            
            double lambda = clamp(0.2 * toolAlpha * dt, 0.0, 0.1);
            Eigen::SparseMatrix<double> smooth =
                id + lambda * mask * laplacian;
            positions = smooth * positions;
            
            mesh.processVecField(&positions, []
                (ModelerVertex &v, const Vec3d &vec) {
                    v.pos = vec;
                });
            mesh.remesh();
            computeAveragedNormals(mesh);
        }));
        
        // Grab
        mouseDrag([this](UI::Buttons button, UI::Modifiers mod)->bool {
            return (button == UI::LEFT && mod == UI::NONE &&
                    currentTool == "Grab");
        }, surfaceGrab(
        [this](ModelerMesh::Isct isct) { // start
            setDistanceFalloff(isct);
            mesh.for_verts([](ModelerVertex &v) {
                v.pastCollisionParity = 0;
            });
        },
        [this](Vec3d, Vec3d displacement) { // update
            Vec3d displace = displacement;
            mesh.for_verts([displace]
              (ModelerVertex &v) {
                double dist = v.distField;
                v.color = brewerRedRamp(dist);
                if(dist >= 1.0f) {
                    v.velocity = Vec3d(0,0,0);
                } else {
                    double frac = (1.0f-dist);
                    v.velocity = frac * displace;
                }
            });
            
            mesh.findCollisions();
            mesh.for_verts([](ModelerVertex &v) {
                v.pos += v.velocity;
            });
            //mesh.previewComponentMajorityVote();
            mesh.previewLinearSolveVote();
            //mesh.previewMinCutVote();
            mesh.for_verts([](ModelerVertex &v) {
                if(v.dead == DEAD)  v.color.w = 0.2f;
            });
            drawMesh->useTransparency = true;
            
            mesh.remesh();
            computeAveragedNormals(mesh);
        },
        [this]() { // finish
            mesh.for_verts([](ModelerVertex &v) {
                v.velocity = Vec3d(0,0,0);
                v.collisions.resize(0);
                v.color = Vec4f(1,1,1,1);
            });
            //mesh.previewComponentMajorityVote();
            mesh.previewLinearSolveVote();
            //mesh.previewMinCutVote();
            mesh.applyDeathField();
            
            mesh.remesh();
            mesh.conformOrientations();
            computeAveragedNormals(mesh);
        }));
        
        // Move Rigid
        mouseDrag([this](UI::Buttons button, UI::Modifiers mod)->bool {
            return (button == UI::LEFT && mod == UI::NONE &&
                    currentTool == "Move Rigid");
        }, surfaceGrab(
        [this](ModelerMesh::Isct isct) { // start
            mesh.existsPathDistances(isct);
            mesh.for_verts([](ModelerVertex &v) {
                v.pastCollisionParity = 0;
            });
        },
        [this](Vec3d, Vec3d displacement) { // update
            mesh.for_verts([displacement]
              (ModelerVertex &v) {
                double dist = v.distField;
                if(dist >= 1.0f) {
                    v.velocity = Vec3d(0,0,0);
                } else {
                    v.velocity = displacement;
                }
            });
            
            mesh.findCollisions();
            mesh.for_verts([](ModelerVertex &v) {
                v.pos += v.velocity;
            });
            //mesh.previewComponentMajorityVote();
            mesh.previewLinearSolveVote();
            //mesh.previewMinCutVote();
            mesh.for_verts([](ModelerVertex &v) {
                if(v.dead == DEAD)  v.color.w = 0.2f;
                else                v.color.w = 1.0f;
            });
            drawMesh->useTransparency = true;
            
            //mesh.remesh();
            //computeAveragedNormals(mesh);
        },
        [this]() { // finish
            mesh.for_verts([](ModelerVertex &v) {
                v.velocity = Vec3d(0,0,0);
                v.collisions.resize(0);
                v.color = Vec4f(1,1,1,1);
            });
            //mesh.previewComponentMajorityVote();
            mesh.previewLinearSolveVote();
            //mesh.previewMinCutVote();
            mesh.applyDeathField();
            
            mesh.remesh();
            mesh.conformOrientations();
            computeAveragedNormals(mesh);
        }));
    }
    virtual ~ModelerCanvas() {}
    
    void draw() {
        glClearColor(0,0,0,1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        viewport.establishPerspective();
        glMatrixMode(GL_MODELVIEW);
        
        wireLightsInScene(scene);
        drawScene(scene);
    }
    
    void postMeshUpdate() {
        drawMesh->postReload(&mesh);
        //mesh.getSurfPoints(&points);
        //drawPoints->postReload(&points);
        postRedisplay();
    }
};




// encapsulate behavior dealing with the file dialog
std::function<void()> loadDialog(
    QMainWindow                     *window,
    string                          dialog_name,
    //string                          default_dir,
    std::function<bool(string)>     action
) {
    return [&]() {
        QString qs = QFileDialog::getOpenFileName(
                        window,
                        dialog_name.c_str(),
                        ".", // directory...
                        "Meshes (*.ifs *.off)");
        string filename = qs.toStdString();
        if(filename != "")
            action(filename);
    };
}
std::function<void()> saveDialog(
    QMainWindow                     *window,
    string                          dialog_name,
    //string                          default_dir,
    std::function<bool(string)>     action
) {
    return [&]() {
        QString qs = QFileDialog::getSaveFileName(
                        window,
                        dialog_name.c_str(),
                        "./untitled.off", // directory...
                        "Meshes (*.ifs *.off)");
        string filename = qs.toStdString();
        if(filename != "")
            action(filename);
    };
}



int main(int argc, char* argv[])
{
    QApplication app(argc, argv);
    initRand();
    
    // set up window...
    QMainWindow     window;
    QWidget         centralWidget;
    window.setCentralWidget(&centralWidget);
    
    // canvas has to be defined a bit early to allow it to be
    // captured in various callbacks...
    ModelerCanvas canvas;
    
    auto loadMesh = [&canvas](string filename) -> bool {
        Files::FileMesh filemesh;
        if(Files::readTriMesh(filename, &filemesh) > 0) {
            ERROR("Unable to load in " + filename);
            return false;
        } else {
            canvas.mesh = file2modeler(filemesh);
            canvas.mesh.remesh(); // disable?
            computeAveragedNormals(canvas.mesh);
            //canvas.mesh.initPointDistribution();
            canvas.postMeshUpdate();
            return true;
        }
    };
    auto saveMesh = [&canvas](string filename) -> bool {
        Files::FileMesh filemesh = modeler2file(canvas.mesh.raw());
        if(Files::writeTriMesh(filename, &filemesh) > 0) {
            ERROR("Unable to write to " + filename);
            return false;
        } else {
            return true;
        }
    };
    auto addMesh = [&canvas](string filename) -> bool {
        Files::FileMesh filemesh;
        if(Files::readTriMesh(filename, &filemesh) > 0) {
            ERROR("Unable to load in " + filename);
            return false;
        } else {
            ModelerMesh addition = file2modeler(filemesh);
            addition.remesh(); // disable?
            computeAveragedNormals(addition);
            canvas.mesh.disjointUnion(addition);
            canvas.postMeshUpdate();
            return true;
        }
    };
    
    QHBoxLayout top_layout(&centralWidget);
    top_layout.setContentsMargins(0,0,0,0);
        // TOOLBOX
        QVBoxLayout toolbox_layout;
        toolbox_layout.setContentsMargins(10,10,0,10);
            vector<string> toolnames;
                toolnames.push_back("Inflate");
                toolnames.push_back("Smooth");
                //toolnames.push_back("Draw +");
                toolnames.push_back("Grab");
                toolnames.push_back("Move Rigid");
            RadioModes tool_mode_radio(toolnames);
            tool_mode_radio.onModeChange([&canvas](string mode) {
                canvas.currentTool = mode;
                cout << "mode changed to " << mode << endl;
            });
        toolbox_layout.addWidget(&tool_mode_radio);
            QFormLayout slider_group;
                Slider radius_slider(2.0, 0.1, 4.0);
                radius_slider.onValueChange([&canvas](double val) {
                    canvas.toolRadius = val;
                    cout << "radius value updated to " << val << endl;
                });
            slider_group.addRow("Radius", &radius_slider);
                Slider alpha_slider(0.15, 0.0, 1.0);
                alpha_slider.onValueChange([&canvas](double val) {
                    canvas.toolAlpha = val;
                    cout << "alpha value updated to " << val << endl;
                });
            slider_group.addRow("Alpha", &alpha_slider);
        toolbox_layout.addLayout(&slider_group);
            QLabel label;
            label.setText(
                "\n"
                "right click to rotate\n"
                "right click + shift to pan\n"
                "scroll wheel to zoom\n"
                "\n"
                "left click to apply tool\n"
                "left + shift to apply \n"
                " opposite tool\n"
                " (only works for inflate)\n"
                "left + cmd (ctrl on linux)"
                " to smooth\n"
            );
        toolbox_layout.addWidget(&label);
        /*    vector<string> toponames(3);
                toponames[0] = "Clay";
                toponames[1] = "Bubble";
                toponames[2] = "None";
            RadioModes topo_mode_radio(toponames);
            topo_mode_radio.onModeChange([&canvas](string mode) {
                canvas.topoMode = mode;
                cout << "Topology Mode Changed To " << mode << endl;
            });
        toolbox_layout.addWidget(&topo_mode_radio);
            CheckBox autoRemeshBox("auto-remeshing");
            autoRemeshBox.onChange([](bool checked) {
                cout << "auto remeshing is checked? " << checked << endl;
            });
        toolbox_layout.addWidget(&autoRemeshBox);*/
        toolbox_layout.addStretch(1); // leave bottom of column blank
    top_layout.addLayout(&toolbox_layout);
        // CANVAS
    top_layout.addWidget(&canvas, 1); // canvas consumes extra space
    
    QMenu *file = window.menuBar()->addMenu("&File");
        MenuItem menuLoad("&Load Mesh", file, &window);
        menuLoad.onCommand(loadDialog(&window,
            "Load Mesh",
            loadMesh // action
        ));
        MenuItem menuAddMesh("&Add Mesh", file, &window);
        menuAddMesh.onCommand(loadDialog(&window,
            "Add Mesh",
            addMesh // action
        ));
        MenuItem menuSave("&Save Mesh", file, &window);
        menuSave.onCommand(saveDialog(&window,
            "Save Mesh",
            saveMesh // action
        ));
        file->addSeparator();
    
    
    string filename = "resources/models/default.off";
    if(argc >= 2) { // if a load file is specified
        filename = argv[1];
    }
    if(!loadMesh(filename)) {
        exit(1);
    }
    
    window.resize(600,500);
    window.show();
    return app.exec();
}




