// +-------------------------------------------------------------------------
// | drawMesh.h
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

#include "glIncludes.h"

#include "drawable.h"
#include "mesh.decl.h"

#include "shader.h"
#include "material.h"

class DrawMesh;
typedef std::shared_ptr<DrawMesh> DrawMeshPtr;

struct DrawMeshVBOEntry {
    Vec3f p;
    Vec3f n;
    Vec4f c;
};


// normal data necessary for drawing
struct RenderVertexData {
    Vec3f normal; // scratch space
    Vec4f color;
};
struct RenderTriangleData {
    Vec3f na, nb, nc; // normal at each vertex
    Vec4f color;
};
//using RenderRawMesh = RawMesh<RenderVertexData, RenderTriangleData>;

// REQUIRES:
//  - MinimalData
//  - RenderData
template<class VertData, class TriData>
void computeFlatNormals(Mesh<VertData, TriData> &mesh);
template<class VertData, class TriData>
void computeAveragedNormals(Mesh<VertData, TriData> &mesh);

class DrawMesh : public Drawable
{
public:
    DrawMesh(); // do not use
    virtual ~DrawMesh();
    
    static DrawMeshPtr create();
    
public: // allow clients to directly manipulate material properties
    PhongMaterial material;
    PhongMaterial wireframeMaterial;
    
public: // manipulate how the mesh gets drawn
    bool showWireframe;
    bool useVertexColor;
    bool useTransparency;
    
public: // load/refresh mesh data to be drawn
    // REQUIRES:
    //  - MinimalData
    //  - NormalData
    template<class VertData, class TriData>
    void postReload(Mesh<VertData, TriData> *mesh);
    //template<class VertData, class TriData>
    //void postReload(RawMesh<VertData, TriData> *mesh);
    // this version loads color data from triangles instead of vertices
    template<class VertData, class TriData>
    void postReloadTriColor(Mesh<VertData, TriData> *mesh);
    
    virtual void draw();

private:
    template<class VertData, class TriData>
    inline void loadTri(TriData &tri, VertData &a, VertData &b, VertData &c);

private:
    std::vector<DrawMeshVBOEntry> data_array;
    uint transparent_offset;
    
private:
    void loadMeshToGPU();
    void initBuffers();
    
    bool initialized;
    bool reload;
    std::function<void()> loadData;
    
    GLuint vboId;
    
    Shader shader;
    Shader depthonly_shader;
    Shader wireframe_shader;
};


#include "drawMesh.tpp"

