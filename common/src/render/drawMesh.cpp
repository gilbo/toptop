// +-------------------------------------------------------------------------
// | drawMesh.cpp
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
#include "drawMesh.h"


DrawMesh::DrawMesh() :
    showWireframe(false),
    useVertexColor(true),
    useTransparency(false),
    initialized(false)
{
    shader.load("common/shaders/standard.vert",
                "common/shaders/standard.frag");
    depthonly_shader.loadDefault();
    wireframe_shader.load("common/shaders/wireframe.vert",
                          "common/shaders/wireframe.frag");
    wireframeMaterial.diffuse = Vec4f(0.8, 0.0, 0.0, 1.0);
}

DrawMesh::~DrawMesh()
{
    if(initialized) {
        glDeleteBuffers(1, &vboId);
    }
}

DrawMeshPtr DrawMesh::create()
{
    DrawMeshPtr dmesh(new DrawMesh);
    return dmesh;
}

void DrawMesh::draw()
{
    initBuffers();
    //reload = true; // dummy for stress testing
    if(reload) {
        loadMeshToGPU();
        reload = false;
    }
    
    glBindBuffer(GL_ARRAY_BUFFER, vboId);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glVertexPointer(3,
                    GL_FLOAT,
                    sizeof(DrawMeshVBOEntry),
                    (void*)(0));
    glNormalPointer(GL_FLOAT,
                    sizeof(DrawMeshVBOEntry),
                    (void*)(sizeof(Vec3f)));
    glColorPointer(4,
                   GL_FLOAT,
                   sizeof(DrawMeshVBOEntry),
                   (void*)(2 * sizeof(Vec3f)));
    
    shader.set();
        material.use();
        
        if(!useTransparency) {
            glDrawArrays(GL_TRIANGLES, 0, data_array.size());
        } else {
            glDrawArrays(GL_TRIANGLES,
                         0,
                         transparent_offset);
            
            // need to render just to depth buffer first
            shader.unset();
            depthonly_shader.set();
                glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
                glDrawArrays(GL_TRIANGLES,
                             transparent_offset,
                             data_array.size() - transparent_offset);
                glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
            depthonly_shader.unset();
            shader.set();
            
            glDepthFunc(GL_LEQUAL);
            glEnable(GL_BLEND);
            glDrawArrays(GL_TRIANGLES,
                         transparent_offset,
                         data_array.size() - transparent_offset);
            glDisable(GL_BLEND);
            glDepthFunc(GL_LESS);
        }
        
        if(showWireframe) {
            shader.unset();
            wireframe_shader.set();
            glLineWidth(1.5);
            wireframeMaterial.use();
            glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
            glDrawArrays(GL_TRIANGLES, 0, data_array.size());
            glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
            wireframe_shader.unset();
            shader.set();
        }
    shader.unset();
    
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
}

void DrawMesh::loadMeshToGPU()
{
    if(!reload)   return;
    
    // copy data into the data array...
    //data_array.clear();
    //mesh->for_tris([&data_array](
    
    
    /*uint numTris = mesh->triangles.size();
    data_array.resize(numTris * 3);
    
    for(uint i=0; i<numTris; i++) {
        const auto &tri = mesh->triangles[i];
        
        data_array[3 * i + 0].p = mesh->vertices[tri.a].pos;
        data_array[3 * i + 0].n = tri.na;
        data_array[3 * i + 1].p = mesh->vertices[tri.b].pos;
        data_array[3 * i + 1].n = tri.nb;
        data_array[3 * i + 2].p = mesh->vertices[tri.c].pos;
        data_array[3 * i + 2].n = tri.nc;
    }*/
    loadData();
    
    glBindBuffer(GL_ARRAY_BUFFER, vboId);
    glBufferData(GL_ARRAY_BUFFER,
                 data_array.size() * sizeof(DrawMeshVBOEntry),
                 &data_array[0],
                 GL_STREAM_DRAW);
}

void DrawMesh::initBuffers()
{
    if(initialized)     return;
    
    glGenBuffers(1, &vboId);
    initialized = reload = true;
}
