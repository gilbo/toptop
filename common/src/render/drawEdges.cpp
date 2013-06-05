// +-------------------------------------------------------------------------
// | drawEdges.cpp
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
#include "drawEdges.h"

DrawEdges::DrawEdges() :
    initialized(false)
{
    shader.load("common/shaders/edges.vert",
                "common/shaders/edges.frag");
    color = Vec4f(1,1,0,1); // yellow
}

DrawEdges::~DrawEdges()
{
    if(initialized) {
        glDeleteBuffers(1, &vboId);
    }
}

DrawEdgesPtr DrawEdges::create()
{
    DrawEdgesPtr dedges(new DrawEdges);
    return dedges;
}

void DrawEdges::postReload(std::vector< std::pair<Vec3d,Vec3d> > *edges)
{
    reload = true;
    
    data_array.resize(edges->size() * 2);
    uint write = 0;
    for(auto edge : (*edges)) {
        data_array[write].p = edge.first;
        write++;
        data_array[write].p = edge.second;
        write++;
    }
}

void DrawEdges::initBuffers()
{
    if(initialized)     return;
    
    glGenBuffers(1, &vboId);
    
    initialized = reload = true;
}

void DrawEdges::loadEdgesToGPU()
{
    if(!reload)         return;
    
    glBindBuffer(GL_ARRAY_BUFFER, vboId);
    glBufferData(GL_ARRAY_BUFFER,
                 data_array.size() * sizeof(DrawEdgesVBOEntry),
                 &data_array[0],
                 GL_STREAM_DRAW);
}

void DrawEdges::draw()
{
    initBuffers();
    if(reload) {
        loadEdgesToGPU();
        reload = false;
    }
    
    glDepthFunc(GL_LEQUAL);
    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, vboId);
    glVertexPointer(3,
                    GL_FLOAT,
                    sizeof(DrawEdgesVBOEntry),
                    (void*)(0));
    
    shader.set();
        glColor4fv(color.v);
        glLineWidth( 3.0 );
        glDrawArrays(GL_LINES, 0, data_array.size());
    shader.unset();
    
    glDisableClientState(GL_VERTEX_ARRAY);
    glDepthFunc(GL_LESS);
}

