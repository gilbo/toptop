// +-------------------------------------------------------------------------
// | drawPoints.cpp
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
#include "drawPoints.h"

DrawPoints::DrawPoints() :
    initialized(false)
{
    shader.load("common/shaders/points.vert",
                "common/shaders/points.frag");
    color = Vec4f(1,1,0,1); // yellow
}

DrawPoints::~DrawPoints()
{
    if(initialized) {
        glDeleteBuffers(1, &vboId);
    }
}

DrawPointsPtr DrawPoints::create()
{
    DrawPointsPtr dpoints(new DrawPoints);
    return dpoints;
}

void DrawPoints::postReload(std::vector<Vec3d> *points)
{
    reload = true;
    
    if(points->size() > 1000000)
        ERROR("trying to show over a million points");
    data_array.resize(points->size());
    uint write = 0;
    for(Vec3d point : (*points)) {
        data_array[write].p = point;
        write++;
    }
}

void DrawPoints::initBuffers()
{
    if(initialized)     return;
    
    glGenBuffers(1, &vboId);
    
    initialized = reload = true;
}

void DrawPoints::loadPointsToGPU()
{
    if(!reload)         return;
    
    glBindBuffer(GL_ARRAY_BUFFER, vboId);
    glBufferData(GL_ARRAY_BUFFER,
                 data_array.size() * sizeof(DrawPointsVBOEntry),
                 &data_array[0],
                 GL_STREAM_DRAW);
}

void DrawPoints::draw()
{
    initBuffers();
    if(reload) {
        loadPointsToGPU();
        reload = false;
    }
    
    glDepthFunc(GL_LEQUAL);
    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, vboId);
    glVertexPointer(3,
                    GL_FLOAT,
                    sizeof(DrawPointsVBOEntry),
                    (void*)(0));
    
    shader.set();
        glColor4fv(color.v);
        glPointSize( 6.0 );
        glDrawArrays(GL_POINTS, 0, data_array.size());
    shader.unset();
    
    glDisableClientState(GL_VERTEX_ARRAY);
    glDepthFunc(GL_LESS);
}




