// +-------------------------------------------------------------------------
// | drawEdges.h
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
#include "shader.h"

#include <vector>


class DrawEdges;
typedef std::shared_ptr<DrawEdges> DrawEdgesPtr;

struct DrawEdgesVBOEntry {
    Vec3f p;
    //Vec3f c;
};

class DrawEdges : public Drawable
{
public:
    DrawEdges(); // do not use
    virtual ~DrawEdges();
    
    static DrawEdgesPtr create();
    
public: // allow the color to be directly manipulated
    Vec4f   color;
    
public:
    void postReload(std::vector< std::pair<Vec3d,Vec3d> > *edges);
    
    virtual void draw();

private:
    std::vector<DrawEdgesVBOEntry> data_array;
    
private:
    void loadEdgesToGPU();
    void initBuffers();
    
    bool initialized;
    bool reload;
    
    GLuint vboId;
    
    Shader shader;
};

