// +-------------------------------------------------------------------------
// | perspectiveViewport.cpp
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
#include "perspectiveViewport.h"

#include "glIncludes.h"

PerspectiveViewport::PerspectiveViewport() :
    fov(75), near(0.5), far(100),
    width(0), height(0)
{}

PerspectiveViewport::~PerspectiveViewport()
{
}

bool PerspectiveViewport::valid() const {
    return (camera != NULL) &&
           (fov < 180.0 && fov > 0.0) &&
           (near > 0.0  && far > near) &&
           (width > 0   && height > 0);
}

void PerspectiveViewport::resize(int w, int h)
{
    width = w;
    height = h;
}

void PerspectiveViewport::establishPerspective()
{
    if(!valid())    return;
    
    // viewport transform
    glViewport(0,0,width,height);
    
    // projection transform
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fov, aspectRatio(), near, far);
    
    // view transform
    glMatrixMode(GL_MODELVIEW);
    Mat4f view = transpose(rootToNode(camera).mat());
    glLoadMatrixf(view.data);
}


Ray3f PerspectiveViewport::pixelToRay(const Vec2f &xy) const
{
    if(!valid())    return Ray3f(Vec3f(0,0,0),Vec3f(0,0,0));
    
    // xy to screen coordinates:
    Vec2f screen = (2.0f*xy - Vec2f(width, height)) / float(height);
    // unproject into 3d eye coordinates
    float tanfov = std::tan(deg2rad(fov) / 2.0f);
    Vec2f eyeXY = (tanfov * near) * screen;
    Vec3f eye(eyeXY.x, eyeXY.y, -near);
    // Now construct a ray in eye space coordinates
    // the ray direction is simply eye - (0,0,0), normalized
    Ray3f eyeRay(eye, normalized(eye));
    // then return this ray in world space coordinates
    return nodeToRoot(camera)(eyeRay);
}
Vec2f PerspectiveViewport::worldToPixel(const Vec3f &p) const
{
    if(!valid())    return Vec2f(0,0);
    
    // first, take the point into eye coordinates
    Vec3f eye = rootToNode(camera)(p);
    // then project into screen coordinates
    if(eye.z == 0.0f) return Vec2f(-1,-1); // guard against z == 0
    float tanfov = std::tan(deg2rad(fov) / 2.0f);
    Vec2f screen = Vec2f(eye.x, eye.y) / (tanfov * -eye.z);
    Vec2d pixel  = (float(height) * screen + Vec2f(width, height)) / 2.0f;
    
    return pixel;
}





