// +-------------------------------------------------------------------------
// | perspectiveViewport.h
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

#include "camera.h"
#include "ray.h"


class PerspectiveViewport {
public:
    PerspectiveViewport();
    ~PerspectiveViewport();

public: // accessors
    int getWidth() const { return width; }
    int getHeight() const { return height; }
    
    CameraPtr   camera;
    float       fov; // angle up to 180 degrees
    float       near, far;  // clip planes for frustum
    
public:
    double aspectRatio() const { return double(width) / double(height); }
    bool valid() const; // checks whether data is sensible

public: // transformations between the screen and 3d space
    Ray3f pixelToRay(const Vec2f &xy) const;
    Vec2f worldToPixel(const Vec3f &p) const;

public: // Drawing related functions
    // The GUI system needs to inform the viewport
    // whenever its size changes.
    void resize(int w, int h);
    // Establish viewport, projection and view transformations
    void establishPerspective();

private:
    int width;
    int height;
};



