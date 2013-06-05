// +-------------------------------------------------------------------------
// | drawable.h
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

#include "scene.h"
#include "light.h"

class Drawable;
typedef std::shared_ptr<Drawable> DrawablePtr;

class Drawable : public SceneNode
{
public:
    virtual ~Drawable() {}
    
    virtual void draw() = 0;
};

// calls the draw routines of the drawables present
// in the specified scene hierarchy.
// No guarantees on draw order are made.
void drawScene(ScnNodePtr sceneRoot);
















