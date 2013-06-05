// +-------------------------------------------------------------------------
// | collide.h
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

#include "vec.h"

struct VertTriCollisionInput
{
    Vec3d   p;
    Vec3d   pv;
    Vec3d   t[3];
    Vec3d   tv[3];
};
struct VertTriCollisionOutput
{
    uint    nCollisions;
    double  time[2];
    Vec3d   bary[2];
};

// returns true if any collisions occurred
bool collide(VertTriCollisionInput *in, VertTriCollisionOutput *out);


















