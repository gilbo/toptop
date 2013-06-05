// +-------------------------------------------------------------------------
// | material.h
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

class PhongMaterial
{
public:
    PhongMaterial();
    PhongMaterial(Vec4f ambientColor,
                  Vec4f diffuseColor,
                  Vec4f specularColor = Vec4f(0,0,0,1),
                  float shininess = 64.0);
    ~PhongMaterial();
    
    void use(); // makes opengl calls
    
public: // data; feel free to change
    Vec4f ambient;
    Vec4f diffuse;
    Vec4f specular;
    float shininess;
};


