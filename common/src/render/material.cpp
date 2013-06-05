// +-------------------------------------------------------------------------
// | material.cpp
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
#include "material.h"

#include "glIncludes.h"


PhongMaterial::PhongMaterial() :
    ambient(0.1,0.1,0.1,1.0),
    diffuse(0.8,0.8,0.8,1.0),
    specular(0.0,0.0,0.0,1.0),
    shininess(64.0)
{}
PhongMaterial::PhongMaterial(
    Vec4f ambientColor,
    Vec4f diffuseColor,
    Vec4f specularColor,
    float shininessExponent
) :
    ambient(ambientColor),
    diffuse(diffuseColor),
    specular(specularColor),
    shininess(shininessExponent)
{}
PhongMaterial::~PhongMaterial()
{}
    
void PhongMaterial::use()
{
    glMaterialfv(GL_FRONT, GL_AMBIENT, ambient.v);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse.v);
    glMaterialfv(GL_FRONT, GL_SPECULAR, specular.v);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
}
