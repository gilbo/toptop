// +-------------------------------------------------------------------------
// | light.cpp
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
#include "light.h"

#include "glIncludes.h"

#include <vector>
using std::vector;
using std::set;


Light::Light() :
    ambient(1,1,1), diffuse(1,1,1), specular(1,1,1)
{
}

// --------------------------------------------------------------------------

PointLight::PointLight() :
    position(0,0,0)
{
    ENSURE(ambient == Vec3f(1,1,1)); // development check; TODO: remove
}

PointLightPtr PointLight::create()
{
    PointLightPtr plight(new PointLight);
    return plight;
}
PointLightPtr PointLight::create(const Vec3f &pos)
{
    PointLightPtr plight(new PointLight);
    plight->trackTo(pos);
    return plight;
}

void PointLight::trackTo(const Vec3f &pos)
{
    position = pos;
}
void PointLight::translateBy(const Vec3f &translate)
{
    position += translate;
}

Vec4f PointLight::glPos()
{
    Vec3f worldPos;
    if(getParent())
        worldPos = nodeToRoot(getParent())(position);
    else
        worldPos = position;
    return toHom(worldPos);
}

// --------------------------------------------------------------------------

DirLight::DirLight() :
    yaw(0), pitch(0), direction(0,0,-1)
{
    updateTransform();
}

DirLightPtr DirLight::create()
{
    DirLightPtr dlight(new DirLight);
    return dlight;
}
DirLightPtr DirLight::create(const Vec3f &dir)
{
    DirLightPtr dlight(new DirLight);
    dlight->setDirection(dir);
    return dlight;
}
DirLightPtr DirLight::create(float yaw, float pitch)
{
    DirLightPtr dlight(new DirLight);
    dlight->setYawPitch(yaw, pitch);
    return dlight;
}


void DirLight::updateTransform()
{
    Mat3f rot = rotByYawPitch(yaw, pitch);
    transform = Transform::rotation(rot, transpose(rot));
}
void DirLight::recomputeDirection()
{
    Vec3f horizontal(-std::sin(yaw), 0, -std::cos(yaw));
    Vec3f vertical(0,1,0);
    direction = std::cos(pitch)*horizontal + std::sin(pitch)*vertical;
}
void DirLight::recomputeYawPitch()
{
    Vec3f xzproj = Vec3f(direction.x, 0, direction.z);
    float horizontal = len(xzproj);
    if(horizontal < 1.0e-6) { // near gimble lock singularity...
        if(direction.y > 0)
            pitch = 90.0;
        else
            pitch = -90.0;
    } else {
        normalize(xzproj);
        // no value wrap necessary b/c atan2 returns val in range (-pi,pi)
        yaw = rad2deg(std::atan2(-xzproj.x, -xzproj.z));
        pitch = rad2deg(std::atan2(direction.y, horizontal));
    }
}


void DirLight::setDirection(const Vec3f &dir)
{
    if(len(dir) < 1.0e-6)   return;
    direction = normalized(dir);
    recomputeYawPitch();
    updateTransform();
}
void DirLight::setYawPitch(float yaw_, float pitch_)
{
    yaw = wrap(yaw_, -180.0f, 180.0f);
    pitch = clamp(pitch_, -90.0f, 90.0f);
    recomputeDirection();
    updateTransform();
}
void DirLight::setYaw(float degrees)
{
    yaw = wrap(degrees, -180.0f, 180.0f);
    recomputeDirection();
    updateTransform();
}
void DirLight::setPitch(float degrees)
{
    pitch = clamp(degrees, -90.0f, 90.0f);
    recomputeDirection();
    updateTransform();
}
void DirLight::panLeft(float degrees)
{
    yaw = wrap(yaw + degrees, -180.0f, 180.0f);
    recomputeDirection();
    updateTransform();
}
void DirLight::panUp(float degrees)
{
    pitch = clamp(pitch + degrees, -90.0f, 90.0f);
    recomputeDirection();
    updateTransform();
}

Vec4f DirLight::glPos()
{
    Transform hereUp = nodeToRoot(shared_from_this());
    Vec3f dir = normalized(hereUp(Vec3f(0,0,1)) - hereUp(Vec3f(0,0,0)));
    return Vec4f(dir.x, dir.y, dir.z, 0);
}

// --------------------------------------------------------------------------

void collectSceneLights(vector<LightPtr> &lights, ScnNodePtr node)
{
    // do a pre-fix traversal of the scene tree
    LightPtr light = std::dynamic_pointer_cast<Light>(node);
    if(light)
        lights.push_back(light);
    
    const set<ScnNodePtr> &children = node->getChildren();
    for(set<ScnNodePtr>::iterator it = children.begin();
        it != children.end(); it++)
    {
        collectSceneLights(lights, *it);
    }
}

void wireUpLights(vector<LightPtr> &lights)
{
    ENSURE(lights.size() <= 3);
    
    int i = 0;
    for(;i<int(lights.size()); i++) {
        LightPtr light = lights[i];
        
        GLenum light_id;
        switch(i) {
            case 0:     light_id = GL_LIGHT0;   break;
            case 1:     light_id = GL_LIGHT1;   break;
            case 2:
            default:    light_id = GL_LIGHT2;   break;
        }
        
        Vec4f pos = light->glPos();
        glLightfv(light_id, GL_POSITION, pos.v);
        Vec4f ambient = toHom(light->ambient);
        glLightfv(light_id, GL_AMBIENT, ambient.v);
        Vec4f diffuse = toHom(light->diffuse);
        glLightfv(light_id, GL_DIFFUSE, diffuse.v);
        Vec4f specular = toHom(light->specular);
        glLightfv(light_id, GL_SPECULAR, specular.v);
        // use this variable to turn on the light
        glLightf(light_id, GL_SPOT_EXPONENT, 1.0);
    }
    // turn off the rest of the lights
    for(;i<3; i++) {
        GLenum light_id;
        switch(i) {
            case 0:     light_id = GL_LIGHT0;   break;
            case 1:     light_id = GL_LIGHT1;   break;
            case 2:
            default:    light_id = GL_LIGHT2;   break;
        }
        // turn off the light!!
        glLightf(light_id, GL_SPOT_EXPONENT, 0.0);
    }
}

void wireLightsInScene(ScnNodePtr sceneRoot)
{
    // get a list of all lights in the scene
    vector<LightPtr> lights;
    collectSceneLights(lights, sceneRoot);
    
    wireUpLights(lights);
}



