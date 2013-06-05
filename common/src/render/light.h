// +-------------------------------------------------------------------------
// | light.h
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

class Light;
class PointLight;
class DirLight;
typedef std::shared_ptr<Light> LightPtr;
typedef std::shared_ptr<PointLight> PointLightPtr;
typedef std::shared_ptr<DirLight> DirLightPtr;

// setup all of the lights with this function
void wireLightsInScene(ScnNodePtr sceneRoot);


class Light : public SceneNode
{
public:
    Light();
    virtual ~Light() {}
    // do not create a Light; create a specific kind of Light
    
    Vec3f ambient;
    Vec3f diffuse;
    Vec3f specular;
    
    // outputs the camera position data in the world coordinate frame
    virtual Vec4f glPos() = 0;
};

class PointLight : public Light
{
public:
    PointLight();
    virtual ~PointLight() {}
    
    static PointLightPtr create();
    static PointLightPtr create(const Vec3f &pos);
    
    void trackTo(const Vec3f &pos);
    void translateBy(const Vec3f &translate);
    
    virtual Vec4f glPos();
private:
    Vec3f position;
};

class DirLight : public Light
{
public:
    DirLight();
    virtual ~DirLight() {}
    
    static DirLightPtr create();
    static DirLightPtr create(const Vec3f &dir);
    static DirLightPtr create(float yaw, float pitch);
    
    void setDirection(const Vec3f &dir);
    void setYawPitch(float yaw, float pitch);
    void setYaw(float degrees);
    void setPitch(float degrees);
    void panLeft(float degrees);
    void panUp(float degrees);
    
    virtual Vec4f glPos();
private:
    void updateTransform();
    void recomputeDirection();
    void recomputeYawPitch();
private:
    float yaw, pitch;
    Vec3f direction;
};



