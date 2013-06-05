// +-------------------------------------------------------------------------
// | standard.vert
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

varying vec3 passNormal, passLightDirection, passEyeDirection;

varying vec4 frontColor, backColor;

void main()
{
    vec3 vertPos = (gl_ModelViewMatrix * gl_Vertex).xyz;
    vec3 eyePos = vec3(0.0,0.0,0.0);
    
    passNormal = normalize(gl_NormalMatrix * gl_Normal);
    passEyeDirection = eyePos - vertPos;
    passLightDirection = gl_LightSource[0].position.xyz;
    if(gl_LightSource[0].position.w > 0.0) // if this is a point light
        passLightDirection -= vertPos;
    
    //if(dot(passNormal,passEyeDirection) < 0.0) // make all surfaces two-sided
    //    passNormal = -passNormal;
    
    vec3 lightDir = normalize(passLightDirection);
    float frontDiffuseCos = max(dot( passNormal, lightDir), 0.0);
    float backDiffuseCos  = max(dot(-passNormal, lightDir), 0.0);
    //float diffuseCos = max(diffuseDot, -diffuseDot);
    
    vec4 diffuse = gl_FrontMaterial.diffuse * gl_Color *
                   gl_LightSource[0].diffuse;
    vec4 ambient = gl_FrontMaterial.ambient * gl_Color *
                   gl_LightSource[0].ambient;
    //vec4 globalAmbient = gl_FrontMaterial.ambient * gl_LightModel.ambient;
    
    // do everything except the specular highlight here
    frontColor = vec4(0,0,0,0);
    backColor  = vec4(0,0,0,0);
    if(gl_LightSource[0].spotExponent > 0.0) { // if light on
        frontColor += ambient + frontDiffuseCos * diffuse;
        backColor  += ambient + backDiffuseCos * diffuse;
    }
    
    gl_Position = ftransform();
}