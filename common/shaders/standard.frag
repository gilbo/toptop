// +-------------------------------------------------------------------------
// | standard.frag
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

vec4 attenuate(vec4 color, float dist)
{
    float factor = 1.0 / (
        gl_LightSource[0].constantAttenuation +
        gl_LightSource[0].linearAttenuation * dist + 
        gl_LightSource[0].quadraticAttenuation * dist * dist
    );
    return color * factor;
}

void main()
{
    vec4 color;
    vec3 normal;
    if(gl_FrontFacing) {
        color = frontColor;
        normal = normalize(passNormal);
    } else {
        color = backColor;
        normal = normalize(-passNormal);
    }
    vec3 lightDir = normalize(passLightDirection);
    vec3 eyeDir = normalize(passEyeDirection);
    
    vec3 refl =  2.0 * dot(lightDir,normal) * normal - lightDir;
    float specularCos = max(dot(refl,eyeDir), 0.0);
    vec4 specular = gl_FrontMaterial.specular * gl_LightSource[0].specular;
    
    if(gl_LightSource[0].spotExponent > 0.0)
        color += pow(specularCos, gl_FrontMaterial.shininess) * specular;
    
    if(gl_LightSource[0].position.w > 0.0)
        gl_FragColor = attenuate(color, length(passLightDirection));
    else
        gl_FragColor = color;
}

