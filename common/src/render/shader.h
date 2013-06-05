// +-------------------------------------------------------------------------
// | shader.h
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

#include "glIncludes.h"

#include <string>
#include "vec.h"

class Shader
{
public:
    Shader();
    Shader(std::string vert_file, std::string frag_file);
    ~Shader();
    
    void load(std::string vert_file, std::string frag_file);
    void loadDefault();
    
    // sets this shader as the current shader with the GL state
    void set();
    void unset();
    
    // returns whether this shader is set to the default
    // can be used to check for a failed load from file;
    // however note that load doesn't actually happen until set()
    bool isDefault() const;

public: // provide access to setting uniform variables
    void setU3f(const std::string &name, const Vec3f &val);

// --------------------------------------------------------------------------
private:
    bool initialized;
    GLuint vert_shader;
    GLuint frag_shader;
    GLuint program;
private: // thunk data
    void resolve_pending_load(); // call to execute lazy thunk
    bool pending; // is a load pending?
    bool pending_default; // is the load pending default or from source?
    std::string pending_vert_src; // shader source text
    std::string pending_frag_src; // "
private:
    void releaseShader();
    // This load function does not perform any automatic cleanup
    bool loadHelper(std::string vert_src, std::string frag_src);
private: // initializing default shader; call setup only once
    static void   setup_default();
    static bool   default_initialized; // track whether we've initialized or not
    static GLuint vert_shader_default;
    static GLuint frag_shader_default;
    static GLuint program_default;
};

