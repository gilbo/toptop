// +-------------------------------------------------------------------------
// | shader.cpp
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
#include "shader.h"

#include "prelude.h"

#include <fstream>
using std::ifstream;


// NOTE: setup_default() should be called before anything
// else in all constructors to ensure that the default
// shader has been set up before we do anything else
Shader::Shader() : initialized(false)
{
    loadDefault();
}

Shader::Shader(std::string vert_file, std::string frag_file)
    : initialized(false)
{
    load(vert_file, frag_file);
}

Shader::~Shader()
{
    releaseShader();
}


void Shader::load(std::string vert_file, std::string frag_file)
{
    // now we load source from file...
    pending = true;
    pending_default = true; // in case of failure, fallback
    
    ifstream vertIn(vert_file.c_str());
    if(!vertIn) {
        ERROR("Failed to open vertex shader file: " + vert_file);
        return;
    }
    ifstream fragIn(frag_file.c_str());
    if(!fragIn) {
        ERROR("Failed to open fragment shader file: " + frag_file);
        return;
    }
    
    // read files into strings
    pending_vert_src = ""; // make sure string is clear...
    while(vertIn) {
        std::string temp;
        getline(vertIn, temp);
        pending_vert_src += temp + '\n';
    }
    pending_frag_src = "";
    while(fragIn) {
        std::string temp;
        getline(fragIn, temp);
        pending_frag_src += temp + '\n';
    }
    
    // success!
    pending_default = false;
}


void Shader::loadDefault() {
    pending = true;
    pending_default = true;
}

void Shader::set()
{
    // evaluate thunk
    if(pending) {
        resolve_pending_load();
    }
    // then set the current program
    glUseProgram(program);
}

void Shader::unset()
{
    glUseProgram(0);
}

bool Shader::isDefault() const
{
    return program == program_default;
}



void Shader::setU3f(const std::string &name, const Vec3f &val)
{
    GLint loc = glGetUniformLocation(program, name.c_str());
    glUniform3f(loc, val.x, val.y, val.z);
}





void Shader::resolve_pending_load()
{
    setup_default(); // ensure that the default shader has been created
    
    // make sure to release any currently held shader before loading
    if(initialized && !isDefault())
        releaseShader();
    
    // we have executed this function at least once
    initialized = true;
    
    bool failed_load = false;
    if(!pending_default) {
        failed_load = !loadHelper(pending_vert_src, pending_frag_src);
        pending_vert_src = "";
        pending_frag_src = "";
    }
    
    // if we actually entered with a default load pending,
    // or if we failed to load shaders from file...
    if(pending_default || failed_load) {
        vert_shader = vert_shader_default;
        frag_shader = frag_shader_default;
        program     = program_default;
    }
    
    pending = false;
}


void Shader::releaseShader()
{
    if(isDefault()) return; // no need to release default
    
    // otherwise, need to actually release...
    glDetachShader(program, vert_shader);
    glDetachShader(program, frag_shader);
    glDeleteProgram(program);
    glDeleteShader(vert_shader);
    glDeleteShader(frag_shader);
}

bool Shader::loadHelper(std::string vert_src, std::string frag_src)
{
    const char *vert_src_c_str = vert_src.c_str();
    const char *frag_src_c_str = frag_src.c_str();
    
    // create and setup the shader objects and final program
    vert_shader = glCreateShader(GL_VERTEX_SHADER);
    frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    program     = glCreateProgram();
    glAttachShader(program, vert_shader);
    glAttachShader(program, frag_shader);
    
    bool failed_compile = false;
    int success;
    // compile vertex shader
    glShaderSource(vert_shader, 1, &vert_src_c_str, NULL);
    glCompileShader(vert_shader);
    glGetShaderiv(vert_shader, GL_COMPILE_STATUS, &success);
    if(!success) {
        int len;
        glGetShaderiv(vert_shader, GL_INFO_LOG_LENGTH, &len);
        char *log = new char[len];
        glGetShaderInfoLog(vert_shader, len, NULL, log);
        ERROR(std::string("vertex shader failed to compile.\n"
                          "The log says:\n") +
                          log);
        delete[] log;
        failed_compile = true;
    }
    
    glShaderSource(frag_shader, 1, &frag_src_c_str, NULL);
    glCompileShader(frag_shader);
    glGetShaderiv(frag_shader, GL_COMPILE_STATUS, &success);
    if(!success) {
        int len;
        glGetShaderiv(frag_shader, GL_INFO_LOG_LENGTH, &len);
        char *log = new char[len];
        glGetShaderInfoLog(frag_shader, len, NULL, log);
        ERROR(std::string("fragment shader failed to compile.\n"
                          "The log says:\n") +
                          log);
        delete[] log;
        failed_compile = true;
    }
    
    if(failed_compile)
        goto shader_load_helper_error;
    
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if(!success) {
        int len;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);
        char *log = new char[len];
        glGetProgramInfoLog(program, len, NULL, log);
        ERROR(std::string("default shader program failed to link.\n"
                          "The log says:\n") +
                          log);
        delete[] log;
        goto shader_load_helper_error;
    }
    
    return true; // success
    
    shader_load_helper_error: // error handling jump point.
    // Make sure the shaders and program are properly freed
    glDetachShader(program, vert_shader);
    glDetachShader(program, frag_shader);
    glDeleteProgram(program);
    glDeleteShader(vert_shader);
    glDeleteShader(frag_shader);
    
    return false;
}


// static members
bool   Shader::default_initialized = false;
GLuint Shader::vert_shader_default;
GLuint Shader::frag_shader_default;
GLuint Shader::program_default;

void Shader::setup_default()
{
    if(default_initialized) return;
    
    std::string vertSrc = "\
    void main()\
    {\
        gl_FrontColor = gl_Color;\
        gl_Position = ftransform();\
    } ";
    std::string fragSrc = "\
    void main()\
    {\
        gl_FragColor = gl_Color;\
    } ";
    const char *vertSrcCstr = vertSrc.c_str();
    const char *fragSrcCstr = fragSrc.c_str();
    
    // create and setup the shader objects and final program
    vert_shader_default = glCreateShader(GL_VERTEX_SHADER);
    frag_shader_default = glCreateShader(GL_FRAGMENT_SHADER);
    program_default     = glCreateProgram();
    glAttachShader(program_default, vert_shader_default);
    glAttachShader(program_default, frag_shader_default);
    
    int success;
    // compile vertex shader
    glShaderSource(vert_shader_default, 1, &vertSrcCstr, NULL);
    glCompileShader(vert_shader_default);
    glGetShaderiv(vert_shader_default, GL_COMPILE_STATUS, &success);
    if(!success) {
        char log[2048]; int len;
        glGetShaderInfoLog(vert_shader_default, 2047, &len, log);
        ERROR(std::string("default vertex shader failed to compile.\n"
                          "The log says:\n") +
                          log);
    }
    glShaderSource(frag_shader_default, 1, &fragSrcCstr, NULL);
    glCompileShader(frag_shader_default);
    glGetShaderiv(frag_shader_default, GL_COMPILE_STATUS, &success);
    if(!success) {
        char log[2048]; int len;
        glGetShaderInfoLog(frag_shader_default, 2047, &len, log);
        ERROR(std::string("default fragment shader failed to compile.\n"
                          "The log says:\n") +
                          log);
    }
    
    glLinkProgram(program_default);
    glGetProgramiv(program_default, GL_LINK_STATUS, &success);
    if(!success) {
        char log[2048]; int len;
        glGetProgramInfoLog(program_default, 2047, &len, log);
        ERROR(std::string("default shader program failed to link.\n"
                          "The log says:\n") +
                          log);
    }
    
    default_initialized = true;
}

