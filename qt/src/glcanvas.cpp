// +-------------------------------------------------------------------------
// | glcanvas.cpp
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

#include "glcanvas.h"
//#include <GL/glew.h>

using std::vector;

UI::Modifiers UI::convertFromQtModifiers(Qt::KeyboardModifiers modifiers)
{
    bool isShift = modifiers & Qt::ShiftModifier;
    bool isCtrl  = modifiers & Qt::ControlModifier;
    bool isAlt   = modifiers & Qt::AltModifier;
    return UI::Modifiers((isShift? SHIFT : 0) |
                         (isCtrl?  CTRL  : 0) |
                         (isAlt?   ALT   : 0));
}

UI::Buttons UI::convertFromQtButtons(Qt::MouseButtons buttons)
{
    bool isLeft     = buttons & Qt::LeftButton;
    bool isRight    = buttons & Qt::RightButton;
    bool isMiddle   = buttons & Qt::MiddleButton;
    return UI::Buttons((isLeft?     LEFT    : 0) |
                       (isRight?    RIGHT   : 0) |
                       (isMiddle?   MIDDLE  : 0));
}

GLCanvas::GLCanvas(QWidget *parent)
    : QGLWidget(parent)
{
}

GLCanvas::~GLCanvas()
{
}

void GLCanvas::postRedisplay()
{
    updateGL();
}

int GLCanvas::flipY(int y)
{
    return height - y;
}


void GLCanvas::initializeGL()
{
    /*glewInit();
    
    if(!glewIsSupported("GL_VERSION_2_0")) {
        ERROR("Need to support at least OpenGL version 2.0");
        std::exit(1);
    }*/
}

void GLCanvas::resizeGL(int w, int h)
{
    width = w;
    height = h;
    resizeMe(w,h);
}

void GLCanvas::paintGL()
{
    draw();
}



void GLCanvas::mousePress(MouseGuard guard, MousePress behavior)
{
    mouse_presses.push_back(std::make_pair(guard, behavior));
}
void GLCanvas::mouseDrag(MouseGuard guard, MouseDrag behavior)
{
    mouse_drags.push_back(DragBundle(guard, behavior));
}
void GLCanvas::keyPress(KeyGuard guard, KeyPress behavior)
{
    key_presses.push_back(std::make_pair(guard, behavior));
}
void GLCanvas::wheelSpin(ModGuard guard, WheelSpin behavior)
{
    wheel_spins.push_back(std::make_pair(guard, behavior));
}


void GLCanvas::mousePressEvent(QMouseEvent *event)
{
    UI::Buttons button = UI::convertFromQtButtons(event->button());
    UI::Modifiers mod  = UI::convertFromQtModifiers(qApp->keyboardModifiers());
    int x = event->pos().x();
    int y = flipY(event->pos().y());
    
    // mouse presses
    for(auto &presser : mouse_presses) {
        if(presser.first(button, mod))
            presser.second(x,y);
    }
    
    
    // starts of mouse drags
    for(auto &dragger : mouse_drags) {
        if(dragger.guard(button, mod)) {
            bool result = dragger.behavior.start(x,y); 
            if(result) { // unless the drag is canceled, make it active
                dragger.up_button = button;
                dragger.active = true;
            }
        }
    }
    
    last_x = x;
    last_y = y;
}

void GLCanvas::mouseMoveEvent(QMouseEvent *event)
{
    /*
    UI::Buttons buttons = UI::convertFromQtButtons(event->buttons());
    UI::Modifiers mod   = UI::convertFromQtModifiers(qApp->keyboardModifiers());
    */
    int x = event->pos().x();
    int y = flipY(event->pos().y());
    int dx = x - last_x;
    int dy = y - last_y;
    
    //std::cout << "button " << buttons << " mods " << mod << std::endl;
    
    for(auto &dragger : mouse_drags) {
        if(!dragger.active)             continue;
        
        bool result = dragger.behavior.drag(x, y, dx, dy);
        
        if(!result) {
            dragger.behavior.finish(x, y);
            dragger.active = false;
        }
    }
    
    last_x = x;
    last_y = y;
}

void GLCanvas::mouseReleaseEvent(QMouseEvent *event)
{
    UI::Buttons button = UI::convertFromQtButtons(event->button());
    //UI::Modifiers mod = UI::convertFromQtModifiers(qApp->keyboardModifiers());
    int x = event->pos().x();
    int y = flipY(event->pos().y());
    
    for(auto &dragger : mouse_drags) {
        if(!dragger.active)              continue;
        if(!dragger.up_button == button) continue;
        
        dragger.behavior.finish(x, y);
        dragger.active = false;
    }
    
    last_x = x;
    last_y = y;
}

void GLCanvas::keyPressEvent(QKeyEvent *event)
{
    char key = event->key();
    UI::Modifiers mod = UI::convertFromQtModifiers(qApp->keyboardModifiers());
    
    // key presses
    for(auto &presser : key_presses) {
        if(presser.first(key, mod))
            presser.second();
    }
}

void GLCanvas::wheelEvent(QWheelEvent *event)
{
    double numDegrees = event->delta() / 8.0;
    UI::Modifiers mod = UI::convertFromQtModifiers(qApp->keyboardModifiers());
    
    for(auto &wheeler : wheel_spins) {
        if(wheeler.first(mod))
            wheeler.second(numDegrees);
    }
}

