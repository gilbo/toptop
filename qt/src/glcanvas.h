// +-------------------------------------------------------------------------
// | glcanvas.h
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

#include <QtOpenGL>

#include "prelude.h"

#include <functional>

class GLCanvas;

namespace UI {
    enum Buttons {
        LEFT = 0x1,
        RIGHT = 0x2,
        MIDDLE = 0x4,
    };
    enum Modifiers {
        NONE = 0x0,
        SHIFT = 0x1,
        CTRL = 0x2,
        ALT = 0x4,
    };
    
    Modifiers convertFromQtModifiers(Qt::KeyboardModifiers modifiers);
    Buttons convertFromQtButtons(Qt::MouseButtons button);
    
    class Handler {
    public:
        byte which;
        Modifiers mod;
        virtual bool match(byte qWhich, Modifiers qMod) const {
            return which == qWhich && mod == qMod;
        };
    };
};

// guard types
using ModGuard      = std::function<bool(UI::Modifiers)>;
using MouseGuard    = std::function<bool(UI::Buttons, UI::Modifiers)>;
using KeyGuard      = std::function<bool(char, UI::Modifiers)>;

// convenience guard constructors
inline
ModGuard modGuard(UI::Modifiers modifiers) {
    return [=](UI::Modifiers m) {
        return m == modifiers;
    };
}
inline
MouseGuard mouseGuard(UI::Buttons buttons, UI::Modifiers modifiers) {
    return [=](UI::Buttons b, UI::Modifiers m) {
        return b == buttons && m == modifiers;
    };
}
inline
KeyGuard keyGuard(char key, UI::Modifiers modifiers) {
    return [=](char k, UI::Modifiers m) {
        return k == key && m == modifiers;
    };
}

// behavior types
using MousePress = std::function<void(int,int)>; // x, y position
using WheelSpin  = std::function<void(double)>; // degrees spun
using KeyPress   = std::function<void()>;
struct MouseDrag {
    // x, y -> false == cancel drag
    std::function<bool(int,int)>            start;
    // x, y, dx, dy -> false == cancel drag
    std::function<bool(int,int,int,int)>    drag;
    // x, y
    std::function<void(int,int)>            finish;
};


// Subclass this class in order to define a new tool/application
class GLCanvas : public QGLWidget
{
    Q_OBJECT
public:
    GLCanvas(QWidget *parent = NULL);
    virtual ~GLCanvas();
    
    void postRedisplay(); // tell screen that it needs to redraw!
    
public: // register event handlers/behaviors
    void mousePress(MouseGuard guard, MousePress behavior);
    void mouseDrag(MouseGuard guard, MouseDrag behavior);
    void keyPress(KeyGuard guard, KeyPress behavior);
    void wheelSpin(ModGuard guard, WheelSpin behavior);
    
protected: // override to implement functionality
    virtual void draw() = 0;
    virtual void resizeMe(int w, int h) = 0; // cannot use name resize

protected: // useful data
    int getWidth() { return width; }
    int getHeight() { return height; }
    
// *         CLIENT          *
// ***************************
// *     IMPLEMENTATION      *
    
private:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void wheelEvent(QWheelEvent *event);
    
private:
    int flipY(int y);
    
private:
    int width;
    int height;
    int last_x;
    int last_y;

private:
    struct DragBundle {
        MouseGuard  guard;
        MouseDrag   behavior;
        UI::Buttons up_button;  // which button to finish drag on release of
        bool        active;     // store whether drag is ongoing
        DragBundle(MouseGuard g, MouseDrag b) :
            guard(g), behavior(b), active(false) {}
    };
private:
    std::vector< std::pair<MouseGuard, MousePress> >    mouse_presses;
    std::vector< DragBundle >                           mouse_drags;
    std::vector< std::pair<KeyGuard, KeyPress> >        key_presses;
    std::vector< std::pair<ModGuard, WheelSpin> >       wheel_spins;
};


