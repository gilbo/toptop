// +-------------------------------------------------------------------------
// | customWidgets.h
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

#include <QtGui>

#include <QGroupBox>
#include <QRadioButton>
#include <QVBoxLayout>
#include <QButtonGroup>
#include <QCheckBox>
#include <QAction>
#include <QMenu>
#include <QSlider>


// collection of radio buttons
class RadioModes : public QGroupBox
{
    Q_OBJECT
public:
    RadioModes(
        std::vector<std::string> modes,
        std::string init="",
        QWidget *parent = 0
    );
    ~RadioModes();
public:
    void setMode(std::string name);
    std::string getMode();
    void onModeChange(std::function<void(std::string)> doThis);
private slots:
    void changeModeTo(int modeid);
private:
    std::vector<std::string>            names;
    
    int                                 current_index;
    
    std::vector<QRadioButton *>         buttons;
    QVBoxLayout                         *layout;
    QButtonGroup                        *group;
    
    std::function<void(std::string)>    callback;
};


class Slider : public QSlider
{
    Q_OBJECT
private:
    void setupSlider(); // helper for constructors
public:
    Slider(double minval, double maxval, QWidget *parent = 0);
    Slider(double initval, double minval, double maxval, QWidget *parent = 0);
    ~Slider();
public:
    void changeVal(double val);
    void onValueChange(std::function<void(double)> doThis);
private slots:
    void updateFromIntValue(int discreteVal);
private:
    double minv;
    double maxv;
    double currv;
    
    std::function<void(double)> callback;
};


class CheckBox : public QCheckBox
{
    Q_OBJECT
public:
    CheckBox(std::string label, QWidget *parent = 0);
    ~CheckBox();
public:
    void onChange(std::function<void(bool)> doThis);
private slots:
    void processChange(int state);
private:
    std::function<void(bool)> callback;
};


class MenuItem : public QAction
{
    Q_OBJECT
public:
    MenuItem(std::string label, QMenu *menu, QWidget *parent);
    ~MenuItem();
public:
    void onCommand(std::function<void()> doThis);
private slots:
    void processTrigger();
private:
    std::function<void()> callback;
};





