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
#include "customWidgets.h"

using std::vector;
using std::string;

RadioModes::RadioModes(
    vector<string> modes,
    string init,
    QWidget *parent
) :
    QGroupBox(parent)
{
    group = new QButtonGroup(this);
    layout = new QVBoxLayout();
    
    current_index = 0;
    
    int N = modes.size();
    names.resize(N);
    buttons.resize(N);
    
    for(int i=0; i<N; i++) {
        names[i] = modes[i];
        buttons[i] = new QRadioButton(names[i].c_str());
        group->addButton(buttons[i], i);
        layout->addWidget(buttons[i]);
        
        if(names[i] == init)
            current_index = i;
    }
    setLayout(layout);
    
    changeModeTo(current_index);
    
    // connect something? register handler?
    connect(group, SIGNAL(buttonClicked(int)),
            this,  SLOT(changeModeTo(int)));
}
RadioModes::~RadioModes()
{
    for(uint i=0; i<buttons.size(); i++)
        delete buttons[i];
    delete layout;
    delete group;
}
void RadioModes::setMode(string name)
{
    // try to find this name
    for(uint i=0; i<names.size(); i++) {
        if(name == names[i]) {
            changeModeTo(i);
            break;
        }
    }
}
string RadioModes::getMode()
{
    return names[current_index];
}
void RadioModes::changeModeTo(int id)
{
    if(id < 0 || id >= int(names.size()))       return;
    
    current_index = id;
    buttons[current_index]->setChecked(true);
    
    if(callback) callback(names[current_index]);
}
void RadioModes::onModeChange(std::function<void(std::string)> doThis)
{
    callback = doThis;
    doThis(names[current_index]);
}






inline
double clamp(double val, double minval, double maxval)
{
    return std::max(minval, std::min(maxval, val));
}

inline
double remap(int val, int intmin, int intmax, double dmin, double dmax)
{
    int intrange = intmax - intmin;
    double drange = dmax - dmin;
    // normed value is in range [0,1]
    double normed_val = double(val-intmin) / double(intrange);
    return drange * normed_val + dmin;
}

inline
int remap(double val, double dmin, double dmax, int intmin, int intmax)
{
    int intrange = intmax - intmin;
    double drange = dmax - dmin;
    double normed_val = val / drange;
    return int(normed_val * intrange + 0.49) + intmin;
}

void Slider::setupSlider()
{
    changeVal(currv);
    
    connect(this, SIGNAL(valueChanged(int)),
            this, SLOT(updateFromIntValue(int)));
}

Slider::Slider(
    double minval,
    double maxval,
    QWidget *parent
) :
    QSlider(Qt::Horizontal, parent), minv(minval), maxv(maxval)
{
    currv = (minv + maxv) / 2.0;
    setupSlider();
}
Slider::Slider(
    double initval,
    double minval,
    double maxval,
    QWidget *parent
) :
    QSlider(Qt::Horizontal, parent), minv(minval), maxv(maxval)
{
    currv = clamp(initval, minv, maxv);
    setupSlider();
}
Slider::~Slider()
{
}
void Slider::changeVal(double val)
{
    val = clamp(val, minv, maxv);
    int intmin = minimum();
    int intmax = maximum();
    setValue(remap(val, minv, maxv, intmin, intmax));
    if(callback) callback(val);
}
void Slider::onValueChange(std::function<void(double)> doThis)
{
    callback = doThis;
    callback(currv);
}
void Slider::updateFromIntValue(int discreteVal)
{
    int intmin = minimum();
    int intmax = maximum();
    double oldval = currv;
    currv = remap(discreteVal, intmin, intmax, minv, maxv);
    if(currv != oldval && callback)
        callback(currv);
}







CheckBox::CheckBox(
    std::string label, QWidget *parent
) :
    QCheckBox(label.c_str(), parent)
{
    connect(this, SIGNAL(stateChanged(int)),
            this, SLOT(processChange(int)));
}
CheckBox::~CheckBox()
{
}
void CheckBox::onChange(std::function<void(bool)> doThis)
{
    callback = doThis;
    callback(isChecked());
}
void CheckBox::processChange(int state)
{
    if(callback) callback(state);
}






MenuItem::MenuItem(
    std::string label,
    QMenu *menu,
    QWidget *parent
) :
    QAction(label.c_str(), parent)
{
    menu->addAction(this);
    connect(this, SIGNAL(triggered()),
            this, SLOT(processTrigger()));
}
MenuItem::~MenuItem()
{
}
void MenuItem::onCommand(std::function<void()> doThis)
{
    callback = doThis;
}
void MenuItem::processTrigger()
{
    if(callback) callback();
}







