// +-------------------------------------------------------------------------
// | scene.h
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

#include "transform.h"

#include <set>
#include <memory>

class SceneNode;
typedef std::shared_ptr<SceneNode> ScnNodePtr;

class SceneNode : public std::enable_shared_from_this<SceneNode>
{
public:
    SceneNode() {}
    virtual ~SceneNode() {}
    // please never copy a node.  child/parent data will be invalid if so!
    static ScnNodePtr create();
    // ScnNodePtr clone()  ...  should be implemented in subclasses
    //  clone should return a copy that is completely detached from the scene.
    
public: // data
    Transform transform; // stored in local->root direction
    
public: // parent/child manipulation and access
    void attach(ScnNodePtr child);
    void detach();
    
    const std::set<ScnNodePtr>  getChildren() const;
    std::set<ScnNodePtr>        getChildren();
    const ScnNodePtr            getParent() const;
    ScnNodePtr                  getParent();
private:
    std::set<ScnNodePtr> childNodes;
    std::weak_ptr<SceneNode> parentNode;
};

// coordinate frame change matrices
Transform rootToNode(ScnNodePtr node);
Transform nodeToRoot(ScnNodePtr node);




