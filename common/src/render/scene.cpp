// +-------------------------------------------------------------------------
// | scene.cpp
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
#include "scene.h"


void SceneNode::attach(ScnNodePtr child)
{
    ENSURE(!(child->parentNode.lock()));
    ScnNodePtr thisNode = shared_from_this();
    child->parentNode = thisNode;
    thisNode->childNodes.insert(child);
}
void SceneNode::detach()
{
    ScnNodePtr thisNode = shared_from_this();
    ScnNodePtr parent = thisNode->parentNode.lock();
    if(!parent)  return; // already detached
    parent->childNodes.erase(thisNode);
    thisNode->parentNode.reset(); // zero out parent pointer
}

ScnNodePtr SceneNode::create()
{
    return std::shared_ptr<SceneNode>(new SceneNode);
}


const std::set<ScnNodePtr>  SceneNode::getChildren() const
{
    return childNodes;
}
std::set<ScnNodePtr>        SceneNode::getChildren()
{
    return childNodes;
}
const ScnNodePtr            SceneNode::getParent() const
{
    return parentNode.lock();
}
ScnNodePtr                  SceneNode::getParent()
{
    return parentNode.lock();
}


Transform rootToNode(ScnNodePtr node)
{
    return nodeToRoot(node).inverse();
}
Transform nodeToRoot(ScnNodePtr node)
{
    ScnNodePtr parent = node->getParent();
    if(parent)
        return nodeToRoot(parent) * node->transform;
    else
        return node->transform;
}

