// +-------------------------------------------------------------------------
// | testObjects.cpp
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
#include "testObjects.h"

// explicitly instantiate the TestObjMesh type
#include "mesh.tpp" // only basic functionality
template class Mesh<TestObjVertex, TestObjTriangle>;

TestObjVertex newVert(float x, float y, float z) {
    TestObjVertex v;
    v.pos = Vec3f(x,y,z);
    return v;
}

TestObjTriangle newTri(int a, int b, int c) {
    TestObjTriangle t;
    t.a = a;
    t.b = b;
    t.c = c;
    return t;
}

/*
    Cube Pattern
    
    IN 3D
    
        5               1               ^ y
           *---------*                  |
          /         /|                  |
         /      0  / |                  +------>  x
     4  *---------*  |                 /
        |  *      |  *  3             /
        |   7     | /                v  z
        |         |/            
        *---------*             
     6               2          
    
    UNWRAPPED (w/ triangulation)
    
               5 *--------* 1
                 | \      |
                 |   \    |
      5        4 |     \  | 0        1
        *--------*--------*--------*
        | \      | \      | \      |
        |   \    |   \    |   \    |
        |     \  |     \  |     \  |
        *--------*--------*--------*
      7        6 | \      | 2        3
                 |   \    |
                 |     \  |
               7 *--------* 3
                 | \      |
                 |   \    |
                 |     \  |
               5 *--------* 1
    
    TRIANGLE LIST (CCW order)
    0,1,5; 0,5,4;
    6,4,5; 6,5,7;
    2,0,4; 2,4,6;
    3,1,0; 3,0,2;
    3,2,6; 3,6,7;
    1,3,7; 1,7,5;
*/
TestObjMesh createCube(double halfSideLength)
{
    double s = halfSideLength;
    
    RawTestObjMesh cube;
    
    cube.vertices.push_back(newVert( s,  s,  s));
    cube.vertices.push_back(newVert( s,  s, -s));
    cube.vertices.push_back(newVert( s, -s,  s));
    cube.vertices.push_back(newVert( s, -s, -s));
    cube.vertices.push_back(newVert(-s,  s,  s));
    cube.vertices.push_back(newVert(-s,  s, -s));
    cube.vertices.push_back(newVert(-s, -s,  s));
    cube.vertices.push_back(newVert(-s, -s, -s));
    
    // now need to layout triangles
    cube.triangles.push_back(newTri(0, 1, 5));
        cube.triangles.push_back(newTri(0, 5, 4));
    cube.triangles.push_back(newTri(6, 4, 5));
        cube.triangles.push_back(newTri(6, 5, 7));
    cube.triangles.push_back(newTri(2, 0, 4));
        cube.triangles.push_back(newTri(2, 4, 6));
    cube.triangles.push_back(newTri(3, 1, 0));
        cube.triangles.push_back(newTri(3, 0, 2));
    cube.triangles.push_back(newTri(3, 2, 6));
        cube.triangles.push_back(newTri(3, 6, 7));
    cube.triangles.push_back(newTri(1, 3, 7));
        cube.triangles.push_back(newTri(1, 7, 5));
    
    // compute normals
    TestObjMesh meshmesh(cube);
    computeFlatNormals(meshmesh);
    return meshmesh;
}










