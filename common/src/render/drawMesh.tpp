// +-------------------------------------------------------------------------
// | drawMesh.tpp
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



template<class VertData, class TriData>
void computeFlatNormals(Mesh<VertData, TriData> &mesh)
{
    mesh.for_tris([](
        TriData &tri,
        VertData &a, VertData &b, VertData &c
    ){
        Vec3f normal = normalized(cross(b.pos-a.pos, c.pos-a.pos));
        
        tri.na = tri.nb = tri.nc = normal;
    });
}

template<class VertData, class TriData>
void computeAveragedNormals(Mesh<VertData, TriData> &mesh)
{
    // zero out the vertex work space
    mesh.for_verts([](VertData &v) {
        v.normal = Vec3f(0,0,0);
    });
    // accumulate into the vertex normals
    mesh.for_tris([](
        TriData &,
        VertData &a, VertData &b, VertData &c
    ){
        Vec3f normal = normalized(cross(b.pos-a.pos, c.pos-a.pos));
        
        a.normal += normal;
        b.normal += normal;
        c.normal += normal;
    });
    // normalize the vertex normals
    mesh.for_verts([](VertData &v) {
        normalize(v.normal);
    });
    // pull back the normals onto the faces
    mesh.for_tris([](
        TriData &tri,
        VertData &a, VertData &b, VertData &c
    ){
        tri.na = a.normal;
        tri.nb = b.normal;
        tri.nc = c.normal;
    });
}


template<class VertData, class TriData>
inline bool transparentTri(
    TriData &tri,
    VertData &a, VertData &b, VertData &c,
    bool useVertexColor
) {
    Vec4f acolor = (useVertexColor)? a.color : tri.color ;
    Vec4f bcolor = (useVertexColor)? b.color : tri.color ;
    Vec4f ccolor = (useVertexColor)? c.color : tri.color ;
    return acolor.w < 1.0f || bcolor.w < 1.0f || ccolor.w < 1.0f;
}

template<class VertData, class TriData>
inline void DrawMesh::loadTri(
    TriData &tri, VertData &a, VertData &b, VertData &c
) {
    Vec4f acolor = (useVertexColor)? a.color : tri.color ;
    Vec4f bcolor = (useVertexColor)? b.color : tri.color ;
    Vec4f ccolor = (useVertexColor)? c.color : tri.color ;
    
    data_array.push_back({a.pos, tri.na, acolor});
    data_array.push_back({b.pos, tri.nb, bcolor});
    data_array.push_back({c.pos, tri.nc, ccolor});
}

template<class VertData, class TriData>
void DrawMesh::postReload(Mesh<VertData, TriData> *mesh) {
    reload = true;
    loadData = [this, mesh]() {
        std::vector<DrawMeshVBOEntry> &data = data_array;
        data_array.clear();
        
        if(!useTransparency) {
            mesh->for_tris([this](
                TriData &tri, VertData &a, VertData &b, VertData &c
            ){ loadTri(tri, a, b, c); });
        } else {
            mesh->for_tris([this](
                TriData &tri,
                VertData &a, VertData &b, VertData &c
            ){
                if(!transparentTri(tri, a,b,c, useVertexColor))
                    loadTri(tri, a, b, c);
            });
            
            transparent_offset = data.size();
            
            mesh->for_tris([this](
                TriData &tri,
                VertData &a, VertData &b, VertData &c
            ){
                if(transparentTri(tri, a,b,c, useVertexColor))
                    loadTri(tri, a, b, c);
            });
        }
    };
}









