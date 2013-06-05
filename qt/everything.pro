
TEMPLATE = app
CONFIG  -= app_bundle

CONFIG  += qt
QT      += opengl

# define the compiler
CLANGXX = clang++
CXX11_FLAGS = -std=c++11 -stdlib=libc++ -Wno-c++11-extensions
QMAKE_CXX = $$CLANGXX
QMAKE_LINK = $$CLANGXX
QMAKE_LFLAGS += $$CXX11_FLAGS -Wl,-no_pie
QMAKE_CXXFLAGS += $$CXX11_FLAGS
QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.7
#QMAKE_CXXFLAGS += -mmacosx-version-min=10.7
# hacking this in one level up directly on generated makefile...

QMAKE_CXXFLAGS += -Winline

DEFINES += DNDEBUG

# need to detect platform in the future...
DEFINES += MACOSX

DIR_COMMON   = ../common
DIR_INCLUDE  = $$DIR_COMMON/include
DIR_MATH     = $$DIR_INCLUDE/math
DIR_UTIL     = $$DIR_INCLUDE/util
DIR_ISCT     = $$DIR_INCLUDE/isct
DIR_GEOM     = $$DIR_INCLUDE/geometry
DIR_RENDER   = $$DIR_INCLUDE/render
DIR_ACCEL    = $$DIR_INCLUDE/accel
#DIR_MAXFLOW  = $$DIR_INCLUDE/maxflow
DIR_FILE     = $$DIR_INCLUDE/file_formats
SUB_DIRS     = $$DIR_FILE/ $$DIR_MATH/ $$DIR_ISCT/  \
               $$DIR_GEOM/ $$DIR_UTIL/ $$DIR_RENDER/ $$DIR_ACCEL/
#               $$DIR_MAXFLOW/

!include( ../makeConstants ) {
    message("oh no failed to include makeConstants")
}

# common lib files
#LIBS        += -L$$GMP_PREFIX/lib -lgmpxx -lgmp
LIBS        += $$GMP_PREFIX/lib/libgmpxx.a $$GMP_PREFIX/lib/libgmp.a
INCLUDEPATH += /opt/local/include/
EIGEN_DIR    = ../
INCLUDEPATH += $$EIGEN_DIR
LIBS        += -L../common/lib/ -ltoptop
POST_TARGETDEPS += $$DIR_COMMON/lib/libtoptop.a
INCLUDEPATH += $$DIR_INCLUDE/ $$SUB_DIRS
# depend path makes sure that qmake itself
# can find these header files in order to process them.
DEPENDPATH  += $$DIR_INCLUDE/ $$SUB_DIRS
HEADERS     += \
    $$DIR_MATH/vec.h \
        $$DIR_MATH/mat.h \
        $$DIR_MATH/bbox.h \
        $$DIR_MATH/transform.h \
        $$DIR_MATH/ray.h \
        $$DIR_MATH/cubicRoots.h \
    $$DIR_UTIL/prelude.h \
        $$DIR_UTIL/memPool.h \
        $$DIR_UTIL/iterPool.h \
        $$DIR_UTIL/maybe.h \
        $$DIR_UTIL/shortVec.h \
        $$DIR_UTIL/unionFind.h \
    $$DIR_ISCT/unsafeRayTriIsct.h \
        $$DIR_ISCT/ext4.h \
        $$DIR_ISCT/fixext4.h \
        $$DIR_ISCT/gmpext4.h \
        $$DIR_ISCT/absext4.h \
        $$DIR_ISCT/quantization.h \
        $$DIR_ISCT/fixint.h \
        $$DIR_ISCT/empty3d.h \
        $$DIR_ISCT/collide.h \
    $$DIR_GEOM/rawMesh.h \
        $$DIR_GEOM/rawMesh.tpp \
        $$DIR_GEOM/mesh.h \
        $$DIR_GEOM/mesh.decl.h \
        $$DIR_GEOM/mesh.tpp \
        $$DIR_GEOM/mesh.topoCache.tpp \
        $$DIR_GEOM/mesh.remesh.tpp \
        $$DIR_GEOM/mesh.isct.tpp \
        $$DIR_GEOM/mesh.collision.tpp \
        $$DIR_GEOM/mesh.distances.tpp \
        $$DIR_GEOM/mesh.operators.tpp \
        $$DIR_GEOM/mesh.topoChange.tpp \
    $$DIR_ACCEL/aabvh.h \
#    $$DIR_MAXFLOW/graph.h \
#        $$DIR_MAXFLOW/block.h \
    $$DIR_RENDER/glIncludes.h \
        $$DIR_RENDER/shader.h \
        $$DIR_RENDER/material.h \
        $$DIR_RENDER/color.h \
        $$DIR_RENDER/scene.h \
        $$DIR_RENDER/camera.h \
        $$DIR_RENDER/light.h \
        $$DIR_RENDER/perspectiveViewport.h \
        $$DIR_RENDER/drawable.h \
        $$DIR_RENDER/drawMesh.h \
        $$DIR_RENDER/drawMesh.tpp \
        $$DIR_RENDER/drawPoints.h \
        $$DIR_RENDER/drawEdges.h \
        $$DIR_RENDER/testObjects.h \
    $$DIR_FILE/files.h

# shared qt files
INCLUDEPATH += src/
HEADERS += src/glcanvas.h src/customWidgets.h \
           src/meshTypes.h
SOURCES += src/glcanvas.cpp src/customWidgets.cpp \
           src/meshTypes.cpp

OBJECTS_DIR  = obj
MOC_DIR      = moc


edit {
    TARGET   = edit
    HEADERS += 
    SOURCES += src/mains/edit.cpp
}

modeler {
    TARGET   = modeler
    HEADERS += 
    SOURCES += src/mains/modeler.cpp
}

# check that all of the files we depend on actually exist
#FILES = $$HEADERS $$SOURCES
#for(file, FILES) {
#    !exists( $$file ) {
#        error( "$$file not found" )
#    }
#}


#CONFIG += warn_on
