# Warning: This makefile rampantly assumes it is being run in UNIX
#   Also, there may be some gnu tool chain assumptions: gcc and gmake?

# ***********
# * PREFACE *
# ***********------+
# | Subdirectories |
# +----------------+
# declare subdirectories of the source directory
SUBDIRECTORIES := util file_formats math isct geometry render accel
#                  maxflow

# make sure the subdirectories are mirrored in
# the obj/ debug/ and depend/ directories
# (HACK: use this dummy variable to get access to shell commands)
SHELL_HACK := $(shell mkdir -p bin lib)
SHELL_HACK := $(shell mkdir -p $(addprefix obj/,$(SUBDIRECTORIES)))
SHELL_HACK := $(shell mkdir -p $(addprefix debug/,$(SUBDIRECTORIES)))
SHELL_HACK := $(shell mkdir -p $(addprefix depend/,$(SUBDIRECTORIES)))
# also make a directory to expose headers in
SHELL_HACK := $(shell mkdir -p $(addprefix include/,$(SUBDIRECTORIES)))

# +----------+
# | Platform |
# +----------+-------------------
# | evaluates to 'Darwin' on mac
# | evaluates to ??? on ???
# +------------------------------
PLATFORM := $(shell uname)

-include ../makeConstants

# *********
# * TOOLS *
# *********--+
# | Programs |
# +----------+
CPP11_FLAGS := -std=c++11 -stdlib=libc++ -Wno-c++11-extensions
#CC  := gcc
#CXX := g++
CC  := clang
CXX := clang++
RM  := rm
CP  := cp

# +--------------+
# | Option Flags |
# +--------------+
# these defines remove extraneous dependencies in inherited code
#CONFIG    := -DDDF_NOQHULL -DMSSTANDALONE
# this line correctly compiles Tetgen as a library
#CONFIG    := $(CONFIG) -DTETLIBRARY

# also let the preprocessor know which platform we're on
ifeq ($(PLATFORM),Darwin)
  CONFIG  := $(CONFIG) -DMACOSX
else
  CONFIG  := $(CONFIG)
endif

# include paths (flattens the hierarchy for include directives)
INC       := -I src/ $(addprefix -I src/,$(SUBDIRECTORIES))
# Place the location of GMP header files here
GMPINC    := -I $(GMP_PREFIX)/include
#ANNINC    := -I ../ann_1.1.2/include
EIGEN_DIR := ../
EIGEN_INC := -I$(EIGEN_DIR)
INC       := $(INC) $(GMPINC) $(ANNINC) $(EIGEN_INC)

# use the second line to disable profiling instrumentation
# PROFILING := -pg
PROFILING :=
CCFLAGS   := -Wall $(INC) $(CONFIG) -O2 -DNDEBUG $(PROFILING) -Winline
CXXFLAGS  := $(CCFLAGS) $(CPP11_FLAGS)
CCDFLAGS  := -Wall $(INC) $(CONFIG) -ggdb
CXXDFLAGS := $(CCDFLAGS)

# Place the location of GMP libraries here
#GMPLD     := -L$(GMP_PREFIX)/lib -lgmpxx -lgmp
# static versionÉ
GMPLD     := $(GMP_PREFIX)/lib/libgmpxx.a $(GMP_PREFIX)/lib/libgmp.a
#ANNLD     := -L../ann_1.1.2/lib -lANN

#GL_LD     := -framework OpenGL -lGLEW
GL_LD     := -framework OpenGL

LINK         := $(CXXFLAGS) $(GMPLD) $(GL_LD) -Wl,-no_pie
LINKD        := $(CXXDFLAGS) $(GMPLD) $(GL_LD) -Wl,-no_pie

# ***********************
# * SOURCE DECLARATIONS *
# ***********************-----------------+
# | SRCS defines a generic bag of sources |
# +---------------------------------------+
MATH_SRCS    := cubicRoots
UTIL_SRCS    := timer log
ISCT_SRCS    := empty3d collide quantization
GEOM_SRCS    := mesh
ACCEL_SRCS   := 
#MAXFLOW_SRCS := graph maxflow
RENDER_SRCS  := shader material \
                scene camera light \
                perspectiveViewport \
                drawable drawMesh drawPoints drawEdges testObjects
FILE_SRCS    := files ifs off
SRCS         := \
    $(addprefix math/,$(MATH_SRCS))\
    $(addprefix util/,$(UTIL_SRCS))\
    $(addprefix isct/,$(ISCT_SRCS))\
    $(addprefix geometry/,$(GEOM_SRCS))\
    $(addprefix accel/,$(ACCEL_SRCS))\
    $(addprefix render/,$(RENDER_SRCS))\
    $(addprefix file_formats/,$(FILE_SRCS))


# +-----------------------------------+
# | HEADERS defines headers to export |
# +-----------------------------------+
MATH_HEADERS      := vec.h mat.h bbox.h transform.h ray.h cubicRoots.h
UTIL_HEADERS      := prelude.h memPool.h iterPool.h maybe.h shortVec.h \
                     unionFind.h
ISCT_HEADERS      := unsafeRayTriIsct.h \
                     ext4.h fixext4.h gmpext4.h absext4.h \
                     quantization.h fixint.h \
                     empty3d.h collide.h \
                     triangle.h
GEOM_HEADERS      := rawMesh.h rawMesh.tpp \
                     mesh.h mesh.decl.h \
                     mesh.tpp mesh.topoCache.tpp \
                     mesh.remesh.tpp mesh.isct.tpp mesh.collision.tpp \
                     mesh.distances.tpp mesh.operators.tpp \
                     mesh.topoChange.tpp
ACCEL_HEADERS     := aabvh.h
#MAXFLOW_HEADERS   := block.h graph.h
RENDER_HEADERS    := glIncludes.h shader.h material.h color.h \
                     scene.h camera.h light.h \
                     perspectiveViewport.h \
                     drawable.h drawMesh.h drawMesh.tpp \
                     drawPoints.h drawEdges.h \
                     testObjects.h
FILE_HEADERS      := files.h
HEADERS           := \
    $(addprefix math/,$(MATH_HEADERS))\
    $(addprefix util/,$(UTIL_HEADERS))\
    $(addprefix isct/,$(ISCT_HEADERS))\
    $(addprefix geometry/,$(GEOM_HEADERS))\
    $(addprefix accel/,$(ACCEL_HEADERS))\
    $(addprefix render/,$(RENDER_HEADERS))\
    $(addprefix file_formats/,$(FILE_HEADERS))
HEADER_COPIES     := $(addprefix include/,$(HEADERS))

# +-----------------------------+
# | stand alone program sources |
# +-----------------------------+
MESHGEN_SRC := \
    $(addprefix math/,$(MATH_SRCS))\
    $(addprefix util/,$(UTIL_SRCS))\
    $(addprefix geometry/,$(GEOM_SRCS))\
    $(addprefix render/,$(RENDER_SRCS))\
    $(addprefix file_formats/,$(FILE_SRCS))\
    meshgen

# +---------------------------------------+
# | all sources for dependency generation |
# +---------------------------------------+
ALL_SRCS     := \
    $(SRCS) meshgen
DEPENDS := $(addprefix depend/,$(addsuffix .d,$(ALL_SRCS)))

# +--------------------------------+
# | Object Aggregation for Targets |
# +--------------------------------+

OBJ               := $(addprefix obj/,$(addsuffix .o,$(SRCS))) \
                     obj/isct/triangle.o
DEBUG             := $(addprefix debug/,$(addsuffix .o,$(SRCS))) \
                     obj/isct/triangle.o

#MESHCHECK_OBJ     := $(addprefix obj/,$(addsuffix .o,$(MESHCHECK_SRC)))
#MESHCHECK_DEBUG   := $(addprefix debug/,$(addsuffix .o,$(MESHCHECK_SRC)))

MESHGEN_OBJ       := $(addprefix obj/,$(addsuffix .o,$(MESHGEN_SRC)))
MESHGEN_DEBUG     := $(addprefix debug/,$(addsuffix .o,$(MESHGEN_SRC)))

LIB_TARGET_NAME   := toptop

# *********
# * RULES *
# *********------+
# | Target Rules |
# +--------------+
all: lib/lib$(LIB_TARGET_NAME).a includes \
     bin/meshgen bin/off2obj
debug: lib/lib$(LIB_TARGET_NAME)debug.a includes

lib/lib$(LIB_TARGET_NAME).a: $(OBJ)
	@echo "Bundling $@"
	@ar rcs $@ $(OBJ)

lib/lib$(LIB_TARGET_NAME)debug.a: $(DEBUG)
	@echo "Bundling $@"
	@ar rcs $@ $(DEBUG)

#meshcheck: $(MESHCHECK_OBJ)
#	@echo "Linking meshcheck"
#	@$(CXX) -o meshcheck $(MESHCHECK_OBJ) $(LINK)

bin/meshgen: $(MESHGEN_OBJ)
	@echo "Linking meshgen"
	@$(CXX) -o bin/meshgen $(MESHGEN_OBJ) $(LINK)

bin/filters: obj/filterBuilder.o
	@echo "Linking filters"
	@$(CXX) -o bin/filters obj/filterBuilder.o $(LINK)

bin/off2obj: obj/off2obj.o
	@echo "Linking off2obj"
	@$(CXX) -o bin/off2obj obj/off2obj.o $(LINK)

assembly: obj/isct/empty3d.s

# +------------------------------+
# | Specialized File Build Rules |
# +------------------------------+

obj/isct/triangle.o: src/isct/triangle.c
	@echo "Compiling the Triangle library"
	@$(CC) -O2 -DNO_TIMER \
               -DREDUCED \
               -DCDT_ONLY -DTRILIBRARY \
               -Wall -DANSI_DECLARATORS \
               -o obj/isct/triangle.o -c src/isct/triangle.c

# generate some assembly for hand inspection
obj/isct/empty3d.s: src/isct/empty3d.cpp
	@echo "Compiling to Readable Assembly $@"
	@$(CXX) $(CXXFLAGS) -S -o obj/isct/empty3d.s -c src/isct/empty3d.cpp

# +------------------------------------+
# | Generic Source->Object Build Rules |
# +------------------------------------+
obj/%.o: src/%.cpp
	@echo "Compiling $@"
	@$(CXX) $(CXXFLAGS) -o $@ -c $<

debug/%.o: src/%.cpp
	@echo "Compiling $@"
	@$(CXX) $(CXXDFLAGS) -o $@ -c $<

# dependency file build rules
depend/%.d: src/%.cpp
	@$(CXX) $(CXXFLAGS) -MM $< | \
        sed -e 's@^\(.*\)\.o:@depend/$*.d debug/$*.o obj/$*.o:@' > $@

# +-------------------+
# | include copy rule |
# +-------------------+---------------------
# | This rule exists to safely propagate
# | header file dependencies to other
# | targets that depend on the common code
# +-----------------------------------------
includes: $(HEADER_COPIES)

include/%.h: src/%.h
	@echo "updating $@"
	@cp $< $@
#also support implementation files
include/%.tpp: src/%.tpp
	@echo "updating $@"
	@cp $< $@

# +---------------+
# | cleaning rule |
# +---------------+
# using /*/* to allow two-deep anonymous pattern matching
clean:
	-@$(RM) -r obj depend debug include bin lib
	-@$(RM) bin/meshgen bin/off2obj
#	-@$(RM) gmon.out
	-@$(RM) lib/lib$(LIB_TARGET_NAME).a
	-@$(RM) lib/lib$(LIB_TARGET_NAME)debug.a

-include $(DEPENDS)
