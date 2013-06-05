# Warning: This makefile rampantly assumes it is being run in UNIX
#   Also, there may be some gnu tool chain assumptions: gcc and gmake?

# ***********
# * PREFACE *
# ***********+
# | Platform |
# +----------+-------------------
# | evaluates to 'Darwin' on mac
# | evaluates to ??? on ???
# +------------------------------
PLATFORM := $(shell uname)

-include makeConstants

SHELL_HACK := $(shell mkdir -p bin)

# *********
# * TOOLS *
# *********--+
# | Programs |
# +----------+
CC  := gcc
CXX := g++
RM  := rm

# also let the preprocessor know which platform we're on
ifeq ($(PLATFORM),Darwin)
  CONFIG   := $(CONFIG) -DMACOSX
else
  CONFIG   := $(CONFIG)
endif

#MAKE_ANN   :=
#ifeq ($(PLATFORM),Darwin)
#  MAKE_ANN := make -C ann_1.1.2 macosx-g++
#else
#  MAKE_ANN := make -C ann_1.1.2 linux-g++
#endif


COMMON_TARGETS        := 
QT_TARGETS            := edit modeler
# *********
# * RULES *
# *********------+
# | Target Rules |
# +--------------+
all: $(COMMON_TARGETS) $(QT_TARGETS)

$(QT_TARGETS): build_common
	@echo "building target: $@"
	@cd qt; \
	qmake -config $@; \
	make
	@cp qt/$@ bin/$@

build_common:
	@echo "building common library"
	@make -C common

build_glut: build_common
	@make -C glut


distribute:
	@echo "packaging distribution"
	@mkdir toptop
	@cp -r Eigen toptop/
	@mkdir -p toptop/common
	@cp -r common/src toptop/common/src
	@cp -r common/shaders toptop/common/shaders
	@cp common/Makefile toptop/common/Makefile
	@mkdir -p toptop/qt
	@cp -r qt/src toptop/qt/src
	@cp -r qt/everything.pro toptop/qt/everything.pro
	@mkdir -p toptop/resources
	@cp -r resources/demos toptop/resources/demos
	@cp -r resources/models toptop/resources/models
	@cp COPYRIGHT toptop/
	@cp Makefile toptop/
	@cp README toptop/
	@cp makeConstants toptop/
	@tar -cvzf toptop.tgz toptop
	@rm -r toptop

#build_ann:
#	@echo "building ANN library"
#	@$(MAKE_ANN)

# +---------------+
# | cleaning rule |
# +---------------+
clean:
	make -C qt clean
	-@$(RM) -r qt/moc qt/obj
	make -C common clean
#	-make -C ann_1.1.2 clean
	-@$(RM) -r bin
#	-@$(RM) gmon.out





