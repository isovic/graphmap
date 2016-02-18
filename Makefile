BIN = ./bin/graphmap-not_release
BIN_DEBUG = ./bin/graphmap-debug
BIN_LINUX = ./bin/Linux-x64/graphmap
BIN_MAC = ./bin/Mac/graphmap
OBJ_TESTING = ./obj_test
OBJ_TESTING_EXT = ./obj_testext
OBJ_DEBUG = ./obj_debug
OBJ_LINUX = ./obj_linux
OBJ_EXTCIGAR = ./obj_extcigar
OBJ_MAC = ./obj_mac
SOURCE = src
CODEBASE = codebase
# This finds all 'src' folders at maximum depth 2 (level one inside each submodule's folder).
CODEBASE_SRC_FOLDERS = $(shell find $(CODEBASE) -maxdepth 2 -type d -name "src" -exec echo "-I"{} \;)
# $(shell find $(CODEBASE) -maxdepth 3 -type d -name "libs" -exec echo "-I"{} \;)
# $(shell find $(CODEBASE) -maxdepth 2 -type d -name "src" -exec echo "-I"{}"/*/" \;)

# ? allows override by user using env var
GCC ?= g++
# define variables for GCC version check here
GCC_MAJOR_VERSION_GE_4 := $(shell expr `$(GCC) -dumpversion | cut -f1 -d.` \>= 4)
GCC_MINOR_VERSION_GE_7 := $(shell expr `$(GCC) -dumpversion | cut -f2 -d.` \>= 7)
GCC_MAC ?= /opt/local/bin/g++-mp-4.8


# CPP_FILES := $(wildcard $(SOURCE)/*/*.cpp) $(wildcard $(SOURCE)/*.cpp) $(wildcard $(SOURCE)/libs/*/*.cpp)
# CC_FILES := $(wildcard $(SOURCE)/*/*.cc) $(wildcard $(SOURCE)/*.cc) $(wildcard $(SOURCE)/libs/*/*.cc)
# H_FILES := $(wildcard $(SOURCE)/*/*.h) $(wildcard $(SOURCE)/*.h) $(wildcard $(SOURCE)/libs/*/*.h)
CPP_FILES :=  $(wildcard $(CODEBASE)/*/src/*.cpp) $(wildcard $(CODEBASE)/*/src/libs/*/*.cpp) $(wildcard $(CODEBASE)/*/src/*/*.cpp) $(wildcard $(SOURCE)/*/*.cpp) $(wildcard $(SOURCE)/*.cpp) $(wildcard $(SOURCE)/libs/*/*.cpp)
CC_FILES :=  $(wildcard $(CODEBASE)/*/src/*.cc) $(wildcard $(CODEBASE)/*/src/libs/*/*.cc) $(wildcard $(CODEBASE)/*/src/*/*.cc) $(wildcard $(SOURCE)/*/*.cc) $(wildcard $(SOURCE)/*.cc) $(wildcard $(SOURCE)/libs/*/*.cc)
H_FILES := $(wildcard $(CODEBASE)/*/src/*.h) $(wildcard $(CODEBASE)/*/src/libs/*/*.h) $(wildcard $(CODEBASE)/*/src/*/*.h) $(wildcard $(SOURCE)/*/*.h) $(wildcard $(SOURCE)/*.h) $(wildcard $(CODEBASE)/*/src/*.hpp) $(wildcard $(CODEBASE)/*/src/*/*.hpp) $(wildcard $(SOURCE)/*/*.hpp) $(wildcard $(SOURCE)/*.hpp) $(wildcard $(SOURCE)/libs/*/*.h)

OBJ_FILES := $(CPP_FILES:.cpp=.o) $(CC_FILES:.cc=.o)
OBJ_FILES_FOLDER_TESTING := $(addprefix $(OBJ_TESTING)/,$(OBJ_FILES))
OBJ_FILES_FOLDER_TESTING_EXT := $(addprefix $(OBJ_TESTING_EXT)/,$(OBJ_FILES))
OBJ_FILES_FOLDER_DEBUG := $(addprefix $(OBJ_DEBUG)/,$(OBJ_FILES))
OBJ_FILES_FOLDER_LINUX := $(addprefix $(OBJ_LINUX)/,$(OBJ_FILES))
OBJ_FILES_FOLDER_EXTCIGAR := $(addprefix $(OBJ_EXTCIGAR)/,$(OBJ_FILES))
OBJ_FILES_FOLDER_MAC := $(addprefix $(OBJ_MAC)/,$(OBJ_FILES))

LIB_DIRS = -L"/usr/local/lib" -L"$(CODEBASE)/seqlib/src/libs/libdivsufsort-2.0.1/build/lib"
CC_LIBS = -static-libgcc -static-libstdc++ -D__cplusplus=201103L
# INCLUDE = -I"./src/" -I"/usr/include/" -I"libs/libdivsufsort-2.0.1/build/include" -I"libs/seqan-library-1.4.2/include"
# INCLUDE = -I"./src/" -I"/usr/include/" -I"src/libs/seqan-library-1.4.2/include"
INCLUDE = -I"./src/" -I"/usr/include/" -I"$(CODEBASE)/seqlib/src/libs/seqan-library-2.0.1/include" -I"$(CODEBASE)/seqlib/src/libs/libdivsufsort-2.0.1-64bit/" $(CODEBASE_SRC_FOLDERS)

CC_FLAGS_DEBUG = -O0 -g -rdynamic -c -fmessage-length=0 -ffreestanding -fopenmp -m64 -std=c++11 -Werror=return-type -pthread -march=native
CC_FLAGS_RELEASE = -DRELEASE_VERSION -O3 -fdata-sections -ffunction-sections -c -fmessage-length=0 -ffreestanding -fopenmp -m64 -std=c++11 -Werror=return-type -pthread # -march=native
CC_FLAGS_EXTCIGAR = -DRELEASE_VERSION -DUSE_EXTENDED_CIGAR_FORMAT -O3 -fdata-sections -ffunction-sections -c -fmessage-length=0 -ffreestanding -fopenmp -m64 -std=c++11 -Werror=return-type -pthread -march=native
CC_FLAGS_NOT_RELEASE = -O3 -fdata-sections -ffunction-sections -c -fmessage-length=0 -ffreestanding -fopenmp -m64 -std=c++11 -Werror=return-type -Wuninitialized -pthread -march=native
CC_FLAGS_NOT_RELEASE_EXT = -O3 -DUSE_EXTENDED_CIGAR_FORMAT -fdata-sections -ffunction-sections -c -fmessage-length=0 -ffreestanding -fopenmp -m64 -std=c++11 -Werror=return-type -Wuninitialized -pthread -march=native
LD_FLAGS = -static-libgcc -static-libstdc++ -m64 -ffreestanding
# LD_LIBS = -lpthread -lgomp -lm -lz -ldivsufsort64
LD_LIBS = -lpthread -lgomp -lm -lz



all: gcc_version_check linux



modules:
#	git submodule update --init --recursive
	git submodule foreach git pull origin master

testing: modules $(OBJ_FILES_FOLDER_TESTING)
	mkdir -p $(dir $(BIN))
	$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN) $(OBJ_FILES_FOLDER_TESTING) $(LD_LIBS)
	
obj_test/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_NOT_RELEASE) -o $@ $<
	
obj_test/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_NOT_RELEASE) -o $@ $<

testingext: modules $(OBJ_FILES_FOLDER_TESTING_EXT)
	mkdir -p $(dir $(BIN))
	$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN) $(OBJ_FILES_FOLDER_TESTING_EXT) $(LD_LIBS)
	
obj_testext/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_NOT_RELEASE_EXT) -o $@ $<
	
obj_testext/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_NOT_RELEASE_EXT) -o $@ $<



gcc_version_check:
ifneq ($(GCC_MAJOR_VERSION_GE_4), 1)
	$(warning "*** WARNING $(GCC) major version <4 ***")
endif	
ifneq ($(GCC_MINOR_VERSION_GE_7), 1)
	$(warning "*** WARNING $(GCC) minor version <7 ***")
endif


debug: modules $(OBJ_FILES_FOLDER_DEBUG)
	mkdir -p $(dir $(BIN_DEBUG))
	$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN_DEBUG) $(OBJ_FILES_FOLDER_DEBUG) $(LD_LIBS)
	
obj_debug/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_DEBUG) -o $@ $<
	
obj_debug/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_DEBUG) -o $@ $<



linux: modules $(OBJ_FILES_FOLDER_LINUX)
	mkdir -p $(dir $(BIN_LINUX))
	$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN_LINUX) $(OBJ_FILES_FOLDER_LINUX) $(LD_LIBS)
	
obj_linux/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<
	
obj_linux/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<



extcigar: modules $(OBJ_FILES_FOLDER_EXTCIGAR)
	mkdir -p $(dir $(BIN_LINUX))
	$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN_LINUX) $(OBJ_FILES_FOLDER_EXTCIGAR) $(LD_LIBS)
	
obj_extcigar/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_EXTCIGAR) -o $@ $<
	
obj_extcigar/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_EXTCIGAR) -o $@ $<



mac: modules $(OBJ_FILES_FOLDER_MAC)
	mkdir -p $(dir $(BIN_MAC))
	$(GCC_MAC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN_MAC) $(OBJ_FILES_FOLDER_MAC) $(LD_LIBS)
	
obj_mac/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC_MAC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<
	
obj_mac/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC_MAC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<



# deps:
# 	cd libs; cd libdivsufsort-2.0.1; make clean; rm -rf build; ./configure; mkdir build ;cd build; cmake -DBUILD_DIVSUFSORT64:BOOL=ON -DCMAKE_BUILD_TYPE="Release" -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="/usr/local" .. ; make


	
clean:
	-rm -rf $(OBJ_LINUX) $(BIN_LINUX)

cleantesting:
	-rm -rf $(OBJ_TESTING) $(BIN)

cleandebug:
	-rm -rf $(OBJ_DEBUG) $(BIN_DEBUG)

cleanlinux:
	-rm -rf $(OBJ_LINUX) $(BIN_LINUX)

cleanextcigar:
	-rm -rf $(OBJ_EXTCIGAR) $(BIN_LINUX)

cleanmac:
	-rm -rf $(OBJ_MAC) $(BIN_MAC)

cleanbin:
	-rm -rf bin/

cleanall: clean cleantest cleandebug cleanmac cleanbin



rebuild: clean all

rebuilddebug: cleandebug debug

rebuildlinux: cleanlinux linux

rebuildtesting: cleantesting testing

rebuildmac: cleanmac mac

# divsufsort:
# 	cd libs; ./build-libdivsufsort.sh

