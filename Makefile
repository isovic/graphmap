BIN = ./bin/graphmap-not_release
BIN_DEBUG = ./bin/graphmap-debug
BIN_LINUX = ./bin/Linux-x64/graphmap
BIN_MAC = ./bin/Mac/graphmap
OBJ = ./obj_test
OBJ_DEBUG = ./obj_debug
OBJ_LINUX = ./obj_linux
OBJ_MAC = ./obj_mac
SOURCE = src

GCC = g++
GCC_MAC = /opt/local/bin/g++-mp-4.8

CPP_FILES := $(wildcard $(SOURCE)/*/*.cpp) $(wildcard $(SOURCE)/*.cpp)
CC_FILES := $(wildcard $(SOURCE)/*/*.cc) $(wildcard $(SOURCE)/*.cc)
H_FILES := $(wildcard $(SOURCE)/*/*.h) $(wildcard $(SOURCE)/*.h)

OBJ_FILES := $(CPP_FILES:.cpp=.o) $(CC_FILES:.cc=.o)
OBJ_FILES_FOLDER := $(addprefix $(OBJ)/,$(OBJ_FILES))
OBJ_FILES_FOLDER_DEBUG := $(addprefix $(OBJ_DEBUG)/,$(OBJ_FILES))
OBJ_FILES_FOLDER_LINUX := $(addprefix $(OBJ_LINUX)/,$(OBJ_FILES))
OBJ_FILES_FOLDER_MAC := $(addprefix $(OBJ_MAC)/,$(OBJ_FILES))

LIB_DIRS = -L"/usr/local/lib" -L"libs/libdivsufsort-2.0.1/build/lib"
CC_LIBS = -static-libgcc -static-libstdc++ -D__cplusplus=201103L
INCLUDE = -I"./src/" -I"/usr/include/" -I"libs/libdivsufsort-2.0.1/build/include" -I"libs/seqan-library-1.4.2/include"

CC_FLAGS_DEBUG = -O0 -g -rdynamic -c -fmessage-length=0 -ffreestanding -fopenmp -m64 -std=c++11 -Werror=return-type -pthread
CC_FLAGS_RELEASE = -DRELEASE_VERSION -O3 -fdata-sections -ffunction-sections -c -fmessage-length=0 -ffreestanding -fopenmp -m64 -std=c++11 -Werror=return-type -pthread
CC_FLAGS_NOT_RELEASE = -O3 -fdata-sections -ffunction-sections -c -fmessage-length=0 -ffreestanding -fopenmp -m64 -std=c++11 -Werror=return-type -Wuninitialized -pthread
LD_FLAGS = -static-libgcc -static-libstdc++ -m64 -ffreestanding
LD_LIBS = -lpthread -lgomp -lm -lz -ldivsufsort64



all: linux



testing: $(OBJ_FILES_FOLDER)
	mkdir -p $(dir $(BIN))
	$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN) $(OBJ_FILES_FOLDER) $(LD_LIBS)
	
obj_test/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_NOT_RELEASE) -o $@ $<
	
obj_test/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_NOT_RELEASE) -o $@ $<



debug: $(OBJ_FILES_FOLDER_DEBUG)
	mkdir -p $(dir $(BIN_DEBUG))
	$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN_DEBUG) $(OBJ_FILES_FOLDER_DEBUG) $(LD_LIBS)
	
obj_debug/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_DEBUG) -o $@ $<
	
obj_debug/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_DEBUG) -o $@ $<



linux: $(OBJ_FILES_FOLDER_LINUX)
	mkdir -p $(dir $(BIN_LINUX))
	$(GCC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN_LINUX) $(OBJ_FILES_FOLDER_LINUX) $(LD_LIBS)
	
obj_linux/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<
	
obj_linux/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<



mac: $(OBJ_FILES_FOLDER_MAC)
	mkdir -p $(dir $(BIN_MAC))
	$(GCC_MAC) $(LD_FLAGS) $(LIB_DIRS) -o $(BIN_MAC) $(OBJ_FILES_FOLDER_MAC) $(LD_LIBS)
	
obj_mac/%.o: %.cc $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC_MAC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<
	
obj_mac/%.o: %.cpp $(H_FILES)
	mkdir -p $(dir $@)
	$(GCC_MAC) $(CC_LIBS) $(INCLUDE) $(CC_FLAGS_RELEASE) -o $@ $<



deps:
	cd libs; cd libdivsufsort-2.0.1; make clean; rm -rf build; ./configure; mkdir build ;cd build; cmake -DBUILD_DIVSUFSORT64:BOOL=ON -DCMAKE_BUILD_TYPE="Release" -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX="/usr/local" .. ; make


	
clean:
	-rm -rf $(OBJ_LINUX) $(BIN_LINUX)

cleantesting:
	-rm -rf $(OBJ) $(BIN)

cleandebug:
	-rm -rf $(OBJ_DEBUG) $(BIN_DEBUG)

cleanlinux:
	-rm -rf $(OBJ_LINUX) $(BIN_LINUX)

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

divsufsort:
	cd libs; ./build-libdivsufsort.sh

