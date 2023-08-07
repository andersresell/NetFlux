
BUILD_DIR_DEBUG = ../build_debug
BUILD_DIR_RELEASE = ../build_release
EXE = NetFlux
EXE_DEBUG = $(BUILD_DIR_DEBUG)/$(EXE)
EXE_RELEASE = $(BUILD_DIR_RELEASE)/$(EXE)

CXX = mpic++ #MPI compiler wrapper (automatically includes headers and libraries)
INC_DIR = ../include

CXX_FLAGS_ALL = -Wall -Wno-unused-function -std=c++20 -fmax-errors=1 #-I$(INC_DIR)
CXX_FLAGS_DEBUG = $(CXX_FLAGS_ALL) -g -D_GLIBCXX_CONCEPT_CHECKS -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
CXX_FLAGS_RELEASE = $(CXX_FLAGS_ALL) -DNDEBUG -O3 -flto 
 
#external libraries
LDLIBS = -lyaml-cpp -lmetis

#add source files here
SRC = NetFlux.cpp \
	  Config.cpp \
	  ConfigParser.cpp \
	  geometry/PrimalGrid.cpp \
	  geometry/FV_Grid.cpp \
	  geometry/GeometryUtilities.cpp \
	  Solver.cpp \
	  Utilities.cpp \
	  Output.cpp \
	  Driver.cpp \
	  SolverData.cpp \
	  Numerics.cpp \
	  parallelization/GridCreator.cpp \
	  parallelization/MPI_Wrapper.cpp \
	  parallelization/MetisWrapper.cpp \
	  parallelization/Communicator.cpp \

OBJ_DEBUG = $(patsubst %.cpp, $(BUILD_DIR_DEBUG)/%.o, $(SRC)) 
OBJ_RELEASE = $(patsubst %.cpp, $(BUILD_DIR_RELEASE)/%.o, $(SRC)) 


all: prep debug 

debug: $(EXE_DEBUG) 

release: $(EXE_RELEASE) 


$(EXE_DEBUG): $(OBJ_DEBUG)
	$(CXX) -o $(EXE_DEBUG) $^ $(LDLIBS) $(CXX_FLAGS_DEBUG)

$(EXE_RELEASE): $(OBJ_RELEASE)
	$(CXX) -o $(EXE_RELEASE) $^ $(LDLIBS) $(CXX_FLAGS_RELEASE)


#The -MMD option creates a list of depenendencies, and it's included later
$(BUILD_DIR_DEBUG)/%.o: %.cpp 
	$(CXX) -MMD -c -o $@ $< $(CXX_FLAGS_DEBUG)

$(BUILD_DIR_RELEASE)/%.o: %.cpp 
	$(CXX) -MMD -c -o $@ $< $(CXX_FLAGS_RELEASE)

-include $(OBJ_DEBUG:.o=.d)
-include $(OBJ_RELEASE:.o=.d)

prep:
	mkdir -p $(BUILD_DIR_DEBUG) 
	mkdir -p $(BUILD_DIR_RELEASE)
	mkdir -p $(BUILD_DIR_DEBUG)/containers
	mkdir -p $(BUILD_DIR_RELEASE)/containers
	mkdir -p $(BUILD_DIR_DEBUG)/geometry
	mkdir -p $(BUILD_DIR_RELEASE)/geometry
	mkdir -p $(BUILD_DIR_DEBUG)/parallelization
	mkdir -p $(BUILD_DIR_RELEASE)/parallelization

.PHONY: clean
clean:
	-rm -r $(BUILD_DIR_DEBUG)/*
	-rm -r $(BUILD_DIR_RELEASE)/*













