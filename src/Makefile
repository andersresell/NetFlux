
CXX = g++

INC_DIR = ../include
BUILD_DIR = ../build

CXX_FLAGS = -g -Wall -I$(INC_DIR)

#external libraries
LDLIBS = -lyaml-cpp


_EXE = C3D
EXE = $(BUILD_DIR)/$(_EXE)

#add source files here
SRC = C3D.cpp \
	  Config.cpp \
	  Grid.cpp \
	  Solver.cpp \
	  Utilities.cpp \
	  Output.cpp \
	  Driver.cpp \
	  SolverData.cpp \
	  Numerics .cpp



INCLUDE = $(patsubst %, $(INC_DIR)/%, $(_INCLUDE))
OBJ = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SRC)) #names of src files put in build dir with .o extension

#target for executable
$(_EXE): $(OBJ)
	$(CXX) -o $(EXE) $^ $(LDLIBS) $(CXX_FLAGS)


#target for object files
#The -MMD option creates a list of depenendencies, and it's included later
$(BUILD_DIR)/%.o: %.cpp 
	$(CXX) -MMD -c -o $@ $< $(CXX_FLAGS)


-include $(OBJ:.o=.d)

.PHONY: clean
clean:
	-rm -f $(BUILD_DIR)/*