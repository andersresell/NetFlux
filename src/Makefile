
CXX = g++

INC_DIR = ../include
BUILD_DIR = ../build

CXX_FLAGS = -g -Wall -I$(INC_DIR)

_EXE = C3D
EXE = $(BUILD_DIR)/$(_EXE)

#add source files here
SRC = C3D.cpp \
	  Config.cpp \
	  Grid.cpp \
	  Solver.cpp \
	  Utilities.cpp

#add header files here
_INCLUDE = includes.hpp \
		   Config.hpp \
           Grid.hpp \
		   Solver.hpp \
		   Utilities.hpp

INCLUDE = $(patsubst %, $(INC_DIR)/%, $(_INCLUDE))
OBJ = $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SRC)) #names of src files put in build dir with .o extension

#target for executable
$(_EXE): $(OBJ)
	$(CXX) -o $(EXE) $^ $(CXX_FLAGS)


#target for object files
$(BUILD_DIR)/%.o: %.cpp $(INCLUDE)
	$(CXX) -c -o $@ $< $(CXX_FLAGS)




.PHONY: clean
clean:
	-rm -f $(BUILD_DIR)/*