# Compiler
CXX = g++-14

# Compiler flags
CXXFLAGS = -O3 -W -Wall -std=c++17 -fopenmp
LDFLAGS = -fopenmp

# Source files
SRC_DIR = src
SRC_FILES = $(SRC_DIR)/main.cpp $(SRC_DIR)/gauss_legendre.cpp
HEADER_FILES = $(wildcard $(SRC_DIR)/*.hpp)

# Object files
OBJ_FILES = $(SRC_FILES:.cpp=.o)

# Executable name
EXEC_NAME = NN-cms.x

# Build rule
$(EXEC_NAME): $(OBJ_FILES)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Compile rule
%.o: %.cpp $(HEADER_FILES)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean rule
clean:
	rm -f $(EXEC_NAME) $(OBJ_FILES)
