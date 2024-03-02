# Compiler
CXX = mpic++ -O3

# Compiler flags
CXXFLAGS = -std=c++17 -fopenmp
LDFLAGS = -fopenmp

# Source files
SRC_DIR = src
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
HEADER_FILES = $(wildcard $(SRC_DIR)/*.hpp)

# Object files
OBJ_FILES = $(SRC_FILES:.cpp=.o)

# Executable name
EXEC_NAME = NNcms.x

# Build rule
$(EXEC_NAME): $(OBJ_FILES)
	$(CXX) $(LDFLAGS) -o $@ $^

# Compile rule
%.o: %.cpp $(HEADER_FILES)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean rule
clean:
	rm -f $(EXEC_NAME) $(OBJ_FILES)

