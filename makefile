CXX = g++
FLAGS = -std=c++20 -O3
LIBS = -llapack

SOURCE = $(wildcard source/*.cpp)
S_OBJ = $(patsubst %.cpp, build/%.o, $(SOURCE))

TESTS = $(wildcard tests/*.cpp)
T_OBJ = $(patsubst %.cpp, build/%.o, $(TESTS))

TIMING = $(wildcard experiments/*.cpp)
TM_OBJ = $(patsubst %.cpp, build/%.o, $(TIMING))

INCLUDE = -I ./
TARGET = test

$(TARGET): $(T_OBJ) $(S_OBJ)
	$(CXX) $(FLAGS) -o $(TARGET) $(T_OBJ) $(S_OBJ) $(LIBS)

timing: $(S_OBJ) $(TM_OBJ)
	$(CXX) $(FLAGS) -o timing $(TM_OBJ) $(S_OBJ) $(LIBS)

build/%.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(FLAGS) -c $< -o $@ $(INCLUDE)

clean:
	rm -rf $(TARGET) build