

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -g

# Source file
SRC = main.cpp

# Output executable
TARGET = main

# Default target
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Clean rule
clean:
	rm -f $(TARGET)

.PHONY: all clean