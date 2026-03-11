CXX := clang++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pedantic
TARGET := assignment04
SRC := src/assignment04.cpp

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
