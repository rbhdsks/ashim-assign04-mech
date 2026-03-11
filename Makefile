CXX := clang++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pedantic

Q1 := question1
Q2 := question2

.PHONY: all run-q1 run-q2 clean

all: $(Q1) $(Q2)

$(Q1): src/question1.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

$(Q2): src/question2.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

run-q1: $(Q1)
	./$(Q1)

run-q2: $(Q2)
	./$(Q2)

clean:
	rm -f $(Q1) $(Q2) assignment04
