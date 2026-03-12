CXX := clang++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pedantic

Q1A := q1a
Q1B := q1b
Q1C := q1c
Q2 := q2

.PHONY: all run-q1a run-q1b run-q1c run-q2 clean

all: $(Q1A) $(Q1B) $(Q1C) $(Q2)

$(Q1A): src/q1a.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

$(Q1B): src/q1b.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

$(Q1C): src/q1c.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

$(Q2): src/question2.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

run-q1a: $(Q1A)
	./$(Q1A)

run-q1b: $(Q1B)
	./$(Q1B)

run-q1c: $(Q1C)
	./$(Q1C)

run-q2: $(Q2)
	./$(Q2)

clean:
	rm -f $(Q1A) $(Q1B) $(Q1C) $(Q2) question1 question2 assignment04
