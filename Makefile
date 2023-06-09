CXX = g++
CXXFLAGS = -Wall -std=c++11
TARGET = main
SRC = main.cpp
GORUN = go run

all: go_run compile_and_run $(TARGET) py

go_run:
	cd input && $(GORUN) main.go

compile_and_run: $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)
	./$(TARGET)

py:
	python3 visualize.py

clean:
	rm -f $(TARGET)
