CXX = g++
CXXFLAGS = -Wall -std=c++11
TARGET = main
SRC = main.cpp
GORUN = go run

all: go_run $(TARGET)
	./$(TARGET)

go_run:
	cd input && $(GORUN) main.go

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)
