CXX = g++
CXXFLAGS = -Wall -std=c++11
TARGET = main
SRC = main.cpp

all: $(TARGET)
	./$(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

clean:
	rm -f $(TARGET)
