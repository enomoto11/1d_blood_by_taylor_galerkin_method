CXX = g++
CXXFLAGS = -Wall -std=c++11
TARGET = main
SRC = main.cpp
GORUN = go run

all: go_run compile_and_run $(TARGET) python_scripts

go_run:
	cd input && $(GORUN) main.go

compile_and_run: $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)
	./$(TARGET)

python_scripts:
	python3 A.py
	python3 Q.py
	python3 U.py

clean:
	rm -f $(TARGET)
