package main

import (
	"fmt"
	"os"
)

const (
	ElementD = 100 // ElementD is the number of elements in the D matrix
)

func outputElement_d() {
	file, err := os.Create("./element_d.dat")
	if err != nil {
		fmt.Println("Error creating file:", err)
		return
	}
	defer file.Close()

	for i := 0; i < ElementD; i++ {
		if i == ElementD-1 {
			fmt.Fprintf(file, "%d %d", i, i+1)
		} else {
			fmt.Fprintf(file, "%d %d\n", i, i+1)
		}
	}
}

func outputNode_d() {
	file, err := os.Create("./node_d.dat")
	if err != nil {
		fmt.Println("Error creating file:", err)
		return
	}
	defer file.Close()

	for i := 0; i <= ElementD; i++ {
		var x float32
		x = float32(i) / float32(ElementD)

		if i == ElementD {
			fmt.Fprintf(file, "%.2f %d %d", x, 0, 0)
		} else {
			fmt.Fprintf(file, "%.2f %d %d\n", x, 0, 0)
		}
	}
}

func main() {
	outputElement_d()
	outputNode_d()
}
