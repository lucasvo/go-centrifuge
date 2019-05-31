package snarks

import (
	//"os"
	"fmt"
	"github.com/stretchr/testify/assert"
	"testing"
	python3 "github.com/DataDog/go-python3"
)

func TestPyCryptoCall(t *testing.T) {
	i, err := python3.Py_Main([]string{"pedersen.py 01238da83"});
	fmt.Println(i, err);
	//i, err = python3.PyRun_AnyFile("/Users/lucasvo/go/src/github.com/centrifuge/go-centrifuge/snarks/wrapper.py");
	//fmt.Println(i, err);
	assert.Nil(t, "not");
}
