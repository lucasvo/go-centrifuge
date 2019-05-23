package snarks

import (
	python "github.com/sbinet/go-python"
	"github.com/stretchr/testify/assert"
	"testing"
)

func TestPyCryptoCall(t *testing.T) {
	python.Initialize()
	defer python.Finalize()

	wrapperModule := python.PyImport_ImportModule("wrapper")
	if wrapperModule == nil {
		panic("Error importing module")
	}

	wrapperFunc := fooModule.GetAttrString("hello")
	if wrapperFunc == nil {
		panic("Error importing function")
	}

	// The Python function takes no params but when using the C api
	// we're required to send (empty) *args and **kwargs anyways.
	helloFunc.Call(python.PyTuple_New(0), python.PyDict_New())
}
