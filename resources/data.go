// Code generated by go-bindata.
// sources:
// ../build/configs/default_config.yaml
// ../build/configs/testing_config.yaml
// DO NOT EDIT!

package resources

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"path/filepath"
	"strings"
	"time"
)

func bindataRead(data []byte, name string) ([]byte, error) {
	gz, err := gzip.NewReader(bytes.NewBuffer(data))
	if err != nil {
		return nil, fmt.Errorf("Read %q: %v", name, err)
	}

	var buf bytes.Buffer
	_, err = io.Copy(&buf, gz)
	clErr := gz.Close()

	if err != nil {
		return nil, fmt.Errorf("Read %q: %v", name, err)
	}
	if clErr != nil {
		return nil, err
	}

	return buf.Bytes(), nil
}

type asset struct {
	bytes []byte
	info  os.FileInfo
}

type bindataFileInfo struct {
	name    string
	size    int64
	mode    os.FileMode
	modTime time.Time
}

func (fi bindataFileInfo) Name() string {
	return fi.name
}
func (fi bindataFileInfo) Size() int64 {
	return fi.size
}
func (fi bindataFileInfo) Mode() os.FileMode {
	return fi.mode
}
func (fi bindataFileInfo) ModTime() time.Time {
	return fi.modTime
}
func (fi bindataFileInfo) IsDir() bool {
	return false
}
func (fi bindataFileInfo) Sys() interface{} {
	return nil
}

var _goCentrifugeBuildConfigsDefault_configYaml = []byte("\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff\xe4\x57\xc9\x72\xdc\x38\x12\xbd\xf3\x2b\x32\xca\x97\x99\x83\x25\x70\x03\xc9\xba\x95\x16\xf7\x62\x59\x23\x95\x64\xab\xa5\x1b\x08\x24\x49\xb4\x58\x00\x0d\x80\xb5\xf4\xd7\x4f\x80\x8b\x2c\xd9\x96\x3a\x62\x26\xe6\x34\x75\x62\x10\x99\x2f\xb7\xf7\xc0\xac\x77\x70\x86\x15\xeb\x5b\x07\x02\xb7\xd8\xea\x6e\x83\xca\x81\x43\xeb\x14\x3a\x60\x35\x93\xca\x3a\x30\x52\x3d\x62\x79\x08\x38\x2a\x67\x64\xd5\xd7\x78\x89\x6e\xa7\xcd\xe3\x12\x4c\x6f\xad\x64\xaa\x91\x6d\x1b\x0c\x60\x52\x21\xb8\x06\x41\x4c\xb8\x6a\xb4\xb4\xe0\x1a\xe6\xe0\xf4\x09\x01\x36\x4c\x2a\xe7\xf1\x83\xd9\x64\x19\x00\xbc\x83\x0b\xcd\x59\x3b\xa4\x20\x55\x0d\x5c\x2b\x67\x18\x77\xc0\x84\x30\x68\x2d\x5a\x50\x88\x02\x9c\x86\x12\xc1\xa2\x83\x9d\x74\x0d\xa0\xda\xc2\x96\x19\xc9\xca\x16\xed\x51\x00\xb3\xbf\x87\x04\x90\x62\x09\x71\x1c\x0f\xcf\xe8\x1a\x34\xd8\x6f\xa6\x0a\x7e\x13\x4b\xc8\xe3\x7c\x3c\x2b\xb5\x76\xd6\x19\xd6\x5d\x21\x1a\x3b\xfa\xbe\x87\xc5\xb1\xec\x92\xe3\x30\xca\x8e\xc8\x11\x39\x0a\x8f\x1d\xef\x8e\xe3\x3c\x22\xd1\xb1\xec\x2a\x7b\x7c\xbd\xb9\xbd\xde\x97\xbb\xc7\xfe\xe1\xfe\xfe\xac\xea\xff\xba\x2d\xf7\xe7\xab\x35\xde\x5e\x9e\x5e\xe8\xbf\x0e\x87\x34\xcd\xb7\xd7\xaa\xfe\xb2\xbd\xfa\xf4\xe7\xc5\xfd\xe3\xe2\x6f\x40\xe3\x19\xf4\x4b\x45\xcf\x2f\xe9\xe6\xf1\xeb\x1d\xfe\x79\xf7\xf1\x2e\xfa\x7a\xd5\x87\xf4\x8f\x4e\xfc\x12\x3f\xfe\xae\xc3\xdb\x78\xd3\xb0\xe6\xea\x24\xbd\xc1\x54\x85\x23\xe8\xdc\xaa\xd5\xdc\xa9\xb1\x00\x5f\x3e\x2a\x27\xdd\xe1\x03\xe3\x4e\x9b\xc3\x12\x16\x8b\xef\x4e\xd6\x58\x4b\xeb\x5e\x1c\x31\xc5\x1b\x6d\xd6\xd8\x69\x2b\xbf\xf3\xea\xd8\xc1\xd3\xe4\x5f\x65\x2b\x6b\xe6\xa4\x56\xc3\xd9\x30\xbc\x4f\x4c\xaa\x9f\x52\x69\x9a\x71\x00\xcf\x19\x33\x26\xf8\x0e\x2e\xfb\x0d\x1a\xc9\xe1\xb7\x33\xd0\xd5\xc0\x9e\x67\x3c\xf9\xe6\x39\x0e\x32\x0d\x27\xaf\x93\x79\x5a\xd0\x4a\xeb\xbc\xa7\xd2\x02\x7f\x24\x5a\x67\xf4\x56\x0e\x07\x7a\xc0\x7e\x96\xc0\x9c\xde\xdf\x4e\x3f\x4e\x8f\xa2\x28\x3d\x8a\x08\x39\x4a\xa2\xef\x19\x10\x46\x67\xf1\x47\xad\xef\x2e\xa4\xe4\xd7\x5f\x76\xb7\xcd\xed\xc9\x3d\xdd\x7f\xe4\x57\xfa\xa2\xa2\xeb\xeb\xfb\xdf\x3f\x74\xbb\x2a\x34\x59\xba\xbb\xd8\x47\x0f\xeb\xb8\x3b\x15\xe1\xe2\x67\xf0\x39\x3d\x8a\x42\xf2\x1a\xfc\xf5\xc3\xa7\x55\xfe\xcb\xd5\xaf\x66\x7b\xfe\x70\x52\xec\xc4\xa3\xfe\xcc\x57\xab\xcd\xe9\xc3\xaf\x5d\x81\x87\xc3\x43\x72\x73\x9e\xd7\x1f\x4c\xdc\xdc\x5e\xfe\xb1\x98\x7a\x74\x3e\xb1\x7d\xee\xa2\x6f\xf1\x7b\x58\x4f\x7a\x7e\x45\x0f\xc9\xe4\x7c\xc1\x7c\x7b\x40\x60\xd7\xea\x03\x0a\xb8\xd9\x30\xe3\xe0\x74\xa2\x99\x85\x4a\x9b\xa1\xa1\xb5\xdc\xa2\x7a\xd1\xca\x1f\xa9\x08\xaf\x72\x91\xec\x0b\x22\xa2\x22\x49\xb3\x10\xb3\x38\x4f\x22\x5a\x64\x8c\xd2\x32\x63\x45\xc1\x48\x21\x04\xe5\x59\x2c\xe2\x94\x8a\x37\x58\x4b\xf6\x05\xa5\x84\x93\xb8\x10\x71\x18\x26\x69\xcc\x2a\x22\xd2\x9c\xa7\x94\xd2\x2c\x8a\x45\xc1\xa3\x8a\x65\x82\x22\x7f\x83\xdf\x64\x9f\x55\x79\x9a\x88\x8a\x15\x39\x09\x23\x91\x55\x2c\x4d\x79\x4e\xe2\xb2\x64\x51\x44\x49\xc9\x05\x62\x52\xa6\x28\xde\x52\x02\xd9\x8b\x92\xa4\x79\xb8\x2a\xe2\x28\xa7\x34\xc9\xd3\x34\x8e\xf2\x95\x38\x2b\xc9\x79\x94\x86\x61\x9e\xd0\x84\x54\x05\xa6\x67\x83\x66\x4a\x34\x8a\xb5\x0d\xca\xba\x71\xf6\x3f\x13\x44\xf4\x5f\x0a\xe2\x45\x0a\xff\xa7\x92\xf8\xa8\xb7\x4c\xbd\x2a\x88\xe8\x7f\xa0\x88\x37\x04\x91\xa7\x65\x1c\x55\x19\x8b\xab\x84\x24\x79\x58\x85\x51\x1c\x27\x24\x09\x69\x46\x78\xce\x4b\x24\x59\x95\x89\xac\xe0\x6f\x0a\x22\x4d\x18\xc6\x59\x5c\x91\x82\x56\xac\x8a\x44\x49\xcb\x9c\x25\x34\x0b\x33\x4e\xca\x22\x47\x5e\x31\x92\xa5\x42\xbc\x29\x88\x24\x49\x2a\x9a\x14\x18\x93\x2c\x49\x22\xcc\x28\xe7\x55\x16\x67\x09\xa5\x98\x46\x55\x48\x49\x51\x16\x79\x44\xc9\xdb\x82\x20\x49\x98\x61\x19\x67\x45\x12\x86\x34\x89\x69\x9e\x90\xf0\x8c\x52\x5a\xe4\x09\x3f\x3f\xcb\x68\x91\xac\x4e\xf8\x49\x19\x2e\x82\xe0\x1d\x78\xaa\xbd\x77\xfa\x7d\x87\x68\x7c\xdb\x2a\x59\xf7\x66\xc0\xb2\x41\x17\x75\xe3\x92\x70\x2b\x37\xa8\x7b\x07\xbb\x06\x15\xe8\x0e\xd5\xb4\x2b\x28\xe4\x83\xa5\xa7\xb6\x07\xb0\x01\xcc\xaf\x27\x97\x25\x2c\x62\x62\x87\x48\xd7\x3d\xf6\xf8\x5d\x88\x61\x84\xcc\x1e\x14\x6f\x8c\x56\xba\xb7\x5e\x2d\x1c\xad\x95\xaa\x0e\xbe\x7a\x87\x31\x81\x71\xd3\xb1\xc3\xb4\x55\xbf\x29\xd1\x78\xbd\x79\xc2\xa0\xb1\xc7\x5c\x2b\xeb\x25\x3c\x69\x6f\xe7\x3f\x35\x25\x02\x6b\x5b\xcd\x99\x43\x01\xcc\x81\x75\xcc\xb8\xbe\x0b\xc0\xfb\xdf\x8d\x8e\x4b\x88\x06\xf4\x0f\x06\xd1\x42\xdf\xc1\xe9\xd5\x67\xe0\x07\xde\xa2\x1d\x4b\x1d\x03\x80\xb4\xb0\x63\x72\x58\x90\x7c\xbe\xb8\x45\xe5\x7c\xa9\xe3\xf1\x1d\x93\x43\xb5\x9f\x6e\x96\x10\xfa\x42\x9f\x28\x6f\x3b\xe4\xb2\x92\xfc\x65\xd1\xc1\x4c\xf9\xb1\xb4\x1b\x6c\xd1\x93\x79\xd7\x48\xde\x3c\xc9\x01\x18\xe7\xba\xf7\x5f\x74\x0d\xbd\xc5\xf9\x5e\xd2\xbe\x09\xd3\x85\x22\x40\xaa\xe1\x25\xef\xad\xd3\x9b\x29\x08\x54\xb2\xc5\x00\xe6\x85\x70\x35\xc2\x5c\xb2\x0d\x2e\x61\xe1\x97\xc0\xc5\xd3\xda\xe7\x93\x99\x81\x9f\xe2\xf2\x56\xfa\x45\xc2\x5f\x65\xf0\x8f\x1d\x82\xc1\xaf\xbd\x34\x08\x3b\x0b\xda\x80\xec\xf8\xb4\x0b\xfa\xd5\xcf\x3f\x72\xe6\x7c\xda\x43\x4b\xfe\xe9\xbb\xab\x05\x7e\x5e\x5f\x2c\xa1\x71\xae\x5b\x1e\x1f\xfb\x11\xb4\x8d\xb6\x6e\x59\xa4\x49\x3a\x0f\x73\xd8\x55\x6b\xe6\x6b\x91\xdc\xa7\x5b\x33\x7b\xe5\x1f\x97\x10\x92\xf9\xf7\x83\x71\x2b\x37\xd2\x8d\xc6\x17\xfe\x71\x09\x49\x16\x46\x71\x9e\xbf\x20\xa9\xd3\xc3\xb4\x46\x6a\xa9\x6f\x95\x39\xc3\x94\x65\x03\x61\xe7\x1a\x84\x18\x77\x5b\x06\x65\xab\xf9\x23\x30\x25\xa6\x52\xc0\x19\x59\xd7\x68\x50\x8c\x94\x76\xb8\x77\xf3\xa0\x47\x5a\x53\xe2\x79\xfd\x5a\x60\x83\x4c\x80\x56\xed\xc1\xcb\x65\x26\xfb\xbc\xe0\xcf\x29\x7d\x83\x5e\x23\x13\x2f\xe1\xc3\x74\x42\xbf\xf4\x93\x78\x9e\x7b\xa7\x75\x0b\x1b\xb6\x07\x83\xce\xc8\xf1\xbb\x62\x51\x09\x60\x2f\xcc\xf4\x16\x4d\x00\xde\x70\x3d\xda\x2d\x21\x9a\x7a\xfa\x73\x48\xa9\x1c\x9a\x2d\x6b\x07\xdc\xc3\x28\x00\xe6\x13\xe4\xbd\x31\xc3\x72\xf9\xcc\xa3\x61\x16\x4a\x44\xbf\x7d\x3a\xe4\x6e\x68\xd3\x0c\xe0\xe3\xf9\xfb\x2c\x9a\x2a\x38\x93\x76\x60\xcb\x80\x68\xf5\xe6\x07\xb6\x59\x10\x1a\x94\x76\x60\xfb\xae\xd3\xc6\x81\xdb\x0f\x19\xb1\x4e\xfa\xff\x16\xfb\x2b\xad\xdb\x15\xf7\xd7\xc2\xb9\xf2\x48\x62\x09\xce\xf4\xe8\xb5\xc6\xd4\x01\x04\x96\x7d\x5d\x4f\x57\x92\x97\xc0\x70\x01\xd4\x1a\x7c\x90\x60\x38\x1d\xa5\xd6\x75\x46\x57\xc3\x78\x9e\x5c\x02\x18\xdf\x2e\xa1\x62\xad\xc5\xe0\xdf\x01\x00\x00\xff\xff\x58\x7c\x2d\x36\xa2\x0d\x00\x00")

func goCentrifugeBuildConfigsDefault_configYamlBytes() ([]byte, error) {
	return bindataRead(
		_goCentrifugeBuildConfigsDefault_configYaml,
		"go-centrifuge/build/configs/default_config.yaml",
	)
}

func goCentrifugeBuildConfigsDefault_configYaml() (*asset, error) {
	bytes, err := goCentrifugeBuildConfigsDefault_configYamlBytes()
	if err != nil {
		return nil, err
	}

	info := bindataFileInfo{name: "go-centrifuge/build/configs/default_config.yaml", size: 3490, mode: os.FileMode(420), modTime: time.Unix(1543934586, 0)}
	a := &asset{bytes: bytes, info: info}
	return a, nil
}

var _goCentrifugeBuildConfigsTesting_configYaml = []byte("\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff\x94\x92\xc9\x6e\xe3\x48\x0c\x86\xef\x7a\x0a\x81\x97\x5c\xbc\xd4\xbe\xbd\xc1\x20\x98\xd3\x0c\x90\x33\xab\xc8\x8a\x05\xdb\xb2\x5a\x4b\x12\x23\xc8\xbb\x37\xe4\x38\x9d\x6b\x1a\xba\x90\x04\x7f\xfe\xa4\xea\xe3\xf9\xc0\x23\x2f\xe7\xd4\xb4\x2d\x96\x72\x59\xfa\x79\x5a\xe3\xb6\x3d\x63\xd7\xa7\xf6\x16\xb6\xed\x91\xaf\xa9\x7d\x78\x07\x24\x1a\x79\x9a\x20\x41\x88\x59\x60\x70\x36\xe8\x62\x8c\x31\x58\x2a\x79\x99\x8d\xd3\x2c\x48\x17\x6b\x91\xa5\x91\x0a\x2d\x6c\xa0\x8c\xd7\x61\xbe\x40\x7a\x87\xd2\x0d\x07\x1e\x21\x01\xf2\xb4\x95\x2a\x6c\xcb\x3c\xae\x0d\xb7\xf2\xcc\x6f\x33\x24\x28\xde\xc7\x1a\xb4\x8f\xe4\xbd\xa0\xa8\x4a\x2d\x92\x88\x0c\x86\xaa\x25\x59\x14\x48\x25\x54\x85\x22\x2b\x94\x46\x48\xed\x05\x69\xa7\x45\xd5\xa1\x88\x12\xf0\xcf\xbc\x01\x47\x3c\x4f\xab\x6d\xf7\x02\x09\xb4\x2b\xd2\x05\xf6\x3a\xd7\x18\x44\x65\x6f\xb3\xf0\xca\xd7\x10\x05\x7a\x89\x04\x1f\x1b\x38\x52\x85\x04\xd3\x6d\x61\xb8\xa5\xdf\x43\xe8\x78\xe2\x1e\x92\x56\x1b\xe8\x21\x29\xa7\xa4\x31\x1b\x18\x20\xc9\x0d\x8c\x90\xc2\x06\x26\x3c\xad\x07\x10\xcb\xcc\xd2\xb1\x2e\x31\xc8\x68\x0c\x49\x2e\xa8\x72\xc8\xca\xb3\x61\xc7\x22\xdb\x5c\xb3\xd1\x99\x85\xf6\x0e\x2d\x85\x10\x62\x45\xe7\x23\xaa\x20\x95\x5a\x17\x39\x63\x59\x7f\x45\x91\x2a\xe4\x20\xad\xb5\x36\xa3\x64\x24\x5f\x90\xa3\x70\x82\x43\x30\x0a\x6b\xc1\xa0\xad\x23\xe1\x8c\xb5\x99\x22\x5a\x6f\x55\x46\x57\x4b\x11\x51\x71\x5d\x27\x75\x04\x09\x8c\x65\xe1\x04\xba\x2d\x29\xe4\xad\xd1\x39\x6c\xa3\x52\x75\x6b\x4c\x50\xd1\xc4\x48\xda\x13\x6c\xe0\x85\xc7\xa9\xbb\xac\x47\x7e\x3c\xdc\x1f\x7e\xc0\x69\x7a\xbd\x8c\x94\xda\x87\xaf\xd2\x9d\x81\xd4\xfe\x14\x81\xa6\xe9\x88\xfb\xb9\x9b\xaf\xff\x50\x6a\x41\xbc\x09\xf9\xfd\x41\xd3\xfc\x5a\x78\xe1\x15\xba\x7e\x39\x3f\x5d\xc6\x23\x8f\x53\x6a\x55\xd3\xb6\xaf\xb7\xe4\x09\xbb\xf9\xff\xee\xcc\xff\xfe\x97\x5a\xd9\x34\x47\xbe\xde\x08\x9d\xba\xe7\xbe\xeb\x9f\x3f\x61\x1d\x96\x7c\xea\xca\xe3\x4a\xe9\x6e\xb7\xdf\xed\xf6\x79\xe9\x4e\xb4\x1f\x79\xba\x2c\x63\xe1\x69\x7f\xef\x7e\xe4\xeb\x6e\x58\xf2\x6e\xe0\xf3\xa7\x6e\xec\x5e\x70\xe6\x9f\x09\x8f\xab\xf8\x26\xe4\xf9\x80\xcb\x7c\xf8\xa1\xf7\xbd\xfb\x2f\x8d\xbf\x54\x5f\xae\xbf\x03\x00\x00\xff\xff\xb0\x1c\xaf\x3f\xaa\x03\x00\x00")

func goCentrifugeBuildConfigsTesting_configYamlBytes() ([]byte, error) {
	return bindataRead(
		_goCentrifugeBuildConfigsTesting_configYaml,
		"go-centrifuge/build/configs/testing_config.yaml",
	)
}

func goCentrifugeBuildConfigsTesting_configYaml() (*asset, error) {
	bytes, err := goCentrifugeBuildConfigsTesting_configYamlBytes()
	if err != nil {
		return nil, err
	}

	info := bindataFileInfo{name: "go-centrifuge/build/configs/testing_config.yaml", size: 938, mode: os.FileMode(420), modTime: time.Unix(1541606744, 0)}
	a := &asset{bytes: bytes, info: info}
	return a, nil
}

// Asset loads and returns the asset for the given name.
// It returns an error if the asset could not be found or
// could not be loaded.
func Asset(name string) ([]byte, error) {
	cannonicalName := strings.Replace(name, "\\", "/", -1)
	if f, ok := _bindata[cannonicalName]; ok {
		a, err := f()
		if err != nil {
			return nil, fmt.Errorf("Asset %s can't read by error: %v", name, err)
		}
		return a.bytes, nil
	}
	return nil, fmt.Errorf("Asset %s not found", name)
}

// MustAsset is like Asset but panics when Asset would return an error.
// It simplifies safe initialization of global variables.
func MustAsset(name string) []byte {
	a, err := Asset(name)
	if err != nil {
		panic("asset: Asset(" + name + "): " + err.Error())
	}

	return a
}

// AssetInfo loads and returns the asset info for the given name.
// It returns an error if the asset could not be found or
// could not be loaded.
func AssetInfo(name string) (os.FileInfo, error) {
	cannonicalName := strings.Replace(name, "\\", "/", -1)
	if f, ok := _bindata[cannonicalName]; ok {
		a, err := f()
		if err != nil {
			return nil, fmt.Errorf("AssetInfo %s can't read by error: %v", name, err)
		}
		return a.info, nil
	}
	return nil, fmt.Errorf("AssetInfo %s not found", name)
}

// AssetNames returns the names of the assets.
func AssetNames() []string {
	names := make([]string, 0, len(_bindata))
	for name := range _bindata {
		names = append(names, name)
	}
	return names
}

// _bindata is a table, holding each asset generator, mapped to its name.
var _bindata = map[string]func() (*asset, error){
	"go-centrifuge/build/configs/default_config.yaml": goCentrifugeBuildConfigsDefault_configYaml,
	"go-centrifuge/build/configs/testing_config.yaml": goCentrifugeBuildConfigsTesting_configYaml,
}

// AssetDir returns the file names below a certain
// directory embedded in the file by go-bindata.
// For example if you run go-bindata on data/... and data contains the
// following hierarchy:
//     data/
//       foo.txt
//       img/
//         a.png
//         b.png
// then AssetDir("data") would return []string{"foo.txt", "img"}
// AssetDir("data/img") would return []string{"a.png", "b.png"}
// AssetDir("foo.txt") and AssetDir("notexist") would return an error
// AssetDir("") will return []string{"data"}.
func AssetDir(name string) ([]string, error) {
	node := _bintree
	if len(name) != 0 {
		cannonicalName := strings.Replace(name, "\\", "/", -1)
		pathList := strings.Split(cannonicalName, "/")
		for _, p := range pathList {
			node = node.Children[p]
			if node == nil {
				return nil, fmt.Errorf("Asset %s not found", name)
			}
		}
	}
	if node.Func != nil {
		return nil, fmt.Errorf("Asset %s not found", name)
	}
	rv := make([]string, 0, len(node.Children))
	for childName := range node.Children {
		rv = append(rv, childName)
	}
	return rv, nil
}

type bintree struct {
	Func     func() (*asset, error)
	Children map[string]*bintree
}
var _bintree = &bintree{nil, map[string]*bintree{
	"go-centrifuge": &bintree{nil, map[string]*bintree{
		"build": &bintree{nil, map[string]*bintree{
			"configs": &bintree{nil, map[string]*bintree{
				"default_config.yaml": &bintree{goCentrifugeBuildConfigsDefault_configYaml, map[string]*bintree{}},
				"testing_config.yaml": &bintree{goCentrifugeBuildConfigsTesting_configYaml, map[string]*bintree{}},
			}},
		}},
	}},
}}

// RestoreAsset restores an asset under the given directory
func RestoreAsset(dir, name string) error {
	data, err := Asset(name)
	if err != nil {
		return err
	}
	info, err := AssetInfo(name)
	if err != nil {
		return err
	}
	err = os.MkdirAll(_filePath(dir, filepath.Dir(name)), os.FileMode(0755))
	if err != nil {
		return err
	}
	err = ioutil.WriteFile(_filePath(dir, name), data, info.Mode())
	if err != nil {
		return err
	}
	err = os.Chtimes(_filePath(dir, name), info.ModTime(), info.ModTime())
	if err != nil {
		return err
	}
	return nil
}

// RestoreAssets restores an asset under the given directory recursively
func RestoreAssets(dir, name string) error {
	children, err := AssetDir(name)
	// File
	if err != nil {
		return RestoreAsset(dir, name)
	}
	// Dir
	for _, child := range children {
		err = RestoreAssets(dir, filepath.Join(name, child))
		if err != nil {
			return err
		}
	}
	return nil
}

func _filePath(dir, name string) string {
	cannonicalName := strings.Replace(name, "\\", "/", -1)
	return filepath.Join(append([]string{dir}, strings.Split(cannonicalName, "/")...)...)
}

