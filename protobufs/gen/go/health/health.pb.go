// Code generated by protoc-gen-go. DO NOT EDIT.
// source: health/health.proto

package healthpb

import proto "github.com/golang/protobuf/proto"
import fmt "fmt"
import math "math"
import empty "github.com/golang/protobuf/ptypes/empty"
import _ "github.com/grpc-ecosystem/grpc-gateway/protoc-gen-swagger/options"
import _ "google.golang.org/genproto/googleapis/api/annotations"

import (
	context "golang.org/x/net/context"
	grpc "google.golang.org/grpc"
)

// Reference imports to suppress errors if they are not otherwise used.
var _ = proto.Marshal
var _ = fmt.Errorf
var _ = math.Inf

// This is a compile-time assertion to ensure that this generated file
// is compatible with the proto package it is being compiled against.
// A compilation error at this line likely means your copy of the
// proto package needs to be updated.
const _ = proto.ProtoPackageIsVersion2 // please upgrade the proto package

// Pong contains basic information about the node
type Pong struct {
	Version              string   `protobuf:"bytes,1,opt,name=version,proto3" json:"version,omitempty"`
	Network              string   `protobuf:"bytes,2,opt,name=network,proto3" json:"network,omitempty"`
	XXX_NoUnkeyedLiteral struct{} `json:"-"`
	XXX_unrecognized     []byte   `json:"-"`
	XXX_sizecache        int32    `json:"-"`
}

func (m *Pong) Reset()         { *m = Pong{} }
func (m *Pong) String() string { return proto.CompactTextString(m) }
func (*Pong) ProtoMessage()    {}
func (*Pong) Descriptor() ([]byte, []int) {
	return fileDescriptor_health_a9f2c3a17a667e6a, []int{0}
}
func (m *Pong) XXX_Unmarshal(b []byte) error {
	return xxx_messageInfo_Pong.Unmarshal(m, b)
}
func (m *Pong) XXX_Marshal(b []byte, deterministic bool) ([]byte, error) {
	return xxx_messageInfo_Pong.Marshal(b, m, deterministic)
}
func (dst *Pong) XXX_Merge(src proto.Message) {
	xxx_messageInfo_Pong.Merge(dst, src)
}
func (m *Pong) XXX_Size() int {
	return xxx_messageInfo_Pong.Size(m)
}
func (m *Pong) XXX_DiscardUnknown() {
	xxx_messageInfo_Pong.DiscardUnknown(m)
}

var xxx_messageInfo_Pong proto.InternalMessageInfo

func (m *Pong) GetVersion() string {
	if m != nil {
		return m.Version
	}
	return ""
}

func (m *Pong) GetNetwork() string {
	if m != nil {
		return m.Network
	}
	return ""
}

func init() {
	proto.RegisterType((*Pong)(nil), "health.Pong")
}

// Reference imports to suppress errors if they are not otherwise used.
var _ context.Context
var _ grpc.ClientConn

// This is a compile-time assertion to ensure that this generated file
// is compatible with the grpc package it is being compiled against.
const _ = grpc.SupportPackageIsVersion4

// HealthCheckServiceClient is the client API for HealthCheckService service.
//
// For semantics around ctx use and closing/ending streaming RPCs, please refer to https://godoc.org/google.golang.org/grpc#ClientConn.NewStream.
type HealthCheckServiceClient interface {
	Ping(ctx context.Context, in *empty.Empty, opts ...grpc.CallOption) (*Pong, error)
}

type healthCheckServiceClient struct {
	cc *grpc.ClientConn
}

func NewHealthCheckServiceClient(cc *grpc.ClientConn) HealthCheckServiceClient {
	return &healthCheckServiceClient{cc}
}

func (c *healthCheckServiceClient) Ping(ctx context.Context, in *empty.Empty, opts ...grpc.CallOption) (*Pong, error) {
	out := new(Pong)
	err := c.cc.Invoke(ctx, "/health.HealthCheckService/Ping", in, out, opts...)
	if err != nil {
		return nil, err
	}
	return out, nil
}

// HealthCheckServiceServer is the server API for HealthCheckService service.
type HealthCheckServiceServer interface {
	Ping(context.Context, *empty.Empty) (*Pong, error)
}

func RegisterHealthCheckServiceServer(s *grpc.Server, srv HealthCheckServiceServer) {
	s.RegisterService(&_HealthCheckService_serviceDesc, srv)
}

func _HealthCheckService_Ping_Handler(srv interface{}, ctx context.Context, dec func(interface{}) error, interceptor grpc.UnaryServerInterceptor) (interface{}, error) {
	in := new(empty.Empty)
	if err := dec(in); err != nil {
		return nil, err
	}
	if interceptor == nil {
		return srv.(HealthCheckServiceServer).Ping(ctx, in)
	}
	info := &grpc.UnaryServerInfo{
		Server:     srv,
		FullMethod: "/health.HealthCheckService/Ping",
	}
	handler := func(ctx context.Context, req interface{}) (interface{}, error) {
		return srv.(HealthCheckServiceServer).Ping(ctx, req.(*empty.Empty))
	}
	return interceptor(ctx, in, info, handler)
}

var _HealthCheckService_serviceDesc = grpc.ServiceDesc{
	ServiceName: "health.HealthCheckService",
	HandlerType: (*HealthCheckServiceServer)(nil),
	Methods: []grpc.MethodDesc{
		{
			MethodName: "Ping",
			Handler:    _HealthCheckService_Ping_Handler,
		},
	},
	Streams:  []grpc.StreamDesc{},
	Metadata: "health/health.proto",
}

func init() { proto.RegisterFile("health/health.proto", fileDescriptor_health_a9f2c3a17a667e6a) }

var fileDescriptor_health_a9f2c3a17a667e6a = []byte{
	// 266 bytes of a gzipped FileDescriptorProto
	0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x02, 0xff, 0x64, 0x8f, 0xcf, 0x4a, 0x33, 0x31,
	0x14, 0xc5, 0x99, 0xd2, 0xaf, 0xfd, 0x8c, 0xae, 0x22, 0xc8, 0x38, 0x75, 0x31, 0x14, 0x04, 0x41,
	0x9b, 0x80, 0xee, 0xdc, 0x59, 0x11, 0x5c, 0xc9, 0xa0, 0x2b, 0xdd, 0x65, 0xe2, 0x6d, 0x26, 0xb4,
	0xcd, 0x0d, 0x99, 0xd8, 0xe2, 0xd6, 0x47, 0xd0, 0x47, 0xf3, 0x15, 0x7c, 0x10, 0xc9, 0x9f, 0xae,
	0x5c, 0x5d, 0x4e, 0x4e, 0xce, 0xb9, 0xbf, 0x4b, 0x0e, 0x3b, 0x10, 0x2b, 0xdf, 0xf1, 0x34, 0x98,
	0x75, 0xe8, 0x91, 0x8e, 0x92, 0xaa, 0x4e, 0x14, 0xa2, 0x5a, 0x01, 0x17, 0x56, 0x73, 0x61, 0x0c,
	0x7a, 0xe1, 0x35, 0x9a, 0x3e, 0xfd, 0xaa, 0x26, 0xd9, 0x8d, 0xaa, 0x7d, 0x5b, 0x70, 0x58, 0x5b,
	0xff, 0x9e, 0xcd, 0x8b, 0x38, 0xe4, 0x4c, 0x81, 0x99, 0xf5, 0x5b, 0xa1, 0x14, 0x38, 0x8e, 0x36,
	0xc6, 0xff, 0x56, 0x4d, 0xaf, 0xc9, 0xb0, 0x41, 0xa3, 0x68, 0x49, 0xc6, 0x1b, 0x70, 0xbd, 0x46,
	0x53, 0x16, 0x75, 0x71, 0xb6, 0xf7, 0xb8, 0x93, 0xc1, 0x31, 0xe0, 0xb7, 0xe8, 0x96, 0xe5, 0x20,
	0x39, 0x59, 0x5e, 0x22, 0xa1, 0xf7, 0x11, 0xf7, 0xb6, 0x03, 0xb9, 0x7c, 0x02, 0xb7, 0xd1, 0x12,
	0xe8, 0x33, 0x19, 0x36, 0xda, 0x28, 0x7a, 0xc4, 0x12, 0x25, 0xdb, 0x51, 0xb2, 0xbb, 0x40, 0x59,
	0x1d, 0xb0, 0x7c, 0x71, 0xd8, 0x3b, 0x3d, 0xff, 0xbc, 0x99, 0x54, 0xc7, 0xa9, 0xa6, 0x96, 0xa1,
	0xa7, 0x5e, 0xa0, 0xab, 0x7d, 0x07, 0xf5, 0x03, 0xbe, 0xc2, 0xc7, 0xf7, 0xcf, 0xd7, 0x60, 0x4c,
	0xff, 0x71, 0xab, 0x8d, 0x9a, 0x9f, 0x12, 0x22, 0x71, 0x9d, 0xf3, 0xf3, 0xfd, 0x94, 0x6a, 0x42,
	0x7d, 0x53, 0xbc, 0xfc, 0x4f, 0xcf, 0xb6, 0x6d, 0x47, 0x71, 0xe3, 0xd5, 0x6f, 0x00, 0x00, 0x00,
	0xff, 0xff, 0x6a, 0xd9, 0xf3, 0x6e, 0x62, 0x01, 0x00, 0x00,
}
