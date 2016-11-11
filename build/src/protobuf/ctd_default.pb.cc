// Generated by the protocol buffer compiler.  DO NOT EDIT!
// source: ctd_default.proto

#define INTERNAL_SUPPRESS_PROTOBUF_FIELD_DEPRECATION
#include "ctd_default.pb.h"

#include <algorithm>

#include <google/protobuf/stubs/common.h>
#include <google/protobuf/stubs/once.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/wire_format_lite_inl.h>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/generated_message_reflection.h>
#include <google/protobuf/reflection_ops.h>
#include <google/protobuf/wire_format.h>
// @@protoc_insertion_point(includes)

namespace {

const ::google::protobuf::Descriptor* CTDMessageDefault_descriptor_ = NULL;
const ::google::protobuf::internal::GeneratedMessageReflection*
  CTDMessageDefault_reflection_ = NULL;

}  // namespace


void protobuf_AssignDesc_ctd_5fdefault_2eproto() {
  protobuf_AddDesc_ctd_5fdefault_2eproto();
  const ::google::protobuf::FileDescriptor* file =
    ::google::protobuf::DescriptorPool::generated_pool()->FindFileByName(
      "ctd_default.proto");
  GOOGLE_CHECK(file != NULL);
  CTDMessageDefault_descriptor_ = file->message_type(0);
  static const int CTDMessageDefault_offsets_[3] = {
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(CTDMessageDefault, depth_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(CTDMessageDefault, temperature_),
    GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(CTDMessageDefault, salinity_),
  };
  CTDMessageDefault_reflection_ =
    new ::google::protobuf::internal::GeneratedMessageReflection(
      CTDMessageDefault_descriptor_,
      CTDMessageDefault::default_instance_,
      CTDMessageDefault_offsets_,
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(CTDMessageDefault, _has_bits_[0]),
      GOOGLE_PROTOBUF_GENERATED_MESSAGE_FIELD_OFFSET(CTDMessageDefault, _unknown_fields_),
      -1,
      ::google::protobuf::DescriptorPool::generated_pool(),
      ::google::protobuf::MessageFactory::generated_factory(),
      sizeof(CTDMessageDefault));
}

namespace {

GOOGLE_PROTOBUF_DECLARE_ONCE(protobuf_AssignDescriptors_once_);
inline void protobuf_AssignDescriptorsOnce() {
  ::google::protobuf::GoogleOnceInit(&protobuf_AssignDescriptors_once_,
                 &protobuf_AssignDesc_ctd_5fdefault_2eproto);
}

void protobuf_RegisterTypes(const ::std::string&) {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedMessage(
    CTDMessageDefault_descriptor_, &CTDMessageDefault::default_instance());
}

}  // namespace

void protobuf_ShutdownFile_ctd_5fdefault_2eproto() {
  delete CTDMessageDefault::default_instance_;
  delete CTDMessageDefault_reflection_;
}

void protobuf_AddDesc_ctd_5fdefault_2eproto() {
  static bool already_here = false;
  if (already_here) return;
  already_here = true;
  GOOGLE_PROTOBUF_VERIFY_VERSION;

  ::dccl::protobuf_AddDesc_goby_2facomms_2fprotobuf_2fdccl_5foption_5fextensions_2eproto();
  ::google::protobuf::DescriptorPool::InternalAddGeneratedFile(
    "\n\021ctd_default.proto\0321goby/acomms/protobu"
    "f/dccl_option_extensions.proto\"\267\001\n\021CTDMe"
    "ssageDefault\022,\n\005depth\030\001 \003(\005B\035\242\?\t1\000\000\000\000\000@\217"
    "@\242\?\t)\000\000\000\000\000\000\000\000\242\?\002P\005\0222\n\013temperature\030\002 \003(\005B"
    "\035\242\?\t1\000\000\000\000\000\0004@\242\?\t)\000\000\000\000\000\000$@\242\?\002P\005\0224\n\010salini"
    "ty\030\003 \003(\001B\"\242\?\t1\000\000\000\000\000\000D@\242\?\t)\000\000\000\000\000\0009@\242\?\002 \002\242"
    "\?\002P\005:\n\242\?\002\010f\242\?\002\020 ", 256);
  ::google::protobuf::MessageFactory::InternalRegisterGeneratedFile(
    "ctd_default.proto", &protobuf_RegisterTypes);
  CTDMessageDefault::default_instance_ = new CTDMessageDefault();
  CTDMessageDefault::default_instance_->InitAsDefaultInstance();
  ::google::protobuf::internal::OnShutdown(&protobuf_ShutdownFile_ctd_5fdefault_2eproto);
}

// Force AddDescriptors() to be called at static initialization time.
struct StaticDescriptorInitializer_ctd_5fdefault_2eproto {
  StaticDescriptorInitializer_ctd_5fdefault_2eproto() {
    protobuf_AddDesc_ctd_5fdefault_2eproto();
  }
} static_descriptor_initializer_ctd_5fdefault_2eproto_;

// ===================================================================

#ifndef _MSC_VER
const int CTDMessageDefault::kDepthFieldNumber;
const int CTDMessageDefault::kTemperatureFieldNumber;
const int CTDMessageDefault::kSalinityFieldNumber;
#endif  // !_MSC_VER

CTDMessageDefault::CTDMessageDefault()
  : ::google::protobuf::Message() {
  SharedCtor();
}

void CTDMessageDefault::InitAsDefaultInstance() {
}

CTDMessageDefault::CTDMessageDefault(const CTDMessageDefault& from)
  : ::google::protobuf::Message() {
  SharedCtor();
  MergeFrom(from);
}

void CTDMessageDefault::SharedCtor() {
  _cached_size_ = 0;
  ::memset(_has_bits_, 0, sizeof(_has_bits_));
}

CTDMessageDefault::~CTDMessageDefault() {
  SharedDtor();
}

void CTDMessageDefault::SharedDtor() {
  if (this != default_instance_) {
  }
}

void CTDMessageDefault::SetCachedSize(int size) const {
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
}
const ::google::protobuf::Descriptor* CTDMessageDefault::descriptor() {
  protobuf_AssignDescriptorsOnce();
  return CTDMessageDefault_descriptor_;
}

const CTDMessageDefault& CTDMessageDefault::default_instance() {
  if (default_instance_ == NULL) protobuf_AddDesc_ctd_5fdefault_2eproto();
  return *default_instance_;
}

CTDMessageDefault* CTDMessageDefault::default_instance_ = NULL;

CTDMessageDefault* CTDMessageDefault::New() const {
  return new CTDMessageDefault;
}

void CTDMessageDefault::Clear() {
  depth_.Clear();
  temperature_.Clear();
  salinity_.Clear();
  ::memset(_has_bits_, 0, sizeof(_has_bits_));
  mutable_unknown_fields()->Clear();
}

bool CTDMessageDefault::MergePartialFromCodedStream(
    ::google::protobuf::io::CodedInputStream* input) {
#define DO_(EXPRESSION) if (!(EXPRESSION)) return false
  ::google::protobuf::uint32 tag;
  while ((tag = input->ReadTag()) != 0) {
    switch (::google::protobuf::internal::WireFormatLite::GetTagFieldNumber(tag)) {
      // repeated int32 depth = 1;
      case 1: {
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_VARINT) {
         parse_depth:
          DO_((::google::protobuf::internal::WireFormatLite::ReadRepeatedPrimitive<
                   ::google::protobuf::int32, ::google::protobuf::internal::WireFormatLite::TYPE_INT32>(
                 1, 8, input, this->mutable_depth())));
        } else if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag)
                   == ::google::protobuf::internal::WireFormatLite::
                      WIRETYPE_LENGTH_DELIMITED) {
          DO_((::google::protobuf::internal::WireFormatLite::ReadPackedPrimitiveNoInline<
                   ::google::protobuf::int32, ::google::protobuf::internal::WireFormatLite::TYPE_INT32>(
                 input, this->mutable_depth())));
        } else {
          goto handle_uninterpreted;
        }
        if (input->ExpectTag(8)) goto parse_depth;
        if (input->ExpectTag(16)) goto parse_temperature;
        break;
      }

      // repeated int32 temperature = 2;
      case 2: {
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_VARINT) {
         parse_temperature:
          DO_((::google::protobuf::internal::WireFormatLite::ReadRepeatedPrimitive<
                   ::google::protobuf::int32, ::google::protobuf::internal::WireFormatLite::TYPE_INT32>(
                 1, 16, input, this->mutable_temperature())));
        } else if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag)
                   == ::google::protobuf::internal::WireFormatLite::
                      WIRETYPE_LENGTH_DELIMITED) {
          DO_((::google::protobuf::internal::WireFormatLite::ReadPackedPrimitiveNoInline<
                   ::google::protobuf::int32, ::google::protobuf::internal::WireFormatLite::TYPE_INT32>(
                 input, this->mutable_temperature())));
        } else {
          goto handle_uninterpreted;
        }
        if (input->ExpectTag(16)) goto parse_temperature;
        if (input->ExpectTag(25)) goto parse_salinity;
        break;
      }

      // repeated double salinity = 3;
      case 3: {
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_FIXED64) {
         parse_salinity:
          DO_((::google::protobuf::internal::WireFormatLite::ReadRepeatedPrimitive<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 1, 25, input, this->mutable_salinity())));
        } else if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag)
                   == ::google::protobuf::internal::WireFormatLite::
                      WIRETYPE_LENGTH_DELIMITED) {
          DO_((::google::protobuf::internal::WireFormatLite::ReadPackedPrimitiveNoInline<
                   double, ::google::protobuf::internal::WireFormatLite::TYPE_DOUBLE>(
                 input, this->mutable_salinity())));
        } else {
          goto handle_uninterpreted;
        }
        if (input->ExpectTag(25)) goto parse_salinity;
        if (input->ExpectAtEnd()) return true;
        break;
      }

      default: {
      handle_uninterpreted:
        if (::google::protobuf::internal::WireFormatLite::GetTagWireType(tag) ==
            ::google::protobuf::internal::WireFormatLite::WIRETYPE_END_GROUP) {
          return true;
        }
        DO_(::google::protobuf::internal::WireFormat::SkipField(
              input, tag, mutable_unknown_fields()));
        break;
      }
    }
  }
  return true;
#undef DO_
}

void CTDMessageDefault::SerializeWithCachedSizes(
    ::google::protobuf::io::CodedOutputStream* output) const {
  // repeated int32 depth = 1;
  for (int i = 0; i < this->depth_size(); i++) {
    ::google::protobuf::internal::WireFormatLite::WriteInt32(
      1, this->depth(i), output);
  }

  // repeated int32 temperature = 2;
  for (int i = 0; i < this->temperature_size(); i++) {
    ::google::protobuf::internal::WireFormatLite::WriteInt32(
      2, this->temperature(i), output);
  }

  // repeated double salinity = 3;
  for (int i = 0; i < this->salinity_size(); i++) {
    ::google::protobuf::internal::WireFormatLite::WriteDouble(
      3, this->salinity(i), output);
  }

  if (!unknown_fields().empty()) {
    ::google::protobuf::internal::WireFormat::SerializeUnknownFields(
        unknown_fields(), output);
  }
}

::google::protobuf::uint8* CTDMessageDefault::SerializeWithCachedSizesToArray(
    ::google::protobuf::uint8* target) const {
  // repeated int32 depth = 1;
  for (int i = 0; i < this->depth_size(); i++) {
    target = ::google::protobuf::internal::WireFormatLite::
      WriteInt32ToArray(1, this->depth(i), target);
  }

  // repeated int32 temperature = 2;
  for (int i = 0; i < this->temperature_size(); i++) {
    target = ::google::protobuf::internal::WireFormatLite::
      WriteInt32ToArray(2, this->temperature(i), target);
  }

  // repeated double salinity = 3;
  for (int i = 0; i < this->salinity_size(); i++) {
    target = ::google::protobuf::internal::WireFormatLite::
      WriteDoubleToArray(3, this->salinity(i), target);
  }

  if (!unknown_fields().empty()) {
    target = ::google::protobuf::internal::WireFormat::SerializeUnknownFieldsToArray(
        unknown_fields(), target);
  }
  return target;
}

int CTDMessageDefault::ByteSize() const {
  int total_size = 0;

  // repeated int32 depth = 1;
  {
    int data_size = 0;
    for (int i = 0; i < this->depth_size(); i++) {
      data_size += ::google::protobuf::internal::WireFormatLite::
        Int32Size(this->depth(i));
    }
    total_size += 1 * this->depth_size() + data_size;
  }

  // repeated int32 temperature = 2;
  {
    int data_size = 0;
    for (int i = 0; i < this->temperature_size(); i++) {
      data_size += ::google::protobuf::internal::WireFormatLite::
        Int32Size(this->temperature(i));
    }
    total_size += 1 * this->temperature_size() + data_size;
  }

  // repeated double salinity = 3;
  {
    int data_size = 0;
    data_size = 8 * this->salinity_size();
    total_size += 1 * this->salinity_size() + data_size;
  }

  if (!unknown_fields().empty()) {
    total_size +=
      ::google::protobuf::internal::WireFormat::ComputeUnknownFieldsSize(
        unknown_fields());
  }
  GOOGLE_SAFE_CONCURRENT_WRITES_BEGIN();
  _cached_size_ = total_size;
  GOOGLE_SAFE_CONCURRENT_WRITES_END();
  return total_size;
}

void CTDMessageDefault::MergeFrom(const ::google::protobuf::Message& from) {
  GOOGLE_CHECK_NE(&from, this);
  const CTDMessageDefault* source =
    ::google::protobuf::internal::dynamic_cast_if_available<const CTDMessageDefault*>(
      &from);
  if (source == NULL) {
    ::google::protobuf::internal::ReflectionOps::Merge(from, this);
  } else {
    MergeFrom(*source);
  }
}

void CTDMessageDefault::MergeFrom(const CTDMessageDefault& from) {
  GOOGLE_CHECK_NE(&from, this);
  depth_.MergeFrom(from.depth_);
  temperature_.MergeFrom(from.temperature_);
  salinity_.MergeFrom(from.salinity_);
  mutable_unknown_fields()->MergeFrom(from.unknown_fields());
}

void CTDMessageDefault::CopyFrom(const ::google::protobuf::Message& from) {
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

void CTDMessageDefault::CopyFrom(const CTDMessageDefault& from) {
  if (&from == this) return;
  Clear();
  MergeFrom(from);
}

bool CTDMessageDefault::IsInitialized() const {

  return true;
}

void CTDMessageDefault::Swap(CTDMessageDefault* other) {
  if (other != this) {
    depth_.Swap(&other->depth_);
    temperature_.Swap(&other->temperature_);
    salinity_.Swap(&other->salinity_);
    std::swap(_has_bits_[0], other->_has_bits_[0]);
    _unknown_fields_.Swap(&other->_unknown_fields_);
    std::swap(_cached_size_, other->_cached_size_);
  }
}

::google::protobuf::Metadata CTDMessageDefault::GetMetadata() const {
  protobuf_AssignDescriptorsOnce();
  ::google::protobuf::Metadata metadata;
  metadata.descriptor = CTDMessageDefault_descriptor_;
  metadata.reflection = CTDMessageDefault_reflection_;
  return metadata;
}


// @@protoc_insertion_point(namespace_scope)

// @@protoc_insertion_point(global_scope)
