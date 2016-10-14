#pragma once

#include <cstdint>
#include <memory>
#include <vector>

#include "drake/common/drake_assert.h"
#include "drake/common/drake_throw.h"
#include "drake/common/drake_export.h"
#include "drake/systems/framework/value.h"

namespace drake {
namespace systems {
namespace lcm {

class DRAKE_EXPORT SerializerInterface {
 public:
  SerializerInterface();
  virtual ~SerializerInterface();

  virtual std::unique_ptr<AbstractValue> CreateDefaultValue() const = 0;

  /**
   * Translates LCM message bytes into a `drake::systems::AbstractValue`
   * object.
   */
  virtual void Deserialize(
      const void* message_bytes, int message_length,
      AbstractValue* abstract_value) const = 0;

  /**
   * Translates a `drake::systems::AbstractValue` object into LCM message
   * bytes.
   */
  virtual void Serialize(const AbstractValue& abstract_value,
                         std::vector<uint8_t>* message_bytes) const = 0;

  // Disable copy and assign.
  SerializerInterface(const SerializerInterface&) = delete;
  SerializerInterface& operator=(const SerializerInterface&) = delete;
};

template <typename LcmtMessage>
class Serializer : public SerializerInterface {
 public:
  Serializer() {}

  std::unique_ptr<AbstractValue> CreateDefaultValue() const override {
    return std::make_unique<Value<LcmtMessage>>(LcmtMessage());
  }

  void Deserialize(
      const void* message_bytes, int message_length,
      AbstractValue* abstract_value) const override {
    DRAKE_DEMAND(abstract_value != nullptr);
    LcmtMessage& message = abstract_value->GetMutableValue<LcmtMessage>();
    int consumed = message.decode(message_bytes, 0, message_length);
    DRAKE_THROW_UNLESS(consumed == message_length);
  }

  void Serialize(const AbstractValue& abstract_value,
                 std::vector<uint8_t>* message_bytes) const override {
    DRAKE_DEMAND(message_bytes != nullptr);
    const LcmtMessage& message = abstract_value.GetValue<LcmtMessage>();
    const int message_length = message.getEncodedSize();
    message_bytes->resize(message_length);
    int consumed = message.encode(message_bytes->data(), 0, message_length);
    DRAKE_THROW_UNLESS(consumed == message_length);
  }

  // Disable copy and assign.
  Serializer(const Serializer&) = delete;
  Serializer& operator=(const Serializer&) = delete;
};

}  // namespace lcm
}  // namespace systems
}  // namespace drake
