#pragma once

#include "drake/common/drake_export.h"
#include "drake/systems/framework/leaf_system.h"
#include "drake/systems/plants/RigidBodyTree.h"
#include "drake/systems/plants/joints/DrakeJoint.h"

namespace drake {
namespace systems {

class DRAKE_EXPORT JointStateEncoder final : public LeafSystem<double> {
 public:
  using PositionMap = std::map<const DrakeJoint*, VectorX<double>>;
  using VelocityMap = std::map<const DrakeJoint*, VectorX<double>>;
  using EffortMap = std::map<const RigidBodyActuator*, double>;

  JointStateEncoder();

  ~JointStateEncoder() override {}

  // Disable copy and assign.
  JointStateEncoder(const JointStateEncoder&) = delete;

  JointStateEncoder& operator=(const JointStateEncoder&) = delete;

  void EvalOutput(const Context<double>& context,
                  SystemOutput<double>* output) const override;

  std::unique_ptr<SystemOutput<double>> AllocateOutput(
      const Context<double>& context) const override;

  const SystemPortDescriptor<double>& get_position_port() const;

  const SystemPortDescriptor<double>& get_velocity_port() const;

  const SystemPortDescriptor<double>& get_effort_port() const;

  const SystemPortDescriptor<double>& get_lcm_message_port() const;

 private:
  // Input ports.
  const SystemPortDescriptor<double>& position_port_;
  const SystemPortDescriptor<double>& velocity_port_;
  const SystemPortDescriptor<double>& effort_port_;

  // Output port.
  const systems::SystemPortDescriptor<double>& message_port_;
};

}  // systems
}  // drake
