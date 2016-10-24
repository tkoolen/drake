#pragma once

#include "drake/common/drake_export.h"
#include "drake/systems/framework/leaf_system.h"
#include "drake/systems/plants/RigidBodyTree.h"

namespace drake {
namespace systems {

// TODO(tkoolen): currently doesn't do anything with the effort part of the
// robot_state_t message.

/**
 * Converts a robot_state_t LCM message into a KinematicsCache object.
 */
class DRAKE_EXPORT RobotStateDecoder : public LeafSystem<double> {
 public:
  RobotStateDecoder(const RigidBodyTree& tree);

  ~RobotStateDecoder() override {}

  // Disable copy and assign.
  RobotStateDecoder(
      const RobotStateDecoder&) = delete;

  RobotStateDecoder& operator=(
      const RobotStateDecoder&) = delete;

  void EvalOutput(const Context<double>& context,
                  SystemOutput<double>* output) const override;

  std::unique_ptr<SystemOutput<double>> AllocateOutput(
      const Context<double>& context) const override;

 private:
  std::map<std::string, const RigidBody *> CreateJointNameToBodyMap(
      const RigidBodyTree &tree);

  const RigidBodyTree& tree_;
  const RigidBody* const floating_body_;
  const SystemPortDescriptor<double>& robot_state_message_port_;
  const SystemPortDescriptor<double>& kinematics_cache_port_;
  // TODO(tkoolen) should map to joint:
  const std::map<std::string, const RigidBody*> joint_name_to_body_;
};

}  // systems
}  // drake
