#pragma once

#include <lcm/lcm-cpp.hpp>

#include "drake/drakeRobotStatePublisher_export.h"
#include "drake/systems/framework/context.h"
#include "drake/systems/framework/leaf_system.h"
#include "drake/systems/plants/RigidbodyTree.h"
#include "drake/systems/robotInterfaces/Side.h"
#include "lcmtypes/bot_core/robot_state_t.hpp"

namespace drake {
namespace examples {  // TODO: right namespace?

class DRAKEROBOTSTATEPUBLISHER_EXPORT RobotStatePublisher
    : public systems::LeafSystem<double> {
  /**
   * Constructor.
   *
   * TODO
   */
  RobotStatePublisher(const RigidBodyTree& tree, int num_actuators,
                      const std::string& channel, lcm::LCM* lcm);

  ~RobotStatePublisher() override;

  // Disable copy and assign.
  RobotStatePublisher(const RobotStatePublisher&) = delete;

  RobotStatePublisher& operator=(const RobotStatePublisher&) = delete;

  std::string get_name() const override;

  /**
   * Takes the VectorBase from the input port of the context and publishes
   * it onto an LCM channel.
   */
  void DoPublish(const systems::ContextBase<double>& context) const override;

  /**
   * This System has no output ports so EvalOutput() does nothing.
   */
  void EvalOutput(const systems::ContextBase<double>& context,
                  systems::SystemOutput<double>* output) const override {}

 private:
  void InitializeMessage() const;

  Eigen::Isometry3d EvalFloatingBodyPose(
      const Eigen::Ref<const Eigen::VectorXd>& q) const;

  TwistVector<double> EvalFloatingBodyTwistInBodyFrame(
      const Eigen::Ref<const Eigen::VectorXd>& q,
      const Eigen::Ref<const Eigen::VectorXd>& v) const;

  TwistVector<double> TransformTwistFromBodyFrameToWorldAlignedBodyFrame(
      Eigen::Isometry3d& floating_body_to_world,
      const TwistVector<double>& twist_in_body) const;

  /// Number of floating joint positions.
  int num_floating_joint_positions() const;

  /// Number of floating joint velocities.
  int num_floating_joint_velocities() const;

  /// Check requirements on the tree, to ensure that robot_state_t can
  /// unambiguously represent its state.
  const RigidBodyTree& CheckPreConditions(const RigidBodyTree& tree);

  /// Tree to which message corresponds.
  const RigidBodyTree& tree_;

  /// Pointer to the body in @p tree_ that is attached to the world with a
  /// floating joint. Null if there is no such body.
  const RigidBody* const floating_body_;

  /// The channel on which to publish LCM messages.
  const std::string channel_;

  /// A pointer to the LCM subsystem.
  ::lcm::LCM* lcm_;

  /// LCM Message.
  mutable bot_core::robot_state_t message_;

  /// Data to send.
  mutable std::vector<uint8_t> message_bytes_;

  /// Input port indices.
  int state_port_index_;
  int torque_port_index_;
  std::map<Side, int> foot_wrench_port_indices_;
};

}  // examples
}  // drake
