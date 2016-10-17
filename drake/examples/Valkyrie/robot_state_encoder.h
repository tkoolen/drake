#pragma once

#include "lcmtypes/bot_core/robot_state_t.hpp"

#include "drake/common/drake_export.h"
#include "drake/systems/framework/context.h"
#include "drake/systems/framework/leaf_system.h"
#include "drake/systems/plants/RigidBodyTree.h"
#include "drake/systems/robotInterfaces/Side.h"

namespace drake {
namespace systems {

// TODO: allocate output
class DRAKE_EXPORT RobotStateEncoder final : public LeafSystem<double> {
 public:
  RobotStateEncoder(const RigidBodyTree& tree);

  ~RobotStateEncoder() override;

  // Disable copy and assign.
  RobotStateEncoder(const RobotStateEncoder&) = delete;

  RobotStateEncoder& operator=(const RobotStateEncoder&) = delete;

  void EvalOutput(const Context<double>& context,
                  SystemOutput<double>* output) const override;

  std::unique_ptr<SystemOutput<double>> AllocateOutput(
      const Context<double>& context) const override;

  const SystemPortDescriptor<double>& get_lcm_message_port() const;

  const SystemPortDescriptor<double>& get_kinematics_results_port() const;

  const SystemPortDescriptor<double>& get_effort_port(
      const RigidBodyActuator& actuator) const;

  const SystemPortDescriptor<double>& get_foot_contact_wrench_port(
      const Side& side) const;

  const SystemPortDescriptor<double>& get_hand_contact_wrench_port(
      const Side& side) const;

 private:
  std::map<const RigidBodyActuator*, const SystemPortDescriptor<double>*>
  DeclareEffortInputPorts();

  std::map<Side, const SystemPortDescriptor<double>*> DeclareWrenchInputPorts();

  void SetStateAndEfforts(bot_core::robot_state_t& message,
                          const Context<double>& context) const;

  void SetForceTorque(bot_core::robot_state_t& message,
                      const Context<double>& context) const;

  // Tree to which message corresponds.
  const RigidBodyTree& tree_;

  // Pointer to the body in @p tree_ that is attached to the world with a
  // floating joint. Null if there is no such body.
  const RigidBody* const floating_body_;

  // Output port.
  const SystemPortDescriptor<double>& lcm_message_port_;

  // Input ports.
  const SystemPortDescriptor<double>& kinematics_results_port_;
  const std::map<const RigidBodyActuator*, const SystemPortDescriptor<double>*>
      effort_ports_;
  const std::map<Side, const SystemPortDescriptor<double>*> foot_wrench_ports_;
  const std::map<Side, const SystemPortDescriptor<double>*> hand_wrench_ports_;
};

}  // systems
}  // drake
