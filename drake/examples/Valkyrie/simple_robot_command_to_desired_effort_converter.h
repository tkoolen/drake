#pragma once

#include <map>

#include "drake/systems/framework/context.h"
#include "drake/systems/framework/leaf_system.h"
#include "drake/systems/plants/rigid_body_actuator.h"

// TODO(tkoolen): add comment: effort only.

// TODO(tkoolen): namespace
// TODO(tkoolen): export
namespace drake {

class SimpleRobotCommandToDesiredEffortConverter
    : public systems::LeafSystem<double> {
 public:
  SimpleRobotCommandToDesiredEffortConverter(
      const std::vector<const RigidBodyActuator*>& actuators_in_order);

  ~SimpleRobotCommandToDesiredEffortConverter() override;

  // Disable copy and assign.
  SimpleRobotCommandToDesiredEffortConverter(
      const SimpleRobotCommandToDesiredEffortConverter&) = delete;

  SimpleRobotCommandToDesiredEffortConverter& operator=(
      const SimpleRobotCommandToDesiredEffortConverter&) = delete;

  std::string get_name() const override;

  void EvalOutput(const systems::Context<double>& context,
                  systems::SystemOutput<double>* output) const override;

 private:
  const systems::SystemPortDescriptor<double>& robot_command_port_;
  const systems::SystemPortDescriptor<double>& desired_effort_port_;
  const std::map<std::string, int> name_to_index_;

  static size_t CalcNumEfforts(const std::vector<const RigidBodyActuator*>& actuators);
  static std::map<std::string, int> CreateNameToIndexMap(const std::vector<const RigidBodyActuator*>& actuators_in_order);
};

}  // drake
