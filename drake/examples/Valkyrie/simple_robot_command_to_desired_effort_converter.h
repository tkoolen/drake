#pragma once

#include "drake/drakeRobotCommandToDesiredEffortConverter_export.h"
#include "drake/systems/framework/context.h"
#include "drake/systems/framework/leaf_system.h"
#include "drake/systems/plants/rigid_body_actuator.h"

// TODO(tkoolen): add comment: effort only.

// TODO(tkoolen): namespace
// TODO(tkoolen): export
namespace drake {

class DRAKEROBOTCOMMANDTODESIREDEFFORTCONVERTER_EXPORT
    SimpleRobotCommandToDesiredEffortConverter
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

  const systems::SystemPortDescriptor<double>& get_robot_command_port_() const;

  const systems::SystemPortDescriptor<double>& get_desired_effort_port_() const;

 private:
  const systems::SystemPortDescriptor<double>& robot_command_port_;
  const systems::SystemPortDescriptor<double>& desired_effort_port_;
  const std::vector<const RigidBodyActuator*> actuators_in_order_;

  static size_t CalcNumEfforts(
      const std::vector<const RigidBodyActuator*>& actuators);
};

}  // drake
