#pragma once

#include <map>

#include "drake/common/drake_export.h"
#include "drake/systems/framework/context.h"
#include "drake/systems/framework/leaf_system.h"
#include "drake/systems/plants/rigid_body_actuator.h"

namespace drake {
namespace systems {

/**
 * Converts an atlas_command_t message into desired efforts, presented on one
 * output port per actuator.
 * Currently just takes the effort part of the
 * atlas_command_t message, without considering the PID control parts.
 */
class DRAKE_EXPORT RobotCommandToDesiredEffortConverter
    : public LeafSystem<double> {
 public:
  RobotCommandToDesiredEffortConverter(
      const std::vector<const RigidBodyActuator*>& actuators);

  ~RobotCommandToDesiredEffortConverter() override {}

  // Disable copy and assign.
  RobotCommandToDesiredEffortConverter(
      const RobotCommandToDesiredEffortConverter&) = delete;

  RobotCommandToDesiredEffortConverter& operator=(
      const RobotCommandToDesiredEffortConverter&) = delete;

  void EvalOutput(const Context<double>& context,
                  SystemOutput<double>* output) const override;

  const SystemPortDescriptor<double>& get_desired_effort_output_port(
      const RigidBodyActuator& actuator) const;

 private:
  const SystemPortDescriptor<double>& robot_command_port_;
  const std::map<const RigidBodyActuator*, const SystemPortDescriptor<double>*>
      desired_effort_ports_;
  const std::map<std::string, const RigidBodyActuator*> name_to_actuator_;

  /// Declare one output port for each RigidBodyActuator and store their
  /// descriptors in a map.
  std::map<const RigidBodyActuator*, const SystemPortDescriptor<double>*>
  DeclareDesiredEffortOutputPorts(
      const std::vector<const RigidBodyActuator*>& actuators);

  /// Map from actuator name to actuator.
  std::map<std::string, const RigidBodyActuator*> CreateNameToActuatorMap(
      const std::vector<const RigidBodyActuator*>& actuators);
};

}  // systems
}  // drake
