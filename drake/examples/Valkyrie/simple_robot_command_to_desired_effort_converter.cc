#include "drake/examples/Valkyrie/simple_robot_command_to_desired_effort_converter.h"
#include "lcmtypes/bot_core/atlas_command_t.hpp"

namespace drake {

using drake::systems::System;
using drake::systems::kAbstractValued;
using drake::systems::kVectorValued;
using drake::systems::kContinuousSampling;

SimpleRobotCommandToDesiredEffortConverter::
    SimpleRobotCommandToDesiredEffortConverter(
        const std::vector<const RigidBodyActuator*>& actuators_in_order)
    : robot_command_port_(
          DeclareInputPort(kAbstractValued, 1, kContinuousSampling)),
      desired_effort_port_(DeclareOutputPort(
          kVectorValued, static_cast<int>(CalcNumEfforts(actuators_in_order)),
          kContinuousSampling)),
      name_to_index_(CreateNameToIndexMap(actuators_in_order)) {};

SimpleRobotCommandToDesiredEffortConverter::
    ~SimpleRobotCommandToDesiredEffortConverter() {
  // empty
}

std::string SimpleRobotCommandToDesiredEffortConverter::get_name() const {
  // TODO
  return System::get_name();
}

void SimpleRobotCommandToDesiredEffortConverter::EvalOutput(
    const systems::Context<double>& context,
    systems::SystemOutput<double>* output) const {
  using systems::AbstractValue;
  using bot_core::atlas_command_t;

  const AbstractValue* input = EvalAbstractInput(context, robot_command_port_.get_index());
  const atlas_command_t& message = input->GetValue<atlas_command_t>();
  auto& efforts =
      *output->GetMutableVectorData(desired_effort_port_.get_index());

  // Currently, all RigidBodyActuators are assumed to be one-dimensional.
  for (size_t i = 0; i < message.joint_names.size(); i++) {
    const auto& joint_name = message.joint_names[i];
    const auto& effort = message.effort[i];
    efforts.SetAtIndex(name_to_index_.at(joint_name), effort);
  }
}

size_t SimpleRobotCommandToDesiredEffortConverter::CalcNumEfforts(
    const std::vector<const RigidBodyActuator*>& actuators) {
  // Currently, all RigidBodyActuators are assumed to be one-dimensional.
  return actuators.size();
}

std::map<std::string, int>
SimpleRobotCommandToDesiredEffortConverter::CreateNameToIndexMap(const std::vector<
    const RigidBodyActuator *> &actuators_in_order) {
  std::map<std::string, int> ret;
  int index = 0;
  for (const auto& actuator_ptr : actuators_in_order) {
    ret[actuator_ptr->name_] = index++;
  }
  return ret;
}

}  // drake
