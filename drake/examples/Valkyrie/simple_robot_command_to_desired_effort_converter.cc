#include "drake/examples/Valkyrie/simple_robot_command_to_desired_effort_converter.h"
#include "drake/examples/Valkyrie/robot_command.h"

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
      actuators_in_order_(actuators_in_order) {}

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
  using drake::systems::AbstractValue;

  const AbstractValue* input = EvalAbstractInput(context, robot_command_port_.get_index());
  const RobotCommand& command = input->GetValue<RobotCommand>();
  auto& efforts =
      *output->GetMutableVectorData(desired_effort_port_.get_index());

  // Currently, all RigidBodyActuators are assumed to be one-dimensional.
  int index = 0;
  for (const auto& actuator_ptr : actuators_in_order_) {
    double effort = command.at(actuator_ptr).effort_;
    efforts.SetAtIndex(index, effort);
    index++;
  }
}

const systems::SystemPortDescriptor<double>&
SimpleRobotCommandToDesiredEffortConverter::get_robot_command_port_() const {
  return robot_command_port_;
}

const systems::SystemPortDescriptor<double>&
SimpleRobotCommandToDesiredEffortConverter::get_desired_effort_port_() const {
  return desired_effort_port_;
}

size_t SimpleRobotCommandToDesiredEffortConverter::CalcNumEfforts(
    const std::vector<const RigidBodyActuator*>& actuators) {
  // Currently, all RigidBodyActuators are assumed to be one-dimensional.
  return actuators.size();
}

}  // drake
