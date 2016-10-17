#include "drake/examples/Valkyrie/robot_command_to_desired_effort_converter.h"
#include "lcmtypes/bot_core/atlas_command_t.hpp"

namespace drake {
namespace systems {

RobotCommandToDesiredEffortConverter::
    RobotCommandToDesiredEffortConverter(
        const std::vector<const RigidBodyActuator*>& actuators)
    : robot_command_port_(
          DeclareInputPort(kAbstractValued, 1, kContinuousSampling)),
      desired_effort_ports_(DeclareDesiredEffortOutputPorts(actuators)),
      name_to_actuator_(CreateNameToActuatorMap(actuators)) {
  set_name("RobotCommandToDesiredEffortConverter");
};

void RobotCommandToDesiredEffortConverter::EvalOutput(
    const systems::Context<double>& context,
    systems::SystemOutput<double>* output) const {
  using bot_core::atlas_command_t;

  const AbstractValue* input =
      EvalAbstractInput(context, robot_command_port_.get_index());
  const atlas_command_t& message = input->GetValue<atlas_command_t>();

  // Currently, all RigidBodyActuators are assumed to be one-dimensional.
  for (size_t i = 0; i < message.joint_names.size(); i++) {
    const std::string& joint_name = message.joint_names[i];
    const double& effort = message.effort[i];
    const RigidBodyActuator* actuator = name_to_actuator_.at(joint_name);
    int port_index = desired_effort_ports_.at(actuator)->get_index();
    output->get_mutable_port(port_index)
        ->GetMutableVectorData<double>()
        ->SetAtIndex(0, effort);
  }
}

const SystemPortDescriptor<double>&
RobotCommandToDesiredEffortConverter::get_desired_effort_output_port(
    const RigidBodyActuator& actuator) const {
  return *desired_effort_ports_.at(&actuator);
}

std::map<const RigidBodyActuator*, const SystemPortDescriptor<double>*>
RobotCommandToDesiredEffortConverter::DeclareDesiredEffortOutputPorts(
    const std::vector<const RigidBodyActuator*>& actuators) {
  // Currently, all RigidBodyActuators are assumed to be one-dimensional.
  const int desired_effort_length = 1;
  std::map<const RigidBodyActuator*, const SystemPortDescriptor<double>*> ret;
  for (const auto& actuator : actuators) {
    ret[actuator] = &DeclareOutputPort(kVectorValued, desired_effort_length,
                                       kContinuousSampling);
  }
  return ret;
}

std::map<std::string, const RigidBodyActuator*>
RobotCommandToDesiredEffortConverter::CreateNameToActuatorMap(
    const std::vector<const RigidBodyActuator*>& actuators) {
  std::map<std::string, const RigidBodyActuator*> ret;
  for (const auto& actuator : actuators) {
    ret[actuator->name_] = actuator;
  }
  return ret;
}

}  // systems
}  // drake
