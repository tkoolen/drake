#include "drake/examples/Valkyrie/robot_command_to_desired_effort_converter.h"
#include "drake/examples/Valkyrie/robot_command.h"

namespace drake {

using drake::systems::System;

RobotCommandToDesiredEffortConverter::RobotCommandToDesiredEffortConverter() {

}

RobotCommandToDesiredEffortConverter::~RobotCommandToDesiredEffortConverter() {
  // empty
}

std::string RobotCommandToDesiredEffortConverter::get_name() const {
  // TODO
  return System::get_name();
}

void RobotCommandToDesiredEffortConverter::EvalOutput(const systems::Context<double> &context,
                                                 systems::SystemOutput<double> *output) const {
  using drake::systems::AbstractValue;

  AbstractValue* data = output->GetMutableData(0);
  RobotCommand& command = data->GetMutableValue<RobotCommand>();


}

} // drake
