#include "robot_command_to_desired_effort_converter.h"

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

void RobotCommandToDesiredEffortConverter::DoPublish(const systems::ContextBase<
    double> &context) const {
  // TODO
}

} // drake
