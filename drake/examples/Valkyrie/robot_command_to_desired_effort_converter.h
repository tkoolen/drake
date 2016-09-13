#pragma once

#include "drake/drakeRobotCommandToDesiredEffortConverter_export.h"
#include "drake/systems/framework/context.h"
#include "drake/systems/framework/leaf_system.h"

// credit to Gregory Izatt for the original implementation of most of this in
// LCM2RosControl

// TODO(tkoolen): namespace
// TODO(tkoolen): export
namespace drake {

class DRAKEROBOTCOMMANDTODESIREDEFFORTCONVERTER_EXPORT
    RobotCommandToDesiredEffortConverter : public systems::LeafSystem<double> {
 public:
  RobotCommandToDesiredEffortConverter();

  ~RobotCommandToDesiredEffortConverter() override;

  // Disable copy and assign.
  RobotCommandToDesiredEffortConverter(
      const RobotCommandToDesiredEffortConverter&) = delete;

  RobotCommandToDesiredEffortConverter& operator=(
      const RobotCommandToDesiredEffortConverter&) = delete;

  std::string get_name() const override;

  void EvalOutput(const systems::Context<double>& context,
                  systems::SystemOutput<double>* output) const override;
};

}  // drake
